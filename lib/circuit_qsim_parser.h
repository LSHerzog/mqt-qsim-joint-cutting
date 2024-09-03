// Copyright 2019 Google LLC. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef CIRCUIT_QSIM_PARSER_H_
#define CIRCUIT_QSIM_PARSER_H_

#include <algorithm>
#include <cctype>
#include <string>
#include <sstream>
#include <vector>
#include <set>

#include "circuit.h"
#include "gates_qsim.h"
#include "util.h"

namespace qsim {

/**
 * Parser for the (deprecated) qsim <a href="https://github.com/quantumlib/qsim/blob/master/docs/input_format.md">file input format</a>.
 * The primary supported interface for designing circuits to simulate with qsim
 * is <a href="https://github.com/quantumlib/Cirq">Cirq</a>, which relies on
 * the Python-based qsimcirq interface. For C++ applications, Cirq gates can be
 * explicitly constructed in code.
 */
template <typename IO>
class CircuitQsimParser final {
 public:
  /**
   * Parses the given input stream into a Circuit object, following the rules
   * defined in "docs/input_format.md".
   * @param maxtime Maximum gate "time" to read operations for (inclusive).
   * @param provider Circuit source; only used for error reporting.
   * @param fs The stream to read the circuit from.
   * @param circuit Output circuit object. If parsing is successful, this will
   *   contain the circuit defined in 'fs'.
   * @return True if parsing succeeds; false otherwise.
   */
  template<typename fp_type, typename Stream>
  static bool FromStream(unsigned maxtime, const std::string& provider,
                        Stream& fs, Circuit<GateQSim<fp_type>>& circuit) {
      circuit.num_qubits = 0;
      circuit.gates.clear();
      circuit.gates.reserve(1024);

      unsigned k = 0;
      std::string line;
      line.reserve(128);

      unsigned max_time = 0;
      unsigned prev_mea_time = 0;

      std::vector<std::stringstream> ss_list; // Temporary vectors for handling blocks
      std::vector<unsigned> time_list;
      std::vector<std::string> gate_names;
      bool in_block = false;                  // Flag to indicate if currently inside a block

      std::vector<unsigned> last_times;

      double t0 = GetTime();

      while (std::getline(fs, line)) {
          ++k;

          if (line.size() == 0 || line[0] == '#')
              continue;

          std::stringstream ss(line);

          if (circuit.num_qubits == 0) {
              ss >> circuit.num_qubits;
              if (circuit.num_qubits == 0) {
                  std::cerr << "Invalid number of qubits in " << provider << " at line " << k << std::endl;
                  return false;
              }

              last_times.resize(circuit.num_qubits, unsigned(-1));
              continue;
          }

          unsigned time;
          std::string gate_name;

          if (line[0] == '[') {
              in_block = true;
              line.erase(0, 1); // Remove the opening bracket
              std::stringstream block_ss(line);
              block_ss >> time >> gate_name;

              time_list.push_back(time);
              gate_names.push_back(gate_name);
              ss_list.push_back(std::move(block_ss));
          } else if (in_block && line.back()!= ']') { //also include lines between lines with brackets
              std::stringstream block_ss(line);
              block_ss >> time >> gate_name;
              time_list.push_back(time);
              gate_names.push_back(gate_name);
              ss_list.push_back(std::move(block_ss));
          } else if (in_block && line.back() == ']') {
              line.pop_back(); // Remove the closing bracket
              std::stringstream block_ss(line);
              block_ss >> time >> gate_name;

              time_list.push_back(time);
              gate_names.push_back(gate_name);
              ss_list.push_back(std::move(block_ss));

              double t0_i = GetTime();
              if (!ParseBlockGates<fp_type, GateQSim<fp_type>>(ss_list, time_list, gate_names, circuit.gates, circuit.num_qubits)) {
                  std::cerr << "ParseBlockGates failed at line " << k << std::endl;
                  InvalidGateError(provider, k);
                  return false;
              }
              double t1_i = GetTime();
              std::cout << "Time for ParseBlockGates: " << t1_i-t0_i << std::endl;

              time_list.clear(); // Clear current block for new block
              gate_names.clear();
              ss_list.clear();

              in_block = false; // Reset in_block flag
          } else {
              ss >> time >> gate_name;
              if (!ss) {
                  std::cerr << "Invalid gate format at line " << k << std::endl;
                  InvalidGateError(provider, k);
                  return false;
              }

              if (time > maxtime)
                  break;

              if (!in_block) {
                  if (gate_name == "c") {
                      if (!ParseGate<fp_type>(ss, time, circuit.num_qubits, gate_name, circuit.gates)) {
                          std::cerr << "ParseControlledGate failed at line " << k << std::endl;
                          InvalidGateError(provider, k);
                          return false;
                      }
                  } else {
                      if (!ParseGate<fp_type>(ss, time, circuit.num_qubits, gate_name, circuit.gates)) {
                          std::cerr << "ParseGate failed at line " << k << std::endl;
                          InvalidGateError(provider, k);
                          return false;
                      }
                  }
              }
          }

          // Additional checks and updates
          const auto& gate = circuit.gates.back();

          if (time < prev_mea_time || (gate.kind == gate::kMeasurement && time < max_time)) {
              std::cerr << "Gate crosses time boundary set by measurement gates at line " << k << " in " << provider << std::endl;
              return false;
          }

          if (gate.kind == gate::kMeasurement) {
              prev_mea_time = time;
          }

          if (GateIsOutOfOrder(time, gate.qubits, last_times)
              || GateIsOutOfOrder(time, gate.controlled_by, last_times)) {
              std::cerr << "Gate is out of time order at line " << k << " in " << provider << std::endl;
              return false;
          }

          if (time > max_time) {
              max_time = time;
          }
      }

      double t1 = GetTime();
      //IO::messagef("Time for Parsing + Block SVDs: %g seconds.\n", t1 - t0);
      std::cout << "Time for Parsing: " << t1-t0 << std::endl;

      return true;
  }

  /**
   * Parses the given file into a Circuit object, following the rules defined
   * in "docs/input_format.md".
   * @param maxtime Maximum gate "time" to read operations for (inclusive).
   * @param file The name of the file to read the circuit from.
   * @param circuit Output circuit object. If parsing is successful, this will
   *   contain the circuit defined in 'file'.
   * @return True if parsing succeeds; false otherwise.
   */
  template <typename fp_type>
  static bool FromFile(unsigned maxtime, const std::string& file,
                       Circuit<GateQSim<fp_type>>& circuit) {
    auto fs = IO::StreamFromFile(file);

    if (!fs) {
      return false;
    } else {
      bool rc = FromStream(maxtime, file, fs, circuit);
      IO::CloseStream(fs);
      return rc;
    }
  }

 private:
  static void InvalidGateError(const std::string& provider, unsigned line) {
    IO::errorf("invalid gate in %s in line %u.\n", provider.c_str(), line);
  }

  /**
   * Checks formatting for a zero-qubit gate parsed from 'ss'.
   * @param ss Input stream containing the gate specification.
   */
  static bool ValidateGate(std::stringstream& ss) {
    return ss && ss.peek() == std::stringstream::traits_type::eof();
  }

  /**
   * Checks formatting for a single-qubit gate parsed from 'ss'.
   * @param ss Input stream containing the gate specification.
   * @param num_qubits Number of qubits, as defined at the start of the file.
   * @param q0 Index of the affected qubit.
   */
  static bool ValidateGate(std::stringstream& ss,
                           unsigned num_qubits, unsigned q0) {
    return ss && ss.peek() == std::stringstream::traits_type::eof()
        && q0 < num_qubits;
  }

  /**
   * Checks formatting for a two-qubit gate parsed from 'ss'.
   * @param ss Input stream containing the gate specification.
   * @param num_qubits Number of qubits, as defined at the start of the file.
   * @param q0 Index of the first affected qubit.
   * @param q1 Index of the second affected qubit.
   */
  static bool ValidateGate(std::stringstream& ss,
                           unsigned num_qubits, unsigned q0, unsigned q1) {
    return ss && ss.peek() == std::stringstream::traits_type::eof()
        && q0 < num_qubits && q1 < num_qubits && q0 != q1;
  }

  /**
   * Checks formatting for a multiqubit gate parsed from 'ss'.
   * @param ss Input stream containing the gate specification.
   * @param num_qubits Number of qubits, as defined at the start of the file.
   * @param qubits Indices of affected qubits.
   */
  static bool ValidateGate(std::stringstream& ss, unsigned num_qubits,
                           const std::vector<unsigned>& qubits) {
    return ss && ValidateQubits(num_qubits, qubits);
  }

  static bool ValidateControlledGate(
      unsigned num_qubits, const std::vector<unsigned>& qubits,
      const std::vector<unsigned>& controlled_by) {
    if (!ValidateQubits(num_qubits, controlled_by)) return false;

    std::size_t i = 0, j = 0;

    while (i < qubits.size() && j < controlled_by.size()) {
      if (qubits[i] == controlled_by[j]) {
        return false;
      } else if (qubits[i] < controlled_by[j]) {
        ++i;
      } else {
        ++j;
      }
    }

    return true;
  }

  static bool ValidateQubits(unsigned num_qubits,
                             const std::vector<unsigned>& qubits) {
    if (qubits.size() == 0 || qubits[0] >= num_qubits) return false;

    // qubits should be sorted.

    for (std::size_t i = 1; i < qubits.size(); ++i) {
      if (qubits[i] >= num_qubits || qubits[i] == qubits[i - 1]) {
        return false;
      }
    }

    return true;
  }

  static bool GateIsOutOfOrder(unsigned time,
                               const std::vector<unsigned>& qubits,
                               std::vector<unsigned>& last_times) {
    for (auto q : qubits) {
      if (last_times[q] != unsigned(-1) && time <= last_times[q]) {
        return true;
      }

      last_times[q] = time;
    }

    return false;
  }

  static bool checkSameSize(const std::vector<std::stringstream>& ss_list, 
                              const std::vector<unsigned>& time_list, 
                              const std::vector<std::string>& gate_names) {
    return ss_list.size() == time_list.size() && time_list.size() == gate_names.size();
  }

  /**
   * Parsing method for a single block of gates
  */
  template <typename fp_type, typename Stream, typename Gate>
  static bool ParseBlockGates(std::vector<std::stringstream>& ss_list, std::vector<unsigned> time_list,
                        std::vector<std::string> gate_names, std::vector<Gate>& gates, unsigned num_qubits) {

      //check that ss_list, times_list and gate_names has the same length
      if (!checkSameSize(ss_list, time_list, gate_names)) return false;
      
      std::vector<Gate> gates_temp; 
      std::set<int> qubit_locs; //gather all qubit_locs in a set to be able to find out the size of the block
      std::vector<std::vector<unsigned int>> qubit_locations_tups; //gather the locations of the qubits to use it for block:create
      // parse each gate of the block and store the created gate into a list
      for (int i = 0; i < ss_list.size(); ++i) {

        std::stringstream& ss = ss_list[i];
        unsigned time = time_list[i];
        std::string gate_name = gate_names[i];

        int q0, q1; 
        q0=-1;
        q1=-1;//initialize -1
        fp_type phi, theta;

        //case distinction, do not inclued measurement gates, unfortunately redundand code
        if (gate_name == "p") {
          ss >> phi;
          if (!ValidateGate(ss)) return false;
          gates_temp.push_back(GateGPh<fp_type>::Create(time, phi));
        } else if (gate_name == "id1") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateId1<fp_type>::Create(time, q0));
        } else if (gate_name == "h") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateHd<fp_type>::Create(time, q0));
        } else if (gate_name == "t") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateT<fp_type>::Create(time, q0));
        } else if (gate_name == "x") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateX<fp_type>::Create(time, q0));
        } else if (gate_name == "y") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateY<fp_type>::Create(time, q0));
        } else if (gate_name == "z") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateZ<fp_type>::Create(time, q0));
        } else if (gate_name == "x_1_2") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateX2<fp_type>::Create(time, q0));
        } else if (gate_name == "y_1_2") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateY2<fp_type>::Create(time, q0));
        } else if (gate_name == "rx") {
          ss >> q0 >> phi;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateRX<fp_type>::Create(time, q0, phi));
        } else if (gate_name == "ry") {
          ss >> q0 >> phi;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateRY<fp_type>::Create(time, q0, phi));
        } else if (gate_name == "rz") {
          ss >> q0 >> phi;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateRZ<fp_type>::Create(time, q0, phi));
        } else if (gate_name == "rxy") {
          ss >> q0 >> theta >> phi;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateRXY<fp_type>::Create(time, q0, theta, phi));
        } else if (gate_name == "hz_1_2") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateHZ2<fp_type>::Create(time, q0));
        } else if (gate_name == "s") {
          ss >> q0;
          if (!ValidateGate(ss, num_qubits, q0)) return false;
          gates_temp.push_back(GateS<fp_type>::Create(time, q0));
        } else if (gate_name == "id2") {
          ss >> q0 >> q1;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          gates_temp.push_back(GateId2<fp_type>::Create(time, q0, q1));
        } else if (gate_name == "cz") {
          ss >> q0 >> q1;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          if (q1 > q0) {
            gates_temp.push_back(GateCZ<fp_type>::Create(time, q0, q1));
          } else {
            gates_temp.push_back(GateCZ<fp_type>::CreateReversed(time, q0, q1));
          }
        } else if (gate_name == "cnot" || gate_name == "cx") {
          ss >> q0 >> q1;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          if (q1 > q0) {
            gates_temp.push_back(GateCNot<fp_type>::Create(time, q0, q1));
          } else {
            gates_temp.push_back(GateCNot<fp_type>::CreateReversed(time, q0, q1));
          }
        } else if (gate_name == "sw") {
          ss >> q0 >> q1;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          if (q1 > q0) {
            gates_temp.push_back(GateSwap<fp_type>::Create(time, q0, q1));
          } else {
            gates_temp.push_back(GateSwap<fp_type>::CreateReversed(time, q0, q1));
          }
        } else if (gate_name == "is") {
          ss >> q0 >> q1;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          if (q1 > q0) {
            gates_temp.push_back(GateIS<fp_type>::Create(time, q0, q1));
          } else {
            gates_temp.push_back(GateIS<fp_type>::CreateReversed(time, q0, q1));
          }
        } else if (gate_name == "fs") {
          ss >> q0 >> q1 >> theta >> phi;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          if (q1 > q0) {
            gates_temp.push_back(GateFS<fp_type>::Create(time, q0, q1, theta, phi));
          } else {
            gates_temp.push_back(GateFS<fp_type>::CreateReversed(time, q0, q1, theta, phi));
          }
        } else if (gate_name == "cp") {
          ss >> q0 >> q1 >> phi;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          if (q1 > q0) {
            gates_temp.push_back(GateCP<fp_type>::Create(time, q0, q1, phi));
          } else {
            gates_temp.push_back(GateCP<fp_type>::CreateReversed(time, q0, q1, phi));
          }
        } else if (gate_name == "cy") {
          ss >> q0 >> q1;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          if (q1 > q0) {
            gates_temp.push_back(GateCY<fp_type>::Create(time, q0, q1));
          } else {
            gates_temp.push_back(GateCY<fp_type>::CreateReversed(time, q0, q1));
          }
        } else if (gate_name == "rzz") {
          ss >> q0 >> q1 >> theta;
          if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
          if (q1 > q0) {
            gates_temp.push_back(GateRZZ<fp_type>::Create(time, q0, q1, theta));
          } else {
            gates_temp.push_back(GateRZZ<fp_type>::CreateReversed(time, q0, q1, theta));
          }
        } else {
          return false;
        }
      qubit_locs.insert(q0);
      if (q1!=-1) { //for single qubit gates q1 is not assigned
        qubit_locs.insert(q1);
        if (q1 > q0) {
          qubit_locations_tups.push_back({static_cast<unsigned int>(q0),static_cast<unsigned int>(q1)});
        } else {
          qubit_locations_tups.push_back({static_cast<unsigned int>(q1),static_cast<unsigned int>(q0)});
        }
        
      } else {
        qubit_locations_tups.push_back({static_cast<unsigned int>(q0)});
      }

      //!todo also include the gate consturctions from parsecontrolledgate
      }
      //afterwards, this list is input for gates_qsim::GateBlocks which does the multiplication and svd and overwrites the "gates" varaible which is input
      unsigned int block_size = *std::max_element(qubit_locs.begin(), qubit_locs.end()) - *std::min_element(qubit_locs.begin(), qubit_locs.end()) + 1;
      GateBlock<fp_type> gateBlock(block_size);
      double t0_mul = GetTime();
      gates.push_back(gateBlock.Create(time_list[0], gates_temp, qubit_locations_tups));
      double t1_mul = GetTime();
      std::cout << "Time for Multiplying Unitaries (gateBlock.Create): " << t1_mul-t0_mul << std::endl;
      return true;
    }


  template <typename fp_type, typename Stream, typename Gate>
  static bool ParseGate(Stream& ss, unsigned time, unsigned num_qubits,
                        const std::string& gate_name,
                        std::vector<Gate>& gates) {
    unsigned q0, q1;
    fp_type phi, theta;

    if (gate_name == "p") {
      ss >> phi;
      if (!ValidateGate(ss)) return false;
      gates.push_back(GateGPh<fp_type>::Create(time, phi));
    } else if (gate_name == "id1") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateId1<fp_type>::Create(time, q0));
    } else if (gate_name == "h") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateHd<fp_type>::Create(time, q0));
    } else if (gate_name == "t") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateT<fp_type>::Create(time, q0));
    } else if (gate_name == "x") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateX<fp_type>::Create(time, q0));
    } else if (gate_name == "y") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateY<fp_type>::Create(time, q0));
    } else if (gate_name == "z") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateZ<fp_type>::Create(time, q0));
    } else if (gate_name == "x_1_2") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateX2<fp_type>::Create(time, q0));
    } else if (gate_name == "y_1_2") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateY2<fp_type>::Create(time, q0));
    } else if (gate_name == "rx") {
      ss >> q0 >> phi;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateRX<fp_type>::Create(time, q0, phi));
    } else if (gate_name == "ry") {
      ss >> q0 >> phi;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateRY<fp_type>::Create(time, q0, phi));
    } else if (gate_name == "rz") {
      ss >> q0 >> phi;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateRZ<fp_type>::Create(time, q0, phi));
    } else if (gate_name == "rxy") {
      ss >> q0 >> theta >> phi;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateRXY<fp_type>::Create(time, q0, theta, phi));
    } else if (gate_name == "hz_1_2") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateHZ2<fp_type>::Create(time, q0));
    } else if (gate_name == "s") {
      ss >> q0;
      if (!ValidateGate(ss, num_qubits, q0)) return false;
      gates.push_back(GateS<fp_type>::Create(time, q0));
    } else if (gate_name == "id2") {
      ss >> q0 >> q1;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateId2<fp_type>::Create(time, q0, q1));
    } else if (gate_name == "cz") {
      ss >> q0 >> q1;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateCZ<fp_type>::Create(time, q0, q1));
    } else if (gate_name == "cnot" || gate_name == "cx") {
      ss >> q0 >> q1;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateCNot<fp_type>::Create(time, q0, q1));
    } else if (gate_name == "sw") {
      ss >> q0 >> q1;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateSwap<fp_type>::Create(time, q0, q1));
    } else if (gate_name == "is") {
      ss >> q0 >> q1;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateIS<fp_type>::Create(time, q0, q1));
    } else if (gate_name == "fs") {
      ss >> q0 >> q1 >> theta >> phi;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateFS<fp_type>::Create(time, q0, q1, theta, phi));
    } else if (gate_name == "cp") {
      ss >> q0 >> q1 >> phi;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateCP<fp_type>::Create(time, q0, q1, phi));
    } else if (gate_name == "cy") {
      ss >> q0 >> q1;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateCY<fp_type>::Create(time, q0, q1));
    } else if (gate_name == "rzz") {
      ss >> q0 >> q1 >> theta;
      if (!ValidateGate(ss, num_qubits, q0, q1)) return false;
      gates.push_back(GateRZZ<fp_type>::Create(time, q0, q1, theta));
    } else if (gate_name == "m") {
      std::vector<unsigned> qubits;
      qubits.reserve(num_qubits);

      while (ss.good()) {
        ss >> q0;
        if (ss) {
          qubits.push_back(q0);
        } else {
          return false;
        }
      }

      gates.push_back(gate::Measurement<GateQSim<fp_type>>::Create(
          time, std::move(qubits)));

      if (!ValidateQubits(num_qubits, gates.back().qubits)) return false;
    } else {
      return false;
    }

    return true;
  }

  template <typename fp_type, typename Stream, typename Gate>
  static bool ParseControlledGate(Stream& ss, unsigned time,
                                  unsigned num_qubits,
                                  std::vector<Gate>& gates) {
    std::vector<unsigned> controlled_by;
    controlled_by.reserve(64);

    std::string gate_name;
    gate_name.reserve(16);

    while (1) {
      while (ss.good()) {
        if (!std::isblank(ss.get())) {
          ss.unget();
          break;
        }
      }

      if (!ss.good()) {
        return false;
      }

      if (!std::isdigit(ss.peek())) {
        break;
      } else {
        unsigned q;
        ss >> q;

        if (!ss.good() || !std::isblank(ss.get())) {
          return false;
        }

        controlled_by.push_back(q);
      }
    }

    if (controlled_by.size() == 0) {
      return false;
    }

    ss >> gate_name;

    if (!ss.good() || !ParseGate<fp_type>(ss, time,
                                          num_qubits, gate_name, gates)) {
      return false;
    }

    gates.back().ControlledBy(std::move(controlled_by));

    if (!ValidateControlledGate(num_qubits, gates.back().qubits,
                                gates.back().controlled_by)) {
      return false;
    }

    return true;
  }
};

}  // namespace qsim

#endif  // CIRCUIT_QSIM_PARSER_H_
