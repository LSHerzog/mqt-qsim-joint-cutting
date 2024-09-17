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

#ifndef HYBRID_H_
#define HYBRID_H_

#include <algorithm>
#include <array>
#include <complex>
#include <vector>

#include "gate.h"
#include "gate_appl.h"

namespace qsim {

/**
 * Hybrid Feynman-Schrodinger simulator.
 */
template <typename IO, typename GateT,
          template <typename, typename> class FuserT, typename For>
struct HybridSimulator final {
 public:
  using Gate = GateT;
  using GateKind = typename Gate::GateKind;
  using fp_type = typename Gate::fp_type;

 private:
  // Note that one can use "struct GateHybrid : public Gate {" in C++17.
  struct GateHybrid {
    using GateKind = HybridSimulator::GateKind;
    using fp_type = HybridSimulator::fp_type;

    GateKind kind;
    unsigned time;
    std::vector<unsigned> qubits;
    std::vector<unsigned> controlled_by;
    uint64_t cmask;
    std::vector<fp_type> params;
    Matrix<fp_type> matrix;
    bool unfusible;
    bool swapped;

    const Gate* parent;
    unsigned id;
  };

  struct GateX {
    GateHybrid* decomposed0;
    GateHybrid* decomposed1;
    schmidt_decomp_type<fp_type> schmidt_decomp;
    unsigned schmidt_bits;
    unsigned swapped;
  };

 public:
  using Fuser = FuserT<IO, GateHybrid>;
  using GateFused = typename Fuser::GateFused;

  /**
   * Contextual data for hybrid simulation.
   */
  struct HybridData {
    /**
     * List of gates on the "0" side of the cut.
     */
    std::vector<GateHybrid> gates0;
    /**
     * List of gates on the "1" side of the cut.
     */
    std::vector<GateHybrid> gates1;
    /**
     * List of gates on the cut.
     */
    std::vector<GateX> gatexs;
    /**
     * Global qubit index to local qubit index map.
     */
    std::vector<unsigned> qubit_map;
    /**
     * Number of qubits on the "0" side of the cut.
     */
    unsigned num_qubits0;
    /**
     * Number of qubits on the "1" side of the cut.
     */
    unsigned num_qubits1;
    /**
     * Number of gates on the cut.
     */
    unsigned num_gatexs;
  };

  /**
   * User-specified parameters for gate fusion and hybrid simulation.
   */
  struct Parameter : public Fuser::Parameter {
    /**
     * Fixed bitstring indicating values to assign to Schmidt decomposition
     * indices of prefix gates.
     */
    uint64_t prefix;
    /**
     * Number of gates on the cut that are part of the prefix. Indices of these
     * gates are assigned the value indicated by `prefix`.
     */
    unsigned num_prefix_gatexs;
    /**
     * Number of gates on the cut that are part of the root. All gates that are
     * not part of the prefix or root are part of the suffix.
     */
    unsigned num_root_gatexs;
    unsigned num_threads;
  };

  template <typename... Args>
  explicit HybridSimulator(Args&&... args) : for_(args...) {}

  /**
   * Splits the lattice into two parts, using Schmidt decomposition for gates
   * on the cut.
   * @param parts Lattice sections to be simulated.
   * @param gates List of all gates in the circuit.
   * @param hd Output data with split parts.
   * @return True if the splitting done successfully; false otherwise.
   */
  static bool SplitLattice(const std::vector<unsigned>& parts,
                           const std::vector<Gate>& gates, HybridData& hd) {
    hd.num_gatexs = 0;
    hd.num_qubits0 = 0;
    hd.num_qubits1 = 0;

    hd.gates0.reserve(gates.size());
    hd.gates1.reserve(gates.size());
    hd.qubit_map.reserve(parts.size());

    unsigned count0 = 0;
    unsigned count1 = 0;

    // Global qubit index to local qubit index map.
    for (std::size_t i = 0; i < parts.size(); ++i) {
      parts[i] == 0 ? ++hd.num_qubits0 : ++hd.num_qubits1;
      hd.qubit_map.push_back(parts[i] == 0 ? count0++ : count1++);
    }

    std::vector<unsigned int> n_top_lst; //only for blocks, list because different value for each block possible
    std::vector<unsigned int> n_bottom_lst;

    // Split the lattice.
    for (const auto& gate : gates) {
      if (gate.kind == gate::kMeasurement) {
        IO::errorf("measurement gates are not suported by qsimh.\n");
        return false;
      }

      if (gate.controlled_by.size() > 0) {
        IO::errorf("controlled gates are not suported by qsimh.\n");
        return false;
      }

      std::vector<int> qubitGroups;
      if (gate.qubits.size()>2) {
        std::vector<int> group0_idx, group1_idx;
        //find qubit indices per group
        for (unsigned int i = 0; i < parts.size(); ++i) {
          if (parts[i] == 0) {
            group0_idx.push_back(i);
          }
          else if (parts[i] == 1) {
            group1_idx.push_back(i);
          }
          else {
            IO::errorf("More than two groups for the cut are not yet implemented");
          }  
        }
        for (int q : gate.qubits) {
            if (find(group0_idx.begin(), group0_idx.end(), q) != group0_idx.end()) {
                qubitGroups.push_back(0);
            } else if (find(group1_idx.begin(), group1_idx.end(), q) != group1_idx.end()) {
                qubitGroups.push_back(1);
            } else {
                IO::errorf("Qubits does not correspond to any group defined by parts.");
            }
        }
      }
      switch (gate.qubits.size()) {
      case 1:  // Single qubit gates.
        switch (parts[gate.qubits[0]]) {
        case 0:
          hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time,
            {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, gate.matrix,
            false, false, nullptr, 0});
          break;
        case 1:
          hd.gates1.emplace_back(GateHybrid{gate.kind, gate.time,
            {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, gate.matrix,
            false, false, nullptr, 0});
          break;
        }
        break;
      case 2:  // Two qubit gates.
        {
          if (gate.kind == 24) { //only add values if block gate (24 for qsim gates, 47 for cirq gates)
            n_top_lst.emplace_back(1); //no other choice for 2 qubit blocks
            n_bottom_lst.emplace_back(1);
          };
          switch ((parts[gate.qubits[1]] << 1) | parts[gate.qubits[0]]) {
          case 0:  // Both qubits in part 0.
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time,
              {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]]},
              {}, 0, gate.params, gate.matrix, false, gate.swapped,
              nullptr, 0});
            break;
          case 1:  // Gate on the cut, qubit 0 in part 1, qubit 1 in part 0.
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time,
              {hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {},
              true, gate.swapped, &gate, hd.num_gatexs}); //! the gate is reference two times, because later they go through all hd.gates0 and hd.gates1 and catch the schmidt
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time,
              {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {},
              true, gate.swapped, &gate, hd.num_gatexs});;

            ++hd.num_gatexs;
            break;
          case 2:  // Gate on the cut, qubit 0 in part 0, qubit 1 in part 1.
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time,
              {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {},
              true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time,
              {hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {},
              true, gate.swapped, &gate, hd.num_gatexs});

            ++hd.num_gatexs;
            break;
          case 3:  // Both qubits in part 1.
            hd.gates1.emplace_back(GateHybrid{gate.kind, gate.time,
              {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]]},
              {}, 0, gate.params, gate.matrix, false, gate.swapped,
              nullptr, 0});
            break;
          default:
            IO::errorf("For the Two-Qubit Case, none of the partitions applies --> check your implementation");
          }
        }
        break;
      case 3: //3 qubit gates
        {
          for (size_t i = 0; i < gate.qubits.size()-1; ++i) {
            if (gate.qubits[i]>gate.qubits[i+1]){
              IO::errorf("A 3-qubit block must be already in proper ordering when passed to SplitLattice");
            }
          }
          
          if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0){
            //all qubits in group 0
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1) {
            //all qubits in group 1
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 1) {
            //qubits[0] in group 0 and qubits[1,2] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 0) {
            //qubits[0] in group 1 and qubits[1,2] in group 0
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 1) {
            //qubits[0,1] in group 0 and qubits[2] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 0) {
            //qubits[0,1] in group 1 and qubits[2] in group 0
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 0) {
            //qubits[0,2] in group 0 and qubits[1] in group 1
            IO::errorf("We assume that a block is cut only once, not twice");
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 1) {
            //qubits[0,2] in group 1 and qubits[1] in group 0
            IO::errorf("We assume that a block is cut only once, not twice");
          }
        }
        break;
      case 4: //4 qubit gates
        {
          for (size_t i = 0; i < gate.qubits.size()-1; ++i) {
            if (gate.qubits[i]>gate.qubits[i+1]){
              IO::errorf("A 4-qubit block must be already in proper ordering when passed to SplitLattice");
            }
          }

          if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0){
            //all qubits in group 0
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1) {
            //all qubits in group 1
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 1) {
            //qubits[0,1,2] in group 0 and qubits[3] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 1 && qubitGroups[3] == 1) {
            //qubits[0,1] in group 0 and qubits[2,3] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1) {
            //qubits[0] in group 0 and qubits[1,2,3] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0) {
            //qubits[0] in group 1 and qubits[1,2,3] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 0 && qubitGroups[3] == 0) {
            //qubits[0,1] in group 1 and qubits[2,3] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 0) {
            //qubits[0,1,2] in group 1 and qubits[3] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(1);
            };
          } else {
            IO::errorf("We assume that a block is cut only once, not twice or three times.");
          }
        }
        break;
      case 5: //5 qubit gates
        {
          for (size_t i = 0; i < gate.qubits.size()-1; ++i) {
            if (gate.qubits[i]>gate.qubits[i+1]){
              IO::errorf("A 5-qubit block must be already in proper ordering when passed to SplitLattice");
            }
          }

          if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0){
            //all qubits in group 0
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1) {
            //all qubits in group 1
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 1) {
            //qubits[0,1,2,3] in group 0 and qubits[4] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 1 && qubitGroups[4] == 1) {
            //qubits[0,1,2] in group 0 and qubits[3,4] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1) {
            //qubits[0,1] in group 0 and qubits[2,3,4] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1) {
            //qubits[0] in group 0 and qubits[1,2,3,4] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 0) {
            //qubits[0,1,2,3] in group 1 and qubits[4] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 0 && qubitGroups[4] == 0) {
            //qubits[0,1,2] in group 1 and qubits[3,4] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0) {
            //qubits[0,1] in group 1 and qubits[2,3,4] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0) {
            //qubits[0] in group 1 and qubits[1,2,3,4] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(4);
            };
          } else {
            IO::errorf("We assume that a block is cut only once, not twice, three or four times.");
          }
        }
        break;
      case 6: //6 qubit gate
        {
          for (size_t i = 0; i < gate.qubits.size()-1; ++i) {
            if (gate.qubits[i]>gate.qubits[i+1]){
              IO::errorf("A 6-qubit block must be already in proper ordering when passed to SplitLattice");
            }
          }

          if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0){
            //all qubits in group 0
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1) {
            //all qubits in group 1
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 1) {
            //qubits[0,1,2,3,4] in group 0 and qubits[5] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 1 && qubitGroups[5] == 1) {
            //qubits[0,1,2,3] in group 0 and qubits[4,5] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1) {
            //qubits[0,1,2] in group 0 and qubits[3,4,5] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1) {
            //qubits[0,1] in group 0 and qubits[2,3,4,5] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1) {
            //qubits[0] in group 0 and qubits[1,2,3,4,5] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0) {
            //qubits[0] in group 1 and qubits[1,2,3,4,5] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0) {
            //qubits[0,1] in group 1 and qubits[2,3,4,5] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0) {
            //qubits[0,1,2] in group 1 and qubits[3,4,5] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 0 && qubitGroups[5] == 0) {
            //qubits[0,1,2,3] in group 1 and qubits[4,5] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 0) {
            //qubits[0,1,2,3,4] in group 1 and qubits[5] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(1);
            };
          }
        }
      break;
      case 7: // 7qubit gates
        {
          for (size_t i = 0; i < gate.qubits.size()-1; ++i) {
            if (gate.qubits[i]>gate.qubits[i+1]){
              IO::errorf("A 7-qubit block must be already in proper ordering when passed to SplitLattice");
            }
          }

          if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0  && qubitGroups[6] == 0){
            //all qubits in group 0
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1) {
            //all qubits in group 1
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 1) {
            //qubits[0,1,2,3,4,5] in group 0 and qubits[6] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(6);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 1 && qubitGroups[6] == 1) {
            //qubits[0,1,2,3,4] in group 0 and qubits[5, 6] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1) {
            //qubits[0,1,2,3] in group 0 and qubits[4,5, 6] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1) {
            //qubits[0,1,2] in group 0 and qubits[3,4,5, 6] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1) {
            //qubits[0,1] in group 0 and qubits[2,3,4,5, 6] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1) {
            //qubits[0] in group 0 and qubits[1,2,3,4,5, 6] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(6);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0) {
            //qubits[0] in group 1 and qubits[1,2,3,4,5,6] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(6);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0) {
            //qubits[0,1] in group 1 and qubits[2,3,4,5,6] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0) {
            //qubits[0,1,2] in group 1 and qubits[3,4,5,6] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0) {
            //qubits[0,1,2,3] in group 1 and qubits[4,5,6] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 0 && qubitGroups[6] == 0) {
            //qubits[0,1,2,3,4] in group 1 and qubits[5,6] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 0) {
            //qubits[0,1,2,3,4,5] in group 1 and qubits[6] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(6);
              n_bottom_lst.emplace_back(1);
            };
          }
        }
      break;
      case 8: // 8 qubit block
        {
          for (size_t i = 0; i < gate.qubits.size()-1; ++i) {
            if (gate.qubits[i]>gate.qubits[i+1]){
              IO::errorf("A 8-qubit block must be already in proper ordering when passed to SplitLattice");
            }
          }

          if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0  && qubitGroups[6] == 0 && qubitGroups[7] == 0){
            //all qubits in group 0
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1, qubitGroups[6] == 1 && qubitGroups[7] == 1) {
            //all qubits in group 1
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 1) {
            //qubits[0,1,2,3,4,5,6] in group 0 and qubits[7] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(7);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 1 && qubitGroups[7] == 1) {
            //qubits[0,1,2,3,4,5] in group 0 and qubits[6,7] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(6);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1) {
            //qubits[0,1,2,3,4] in group 0 and qubits[5,6,7] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1) {
            //qubits[0,1,2,3] in group 0 and qubits[4,5,6,7] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1) {
            //qubits[0,1,2] in group 0 and qubits[3,4,5,6,7] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1) {
            //qubits[0,1] in group 0 and qubits[2,3,4,5,6,7] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(6);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1) {
            //qubits[0] in group 0 and qubits[1,2,3,4,5,6,7] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(7);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0) {
            //qubits[0] in group 1 and qubits[1,2,3,4,5,6,7] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(7);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0) {
            //qubits[0,1] in group 1 and qubits[2,3,4,5,6,7] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(6);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0) {
            //qubits[0,1,2] in group 1 and qubits[3,4,5,6,7] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0) {
            //qubits[0,1,2,3] in group 1 and qubits[4,5,6,7] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0) {
            //qubits[0,1,2,3,4] in group 1 and qubits[5,6,7] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 0 && qubitGroups[7] == 0) {
            //qubits[0,1,2,3,4,5] in group 1 and qubits[6,7] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(6);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 0) {
            //qubits[0,1,2,3,4,5,6] in group 1 and qubits[7] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(7);
              n_bottom_lst.emplace_back(1);
            };
          } 
        }
      break;
      case 9 : //9 qubit blocks
        {
          for (size_t i = 0; i < gate.qubits.size()-1; ++i) {
            if (gate.qubits[i]>gate.qubits[i+1]){
              IO::errorf("A 9-qubit block must be already in proper ordering when passed to SplitLattice");
            }
          }

          if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0  && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0){
            //all qubits in group 0
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]], hd.qubit_map[gate.qubits[8]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1, qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1) {
            //all qubits in group 1
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]], hd.qubit_map[gate.qubits[8]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 1) {
            //qubits[0,1,2,3,4,5,6,7] in group 0 and qubits[8] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(8);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 1 && qubitGroups[8] == 1) {
            //qubits[0,1,2,3,4,5,6] in group 0 and qubits[7,8] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(7);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1) {
            //qubits[0,1,2,3,4,5] in group 0 and qubits[6,7,8] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(6);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1) {
            //qubits[0,1,2,3,4] in group 0 and qubits[5,6,7,8] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1) {
            //qubits[0,1,2,3] in group 0 and qubits[4,5,6,7,8] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1) {
            //qubits[0,1,2] in group 0 and qubits[3,4,5,6,7,8] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(6);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1) {
            //qubits[0,1] in group 0 and qubits[2,3,4,5,6,7,8] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(7);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1) {
            //qubits[0] in group 0 and qubits[1,2,3,4,5,6,7,8] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(8);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0) {
            //qubits[0] in group 1 and qubits[1,2,3,4,5,6,7,8] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(8);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0) {
            //qubits[0,1] in group 1 and qubits[2,3,4,5,6,7,8] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(7);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0) {
            //qubits[0,1,2] in group 1 and qubits[3,4,5,6,7,8] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(6);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0) {
            //qubits[0,1,2,3] in group 1 and qubits[4,5,6,7,8] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0) {
            //qubits[0,1,2,3,4] in group 1 and qubits[5,6,7,8] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0) {
            //qubits[0,1,2,3,4,5] in group 1 and qubits[6,7,8] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(6);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 0 && qubitGroups[8] == 0) {
            //qubits[0,1,2,3,4,5,6] in group 1 and qubits[7,8] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(7);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 0) {
            //qubits[0,1,2,3,4,5,6,7] in group 1 and qubits[8] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(8);
              n_bottom_lst.emplace_back(1);
            };
          }
        }
      break;
      case 10 : //10 qubit block
        {
          for (size_t i = 0; i < gate.qubits.size()-1; ++i) {
            if (gate.qubits[i]>gate.qubits[i+1]){
              IO::errorf("A 10-qubit block must be already in proper ordering when passed to SplitLattice");
            }
          }

          if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0  && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0 && qubitGroups[9] == 0){
            //all qubits in group 0
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]], hd.qubit_map[gate.qubits[8]], hd.qubit_map[gate.qubits[9]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1, qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1 && qubitGroups[9] == 1) {
            //all qubits in group 1
            hd.gates0.emplace_back(GateHybrid{gate.kind, gate.time, {hd.qubit_map[gate.qubits[0]], hd.qubit_map[gate.qubits[1]], hd.qubit_map[gate.qubits[2]], hd.qubit_map[gate.qubits[3]], hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]], hd.qubit_map[gate.qubits[8]], hd.qubit_map[gate.qubits[9]]}, {} ,0, gate.params, gate.matrix, false, gate.swapped, nullptr, 0});
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0  && qubitGroups[9] == 1) {
            //qubits[0,1,2,3,4,5,6,7,8] in group 0 and qubits[9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]], hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(9);
              n_bottom_lst.emplace_back(1);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 1  && qubitGroups[9] == 1) {
            //qubits[0,1,2,3,4,5,6,7] in group 0 and qubits[8,9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]], hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(8);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 1 && qubitGroups[8] == 1  && qubitGroups[9] == 1) {
            //qubits[0,1,2,3,4,5,6] in group 0 and qubits[7,8,9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]], hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(7);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1  && qubitGroups[9] == 1) {
            //qubits[0,1,2,3,4,5] in group 0 and qubits[6,7,8,9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]], hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(6);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1  && qubitGroups[9] == 1) {
            //qubits[0,1,2,3,4] in group 0 and qubits[5,6,7,8,9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1  && qubitGroups[9] == 1) {
            //qubits[0,1,2,3] in group 0 and qubits[4,5,6,7,8,9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(6);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1  && qubitGroups[9] == 1) {
            //qubits[0,1,2] in group 0 and qubits[3,4,5,6,7,8,9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(7);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 0 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1  && qubitGroups[9] == 1) {
            //qubits[0,1] in group 0 and qubits[2,3,4,5,6,7,8,9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(8);
            };
          } else if (qubitGroups[0] == 0 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1  && qubitGroups[9] == 1) {
            //qubits[0] in group 0 and qubits[1,2,3,4,5,6,7,8,9] in group 1
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(9);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 0 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0  && qubitGroups[9] == 0) {
            //qubits[0] in group 1 and qubits[1,2,3,4,5,6,7,8,9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(1);
              n_bottom_lst.emplace_back(9);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 0 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0  && qubitGroups[9] == 0) {
            //qubits[0,1] in group 1 and qubits[2,3,4,5,6,7,8,9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(2);
              n_bottom_lst.emplace_back(8);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 0 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0  && qubitGroups[9] == 0) {
            //qubits[0,1,2] in group 1 and qubits[3,4,5,6,7,8,9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(3);
              n_bottom_lst.emplace_back(7);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 0 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0  && qubitGroups[9] == 0) {
            //qubits[0,1,2,3] in group 1 and qubits[4,5,6,7,8,9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(4);
              n_bottom_lst.emplace_back(6);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 0 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0  && qubitGroups[9] == 0) {
            //qubits[0,1,2,3,4] in group 1 and qubits[5,6,7,8,9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(5);
              n_bottom_lst.emplace_back(5);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 0 && qubitGroups[7] == 0 && qubitGroups[8] == 0  && qubitGroups[9] == 0) {
            //qubits[0,1,2,3,4,5] in group 1 and qubits[6,7,8,9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(6);
              n_bottom_lst.emplace_back(4);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 0 && qubitGroups[8] == 0  && qubitGroups[9] == 0) {
            //qubits[0,1,2,3,4,5,6] in group 1 and qubits[7,8,9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(7);
              n_bottom_lst.emplace_back(3);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 0  && qubitGroups[9] == 0) {
            //qubits[0,1,2,3,4,5,6,7] in group 1 and qubits[8,9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[8]],hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(8);
              n_bottom_lst.emplace_back(2);
            };
          } else if (qubitGroups[0] == 1 && qubitGroups[1] == 1 && qubitGroups[2] == 1 && qubitGroups[3] == 1 && qubitGroups[4] == 1 && qubitGroups[5] == 1 && qubitGroups[6] == 1 && qubitGroups[7] == 1 && qubitGroups[8] == 1  && qubitGroups[9] == 0) {
            //qubits[0,1,2,3,4,5,6,7,8] in group 1 and qubits[9] in group 0
            hd.gates0.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[9]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            hd.gates1.emplace_back(GateHybrid{GateKind::kDecomp, gate.time, {hd.qubit_map[gate.qubits[0]],hd.qubit_map[gate.qubits[1]],hd.qubit_map[gate.qubits[2]],hd.qubit_map[gate.qubits[3]],hd.qubit_map[gate.qubits[4]],hd.qubit_map[gate.qubits[5]],hd.qubit_map[gate.qubits[6]],hd.qubit_map[gate.qubits[7]],hd.qubit_map[gate.qubits[8]]}, {}, 0, gate.params, {}, true, gate.swapped, &gate, hd.num_gatexs});
            ++hd.num_gatexs;
            if (gate.kind == 24 ) {
              n_top_lst.emplace_back(9);
              n_bottom_lst.emplace_back(1);
            };
          }
        }
      break;
      default:
        IO::errorf("multi-qubit gates (>10) are not suported by qsimh.\n");
        return false;
      }
    }

    auto compare = [](const GateHybrid& l, const GateHybrid& r) -> bool {
      return l.time < r.time || (l.time == r.time &&
          (l.parent < r.parent || (l.parent == r.parent && l.id < r.id)));
    };

    // Sort gates.
    std::sort(hd.gates0.begin(), hd.gates0.end(), compare);
    std::sort(hd.gates1.begin(), hd.gates1.end(), compare);

    hd.gatexs.reserve(hd.num_gatexs); 

    // Get Schmidt matrices.
    unsigned int counter_ntop_nbottom = 0;
    for (auto& gate0 : hd.gates0) {
      if (gate0.parent != nullptr) { //condition met only for cut gates
        unsigned int n_top_temp;
        unsigned int n_bottom_temp;
        if (n_top_lst.size()!=0){
          n_top_temp = n_top_lst[counter_ntop_nbottom];
          n_bottom_temp = n_bottom_lst[counter_ntop_nbottom];
          if (counter_ntop_nbottom >= n_top_lst.size()+1) {
            throw std::runtime_error("Unexpected Behaviour will arise if your n_top_lst / n_bottom_lst have fewer elements than counter_ntop_nbottom");
          } 
        }
        double t0_sch = GetTime();
        auto d = GetSchmidtDecomp(gate0.parent->kind, gate0.parent->params, gate0.parent->matrix, gate0.parent->qubits, gate0.parent->flattened_qubits, n_top_temp, n_bottom_temp);
        double t1_sch = GetTime();

        if (gate0.parent->kind == kBlock){ //counter should increase only for kBlock, otherwise mixed circuits with kBlock and other cut gates not possible
          counter_ntop_nbottom += 1;
        }
        
        if (d.size() == 0) {
          IO::errorf("no Schmidt decomposition for gate kind %u.\n",
                     gate0.parent->kind);
          return false;
        }

        //adapt here for multiqubit gates
        unsigned schmidt_bits = SchmidtBits(d.size());
        if (schmidt_bits > 10) {
          std::cout << d.size() << d.size() << std::endl;
          std::cout << "schmidt rank "<< schmidt_bits<< std::endl;
          IO::errorf("Schmidt rank is too large for gate kind %u.\n",
                     gate0.parent->kind);
          return false;
        }
        
        unsigned swapped = parts[gate0.parent->qubits[0]];
        if (gate0.parent->swapped) swapped = 1 - swapped;
        hd.gatexs.emplace_back(GateX{&gate0, nullptr, std::move(d),
                                     schmidt_bits, swapped});
      }
    }

    unsigned count = 0;
    for (auto& gate1 : hd.gates1) {
      if (gate1.parent != nullptr) {
        hd.gatexs[count++].decomposed1 = &gate1;
      }
    }

    for (auto& gatex : hd.gatexs) {
      if (gatex.schmidt_decomp.size() == 1) {
        FillSchmidtMatrices(0, gatex);
      }
    }

    return true;
  }

  /**
   * Runs the hybrid simulator on a sectioned lattice.
   * @param param Options for parallelism and logging. Also specifies the size
   *   of the 'prefix' and 'root' sections of the lattice.
   * @param factory Object to create simulators and state spaces.
   * @param hd Container object for gates on the boundary between lattice
   *   sections.
   * @param parts Lattice sections to be simulated.
   * @param fgates0 List of gates from one section of the lattice.
   * @param fgates1 List of gates from the other section of the lattice.
   * @param bitstrings List of output states to simulate, as bitstrings.
   * @param results Output vector of amplitudes. After a successful run, this
   *   will be populated with amplitudes for each state in 'bitstrings'.
   * @return True if the simulation completed successfully; false otherwise.
   */
  template <typename Factory, typename Results>
  bool Run(const Parameter& param, const Factory& factory,
           HybridData& hd, const std::vector<unsigned>& parts,
           const std::vector<GateFused>& fgates0,
           const std::vector<GateFused>& fgates1,
           const std::vector<uint64_t>& bitstrings, Results& results) const {
    using Simulator = typename Factory::Simulator;
    using StateSpace = typename Simulator::StateSpace;
    using State = typename StateSpace::State;

    unsigned num_p_gates = param.num_prefix_gatexs;
    unsigned num_pr_gates = num_p_gates + param.num_root_gatexs;

    auto bits = CountSchmidtBits(param, hd.gatexs);

    uint64_t rmax = uint64_t{1} << bits.num_r_bits;
    uint64_t smax = uint64_t{1} << bits.num_s_bits;

    std::cout << "rmax: " << rmax << std::endl;
    std::cout << "smax: " << smax << std::endl;

    auto loc0 = CheckpointLocations(param, fgates0);
    auto loc1 = CheckpointLocations(param, fgates1);

    struct Index {
      unsigned i0;
      unsigned i1;
    };

    std::vector<Index> indices;
    indices.reserve(bitstrings.size());

    //count the number of paths
    unsigned int num_paths = 1; //neutral element of mult
    for (auto& gate_obj : hd.gatexs) {
        num_paths = num_paths * gate_obj.schmidt_decomp.size();
    }
    std::cout << "Number of Feynman Paths: " << num_paths << std::endl;

    // Bitstring indices for part 0 and part 1. TODO: optimize.
    for (const auto& bitstring : bitstrings) {
      Index index{0, 0};

      for (uint64_t i = 0; i < hd.qubit_map.size(); ++i) {
        unsigned m = ((bitstring >> i) & 1) << hd.qubit_map[i];
        parts[i] ? index.i1 |= m : index.i0 |= m;
      }

      indices.push_back(index);
    }

    StateSpace state_space = factory.CreateStateSpace();

    State* rstate0;
    State* rstate1;

    State state0p = state_space.Null();
    State state1p = state_space.Null();
    State state0r = state_space.Null();
    State state1r = state_space.Null();
    State state0s = state_space.Null();
    State state1s = state_space.Null();

    // Create states.

    if (!CreateStates(hd.num_qubits0, hd.num_qubits1, state_space, true,
                      state0p, state1p, rstate0, rstate1)) {
      return false;
    }

    if (!CreateStates(hd.num_qubits0, hd.num_qubits1, state_space, rmax > 1,
                      state0r, state1r, rstate0, rstate1)) {
      return false;
    }

    if (!CreateStates(hd.num_qubits0, hd.num_qubits1, state_space, smax > 1,
                      state0s, state1s, rstate0, rstate1)) {
      return false;
    }

    state_space.SetStateZero(state0p);
    state_space.SetStateZero(state1p);

    Simulator simulator = factory.CreateSimulator();

    std::vector<unsigned> prev(hd.num_gatexs, unsigned(-1));

    double t0_prefix = GetTime();
    // param.prefix encodes the prefix path.
    unsigned gatex_index = SetSchmidtMatrices(
        0, num_p_gates, param.prefix, prev, hd.gatexs);

    if (gatex_index == 0) {
      // Apply gates before the first checkpoint.
      ApplyGates(fgates0, 0, loc0[0], simulator, state0p);
      ApplyGates(fgates1, 0, loc1[0], simulator, state1p);
    } else {
      IO::errorf("invalid prefix %lu for prefix gate index %u.\n",
                 param.prefix, gatex_index - 1);
      return false;
    }
    double t1_prefix = GetTime();

    double t0_root = GetTime();
    // Branch over root gates on the cut. r encodes the root path.
    for (uint64_t r = 0; r < rmax; ++r) {
      if (rmax > 1) {
        state_space.Copy(state0p, state0r);
        state_space.Copy(state1p, state1r);
      }

      if (SetSchmidtMatrices(num_p_gates, num_pr_gates,
                             r, prev, hd.gatexs) == 0) {
        // Apply gates before the second checkpoint.
        ApplyGates(fgates0, loc0[0], loc0[1], simulator, state0r);
        ApplyGates(fgates1, loc1[0], loc1[1], simulator, state1r);
      } else {
        continue;
      }
      double t1_root = GetTime();

      double t0_suffix = GetTime();
      // Branch over suffix gates on the cut. s encodes the suffix path.
      for (uint64_t s = 0; s < smax; ++s) {
        if (smax > 1) {
          state_space.Copy(rmax > 1 ? state0r : state0p, state0s);
          state_space.Copy(rmax > 1 ? state1r : state1p, state1s);
        }

        double t0_temp = GetTime();
        if (SetSchmidtMatrices(num_pr_gates, hd.num_gatexs,
                               s, prev, hd.gatexs) == 0) {
          // Apply the rest of the gates.
          ApplyGates(fgates0, loc0[1], fgates0.size(), simulator, state0s);
          ApplyGates(fgates1, loc1[1], fgates1.size(), simulator, state1s);
        } else {
          continue;
        }

        t0_temp = GetTime();
        auto f = [](unsigned n, unsigned m, uint64_t i,
                    const StateSpace& state_space,
                    const State& state0, const State& state1,
                    const std::vector<Index>& indices, Results& results) {
          // TODO: make it faster for the CUDA state space.
          auto a0 = state_space.GetAmpl(state0, indices[i].i0);
          auto a1 = state_space.GetAmpl(state1, indices[i].i1);
          results[i] += a0 * a1;
        };

        // Collect results.
        t0_temp = GetTime();
        for_.Run(results.size(), f,
                 state_space, *rstate0, *rstate1, indices, results);
      }
      double t1_suffix = GetTime();
    }

    return true;
  }

 private:
  /**
   * Identifies when to save "checkpoints" of the simulation state. These allow
   * runs with different cut-index values to reuse parts of the simulation.
   * @param param Options for parallelism and logging. Also specifies the size
   *   of the 'prefix' and 'root' sections of the lattice.
   * @param fgates Set of gates for which to find checkpoint locations.
   * @return A pair of numbers specifying how many gates to apply before the
   *   first and second checkpoints, respectively.
   */
  static std::array<unsigned, 2> CheckpointLocations(
      const Parameter& param, const std::vector<GateFused>& fgates) {
    std::array<unsigned, 2> loc{0, 0};

    unsigned num_decomposed = 0;
    unsigned num_p_gates = param.num_prefix_gatexs;
    unsigned num_pr_gates = num_p_gates + param.num_root_gatexs;

    for (std::size_t i = 0; i < fgates.size(); ++i) {
      for (auto gate: fgates[i].gates) {
        if (gate->parent != nullptr) {
          ++num_decomposed;
          // There should be only one decomposed gate in fused gate.
          break;
        }
      }

      if (num_decomposed <= num_p_gates) {
        loc[0] = i + 1;
      }

      if (num_decomposed <= num_pr_gates) {
        loc[1] = i + 1;
      }
    }

    return loc;
  }

  struct Bits {
    unsigned num_p_bits;
    unsigned num_r_bits;
    unsigned num_s_bits;
  };

  static Bits CountSchmidtBits(
      const Parameter& param, const std::vector<GateX>& gatexs) {
    Bits bits{0, 0, 0};

    unsigned num_p_gates = param.num_prefix_gatexs;
    unsigned num_pr_gates = num_p_gates + param.num_root_gatexs;

    for (std::size_t i = 0; i < gatexs.size(); ++i) {
      const auto& gatex = gatexs[i];
      if (i < num_p_gates) {
        bits.num_p_bits += gatex.schmidt_bits;
      } else if (i < num_pr_gates) {
        bits.num_r_bits += gatex.schmidt_bits;
      } else {
        bits.num_s_bits += gatex.schmidt_bits;
      }
    }

    return bits;
  }

  static unsigned SetSchmidtMatrices(std::size_t i0, std::size_t i1,
                                     uint64_t path,
                                     std::vector<unsigned>& prev_k,
                                     std::vector<GateX>& gatexs) {
    unsigned shift_length = 0;

    for (std::size_t i = i0; i < i1; ++i) {
      const auto& gatex = gatexs[i];

      if (gatex.schmidt_bits == 0) {
        // Continue if gatex has Schmidt rank 1.
        continue;
      }

      unsigned k = (path >> shift_length) & ((1 << gatex.schmidt_bits) - 1);
      shift_length += gatex.schmidt_bits;

      if (k != prev_k[i]) {
        if (k >= gatex.schmidt_decomp.size()) {
          // Invalid path. Returns gatex index plus one to report error in case
          // of invalid prefix.
          return i + 1;
        }

        FillSchmidtMatrices(k, gatex);

        prev_k[i] = k;
      }
    }

    return 0;
  }

  static void FillSchmidtMatrices(unsigned k, const GateX& gatex) {
    unsigned part0 = gatex.swapped;
    unsigned part1 = 1 - part0;
    {
      gatex.decomposed0->matrix.resize(gatex.schmidt_decomp[k][part0].size());
      auto begin = gatex.schmidt_decomp[k][part0].begin();
      auto end = gatex.schmidt_decomp[k][part0].end();
      std::copy(begin, end, gatex.decomposed0->matrix.begin());
    }
    {
      gatex.decomposed1->matrix.resize(gatex.schmidt_decomp[k][part1].size());
      auto begin = gatex.schmidt_decomp[k][part1].begin();
      auto end = gatex.schmidt_decomp[k][part1].end();
      std::copy(begin, end, gatex.decomposed1->matrix.begin());
    }
  }

  template <typename Simulator>
  static void ApplyGates(const std::vector<GateFused>& gates,
                         std::size_t i0, std::size_t i1,
                         const Simulator& simulator,
                         typename Simulator::State& state) {
    for (std::size_t i = i0; i < i1; ++i) {
      if (gates[i].matrix.size() > 0) {
        ApplyFusedGate(simulator, gates[i], state);
      } else {
        auto gate = gates[i];
        CalculateFusedMatrix(gate);
        ApplyFusedGate(simulator, gate, state);
      }
    }
  }

  static unsigned SchmidtBits(unsigned size) {
    if (16 < size && size <= 32 ) {
      return 5;
    } else if (32 < size && size <= 64) {
      return 6;
    } else if (64 < size && size <= 128) {
      return 7;
    } else if (128 < size && size <= 256) {
      return 8;
    } else if (256 < size && size <= 512) {
      return 9;
    } else if (512 < size && size <= 1024) {
      return 10;
    }
    switch (size) {
    case 1:
      return 0;
    case 2:
      return 1;
    case 3:
      return 2;
    case 4:
      return 2;
    case 5:
      return 3; //not fully sure whether these assignemnts are right
    case 6:
      return 3;
    case 7:
      return 3;
    case 8:
      return 3;
    case 9: //return the number of qubits involved depending on the length of feynmanPaths, but i doubt that this is correct
      return 4;
    case 10:
      return 4;
    case 11:
      return 4;
    case 12:
      return 4;
    case 13:
      return 4;
    case 14:
      return 4;
    case 15:
      return 4;
    case 16:
      return 4;
    default:
      // Not supported.
      return 42;
    }
  }

  template <typename StateSpace>
  static bool CreateStates(unsigned num_qubits0,unsigned num_qubits1,
                           const StateSpace& state_space, bool create,
                           typename StateSpace::State& state0,
                           typename StateSpace::State& state1,
                           typename StateSpace::State* (&rstate0),
                           typename StateSpace::State* (&rstate1)) {
    if (create) {
      state0 = state_space.Create(num_qubits0);
      state1 = state_space.Create(num_qubits1);

      if (state_space.IsNull(state0) || state_space.IsNull(state1)) {
        IO::errorf("not enough memory: is the number of qubits too large?\n");
        return false;
      }

      rstate0 = &state0;
      rstate1 = &state1;
    }

    return true;
  }

  For for_;
};

}  // namespace qsim

#endif  // HYBRID_H_
