import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Union
from cycler import cycler
import qiskit as qk
import itertools
import sys
import matplotlib.font_manager as fm

try: 
    font_path = '/mnt/c/Windows/Fonts/TIMES.TTF'
    fm.fontManager.addfont(font_path)
    prop = fm.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = prop.get_name()

    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = prop.get_name()  # Roman (normal) font
    
    plt.rcParams['font.family'] = 'Times New Roman'

    colorlen = 10 # number of colors
    plt.rcParams["axes.prop_cycle"] = cycler('color', plt.cm.rainbow(np.linspace(0,1,colorlen)) )
except:
    pass


def readout_cut(filename:str) -> dict:
    """reads out the data for block."""

    extracted_data = {}

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('Number of Feynman Paths'):
                # Extract number of Feynman Paths
                paths = int(line.split(':')[-1].strip())
                extracted_data['Paths'] = paths
            elif line.startswith('time elapsed'):
                # Extract time elapsed
                time_elapsed = float(line.split()[2])
                extracted_data['Time Elapsed (seconds)'] = time_elapsed
            elif line.startswith("Time for Run: "):
                time_sim = float(line.split()[3])
                extracted_data['time sim'] = time_sim

    return extracted_data

def readout_nocut(file_path:str) -> dict:
    """reads out the times without cut, with qsim"""

    extracted_data = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace
            
            if line.startswith('time elapsed'):
                # Extract the time value
                time_value = float(line.split()[2])
                extracted_data['Time Elapsed (seconds)']=time_value
            elif line.startswith('Time for ApplyFusedGates: '):
                time_value = float(line.split()[3])
                extracted_data["ApplyFusedGates"]=time_value

    return extracted_data
        
def read_file(file_path):
    """should be used only for the amplitude output files"""
    with open(file_path, 'r') as file:
        data = []
        for line in file:
            # Split the line into two columns and convert to float
            columns = line.split()
            if len(columns) == 2:
                data.append([float(columns[0]), float(columns[1])])
        return np.array(data)

def compare_files(file1, file2, tolerance=1e-7):
    """comparison of amplitude output files"""
    data1 = read_file(file1)
    data2 = read_file(file2)

    if data1.shape != data2.shape:
        print("Files have different shapes.")
        return False

    counter_falses = 0

    all_fine = True
    for i, j in itertools.product(range(data1.shape[0]), range(data1.shape[1])):
        diff = abs(data1[i, j] - data2[i, j])
        if diff > tolerance:
            print(f"Issue occurred at element ({i}, {j}) with difference {diff}")
            all_fine = False
            counter_falses += 1
            if counter_falses == 100:
                break

    if all_fine:
        print(f"All fine within the given tolerance {tolerance}")

    return all_fine


def visualize_circ(path:str):
    """visualizes a circuit from a qsimh input file. this means we turn the contents of the file into a qiskit object and plot it"""
    
    inblock=False
    with open(path, 'r') as file:
        lines = file.readlines()
        for i,line in enumerate(lines):
            if i==0:
                num_qubits = int(line.strip())
                circ = qk.QuantumCircuit(num_qubits)
            else:
                if line.strip():
                    if ']' in line:
                        inblock=True
                    line = line.replace('[', '')
                    line = line.replace(']', '')
                    parts = line.split()
                    gate = parts[1]
                    
                    #turn "gate" into something qiskit knows
                    if gate=="cnot":
                        gate="cx"
                    elif gate=="sw":
                        gate="swap"
                    elif gate=="is":
                        gate="iswap"
                    elif gate=="fs":#even if this is wrong, just do antoher gate to display the structure
                        gate="cry"

                    print(parts)

                    #add to circ
                    if len(parts)-2 == 1: #single qubit gate
                        getattr(circ, gate)(int(parts[2]))
                    elif len(parts)-2 == 2: #two qubit gate without angles or single qubit with angle
                        try:
                            getattr(circ, gate)(int(parts[2]), int(parts[3])) #two qubit no angles
                        except:
                            getattr(circ, gate)(float(parts[3]), int(parts[2])) #single qubit angle
                    elif len(parts)-2 == 3: #2 gate wiht angle
                        getattr(circ, gate)(float(parts[4]), int(parts[2]), int(parts[3]))
                    elif len(parts)-2 == 4: #2 gate wiht angles
                        getattr(circ, gate)(float(parts[4]), int(parts[2]), int(parts[3])) #this case only happens for the fsim gate which appears not to be present in qiskit,hence for illustration just take another gate with a single angel
                    
                    if inblock==True:
                        circ.barrier() 
                        inblock=False  

    circuit_diagram = circ.draw(output='mpl')
    circuit_diagram.savefig(path + "diagram.pdf")      
    print("stored: "+ path + "diagram.pdf")     

def print_diffs(reps, path_nocut, path_block, path_noblock, suffix=""):
    """prints the results of the runs of a single circuit.
    note that the paths only include everything up to `_rep{rep}.tlog` because we need to loop over this in this functio here
    in case your files have a suffix, include the suffix"""
    full_nocut = []
    full_block = []
    full_noblock = []

    sim_nocut = []
    sim_block = []
    sim_noblock = []

    paths_block = 0
    paths_noblock = 0


    for rep in range(reps):
        times_nocut = path_nocut + f"_rep{rep}" + suffix + ".tlog"
        times_block = path_block + f"_rep{rep}" + suffix + ".tlog"
        times_noblock = path_noblock + f"_rep{rep}" + suffix + ".tlog"
        dct = readout_nocut(times_nocut)
        full_nocut.append(dct['Time Elapsed (seconds)'])
        sim_nocut.append(dct["ApplyFusedGates"])

        dct = readout_cut(times_block)
        full_block.append(dct['Time Elapsed (seconds)'])
        sim_block.append(dct['time sim'])

        paths_block = dct["Paths"]

        flag_noblock = False
        try:
            dct = readout_cut(times_noblock)
            full_noblock.append(dct['Time Elapsed (seconds)'])
            sim_noblock.append(dct['time sim'])

            paths_noblock = dct["Paths"]
        except:
            flag_noblock = True
            print("No block was obviously timed out.")

    print("--------Full times---------")
    print(f"No cut: Mean ={np.mean(full_nocut)}, std={np.std(full_nocut)}")
    print(f"Block: Mean ={np.mean(full_block)}, std={np.std(full_block)}")
    if flag_noblock == False:
        print(f"No Block: Mean ={np.mean(full_noblock)}, std={np.std(full_noblock)}")

    print("--------Sim times---------")
    print(f"No cut: Mean ={np.mean(sim_nocut)}, std={np.std(sim_nocut)}")
    print(f"Block: Mean ={np.mean(sim_block)}, std={np.std(sim_block)}")
    if flag_noblock == False:
        print(f"No Block: Mean ={np.mean(sim_noblock)}, std={np.std(sim_noblock)}")

    print("-----Paths-----")
    print(f"Block: {paths_block}")
    if flag_noblock == False:
        print(f"No Block: {paths_noblock}")


    print("--------Ratios--------")
    print(f"S/J = {np.mean(full_nocut) / np.mean(full_block)}")
    if flag_noblock == False:
        print(f"T/J = {np.mean(full_noblock) / np.mean(full_block)}")
    else:
        print(f"T/J >= {3600 / np.mean(full_block)}")


