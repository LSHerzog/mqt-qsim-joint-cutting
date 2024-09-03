# Joint Cutting for Hybrid Schrödinger Feynman Simulation (`qsimh`)

Everything one needs to know about using the qsimh package with joint cutting is described in this file.

## Additional Dependencies

In order to generally perform the Schmidt Decomposition on blocks of gates (i.e. doing the cut fully automatically, without having to consider analytical forms beforehand) one has to perform a SVD. For this one needs [`xtensor`](https://xtensor.readthedocs.io/en/latest/installation.html) as well as [`xtensor-blas`](https://xtensor-blas.readthedocs.io/en/latest/installation.html). One has to ensure that those dependencies are properly recognized by `qsim`.

## Usage and Examples

In order to perform joint cutting, one has to gather the corresponding gates into blocks in the input file. Consider running `qsimh` with the following command:

```
./qsimh_amplitudes.x -c../circuits/circuit_q18 -d 1000 -k 0,1,2,3,4,5,6,7,8 -t 8 -w 0 -p 0 -r1 -i ../circuits/bitstrings_q18 -o ../circuits/ampl_q18_block -v 1 
```

This means that there are 18 qubits in the circuit and the cut is performed after qubit 8. Standard Cutting would be performed if circuit_q18 was chosen as

```
18
0 h 0
0 h 1
0 h 2
0 h 3
0 h 4
0 h 5
0 h 6
0 h 7
0 h 8
0 h 9
0 h 10
0 h 11
0 h 12
0 h 13
0 h 14
0 h 15
0 h 16
0 h 17
1 rzz 8 9 0.2
2 rzz 8 10 0.5
3 rzz 8 11 2.4
4 rzz 8 14 2.5
5 rzz 8 17 1.9
```

However, it would be suitable to summarize the RZZ gates into a block. The parsing is adapted such that one just has to put brackets around the gates in the block and voilá, qsimh cuts those gates jointly now.

```
18
0 h 0
0 h 1
0 h 2
0 h 3
0 h 4
0 h 5
0 h 6
0 h 7
0 h 8
0 h 9
0 h 10
0 h 11
0 h 12
0 h 13
0 h 14
0 h 15
0 h 16
0 h 17
[1 rzz 8 9 0.2
2 rzz 8 10 0.5
3 rzz 8 11 2.4
4 rzz 8 14 2.5
5 rzz 8 17 1.9]
```

> ⚠️ **Warning:** In general, be aware that `qsim` has issues with spaces after the parameters of a gate per line. This will cause issues during the parsing. Specifically regarding the joint cutting, note that the brackets must be placed as in the above example. Something like e.g.:
> 
> ```
> [
> 3 rzz 8 11 2.4
> 4 rzz 8 14 2.5
> ]
> ```
> 
> cannot be parsed properly.


## Limitations

It is important to know that one cannot have too many qubits within a partition of a block. A hard limit is set such that maximally 6 qubits can be in a partition, i.e., a block can contain maximally 12 qubits in principle. However, in practice, it turned out that 5 qubits per partition of a block, i.e. maximally 10 qubits in a block is more stable. Thus, do not go beyond this value as otherwise your amplitudes might come out incorrect.

So far, the `cirq` Python Bindings were not adapted yet. Thus, joint cutting can only be performed on the C++ Level.

## Timestamps 

In order to check the performance of qsimh with and without joint cutting, additional timestamps were added. The most important time stamps can be read out using ``/benchmarks_join_cutting/readout.py`. This uses the output files from the qsim(h) computation of multiple runs and yields an output like this

```
--------Full times---------
No cut: Mean =43.17724, std=0.14887001847249026
Block: Mean =0.6180796000000001, std=0.02119766231073605
No Block: Mean =118.02239999999999, std=2.040881241032902
--------Sim times---------
No cut: Mean =42.91472, std=0.15012290165061354
Block: Mean =0.6176368, std=0.021206826772527762
No Block: Mean =118.02239999999999, std=2.040881241032902
-----Paths-----
Block: 2048
No Block: 131072
--------Ratios--------
S/J = 69.85708636881074
T/J = 190.9501624062661
```

Here, `no cut` represents the Schrödinger Simulation with `qsim`. `Block` denotes the joint cutting for `qsimh` and `No Block` the traditional, single cutting. Everything below `Sim times` shows the runtimes for the simulation only - meaning the time for executing `HybridSimulator().Run`, excluding the parsing, SVDs and fusing. It only includes allocating the cut partitions and building the paths as well as executing them and adding them together. `Full times` includes both the latter as well as the first. Also for the Schrödinger simulation `sim times` and `full times` is distinguished as the first only includes only the execution time of all `ApplyFusedGate` runs. It excludes the parsing and the fusing of the gates, the latter includes all of those runtimes.

As the output files of multiple runs are given to this function, it computes the mean as well as the standard deviation for each variant.

In more detail, one can find a more verbose variety of timestamps in the output files directly. For qsim (i.e. Schrödinger Simulation):

| Timestamp   | Description   | source |
|-------------|---------------|--------|
| Time for Parsing  | Complete parsing time of single gates | `circuit_qsim_parser.h`|
| init time is | initialization of `CreateStateSpace` etc. | `run_qsim.h`|
| fuse time is | fusing gates | `run_qsim.h`|
| Time for ApplyFusedGates | Simulation runtime as described above | `run_qsim.h`|
| time is | full time without init time | `run_qsim.h`|
| time elapsed | full time as described above | `run_qsim.h`|

and for qsimh:

| Timestamp in output   | Description   | source |
|-----------------------|---------------|--------|
| Time for Multiplying Unitaries (gateBlock.Create) | Time for creating the joint gates of a block by multiplying and removing empty qubits | `circuit_qsim_parser.h` |
| Time for ParseBlockGates | Reading the gates from the file and gathering the single gates in vectors | `circuit_qsim_parser.h` |
| Time for Parsing | Complete parsing time (ParseBlockGates inclusive)| `circuit_qsim_parser.h`|
| Time for SchmidtDecomposition | Time for performing the Schmidt Decomposition, i.e. reshaping legs as well as SVD via `xtensor-blas`, output occurs for each joint cut | `hybrid.h` |
| Time for SplitLattice | allocates and constructs the different paths (not run yet) | `run_qsimh.h`|
| Time for Fuser | Fusing the gates (both cut and uncut) | `run_qsimh.h`|
| Time for Run | simulation time as described above | `run_qsimh.h`|
| time elapsed | full runtime as described above | `run_qsimh.h`|


## Benchmarks and Tutorial Notebooks

### Tutorial

A Jupyter Notebook demonstrating the basic usage of qsimh with joint cutting can be found in `benchmarks_joint_cutting/toy_examples`.

### Benchmarks

TODO: Include Benchmarks for weighted Max Cut (+ helper functions) as well as random Boixo circuits.

## Summary of Adaptions in `qsimh`

todo