# spectral-mapping
Spectral Tools for Mapping Logical Qubits of Circuits onto Connectivity-Constrained Physical Architectures

## Requirements
This code has been tested using Python 3.6.5 on MacOS.

The code makes use of `networkx`, `numpy`, `scipy`, and `matplotlib`. Quick installation of these packages can be performed by running

```pip install -r requirements.txt```

## Usage
To generate an equivalent SWAP-compliant circuit for a given example circuit, run the command
```
python3 main.py <filename> <benchmark_folder> <result_folder>
```
For example, if you have a benchmark in `benchmarks/example.qasm` and wish to output the results into a folder `results`, run
```
python3 main.py example.qasm benchmarks results
```

Note that the input file is to be given in OpenQASM format and should only contain single- or two-qubit gates. 

The connectivity compliant circuits generated this algorithm are written to the `<result_folder>` provided. 
