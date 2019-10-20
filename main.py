import sys
import time
import os
from qasm_io import *
from mapper import generate_compliant_circuit

# This file is run as python3 main.py filename benchmark_folder result_folder
# which get path compliant circuits using a variety of spectral options
# and also prints out total cpu time used
if __name__ == '__main__':
    start = time.process_time()

    filename = sys.argv[1]
    benchmark_folder = sys.argv[2]
    path = f'{benchmark_folder}/{filename}'

    interactions, gate_blocks, nqubits = qasm_to_2Q_blocks(path)
    num_original_gates = sum(len(block) for block in gate_blocks)
    forward_only = True
    fallback_only = False
    m, n = nqubits, 1
    best_swaps = float('inf')

    configs = [(5, 1), (8, 6), (7, 1), (2, 3), (9, 9), (4, 1), (8, 1), (3, 4), (5, 6), (8, 2)]

    result_folder = sys.argv[3]
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)

    root, ext = os.path.splitext(filename)
    for i, j in configs:
        start_perm, gates = generate_compliant_circuit(interactions, gate_blocks, nqubits, m, n, i / 10, j / 10, forward_only, fallback_only)
        result_path = f'{result_folder}/{root}_{forward_only}_{fallback_only}_{i}_{j}.qasm'
        write_compliant_circuit_to_open_qasm(gates, result_path, nqubits)
        nswaps = len(gates) - num_original_gates
        best_swaps = min(nswaps, best_swaps)
    cputime = time.process_time() - start
    print(f'{filename},{best_swaps},{cputime}')