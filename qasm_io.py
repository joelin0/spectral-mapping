import re


def qasm_to_2Q_blocks(filename):
    # absorbs 1Q gates in with closest previous 2Q gate
    # except for left-most layer 1Q gates, which are absorbed into closest next 2Q gate
    with open(filename) as f:
        lines = f.readlines()
    qubit_label = re.compile("qreg (.*)\[(.*)\]").findall(lines[2])[0][0]
    interactions = []
    gate_blocks = []
    latest_qubit_interaction_index = {}
    qubit_to_unassigned_gate = {}
    unique_qubits = {}
    nqubits = 0
    lines = lines[4:]
    ninteractions = 0
    for i in range(len(lines)):
        res = re.compile("(.*) " + qubit_label + "\[(.*)\]," + qubit_label + "\[(.*)\]").findall(lines[i])
        if not res:
            res = re.compile("(.*) " + qubit_label + "\[(.*)\]").findall(lines[i])
            gate, qlabel = res[0]
            if qlabel not in unique_qubits:
                unique_qubits[qlabel] = nqubits
                nqubits += 1
            q = unique_qubits[qlabel]
            if q not in latest_qubit_interaction_index:
                qubit_to_unassigned_gate.setdefault(q, []).append((gate, q))
            else:
                gate_blocks[latest_qubit_interaction_index[q]].append((gate, q))
        else:
            gate, qlabel1, qlabel2 = res[0]
            if qlabel1 not in unique_qubits:
                unique_qubits[qlabel1] = nqubits
                nqubits += 1
            if qlabel2 not in unique_qubits:
                unique_qubits[qlabel2] = nqubits
                nqubits += 1

            q1, q2 = unique_qubits[qlabel1], unique_qubits[qlabel2]
            if q1 in latest_qubit_interaction_index and q2 in latest_qubit_interaction_index \
                    and latest_qubit_interaction_index[q1] == latest_qubit_interaction_index[q2]:
                gate_blocks[latest_qubit_interaction_index[q1]].append((gate, q1, q2))
            else:
                block = qubit_to_unassigned_gate.pop(q1, []) + qubit_to_unassigned_gate.pop(q2, [])
                block.append((gate, q1, q2))
                gate_blocks.append(block)
                interactions.append((q1, q2))
                latest_qubit_interaction_index[q1] = ninteractions
                latest_qubit_interaction_index[q2] = ninteractions
                ninteractions += 1
    if qubit_to_unassigned_gate:
        raise RuntimeError("Single qubit gate unabsorbed")
    return interactions, gate_blocks, len(unique_qubits)


def write_compliant_circuit_to_open_qasm(gates, filename, nqubits):
    with open(filename, 'w') as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{nqubits}];\ncreg c[{nqubits}];\n")
        for gate, *qs in gates:
            line = gate + ' '
            for q in qs:
                line += f"q[{q}],"
            line = line[:-1] + ';\n'
            f.write(line)


def quipper_to_qasm(filename):
    # can use
    # tdg for t-dagger
    # rz(phi) for z rotation
    # x for not
    # cx for cnot
    # sdg for s-dagger

    gate_translate = {"not": "x",
                      "H": "h",
                      "T": "t",
                      "S": "s",
                      "Z": "z"}
    dagger_str = "dg"
    control_str = "c"

    with open(f"../optimizer/Arithmetic_and_Toffoli/{filename}") as f:
        in_lines = f.readlines()

    qubits = set()
    out_lines = []

    for i in range(1, len(in_lines) -1):
        line = in_lines[i]
        res = re.compile("QGate\\[\"(.*?)\"\\]([*]?)\((.*?)\) (with controls=\[\+(.*?)\] )?with nocontrol").findall(line)[0]
        in_gate = res[0]
        out_gate = gate_translate[in_gate]
        conjugated = res[1] == '*'
        target = int(res[2])
        control = None if res[4] == '' else int(res[4])
        if conjugated:
            out_gate += dagger_str
        if control is not None:
            out_gate = control_str + out_gate

        qubits.add(target)
        if control is None:
            out_lines.append(f'{out_gate} q[{target}];\n')
        else:
            out_lines.append(f'{out_gate} q[{control}],q[{target}];\n')
            qubits.add(control)

    nqubits = len(qubits)
    with open(f'examples/{filename}.qasm', 'w') as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{nqubits}];\ncreg c[{nqubits}];\n")
        f.writelines(out_lines)