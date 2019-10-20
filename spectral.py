import networkx as nx
from scipy.linalg import eigh
import numpy as np
from random import random


def labeling_from_graph_laplacian(G, m, n, nqubits, previous_mapping):
    # use generate labelings and then get the one closest to previous_mappings
    # if previous_mapping is None, just pick one
    mapping_generator = generate_equivalent_labelings_from_graph_laplacian(G, m, n, nqubits)
    return get_best_orientation(previous_mapping, mapping_generator)


# TODO: simpler may be to use networkx.linalg.algebraicconnectivity.fiedler_vector
# rather than generate it on my own
def generate_equivalent_labelings_from_graph_laplacian(G, m, n, nqubits):
    if m < n:
        m, n = n, m

    L = nx.laplacian_matrix(G,range(nqubits)).todense()

    # eigvals, eigvecs = eig(L)
    # idx = eigvals.argsort()  # order eigenvalues smallest to largest
    # eigvals = eigvals[idx]  # order eigenvalues smallest to largest
    # eigvecs = eigvecs[:,idx]  # order eigenvectors, stored in each columns, according to eigenvalues
    # v1, v2 = eigvecs[:,1], eigvecs[:,2]

    eigvals, eigvecs = eigh(L, overwrite_a=True, overwrite_b=True, check_finite=False, eigvals=(1, 2))
    v1, v2 = eigvecs[:,0], eigvecs[:,1]
    yield labeling_from_vectors(v1, v2, m, n, nqubits)
    yield labeling_from_vectors(-v1, v2, m, n, nqubits)
    if n > 1:
        yield labeling_from_vectors(v1, -v2, m, n, nqubits)
        yield labeling_from_vectors(-v1, -v2, m, n, nqubits)
        if m == n:
            yield labeling_from_vectors(v2, v1, m, n, nqubits)
            yield labeling_from_vectors(-v2, v1, m, n, nqubits)
            yield labeling_from_vectors(v2, -v1, m, n, nqubits)
            yield labeling_from_vectors(-v2, -v1, m, n, nqubits)


def labeling_from_vectors(v1, v2, m, n, nqubits):
    order1 = v1.argsort()
    mappings = [None] * nqubits
    for i in range(0, len(order1), n):
        row = order1[i:i+n]  # take n of the components at a time, in order, for that given row
        row_vec_comp = v2[row]  # now take their corresponding v2 components to order along a column
        order2 = row_vec_comp.argsort()
        for j, row_idx in enumerate(order2):
            logical_qubit = row[row_idx]
            physical_location = (i // n, j)
            mappings[logical_qubit] = physical_location

    return mappings


def get_best_orientation(previous_mapping, mapping_generator):
    # given a sample mapping, look at all possible symmetries and select the one with
    # the closest total L1 distance from the previous mapping
    # if previous_mapping is None, just use sample_mapping
    closest = float('inf')
    best_mapping = None
    for mapping in mapping_generator:
        if previous_mapping is None:
            best_mapping = mapping
            break
        distance = 0
        for pos1, pos2 in zip(mapping, previous_mapping):
            distance += l1_norm(pos1, pos2)
        if distance < closest:
            closest = distance
            best_mapping = mapping

    return best_mapping, {best_mapping[key]: key for key in range(len(best_mapping))}


def generate_equivalent_labelings(mappings, m, n):
    # originally generated mapping
    yield mappings

    # calculate horizontal mirror flip
    horizontal_flipped_mappings = [(m - 1 - i, j) for (i, j) in mappings]
    yield horizontal_flipped_mappings

    # 2D case
    if n > 1:
        # vertical flip
        vertical_flipped_mappings = [(i, n - 1 - j) for (i, j) in mappings]
        yield vertical_flipped_mappings

        # horiz + vert flip = 180 degree rotation
        vertical_horizontal_flipped_mappings = [(i, n - 1 - j) for (i, j) in horizontal_flipped_mappings]
        yield vertical_horizontal_flipped_mappings

        # square case
        if m == n:
            # for each of the four above cases, just switch the x and y axes
            swapped_mappings = [(j, i) for (i, j) in mappings]
            yield swapped_mappings

            swapped_horizontal_flipped_mappings = [(j, i) for (i, j) in horizontal_flipped_mappings]
            yield swapped_horizontal_flipped_mappings

            swapped_vertical_flipped_mappings = [(j, i) for (i, j) in vertical_flipped_mappings]
            yield swapped_vertical_flipped_mappings

            swapped_vertical_horizontal_flipped_mappings = [(j, i) for (i, j) in vertical_horizontal_flipped_mappings]
            yield swapped_vertical_horizontal_flipped_mappings


def generate_labeling_from_graph_laplacian(G, m, n, nqubits):
    if m < n:
        m, n = n, m

    L = nx.laplacian_matrix(G,range(nqubits)).todense()

    # eigvals, eigvecs = eig(L)
    # idx = eigvals.argsort()  # order eigenvalues smallest to largest
    # eigvals = eigvals[idx]  # order eigenvalues smallest to largest
    # eigvecs = eigvecs[:,idx]  # order eigenvectors, stored in each columns, according to eigenvalues
    # v1, v2 = eigvecs[:,1], eigvecs[:,2]

    eigvals, eigvecs = eigh(L, overwrite_a=True, overwrite_b=True, check_finite=False, eigvals=(1, 2))
    v1, v2 = eigvecs[:,0], eigvecs[:,1]
    # TODO: make this bottom part its own function
    # then to generate all mappings, including symmetry
    # just pass in (v1, v2), (-v1, v2), etc.
    # rather than performing symmetry on the final product
    # this way the "holes" are guaranteed to be in the bottom right corner
    # for any dimensions
    order1 = v1.argsort()
    mappings = [None] * len(L)
    for i in range(0, len(order1), n):
        row = order1[i:i+n]  # take n of the components at a time, in order, for that given row
        row_vec_comp = v2[row]  # now take their corresponding v2 components to order along a column
        order2 = row_vec_comp.argsort()
        for j, row_idx in enumerate(order2):
            logical_qubit = row[row_idx]
            physical_location = (i // n, j)
            mappings[logical_qubit] = physical_location

    return mappings


def lazy_mapping_from_node_coordinates(coords, node_to_qubits, qubit_to_node, nqubits, m, n):
    # just go back from nodes to qubits and map accordingly
    noise1 = [10**-6 * random() for node in node_to_qubits] # this is meant to perturb nodes off of each other randomly
    noise2 = [10 ** -6 * random() for node in
              node_to_qubits]  # this is meant to perturb nodes off of each other randomly
    v1 = np.array([coords[qubit_to_node[q]][0] + noise1[qubit_to_node[q]] for q in range(nqubits)])
    v2 = np.array([coords[qubit_to_node[q]][1] + noise2[qubit_to_node[q]] for q in range(nqubits)])
    order1 = v1.argsort()
    mappings = [None] * nqubits
    for i in range(0, len(order1), n):
        # check node pair overflow, only for 2D case
        if n > 1 and i + n < nqubits:
            q1, q2 = order1[i+n-1], order1[i+n]
            if qubit_to_node[q1] == qubit_to_node[q2]:
                swap_idx = i + n + 1
                while len(node_to_qubits[qubit_to_node[order1[swap_idx]]]) > 1:
                    swap_idx += 1
                q3 = order1[swap_idx]
                order1 = np.delete(order1, swap_idx)
                order1 = np.insert(order1, i + n - 1, q3)
        row = order1[i:i + n]  # take n of the components at a time, in order, for that given row

        row_vec_comp = v2[row]  # now take their corresponding v2 components to order along a column
        order2 = row_vec_comp.argsort()
        for j, row_idx in enumerate(order2):
            logical_qubit = row[row_idx]
            physical_location = (i // n, j)
            mappings[logical_qubit] = physical_location

    # mappings is logical --> physical
    return mappings



def graph_from_rev_perm(rev_perm, m, n, perm_weight):
    # rev_perm maps physical location to logical qubit
    # return laplacian for logical qubit
    G = nx.Graph()
    for i in range(m):
        for j in range(n):
            q1 = rev_perm.get((i, j))
            if q1 is None:
                continue
            right = (i + 1, j)
            q2 = rev_perm.get(right)
            if q2 is not None:
                G.add_edge(q1, q2, weight=perm_weight)
            down = (i, j + 1)
            q3 = rev_perm.get(down)
            if q3 is not None:
                G.add_edge(q1, q3, weight=perm_weight)
    return G


def l1_norm(pos1, pos2):
    norm = 0
    for i, j in zip(pos1, pos2):
        norm += abs(i - j)
    return norm


def increment_edge_weight(G, i, j, weight):
    if G.has_edge(i, j):
        G[i][j]['weight'] += weight
    else:
        G.add_edge(i, j, weight=weight)


def node_coordinates_from_laplacian(G, m, n):
    if m < n:
        m, n = n, m

    nnodes = len(G.nodes)
    L = nx.laplacian_matrix(G,range(nnodes)).todense()
    if n == 1:
        eigvals, eigvecs = eigh(L, overwrite_a=True, overwrite_b=True, check_finite=False, eigvals=(1, 1))
        v1 = eigvecs[:,0]
        return [(v1[i], 0) for i in range(nnodes)]
    else:
        eigvals, eigvecs = eigh(L, overwrite_a=True, overwrite_b=True, check_finite=False, eigvals=(1, 2))
        v1, v2 = eigvecs[:,0], eigvecs[:,1]
        return [(v1[i], v2[i]) for i in range(nnodes)]


def location_to_label(i, j, m, n):
    return i + m * j
