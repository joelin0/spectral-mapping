import networkx as nx
from token_swapping import swap_between_perm
from spectral import *

# TODO: refactor all mentions of CNOT (and cnot variables) to interactions instead
class CNOTContainer(object):
    def __init__(self, cnots):
        """

        :param cnots: a list of CNOTs
        :return:
        dependencies: a dictionary, keyed by cnot index. values are another dictionary, including a "left" set of
        cnots that are left adjacent, a "right" set of cnots that are right adjacent, a "left_layer" and a "right_layer"
        left_layers_to_cnot: dictionary mapping each left layer to the cnot indices within them
        right_layers_to_cnot: dictionary mapping each right layer to the cnot indices within them
        """
        self.dependencies = {}
        self.left_layers_to_cnot = {}
        self.right_layers_to_cnot = {}
        self.applied_indices = set()
        self.cnots = tuple(cnots)  # frozen
        self.last_layer = 0

        qubit_to_latest_cnot = {}
        for i, (q1, q2) in enumerate(cnots):
            left_q1 = qubit_to_latest_cnot.get(q1, None)
            left_q2 = qubit_to_latest_cnot.get(q2, None)
            lefts = set()
            if left_q1 is not None:
                lefts.add(left_q1)
                self.dependencies[left_q1]['right'].add(i)
            if left_q2 is not None:
                lefts.add(left_q2)
                self.dependencies[left_q2]['right'].add(i)
            self.dependencies[i] = {}
            self.dependencies[i]['left'] = lefts
            self.dependencies[i]['right'] = set()
            qubit_to_latest_cnot[q1] = i
            qubit_to_latest_cnot[q2] = i
        ncnots = len(cnots)

        for i in range(ncnots):
            lefts = self.dependencies[i]['left']
            layer = 0
            for left in lefts:
                layer = max(self.dependencies[left]['left_layer'] + 1, layer)
            self.dependencies[i]['left_layer'] = layer
            self.last_layer = max(self.last_layer, layer)
            self.left_layers_to_cnot.setdefault(layer, set()).add(i)

        for j in range(ncnots):
            i = ncnots - j - 1
            rights = self.dependencies[i]['right']
            layer = 0
            for right in rights:
                layer = max(self.dependencies[right]['right_layer'] + 1, layer)
            self.dependencies[i]['right_layer'] = layer
            self.right_layers_to_cnot.setdefault(layer, set()).add(i)

    def apply_permutation(self, mappings, left):
        # mappings goes from a logical qubit to a location on the grid
        # left specifies whether we are looking at applying leftmost or rightmost cnots

        applied = set()
        # finished applying everything
        if self.last_layer < 0:
            return applied

        i = 0
        to_apply = list(self.left_layers_to_cnot[0]) if left else list(self.right_layers_to_cnot[0])

        last_layer_touched = 0

        while i < len(to_apply):
            cnot_idx = to_apply[i]
            q1, q2 = self.cnots[cnot_idx]
            pos1, pos2 = mappings[q1], mappings[q2]
            if l1_norm(pos1, pos2) <= 1.05:
                applied.add(cnot_idx)
                self.applied_indices.add(cnot_idx)

                info = self.dependencies[cnot_idx]
                del self.dependencies[cnot_idx]
                self.left_layers_to_cnot[info['left_layer']].remove(cnot_idx)
                self.right_layers_to_cnot[info['right_layer']].remove(cnot_idx)

                # on left updates, the children are to the right
                # on right updates, the children are to the left
                adj = info['right'] if left else info['left']
                for child in adj:
                    child_info = self.dependencies[child]
                    adj_parents = child_info['left'] if left else child_info['right']
                    adj_parents.remove(cnot_idx)
                    if len(adj_parents) == 0:
                        to_apply.append(child)
                        last_layer_touched = max(last_layer_touched,
                                                 child_info['left_layer'] if left else child_info['right_layer'])
            i += 1

        # update layer data
        # if it was a left update, only left layer data needs to be updated
        # as no rightward dependencies change
        # and vice versa
        layer_data = self.left_layers_to_cnot if left else self.right_layers_to_cnot
        layer = 0
        while True:
            changed = False
            for cnot_idx in set(layer_data.get(layer, set())):
                new_layer = 0
                info = self.dependencies[cnot_idx]
                parents = info['left'] if left else info['right']
                for parent in parents:
                    parent_info = self.dependencies[parent]
                    new_layer = max(new_layer, (parent_info['left_layer'] if left else parent_info['right_layer']) + 1)
                if new_layer != layer:
                    layer_data[layer].remove(cnot_idx)
                    layer_data[new_layer].add(cnot_idx)
                    info['left_layer' if left else 'right_layer'] = new_layer
                    changed = True
            if not changed and layer > last_layer_touched:
                break
            layer += 1

        for layer in range(self.last_layer, -1, -1):
            if len(self.left_layers_to_cnot[layer]) == 0:
                del self.left_layers_to_cnot[layer]
                del self.right_layers_to_cnot[layer]
                self.last_layer -= 1

        return applied

    def __str__(self):
        return str(self.dependencies)


class Mapper(object):
    def __init__(self, container):
        self.container = container

    def mapper(self, m, n, nqubits, left):
        raise NotImplemented()

# Note: change as of 5/8/2019 to try Greedy Layer choices but with Lazy Layer discounting
class GreedyLayerLimitedMapper(Mapper):
    def __init__(self, container, nlayers):
        Mapper.__init__(self, container)
        self.nlayers = nlayers

    def mapper(self, m, n, nqubits, left, prev_perm=None, discount=0.5, prev_perm_weight=0.5):
        G = nx.Graph() if prev_perm is None else graph_from_rev_perm(prev_perm, m, n, prev_perm_weight)
        layer_info = self.container.left_layers_to_cnot if left else self.container.right_layers_to_cnot
        lazy_layer = 'right_layer' if left else 'left_layer'
        for q in range(nqubits):
            G.add_node(q)
        for i in range(self.nlayers):
            if i not in layer_info:
                break
            for cnot in layer_info[i]:
                weight = discount ** (self.container.last_layer - self.container.dependencies[cnot][lazy_layer])
                q1, q2 = self.container.cnots[cnot]
                increment_edge_weight(G, q1, q2, weight)

        prev_mapping = None
        if prev_perm is not None:
            prev_mapping = [None] * nqubits
            for key in prev_perm:
                prev_mapping[prev_perm[key]] = key

        return labeling_from_graph_laplacian(G, m, n, nqubits, prev_mapping)


# Note: change as of 5/8/2019 to try Greedy Layer choices but with Lazy Layer discounting
class GreedyForcedPairMapper(Mapper):
    def __init__(self, container, nlayers):
        Mapper.__init__(self, container)
        self.nlayers = nlayers

    def mapper(self, m, n, nqubits, left, prev_perm=None, discount=0.5, prev_perm_weight=0.5):
        node_to_qubits, qubit_to_node = self.node_mapping_with_forced_pairs(nqubits, left)
        G = self.generate_connectivity_with_forced_pairs(m, n, left, prev_perm, discount,
                                                         prev_perm_weight, node_to_qubits, qubit_to_node)
        coords = node_coordinates_from_laplacian(G, m, n)

        prev_mapping = None
        if prev_perm is not None:
            prev_mapping = [None] * nqubits
            for key in prev_perm:
                prev_mapping[prev_perm[key]] = key

        return get_best_orientation(prev_mapping,
                                    self.mapping_generator(coords, node_to_qubits, qubit_to_node, nqubits, m, n))

    def mapping_generator(self, coords, node_to_qubits, qubit_to_node, nqubits, m, n):
        yield lazy_mapping_from_node_coordinates(coords, node_to_qubits, qubit_to_node, nqubits, m, n)
        yield lazy_mapping_from_node_coordinates([(-c[0], c[1]) for c in coords], node_to_qubits, qubit_to_node, nqubits, m, n)
        if n > 1:
            yield lazy_mapping_from_node_coordinates([(c[0], -c[1]) for c in coords], node_to_qubits, qubit_to_node,
                                                     nqubits, m, n)
            yield lazy_mapping_from_node_coordinates([(-c[0], -c[1]) for c in coords], node_to_qubits, qubit_to_node,
                                                     nqubits, m, n)
            if m == n:
                yield lazy_mapping_from_node_coordinates([(c[1], c[0]) for c in coords], node_to_qubits, qubit_to_node, nqubits, m, n)
                yield lazy_mapping_from_node_coordinates([(-c[1], c[0]) for c in coords], node_to_qubits, qubit_to_node,
                                                         nqubits, m, n)
                yield lazy_mapping_from_node_coordinates([(c[1], -c[0]) for c in coords], node_to_qubits, qubit_to_node,
                                                         nqubits, m, n)
                yield lazy_mapping_from_node_coordinates([(-c[1], -c[0]) for c in coords], node_to_qubits, qubit_to_node,
                                                         nqubits, m, n)

    def node_mapping_with_forced_pairs(self, nqubits, left):
        # forces the first layer cnots to be together
        # laziness means first layer is determined by highest right-layer
        # generate a new mapping between qubits and connectivity nodes
        node = 0
        node_to_qubits = {}
        qubit_to_node = {}
        # NOTE this may seem backwards, but it's due to the LAZY constraint: if we want to
        # apply leftward, we want to "push" the CNOTs right then take the leftmost layer
        # and vice versa
        layer_info = self.container.right_layers_to_cnot if left else self.container.left_layers_to_cnot
        for i in layer_info[self.container.last_layer]:
            q1, q2 = self.container.cnots[i]
            if q1 in qubit_to_node or q2 in qubit_to_node:
                continue
            else:
                node_to_qubits[node] = {q1, q2}
                qubit_to_node[q1] = node
                qubit_to_node[q2] = node
                node += 1
        for q in range(nqubits):
            if q not in qubit_to_node:
                node_to_qubits[node] = {q}
                qubit_to_node[q] = node
                node += 1

        return node_to_qubits, qubit_to_node

    def generate_connectivity_with_forced_pairs(self, m, n, left, prev_perm, discount, prev_perm_weight, node_to_qubits, qubit_to_node):
        # based on naive exponential decay
        # will be setting node_to_qubits, qubit_to_node = node_mapping_with_forced_pairs(cnots, m * n)
        nnodes = len(node_to_qubits)
        G = nx.Graph()
        if prev_perm is not None:
            prev_qubit_G = graph_from_rev_perm(prev_perm, m, n, prev_perm_weight)
            for q1, q2 in prev_qubit_G.edges:
                n1, n2 = qubit_to_node[q1], qubit_to_node[q2]
                if n1 == n2:
                    continue
                G.add_edge(n1, n2, weight=prev_perm_weight)
        for q in range(nnodes):
            G.add_node(q)

        layer_info = self.container.left_layers_to_cnot if left else self.container.right_layers_to_cnot
        lazy_layer = 'right_layer' if left else 'left_layer'
        for j in range(self.nlayers):
            if j not in layer_info:
                break
            for cnot in layer_info[j]:
                q1, q2 = self.container.cnots[cnot]
                n1, n2 = qubit_to_node[q1], qubit_to_node[q2]
                if n1 == n2:
                    continue
                weight = discount ** (self.container.last_layer - self.container.dependencies[cnot][lazy_layer])
                increment_edge_weight(G, n1, n2, weight)
        return G


class LazyLayerLimitedMapper(Mapper):
    def __init__(self, container, nlayers):
        Mapper.__init__(self, container)
        self.nlayers = nlayers

    def mapper(self, m, n, nqubits, left, prev_perm=None, discount=0.5, prev_perm_weight=0.5):
        G = nx.Graph() if prev_perm is None else graph_from_rev_perm(prev_perm, m, n, prev_perm_weight)
        layer_info = self.container.right_layers_to_cnot if left else self.container.left_layers_to_cnot
        for q in range(nqubits):
            G.add_node(q)
        for j in range(self.nlayers):
            i = self.container.last_layer - j
            if i not in layer_info:
                break
            for cnot in layer_info[i]:
                weight = discount ** j
                q1, q2 = self.container.cnots[cnot]
                increment_edge_weight(G, q1, q2, weight)

        prev_mapping = None
        if prev_perm is not None:
            prev_mapping = [None] * nqubits
            for key in prev_perm:
                prev_mapping[prev_perm[key]] = key

        return labeling_from_graph_laplacian(G, m, n, nqubits, prev_mapping)


class LazyForcedPairMapper(Mapper):
    def __init__(self, container, nlayers):
        Mapper.__init__(self, container)
        self.nlayers = nlayers

    def mapper(self, m, n, nqubits, left, prev_perm=None, discount=0.5, prev_perm_weight=0.5):
        node_to_qubits, qubit_to_node = self.node_mapping_with_forced_pairs(nqubits, left)
        G = self.generate_connectivity_with_forced_pairs(m, n, left, prev_perm, discount,
                                                         prev_perm_weight, node_to_qubits, qubit_to_node)
        coords = node_coordinates_from_laplacian(G, m, n)

        prev_mapping = None
        if prev_perm is not None:
            prev_mapping = [None] * nqubits
            for key in prev_perm:
                prev_mapping[prev_perm[key]] = key

        return get_best_orientation(prev_mapping,
                                    self.mapping_generator(coords, node_to_qubits, qubit_to_node, nqubits, m, n))

    def mapping_generator(self, coords, node_to_qubits, qubit_to_node, nqubits, m, n):
        yield lazy_mapping_from_node_coordinates(coords, node_to_qubits, qubit_to_node, nqubits, m, n)
        yield lazy_mapping_from_node_coordinates([(-c[0], c[1]) for c in coords], node_to_qubits, qubit_to_node, nqubits, m, n)
        if n > 1:
            yield lazy_mapping_from_node_coordinates([(c[0], -c[1]) for c in coords], node_to_qubits, qubit_to_node,
                                                     nqubits, m, n)
            yield lazy_mapping_from_node_coordinates([(-c[0], -c[1]) for c in coords], node_to_qubits, qubit_to_node,
                                                     nqubits, m, n)
            if m == n:
                yield lazy_mapping_from_node_coordinates([(c[1], c[0]) for c in coords], node_to_qubits, qubit_to_node, nqubits, m, n)
                yield lazy_mapping_from_node_coordinates([(-c[1], c[0]) for c in coords], node_to_qubits, qubit_to_node,
                                                         nqubits, m, n)
                yield lazy_mapping_from_node_coordinates([(c[1], -c[0]) for c in coords], node_to_qubits, qubit_to_node,
                                                         nqubits, m, n)
                yield lazy_mapping_from_node_coordinates([(-c[1], -c[0]) for c in coords], node_to_qubits, qubit_to_node,
                                                         nqubits, m, n)

    def node_mapping_with_forced_pairs(self, nqubits, left):
        # forces the first layer cnots to be together
        # laziness means first layer is determined by highest right-layer
        # generate a new mapping between qubits and connectivity nodes
        node = 0
        node_to_qubits = {}
        qubit_to_node = {}
        # NOTE this may seem backwards, but it's due to the LAZY constraint: if we want to
        # apply leftward, we want to "push" the CNOTs right then take the leftmost layer
        # and vice versa
        layer_info = self.container.right_layers_to_cnot if left else self.container.left_layers_to_cnot
        for i in layer_info[self.container.last_layer]:
            q1, q2 = self.container.cnots[i]
            if q1 in qubit_to_node or q2 in qubit_to_node:
                continue
            else:
                node_to_qubits[node] = {q1, q2}
                qubit_to_node[q1] = node
                qubit_to_node[q2] = node
                node += 1
        for q in range(nqubits):
            if q not in qubit_to_node:
                node_to_qubits[node] = {q}
                qubit_to_node[q] = node
                node += 1

        return node_to_qubits, qubit_to_node

    def generate_connectivity_with_forced_pairs(self, m, n, left, prev_perm, discount, prev_perm_weight, node_to_qubits, qubit_to_node):
        # based on naive exponential decay
        # will be setting node_to_qubits, qubit_to_node = node_mapping_with_forced_pairs(cnots, m * n)
        nnodes = len(node_to_qubits)
        G = nx.Graph()
        if prev_perm is not None:
            prev_qubit_G = graph_from_rev_perm(prev_perm, m, n, prev_perm_weight)
            for q1, q2 in prev_qubit_G.edges:
                n1, n2 = qubit_to_node[q1], qubit_to_node[q2]
                if n1 == n2:
                    continue
                G.add_edge(n1, n2, weight=prev_perm_weight)
        for q in range(nnodes):
            G.add_node(q)

        layer_info = self.container.right_layers_to_cnot if left else self.container.left_layers_to_cnot

        for j in range(self.nlayers):
            i = self.container.last_layer - j
            if i not in layer_info:
                break
            for cnot in layer_info[i]:
                q1, q2 = self.container.cnots[cnot]
                n1, n2 = qubit_to_node[q1], qubit_to_node[q2]
                if n1 == n2:
                    continue
                weight = discount ** j
                increment_edge_weight(G, n1, n2, weight)
        return G


def generate_permutations(cnots, nqubits, m, n, discount=0.5, prev_perm_weight=0.5, forward_only=False, fallback_only=False):
    # plan
    # do both forward and backward, meet in the middle (idea is that forward/backward symmetry, as well as that
    # the first/last permutations should insist on something from first layer being applied)
    # for each gate in the current direction, assign with a weight equal to discount^layer, where
    # the layer number for a cnot(q1, q2) is the max(#cnots ahead with q1, #cnots ahead with q2)
    # finally, add as a "constant layer" a prev_perm_weight edge-weight grid that gives some memory
    # of the previous permutation.
    # for each permutation generated, greedily apply as many gates as possible.
    container = CNOTContainer(cnots)
    ncnots = len(cnots)
    # mapper = LazyLayerLimitedMapper(container, nqubits)
    # fallback_mapper = LazyForcedPairMapper(container, 4 * nqubits)

    mapper = GreedyLayerLimitedMapper(container, nqubits)
    fallback_mapper = GreedyForcedPairMapper(container, 4 * nqubits)

    forward_prev_perm = None
    backward_prev_perm = None

    forward_permutations = []
    forward_index_layer = []

    backward_permutations = []
    backward_index_layer = []

    while len(container.applied_indices) < ncnots:

        applied_before_generation = len(container.applied_indices)

        if not fallback_only:
            forward_mappings, forward_rev_mapping = mapper.mapper(m, n, nqubits, True, forward_prev_perm, discount, prev_perm_weight)
            if not forward_only:
                backward_mappings, backward_rev_mapping = mapper.mapper(m, n, nqubits, False, backward_prev_perm, discount, prev_perm_weight)

            forward_applied_indices = container.apply_permutation(forward_mappings, True)
            if not forward_only:
                backward_applied_indices = container.apply_permutation(backward_mappings, False)

        applied_after_generation = len(container.applied_indices)
        if applied_after_generation == applied_before_generation:
            # fallback, currently forcing first layer pairs together
            forward_mappings, forward_rev_mapping = fallback_mapper.mapper(m, n, nqubits, True, forward_prev_perm, discount, prev_perm_weight)
            if not forward_only:
                backward_mappings, backward_rev_mapping = fallback_mapper.mapper(m, n, nqubits, False, backward_prev_perm, discount, prev_perm_weight)
            forward_applied_indices = container.apply_permutation(forward_mappings, True)
            if not forward_only:
                backward_applied_indices = container.apply_permutation(backward_mappings, False)


            applied_after_fallback = len(container.applied_indices)

            if applied_after_fallback == applied_before_generation:
                raise RuntimeError("Fails even with fallback")

        if len(forward_applied_indices) > 0:
            forward_permutations.append(forward_mappings)
            forward_index_layer.append(forward_applied_indices)

        forward_prev_perm = forward_rev_mapping

        if not forward_only:
            if len(backward_applied_indices) > 0:
                backward_permutations.append(backward_mappings)
                backward_index_layer.append(backward_applied_indices)

            backward_prev_perm = backward_rev_mapping

    return forward_permutations, forward_index_layer, \
           backward_permutations, backward_index_layer


def generate_compliant_circuit(interactions, gate_blocks, nqubits, m, n, discount=0.5, prev_perm_weight=0.5, forward_only=False, fallback_only=False):
    # physical location label (i,j) --> i + m * j
    forward_permutations, forward_index_layer, \
    backward_permutations, backward_index_layer = generate_permutations(interactions, nqubits, m, n, discount, prev_perm_weight, forward_only, fallback_only)

    permutations = forward_permutations + backward_permutations[::-1]
    index_layers = forward_index_layer + backward_index_layer[::-1]

    gates = []

    start_perm = permutations[0]
    current_perm = permutations[0]
    i = 0
    while True:
        # apply current CNOTs to the current permutation
        indices_to_apply = index_layers[i]
        for j in indices_to_apply:
            for gate, *qs in gate_blocks[j]:
                locs = [current_perm[q] for q in qs]
                physes = [location_to_label(*loc, m, n) for loc in locs]
                gates.append((gate, *physes))
        i += 1
        if i >= len(permutations):
            break
        # swap to next permutation
        old_perm = current_perm
        current_perm = permutations[i]
        perm1 = {key: old_perm[key] for key in range(len(old_perm))}
        perm2 = {key: current_perm[key] for key in range(len(current_perm))}

        # Approximate
        for loc1, loc2 in swap_between_perm(perm1, perm2):
            phys1, phys2 = location_to_label(*loc1, m, n), location_to_label(*loc2, m, n)
            gates.append(('swap', phys1, phys2))

        # # Brute force optimal
        # for phys1, phys2 in brute_force_permutation(perm1, perm2, m, n):
        #     gates.append(('swap', phys1, phys2))

    return start_perm, gates