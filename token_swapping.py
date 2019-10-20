import random
import numpy as np

# based on https://arxiv.org/pdf/1602.05150.pdf


def assert_consistent(perm, rev_perm):
    for k in perm:
        if rev_perm[perm[k]] != k:
            return False
    return True


def swap_between_perm(perm1, perm2):
    perm1 = {key: perm1[key] for key in perm1}  # copy perm1
    rev_perm1 = {perm1[key]: key for key in perm1}  # reverse perm1

    out_of_position = {key for key in perm1 if perm1[key] != perm2[key]}
    while len(out_of_position) > 0:
        current_token = random.choice(tuple(out_of_position))  # pick a random starting point
        current_location = perm1[current_token]
        visited = set()
        path = []
        happy = False
        while True:
            path.append(current_token)
            if current_location in visited:
                # happy swap chain detected
                happy = True
                break
            visited.add(current_location)
            desired_location = perm2[current_token]
            if desired_location == current_location:
                # unhappy swap detected
                break

            dx = np.sign(desired_location[0] - current_location[0])
            dy = np.sign(desired_location[1] - current_location[1])
            if dx == 0:
                go_x_dir = False
            elif dy == 0:
                go_x_dir = True
            else:
                # randomly choose which way to go
                go_x_dir = bool(random.getrandbits(1))
            if go_x_dir:
                current_location = (current_location[0] + dx, current_location[1])
            else:
                current_location = (current_location[0], current_location[1] + dy)
            current_token = rev_perm1[current_location]

        if happy:
            end_token = path[-1]  # where the cycle ended
            token1 = path[-2]  # this token will be swapped all the way around
            i = -3
            while True:
                loc1 = perm1[token1]
                token2 = path[i]
                loc2 = perm1[token2]
                yield (loc1, loc2)
                perm1[token1], perm1[token2] = perm1[token2], perm1[token1]  # tokens swap locations
                rev_perm1[loc1], rev_perm1[loc2] = rev_perm1[loc2], rev_perm1[loc1]  # locations swap tokens
                if perm1[token2] == perm2[token2]:
                    out_of_position.remove(token2)
                if token2 == end_token:
                    break
                i -= 1
            if perm1[token1] == perm2[token1]:
                out_of_position.remove(token1)
        else:
            token1 = path[-1]
            token2 = path[-2]
            loc1 = perm1[token1]
            loc2 = perm1[token2]
            yield (loc1, loc2)
            # this unhappy swap moves token1 out of position, but clearly does not move token2 into position
            out_of_position.add(token1)
            perm1[token1], perm1[token2] = perm1[token2], perm1[token1]
            rev_perm1[loc1], rev_perm1[loc2] = rev_perm1[loc2], rev_perm1[loc1]


def random_perm(m, n):
    perm = {}
    order = np.random.permutation(range(m * n))
    for idx, label in enumerate(order):
        perm[label] = (idx//m, idx % m)
    return perm


def brute_force_permutation(initial_perm, final_perm, m, n):
    # naive BFS search
    # perms given as logical --> physical location

    # want initial, final to be physical label --> logical label
    initial = [None] * m * n
    final = [None] * m * n
    for logical in initial_perm:
        loc = initial_perm[logical]
        physical = loc[0] + loc[1] * m
        initial[physical] = logical

    for logical in final_perm:
        loc = final_perm[logical]
        physical = loc[0] + loc[1] * m
        final[physical] = logical

    initial = tuple(initial)
    final = tuple(final)

    path = {final: {"parent": None}}
    queue = [final]
    idx = 0
    while idx < len(queue):
        current = queue[idx]
        for i in range(m - 1):
            for j in range(n):
                adj = list(current)
                loc1 = i + m * j
                loc2 = (i + 1) + m * j
                adj[loc1], adj[loc2] = adj[loc2], adj[loc1]
                adj = tuple(adj)
                if adj not in path:
                    path[adj] = {"parent": current, "swap": (loc1, loc2)}
                    if adj == initial:
                        break
                    queue.append(adj)
        idx += 1
    final_path = []
    current = initial
    while current != final:
        final_path.append(path[current]['swap'])
        current = path[current]['parent']
    return final_path