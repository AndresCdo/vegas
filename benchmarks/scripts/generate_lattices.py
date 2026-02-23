#!/usr/bin/env python3
"""Generate lattice files for VEGAS benchmarks."""

import math


def generate_ising_chain_100():
    """Generate 1D Ising chain with 100 atoms, periodic boundary."""
    n_atoms = 100

    # Header: num_ions num_interactions num_types
    # Bidirectional bonds: each of 100 bonds stored as i→j AND j→i
    num_interactions = 200  # 100 bonds * 2 directions

    lines = []
    lines.append(f"{n_atoms} {num_interactions} 1")
    lines.append("Fe")  # Atom type

    # Atoms: index px py pz spinNorm hx hy hz type model
    # Position on a line, spin normalization 1.0, no external field
    for i in range(n_atoms):
        lines.append(f"{i} {i} 0 0 1.0 0 0 0 Fe ising")

    # Interactions: index neighbor exchange
    # J = 1.0 for Ising chain (ferromagnetic in VEGAS convention)
    # Bidirectional: i→i+1 and i+1→i
    J = 1.0
    for i in range(n_atoms):
        j = (i + 1) % n_atoms  # Periodic
        lines.append(f"{i} {j} {J}")  # i → j
        lines.append(f"{j} {i} {J}")  # j → i

    with open("benchmarks/lattices/ising_chain_100.txt", "w") as f:
        f.write("\n".join(lines))
    print(
        f"Created ising_chain_100.txt: {n_atoms} atoms, {num_interactions} interactions"
    )


def generate_ising_square_20x20():
    """Generate 2D square lattice 20x20 with periodic boundary conditions."""
    L = 20
    n_atoms = L * L

    # Each interior atom has 4 neighbors, edges have fewer with PBC
    # Total bonds: 2 * L * L (horizontal + vertical), each stored bidirectionally
    num_interactions = 2 * n_atoms * 2  # 2 directions per bond

    lines = []
    lines.append(f"{n_atoms} {num_interactions} 1")
    lines.append("Fe")

    # Atoms: index = y * L + x
    for y in range(L):
        for x in range(L):
            idx = y * L + x
            lines.append(f"{idx} {x} {y} 0 1.0 0 0 0 Fe ising")

    # Interactions: periodic boundary conditions
    J = 1.0
    for y in range(L):
        for x in range(L):
            idx = y * L + x

            # Right neighbor
            xr = (x + 1) % L
            idxr = y * L + xr
            lines.append(f"{idx} {idxr} {J}")
            lines.append(f"{idxr} {idx} {J}")

            # Up neighbor
            yu = (y + 1) % L
            idxu = yu * L + x
            lines.append(f"{idx} {idxu} {J}")
            lines.append(f"{idxu} {idx} {J}")

    with open("benchmarks/lattices/ising_square_20x20.txt", "w") as f:
        f.write("\n".join(lines))
    print(
        f"Created ising_square_20x20.txt: {n_atoms} atoms, {num_interactions} interactions"
    )


def generate_simple_10x10():
    """Generate simple 10x10 lattice for B5 and B6."""
    L = 10
    n_atoms = L * L
    num_interactions = 2 * n_atoms * 2  # Bidirectional

    lines = []
    lines.append(f"{n_atoms} {num_interactions} 1")
    lines.append("Fe")

    for y in range(L):
        for x in range(L):
            idx = y * L + x
            lines.append(f"{idx} {x} {y} 0 1.0 0 0 0 Fe heisenberg")

    J = 1.0
    for y in range(L):
        for x in range(L):
            idx = y * L + x

            xr = (x + 1) % L
            idxr = y * L + xr
            lines.append(f"{idx} {idxr} {J}")
            lines.append(f"{idxr} {idx} {J}")

            yu = (y + 1) % L
            idxu = yu * L + x
            lines.append(f"{idx} {idxu} {J}")
            lines.append(f"{idxu} {idx} {J}")

    with open("benchmarks/lattices/simple_10x10.txt", "w") as f:
        f.write("\n".join(lines))
    print(f"Created simple_10x10.txt: {n_atoms} atoms, {num_interactions} interactions")


def generate_initial_states():
    """Generate initial state files for ergodicity test."""
    # All spins up (+z direction)
    lines_up = []
    for i in range(400):  # Match 20x20 lattice
        lines_up.append("0 0 1")
    with open("benchmarks/lattices/all_spins_up.txt", "w") as f:
        f.write("\n".join(lines_up))
    print("Created all_spins_up.txt")

    # All spins down (-z direction)
    lines_down = []
    for i in range(400):
        lines_down.append("0 0 -1")
    with open("benchmarks/lattices/all_spins_down.txt", "w") as f:
        f.write("\n".join(lines_down))
    print("Created all_spins_down.txt")


if __name__ == "__main__":
    generate_ising_chain_100()
    generate_ising_square_20x20()
    generate_simple_10x10()
    generate_initial_states()
    print("All lattice files created!")
