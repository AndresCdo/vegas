#!/usr/bin/env python3
"""
Generate a small triangular lattice for testing graph coloring.
Triangular lattice requires 3 colors (non-bipartite).
Creates a 3x3 lattice with periodic boundary conditions.
"""

import os


def write_triangular_lattice_3x3(output_path):
    """
    Generate a 3x3 triangular lattice (9 atoms) with periodic boundaries.
    Each atom has 6 neighbors (hexagonal).
    Positions: (x, y) with basis vectors:
        a1 = (1, 0)
        a2 = (0.5, sqrt(3)/2) ≈ (0.5, 0.866)
    """
    L = 3  # 3x3 unit cells
    atoms = []
    interactions = []

    # Generate atoms
    idx = 0
    for ix in range(L):
        for iy in range(L):
            # Position in triangular lattice coordinates
            x = ix + 0.5 * iy
            y = iy * 0.8660254  # sqrt(3)/2

            atoms.append(
                {
                    "idx": idx,
                    "px": float(x),
                    "py": float(y),
                    "pz": 0.0,
                    "spin_norm": 1.0,
                    "hx": 0.0,
                    "hy": 0.0,
                    "hz": 1.0,
                    "type_name": "Fe",
                    "model_name": "heisenberg",  # model doesn't matter for coloring
                }
            )
            idx += 1

    # Neighbor vectors for triangular lattice (6 neighbors)
    # Using integer offsets in the two basis directions
    neighbors = [
        (1, 0),  # right
        (0, 1),  # up-right
        (-1, 1),  # up-left
        (-1, 0),  # left
        (0, -1),  # down-left
        (1, -1),  # down-right
    ]

    J = 1.0  # exchange constant

    # Create interactions with periodic boundary conditions
    for i in range(L * L):
        ix1 = i // L
        iy1 = i % L
        for dx, dy in neighbors:
            ix2 = (ix1 + dx) % L
            iy2 = (iy1 + dy) % L
            j = ix2 * L + iy2

            # Add bidirectional bond (i -> j and j -> i)
            # Avoid duplicates by only adding if i < j
            if i < j:
                interactions.append((i, j, J))
                interactions.append((j, i, J))

    # Write lattice file
    n_atoms = len(atoms)
    n_inter = len(interactions)
    n_types = 1

    with open(output_path, "w") as f:
        f.write(f"{n_atoms} {n_inter} {n_types}\n")
        f.write("Fe\n")
        for a in atoms:
            f.write(
                f"{a['idx']} "
                f"{a['px']:.6f} {a['py']:.6f} {a['pz']:.6f} "
                f"{a['spin_norm']:.6f} "
                f"{a['hx']:.6f} {a['hy']:.6f} {a['hz']:.6f} "
                f"{a['type_name']} {a['model_name']}\n"
            )
        for i, j, Jval in interactions:
            f.write(f"{i} {j} {Jval:.6f}\n")

    print(f"Generated triangular lattice: {output_path}")
    print(f"  Atoms: {n_atoms}, Interactions: {n_inter}")
    print(f"  Expected colors: 3 (triangular lattice is non-bipartite)")


if __name__ == "__main__":
    output_path = os.path.join(os.path.dirname(__file__), "triangular_3x3_pbc.txt")
    write_triangular_lattice_3x3(output_path)
