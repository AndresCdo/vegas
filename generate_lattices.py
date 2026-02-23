#!/usr/bin/env python3
"""
VEGAS Benchmark Lattice Generator
Generates all lattice files needed for the physical validation benchmark suite.

Usage:
    python generate_lattices.py [--output-dir benchmarks/lattices]

All lattices use VEGAS format:
    Line 1:  num_atoms num_interactions num_types
    Line 2:  type_name_1 [type_name_2 ...]
    Lines 3..N+2:  idx px py pz spin_norm hx hy hz type_name model_name
    Lines N+3..:   atom_i atom_j exchange_J
"""

import argparse
import math
import os


def write_lattice(path, atoms, interactions, type_names):
    """
    atoms: list of dicts with keys:
        idx, px, py, pz, spin_norm, hx, hy, hz, type_name, model_name
    interactions: list of (i, j, J)
    type_names: ordered list of unique type name strings
    """
    n_atoms = len(atoms)
    n_inter = len(interactions)
    n_types = len(type_names)

    with open(path, "w") as f:
        f.write(f"{n_atoms} {n_inter} {n_types}\n")
        f.write(" ".join(type_names) + "\n")
        for a in atoms:
            f.write(
                f"{a['idx']} "
                f"{a['px']:.6f} {a['py']:.6f} {a['pz']:.6f} "
                f"{a['spin_norm']:.6f} "
                f"{a['hx']:.6f} {a['hy']:.6f} {a['hz']:.6f} "
                f"{a['type_name']} {a['model_name']}\n"
            )
        for i, j, J in interactions:
            f.write(f"{i} {j} {J:.6f}\n")

    print(f"  Written: {path}  ({n_atoms} atoms, {n_inter} interactions)")


def make_ising_chain_100(out_dir):
    """
    1D open chain, 100 atoms, nearest-neighbour exchange J = 1.0.
    Open boundaries (no wrap): 198 interactions (bidirectional).

    Exchange convention (standard physics): E = -J * S_i · S_j
    - J > 0: Ferromagnetic (favors alignment)
    - J < 0: Antiferromagnetic (favors anti-alignment)
    """

    N = 100
    J = 1.0  # Ferromagnetic coupling

    atoms = [
        {
            "idx": i,
            "px": float(i),
            "py": 0.0,
            "pz": 0.0,
            "spin_norm": 1.0,
            "hx": 0.0,
            "hy": 0.0,
            "hz": 1.0,
            "type_name": "Fe",
            "model_name": "ising",
        }
        for i in range(N)
    ]

    # Store each bond in BOTH directions
    interactions = []
    for i in range(N - 1):
        interactions.append((i, i + 1, J))  # i -> i+1
        interactions.append((i + 1, i, J))  # i+1 -> i

    write_lattice(
        os.path.join(out_dir, "ising_chain_100.txt"), atoms, interactions, ["Fe"]
    )


def make_ising_square_20x20_pbc(out_dir):
    """
    20x20 square lattice, J = +1.0 (ferromagnetic in VEGAS convention), with periodic boundary conditions.
    Each atom connects to right, up, left, down neighbours (bidirectional).

    NOTE: For critical temperature test, we need FERROMAGNETIC coupling to observe
    spontaneous magnetization below Tc. VEGAS uses INVERTED convention:
    - J = +1 gives ferromagnetic behavior (standard physics J < 0)

    CRITICAL: Each bond must be stored for BOTH atoms for proper Metropolis dynamics!
    """
    L = 20
    J = 1.0  # Ferromagnetic in VEGAS convention

    def idx(i, j):
        return i * L + j

    atoms = [
        {
            "idx": idx(i, j),
            "px": float(j),
            "py": float(i),
            "pz": 0.0,
            "spin_norm": 1.0,
            "hx": 0.0,
            "hy": 0.0,
            "hz": 1.0,
            "type_name": "Fe",
            "model_name": "ising",
        }
        for i in range(L)
        for j in range(L)
    ]

    interactions = []
    for i in range(L):
        for j in range(L):
            # Store each bond in BOTH directions
            # Right neighbor
            interactions.append((idx(i, j), idx(i, (j + 1) % L), J))
            interactions.append((idx(i, (j + 1) % L), idx(i, j), J))
            # Up neighbor
            interactions.append((idx(i, j), idx((i + 1) % L, j), J))
            interactions.append((idx((i + 1) % L, j), idx(i, j), J))

    write_lattice(
        os.path.join(out_dir, "ising_square_20x20_pbc.txt"), atoms, interactions, ["Fe"]
    )


def make_ferromagnet_square_20x20_pbc(out_dir):
    """
    Same geometry as B2 but with J = +1.0 (ferromagnetic coupling).

    Ground state energy per atom:
        e_0 = -|J| * z / 2 = -1.0 * 4 / 2 = -2.0  (z=4 for square lattice)

    Exchange convention (standard physics): E = -J * S_i · S_j
    - J > 0: Ferromagnetic (favors alignment)
    - J < 0: Antiferromagnetic (favors anti-alignment)

    CRITICAL: Each bond must be stored for BOTH atoms for proper Metropolis dynamics!
    """
    L = 20
    J = 1.0  # Ferromagnetic coupling

    def idx(i, j):
        return i * L + j

    atoms = [
        {
            "idx": idx(i, j),
            "px": float(j),
            "py": float(i),
            "pz": 0.0,
            "spin_norm": 1.0,
            "hx": 0.0,
            "hy": 0.0,
            "hz": 1.0,
            "type_name": "Fe",
            "model_name": "ising",
        }
        for i in range(L)
        for j in range(L)
    ]

    interactions = []
    for i in range(L):
        for j in range(L):
            # Store each bond in BOTH directions
            # Right neighbor
            interactions.append((idx(i, j), idx(i, (j + 1) % L), J))
            interactions.append((idx(i, (j + 1) % L), idx(i, j), J))
            # Up neighbor
            interactions.append((idx(i, j), idx((i + 1) % L, j), J))
            interactions.append((idx((i + 1) % L, j), idx(i, j), J))

    write_lattice(
        os.path.join(out_dir, "ferromagnet_square_20x20_pbc.txt"),
        atoms,
        interactions,
        ["Fe"],
    )


def make_initial_states(out_dir, n_atoms=400):
    """
    Initial state files for ergodicity test (B4).
    """
    up_path = os.path.join(out_dir, "all_spins_up.txt")
    down_path = os.path.join(out_dir, "all_spins_down.txt")

    with open(up_path, "w") as f:
        for _ in range(n_atoms):
            f.write("0.0 0.0 1.0\n")
    print(f"  Written: {up_path}  ({n_atoms} spins pointing +z)")

    with open(down_path, "w") as f:
        for _ in range(n_atoms):
            f.write("0.0 0.0 -1.0\n")
    print(f"  Written: {down_path}  ({n_atoms} spins pointing -z)")


def make_heisenberg_chain_100(out_dir):
    """
    1D Heisenberg chain for isotropy test at T -> infinity (B6).

    At T = 1000 >> J = 1, the spin distribution must be isotropic:
        <Sx> = <Sy> = <Sz> = 0
        <Sx²> = <Sy²> = <Sz²> = 1/3

    Uses 'heisenberg' model so the v2.3.0 Gaussian perturbation proposal is exercised.

    CRITICAL: Each bond must be stored for BOTH atoms for proper Metropolis dynamics!
    """
    N = 100
    J = 1.0

    atoms = [
        {
            "idx": i,
            "px": float(i),
            "py": 0.0,
            "pz": 0.0,
            "spin_norm": 1.0,
            "hx": 0.0,
            "hy": 0.0,
            "hz": 1.0,
            "type_name": "Fe",
            "model_name": "heisenberg",
        }
        for i in range(N)
    ]

    # Store each bond in BOTH directions
    interactions = []
    for i in range(N - 1):
        interactions.append((i, i + 1, J))
        interactions.append((i + 1, i, J))

    write_lattice(
        os.path.join(out_dir, "heisenberg_chain_100.txt"), atoms, interactions, ["Fe"]
    )


def make_small_ising_for_sigma(out_dir):
    """
    Small 5x5 Ising lattice for B5 (sigma freeze test).
    Small size means each benchmark run is fast, suitable for repeated testing.

    CRITICAL: Each bond must be stored for BOTH atoms for proper Metropolis dynamics!
    """
    L = 5
    J = 1.0

    def idx(i, j):
        return i * L + j

    atoms = [
        {
            "idx": idx(i, j),
            "px": float(j),
            "py": float(i),
            "pz": 0.0,
            "spin_norm": 1.0,
            "hx": 0.0,
            "hy": 0.0,
            "hz": 1.0,
            "type_name": "Fe",
            "model_name": "ising",
        }
        for i in range(L)
        for j in range(L)
    ]

    interactions = []
    for i in range(L):
        for j in range(L):
            # Store each bond in BOTH directions
            # Right neighbor
            interactions.append((idx(i, j), idx(i, (j + 1) % L), J))
            interactions.append((idx(i, (j + 1) % L), idx(i, j), J))
            # Up neighbor
            interactions.append((idx(i, j), idx((i + 1) % L, j), J))
            interactions.append((idx((i + 1) % L, j), idx(i, j), J))

    write_lattice(
        os.path.join(out_dir, "ising_small_5x5_pbc.txt"), atoms, interactions, ["Fe"]
    )


def verify_file(path):
    with open(path) as f:
        lines = f.readlines()

    header = lines[0].split()
    n_atoms = int(header[0])
    n_inter = int(header[1])
    n_types = int(header[2])

    atom_lines = lines[2 : 2 + n_atoms]
    inter_lines = lines[2 + n_atoms : 2 + n_atoms + n_inter]

    J_values = [float(l.split()[2]) for l in inter_lines]
    J_unique = set(J_values)

    models = set(l.split()[9] for l in atom_lines)
    types = set(l.split()[8] for l in atom_lines)

    coord = {}
    for l in inter_lines:
        a, b = int(l.split()[0]), int(l.split()[1])
        coord[a] = coord.get(a, 0) + 1
        coord[b] = coord.get(b, 0) + 1
    coord_values = list(coord.values())
    min_z = min(coord_values)
    max_z = max(coord_values)
    mean_z = sum(coord_values) / len(coord_values)

    print(f"  {os.path.basename(path)}:")
    print(f"    atoms={n_atoms}, interactions={n_inter}, types={n_types}")
    print(f"    models={models}, atom_types={types}")
    print(f"    J_values={J_unique}")
    print(f"    coordination: min={min_z}, max={max_z}, mean={mean_z:.2f}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate VEGAS benchmark lattice files"
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        default="benchmarks/lattices",
        help="Output directory for lattice files (default: benchmarks/lattices)",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Print summary statistics for each generated file",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nGenerating lattice files in: {args.output_dir}\n")

    generators = [
        ("B1 - 1D Ising chain (100 atoms, open, J=+1)", make_ising_chain_100),
        ("B2 - 2D Ising square 20x20 (PBC, J=+1)", make_ising_square_20x20_pbc),
        (
            "B3 - 2D Ferromagnet square 20x20 (PBC, J=-1)",
            make_ferromagnet_square_20x20_pbc,
        ),
        (
            "B4 - Initial state files (up / down, 400 spins)",
            lambda d: make_initial_states(d, n_atoms=400),
        ),
        ("B5 - Small Ising 5x5 (PBC, J=+1, sigma test)", make_small_ising_for_sigma),
        ("B6 - 1D Heisenberg chain (100 atoms, open, J=+1)", make_heisenberg_chain_100),
    ]

    for label, fn in generators:
        print(f"{label}")
        fn(args.output_dir)
        print()

    if args.verify:
        print("\n--- Verification ---\n")
        lattice_files = [
            "ising_chain_100.txt",
            "ising_square_20x20_pbc.txt",
            "ferromagnet_square_20x20_pbc.txt",
            "ising_small_5x5_pbc.txt",
            "heisenberg_chain_100.txt",
        ]
        for fname in lattice_files:
            path = os.path.join(args.output_dir, fname)
            if os.path.exists(path):
                verify_file(path)
            else:
                print(f"  MISSING: {path}")

    print("\nDone. Lattice files generated successfully.")
    print(
        "B3 ground state energy: e_expected = -2.0 per atom (0.5x correction applied in totalEnergy)"
    )


if __name__ == "__main__":
    main()
