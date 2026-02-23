#!/usr/bin/env python3
"""Extract observables from VEGAS HDF5 output files."""

import sys
import h5py
import numpy as np
import glob


def extract_observable(h5_file, observable):
    """Extract a scalar observable from HDF5 output."""

    with h5py.File(h5_file, "r") as f:
        temps = f["temperature"][:]
        kb = f.attrs["kb"]

        if observable == "mean_energy_per_atom":
            energy = f["energy"][:]
            n_atoms = f["positions"].shape[0]
            return np.mean(energy[0]) / n_atoms

        elif observable == "energy_at_temp":
            temp_idx = int(sys.argv[3]) if len(sys.argv) > 3 else 0
            energy = f["energy"][:]
            n_atoms = f["positions"].shape[0]
            return np.mean(energy[temp_idx]) / n_atoms

        elif observable == "magnetization_at_low_T":
            mx = f["magnetization_x"][:]
            my = f["magnetization_y"][:]
            mz = f["magnetization_z"][:]
            n_atoms = f["positions"].shape[0]
            M = (
                np.sqrt(
                    np.mean(mx, axis=1) ** 2
                    + np.mean(my, axis=1) ** 2
                    + np.mean(mz, axis=1) ** 2
                )
                / n_atoms
            )
            return M[0]

        elif observable == "magnetization_at_high_T":
            mx = f["magnetization_x"][:]
            my = f["magnetization_y"][:]
            mz = f["magnetization_z"][:]
            n_atoms = f["positions"].shape[0]
            M = (
                np.sqrt(
                    np.mean(mx, axis=1) ** 2
                    + np.mean(my, axis=1) ** 2
                    + np.mean(mz, axis=1) ** 2
                )
                / n_atoms
            )
            return M[-1]

        elif observable == "magnetization_all":
            mx = f["magnetization_x"][:]
            my = f["magnetization_y"][:]
            mz = f["magnetization_z"][:]
            n_atoms = f["positions"].shape[0]
            M = (
                np.sqrt(
                    np.mean(mx, axis=1) ** 2
                    + np.mean(my, axis=1) ** 2
                    + np.mean(mz, axis=1) ** 2
                )
                / n_atoms
            )
            return M

        elif observable == "specific_heat_peak_location":
            energy = f["energy"][:]
            n_atoms = f["positions"].shape[0]
            cv = np.var(energy, axis=1) / (kb * temps**2) / n_atoms
            return temps[np.argmax(cv)]

        elif observable == "susceptibility_peak_location":
            mx = f["magnetization_x"][:]
            n_atoms = f["positions"].shape[0]
            chi = (
                (np.mean(mx**2, axis=1) - np.mean(mx, axis=1) ** 2)
                / (kb * temps)
                / n_atoms
            )
            return temps[np.argmax(chi)]

        elif observable == "susceptibility":
            mx = f["magnetization_x"][:]
            n_atoms = f["positions"].shape[0]
            chi = (
                (np.mean(mx**2, axis=1) - np.mean(mx, axis=1) ** 2)
                / (kb * temps)
                / n_atoms
            )
            return chi

        elif observable == "specific_heat":
            energy = f["energy"][:]
            n_atoms = f["positions"].shape[0]
            cv = np.var(energy, axis=1) / (kb * temps**2) / n_atoms
            return cv

        elif observable == "mean_magnetization_per_atom":
            mx = f["magnetization_x"][:]
            my = f["magnetization_y"][:]
            mz = f["magnetization_z"][:]
            n_atoms = f["positions"].shape[0]
            M = (
                np.sqrt(
                    np.mean(mx, axis=1) ** 2
                    + np.mean(my, axis=1) ** 2
                    + np.mean(mz, axis=1) ** 2
                )
                / n_atoms
            )
            return M

        elif observable == "mean_magnetization_components":
            mx = f["magnetization_x"][:]
            my = f["magnetization_y"][:]
            mz = f["magnetization_z"][:]
            n_atoms = f["positions"].shape[0]
            return (
                np.mean(mx, axis=1) / n_atoms,
                np.mean(my, axis=1) / n_atoms,
                np.mean(mz, axis=1) / n_atoms,
            )

        elif observable == "second_moment_magnetization":
            mx = f["magnetization_x"][:]
            my = f["magnetization_y"][:]
            mz = f["magnetization_z"][:]
            n_atoms = f["positions"].shape[0]
            return (
                np.mean(mx**2, axis=1) / n_atoms,
                np.mean(my**2, axis=1) / n_atoms,
                np.mean(mz**2, axis=1) / n_atoms,
            )

        elif observable == "energy_std":
            energy = f["energy"][:]
            return np.std(energy[0])

        else:
            raise ValueError(f"Unknown observable: {observable}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: extract_observable.py <h5_file> <observable> [args...]")
        print(
            "Observables: mean_energy_per_atom, magnetization_at_low_T, magnetization_at_high_T,"
        )
        print("             specific_heat_peak_location, susceptibility_peak_location,")
        print("             mean_magnetization_per_atom, second_moment_magnetization")
        sys.exit(1)

    h5_file = sys.argv[1]
    observable = sys.argv[2]

    result = extract_observable(h5_file, observable)

    if isinstance(result, tuple):
        print(" ".join(str(x) for x in result))
    else:
        print(result)
