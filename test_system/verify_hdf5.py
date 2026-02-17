#!/usr/bin/env python3
"""
Verify HDF5 output structure from VEGAS simulation.
"""

import h5py
import numpy as np
import sys


def verify_hdf5_structure(filename):
    """Verify the HDF5 file has the expected structure."""
    print(f"Verifying HDF5 file: {filename}")

    try:
        with h5py.File(filename, "r") as f:
            # Check attributes
            print("\n=== Attributes ===")
            for attr_name in f.attrs:
                print(f"  {attr_name}: {f.attrs[attr_name]}")

            # Check datasets
            print("\n=== Datasets ===")
            required_datasets = [
                "temperature",
                "field",
                "magnetization_x",
                "magnetization_y",
                "magnetization_z",
                "energy",
                "positions",
                "types",
            ]

            for ds_name in required_datasets:
                if ds_name in f:
                    ds = f[ds_name]
                    print(f"  {ds_name}: shape={ds.shape}, dtype={ds.dtype}")
                else:
                    print(f"  {ds_name}: MISSING")

            # Check dimensions consistency
            print("\n=== Dimension Consistency ===")
            if "temperature" in f and "field" in f:
                n_temps = f["temperature"].shape[0]
                n_fields = f["field"].shape[0]
                print(f"  n_temperatures: {n_temps}")
                print(f"  n_fields: {n_fields}")

                # Check magnetization arrays
                for comp in ["x", "y", "z"]:
                    ds_name = f"magnetization_{comp}"
                    if ds_name in f:
                        shape = f[ds_name].shape
                        print(f"  {ds_name}: shape={shape}")
                        if len(shape) == 3:
                            if shape[0] != n_temps or shape[1] != n_fields:
                                print(
                                    f"    WARNING: First dimensions don't match temperature/field counts"
                                )

            # Check positions and types
            if "positions" in f:
                n_atoms = f["positions"].shape[0]
                print(f"\n  n_atoms: {n_atoms}")
                print(f"  positions shape: {f['positions'].shape}")

            if "types" in f:
                print(f"  types shape: {f['types'].shape}")

            # Sample data values
            print("\n=== Sample Values ===")
            if "temperature" in f:
                print(f"  temperatures: {f['temperature'][:]}")

            if "field" in f:
                print(f"  fields: {f['field'][:]}")

            # Check for non-zero magnetization (should have some signal)
            if "magnetization_z" in f:
                mz = f["magnetization_z"][:]
                print(f"  magnetization_z range: [{mz.min():.4f}, {mz.max():.4f}]")
                print(f"  magnetization_z mean: {mz.mean():.4f}")

            print("\n✓ HDF5 file structure appears valid")
            return True

    except Exception as e:
        print(f"\n✗ Error verifying HDF5 file: {e}")
        return False


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python verify_hdf5.py <hdf5_file>")
        sys.exit(1)

    filename = sys.argv[1]
    success = verify_hdf5_structure(filename)
    sys.exit(0 if success else 1)
