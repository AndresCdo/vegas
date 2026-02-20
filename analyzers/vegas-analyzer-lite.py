"""
VEGAS Lite Analyzer - Simple magnetization analysis for HDF5 output files.
"""

from typing import Optional
import h5py
import matplotlib

matplotlib.use("pdf")
from matplotlib import pyplot
import click
import numpy as np


@click.command()
@click.argument("file", type=click.Path(exists=True))
def main(file: str) -> None:
    """Analyze VEGAS simulation output and generate plots."""
    with h5py.File(file, mode="r") as dataset:
        mcs: int = dataset.attrs["mcs"]
        seed: int = dataset.attrs["seed"]
        tau: int = mcs // 5
        num_ions: int = len(dataset.get("positions"))
        temps: np.ndarray = dataset.get("temperature")[:]
        fields: np.ndarray = dataset.get("field")[:]
        Mz: np.ndarray = dataset.get("magnetization_z")[:, tau:] / num_ions
        Mz_mean: np.ndarray = np.mean(Mz, axis=1)
        zeros: np.ndarray = np.zeros_like(Mz_mean)

    with open(file.replace(".h5", ".mean"), mode="w") as file_:
        file_.write("# seed = {}\n".format(seed))
        file_.write("#\tT\tH\tE\tCv\tM\tMz\tX\n")
        for i, T in enumerate(temps):
            file_.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    T, fields[i], zeros[i], zeros[i], zeros[i], Mz_mean[i], zeros[i]
                )
            )

    pyplot.figure()
    pyplot.plot(fields, Mz_mean, "-o", label="$M_{z}$")
    pyplot.grid()
    pyplot.xlabel("$H$", fontsize=20)
    pyplot.ylabel("$M_{z}$", fontsize=20)
    pyplot.tight_layout()
    pyplot.savefig(file.replace(".h5", ".pdf"))
    pyplot.close()


if __name__ == "__main__":
    main()
