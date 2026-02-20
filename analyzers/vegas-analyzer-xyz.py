"""
VEGAS XYZ Analyzer - Extended magnetization analysis for HDF5 output files.
"""

from typing import List, Tuple
import numpy as np
import matplotlib

matplotlib.use("pdf")
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import click
import h5py


@click.command()
@click.argument("file", type=click.Path(exists=True))
def main(file: str) -> None:
    """Analyze VEGAS simulation output and generate comprehensive plots."""
    pp = PdfPages(file.replace(".h5", "_xyz.pdf"))

    with h5py.File(file, mode="r") as dataset:
        mcs = int(dataset.attrs.get("mcs", 0))
        kb = float(dataset.attrs.get("kb", 1.0))
        seed = dataset.attrs.get("seed")

        tau = mcs // 2

        temps = dataset.get("temperature")[:]
        fields = dataset.get("field")[:]

        energy = dataset.get("energy")[:, tau:]
        magx = dataset.get("magnetization_x")[:, tau:]
        magy = dataset.get("magnetization_y")[:, tau:]
        magz = dataset.get("magnetization_z")[:, tau:]

    magx_mean = np.mean(magx, axis=1)
    magx_abs_mean = np.abs(np.mean(magx, axis=1))

    magy_mean = np.mean(magy, axis=1)
    magy_abs_mean = np.abs(np.mean(magy, axis=1))

    magz_mean = np.mean(magz, axis=1)
    magz_abs_mean = np.abs(np.mean(magz, axis=1))

    mag = np.array([magx, magy, magz])
    mag_mean = np.mean(np.linalg.norm(mag, axis=0), axis=1)

    energy_mean = np.mean(energy, axis=1)

    # Safe division with zero temperature check
    temps_safe = np.where(temps > 0, temps, 1.0)  # Avoid division by zero
    chix = np.std(magx, axis=1) ** 2 / (kb * temps_safe)
    chiy = np.std(magy, axis=1) ** 2 / (kb * temps_safe)
    chiz = np.std(magz, axis=1) ** 2 / (kb * temps_safe)

    chi = np.std(np.linalg.norm(mag, axis=0), axis=1) ** 2 / (kb * temps_safe)

    cv = np.std(energy, axis=1) ** 2 / (kb * temps_safe**2)

    for xlabel, xarr in zip(["T", "H"], [temps, fields]):
        fig = pyplot.figure(figsize=(20, 10))
        ax1 = pyplot.subplot2grid((2, 5), (0, 0))
        ax2 = pyplot.subplot2grid((2, 5), (0, 1))
        ax3 = pyplot.subplot2grid((2, 5), (0, 2))
        ax4 = pyplot.subplot2grid((2, 5), (0, 3))
        ax5 = pyplot.subplot2grid((2, 5), (0, 4))
        ax6 = pyplot.subplot2grid((2, 5), (1, 0))
        ax7 = pyplot.subplot2grid((2, 5), (1, 1))
        ax8 = pyplot.subplot2grid((2, 5), (1, 2))
        ax9 = pyplot.subplot2grid((2, 5), (1, 3))
        ax10 = pyplot.subplot2grid((2, 5), (1, 4))

        ax1.plot(xarr, magx_mean, "-g", label=r"$\left< M_{x} \right>$")
        ax1.plot(
            xarr, magx_abs_mean, "-r", label=r"$ \left | \left< M_{x} \right> \right |$"
        )

        ax2.plot(xarr, magy_mean, "-g", label=r"$\left< M_{y} \right>$")
        ax2.plot(
            xarr, magy_abs_mean, "-r", label=r"$ \left | \left< M_{y} \right> \right |$"
        )

        ax3.plot(xarr, magz_mean, "-g", label=r"$\left< M_{z} \right>$")
        ax3.plot(
            xarr, magz_abs_mean, "-r", label=r"$ \left | \left< M_{z} \right> \right |$"
        )

        ax4.plot(xarr, mag_mean, "-g", label=r"$\left< M \right>$")

        ax5.plot(xarr, energy_mean, "-g", label=r"$\left< E \right>$")

        ax6.plot(xarr, chix, "-g", label=r"$\chi_{x}$")
        ax7.plot(xarr, chiy, "-g", label=r"$\chi_{y}$")
        ax8.plot(xarr, chiz, "-g", label=r"$\chi_{z}$")

        ax9.plot(xarr, chi, "-g", label=r"$\chi$")
        ax10.plot(xarr, cv, "-g", label=r"$C_{v}$")

        for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]:
            ax.grid()
            ax.legend(loc="best", fontsize=20)
            ax.set_xlabel(r"$%s$" % xlabel, fontsize=20)

        pyplot.tight_layout()
        pyplot.savefig(pp, format="pdf")
        pyplot.close()

    pp.close()

    with open(file.replace(".h5", "_xyz.mean"), mode="w") as outdata:
        outdata.write("# seed = %s\n" % seed)
        outdata.write("# T H E Cv Mx My Mz M Chix Chiy Chiz Chi\n")
        for i in range(len(temps)):
            outdata.write(
                "{} {} {} {} {} {} {} {} {} {} {} {}\n".format(
                    temps[i],
                    fields[i],
                    energy_mean[i],
                    cv[i],
                    magx_mean[i],
                    magy_mean[i],
                    magz_mean[i],
                    mag_mean[i],
                    chix[i],
                    chiy[i],
                    chiz[i],
                    chi[i],
                )
            )


if __name__ == "__main__":
    main()
