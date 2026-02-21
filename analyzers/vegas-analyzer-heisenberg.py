"""
VEGAS Heisenberg Analyzer - Comprehensive analysis for Heisenberg model simulations.
"""

from typing import Dict, List, Tuple, Set
import h5py
import numpy as np
import os
import matplotlib

matplotlib.use("PDF")
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import click


def get_equals(lista: List[float]) -> List[Tuple[int, ...]]:
    """Find groups of equal consecutive values in a list."""
    elements: List[Set[int]] = [{i} for i in range(len(lista))]
    for i, val1 in enumerate(lista):
        for j in range(i + 1, len(lista)):
            if val1 == lista[j]:
                elements[i].add(j)
            else:
                break

    lista2 = list(reversed(lista))
    elements2: List[Set[int]] = [{len(lista) - 1 - i} for i in range(len(lista))]
    for i, val1 in enumerate(lista2):
        for j in range(i + 1, len(lista2)):
            if val1 == lista2[j]:
                elements2[i].add(len(lista) - 1 - j)
            else:
                break
    elements2 = list(reversed(elements2))
    final_elements = [
        tuple(sorted(elements[i] | elements2[i])) for i in range(len(lista))
    ]
    return sorted(set(final_elements))


@click.command()
@click.argument("file", type=click.Path(exists=True))
def main(file: str) -> None:
    """Analyze VEGAS Heisenberg simulation output."""
    out = os.path.abspath(os.path.join(file, os.pardir))

    with h5py.File(file, "r") as data:
        mcs = int(data.attrs["mcs"])
        seed = data.attrs["seed"]
        kb = float(data.attrs.get("kb", 1.0))
        temps_full = data.get("temperature")[:]
        fields_full = data.get("field")[:]
        num_sites = len(data.get("positions"))

        types = np.unique(data.get("types")[:])
        types = [t.decode("utf-8") for t in types]

        num_sites_types: Dict[str, int] = {
            t: len(data.get("types")[data.get("types")[:] == t.encode()]) for t in types
        }

        join = [(temps_full[i], fields_full[i]) for i in range(len(temps_full))]
        equals = get_equals([(float(t), float(f)) for t, f in join])

        temps: List[float] = []
        fields: List[float] = []
        for i in equals:
            temps.append(float(temps_full[i[0]]))
            fields.append(float(fields_full[i[0]]))
            for j in i:
                for k in i:
                    assert (temps_full[j], fields_full[j]) == (
                        temps_full[k],
                        fields_full[k],
                    )

        temps_arr = np.array(temps)
        fields_arr = np.array(fields)

        mags_x: List[List[float]] = [[] for _ in range(len(equals))]
        mags_y: List[List[float]] = [[] for _ in range(len(equals))]
        mags_z: List[List[float]] = [[] for _ in range(len(equals))]
        mags_types_x: Dict[str, List[List[float]]] = {
            t: [[] for _ in range(len(equals))] for t in types
        }
        mags_types_y: Dict[str, List[List[float]]] = {
            t: [[] for _ in range(len(equals))] for t in types
        }
        mags_types_z: Dict[str, List[List[float]]] = {
            t: [[] for _ in range(len(equals))] for t in types
        }
        energy: List[List[float]] = [[] for _ in range(len(equals))]

        for i, eq in enumerate(equals):
            for j in eq:
                mags_x[i] += list(data.get("magnetization_x")[j, :])
                mags_y[i] += list(data.get("magnetization_y")[j, :])
                mags_z[i] += list(data.get("magnetization_z")[j, :])
                energy[i] += list(data.get("energy")[j, :])
                for t in types:
                    mags_types_x[t][i] += list(data.get("%s_x" % t)[j, :])
                    mags_types_y[t][i] += list(data.get("%s_y" % t)[j, :])
                    mags_types_z[t][i] += list(data.get("%s_z" % t)[j, :])

        mags_x = [np.array(M) for M in mags_x]
        mags_y = [np.array(M) for M in mags_y]
        mags_z = [np.array(M) for M in mags_z]
        energy = [np.array(M) for M in energy]
        for t in types:
            mags_types_x[t] = [np.array(M) for M in mags_types_x[t]]
            mags_types_y[t] = [np.array(M) for M in mags_types_y[t]]
            mags_types_z[t] = [np.array(M) for M in mags_types_z[t]]

        tau = [len(M) // 5 for M in energy]

        mean_mags = np.array(
            [
                np.mean(
                    np.linalg.norm([mags_x[i], mags_y[i], mags_z[i]], axis=0)[tau[i] :]
                )
                / num_sites
                for i, _ in enumerate(mags_x)
            ]
        )

        mean_mags_types: Dict[str, np.ndarray] = {
            t: np.array(
                [
                    np.mean(
                        np.linalg.norm(
                            [
                                mags_types_x[t][i],
                                mags_types_y[t][i],
                                mags_types_z[t][i],
                            ],
                            axis=0,
                        )[tau[i] :]
                    )
                    / num_sites_types[t]
                    for i, _ in enumerate(mags_types_x[t])
                ]
            )
            for t in types
        }

        mean_mags_z = np.array(
            [np.mean(mags_z[i][tau[i] :]) / num_sites for i, _ in enumerate(mags_z)]
        )
        mean_mags_types_z: Dict[str, np.ndarray] = {
            t: np.array(
                [
                    np.mean(mags_types_z[t][i][tau[i] :]) / num_sites_types[t]
                    for i, _ in enumerate(mags_types_z[t])
                ]
            )
            for t in types
        }

        mean_ene = np.array(
            [np.mean(energy[i][tau[i] :]) / num_sites for i, _ in enumerate(energy)]
        )

        temps_safe = np.where(temps_arr > 0, temps_arr, 1.0)

        susceptibility = np.array(
            [
                (
                    np.mean(
                        mags_x[i][tau[i] :] ** 2
                        + mags_y[i][tau[i] :] ** 2
                        + mags_z[i][tau[i] :] ** 2
                    )
                    - np.mean(
                        np.sqrt(
                            mags_x[i][tau[i] :] ** 2
                            + mags_y[i][tau[i] :] ** 2
                            + mags_z[i][tau[i] :] ** 2
                        )
                    )
                    ** 2
                )
                / (kb * temps_safe[i])
                / num_sites
                for i, _ in enumerate(mags_x)
            ]
        )

        susceptibility_types: Dict[str, np.ndarray] = {
            t: np.array(
                [
                    (
                        np.mean(
                            mags_types_x[t][i][tau[i] :] ** 2
                            + mags_types_y[t][i][tau[i] :] ** 2
                            + mags_types_z[t][i][tau[i] :] ** 2
                        )
                        - np.mean(
                            np.sqrt(
                                mags_types_x[t][i][tau[i] :] ** 2
                                + mags_types_y[t][i][tau[i] :] ** 2
                                + mags_types_z[t][i][tau[i] :] ** 2
                            )
                        )
                        ** 2
                    )
                    / (kb * temps_safe[i])
                    / num_sites_types[t]
                    for i, _ in enumerate(mags_x)
                ]
            )
            for t in types
        }

        specific_heat = np.array(
            [
                np.std(energy[i][tau[i] :]) ** 2 / (kb * temps_safe[i] ** 2) / num_sites
                for i, _ in enumerate(energy)
            ]
        )

        pp = PdfPages(file.replace(".h5", ".pdf"))

        fig = pyplot.figure(figsize=(16, 9))
        ax = fig.add_subplot(221)
        ax.plot(temps_arr, mean_mags, label=r"$M_{T}$")
        for t in types:
            ax.plot(temps_arr, mean_mags_types[t], label=r"$M_{%s}$" % t)
        ax.grid()
        ax.legend(loc="best", fontsize=4)
        ax.set_xlabel(r"$T$", fontsize=30)
        ax.set_ylabel(r"$M$", fontsize=30)
        ax.set_title("Magnetization")

        text = "Seed = {}\nMCS = {}".format(seed, mcs)
        props = dict(boxstyle="round", facecolor="yellowgreen", alpha=1.0)
        ax.text(
            pyplot.xlim()[-1],
            pyplot.ylim()[-1],
            text,
            fontsize=4,
            horizontalalignment="right",
            verticalalignment="top",
            bbox=props,
            fontweight="bold",
        )

        ax = fig.add_subplot(222)
        ax.plot(temps_arr, mean_ene)
        ax.grid()
        ax.set_xlabel(r"$T$", fontsize=30)
        ax.set_ylabel(r"$E$", fontsize=30)
        ax.set_title("Energy")

        ax = fig.add_subplot(223)
        ax.plot(temps_arr, susceptibility, label=r"$\chi_{T}$")
        for t in types:
            ax.plot(temps_arr, susceptibility_types[t], label=r"$\chi_{%s}$" % t)
        ax.grid()
        ax.legend(loc="best", fontsize=4)
        ax.set_xlabel(r"$T$", fontsize=30)
        ax.set_ylabel(r"$\chi$", fontsize=30)
        ax.set_title("Susceptibility")

        ax = fig.add_subplot(224)
        ax.plot(temps_arr, specific_heat)
        ax.grid()
        ax.set_xlabel(r"$T$", fontsize=30)
        ax.set_ylabel(r"$C_V$", fontsize=30)
        ax.set_title("Specific heat")

        fig.subplots_adjust(hspace=0.5)
        fig.subplots_adjust(wspace=0.5)

        pyplot.savefig(pp, format="pdf")

        fig = pyplot.figure(figsize=(16, 9))
        ax = fig.add_subplot(221)
        ax.plot(fields_arr, mean_mags_z, label=r"$M_{T}$")
        for t in types:
            ax.plot(fields_arr, mean_mags_types_z[t], label=r"$M_{%s}$" % t)
        ax.grid()
        ax.legend(loc="best", fontsize=4)
        ax.set_xlabel(r"$H$", fontsize=30)
        ax.set_ylabel(r"$M$", fontsize=30)
        ax.set_title("Magnetization")

        ax = fig.add_subplot(222)
        ax.plot(fields_arr, mean_ene)
        ax.grid()
        ax.set_xlabel(r"$H$", fontsize=30)
        ax.set_ylabel(r"$E$", fontsize=30)
        ax.set_title("Energy")

        ax = fig.add_subplot(223)
        ax.plot(fields_arr, susceptibility, label=r"$\chi_{T}$")
        for t in types:
            ax.plot(fields_arr, susceptibility_types[t], label=r"$\chi_{%s}$" % t)
        ax.grid()
        ax.legend(loc="best", fontsize=4)
        ax.set_xlabel(r"$H$", fontsize=30)
        ax.set_ylabel(r"$\chi$", fontsize=30)
        ax.set_title("Susceptibility")

        ax = fig.add_subplot(224)
        ax.plot(fields_arr, specific_heat)
        ax.grid()
        ax.set_xlabel(r"$H$", fontsize=30)
        ax.set_ylabel(r"$C_V$", fontsize=30)
        ax.set_title("Specific heat")

        fig.subplots_adjust(hspace=0.5)
        fig.subplots_adjust(wspace=0.5)

        pyplot.savefig(pp, format="pdf")

        colors = [
            "purple",
            "pink",
            "gold",
            "black",
            "skyblue",
            "crimson",
            "green",
            "blue",
            "orange",
        ]
        colors_types: Dict[str, str] = {t: colors.pop() for t in types}

        fig = pyplot.figure(figsize=(16, 16))
        pos = data.get("positions")[data.get("positions")[:, 2] == 0, :]
        sub_types = data.get("types")[data.get("positions")[:, 2] == 0]
        colors_list = [colors_types[t.decode("utf-8")] for t in sub_types]
        pyplot.scatter(pos[:, 0], pos[:, 1], c=colors_list, s=500)
        pyplot.xlabel(r"$x$", fontsize=4)
        pyplot.ylabel(r"$y$", fontsize=4)
        pyplot.savefig(pp, format="pdf")

        def get_indices(
            temps_arr: np.ndarray, fields_arr: np.ndarray
        ) -> Tuple[List[int], List[str]]:
            index: List[int] = []
            titles: List[str] = []
            for i in range(1, len(temps_arr)):
                if temps_arr[i] == temps_arr[i - 1]:
                    index.append(i - 1)
                    titles.append("Start")
                    break

            index.append(int(np.argmin(fields_arr)))
            titles.append("")

            index.append(-1)
            titles.append("Final")

            return index, titles

        finalstates = np.array([data.get("finalstates")[i[-1], :] for i in equals])
        indices, titles = get_indices(temps_arr, fields_arr)

        for n, i in enumerate(indices):
            fig = pyplot.figure(figsize=(32, 16))
            fig.suptitle(
                r"$T=%f \hspace{2} H=%f \hspace{2} %s$"
                % (temps_arr[i], fields_arr[i], titles[n]),
                fontsize=50,
            )

            ax = fig.add_subplot(121)
            pos = data.get("positions")[data.get("positions")[:, 1] == 0, :]
            sub_types = data.get("types")[data.get("positions")[:, 1] == 0]
            colors_list = [colors_types[t.decode("utf-8")] for t in sub_types]
            spins = finalstates[i, data.get("positions")[:, 1] == 0, :]

            ax.quiver(
                pos[:, 0],
                pos[:, 2],
                spins[:, 0],
                spins[:, 2],
                pivot="middle",
                color=colors_list,
            )

            ax.set_xlabel(r"$x$", fontsize=30)
            ax.set_ylabel(r"$z$", fontsize=30)
            ax.set_title("XZ view for y = 0")
            ax.set_aspect("equal")

            ax = fig.add_subplot(122)
            pos = data.get("positions")[data.get("positions")[:, 0] == 0, :]
            sub_types = data.get("types")[data.get("positions")[:, 0] == 0]
            colors_list = [colors_types[t.decode("utf-8")] for t in sub_types]
            spins = finalstates[i, data.get("positions")[:, 0] == 0, :]

            ax.quiver(
                pos[:, 1],
                pos[:, 2],
                spins[:, 1],
                spins[:, 2],
                pivot="middle",
                color=colors_list,
            )

            ax.set_xlabel(r"$y$", fontsize=30)
            ax.set_ylabel(r"$z$", fontsize=30)
            ax.set_title("YZ view for x = 0")
            ax.set_aspect("equal")

            pyplot.savefig(pp, format="pdf")

        pp.close()

    with open(file.replace(".h5", ".mean"), mode="w") as file_:
        file_.write("# seed = {}\n".format(seed))
        file_.write("#\tT\tH\tE\tCv\tM\tMz\tX\t")
        for t in types:
            file_.write("%s\t%sz\tX_%s\t" % (t, t, t))
        file_.write("\n")
        for i, T in enumerate(temps_arr):
            file_.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                    T,
                    fields_arr[i],
                    mean_ene[i],
                    specific_heat[i],
                    mean_mags[i],
                    mean_mags_z[i],
                    susceptibility[i],
                )
            )
            for t in types:
                file_.write(
                    "{}\t{}\t{}\t".format(
                        mean_mags_types[t][i],
                        mean_mags_types_z[t][i],
                        susceptibility_types[t][i],
                    )
                )
            file_.write("\n")


if __name__ == "__main__":
    main()
