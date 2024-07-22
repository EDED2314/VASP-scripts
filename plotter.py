# plot_potential_surface() repurposed from https://corinwagen.github.io/public/blog/20220905_pypes.html
# plot_reaction_coord repurposed from https://github.com/Paci-Group/Reaction-Coordinate-Plotter/blob/main/plot_rxn_coords.py

# get packages
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from PIL import Image


def plot_potential_surface(
    X,
    Y,
    labels,
    xlabel: str = "Reaction Coordinate",
    ylabel: str = "Potential Energy (eV)",
):
    """
    x and y positions. y in any units you want, if you want, and x in the range [0,1].

    Y = [2.49, 3.5, 0, 20.2, 19, 21.5, 20, 20.3, -5]

    X = [0, 0.15, 0.3, 0.48, 0.55, 0.63, 0.70, 0.78, 1]

    labels for points. False if you don't want a label

    label = [
        "label1",
        False,
        "label2",
        "label3",
        "label4",
        "label5",
        "label6",
        "label7",
        "label8",
    ]
    """
    label = labels

    # make matplotlib look good
    plt.rc("font", size=11, family="serif")
    plt.rc("axes", titlesize=12, labelsize=12)
    plt.rc(["xtick", "ytick"], labelsize=11)
    plt.rc("legend", fontsize=12)
    plt.rc("figure", titlesize=14)

    TS = []
    for idx in range(len(Y)):
        if idx == 0 or idx == len(Y) - 1:
            TS.append(False)
        else:
            TS.append((Y[idx] > Y[idx + 1]) and (Y[idx] > Y[idx - 1]))

    # sanity checks
    assert len(X) == len(Y), "need X and Y to match length"
    assert len(X) == len(label), "need right number of labels"

    # now we start building the figure, axes first
    f = plt.figure(figsize=(8, 8))
    ax = f.gca()
    xgrid = np.linspace(0, 1, 1000)
    ax.spines[["right", "bottom", "top"]].set_visible(False)

    YMAX = 1.1 * max(Y) - 0.1 * min(Y)
    YMIN = 1.1 * min(Y) - 0.1 * max(Y)

    plt.xlim(-0.1, 1.1)
    plt.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    plt.ylim(bottom=YMIN, top=YMAX)
    ax.plot(-0.1, YMAX, "^k", clip_on=False)

    # label axes
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # plot the points
    plt.plot(X, Y, "o", markersize=7, c="black")

    # add labels
    for i in range(len(X)):
        if label[i]:
            delta_y = 0.6 if TS[i] else -1.2
            plt.annotate(
                label[i],
                (X[i], Y[i] + delta_y),
                fontsize=12,
                fontweight="normal",
                ha="center",
            )

    # add connecting lines
    for i in range(len(X) - 1):
        idxs = np.where(np.logical_and(xgrid >= X[i], xgrid <= X[i + 1]))
        smoother = interp.BPoly.from_derivatives(
            [X[i], X[i + 1]], [[y, 0] for y in [Y[i], Y[i + 1]]]
        )
        plt.plot(xgrid[idxs], smoother(xgrid[idxs]), ls="-", c="black", lw=2)

    # finish up!
    plt.tight_layout()
    plt.show()


def plotOnlyPoints(X, Y, xlabel, ylabel):
    plt.rc("font", size=11, family="serif")
    plt.rc("axes", titlesize=12, labelsize=12)
    plt.rc(["xtick", "ytick"], labelsize=11)
    plt.rc("legend", fontsize=12)
    plt.rc("figure", titlesize=14)

    # sanity checks
    assert len(X) == len(Y), "need X and Y to match length"

    # now we start building the figure, axes first
    f = plt.figure(figsize=(8, 8))
    ax = f.gca()
    xgrid = np.linspace(0, 1, 1000)
    ax.spines[["right", "bottom", "top"]].set_visible(False)

    YMAX = 1.1 * max(Y) - 0.1 * min(Y)
    YMIN = 1.1 * min(Y) - 0.1 * max(Y)

    plt.xlim(-0.1, 1.1)
    plt.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    plt.ylim(bottom=YMIN, top=YMAX)
    ax.plot(-0.1, YMAX, "^k", clip_on=False)

    # label axes
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # plot the points
    plt.plot(X, Y, "o", markersize=7, c="black")

    plt.tight_layout()
    plt.show()


def open_image_local(path_to_image):
    image = Image.open(path_to_image)  # Open the image
    image_array = np.array(image)  # Convert to a numpy array
    return image_array  # Output


def customPlot(
    X,
    Y1,
    Y2,
    labels1,
    labels2,
    images,
    xlabel: str = "Reaction Coordinate",
    ylabel: str = "Potential Energy (eV)",
    image_width=0.15,
    image_height=0.15,
):
    """
    x and y positions. y in any units you want, if you want, and x in the range [0,1].

    Y = [2.49, 3.5, 0, 20.2, 19, 21.5, 20, 20.3, -5]

    X = [0, 0.15, 0.3, 0.48, 0.55, 0.63, 0.70, 0.78, 1]

    labels for points. False if you don't want a label

    labelsh = [
        {"label": "WO3 + 2H", "pos": "B"},
        {"label": "WO3--H + H", "pos": "B"},
        {"label": "WO3--H2", "pos": "T"},
        {"label": "WO3 (vac) + H2O", "pos": "B"},
    ]
    """

    # make matplotlib look good
    plt.rc("font", size=11, family="serif")
    plt.rc("axes", titlesize=12, labelsize=12)
    plt.rc(["xtick", "ytick"], labelsize=11)
    plt.rc("legend", fontsize=12)
    plt.rc("figure", titlesize=14)

    # sanity checks
    assert len(X) == len(Y1), "need X and Y to match length"
    assert len(X) == len(Y2)
    assert len(X) == len(labels1), "need right number of labels"
    assert len(X) == len(labels2)

    # create figure and axis
    fig, ax = plt.subplots(figsize=(8, 8))
    xgrid = np.linspace(0, 1, 1000)
    ax.spines[["right", "bottom", "top"]].set_visible(False)

    YMAX = max(1.1 * max(Y1) - 0.1 * min(Y1), 1.1 * max(Y2) - 0.1 * min(Y2))
    YMIN = min(1.1 * min(Y1) - 0.1 * max(Y1), 1.1 * min(Y2) - 0.1 * max(Y2)) - 0.6

    ax.set_xlim(-0.1, 1.1)
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    ax.set_ylim(bottom=YMIN, top=YMAX)
    ax.plot(-0.1, YMAX, "^k", clip_on=False)

    # label axes
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # plot the points
    ax.plot(X, Y1, "o", markersize=7, c="black")
    ax.plot(X, Y2, "o", markersize=7, c="blue")

    # add labels
    for i in range(len(X)):
        if labels1[i]:
            # delta_y = 0.6 if TS1[i] else -1.2
            delta_y = 0.6 if labels1[i]["pos"] == "T" else -1.2
            ax.annotate(
                labels1[i]["label"],
                (X[i], Y1[i] + delta_y),
                fontsize=12,
                fontweight="normal",
                ha="center",
            )
        if labels2[i]:
            # delta_y = 0.6 if TS2[i] else -1.2
            delta_y = 0.6 if labels1[i]["pos"] == "T" else -1.2
            ax.annotate(
                labels2[i]["label"],
                (X[i], Y2[i] + delta_y),
                fontsize=12,
                fontweight="normal",
                ha="center",
            )

    # put images on graph
    for i in range(len(X)):
        delta_y = 2 if images[i]["pos"] == "T" else -2
        if images[i]["ref"] == 0:
            pos = (X[i] + images[i]["dis_x"], Y1[i] + delta_y)
        else:
            pos = (X[i] + images[i]["dis_x"], Y2[i] + delta_y)

        image = open_image_local(images[i]["img"])

        display_coords = ax.transData.transform(pos)
        figure_coords = fig.transFigure.inverted().transform(display_coords)
        image_xaxis, image_yaxis = figure_coords

        ax_image = fig.add_axes(
            [
                image_xaxis - image_width / 2,
                image_yaxis - image_height / 2,
                image_width,
                image_height,
            ]
        )
        # ax_image = fig.add_axes([image_xaxis, image_yaxis, image_width, image_height])
        ax_image.imshow(image)
        ax_image.axis("off")

    # add connecting lines
    for i in range(len(X) - 1):
        idxs = np.where(np.logical_and(xgrid >= X[i], xgrid <= X[i + 1]))
        smoother1 = interp.BPoly.from_derivatives(
            [X[i], X[i + 1]], [[y, 0] for y in [Y1[i], Y1[i + 1]]]
        )
        ax.plot(xgrid[idxs], smoother1(xgrid[idxs]), ls="-", c="black", lw=2)
        smoother2 = interp.BPoly.from_derivatives(
            [X[i], X[i + 1]], [[y, 0] for y in [Y2[i], Y2[i + 1]]]
        )
        ax.plot(xgrid[idxs], smoother2(xgrid[idxs]), ls="-", c="blue", lw=2)

    # finish up!
    fig.tight_layout()
    plt.show()


import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
from matplotlib import colormaps as cmap
import matplotlib
import numpy as np


def plot_rxn_coords(ax, energies, color, label, linewidth=2, scale=0.32):
    energies = np.array(energies)
    for j, l in enumerate(energies):
        ax.plot(
            [(j + 1 - scale), (j + 1 + scale)], [l, l], color=color, linewidth=linewidth
        )
        if j <= len(energies) - 2:
            ax.plot(
                [(j + 1 + scale), (j + 2 - scale)],
                [l, energies[j + 1]],
                linestyle=":",
                color=color,
                linewidth=linewidth,
            )
    ax.plot([], [], color=color, linewidth=linewidth, label=label)


def plot_rxn_delta_es(
    ax,
    energies,
    color,
    label,
    width=0.25,
    shift=0,
    add_labels=False,
    fmt="%0.2f",
    label_padding=3,
    add_zero_line=True,
    linewidth=1,
    linestyle=":",
):
    x = np.arange(len(energies) - 1) + 1
    deltas = energies[1:] - energies[:-1]
    rects = ax.bar(x + shift, deltas, width, label=label, color=color)
    if add_labels:
        ax.bar_label(rects, padding=label_padding, fmt=fmt)
    if add_zero_line:
        ax.axhline(linewidth=linewidth, linestyle=linestyle, color="k")


def plot_rxn_coord_custom(energies1, name, energies2, name2, index="1"):
    energies1 = np.array(energies1)
    energies2 = np.array(energies2)
    # energies = np.array([6.0, 5.5, 5.7, 5.0, 4.3, 3.4, 3.0])
    deltaLabels = [1, 2, 3], ["1-2", "2-3", "3-4"]
    standardLabels = [1, 2, 3, 4], ["1", "2", "3", "4"]

    # ##########################################################################################
    # # Plot Delta E's as Bar Grap
    # ##########################################################################################
    # fig = plt.figure(figsize=(8, 6))
    # ax = fig.add_subplot(111)
    # width = 0.5
    # plot_rxn_delta_es(
    #     ax,
    #     energies1,
    #     "red",
    #     "",
    #     width=width,
    #     shift=0,
    #     label_padding=5,
    #     add_labels=True,
    # )

    # ax.set_ylabel("$\Delta$ E (eV)", fontsize=16)
    # ax.set_xlabel("Step", fontsize=16)
    # x, labels = deltaLabels
    # ax.set_xticks(x, labels)
    # ax.tick_params(labelsize=14)
    # ax.legend(loc="lower left", frameon=False, fontsize=14)
    # plt.tight_layout()

    # fig.savefig(f"data/delta_es_{index}.png", dpi=300)
    # plt.close(fig)

    # ##########################################################################################
    # # Plot Standard Rxn Coordinate Diagram
    # ##########################################################################################
    # fig2 = plt.figure(figsize=(8, 6))
    # ax2 = fig2.add_subplot(111)
    # width = 0.5
    # plot_rxn_coords(ax2, energies1, "red", name, zero=0, linewidth=2, scale=0.32)

    # ax2.set_ylabel("E (eV)", fontsize=16)
    # ax2.set_xlabel("Reaction Coordinate", fontsize=16)
    # x, labels = standardLabels
    # ax2.set_xticks(x, labels)
    # ax2.tick_params(labelsize=14)
    # ax2.legend(loc="lower left", frameon=False, fontsize=14)
    # plt.tight_layout()

    # fig2.savefig(f"data/es_{index}.png", dpi=300)
    # plt.close(fig2)

    ##########################################################################################
    # Plot Rxn Coordinate Diagram with Delta E Inset
    ##########################################################################################
    fig3 = plt.figure(figsize=(8, 6))
    ax3 = fig3.add_subplot(111)
    width = 0.5
    plot_rxn_coords(
        ax3,
        np.array(energies1),
        "red",
        name,
        linewidth=2,
        scale=0.32,
    )
    plot_rxn_coords(
        ax3,
        np.array(energies2),
        "blue",
        name2,
        linewidth=2,
        scale=0.32,
    )

    ax3.set_ylabel("E (eV)", fontsize=16)
    ax3.set_xlabel("Reaction Coordinate", fontsize=16)
    x, labels = standardLabels
    ax3.set_xticks(x, labels)
    ax3.tick_params(labelsize=14)
    ax3.legend(loc="lower left", frameon=False, fontsize=14)

    ax4 = ax3.inset_axes([0.6, 0.7, 0.35, 0.25])
    plot_rxn_delta_es(
        ax4,
        energies1,
        "red",
        name,
        width=width,
        shift=0,
        label_padding=5,
        add_labels=False,
    )
    ax5 = ax3.inset_axes([0.6, 0.3, 0.35, 0.25])
    plot_rxn_delta_es(
        ax5,
        energies2,
        "blue",
        name,
        width=width,
        shift=0,
        label_padding=5,
        add_labels=False,
    )
    ax4.tick_params(labelsize=12)
    ax4.locator_params(axis="y", nbins=6)
    ax4.set_xticks(deltaLabels[0], deltaLabels[1])
    ax4.set_ylabel("$\Delta$ E (eV)", fontsize=12)
    ax4.set_xlabel("Step", fontsize=12)
    
    ax5.tick_params(labelsize=12)
    ax5.locator_params(axis="y", nbins=6)
    ax5.set_xticks(deltaLabels[0], deltaLabels[1])
    ax5.set_ylabel("$\Delta$ E (eV)", fontsize=12)
    ax5.set_xlabel("Step", fontsize=12)

    plt.tight_layout()

    fig3.savefig(f"data/es_delta_es_inset_{index}.png", dpi=300)
    plt.close(fig3)
