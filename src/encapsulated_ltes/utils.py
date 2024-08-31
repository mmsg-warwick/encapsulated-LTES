#
# Auxiliary functions for the project
#
from pathlib import Path

import matplotlib as mpl


def root_dir():
    return Path(__file__).resolve().parents[2]


def get_interface_position(solution):
    pass


def set_plotting_format(mode="presentation"):
    mpl.style.use(["science", "vibrant"])

    if mode == "presentation":
        mpl.rcParams.update(
            {
                "font.family": "sans-serif",
                "text.usetex": False,
                "font.size": 10,
                "axes.labelsize": 12,
                "lines.linewidth": 2,
            }
        )

    elif mode == "paper":
        mpl.rcParams.update(
            {
                "font.family": "sans-serif",
                "text.usetex": False,
                "font.size": 6,
                "axes.labelsize": 8,
            }
        )

    else:
        raise KeyError("Mode should be either 'presentation' or 'paper'.")
