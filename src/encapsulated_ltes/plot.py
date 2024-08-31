#
# Plotting methods
#
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from .utils import set_plotting_format


def plot_comparison_data(simulation, datasets, xs=None, plotting_format="paper"):
    set_plotting_format(plotting_format)

    if xs is None:
        xs = [0.25, 0.5, 0.75, 1]
    fig, axes = plt.subplots(2, 2, figsize=(5.5, 3.5), sharey=True, sharex=True)

    L = simulation.parameter_values["Pipe length [m]"]
    R = simulation.parameter_values["Capsule radius [m]"]

    solution = simulation.solution

    for data, x in zip(datasets["PCM"], xs):
        time = solution["Time [s]"].entries
        T_HTF = solution["Phase-change material temperature [degC]"](t=time, x=x*L, r=0.8*R)
        axes[1, 0].plot(time, T_HTF)
        axes[1, 1].plot(data["Time [min]"] * 60, data["PCM Temperature [degC]"], '.-', label=f"{x:.2f}L")

    for data, x in zip(datasets["HTF"], xs):
        time = solution["Time [s]"].entries
        T_HTF = solution["Heat transfer fluid temperature [degC]"](t=time, x=x*L)
        axes[0, 0].plot(time, T_HTF, label=f"{x:.2f}L")
        axes[0, 1].plot(data["Time [min]"] * 60, data["HTF Temperature [degC]"], '.-', label=f"{x:.2f}L")

    for ax in axes[1, :]:
        ax.set_xlabel("Time [s]")

    for (ax, letter) in zip(axes.flatten(), ["a", "b", "c", "d"]):
        ax.text(0, 70, f"({letter})", horizontalalignment='left', verticalalignment='top',)

    axes[1, 0].set_ylabel("PCM temperature [°C]")
    axes[0, 0].set_ylabel("HTF temperature [°C]")
    axes[0, 0].set_title("Full model")
    axes[0, 1].set_title("Experimental data")
    axes[0, 0].legend()

    return fig, axes

def compare_0D_variables(simulations, output_variables=None, variable_names=None, plotting_format="paper"):
    if output_variables is None:
        output_variables = [
            "Outlet temperature [degC]",
            "X-averaged state of charge",
            "Stored energy per unit area [J.m-2]",
            "Relative error in energy conservation [%]",
        ]
        if variable_names is None:
            variable_names = [
                "Outlet temperature [°C]",
                "State of charge",
                "Stored energy per\n unit area [J m${}^{-2}$]",
                "Relative error in\n energy conservation [%]",
            ]
    if variable_names is None:
        variable_names = output_variables

    if not isinstance(simulations, list):
        simulations = [simulations]

    set_plotting_format(plotting_format)
    N_rows = math.ceil(len(output_variables)/ 2)
    fig, axes = plt.subplots(N_rows, 2, figsize=(5.5, 0.5 + 1.5 * N_rows), sharey=False, sharex=True)

    if len(simulations) > 1:
        for i, var_name in enumerate(output_variables):
            ax = axes.flat[i]
            time = simulations[1].solution["Time [s]"].data
            var = simulations[1].solution[var_name].data
            ax.plot(time, var, "k--", label=simulations[1].model.name)

    for i, var_name in enumerate(output_variables):
        ax = axes.flat[i]
        time = simulations[0].solution["Time [s]"].data
        var = simulations[0].solution[var_name].data
        ax.plot(time, var, label=simulations[0].model.name)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(variable_names[i])

    axes[0, 0].legend()
    fig.tight_layout()
    return fig, axes

def compare_1D_variables(simulations, output_variables=None, variable_names=None, times=5, plotting_format="paper"):
    if output_variables is None:
        output_variables = [
            "Heat transfer fluid temperature [degC]",
            "Phase-change material surface temperature [degC]",
        ]
        if variable_names is None:
            variable_names = [
                "HTF temperature [°C]",
                "PCM surface\n temperature [°C]",
            ]
    if variable_names is None:
        variable_names = output_variables

    if not isinstance(simulations, list):
        simulations = [simulations]

    if not isinstance(times, list):
        if isinstance(times, int):
            end_time = simulations[0].solution["Time [s]"].data[-1]
            times = [end_time * i / (times - 1) for i in range(times)]
        else:
            raise ValueError("times must be an integer or a list")


    set_plotting_format(plotting_format)
    N_rows = math.ceil(len(output_variables)/ 2)
    fig, axes = plt.subplots(N_rows, 2, figsize=(5.5, 0.5 + 1.5 * N_rows), sharey=False, sharex=True)

    viridis = cm.get_cmap('viridis')
    colours = viridis(np.linspace(0, 0.9, len(times)))

    if len(simulations) > 1:
        for i, var_name in enumerate(output_variables):
            label = simulations[1].model.name
            for t in times:
                ax = axes.flat[i]
                time = simulations[1].solution["x [m]"](t=t)
                var = simulations[1].solution[var_name](t=t)
                ax.plot(time, var, "k--", label=label)
                label = None

    for i, var_name in enumerate(output_variables):
        label = simulations[0].model.name
        for t, c in zip(times, colours):
            ax = axes.flat[i]
            time = simulations[0].solution["x [m]"](t=t)
            var = simulations[0].solution[var_name](t=t)
            ax.plot(time, var, color=c, label=label)
            label = None

            ax.set_xlabel("z [m]")
            ax.set_ylabel(variable_names[i])

    axes.flat[0].legend()
    fig.tight_layout()
    return fig, axes

def compare_2D_variables(simulations, output_variable=None, variable_name=None, times=4, xs=5, plotting_format="paper"):
    if output_variable is None:
        output_variable = "Phase-change material enthalpy [J.m-3]"
        if variable_name is None:
            variable_name = "Phase-change material\n enthalpy [J m${}^{-3}$]"
    if variable_name is None:
        variable_name = output_variable

    if not isinstance(times, list):
        if isinstance(times, int):
            end_time = simulations[0].solution["Time [s]"].data[-1]
            times = [end_time * i / (times - 1) for i in range(times)]
        else:
            raise ValueError("times must be an integer or a list")

    if not isinstance(xs, list):
        if isinstance(xs, int):
            Z = simulations[0].parameter_values["Pipe length [m]"]
            xs = [Z * i / (xs - 1) for i in range(xs)]
        else:
            raise ValueError("xs must be an integer or a list")

    set_plotting_format(plotting_format)
    N_rows = math.ceil(len(times)/ 2)
    fig, axes = plt.subplots(N_rows, 2, figsize=(5.5, 0.5 + 1.5 * N_rows), sharey=True, sharex=True)

    for t, ax in zip(times, axes.flat):
        if len(simulations) > 1:
            label = simulations[1].model.name
            for x in xs:
                r = simulations[1].solution["r [mm]"](t=t, x=Z/2)
                var = simulations[1].solution[output_variable](t=t, x=x)
                ax.plot(r, var, "lightgray", label=label)
                label = None

        r = simulations[0].solution["r [mm]"](t=t, x=Z/2)
        var = simulations[0].solution[output_variable](t=t, x=Z/2)
        ax.plot(r, var, label=simulations[0].model.name)
        ax.set_title(f"t = {t:.0f} s")

        ax.set_xlabel("r [mm]")
        ax.set_ylabel(variable_name)

    axes.flat[0].legend()
    fig.tight_layout()
    return fig, axes

def draw_loglog_slope(
    fig,
    ax,
    origin,
    width_inches,
    slope,
    inverted=False,
    color=None,
    polygon_kwargs=None,
    labelcolor=None,
    label_kwargs=None,
    zorder=None,
):
    """
    Adapted from https://gist.github.com/w1th0utnam3/a0189dc8a2c067ccb56e1de8c317b190

    This function draws slopes or "convergence triangles" into loglog plots.
    @param fig: The figure
    @param ax: The axes object to draw to
    @param origin: The 2D origin (usually lower-left corner) coordinate of the triangle
    @param width_inches: The width in inches of the triangle
    @param slope: The slope of the triangle, i.e. order of convergence
    @param inverted: Whether to mirror the triangle around the origin, i.e. whether
        it indicates the slope towards the lower left instead of upper right (defaults to false)
    @param color: The color of the of the triangle edges (defaults to default color)
    @param polygon_kwargs: Additional kwargs to the Polygon draw call that creates the slope
    @param label: Whether to enable labeling the slope (defaults to true)
    @param labelcolor: The color of the slope labels (defaults to the edge color)
    @param label_kwargs: Additional kwargs to the Annotation draw call that creates the labels
    @param zorder: The z-order value of the triangle and labels, defaults to a high value
    """

    if polygon_kwargs is None:
        polygon_kwargs = {}
    if label_kwargs is None:
        label_kwargs = {}

    if color is not None:
        polygon_kwargs["color"] = color
    if "linewidth" not in polygon_kwargs:
        polygon_kwargs["linewidth"] = 0.75 * mpl.rcParams["lines.linewidth"]
    if labelcolor is not None:
        label_kwargs["color"] = labelcolor
    if "color" not in label_kwargs:
        label_kwargs["color"] = polygon_kwargs["color"]
    if "fontsize" not in label_kwargs:
        label_kwargs["fontsize"] = 0.8 * mpl.rcParams["font.size"]

    if inverted:
        width_inches = -width_inches
    if zorder is None:
        zorder = 10

    # For more information on coordinate transformations in Matplotlib see
    # https://matplotlib.org/3.1.1/tutorials/advanced/transforms_tutorial.html

    # Convert the origin into figure coordinates in inches
    origin_disp = ax.transData.transform(origin)
    origin_dpi = fig.dpi_scale_trans.inverted().transform(origin_disp)

    # Obtain the bottom-right corner in data coordinates
    corner_dpi = origin_dpi + width_inches * np.array([1.0, 0.0])
    corner_disp = fig.dpi_scale_trans.transform(corner_dpi)
    corner = ax.transData.inverted().transform(corner_disp)

    (x1, y1) = (origin[0], origin[1])
    x2 = corner[0]

    # The width of the triangle in data coordinates
    width = x2 - x1
    # Compute offset of the slope
    log_offset = y1 / (x1**slope)

    y2 = log_offset * ((x1 + width) ** slope)

    # The vertices of the slope
    a = origin
    b = corner
    c = [x2, y2]

    # Draw the slope triangle
    X = np.array([a, b, c])
    triangle = plt.Polygon(X[:3, :], fill=True, zorder=zorder, **polygon_kwargs)
    ax.add_patch(triangle)

    # Convert vertices into display space
    a_disp = ax.transData.transform(a)
    b_disp = ax.transData.transform(b)
    c_disp = ax.transData.transform(c)

    # Figure out the center of the triangle sides in display space
    bottom_center_disp = a_disp + 0.5 * (b_disp - a_disp)
    bottom_center = ax.transData.inverted().transform(bottom_center_disp)

    right_center_disp = b_disp + 0.5 * (c_disp - b_disp)
    right_center = ax.transData.inverted().transform(right_center_disp)

    # Label alignment depending on inversion parameter
    va_xlabel = "top" if (not inverted and slope > 0) else "bottom"
    ha_ylabel = "left" if not inverted else "right"

    # Label offset depending on inversion parameter
    offset_xlabel = (
        [0.0, -0.33 * label_kwargs["fontsize"]]
        if (not inverted and slope > 0)
        else [0.0, 0.33 * label_kwargs["fontsize"]]
    )
    offset_ylabel = (
        [0.33 * label_kwargs["fontsize"], 0.0]
        if not inverted
        else [-0.33 * label_kwargs["fontsize"], 0.0]
    )

    # Draw the slope labels
    ax.annotate(
        "$1$",
        bottom_center,
        xytext=offset_xlabel,
        textcoords="offset points",
        ha="center",
        va=va_xlabel,
        zorder=zorder,
        **label_kwargs,
    )
    ax.annotate(
        f"${slope}$",
        right_center,
        xytext=offset_ylabel,
        textcoords="offset points",
        ha=ha_ylabel,
        va="center",
        zorder=zorder,
        **label_kwargs,
    )
