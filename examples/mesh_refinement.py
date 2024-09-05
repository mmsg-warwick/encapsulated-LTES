import gc
import math

import matplotlib.pyplot as plt
import numpy as np
import pybamm

import encapsulated_ltes as ltes

ltes.set_plotting_format("paper")

models = [
    ltes.ReducedModel(),
    ltes.FullModel(),
]

min_mesh = 0
max_mesh = 4

level_range = range(min_mesh, max_mesh + 1)


param = ltes.get_parameter_values("Nallusamy2007")
param["Heat transfer coefficient [W.m-2.K-1]"] = 1000

solutions = [[], []]
times = [[], []]

plt.figure(figsize=(2.75, 2))

for l in range(min_mesh, max_mesh + 1):
    for model, solution, time in zip(models, solutions, times):
        r = model.variables["r [m]"]
        x = model.variables["x [m]"]
        var_pts = {r: math.floor(10 * 2**l), x: math.floor(20 * 2**l)}
        sim = pybamm.Simulation(model, parameter_values=param, var_pts=var_pts)
        # sol = sim.solve(np.linspace(0, 10000, 2000))
        sol = sim.solve([0, 10000])

        del sim
        gc.collect()

        solution.append(sol)
        time.append(sol.solve_time.value)

        print(f"{model.name}: {sol.solve_time}")

        del sol
        gc.collect()


# Plot convergence of energy conservation
print("Generating energy conservation plot")
errors = [[], []]
for solution, error in zip(solutions, errors):
    for sol in solution:
        err = sol["Relative error in energy conservation [%]"].data.mean()
        error.append(err)

fig, axes = plt.subplots(1, 2, figsize=(5.5, 2))

mesh = [2**i for i in level_range]
for model, error, time in zip(models, errors, times):
    axes[0].loglog(mesh, error, ".-", label=model.name)
    axes[1].loglog(mesh, time, ".-", label=model.name)

axes[0].legend()
for ax in axes:
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Mesh refinement factor")
axes[0].set_ylabel("Relative error in\n energy conservation [%]")
axes[1].set_ylabel("Solve time [s]")
i_mid = (len(level_range) - 1) // 2
ltes.draw_loglog_slope(
    fig,
    axes[0],
    [mesh[i_mid], 1.2 * np.max([error[i_mid] for error in errors])],
    0.2,
    -1,
    inverted=False,
    color="lightgray",
    labelcolor="black",
)
ltes.draw_loglog_slope(
    fig,
    axes[1],
    [mesh[i_mid], 0.6 * times[0][i_mid]],
    0.2,
    1,
    inverted=False,
    color="lightgray",
    labelcolor="black",
)
ltes.draw_loglog_slope(
    fig,
    axes[1],
    [mesh[i_mid], 0.6 * times[1][i_mid]],
    0.2,
    3,
    inverted=False,
    color="lightgray",
    labelcolor="black",
)
fig.tight_layout()

fig.savefig(ltes.root_dir() / "figures" / "convergence_conservation_error.png", dpi=300)

print("Energy conservation plot saved")

gc.collect()

## Plot convergence of variables
x = np.linspace(0, param["Pipe length [m]"], 50)
r = np.linspace(0, param["Capsule radius [m]"], 50)
t = np.linspace(0, 10000, 100)
data = {}

# process HTF temperature
print("Generating HTF temperature plot")
errors = [[], []]
for solution, error in zip(solutions, errors):
    benchmark = solution[-1]["Heat transfer fluid temperature [K]"](t=t, x=x)
    for sol in solution[:-1]:
        result = sol["Heat transfer fluid temperature [K]"](t=t, x=x)
        err = np.sqrt(
            (((result - benchmark) ** 2).mean()) / ((benchmark**2).mean())
        )
        error.append(err)

data["HTF temperature [K]"] = errors

del result
gc.collect()

# process PCM temperature
print("Generating PCM temperature plot")
errors = [[], []]
for solution, error in zip(solutions, errors):
    benchmark = solution[-1]["Phase-change material temperature [K]"](
        t=t, x=x, r=r
    )
    for sol in solution[:-1]:
        result = sol["Phase-change material temperature [K]"](t=t, x=x, r=r)
        err = np.sqrt(((result - benchmark) ** 2).mean()) / np.sqrt(
            (benchmark**2).mean()
        )
        error.append(err)

data["PCM temperature [K]"] = errors

fig, axes = plt.subplots(1, 2, figsize=(5.5, 2), sharey=True)
for (var, errors), ax in zip(data.items(), axes):
    for model, error in zip(models, errors):
        ax.loglog([2**i for i in level_range[:-1]], error, ".-", label=model.name)
        ax.set_xscale("log", base=2)
        ax.set_xlabel("Mesh refinement factor")
        ax.set_ylabel("Relative error")
        ax.set_title(var)

i_mid = (len(level_range) - 2) // 2
ltes.draw_loglog_slope(
    fig,
    axes[0],
    [
        mesh[i_mid],
        1.2 * np.max([error[i_mid] for error in data["HTF temperature [K]"]]),
    ],
    0.2,
    -1,
    inverted=False,
    color="lightgray",
    labelcolor="black",
)
ltes.draw_loglog_slope(
    fig,
    axes[1],
    [
        mesh[i_mid],
        1.2 * np.max([error[i_mid] for error in data["PCM temperature [K]"]]),
    ],
    0.2,
    -1,
    inverted=False,
    color="lightgray",
    labelcolor="black",
)


axes[0].legend()
fig.tight_layout()

gc.collect()

print("Saving to file")

fig.savefig(ltes.root_dir() / "figures" / "convergence_variables.png", dpi=300)

print("Plots saved")
