#
# Full model for encapuslated LTES
#
import pybamm
from .base_LTES_model import BaseLTESModel


class ReducedModel(BaseLTESModel):
    def __init__(self, name="Reduced model"):
        super().__init__(name=name)

        param = self.param

        ######################
        # Variables
        ######################
        H_av = pybamm.Variable(
            "X-averaged phase-change material enthalpy [J.m-3]",
            domains={"primary": "capsule"},
            scale=param.T2H(param.T_0),
        )
        T_c_av = pybamm.Variable(
            "X-averaged phase-change material temperature [K]",
            domains={"primary": "capsule"},
            scale=abs(param.T_m - param.T_in),
            reference=param.T_m,
        )
        H = pybamm.SecondaryBroadcast(H_av, "pipe")
        T_c = pybamm.SecondaryBroadcast(T_c_av, "pipe")

        ######################
        # Governing equations
        ######################
        self.rhs = {H_av: pybamm.div(param.k(H_av) * pybamm.grad(T_c_av))}
        self.algebraic = {T_c_av: T_c_av - param.H2T(H_av)}

        self.boundary_conditions = {
            T_c_av: {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (
                    -param.h
                    / param.k(pybamm.surf(H_av))
                    * (pybamm.surf(T_c_av) - param.T_in),
                    "Neumann",
                ),
            },
        }

        self.initial_conditions = {H_av: param.T2H(param.T_0), T_c_av: param.T_0}

        ######################
        # (Some) variables
        ######################
        q = -param.k(pybamm.surf(H)) * pybamm.BoundaryGradient(T_c, "right")
        T_f = (
            param.T_in
            + (
                param.a * q / (param.epsilon * param.rho_f * param.c_p_f)
                + param.dT_indt
            )
            * self.x
            / param.u
        )
        T_c_surf = pybamm.surf(T_c)
        T_c_surf_av = pybamm.Integral(T_c_surf, self.x) / param.Z
        self.variables.update(
            {
                "Heat transfer fluid temperature [K]": T_f,
                "Phase-change material temperature [K]": T_c,
                "X-averaged phase-change material temperature [K]": T_c_av,
                "Phase-change material surface temperature [K]": T_c_surf,
                "X-averaged phase-change material surface temperature [K]": T_c_surf_av,
                "Phase-change material enthalpy [J.m-3]": H,
                "X-averaged phase-change material enthalpy [J.m-3]": H_av,
                "Phase": (
                    H
                    >= param.rho_s * param.c_p_s * param.T_m + param.rho_s * param.L / 2
                ),
                "X-averaged phase": (
                    H_av
                    >= param.rho_s * param.c_p_s * param.T_m + param.rho_s * param.L / 2
                ),
            }
        )
