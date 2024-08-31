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
        scale = 200
        T_f = pybamm.Variable(
            "Heat transfer fluid temperature [K]",
            domains={"primary": "pipe"},
            scale=scale,
            reference=param.T_m,
        )
        H_av = pybamm.Variable(
            "X-averaged phase-change material enthalpy [J.m-3]",
            domains={"primary": "capsule"},
            scale=param.T2H(param.T_0),
        )
        T_c_av = pybamm.Variable(
            "X-averaged phase-change material temperature [K]",
            domains={"primary": "capsule"},
            scale=scale,
            reference=param.T_m,
        )
        H = pybamm.SecondaryBroadcast(H_av, "pipe")
        T_c = pybamm.SecondaryBroadcast(T_c_av, "pipe")
        Q = pybamm.Variable(
            "Stored energy per unit area [J.m-2]",
            scale=param.Z * (1 - param.epsilon) * param.L * param.rho_s,
        )

        ######################
        # Governing equations
        ######################
        u_edge = pybamm.PrimaryBroadcastToEdges(param.u, ["pipe"])
        H_surf = pybamm.surf(H)
        H_av_surf = pybamm.surf(H_av)
        q = -param.k(H_surf) * pybamm.BoundaryGradient(T_c, "right")
        dTfdt = (
            -pybamm.div(
                u_edge * pybamm.upwind(T_f)
                # - param.k_f / (param.rho_f * param.c_p_f) * pybamm.grad(T_f)
            )
            + param.a / (param.epsilon * param.rho_f * param.c_p_f) * q
        )
        dHavdt = pybamm.div(param.k(H_av) * pybamm.grad(T_c_av))
        T_in = pybamm.boundary_value(T_f, "left")
        T_out = pybamm.boundary_value(T_f, "right")
        dQdt = param.epsilon * param.rho_f * param.c_p_f * param.u * (T_in - T_out)
        self.rhs = {T_f: dTfdt, H_av: dHavdt, Q: dQdt}
        self.algebraic = {T_c_av: T_c_av - param.H2T(H_av)}

        T_f_av = pybamm.Integral(T_f, self.x) / param.Z

        self.boundary_conditions = {
            T_f: {
                "left": (param.T_in, "Dirichlet"),
                # "right": (pybamm.Scalar(0), "Neumann"),
            },
            T_c_av: {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (
                    -param.h / param.k(H_av_surf) * (pybamm.surf(T_c_av) - T_f_av),
                    "Neumann",
                ),
            },
        }

        self.initial_conditions = {
            T_f: param.T_0,
            H_av: param.T2H(param.T_0),
            T_c_av: param.T_0,
            Q: pybamm.Scalar(0),
        }

        ######################
        # Set variables
        ######################
        self._set_output_variables(T_f, T_c, H, Q)
