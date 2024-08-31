#
# Full model for encapuslated LTES
#
import pybamm
from .base_LTES_model import BaseLTESModel


class FullModel(BaseLTESModel):
    def __init__(self, name="Full model"):
        super().__init__(name=name)

        param = self.param

        ######################
        # Variables
        ######################
        scale = abs(param.T_m - param.T_in)
        scale = 200
        T_f = pybamm.Variable(
            "Heat transfer fluid temperature [K]",
            domains={"primary": "pipe"},
            scale=scale,
            reference=param.T_m,
        )
        H = pybamm.Variable(
            "Phase-change material enthalpy [J.m-3]",
            domains={"primary": "capsule", "secondary": "pipe"},
            scale=param.T2H(param.T_0),
        )
        T_c = pybamm.Variable(
            "Phase-change material temperature [K]",
            domains={"primary": "capsule", "secondary": "pipe"},
            scale=scale,
            reference=param.T_m,
        )
        Q = pybamm.Variable(
            "Stored energy per unit area [J.m-2]",
            scale=param.Z * (1 - param.epsilon) * param.L * param.rho_s,
        )

        ######################
        # Governing equations
        ######################
        u_edge = pybamm.PrimaryBroadcastToEdges(param.u, ["pipe"])
        H_surf = pybamm.surf(H)
        q = -param.k(H_surf) * pybamm.BoundaryGradient(T_c, "right")
        dTfdt = (
            -pybamm.div(
                (pybamm.upwind(T_f) * u_edge)
                # - param.k_f / (param.rho_f * param.c_p_f) * pybamm.grad(T_f)
            )
            + param.a / (param.epsilon * param.rho_f * param.c_p_f) * q
        )
        dHdt = pybamm.div(param.k(H) * pybamm.grad(T_c))
        T_in = pybamm.boundary_value(T_f, "left")
        T_out = pybamm.boundary_value(T_f, "right")
        dQdt = param.epsilon * param.rho_f * param.c_p_f * param.u * (T_in - T_out)
        self.rhs = {T_f: dTfdt, H: dHdt, Q: dQdt}
        self.algebraic = {T_c: T_c - param.H2T(H)}

        self.boundary_conditions = {
            T_f: {
                "left": (param.T_in, "Dirichlet"),
                # "right": (pybamm.Scalar(0), "Neumann"),
            },
            T_c: {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (
                    -param.h / param.k(H_surf) * (pybamm.surf(T_c) - T_f),
                    "Neumann",
                ),
            },
        }

        self.initial_conditions = {
            T_f: param.T_0,
            H: param.T2H(param.T_0),
            T_c: param.T_0,
            Q: pybamm.Scalar(0),
        }

        ######################
        # Set variables
        ######################
        self._set_output_variables(T_f, T_c, H, Q)

    # @property
    # def default_var_pts(self):
    #     return {self.r: 50, self.x: 100}