#
# Standard parameters for encapsulated LTES models
#
import pybamm


class EncapsulatedLTESParameters:
    """
    Standard parameters for encapsulated LTES models
    """

    def __init__(self):
        # Set parameters
        self._set_parameters()

    def _set_parameters(self):
        """Defines the dimensional parameters"""

        # HTF parameters
        self.T_0 = pybamm.Parameter("Initial temperature [K]")
        self.T_in = pybamm.FunctionParameter(
            "Inlet temperature [K]", {"Time [s]": pybamm.t}
        )
        self.dT_indt = pybamm.FunctionParameter(
            "Inlet temperature [K]",
            {"Time [s]": pybamm.t},
            diff_variable=pybamm.t,
        )
        self.u = pybamm.Parameter("Inlet velocity [m.s-1]")
        self.rho_f = pybamm.Parameter("Heat transfer fluid density [kg.m-3]")
        self.c_p_f = pybamm.Parameter(
            "Heat transfer fluid specific heat capacity [J.kg-1.K-1]"
        )
        self.k_f = pybamm.Parameter("Heat transfer fluid conductivity [W.m-1.K-1]")
        self.h = pybamm.Parameter("Heat transfer coefficient [W.m-2.K-1]")

        # PCM parameters
        self.k_s = pybamm.Parameter("Solid phase conductivity [W.m-1.K-1]")
        self.k_l = pybamm.Parameter("Liquid phase conductivity [W.m-1.K-1]")
        self.rho_s = pybamm.Parameter("Solid phase density [kg.m-3]")
        self.rho_l = pybamm.Parameter("Liquid phase density [kg.m-3]")
        self.c_p_s = pybamm.Parameter("Solid phase specific heat capacity [J.kg-1.K-1]")
        self.c_p_l = pybamm.Parameter(
            "Liquid phase specific heat capacity [J.kg-1.K-1]"
        )
        self.L = pybamm.Parameter("Latent heat [J.kg-1]")
        self.T_m = pybamm.Parameter("Melting temperature [K]")

        # Geometric parameters
        self.R = pybamm.Parameter("Capsule radius [m]")
        self.Z = pybamm.Parameter("Pipe length [m]")
        self.epsilon = pybamm.Parameter("Porosity")
        self.a = 3 * (1 - self.epsilon) / self.R

    def H2T(self, H):
        """Convert enthalpy to temperature"""
        solid = H / (self.rho_s * self.c_p_s)
        liquid = self.T_m + (H - self.rho_s * (self.c_p_s * self.T_m + self.L)) / (
            self.rho_l * self.c_p_l
        )
        return (
            pybamm.minimum(solid, self.T_m)
            + pybamm.maximum(liquid, self.T_m)
            - self.T_m
        )

    def T2H(self, T):
        """Convert temperature to enthalpy"""
        if T == self.T_m:
            raise ValueError("Enthalpy is not uniquely defined at melting temperature")
        H_l = self.rho_l * self.c_p_l * (T - self.T_m) + self.rho_s * (
            self.c_p_s * self.T_m + self.L
        )
        H_s = self.rho_s * self.c_p_s * T
        return H_l * (T > self.T_m) + H_s * (T < self.T_m)

    def k(self, H):
        """Effective conductivity as a function of the enthalpy"""
        H_s = self.rho_s * self.c_p_s * self.T_m
        H_l = self.rho_s * (self.c_p_s * self.T_m + self.L)
        k_m = (self.k_l - self.k_s) / (self.rho_s * self.L) * (H - H_s) + self.k_s
        return (
            self.k_s * (H <= H_s) + k_m * (H > H_s) * (H < H_l) + self.k_l * (H >= H_l)
        )
