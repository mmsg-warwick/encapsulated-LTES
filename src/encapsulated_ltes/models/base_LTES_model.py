#
# Base model for encapuslated LTES
#
import numpy as np
import pybamm

from ..parameter_values import get_parameter_values
from ..parameters import EncapsulatedLTESParameters


class BaseLTESModel(pybamm.models.base_model.BaseModel):
    def __init__(self, name="Unnamed LTES model"):
        super().__init__(name=name)
        ######################
        # Get parameters
        ######################
        self.param = EncapsulatedLTESParameters()

        ######################
        # Define spatial variables
        ######################
        self.r = pybamm.SpatialVariable(
            "r",
            domains={"primary": "capsule", "secondary": "pipe"},
            coord_sys="spherical polar",
        )
        self.x = pybamm.SpatialVariable(
            "x", domains={"primary": "pipe"}, coord_sys="cartesian"
        )

        ######################
        # (Some) variables
        ######################
        self.variables = {
            "Time [s]": pybamm.t,
            "Time [min]": pybamm.t / 60,
            "Time [h]": pybamm.t / 3600,
            "r [m]": self.r,
            "r [mm]": self.r * 1000,
            "x [m]": self.x,
            "x [mm]": self.x * 1000,
        }

    def _set_output_variables(self, T_f, T_c, H, Q):
        param = self.param
        T_c_av = pybamm.Integral(T_c, self.x) / param.Z
        H_av = pybamm.Integral(H, self.x) / param.Z
        H_surf = pybamm.surf(H)
        q = -param.k(H_surf) * pybamm.BoundaryGradient(T_c, "right")
        q_av = pybamm.Integral(q, self.x) / param.Z
        T_c_surf = pybamm.surf(T_c)
        T_c_surf_av = pybamm.Integral(T_c_surf, self.x) / param.Z
        H_f = param.epsilon * pybamm.Integral(param.rho_f * param.c_p_f * T_f, self.x)
        H_f0 = param.epsilon * param.rho_f * param.c_p_f * param.T_0 * param.Z
        H_PCM = (1 - param.epsilon) * pybamm.Integral(pybamm.Integral(H, self.r) * 3 / (param.R ** 3 * 4 * np.pi), self.x)
        H_PCM0 = (1 - param.epsilon) * param.T2H(param.T_0) * param.Z
        H_tot = H_f + H_PCM
        H_tot0 = H_f0 + H_PCM0
        T_in = pybamm.boundary_value(T_f, "left")
        T_out = pybamm.boundary_value(T_f, "right")

        phase = (param.rho_s * param.c_p_s * param.T_m + param.rho_s * param.L / 2 <= H)
        phase_av = (H_av >= param.rho_s * param.c_p_s * param.T_m + param.rho_s * param.L / 2)

        ones_xr = pybamm.FullBroadcast(pybamm.Scalar(1), broadcast_domains={"primary": "capsule", "secondary": "pipe"})
        ones_r = pybamm.FullBroadcast(pybamm.Scalar(1), broadcast_domains={"primary": "capsule"})
        SoC = pybamm.Integral(phase, self.r) / pybamm.Integral(ones_xr, self.r)
        SoC_av = pybamm.Integral(phase_av, self.r) / pybamm.Integral(ones_r, self.r)

        self.variables.update(
            {
                "Heat transfer fluid temperature [K]": T_f,
                "Heat transfer fluid temperature [degC]": T_f - pybamm.Scalar(273.15),
                "Phase-change material temperature [K]": T_c,
                "Phase-change material temperature [degC]": T_c - pybamm.Scalar(273.15),
                "X-averaged phase-change material temperature [K]": T_c_av,
                "X-averaged phase-change material temperature [degC]": T_c_av - pybamm.Scalar(273.15),
                "Phase-change material surface temperature [K]": T_c_surf,
                "Phase-change material surface temperature [degC]": T_c_surf - pybamm.Scalar(273.15),
                "X-averaged phase-change material surface temperature [K]": T_c_surf_av,
                "X-averaged phase-change material surface temperature [degC]": T_c_surf_av - pybamm.Scalar(273.15),
                "Phase-change material enthalpy [J.m-3]": H,
                "X-averaged phase-change material enthalpy [J.m-3]": H_av,
                "Phase": phase,
                "X-averaged phase": phase_av,
                "State of charge": SoC,
                "X-averaged state of charge": SoC_av,
                "Inlet temperature [K]": T_in,
                "Outlet temperature [K]": T_out,
                "Inlet temperature [degC]": T_in - pybamm.Scalar(273.15),
                "Outlet temperature [degC]": T_out - pybamm.Scalar(273.15),
                "Total enthalpy of phase-change material per unit area [J.m-2]": H_PCM,
                "Total enthalpy of heat transfer fluid per unit area [J.m-2]": H_f,
                "Total enthalpy per unit area [J.m-2]": H_tot,
                "Variation in total enthalpy per unit area [J.m-2]": H_tot - H_tot0,
                "Stored energy per unit area [J.m-2]": Q,
                "Flux into phase-change material [W.m-2]": q,
                "X-averaged flux into phase-change material [W.m-2]": q_av,
                "Error in energy conservation [J.m-2]": H_tot - H_tot0 - Q,
                "Relative error in energy conservation [%]": (H_tot - H_tot0 - Q) / H_tot * 100,
            }
        )

    @property
    def default_geometry(self):
        return pybamm.Geometry(
            {
                "capsule": {self.r: {"min": pybamm.Scalar(0), "max": self.param.R}},
                "pipe": {self.x: {"min": pybamm.Scalar(0), "max": self.param.Z}},
            }
        )

    @property
    def default_submesh_types(self):
        return {
            "capsule": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh),
            "pipe": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh),
        }

    @property
    def default_var_pts(self):
        return {self.r: 40, self.x: 80}

    @property
    def default_spatial_methods(self):
        return {
            "capsule": pybamm.FiniteVolume(),
            "pipe": pybamm.FiniteVolume(),
        }

    @property
    def default_solver(self):
        return pybamm.IDAKLUSolver()

    @property
    def default_quick_plot_variables(self):
        return [
            "Heat transfer fluid temperature [K]",
            "Phase-change material temperature [K]",
            "Phase-change material enthalpy [J.m-3]",
            "Phase",
            "Stored energy per unit area [J.m-2]",
        ]

    @property
    def default_parameter_values(self):
        return get_parameter_values()
