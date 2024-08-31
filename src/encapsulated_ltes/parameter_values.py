#
# Parameter values for encapsulated LTES models
#
import pybamm


def get_parameter_values(parameter_set="Raul2018"):
    """
    Returns the parameters for the encapsulated LTES models

    Parameters
    ----------
    parameter_set : str, optional
        The name of the parameter set to use (default is "Raul2018")

    Returns
    -------
    parameters : parameter_values.ParameterValues
        The parameters for the encapsulated LTES model
    """
    if parameter_set == "Raul2018 enthalpy":
        parameter_values = {
            "Solid phase conductivity [W.m-1.K-1]": 0.45,
            "Liquid phase conductivity [W.m-1.K-1]": 0.45,
            "Solid phase density [kg.m-3]": 1500,
            "Liquid phase density [kg.m-3]": 1500,
            "Solid phase specific heat capacity [J.kg-1.K-1]": 2013.3,
            "Liquid phase specific heat capacity [J.kg-1.K-1]": 2013.3,
            "Latent heat [J.kg-1]": 2.497e5,
            "Melting temperature [K]": 441.85,
            "Boundary temperature [K]": 453.15,
            "Radius [m]": 1.55e-2,
            "Initial enthalpy [J.m-3]": 1500 * 2013.3 * 400,
        }
    elif parameter_set == "Raul2018":
        parameter_values = {
            "Solid phase conductivity [W.m-1.K-1]": 0.45,
            "Liquid phase conductivity [W.m-1.K-1]": 0.45,
            "Solid phase density [kg.m-3]": 1500,
            "Liquid phase density [kg.m-3]": 1500,
            "Solid phase specific heat capacity [J.kg-1.K-1]": 2013.3,
            "Liquid phase specific heat capacity [J.kg-1.K-1]": 2013.3,
            "Latent heat [J.kg-1]": 2.497e5,
            "Melting temperature [K]": 441.85,
            "Initial temperature [K]": 400,
            "Inlet temperature [K]": 453.15,
            "Capsule radius [m]": 1.55e-2,
            "Pipe length [m]": 0.36,
            "Inlet velocity [m.s-1]": 1.0e-2,
            "Porosity": 0.6,
            "Heat transfer coefficient [W.m-2.K-1]": 100,
            "Heat transfer fluid density [kg.m-3]": 720.9,
            "Heat transfer fluid specific heat capacity [J.kg-1.K-1]": 3097.4,
            "Heat transfer fluid conductivity [W.m-1.K-1]": 0.116,
        }
    elif parameter_set == "Nallusamy2007":
        parameter_values = {
            "Solid phase conductivity [W.m-1.K-1]": 0.4,
            "Liquid phase conductivity [W.m-1.K-1]": 0.15,
            "Solid phase density [kg.m-3]": 861,
            "Liquid phase density [kg.m-3]": 778,
            "Solid phase specific heat capacity [J.kg-1.K-1]": 1850,
            "Liquid phase specific heat capacity [J.kg-1.K-1]": 2384,
            "Latent heat [J.kg-1]": 213000,
            "Melting temperature [K]": 333.15,
            "Initial temperature [K]": 305.15,
            "Inlet temperature [K]": 343.15,
            "Capsule radius [m]": 27.5e-3,
            "Pipe length [m]": 0.46,
            "Inlet velocity [m.s-1]": 6.5e-4,
            "Porosity": 0.5,
            "Heat transfer coefficient [W.m-2.K-1]": 100,
            "Heat transfer fluid density [kg.m-3]": 1000,
            "Heat transfer fluid specific heat capacity [J.kg-1.K-1]": 4186,
            "Heat transfer fluid conductivity [W.m-1.K-1]": 0.6,
        }
    else:
        msg = f"Parameter set '{parameter_set}' not recognised"
        raise ValueError(msg)

    return pybamm.ParameterValues(parameter_values)
