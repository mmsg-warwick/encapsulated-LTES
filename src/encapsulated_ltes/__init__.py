"""
Copyright (c) 2024 Ferran Brosa Planella. All rights reserved.

encapsulated-LTES: A project to simulate encapsulated latent thermal energy
"""
__version__ = "0.1.0"

import pybamm

from .models import *
from .parameter_values import get_parameter_values
from .parameters import EncapsulatedLTESParameters
from .utils import get_interface_position, root_dir
from .plot import *

__all__ = [
    "__version__",
    "pybamm",
    "parameter_sets",
    "Model",
    "models",
]