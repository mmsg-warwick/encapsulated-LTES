# Encapsulated LTES models

[![Tests](https://github.com/mmsg-warwick/encapsulated-LTES/actions/workflows/periodic_tests.yml/badge.svg?branch=main)](https://github.com/mmsg-warwick/encapsulated-LTES/actions/workflows/periodic_tests.yml)

<!-- [![Documentation Status][rtd-badge]][rtd-link]

[![PyPI version][pypi-version]][pypi-link]
[![Conda-Forge][conda-badge]][conda-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

[![GitHub Discussions][github-discussions-badge]][github-discussions-link] -->

<!-- SPHINX-START -->

<!-- prettier-ignore-start -->

[actions-badge]:            https://github.com/mmsg-warwick/supercapacitors/workflows/CI/badge.svg
[actions-link]:             https://github.com/mmsg-warwick/supercapacitors/actions
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/supercapacitors
[conda-link]:               https://github.com/conda-forge/supercapacitors-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/mmsg-warwick/supercapacitors/discussions
[pypi-link]:                https://pypi.org/project/supercapacitors/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/supercapacitors
[pypi-version]:             https://img.shields.io/pypi/v/supercapacitors
[rtd-badge]:                https://readthedocs.org/projects/supercapacitors/badge/?version=latest
[rtd-link]:                 https://supercapacitors.readthedocs.io/en/latest/?badge=latest

<!-- prettier-ignore-end -->

This repository contains the models for encapsulated latent energy storage (LTES) from the article:
> F. Brosa Planella, A simple model for latent thermal energy storage systems with encapsulated phase-change material, Submitted for publication (2024).

## 🚀 Installing the package
The package is not yet available on PyPI so it needs to be installed from the source code. These instructions assume that you have a compatible Python version installed (between 3.9 and 3.12).

### Linux and macOS
First clone the repository, either from the command line or using a Git client:

```bash
git clone git@github.com:mmsg-warwick/encapsulated-LTES.git
```

If you do not have nox installed, install it with

```bash
python3 -m pip install nox
```

Then, navigate to the repository you just cloned and run

```bash
nox -s dev
```

This will create a virtual environment called `venv` in your current directory and install the package in editable mode with all the development dependencies. To activate the virtual environment, run

```bash
source env/bin/activate
```

You can now run the examples in the `examples` directory.

If needed, you can deactivate your virtual environment with

```bash
deactivate
```

### Windows
First clone the repository, either from the command line or using a Git client:

```bash
git clone git@github.com:mmsg-warwick/encapsulated-LTES.git
```

If you do not have nox installed, install it with

```bash
python3 -m pip install nox
```

Then, navigate to the repository you just cloned and run

```bash
nox -s dev
```

This will create a virtual environment called `venv` in your current directory and install the package in editable mode with all the development dependencies. To activate the virtual environment, run

```bash
venv\Scripts\activate.bat
```
if you are using Command Prompt, or
```bash
venv\Scripts\Activate.ps1
```
if you are using PowerShell.


You can now run the examples in the `examples` directory.

If needed, you can deactivate your virtual environment with

```bash
deactivate
```
