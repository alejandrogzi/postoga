#!/usr/bin/env python3


""" A moduule to install necessary packages. """

import subprocess


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.4.0-devel"


packages = ["pandas==2.0.2", "numpy==1.24.3"]

for package in packages:
    try:
        subprocess.check_call(["pip", "install", package])
        print(f"Installed {package} successfully.")
    except subprocess.CalledProcessError:
        print(f"Failed to install {package}.")
