#!/usr/bin/env python3

"""Version handler for postoga."""

import os


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__credits__ = ["Bogdan M. Kirilenko"]


class Version:
    def __init__(self, major, minor, patch, dev=False):
        self.major = major
        self.minor = minor
        self.patch = patch
        self.dev = dev
        self.color = "blue"
        self.version_repr = f"{major}.{minor}.{patch}"
        self.readme_repr = self.version_repr
        if self.dev:
            self.version_repr += f"-devel"
            self.readme_repr += f"--devel"
            self.color = "orange"

    def update_readme(self, filename="README.md"):
        with open(filename, "r") as f:
            lines = f.readlines()

        with open(filename, "w") as f:
            for line in lines:
                if "img.shields.io/badge/version-" in line:
                    line = f"![version](https://img.shields.io/badge/version-{self.readme_repr}-{self.color})\n"
                elif "## What's new" in line:
                    line = f"## What's new on version {self.version_repr}\n\n"
                f.write(line)

    def get_py_scripts(self):
        scripts = []
        for root, dirs, files in os.walk("."):
            for file in files:
                if file.endswith(".py"):
                    if file != "version.py":
                        scripts.append(os.path.join(root, file))
                elif file.endswith(".sh"):
                    scripts.append(os.path.join(root, file))
                elif file.endswith("Makefile"):
                    scripts.append(os.path.join(root, file))
        return scripts

    def check_uncomitted(self, file):
        cmd = f"git diff --exit-code {file}"
        return os.system(cmd) != 0

    def update_scripts(self):
        scripts = self.get_py_scripts()
        for script in scripts:
            if self.check_uncomitted(script) or script == "logger.py":
                with open(script, "r") as f:
                    lines = f.readlines()
                with open(script, "w") as f:
                    for line in lines:
                        if "__version__ =" in line:
                            line = f'__version__ = "{self.version_repr}"\n'
                        elif "# version:" in line:
                            line = f"# version: {self.version_repr}\n"
                        elif "VERSION=" in line:
                            line = f'VERSION="{self.version_repr}"\n'
                        f.write(line)

    def __repr__(self):
        return self.version_repr

    def to_string(self):
        return self.version_repr


__version__ = Version(0, 9, 3, dev=True)

if __name__ == "__main__":
    print(f"postoga v.{__version__}")
    __version__.update_readme()
    __version__.update_scripts()
