#!/usr/bin/env python3

"""Version handler for postoga."""


__author__ = "Alejandro Gonzales-Irribarren"
__credits__ = "Bogdan M. Kirilenko"


class Version:
    def __init__(self, major, minor, patch, dev=False):
        self.major = major
        self.minor = minor
        self.patch = patch
        self.dev = dev
        self.color = "blue"
        self.version_repr = f"{major}.{minor}.{patch}"
        if self.dev:
            self.version_repr += f"--devel"
            self.color = "orange"

    def update_readme(self, filename="README.md"):
        with open(filename, "r") as f:
            lines = f.readlines()

        with open(filename, "w") as f:
            for line in lines:
                if "img.shields.io/badge/version-" in line:
                    line = f"![version](https://img.shields.io/badge/version-{self.version_repr}-{self.color})\n"
                f.write(line)

    def __repr__(self):
        return self.version_repr

    def to_string(self):
        return self.version_repr


__version__ = Version(0, 2, 0, dev=True)

if __name__ == "__main__":
    print(f"postoga v.{__version__}")
    __version__.update_readme()
