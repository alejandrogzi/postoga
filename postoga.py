#!/usr/bin/env python3

"""
postoga.py

The post-processing TOGA pipeline.
At its core, this tool takes a TOGA results directory and produces
a series of steps to reduce the amount of manual work required
to pre-process files for downstream analysis.
"""

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def main():
    from run import TogaDir, parser

    args = parser()
    TogaDir(args).run()


if __name__ == "__main__":
    main()
