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
__version__ = "0.10.1"


def main():
    try:
        from .run import TogaDir, parse_args
    except ImportError:
        from run import TogaDir, parse_args

    args = parse_args()
    TogaDir(args).run()


if __name__ == "__main__":
    main()
