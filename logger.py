#!/usr/env/bin python3


"""
Logger class for postoga.

Manages the logging of postoga, recording
all the steps and errors that occur during
the execution of the tool.

It is recommended to use the `Log` class
in the following way:

``` python
from logger import Log

log = Log("path/to/log", "name")
log.start()
log.record("This is a message.")
log.close()
```

The above code will create a log file in
`path/to/log` named `name.log` with the
following contents:

```
    Log file for name.

    This is a message.
```

This class have been included in postoga
to better manage the logging of the tool.
"""


import os
import datetime
from version import __version__
from constants import Constants
from modules.utils import shell
import logging


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.4.0-devel"


class Log:
    """Logger class for postoga."""

    def __init__(self, path: str, log_file: str):
        self.log_file = os.path.join(path, log_file)
        self.version = __version__
        self.commit = shell(Constants.Commands.COMMIT)
        self.branch = shell(Constants.Commands.BRANCH)

    def start(self):
        logging.basicConfig(
            filename=self.log_file,
            level=logging.INFO,
            format="[%(asctime)s] - %(levelname)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    def intro(self):
        start_message = f"{'#'*36}\npostoga: the post-TOGA processing pipeline"
        version = f"version: {self.version}"
        commit = f"commit: {self.commit}"
        branch = f"branch: {self.branch}\n\n"
        metadata = "\n".join([version, commit, branch])

        with open(self.log_file, "w") as log:
            log.write("\n".join([start_message, metadata]))

    def record(self, message, timestamp=True):
        logging.info(message)

    @classmethod
    def connect(cls, path, log_file):
        log = cls(path, log_file)
        log.start()
        return log

    def close(self):
        end_message = f"postoga finished!\n{'#'*36}"
        logging.info(end_message)


if __name__ == "__main__":
    # To create a new log file
    log = Log("path/to/log", "name")
    log.start()
    log.record("This is a message.")
    log.record("Another message.")

    # To use an existing log file
    log = Log.connect("my_existing_log.log")
    log.record("Log this message in the existing log.")
