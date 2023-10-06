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
log.write("This is a message.")
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


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"


import os
import datetime
from version import __version__
from constants import Constants
import subprocess


class Log:
    """Logger class for postoga."""

    def __init__(self, path, name):
        self.path = path
        self.name = name
        self.log_file = os.path.join(path, name)
        self.version = __version__
        self.commit = self.shell(Constants.Commands.COMMIT)
        self.branch = self.shell(Constants.Commands.BRANCH)

    def start(self):
        start_message = f"{'#'*36}\npostoga: the post-TOGA processing pipeline"
        version = f"version: {self.version}"
        commit = f"commit: {self.commit}"
        branch = f"branch: {self.branch}"
        metadata = "\n".join([version, commit, branch])
        with open(self.log_file, "w") as log:
            log.write(f"{start_message}\n\n{metadata}\n\n")

    def record(self, message, timestamp=True):
        formatted_message = message
        if timestamp:
            timestamp = datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
            formatted_message = f"{timestamp} - {message}"

        with open(self.log_file, "a") as log:
            log.write(f"{formatted_message}\n")

    def close(self):
        end_message = f"postoga finished!\n{'#'*36}"
        with open(self.log_file, "a") as log:
            log.write(f"\n{end_message}\n")

    def shell(self, cmd):
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.stdout.strip()


if __name__ == "__main__":
    log = Log(".", "test")
    log.start()
    log.write("This is a message.")
    log.close()
