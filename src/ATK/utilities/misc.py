import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from importlib.metadata import version


@contextmanager
def suppress_stdout() -> None:
    """
    Suppress stdout temporarily, should be used inside with block
    """

    with open(os.devnull, "w") as fnull, redirect_stdout(fnull), redirect_stderr(fnull):
        yield


def get_package_version():
    return version("AstroToolkit")


def angle_to_quantity(angle, unit):
    return angle.to_value(unit) * unit
