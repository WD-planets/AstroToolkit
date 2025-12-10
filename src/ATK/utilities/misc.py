import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout


@contextmanager
def suppress_stdout() -> None:
    """
    Suppress stdout temporarily, should be used inside with block
    """

    with open(os.devnull, "w") as fnull, redirect_stdout(fnull), redirect_stderr(fnull):
        yield
