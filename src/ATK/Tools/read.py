from pathlib import Path

from ..io.files.read import read_local


def read(path: str | Path):
    """
    Reads a local ATK fits file to recreate the original data structure
    """

    return read_local(path)
