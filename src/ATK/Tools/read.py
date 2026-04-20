from pathlib import Path

from ..io.files.read import read_local


def read(path: str | Path):
    """
    Reads a local ATK FITS file to recreate the original :class:`~ATK.Models.DataSet`.

    Parameters
    ----------
    path : str or :class:`~pathlib.Path`
        Path to file.

    See Also
    --------
    :class:`~ATK.Models.DataSet` : Data structure from which a local ATK fits file can be generated.
    :meth:`~ATK.Models.DataSet.store` : :class:`~ATK.Models.DataSet` method from which an ATK FITS file can be generated.
    """

    return read_local(path)
