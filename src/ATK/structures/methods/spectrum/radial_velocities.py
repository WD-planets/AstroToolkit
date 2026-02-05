import astropy.units as u
import numpy as np
from astropy.units import Quantity

C_KMS = 299792.458


def get_velocities(wav: Quantity, wav_ref: Quantity) -> np.ndarray:
    return (wav - wav_ref) / wav_ref * C_KMS * u.Unit("km s-1")
