def identifier_from_pos(pos, kind="identifier"):
    import numpy as np
    from astropy import units as u
    from astropy.coordinates import Angle

    # Check for a negative declination, only used to later make sure correct signs are shown in the resulting identifier
    ra, dec = pos[0], pos[1]
    if dec < 0:
        negativeDec = True
    else:
        negativeDec = False

    ra = Angle(ra, u.degree)
    dec = Angle(dec, u.degree)

    # Do unit conversion from deg --> hms/dms
    ra = ra.hms
    dec = dec.dms

    # Create ra_arr and dec_arr containing [H,M,S] and [D,M,S]
    ra_arr = np.array([0, 0, 0], dtype=float)
    dec_arr = np.array([0, 0, 0], dtype=float)

    ra_arr = [x for x in ra]
    dec_arr = [x for x in dec]

    def format_arr(arr):
        arr = [abs(x) for x in arr]
        remainders = [int(x) - x for x in arr]
        arr = [str(int(x)).zfill(2) for x in arr]

        return arr, remainders[-1]

    ra_str_arr, ra_remainder = format_arr(ra_arr)
    dec_str_arr, dec_remainder = format_arr(dec_arr)

    if kind == "identifier":
        # Format remainder: force 2 decimal places, round to 2 decimal places and remove '0.'
        ra_str_arr[2] += str("{:.2f}".format(round(ra_remainder, 2))[2:])
        dec_str_arr[2] += str("{:.2f}".format(round(dec_remainder, 2))[2:])
        prefix = "J"
    elif kind == "conversion":
        ra_str_arr[2] += str(ra_remainder)[2:]
        dec_str_arr[2] += str(dec_remainder)[2:]
        prefix = ""
    else:
        raise Exception("Invalid type.")

    # Writes final identifier and takes account of negative decs.
    if negativeDec:
        objRef = f"{prefix}{ra_str_arr[0]}{ra_str_arr[1]}{ra_str_arr[2]}-{dec_str_arr[0]}{dec_str_arr[1]}{dec_str_arr[2]}"
    else:
        objRef = f"{prefix}{ra_str_arr[0]}{ra_str_arr[1]}{ra_str_arr[2]}+{dec_str_arr[0]}{dec_str_arr[1]}{dec_str_arr[2]}"

    return objRef
