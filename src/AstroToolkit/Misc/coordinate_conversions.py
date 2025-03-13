def conv_hms_to_deg(pos):
    if "+" in pos:
        separator = "+"
    elif "-" in pos:
        separator = "-"

    pos = pos.replace(" ", "")
    pos = pos.lower()
    if pos.startswith("j"):
        pos = pos.lstrip("j")

    ra, dec = pos.split(separator)[0], pos.split(separator)[1]

    h, m, s = float(ra[0:2]), float(ra[2:4]), float(ra[4:])
    ra_deg = 15 * h + 15 * m / 60 + 15 * s / 3600

    d, m, s = float(dec[0:2]), float(dec[2:4]), float(dec[4:])
    dec_deg = d + m / 60 + s / 3600

    if separator == "-":
        dec_deg *= -1

    return [ra_deg, dec_deg]
