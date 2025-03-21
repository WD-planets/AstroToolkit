import argparse
from argparse import RawTextHelpFormatter

from AstroToolkit.Tools import deg2hms, hms2deg


def main():
    parser = argparse.ArgumentParser(
        description="Converts a position in degrees to a position in HHMMSS.SS±DDMMSS.SS format and vice versa.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "Position",
        nargs="+",
        help="2 Arguments:\n  [RA] [DEC] in degrees\n1 Argument:\n  [POSITION] in HHMMSS.SS±DDMMSS.SS format",
    )

    args = parser.parse_args()

    if len(args.input) == 1:
        print(hms2deg(args.input[0]))
    elif len(args.input) == 2:
        pos = [float(x) for x in args.input]
        print(deg2hms(pos))
    else:
        raise ValueError("Invalid input.")
