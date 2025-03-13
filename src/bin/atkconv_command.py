import argparse

from AstroToolkit.Tools import deg2hms, hms2deg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="+")

    args = parser.parse_args()

    if len(args.input) == 1:
        print(hms2deg(args.input[0]))
    elif len(args.input) == 2:
        pos = [float(x) for x in args.input]
        print(deg2hms(pos))
    else:
        raise ValueError("Invalid input.")
