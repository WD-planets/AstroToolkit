import argparse

from AstroToolkit.Tools import readdata

from .atkquery_command import jobs_ui


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", type=str)

    args = parser.parse_args()

    data = readdata(fname=args.fname)

    data_exists = False
    if data.kind in ["lightcurve"]:
        for band in data.data:
            if band["mag"]:
                data_exists = True
    else:
        if data.data:
            data_exists = True

    if data_exists:
        jobs_ui(data)
