import argparse

from AstroToolkit.Tools import readdata
from AstroToolkit.Utility import openFileDialogue

from .atkquery_command import jobs_ui


def main():
    parser = argparse.ArgumentParser(
        description="Reads a local ATK file and provides a set of jobs to perform on the returned data structure."
    )
    parser.add_argument("fname", type=str, help="File path", nargs="?")

    args = parser.parse_args()

    if not args.fname:
        args.fname = openFileDialogue()

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
