import argparse
import inspect
import re

from AstroToolkit.Configuration.baseconfig import ConfigStruct
from AstroToolkit.Tools import query

config = ConfigStruct()

newline = "\n"


def jobs_ui(data):
    while True:
        if data.kind in ["data"]:
            accepted_jobs = ["showdata", "savedata", "exit"]
            accepted_jobs_str = "showdata, savedata <filename: str, optional>, exit"
        else:
            accepted_jobs = ["showdata", "savedata", "showplot", "saveplot", "exit"]
            accepted_jobs_str = "showdata, savedata <filename: str, optional> showplot <filename: str, optional>, saveplot <filename: str, optional>, exit"

        print(f"Available Jobs: {accepted_jobs_str}{newline}")

        job = str(input("Job? "))
        job = [a for a in re.split(r"(\s|\,)", job.strip()) if a]
        job = [x for x in job if x != " " and x != ","]

        if len(job) > 1:
            fname = job[1]
            job = job[0]
        else:
            job = job[0]
            fname = None

        if job == "showdata" and fname:
            raise Exception("fname provided for showdata job.")

        if job not in accepted_jobs:
            print(f"\nInvalid job. Accepted jobs: {accepted_jobs_str}\n")
            continue

        if job == "showdata":
            data.showdata(print_methods=False)
            print()
        elif job == "savedata":
            data.savedata(fname=fname)
        elif job == "showplot":
            if data.kind == "lightcurve":
                while True:
                    plot_kind = str(input("Plot Type? "))
                    if plot_kind in ["lightcurve", "phasefold", "powspec"]:
                        data.plot(kind=plot_kind).showplot(fname=fname)
                        break
                    else:
                        print("Invalid plot type. Accepted plot types: lightcurve,phasefold,powspec")
            else:
                data.plot().showplot(fname=fname)
        elif job == "saveplot":
            if data.kind == "lightcurve":
                while True:
                    plot_kind = str(input("Plot Type? "))
                    if plot_kind in ["lightcurve", "phasefold", "powspec"]:
                        data.plot(kind=plot_kind).saveplot(fname=fname)
                        break
                    else:
                        print("Invalid plot type. Accepted plot types: lightcurve,phasefold,powspec")
            else:
                data.plot().saveplot(fname=fname)
        elif job == "exit":
            break


def main():
    config.read_config()

    params = {
        "data": ["target", "survey", "r"],
        "reddening": ["target", "survey", "r"],
        "bulkdata": ["target", "r"],
        "image": ["target", "survey", "s"],
        "lightcurve": ["target", "survey", "r", "username", "password"],
        "hrd": ["sources"],
        "sed": ["target", "r"],
        "spectrum": ["target", "survey", "r"],
    }

    all_params = ["target", "survey", "r", "s", "username", "password", "sources"]

    parser = argparse.ArgumentParser(description="Fetches data for a target from a given survey")
    sub_parsers = parser.add_subparsers(dest="kind")

    for kind in params:
        sub_parser = sub_parsers.add_parser(kind, help=f"Performs a {kind} query")

        if "survey" in params[kind]:
            survey_msg = "Target Survey"
            if kind == "data":
                survey_msg += " or Vizier catalogue alias"
            sub_parser.add_argument("survey", help="Target survey")

        if "username" in params[kind]:
            sub_parser.add_argument(
                "--username", nargs=1, type=str, help="ATLAS username (only needed in ATLAS queries)"
            )
        if "password" in params[kind]:
            sub_parser.add_argument(
                "--password", nargs=1, type=str, help="ATLAS password (only needed in ATLAS queries)"
            )

        if "target" in params[kind]:
            group = sub_parser.add_mutually_exclusive_group(required=True)
            group.add_argument("--source", nargs=1, type=int, help="Gaia DR3 Source ID")
            group.add_argument("--pos", nargs=2, type=float, help="Position in degrees", metavar=("RA", "DEC"))

        if "sources" in params[kind]:
            sub_parser.add_argument("sources", type=int, nargs="+", help="Sources to overlay")

        if "r" in params[kind]:
            sub_parser.add_argument("-r", help="Radius of search in arcseconds", type=float, metavar=("RADIUS"))

        elif "s" in params[kind]:
            sub_parser.add_argument("-s", help="Size of image in arcseconds", type=float, metavar=("SIZE"))

    args = parser.parse_args()

    if hasattr(args, "source") and isinstance(args.source, list):
        args.source = args.source[0]
    if hasattr(args, "password") and isinstance(args.password, list):
        args.password = args.password[0]
    if hasattr(args, "username") and isinstance(args.username, list):
        args.username = args.username[0]

    for arg in vars(args):
        if arg not in params[args.kind] and arg not in ["kind", "source", "pos"]:
            setattr(args, arg, None)
    if "target" not in params[args.kind]:
        args.source = None
        args.pos = None

    for param in all_params:
        if not hasattr(args, param):
            setattr(args, param, None)

    if hasattr(args, "s") and not args.s:
        args.s = config.query_image_size
    if hasattr(args, "r") and not args.r:
        args.r = getattr(config, f"query_{kind}_radius")

    data = query(
        kind=args.kind,
        survey=args.survey,
        pos=args.pos,
        source=args.source,
        sources=args.sources,
        radius=args.r,
        size=args.s,
        username=args.username,
        password=args.password,
    )

    data_exists = False
    if args.kind in ["lightcurve"]:
        for band in data.data:
            if band["mag"]:
                data_exists = True
    else:
        if data.data:
            data_exists = True

    if data_exists:
        jobs_ui(data)
