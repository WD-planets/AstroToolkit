import argparse
import re

from AstroToolkit.Configuration.baseconfig import ConfigStruct
from AstroToolkit.Tools import query

config = ConfigStruct()


def jobs_ui(data):
    while True:
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

        if data.kind in ["data"]:
            accepted_jobs = ["showdata", "savedata", "exit"]
            accepted_jobs_str = "showdata, savedata <filename (optional)>, exit"
        else:
            accepted_jobs = ["showdata", "savedata", "showplot", "saveplot", "exit"]
            accepted_jobs_str = "showdata, savedata <filename (optional)> showplot <filename (optional)>, saveplot<filename (optional)>, exit"

        if job not in accepted_jobs:
            print(f"\nInvalid job. Accepted jobs: {accepted_jobs_str}\n")
            continue

        if job == "showdata":
            data.showdata()
            print()
        elif job == "savedata":
            data.savedata(fname=fname)
        elif job == "showplot":
            if data.kind == "lightcurve":
                while True:
                    plot_kind = str(input("Plot Type? "))
                    if plot_kind in [
                        "lightcurve",
                        "phasefold",
                        "powspec",
                        "phase",
                        "fold",
                    ]:
                        data.plot(kind=plot_kind).showplot(fname=fname)
                        break
                    else:
                        print(
                            "Invalid plot type. Accepted plot types: lightcurve,phasefold,powspec"
                        )
            else:
                data.plot().showplot(fname=fname)
        elif job == "saveplot":
            if data.kind == "lightcurve":
                while True:
                    plot_kind = str(input("Plot Type? "))
                    if plot_kind in [
                        "lightcurve",
                        "phasefold",
                        "powspec",
                        "phase",
                        "fold",
                    ]:
                        data.plot(kind=plot_kind).saveplot()
                        break
                    else:
                        print(
                            "Invalid plot type. Accepted plot types: lightcurve,phasefold,powspec"
                        )
            else:
                data.plot().saveplot(fname=fname)
        elif job == "exit":
            break


def main():
    config.read_config()

    parser = argparse.ArgumentParser()
    parser.add_argument("kind", type=str)
    parser.add_argument("survey", type=str)
    parser.add_argument("target", nargs="+")

    args = parser.parse_args()

    if len(args.target) > 2:
        pos = [float(args.target[0]), float(args.target[1])]
        source = None
        radius = float(args.target[2])

    elif len(args.target) > 1:
        if float(args.target[0]) > pow(10, 10):
            source = args.target[0]
            radius = args.target[1]
            pos = None
        else:
            pos = [float(args.target[0]), float(args.target[1])]
            source = None
            radius = None
    else:
        source = int(args.target[0])
        pos = None
        radius = None

    if not radius:
        if args.kind not in ["image"]:
            radius = getattr(config, f"query_{args.kind}_radius")
        else:
            size = config.query_image_size
            radius = None

    if args.kind != "image":
        size = None

    if args.survey == "atlas":
        atlas_login = str(input("ATLAS username and password: "))
        atlas_login = [a for a in re.split(r"(\s|\,)", atlas_login.strip()) if a]
        atlas_login = [x for x in atlas_login if x != " " and x != ","]

        if len(atlas_login) < 2:
            raise Exception("Only one argument for ATLAS login provided.")
        else:
            atlas_username, atlas_password = atlas_login[0], atlas_login[1]
    else:
        atlas_username, atlas_password = None, None

    data = query(
        kind=args.kind,
        survey=args.survey,
        pos=pos,
        source=source,
        radius=radius,
        username=atlas_username,
        password=atlas_password,
        size=size,
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
