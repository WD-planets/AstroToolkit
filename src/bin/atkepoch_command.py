import argparse

from AstroToolkit.Epochs import (delEpoch, openEpochs, resetEpochs, setEpoch,
                                 showEpochs)


def main():
    parser = argparse.ArgumentParser(description="Facilitates the viewing and editing of the ATK epochs file.")
    sub_parsers = parser.add_subparsers(dest="job")

    sub_parsers.add_parser("reset", help="Resets the epoch list to its default state")
    sub_parsers.add_parser("show", help="Outputs the epoch list to stdout")
    sub_parsers.add_parser("open", help="Opens the epoch list in the default text editor")

    set_parser = sub_parsers.add_parser("set", help="Sets the epoch of a survey or Vizier catalogue alias")
    set_parser.add_argument(
        "section", type=str, help="Section of the epochs file in which survey or Vizier catalogue alias is found"
    )
    set_parser.add_argument("survey", type=str, help="Survey for which epoch should be set")
    set_parser.add_argument("epoch", type=str, help="Epoch to set (e.g. 2016,0] for Jan 2016)")

    del_parser = sub_parsers.add_parser("del", help="Deletes an existing epoch definition for a Vizier catalogue alias")
    del_parser.add_argument("alias", type=str, help="Name of alias for which an epoch definition should be deleted")

    args = parser.parse_args()

    if args.job == "reset":
        resetEpochs()
    elif args.job == "show":
        showEpochs()
    elif args.job == "open":
        openEpochs()
    elif args.job == "set":
        epoch = args.epoch.split(",")
        epoch = [int(x) for x in epoch]
        setEpoch(section=args.section, survey=args.survey, epoch=epoch)
    elif args.job == "del":
        delEpoch(alias=args.alias)
