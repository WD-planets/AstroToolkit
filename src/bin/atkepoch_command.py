import argparse

from ATK.configuration.epoch_config import EPOCH_CONFIG


def main():
    parser = argparse.ArgumentParser(description="Facilitates the viewing and editing of the ATK epochs file.")
    sub_parsers = parser.add_subparsers(dest="job")

    sub_parsers.add_parser("reset", help="Resets the epoch list to its default state")
    sub_parsers.add_parser("show", help="Outputs the epoch list to stdout")
    sub_parsers.add_parser("open", help="Opens the epoch list in the default text editor")

    set_parser = sub_parsers.add_parser("set", help="Sets the epoch of a survey or Vizier catalogue alias")
    set_parser.add_argument("section", type=str, help="Section of the epochs file in which the survey or Vizier catalogue alias is found")
    set_parser.add_argument("survey", type=str, help="Survey for which epoch should be set")
    set_parser.add_argument("epoch", type=str, help="Epoch to set (e.g. 2016,0] for Jan 2016)")

    del_parser = sub_parsers.add_parser("del", help="Deletes an existing epoch definition for a Vizier catalogue alias")
    del_parser.add_argument("section", type=str, help="Section of the epochs file in which the survey or Vizier catalogue alias is found")
    del_parser.add_argument("alias", type=str, help="Name of alias for which an epoch definition should be deleted")

    args = parser.parse_args()

    if args.job == "reset":
        EPOCH_CONFIG._reset()
    elif args.job == "show":
        EPOCH_CONFIG._show()
    elif args.job == "open":
        EPOCH_CONFIG._open()
    elif args.job == "set":
        EPOCH_CONFIG._set(args.section, args.survey, args.epoch)
    elif args.job == "del":
        EPOCH_CONFIG._del(args.section, args.alias)
