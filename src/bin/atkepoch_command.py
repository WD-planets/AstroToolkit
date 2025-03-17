import argparse

from AstroToolkit.Epochs import (delEpoch, openEpochs, resetEpochs, setEpoch,
                                 showEpochs)


def main():
    parser = argparse.ArgumentParser()
    sub_parsers = parser.add_subparsers(dest="job")

    sub_parsers.add_parser("reset")
    sub_parsers.add_parser("show")
    sub_parsers.add_parser("open")

    set_parser = sub_parsers.add_parser("set")
    set_parser.add_argument("section", type=str)
    set_parser.add_argument("survey", type=str)
    set_parser.add_argument("epoch", type=str)

    del_parser = sub_parsers.add_parser("del")
    del_parser.add_argument("alias", type=str)

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
