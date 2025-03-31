import argparse

from AstroToolkit.Overlays import (addOverlay, delOverlay, openOverlays,
                                   resetOverlays)


def main():
    parser = argparse.ArgumentParser(description="Facilitates the viewing and editing of the ATK overlays file.")

    sub_parsers = parser.add_subparsers(dest="job")

    add_parser = sub_parsers.add_parser("add", help="Adds an overlay definition for a given survey.")
    add_parser.add_argument("alias", type=str, nargs=1, help="Overlay alias")
    add_parser.add_argument(
        "--ra", type=str, nargs=1, required=True, metavar=("RA_NAME"), help="Name of right ascension column in Vizier"
    )
    add_parser.add_argument(
        "--dec", type=str, nargs=1, required=True, metavar=("DEC_NAME"), help="Name of declination column in Vizier"
    )
    add_parser.add_argument(
        "--id", type=str, nargs=1, required=True, metavar=("ID_NAME"), help="Name of ID column in Vizier"
    )
    add_parser.add_argument(
        "--mags",
        type=str,
        nargs="+",
        required=False,
        metavar=("MAG_NAMES"),
        help="Names of magnitude columns in Vizier. If provided, a scaled_detection overlay defintion is generated. Otherwise, a detection overlay is generated.",
    )

    del_parser = sub_parsers.add_parser(
        "del", help="Deletes an existing overlay definition of a given kind for a given survey."
    )
    del_parser.add_argument(
        "kind",
        type=str,
        nargs=1,
        choices=["scaled_detection", "detection"],
        help="Kind of overlay definition to delete",
    )
    del_parser.add_argument("alias", type=str, nargs=1, help="Overlay alias")

    sub_parsers.add_parser("reset", help="Resets the list of overlay definitions to its default state")
    sub_parsers.add_parser("open", help="Opens the list of overlay definitions in the default text editor")

    args = parser.parse_args()

    if args.job == "add":
        addOverlay(args.alias[0], args.ra[0], args.dec[0], args.id[0], args.mags)
    elif args.job == "del":
        delOverlay(args.kind[0], args.alias[0])
    elif args.job == "reset":
        resetOverlays()
    elif args.job == "open":
        openOverlays()
