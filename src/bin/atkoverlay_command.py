import argparse

from ATK.configuration.alias_config import ALIAS_CONFIG
from ATK.configuration.overlay_config import OVERLAY_CONFIG

sections = ["photometric", "positional"]
aliases = list(ALIAS_CONFIG.as_dict()["vizier_aliases"].keys())


def validate_common_args(parser, section, args):
    """
    Ensure lat/lon/ID/frame args existfor all overlay sections
    """

    missing = [f"--{name}" for name in ("lat", "lon", "id", "frame") if getattr(args, name) is None]

    if missing:
        parser.error(f"For section '{section}', the following arguments are required: {', '.join(missing)}")


def validate_specific_args(parser, section, args):
    if section == "photometric":
        if args.mags is None or args.errors is None:
            parser.error("Photometric overlays require both --mags and --errors.")
    elif section == "positional":
        if args.mags or args.errors:
            parser.error("Positional overlays do not use --mags or --errors.")
    else:
        parser.error(f"Unexpected section '{section}'. Section should be one of 'photometric', 'positional'.")


def handle_set(parser, args):
    section = args.type.lower()

    validate_common_args(parser, section, args)
    validate_specific_args(parser, section, args)

    OVERLAY_CONFIG._set(
        section,
        args.alias,
        lon_column=args.lon,
        lat_column=args.lat,
        frame=args.frame,
        id_column=args.id,
        mags=args.mags,
        errors=args.errors,
    )


def handle_del(args):
    section = args.type.lower()

    OVERLAY_CONFIG._del(section, args.alias)


def handle_reset(args):
    OVERLAY_CONFIG._reset()


def handle_open(args):
    OVERLAY_CONFIG._open()


def handle_show(args):
    OVERLAY_CONFIG._show()


def main():
    parser = argparse.ArgumentParser(description="Facilitates viewing and editing of the ATK overlays file.")
    subparsers = parser.add_subparsers(dest="job", required=True)

    # show()
    p_show = subparsers.add_parser("show", help="Prints the ATK overlays file to stdout.")
    p_show.set_defaults(func=handle_show)

    # reset()
    p_reset = subparsers.add_parser("reset", help="Reset the ATK overlays file to its default state.")
    p_reset.set_defaults(func=handle_reset)

    # open()
    p_open = subparsers.add_parser("open", help="Opens the ATK overlays file in the default editor.")
    p_open.set_defaults(func=handle_open)

    # set()
    p_set = subparsers.add_parser("set", help="Sets an overlay definition in the ATK overlays file.")
    p_set.add_argument(
        "type",
        metavar="<TYPE>",
        help=f"Overlay type, from: {', '.join(sections)}. Photometric overlays use magnitudes to scale markers.",
        choices=sections,
    )
    p_set.add_argument(
        "alias",
        type=str,
        metavar="<ALIAS>",
        help=f"Corresponding alias in the ATK alias file, from: {', '.join(aliases)}.",
        choices=aliases,
    )

    # common args
    p_set.add_argument("--lon", help="Name of longitudinal column in Vizier table (e.g. RA, GLON).", required=True)
    p_set.add_argument("--lat", help="Name of latitudinal column in Vizier table (e.g. DEC, GLAT).", required=True)
    p_set.add_argument("--id", help="Survey-specific ID column in Vizier table.", required=True)
    p_set.add_argument("--frame", help="Frame of coordinates (e.g. icrs, galactic).", required=True)

    # photometric only
    list_group = p_set.add_argument_group("Only required for photometric overlays")
    list_group.add_argument("--mags", nargs="+", metavar="MAGS", help="Names of magnitude columns in Vizier table.")
    list_group.add_argument(
        "--errors", nargs="+", metavar="ERRORS", help="Names of magnitude error columns in Vizier table."
    )

    p_set.set_defaults(func=lambda args: handle_set(p_set, args))

    # del()
    p_del = subparsers.add_parser("del", help="Delete an existing overlay definition in the ATK overlays file.")
    p_del.add_argument(
        "type",
        metavar="<TYPE>",
        help=f"Overlay type, from: {', '.join(sections)}. Photometric overlays use magnitudes to scale markers.",
        choices=sections,
    )
    p_del.add_argument(
        "alias",
        type=str,
        metavar="<ALIAS>",
        help="Name of survey or alias for which an overlay definition should be deleted.",
    )
    p_del.set_defaults(func=handle_del)

    # parse and dispatch
    args = parser.parse_args()
    args.func(args)
