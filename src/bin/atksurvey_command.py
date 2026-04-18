import argparse
import warnings

from ATK.configuration.survey_config import SURVEY_CONFIG


def validate_specific_args(parser, section, args):
    ignore_args = ["job", "section", "name", "func", "epoch"]

    if section == "vizier":
        pass
    else:
        additional_args = list(arg for arg in args.__dict__.keys() if arg not in ignore_args and args.__dict__[arg] is not None)
        for arg in additional_args:
            warnings.warn(f"Argument '--{arg}' is invalid in given section, and has been ignored.")
            setattr(args, arg, None)


def handle_set(parser, args):
    section = args.section.lower()

    validate_specific_args(parser, section, args)

    SURVEY_CONFIG._set(
        section,
        args.name,
        id=args.id,
        epoch=args.epoch,
        mags=args.mags,
        errors=args.errors,
        lon=args.lon,
        lat=args.lat,
        frame=args.frame,
    )


def handle_del(args):
    section = args.type.lower()

    SURVEY_CONFIG._del(section, args.alias)


def handle_reset(args):
    SURVEY_CONFIG.reset()


def handle_open(args):
    SURVEY_CONFIG.open()


def handle_show(args):
    SURVEY_CONFIG.show()


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
    p_set = subparsers.add_parser("set", help="Sets a survey definition in the ATK survey file.")
    p_set.add_argument("section", type=str, metavar="<SECTION>", help="Entry kind, e.g. 'vizier', 'lightcurve', 'image'.")
    p_set.add_argument("name", type=str, metavar="<NAME>", help="Name of survey or Vizier catalogue alias.")

    # common args
    p_set.add_argument("--epoch", help="Epoch of survey or catalogue alias.")

    # specific args
    vizier_args = p_set.add_argument_group("Required for vizier entries")
    vizier_args.add_argument("--id", help="Vizier catalogue ID")

    overlay_args = p_set.add_argument_group("Overlay parameters (only allowed in vizier entries)")
    overlay_args.add_argument("--lon", help="Name of longitudinal column in Vizier table (e.g. RA, GLON).")
    overlay_args.add_argument("--lat", help="Name of latitudinal column in Vizier table (e.g. DEC, GLAT).")
    overlay_args.add_argument("--frame", help="Frame of coordinates (e.g. icrs, galactic).")
    overlay_args.add_argument("--mags", nargs="+", metavar="MAGS", help="Names of magnitude columns in Vizier table. If omitted, overlay will be positional.")
    overlay_args.add_argument(
        "--errors",
        nargs="+",
        metavar="ERRORS",
        help="Names of magnitude error columns in Vizier table. If omitted, overaly will be positional.",
    )

    p_set.set_defaults(func=lambda args: handle_set(p_set, args))

    # del()
    p_del = subparsers.add_parser("del", help="Delete an existing survey from the ATK surveys file.")
    p_del.add_argument("survey", type=str, metavar="<ALIAS>", help="Name of survey or Vizier catalogue alias.")
    p_del.set_defaults(func=handle_del)

    # parse and dispatch
    args = parser.parse_args()
    args.func(args)
