import argparse

from ATK.configuration.epoch_config import EPOCH_CONFIG

sections = list(EPOCH_CONFIG.as_dict().keys())


def handle_set(args):
    EPOCH_CONFIG._set("vizier_aliases", args.alias.lower(), args.epoch)


def handle_del(args):
    EPOCH_CONFIG._del("vizier_aliases", args.alias.lower())


def handle_reset(args):
    EPOCH_CONFIG._reset()


def handle_open(args):
    EPOCH_CONFIG._open()


def handle_show(args):
    EPOCH_CONFIG._show()


def main():
    parser = argparse.ArgumentParser(description="Facilitates the viewing and editing of the ATK epochs file.")
    subparsers = parser.add_subparsers(dest="job")

    # show()
    p_show = subparsers.add_parser("show", help="Prints the ATK epoch file to stdout.")
    p_show.set_defaults(func=handle_show)

    # reset()
    p_reset = subparsers.add_parser("reset", help="Reset the ATK epoch file to its default state.")
    p_reset.set_defaults(func=handle_reset)

    # open()
    p_open = subparsers.add_parser("open", help="Opens the ATK epoch file in the default editor.")
    p_open.set_defaults(func=handle_open)

    # set()
    p_set = subparsers.add_parser("set", help="Sets the epoch of a survey or Vizier table alias.")
    p_set.add_argument("alias", type=str, metavar="<ALIAS>", help="Survey or alias of new epoch definition.")
    p_set.add_argument(
        "epoch", type=str, metavar="<EPOCH>", help="Epoch of survey in ISOT format (YYYY-MM-DDTHH:MM:SS.SSS)."
    )
    p_set.set_defaults(func=handle_set)

    # del()
    del_parser = subparsers.add_parser("del", help="Deletes an existing epoch definition for a Vizier catalogue alias.")
    del_parser.add_argument(
        "alias", type=str, metavar="<ALIAS>", help="Name of alias for which an epoch definition should be deleted."
    )

    # parse and dispatch
    args = parser.parse_args()
    args.func(args)
