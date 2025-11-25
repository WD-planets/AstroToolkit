import argparse

from ATK.configuration.alias_config import ALIAS_CONFIG


def handle_set(args):
    ALIAS_CONFIG._set("vizier_aliases", args.alias.lower(), args.table_id)


def handle_del(args):
    ALIAS_CONFIG._del("vizier_aliases", args.alias.lower())


def handle_reset(args):
    ALIAS_CONFIG._reset()


def handle_open(args):
    ALIAS_CONFIG._open()


def handle_show(args):
    ALIAS_CONFIG._show()


def main():
    parser = argparse.ArgumentParser(description="Facilitates viewing and editing of the ATK aliases file.")
    subparsers = parser.add_subparsers(dest="job", required=True)

    # show()
    p_show = subparsers.add_parser("show", help="Prints the ATK alias file to stdout.")
    p_show.set_defaults(func=handle_show)

    # reset()
    p_reset = subparsers.add_parser("reset", help="Reset the ATK alias file to its default state.")
    p_reset.set_defaults(func=handle_reset)

    # open()
    p_open = subparsers.add_parser("open", help="Opens the ATK alias file in the default editor.")
    p_open.set_defaults(func=handle_open)

    # set()
    p_set = subparsers.add_parser("set", help="Sets an alias to a Vizier table in the ATK alias file.")
    p_set.add_argument("alias", type=str, metavar="<ALIAS>", help="Name of alias to a Vizier table.")
    p_set.add_argument("table_id", type=str, metavar="<TABLE_ID>", help="Vizier table ID (e.g. I/355/gaiadr3).")
    p_set.set_defaults(func=handle_set)

    # del()
    del_parser = subparsers.add_parser("del", help="Deletes an existing alias in the ATK alias file.")
    del_parser.add_argument("alias", type=str, metavar="<ALIAS>", help="Name of alias to delete.")

    # parse and dispatch
    args = parser.parse_args()
    args.func(args)
