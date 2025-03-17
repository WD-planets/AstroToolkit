import argparse

from AstroToolkit.Aliases import (addAlias, delAlias, openAliases,
                                  resetAliases, showAliases)


def main():
    parser = argparse.ArgumentParser()
    sub_parsers = parser.add_subparsers(dest="job")

    sub_parsers.add_parser("reset")
    sub_parsers.add_parser("show")
    sub_parsers.add_parser("open")

    add_parser = sub_parsers.add_parser("add")
    add_parser.add_argument("alias", type=str)
    add_parser.add_argument("id", type=str)

    del_parser = sub_parsers.add_parser("del")
    del_parser.add_argument("alias", type=str)

    args = parser.parse_args()

    if args.job == "reset":
        resetAliases()
    elif args.job == "show":
        showAliases()
    elif args.job == "open":
        openAliases()
    elif args.job == "add":
        addAlias(name=args.alias, id=args.id)
    elif args.job == "del":
        delAlias(name=args.alias)
