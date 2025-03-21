import argparse

from AstroToolkit.Aliases import (addAlias, delAlias, openAliases,
                                  resetAliases, showAliases)


def main():
    parser = argparse.ArgumentParser(description="Facilitates the viewing and editing of the ATK aliases file.")
    sub_parsers = parser.add_subparsers(dest="job")

    sub_parsers.add_parser("reset", help="Resets the alias list to its default state")
    sub_parsers.add_parser("show", help="Outputs the alias list to stdout")
    sub_parsers.add_parser("open", help="Opens the alias list in the default text editor")

    add_parser = sub_parsers.add_parser("add", help="Adds an alias for a Vizier Catalogue ID")
    add_parser.add_argument("alias", type=str, help="Alias name")
    add_parser.add_argument("id", type=str, help="Vizier Catalogue ID")

    del_parser = sub_parsers.add_parser("del", help="Deletes an existing alias")
    del_parser.add_argument("alias", type=str, help="Alias name")

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
