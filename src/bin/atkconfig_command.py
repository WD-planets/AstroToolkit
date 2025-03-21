import argparse

from AstroToolkit.Config import editConfig, openConfig, resetConfig, showConfig


def main():
    parser = argparse.ArgumentParser(description="Facilitates the viewing and editing of the ATK config file.")
    sub_parsers = parser.add_subparsers(dest="job")

    sub_parsers.add_parser("reset", help="Resets the config to default values")
    sub_parsers.add_parser("show", help="Outputs the config to stdout")
    sub_parsers.add_parser("open", help="Opens the config in the default text editor")

    edit_parser = sub_parsers.add_parser("edit", help="Sets the value of a given config key")
    edit_parser.add_argument("key", type=str, help="Name of config key to edit")
    edit_parser.add_argument("value", type=str, help="New value of config key")
    args = parser.parse_args()

    if args.job == "reset":
        resetConfig()
    elif args.job == "show":
        showConfig()
    elif args.job == "open":
        openConfig()
    elif args.job == "edit":
        editConfig(key=args.key, value=args.value)
