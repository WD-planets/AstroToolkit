import argparse

from AstroToolkit.Config import editConfig, openConfig, resetConfig, showConfig


def main():
    parser = argparse.ArgumentParser()
    sub_parsers = parser.add_subparsers(dest="job")

    sub_parsers.add_parser("reset")
    sub_parsers.add_parser("show")
    sub_parsers.add_parser("open")

    edit_parser = sub_parsers.add_parser("edit")
    edit_parser.add_argument("key", type=str)
    edit_parser.add_argument("value", type=str)

    args = parser.parse_args()

    if args.job == "reset":
        resetConfig()
    elif args.job == "show":
        showConfig()
    elif args.job == "open":
        openConfig()
    elif args.job == "edit":
        editConfig(key=args.key, value=args.value)
