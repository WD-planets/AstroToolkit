import argparse

from ATK.configuration.base_config import BASE_CONFIG


def main():
    parser = argparse.ArgumentParser(description="Facilitates the viewing and editing of the ATK config file.")
    sub_parsers = parser.add_subparsers(dest="job")

    sub_parsers.add_parser("reset", help="Resets the config to default values")
    sub_parsers.add_parser("show", help="Outputs the config to stdout")
    sub_parsers.add_parser("open", help="Opens the config in the default text editor")

    edit_parser = sub_parsers.add_parser("set", help="Sets the value of a given config key in a given section")
    edit_parser.add_argument("section", type=str, help="Name of section in which a key should be edited")
    edit_parser.add_argument("key", type=str, help="Name of config key to edit")
    edit_parser.add_argument("value", type=str, help="New value of config key")
    args = parser.parse_args()

    if args.job == "reset":
        BASE_CONFIG._reset()
    elif args.job == "show":
        BASE_CONFIG._show()
    elif args.job == "open":
        BASE_CONFIG._open()
    elif args.job == "set":
        BASE_CONFIG._set(args.section, args.key, args.value)
