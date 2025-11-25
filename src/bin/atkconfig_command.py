import argparse

from ATK.configuration.base_config import BASE_CONFIG

sections = list(BASE_CONFIG.as_dict().keys())


def handle_set(args):
    BASE_CONFIG._set(args.section.lower(), args.key.lower(), args.value)


def handle_reset(args):
    BASE_CONFIG._reset()


def handle_open(args):
    BASE_CONFIG._open()


def handle_show(args):
    BASE_CONFIG._show()


def main():
    parser = argparse.ArgumentParser(description="Facilitates viewing and editing of the ATK config file.")
    subparsers = parser.add_subparsers(dest="job", required=True)

    # show()
    p_show = subparsers.add_parser("show", help="Prints the ATK config file to stdout.")
    p_show.set_defaults(func=handle_show)

    # reset()
    p_reset = subparsers.add_parser("reset", help="Reset the ATK config file to its default state.")
    p_reset.set_defaults(func=handle_reset)

    # open()
    p_open = subparsers.add_parser("open", help="Opens the ATK config file in the default editor.")
    p_open.set_defaults(func=handle_open)

    # set()
    p_set = subparsers.add_parser("set", help="Sets the value of a key in a given section of the ATK config file.")
    p_set.add_argument(
        "section",
        type=str,
        metavar="<SECTION>",
        help=f"Name of config section in which key is located, from: {', '.join(sections)}.",
        choices=sections,
    )
    p_set.add_argument("key", type=str, metavar="<KEY>", help="Name of config key in given section.")
    p_set.add_argument("value", type=str, metavar="<VALUE>", help="New value of config key.")
    p_set.set_defaults(func=handle_set)

    args = parser.parse_args()
    args.func(args)
