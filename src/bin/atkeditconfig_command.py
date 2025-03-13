import argparse

from AstroToolkit.Config import editconfig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("key", type=str)
    parser.add_argument("value", type=str)

    args = parser.parse_args()

    editconfig(key=args.key, value=args.value)
