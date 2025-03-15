import argparse

from AstroToolkit.Aliases import addAlias


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("key", type=str)
    parser.add_argument("value", type=str)

    args = parser.parse_args()

    addAlias(name=args.key, id=args.value)
