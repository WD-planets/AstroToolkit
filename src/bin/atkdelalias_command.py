import argparse

from AstroToolkit.Aliases import delAlias


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("key", type=str)

    args = parser.parse_args()

    delAlias(name=args.key)
