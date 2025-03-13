import argparse

from AstroToolkit.Configuration.baseconfig import ConfigStruct
from AstroToolkit.Tools import search


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", type=str)
    parser.add_argument("target", nargs="+")

    args = parser.parse_args()

    if len(args.target) > 2:
        pos = [float(args.target[0]), float(args.target[1])]
        source = None
        radius = float(args.target[2])

    elif len(args.target) > 1:
        if float(args.target[0]) > pow(10, 10):
            source = args.target[0]
            radius = args.target[1]
            pos = None
        else:
            pos = [float(args.target[0]), float(args.target[1])]
            source = None
            radius = None
    else:
        source = int(args.target[0])
        pos = None
        radius = None

    if not radius:
        config = ConfigStruct()
        config.read_config()
        radius = config.search_radius

    search(args.kind, pos=pos, source=source, radius=radius)
