import argparse

from AstroToolkit.Configuration.baseconfig import ConfigStruct
from AstroToolkit.Tools import search


def main():
    parser = argparse.ArgumentParser(
        description="Performs a search for an object in SIMBAD or Vizier, given its position or Gaia Source ID. The result is opened in the default browser."
    )

    parser.add_argument("kind", choices=["simbad", "vizier"], help="Where to perform search", type=str.lower)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--source", nargs=1, type=int, help="Gaia DR3 Source ID")
    group.add_argument("--pos", nargs=2, type=float, help="Position in degrees", metavar=("RA", "DEC"))
    parser.add_argument("-r", help="Radius of search in arcseconds", type=float, metavar=("RADIUS"))

    args = parser.parse_args()

    if isinstance(args.source, list):
        args.source = args.source[0]

    if not args.r:
        config = ConfigStruct()
        config.read_config()
        args.radius = config.search_radius

    search(kind=args.kind, pos=args.pos, source=args.source, radius=args.r)


def main_new1():
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=["simbad", "vizier"], help="Where to perform search")
    sub_parsers = parser.add_subparsers(dest="target_kind", help="Format of target")

    pos_parser = sub_parsers.add_parser("pos")
    pos_parser.add_argument("ra", type=float)
    pos_parser.add_argument("dec", type=float)
    pos_parser.add_argument("radius", nargs="?", type=float)

    source_parser = sub_parsers.add_parser("source")
    source_parser.add_argument("source", type=int)
    source_parser.add_argument("radius", nargs="?", type=float)

    args = parser.parse_args()

    if not args.radius:
        config = ConfigStruct()
        config.read_config()
        args.radius = config.search_radius

    search(kind=args.kind, pos=args.pos, source=args.source, radius=args.radius)


"""
def main_old():
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
"""
