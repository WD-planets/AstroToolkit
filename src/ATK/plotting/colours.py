import matplotlib.colors as mcolors
import numpy as np
from bokeh.palettes import Category10, Category20
from matplotlib.colors import hsv_to_rgb

GRADIENT_MAPS = {
    "green": ("forestgreen", "greenyellow"),
    "red": ("red", "yellow"),
    "blue": ("royalblue", "aqua"),
    "black": ("black", "lightgray"),
    "orange": ("orange", "gold"),
    "purple": ("darkviolet", "orchid"),
}


def generate_palette(n: int):
    """
    Generates a pallete with a given number of colours
    """

    hues = np.linspace(0, 1, n, endpoint=False)
    colours = hsv_to_rgb(np.column_stack([hues, np.full(n, 0.7), np.full(n, 0.9)]))

    return ["#{:02x}{:02x}{:02x}".format(*(int(c * 255) for c in colour)) for colour in colours]


def get_palette(n: int, shift: int = 0):
    """
    Generates a palette of colours with a given shift (used to de-prioritise blue markers as these are hard to see with viridis colour map).
    1-20 colours are generated from Bokeh presets, anything above this is done manually
    """

    if n <= 10:
        colours = Category10[10]
    elif n <= 20:
        colours = Category20[20]
    else:
        return generate_palette(n)

    # apply shift
    if shift:
        colours = colours[shift:] + colours[:shift]

    return colours[:n]


def get_gradient(colour: str, n: int, cycles: int = 1, reverse: bool = False) -> list[str] | str:
    low, high = GRADIENT_MAPS[colour]
    if cycles == 0:
        return mcolors.to_hex(low)

    if reverse:
        low, high = high, low

    points = [low, high]
    for i in range(cycles - 1):
        points.append(low if i % 2 == 0 else high)

    cmap = mcolors.LinearSegmentedColormap.from_list("", points)
    palette = [mcolors.rgb2hex(c) for c in cmap(np.linspace(0, 1, n))]

    return palette


def assign_gradient_palettes(
    n: int, colours: list[str] | None = None, gradient_size: int = 256, cycles: int = 1, reverse: bool = False
) -> list[list[str]]:
    cycle = [c for c in GRADIENT_MAPS if c != "black"]

    if not colours:
        colour_names = [cycle[i % len(cycle)] for i in range(n)]
    else:
        colours = list(colours)
        if len(colours) >= n:
            colour_names = colours[:n]
        else:
            used = list(dict.fromkeys(colours))
            remaining = [c for c in cycle if c not in used]
            result = list(colours)
            for c in remaining:
                if len(result) >= n:
                    break
                result.append(c)
            i = 0
            while len(result) < n:
                result.append(cycle[i % len(cycle)])
                i += 1
            colour_names = result

    palettes = [get_gradient(c, gradient_size, cycles=cycles, reverse=reverse) for c in colour_names]

    return palettes
