import numpy as np
from bokeh.palettes import Category10, Category20
from matplotlib.colors import hsv_to_rgb


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
