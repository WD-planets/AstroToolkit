import numpy as np
from bokeh.palettes import Category10, Category20
from matplotlib.colors import hsv_to_rgb


def generate_palette(n: int):
    hues = np.linspace(0, 1, n, endpoint=False)
    colours = hsv_to_rgb(np.column_stack([hues, np.full(n, 0.7), np.full(n, 0.9)]))

    return ["#{:02x}{:02x}{:02x}".format(*(int(c * 255) for c in colour)) for colour in colours]


def get_palette(n: int, shift: int = 0):
    if n <= 10:
        colours = Category10[10]
    elif n <= 20:
        colours = Category20[20]
    else:
        return generate_palette(n)

    if shift:
        colours = colours[shift:] + colours[:shift]

    return colours[:n]
