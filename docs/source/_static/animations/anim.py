import re

from manim import *

config.background_color = "#ffffff"
config.max_files_cached = 500
config.disable_caching = True
config.pixel_width = 1080
config.pixel_height = 1080

MIN_WIDTH = 4
SCALE = 0.6
DELAY_MULT = 0.75
GRID_COLS = 8

TERMINAL = False


def get_char_speed(char):
    if char == "\n":
        t = 0.25 * DELAY_MULT
    elif char in "().":
        t = 0.15 * DELAY_MULT
    else:
        t = 0.1 * DELAY_MULT

    return t


class Anim(MovingCameraScene):
    def render_frame(self, cursor_pos=None):
        code_str = self.displayed if self.displayed.strip() else ""
        displayed_code = Code(code_string=code_str, language="python", background="rectangle", tab_width=4, paragraph_config={"font": "Ubuntu Sans Mono"})
        displayed_code.scale(SCALE)

        bg = displayed_code.submobjects[0]
        if bg.width < MIN_WIDTH:
            new_bg = RoundedRectangle(
                width=MIN_WIDTH,
                height=bg.height,
                fill_color=bg.get_fill_color(),
                fill_opacity=bg.get_fill_opacity(),
                stroke_color=GREY,
                stroke_width=bg.get_stroke_width(),
                corner_radius=0.1,
            )
            new_bg.move_to(bg.get_center())
            content_offset = new_bg.get_left()[0] - bg.get_left()[0]
            displayed_code.submobjects[1].shift(RIGHT * content_offset)
            displayed_code.submobjects[2].shift(RIGHT * content_offset)
            displayed_code.submobjects[0].become(new_bg)

        # centre align
        displayed_code.set_x(0)
        frame_top = self.camera.frame.get_top()[1]
        displayed_code.set_y(frame_top - displayed_code.height / 2 - 0.3)
        if self.terminal:
            self.terminal.next_to(self.initial_code, DOWN, buff=0.8)
            self.terminal.set_x(0)

        self.remove(self.initial_code)
        self.remove(self.cursor)
        self.add(displayed_code)
        self.add(self.cursor)

        # Render a trimmed copy just to find cursor position
        cursor_str = self.displayed[:cursor_pos] if cursor_pos is not None else self.displayed
        cursor_code = Code(
            code_string=cursor_str if cursor_str.strip() else "",
            language="python",
            background="rectangle",
            tab_width=4,
            paragraph_config={"font": "Ubuntu Sans Mono"},
        )
        cursor_code.scale(SCALE)
        cursor_code.align_to(displayed_code, UP + LEFT)

        code_paragraph = cursor_code.submobjects[2]
        for line_group in reversed(code_paragraph.submobjects):
            chars = [m for m in line_group.submobjects if m.width > 0]
            if chars:
                self.cursor.next_to(chars[-1], RIGHT, buff=0.04)
                self.cursor.set_y(line_group.get_center()[1])  # use line, not char

                break

        self.initial_code = displayed_code

    @property
    def grid_unit(self):
        return self.camera.frame.width / GRID_COLS

    def type_code(self, text):
        for char in text:
            self.displayed += char
            self.render_frame()
            self.wait(get_char_speed(char))

    def delete_chars(self, n):
        for _ in range(n):
            self.displayed = self.displayed[:-1]
            self.render_frame()
            self.wait(0.05)

    def replace_code(self, old, new):
        pos = self.displayed.find(old)
        if pos == -1:
            return

        # backspace through old right to left, cursor follows
        for i in range(len(old)):
            remove_pos = pos + len(old) - i - 1
            self.displayed = self.displayed[:remove_pos] + self.displayed[remove_pos + 1 :]
            self.render_frame(cursor_pos=remove_pos)
            self.wait(0.05)

        # type new at pos, cursor follows
        for i, char in enumerate(new):
            self.displayed = self.displayed[: pos + i] + char + self.displayed[pos + i :]
            self.render_frame(cursor_pos=pos + i + 1)
            self.wait(get_char_speed(char))

    def delete_line(self, line_index, animate=False):
        lines = self.displayed.split("\n")

        if line_index < 0:
            line_index += len(lines)
        if line_index < 0 or line_index >= len(lines):
            return

        # Compute cursor position at start of the line
        cursor_pos = sum(len(line) + 1 for line in lines[:line_index])

        if animate:
            self.render_frame(cursor_pos=cursor_pos)
            self.wait(0.1)

        # Remove the line
        del lines[line_index]

        # Rebuild string
        self.displayed = "\n".join(lines)

        # Render once after deletion
        self.render_frame(cursor_pos=cursor_pos)
        self.cursor_blink(0.25)

    def show_output(self, output: str, pause_before=0.5):
        global TERMINAL

        self.run_code(pause_before)
        terminal = Code(
            code_string=output,
            language="text",
            background="window",
            tab_width=4,
            background_config={"stroke_color": GREY},
            paragraph_config={"font": "Ubuntu Sans Mono"},
            add_line_numbers=False,
        )
        terminal.scale(SCALE)
        terminal.next_to(self.initial_code, DOWN, buff=0.8)

        if TERMINAL:
            self.remove(self.terminal)
            self.add(terminal)
        else:
            self.play(FadeIn(terminal))

        self.terminal = terminal
        TERMINAL = True

        self.cursor_blink(3.0)

    def show_png(self, path, grid_w, pause_before=0.5):
        global TERMINAL

        self.run_code(0.5)

        img = ImageMobject(path)
        img.scale_to_fit_width(grid_w * self.grid_unit)
        img.next_to(self.initial_code, DOWN, buff=0.8)

        # frame = RoundedRectangle(width=img.width + 0.3, height=img.height + 0.3, corner_radius=0.1, stroke_color=GREY, stroke_width=1)
        # frame.move_to(img)
        # group = Group(frame, img)

        if TERMINAL:
            self.remove(self.terminal)
        self.add(img)

        self.terminal = img
        TERMINAL = True

        self.cursor_blink(3.0)

    def show_svg(self, path, grid_w, pause_before=0.5):
        global TERMINAL

        self.run_code(0.5)

        img = SVGMobject(path)
        img.scale_to_fit_width(grid_w * self.grid_unit)
        img.next_to(self.initial_code, DOWN, buff=0.8)

        if TERMINAL:
            self.remove(self.terminal)
        self.add(img)

        self.terminal = img
        TERMINAL = True

        self.cursor_blink(3.0)

    def run_code(self, runtime=1.5):
        spinner = Arc(radius=0.14, angle=TAU * 0.75, color=GREY, stroke_width=2.5)
        spinner.next_to(self.initial_code, DOWN, buff=0.3)
        spinner.set_x(0)

        self.add(spinner)
        self.start_cursor_blink()

        self.play(Rotate(spinner, TAU * 2, run_time=runtime, rate_func=linear))

        self.stop_cursor_blink()
        self.remove(spinner)

    def start_cursor_blink(self, rate=2):
        def updater(m, dt):
            m.time = getattr(m, "time", 0) + dt
            m.set_opacity(0 if int(m.time * rate) % 2 == 0 else 1)

        self.cursor.add_updater(updater)

    def stop_cursor_blink(self):
        self.cursor.clear_updaters()
        self.cursor.set_opacity(1)

    def cursor_blink(self, runtime, rate=2):
        self.start_cursor_blink(rate)
        self.wait(runtime)
        self.stop_cursor_blink()

    def construct(self):
        self.terminal = None
        self.displayed = ""
        self.initial_code = Code(
            code_string="",
            language="python",
            background="rectangle",
            background_config={"stroke_color": "maroon"},
            paragraph_config={"font": "Ubuntu Sans Mono"},
        )
        self.initial_code.scale(SCALE)
        self.initial_code.to_edge(UP, buff=0.2)
        self.cursor = Rectangle(width=0.05, height=0.25, fill_opacity=1, fill_color=WHITE, stroke_width=0)
        self.cursor.scale(SCALE * 1.4)
        self.add(self.cursor)

        self.type_code("from ATK import query\n")
        self.type_code("from astropy.coordinates import SkyCoord\n")
        self.type_code("\n")
        self.type_code("vMa2 = SkyCoord(12.291, 5.389, unit='deg', frame='icrs')\n")
        self.type_code("data = query('vizier', targets=vMa2, survey='galex')\n")
        self.type_code("data.show()")

        output = """
        <vizier DataSet>

        .kind:      vizier
        .targets:   12.291° 5.389° (icrs, 2000-01-01T00:00:00.000, 3.0″)
        .exception: False
        .data:      <empty list>
        """

        self.show_output(output)
        self.replace_code("SkyCoord(12.291, 5.389, unit='deg', frame='icrs')", "2552928187080872832")

        output = """
        <vizier DataSet>

        .kind:      vizier
        .targets:   2552928187080872832 | 12.297° 5.377° (icrs, 2016-01-01T00:00:00.000, 3.0″)
        .exception: False
        .data:                
              <galex (II/335/galex_ais) Record>
                    survey:     galex
                    catalogue:  II/335/galex_ais
                    correction: full
                    search_pos: 12.294° 5.384° (icrs, 2006-08-01T00:00:00.000)
                    table:               
                          (astropy.Table)
                                RAJ2000:  [12.293] °
                                DEJ2000:  [5.384] °
                                Name:     [GALEX J004910.4+052300]
                                objid:    [6380239640935270574]
                                ...
        """

        self.show_output(output)
        self.replace_code("2552928187080872832", "SkyCoord(12.291, 5.389, unit='deg', frame='icrs')")
        self.replace_code("data", "img")
        self.replace_code("vizier", "image")
        self.replace_code("'galex'", "'panstarrs', band='y', size=2 * u.arcmin, overlays=['galex']")
        self.replace_code("data.show()", "img.show()")

        output = """
        <image DataSet>

        .kind:      image
        .targets:   12.291° 5.389° (icrs, 2000-01-01T00:00:00.000)
        .exception: False
        .data:                
              <panstarrs g-band Image>
                    survey:     panstarrs
                    correction: none
                    search_pos: 12.291° 5.389° (icrs, 2011-12-16T21:11:06.224)
                    band:       g
                    size:       30.0″
                    epoch:      2011-12-16T21:11:06.224
                    hdu:        <ImageHDU>
                    wcs:        <WCS>
                """

        self.show_output(output)
        self.replace_code(".show()", ".open()")
        self.show_png("uncorrected_image.png", 4.5)
        self.replace_code("SkyCoord(12.291, 5.389, unit='deg', frame='icrs')", "2552928187080872832")
        self.show_png("corrected_image.png", 4.5)
        self.delete_line(-1)
        self.delete_line(-1)
        self.replace_code("vMa2 = 2552928187080872832", "ARSco = 6050296829033196032")
        self.type_code("\nlc = query('lightcurve', targets=ARSco, survey='asassn')")
        self.type_code("\nlc.open()")
        self.show_png("lightcurve.png", 4.5)
        self.replace_code("lc.open()", "pspec = lc.apply('pspec', fmin=0, fmax=10, samples=100000, inplace=False)")
        self.type_code("\npspec.open()")
        self.show_png("pspec.png", 4.5)
        self.replace_code("pspec", "fold")
        self.replace_code("'pspec'", "'fold'")
        self.replace_code("pspec", "fold")
        self.show_png("folded_lightcurve.png", 4.5)
        self.delete_line(-1)
        self.delete_line(-1)
        self.delete_line(-1)
        self.replace_code("ARSco = 6050296829033196032", "HuLeo = 587316166180416640")
        self.type_code("\nspec = query('spectrum', targets=HuLeo, survey='sdss')")
        self.type_code("\nspec.open()")
        self.show_png("spec.png", 4.5)
        self.replace_code("spec.open()", "fitted = spec.apply('fit', inplace=False)")
        self.type_code("\nfitted.open()")
        self.show_png("fitted_spec.png", 4.5)
        self.delete_line(-1)
        self.delete_line(-1)
        self.delete_line(-1)
        self.type_code("\nsed = query('sed', targets=HuLeo)")
        self.type_code("\nsed.open()")
        self.show_png("sed.png", 4.5)
        self.delete_line(-1)
        self.delete_line(-1)
        self.type_code("\nhrd = query('hrd', targets=HuLeo)")
        self.type_code("\nhrd.open()")
        self.show_png("hrd.png", 4.5)

        self.stop_cursor_blink()
        mobjects_to_fade = [self.initial_code, self.cursor]
        if self.terminal:
            mobjects_to_fade.append(self.terminal)
        self.play(*[FadeOut(m) for m in mobjects_to_fade], run_time=2.0)


"""
add a small buffer in height so that entering e.g. a "g" or a "y" doesn't change the size?
stop cursor from going up a line when at the zero position?
change spectral fitting params to fit H-alpha
fix hrd
look into file size of mp4 file
"""
