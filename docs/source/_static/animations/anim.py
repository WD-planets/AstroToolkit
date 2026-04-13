from manim import *

config.background_color = "#0c0c0c"
code = """from ATK import query
from astropy.coordinates import SkyCoord
vMa2 = SkyCoord(12.291, 5.389, unit="deg", frame="icrs")
data = query("vizier", target=vMa2, survey="galex")
data.show()
"""


class Anim(Scene):
    def construct(self):
        initial_code = Code(code_string="", language="python", background="window", background_config={"stroke_color": "maroon"})

        self.play(FadeIn(initial_code))

        cursor = Rectangle(width=0.06, height=0.35, fill_opacity=1, fill_color=WHITE, stroke_width=0)
        self.add(cursor)

        displayed = ""
        for char in code:
            displayed += char
            displayed_code = Code(code_string=displayed, language="python", background="window", tab_width=4)
            displayed_code.scale(0.7)
            displayed_code.move_to(initial_code)
            # typing speed logic
            if char == "\n":
                t = 0.25
            elif char in "().":
                t = 0.08
            else:
                t = 0.025

            self.play(Transform(initial_code, displayed_code), run_time=t, rate_func=linear)

            # move cursor to end of code block
            cursor.next_to(displayed_code, RIGHT, buff=0.05)

        self.play(cursor.animate.set_opacity(0), rate_func=there_and_back, run_time=0.5)

