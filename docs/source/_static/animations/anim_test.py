from manim import *

config.background_color = "#0c0c0c"

CODE_STYLE = dict(language="python", background="window", background_config={"stroke_color": "maroon"}, formatter_style="vim")


def char_speed(ch: str) -> float:
    if ch == "\n":
        return 0.20
    if ch in "()[]{}.,=:\"'":
        return 0.07
    return 0.03


def make_code(text: str, scale: float, center) -> Code:
    mob = Code(code_string=text, **CODE_STYLE)
    mob.scale(scale)
    mob.move_to(center)
    return mob


def cursor_anchor(code_mob: Code) -> np.ndarray:
    para = code_mob.submobjects[2]
    for line in reversed(para.submobjects):
        if len(line):
            g = line[-1]
            return np.array([g.get_right()[0], line.get_center()[1], 0.0])
    return para.get_left() + LEFT * 0.05


def type_string(scene, mob, cursor, text: str, scale: float, center):
    displayed = getattr(mob, "_displayed", "")
    for ch in text:
        displayed += ch
        target = make_code(displayed, scale, center)
        scene.play(
            Transform(mob, target),
            cursor.animate.move_to(cursor_anchor(target) + RIGHT * cursor.width / 2),
            run_time=char_speed(ch),
            rate_func=linear,
        )
    mob._displayed = displayed
    return displayed


def delete_chars(scene, mob, cursor, n: int, scale: float, center):
    displayed = getattr(mob, "_displayed", "")
    for _ in range(n):
        if not displayed:
            break
        ch = displayed[-1]
        displayed = displayed[:-1]
        target = make_code(displayed or " ", scale, center)
        scene.play(
            Transform(mob, target),
            cursor.animate.move_to(cursor_anchor(target) + RIGHT * cursor.width / 2),
            run_time=char_speed(ch),
            rate_func=linear,
        )
    mob._displayed = displayed
    return displayed


def blink(scene, cursor, run_time: float = 0.5):
    scene.play(cursor.animate.set_opacity(0), rate_func=there_and_back, run_time=run_time)


CODE = """\
from ATK import query
from astropy.coordinates import SkyCoord
vMa2 = SkyCoord(12.291, 5.389, unit="deg", frame="icrs")
data = query("vizier", target=vMa2, survey="galex")
data.show()
"""


class Anim(Scene):
    def construct(self):
        _final = make_code(CODE, 1.0, ORIGIN)
        scale = config.frame_width * 0.95 / _final.width
        center = ORIGIN

        # Initial empty window
        mob = make_code("", scale, center)
        mob._displayed = ""
        self.play(FadeIn(mob))

        _chars = [ch for line in _final.submobjects[2] for ch in line]
        gh = (_chars[0].height * scale) if _chars else 0.28
        cursor = Rectangle(width=gh * 0.18, height=gh * 1.1, fill_opacity=1, fill_color=WHITE, stroke_width=0)
        cursor.move_to(mob.submobjects[2].get_left() + LEFT * cursor.width)
        self.add(cursor)

        type_string(self, mob, cursor, CODE, scale, center)
        blink(self, cursor)
        self.wait(1)

        delete_chars(self, mob, cursor, 12, scale, center)
        self.wait(0.3)
        type_string(self, mob, cursor, "data.show()\n", scale, center)
        blink(self, cursor)
        self.wait(0.5)

        delete_chars(self, mob, cursor, 15, scale, center)
        self.wait(0.2)
        type_string(self, mob, cursor, '"gaia")\ndata.show()\n', scale, center)
        blink(self, cursor)
        self.wait(1)
