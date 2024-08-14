from bokeh.layouts import column
from bokeh.models import Button, CustomJS, InlineStyleSheet

def get_search_buttons(radius, grid_size, source=None, pos=None):
    from ..Misc.search import DatapageButtons

    query_object = DatapageButtons(radius=radius, source=source, pos=pos)
    query_object.get_params()
    simbad_url, vizier_url = query_object.get_urls()

    button_width, button_height = (
        int(int(grid_size) / 2),
        int(int(grid_size) / 4),
    )

    button_style = InlineStyleSheet(
        css=".bk-btn {--base-font:var(--bokeh-base-font, Helvetica, Arial, sans-serif);--mono-font:var(--bokeh-mono-font, monospace);--font-size:var(--bokeh-font-size, 25px);--line-height:calc(20 / 14);--line-height-computed:calc(var(--font-size) * var(--line-height));--border-radius:4px;--padding-vertical:6px;--padding-horizontal:12px;--bokeh-top-level:1000;}:host{box-sizing:border-box;font-family:var(--base-font);font-size:var(--font-size);line-height:var(--line-height);}*,*:before,*:after{box-sizing:inherit;font-family:inherit;}pre,code{font-family:var(--mono-font);margin:0;}"
    )

    simbad_button = Button(
        label="SIMBAD Search",
        button_type="light",
        height=button_height,
        width=button_width,
        sizing_mode="stretch_width",
        align="end",
        stylesheets=[button_style]
    )
    vizier_button = Button(
        label="Vizier Search",
        button_type="light",
        height=button_height,
        width=button_width,
        sizing_mode="stretch_width",
        align="end",
        stylesheets=[button_style]
    )

    simbad_button_js = CustomJS(
        args=dict(url=simbad_url),
        code="""
        window.open(url)
    """,
    )
    simbad_button.js_on_event("button_click", simbad_button_js)

    vizier_button_js = CustomJS(
        args=dict(url=vizier_url),
        code="""
        window.open(url)
    """,
    )
    vizier_button.js_on_event("button_click", vizier_button_js)

    simbad_button.margin = [
        round(button_height),
        round(button_height),
        round(button_height / 4),
        round(button_height),
    ]
    vizier_button.margin = [
        0,
        round(button_height),
        round(button_height / 4),
        round(button_height),
    ]

    buttons = column(simbad_button, vizier_button, align="center")

    return buttons
