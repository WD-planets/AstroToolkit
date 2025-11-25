import yaml

# class CustomDumper(yaml.SafeDumper):
#     def write_line_break(self, data=None):
#         super().write_line_break(data)
#         if not self.indent:
#             super().write_line_break()


class CustomDumper(yaml.SafeDumper):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Track previous indent and line count
        self._prev_indent = 0
        self._first_survey = True

    def write_line_break(self, data=None):
        # Normal line break
        super().write_line_break(data)

        # Current indent
        indent = getattr(self, "indent", None)

        # --- Insert blank lines only when appropriate ---
        # 1) Before a new top-level section (indent == 0)
        if indent == 0 and self._prev_indent != 0:
            super().write_line_break(data)

        # 2) Between surveys (first-level mapping inside section: indent == 2)
        if indent == 4:
            if not self._first_survey:
                super().write_line_break(data)
            else:
                # Skip extra newline for first survey
                self._first_survey = False

        # Update previous indent
        self._prev_indent = indent

    # Reset first survey flag when starting a new mapping
    def increase_indent(self, flow=False, indentless=False):
        indent = super().increase_indent(flow, indentless)
        # When indent goes from 0 to 2, we are entering first survey
        if getattr(self, "indent", 0) == 0:
            self._first_survey = True
        return indent


def default_printer(data: dict):
    for section, values in data.items():
        print(f"[{section}]")
        if isinstance(values, dict):
            for key, value in values.items():
                print(f"{key} = {value}")
        else:
            print(values)
        print()
