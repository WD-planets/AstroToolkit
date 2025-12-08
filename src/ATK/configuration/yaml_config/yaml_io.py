import yaml


class CustomDumper(yaml.SafeDumper):
    """
    Overrides the default yaml dumper to add additional newlines for readability
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._prev_indent = 0

        from .YAMLConfig import YAML_INDENT

        self.INDENT = YAML_INDENT

    def write_line_break(self, data=None):
        # normal line break
        super().write_line_break(data)

        # get current indent
        current_indent = getattr(self, "indent", None)

        if self._prev_indent is None or current_indent is None:
            pass

        # print newlines between surveys and sections
        elif self._prev_indent > current_indent:
            # print extra newline before new section
            if not current_indent:
                super().write_line_break()
            super().write_line_break(data)

        # update previous indent
        self._prev_indent = current_indent


def default_printer(data: dict):
    """
    Default printer for printing yaml files to stdout
    """

    for section, values in data.items():
        print(f"[{section}]")
        if isinstance(values, dict):
            for key, value in values.items():
                print(f"{key} = {value}")
        else:
            print(values)
        print()


def seq_representer(dumper, seq):
    """
    Overrides CustomDumper to print only lists of values in flow style (i.e. [one,two,three])
    """

    return dumper.represent_sequence("tag:yaml.org,2002:seq", seq, flow_style=True)


CustomDumper.add_representer(list, seq_representer)
