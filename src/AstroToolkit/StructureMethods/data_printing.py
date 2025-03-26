def print_data(struct, pprint, print_methods):
    import numpy as np

    print()

    methods = ["savedata", "showdata", "saveplot", "showplot"]
    headings = ["band", "survey"]
    newline = "\n"
    list_exceptions = ["pos", "sources", "identifiers", "positions"]
    header_exceptions = {
        "image_data": "<Image Data>",
        "image_header": "<Image Header>",
        "wcs": "<WCS Object>",
        "overlay": "<Overlay Data>",
    }
    ignore_methods = ["__init__", "__str__"]

    np.set_printoptions(threshold=10, linewidth=1000)

    var_len_list = []
    var_len_group = []

    vars_list = [x for x in vars(struct) if x not in methods]
    if "data" in vars_list:
        vars_list.remove("data")
        vars_list.append("data")

    for var in vars_list:
        if isinstance(getattr(struct, var), (list, dict)) and var not in list_exceptions:
            var_len_list.append(max(var_len_group))
            var_len_group = []
        else:
            var_len_group.append(len(var))
            if var == vars_list[-1]:
                var_len_list.append(max(var_len_group))

    indent_width = 4

    pad_counter = 0
    for var in vars_list:
        pad_length = var_len_list[pad_counter] + 1
        if isinstance(getattr(struct, var), (list, dict)) and var not in list_exceptions:
            pad_counter += 1

        if var not in methods:
            if isinstance(getattr(struct, var), dict):
                print(f"{newline}.{var}:")

                def print_dict(dictionary, indent_index=1):
                    header_len_list = []
                    for key in dictionary:
                        header_len_list.append(len(key))
                    header_pad_length = max(header_len_list) + 1

                    for key, item in dictionary.items():
                        if key not in header_exceptions:
                            if not pprint:
                                print(
                                    f"{' ' * indent_index * indent_width}{(key + ':').ljust(header_pad_length)} {item}"
                                )
                            else:
                                print(
                                    f"{' ' * indent_index * indent_width}{(key + ':').ljust(header_pad_length)} {np.asarray(item)}"
                                )
                        else:
                            print(
                                f"{' ' * indent_index * indent_width}{(key + ':').ljust(header_pad_length)} {header_exceptions[key]}"
                            )
                    if indent_index == 1:
                        print()

                completed = False
                for key, item in getattr(struct, var).items():
                    if isinstance(item, dict):
                        print(f"{' ' * indent_width}{key}:")
                        print_dict(item, indent_index=2)
                        print()
                        completed = True

                if not completed:
                    print_dict(getattr(struct, var))

            elif isinstance(getattr(struct, var), list):
                if var not in list_exceptions:
                    print(f"{newline}.{var}:")
                    for element in getattr(struct, var):
                        if isinstance(element, dict):
                            header_len_list = []
                            for key in element:
                                header_len_list.append(len(key))
                            header_pad_length = max(header_len_list) + 1

                            for key, item in element.items():
                                if key in headings:
                                    print(f"{' ' * indent_width}{key}: {item}")
                                else:
                                    if not pprint:
                                        print(f"{' ' * 2 * indent_width}{(key + ':').ljust(header_pad_length)} {item}")
                                    else:
                                        print(
                                            f"{' ' * 2 * indent_width}{(key + ':').ljust(header_pad_length)} {np.asarray(item)}"
                                        )
                        print()
                else:
                    print(f".{(var + ':').ljust(pad_length)} {getattr(struct, var)}")
            else:
                if var in ["trace", "traces"]:
                    if var == "trace":
                        shift = 1
                    else:
                        shift = 0
                    if not getattr(struct, var):
                        print(f".{(var + ':').ljust(pad_length)} {getattr(struct, var)}")
                        continue

                    split = getattr(struct, var).split("|")
                    if len(split) > 1:
                        print(f".{(var + ':').ljust(pad_length)} {split[0]}")
                        for element in split[1:-1]:
                            print(f"{' ' * (len(var) + shift + pad_length)} {element}")
                        print(f"{' ' * pad_length} {split[-1]}")
                    else:
                        print(f".{(var + ':').ljust(pad_length)} {getattr(struct, var)}")

                else:
                    print(f".{(var + ':').ljust(pad_length)} {getattr(struct, var)}")

    if print_methods:
        import inspect

        methods = [
            method[0]
            for method in inspect.getmembers(struct, predicate=inspect.ismethod)
            if method[0] not in ignore_methods
        ]
        method_str = "Available Methods: "
        for method in methods:
            method_str += f".{method}(), "
        method_str = method_str[:-2]
        print(method_str)
