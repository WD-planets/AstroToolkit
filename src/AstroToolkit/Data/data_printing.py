def print_data(struct, raw=False):
    import numpy as np

    print()

    methods = ["savedata", "showdata", "saveplot", "showplot"]
    headings = ["band", "survey"]
    newline = "\n"
    list_exceptions = ["pos", "sources", "identifiers"]
    header_exceptions = {
        "image_data": "<Image Data>",
        "image_header": "<Image Header>",
        "wcs": "<WCS Object>",
        "overlay": "<Overlay Data>",
    }

    np.set_printoptions(threshold=10, linewidth=1000)

    var_len_list = []
    var_len_group = []

    vars_list = [x for x in vars(struct) if x not in methods]
    for var in vars_list:
        if (
            isinstance(getattr(struct, var), (list, dict))
            and var not in list_exceptions
        ):
            var_len_list.append(max(var_len_group))
            var_len_group = []
        else:
            var_len_group.append(len(var))
            if var == vars_list[-1]:
                var_len_list.append(max(var_len_group))

    indent_width = 4

    pad_counter = 0
    for var in vars(struct):
        pad_length = var_len_list[pad_counter] + 1
        if (
            isinstance(getattr(struct, var), (list, dict))
            and var not in list_exceptions
        ):
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
                            if raw:
                                print(
                                    f"{' ' * indent_index * indent_width}{(key+':').ljust(header_pad_length)} {item}"
                                )
                            else:
                                print(
                                    f"{' ' * indent_index * indent_width}{(key+':').ljust(header_pad_length)} {np.asarray(item)}"
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
                                    print(f"{" " * indent_width}{key}: {item}")
                                else:
                                    if raw:
                                        print(
                                            f"{" " * 2*indent_width}{(key+':').ljust(header_pad_length)} {item}"
                                        )
                                    else:
                                        print(
                                            f"{" " * 2*indent_width}{(key+':').ljust(header_pad_length)} {np.asarray(item)}"
                                        )
                        print()
                else:
                    print(f".{(var+':').ljust(pad_length)} {getattr(struct, var)}")
            else:
                print(f".{(var+':').ljust(pad_length)} {getattr(struct, var)}")
