from ..Configuration.baseconfig import ConfigStruct

config = ConfigStruct()
config.read_config()
newline = "\n"


def savedata(struct, fname):
    if fname:
        if not fname.endswith(".fits"):
            fname += ".fits"
    else:
        fname = struct.dataname

    from ..FileHandling.file_writing import generate_local_file

    success = generate_local_file(struct, fname)

    if success:
        config.read_config()
        if config.enable_notifications:
            print(f"Saving data to local storage: {fname}{newline}")

    return fname
