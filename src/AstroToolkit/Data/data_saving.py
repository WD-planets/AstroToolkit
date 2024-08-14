from ..Configuration.baseconfig import ConfigStruct

config = ConfigStruct()
config.read_config()
newline = "\n"


def savedata(self, fname=None):
    if fname:
        if self.kind != "image":
            if not fname.endswith(".csv"):
                fname += ".csv"
        else:
            if not fname.endswith(".fits"):
                fname += ".fits"
    else:
        fname = self.dataname

    from ..FileHandling.file_writing import generate_local_file

    success = generate_local_file(self, fname)

    if success:
        config.read_config()
        if config.enable_notifications == "True":
            print(f"Saving data to local storage: {fname}{newline}")

    return fname
