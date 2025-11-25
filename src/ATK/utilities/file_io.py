import os
import shutil
import subprocess
import sys
from pathlib import Path

# possible ways of opening files across linux distributions, in order of priority
LINUX_OPENERS = {
    "xdg-open": ["xdg-open"],
    "vim": ["vim"],
    "vi": ["vi"],
    "nano": ["nano"],
    "gio": ["gio", "open"],
    "gvfs-open": ["gvfs-open"],
}


def open_file(path: Path):
    # Windows
    if sys.platform.startswith("win"):
        os.startfile(path)

    # MacOS
    elif sys.platform == "darwin":
        subprocess.run(["open", path], check=False)

    # Linux
    else:
        opener = None
        for tool, command in LINUX_OPENERS.items():
            if shutil.which(tool):
                opener = command
                break

        if opener:
            print([opener, path])
            subprocess.run([*opener, path], check=False)
        else:
            print(f"No system opener found; file is located at {path}")
