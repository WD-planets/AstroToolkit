"""
This module assists in setting up the PyAOV time series analysis routines.
"""

newline = "\n"


def tsguide(platform="auto"):
    """
    Prints instructions for installing the PyAOV routines on the current platform to stdout. This function can also be called directly from the :ref:`Command-Line`.

    Your platform should be automatically detected. If it isn't, this can be overriden with:

    :param str platform: The platform for which instructions should be output, from:

        - windows
        - linux

    :return: None

    |

    """

    import os
    import platform as platform_get
    from pathlib import Path

    from .Timeseries.pyaov import pyaov

    if platform == "auto":
        plat = platform_get.system().lower()
    else:
        plat = platform

    os.chdir(os.path.join(Path(pyaov.__file__).parent.absolute(), "README"))

    if plat == "windows":
        print(newline * 5)
        plat_str = "Windows"
        with open("Windows_README.txt", "r") as f:
            print(f.read())

    elif plat in ["posix", "linux"]:
        print(newline * 5)
        plat_str = "Linux"
        with open("Linux_README.txt", "r") as f:
            print(f.read())

    else:
        if platform == "auto":
            raise Exception(
                f"Could not detect platform. Instructions may be found for your platform via {os.getcwd()}"
            )
        else:
            raise Exception(
                f"Unsupported platform. Instructions may be found for Windows/Linux via {os.getcwd()}"
            )

    if platform == "auto":
        print(
            f"{newline}{newline}Note: These instructions are for {plat_str} as per your detected platform."
        )
    else:
        pass


def tsbuild(platform="auto"):
    """
    Builds PyAOV time series analysis routines according to current platform. This function can also be called directly from the :ref:`Command-Line`.

    Your platform should be automatically detected. If it isn't, this can be overriden with:

    :param str platform: The current platform, from:

        - windows
        - linux

    :return: None

    |

    """
    import os
    import platform as platform_get
    import subprocess
    import time
    from pathlib import Path

    retry_count = 1

    def do_build(retry_count):
        if platform == "auto":
            plat = platform_get.system().lower()
        else:
            plat = platform

        if plat == "windows":
            from .Timeseries.pyaov import pyaov

            path = os.path.join(Path(pyaov.__file__).parent.absolute(), "buildwin.bat")
            os.chdir(Path(pyaov.__file__).parent.absolute())
            p = subprocess.Popen(
                ["cmd", "/c", str(path)],
                stdout=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = p.communicate()
            lines = stdout.splitlines()

            success = True if lines[-3] == "Compilation Successful." else False
            if not success and retry_count == 0:
                retry_count += 1
                print(f"{newline}Compilation unsuccessful, retrying...{newline}")
                time.sleep(1)
                do_build(retry_count)
            elif success:
                for line in lines:
                    print(line)
            else:
                print(f"{newline}Compilation unsuccessful. Please retry installation.")

        elif plat in ["posix", "linux"]:
            from .Timeseries.pyaov import pyaov

            path = os.path.join(Path(pyaov.__file__).parent.absolute(), "build.sh")
            os.chdir(Path(pyaov.__file__).parent.absolute())
            subprocess.run(["chmod", "+x", str(path)])
            p = subprocess.Popen(
                ["sh", str(path)],
                stdout=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = p.communicate()
            lines = stdout.splitlines()

            success = True if lines[-3] == "Compilation Successful." else False
            if not success and retry_count == 0:
                retry_count += 1
                print(f"{newline}Compilation unsuccessful, retrying...{newline}")
                time.sleep(1)
                do_build(retry_count)
            elif success:
                for line in lines:
                    print(line)
            else:
                print(f"{newline}Compilation unsuccessful. Please retry installation.")

        else:
            if platform == "auto":
                from .Timeseries.pyaov import pyaov

                os.chdir(Path(pyaov.__file__).parent.absolute())
                raise Exception(
                    f"Could not detect platform. Manual compilation may be possible via {os.getcwd()}."
                )
            else:
                raise Exception(
                    f"Unsupported Platform. Manual compilation may be possible via {os.getcwd()}."
                )

    do_build(retry_count)
