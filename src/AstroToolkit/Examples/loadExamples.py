import glob
import os
from pathlib import Path

file_path = Path(__file__).parent.absolute()

scripts = glob.glob(os.path.join(file_path, "*.py"))
scripts = [Path(x).stem for x in scripts]

for script in scripts:
    __import__(f"AstroToolkit.Examples.{script}")
