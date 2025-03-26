def openFileDialogue():
    import os

    from PyQt5.QtWidgets import QApplication, QFileDialog

    app = QApplication([])
    fname, _ = QFileDialog.getOpenFileName(None, "Open File", str(os.getcwd()), "ATK Data File (*.fits)")

    return fname
