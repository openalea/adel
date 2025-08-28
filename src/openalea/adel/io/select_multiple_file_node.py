from openalea.core import path as path, Node
from openalea.core.observer import lock_notify
from PyQt4.QtGui import *
from PyQt4.QtCore import *


class SelectMultipleFiles(NodeWidget, QDialog):
    def __init__(self, node, parent):
        QDialog.__init__(self, parent)
        NodeWidget.__init__(self, node)

        self.gridlayout = QGridLayout(self)
        self.gridlayout.setMargin(3)
        self.gridlayout.setSpacing(5)

        self.package_lineedit_label = QLabel("1. Set the package", self)
        self.gridlayout.addWidget(self.package_lineedit_label, 0, 0)

        self.package_lineedit = QLineEdit(self)
        self.cursorPosition_package_lineedit = self.package_lineedit.cursorPosition()
        self.gridlayout.addWidget(self.package_lineedit, 0, 1, 1, 4)
        self.connect(
            self.package_lineedit, SIGNAL("textChanged(QString)"), self.package_changed
        )

        self.select_files_label = QLabel(
            "2. Set the patterns and select the \
files",
            self,
        )
        self.gridlayout.addWidget(self.select_files_label, 1, 0, 1, 4)

        self.nb_of_inputs = node.get_nb_input()
        self.nb_of_files = (self.nb_of_inputs - 1) / 2

        self.methods_mapping = {}
        self.files_mapping = {}

        for i in range(self.nb_of_files):
            pattern_input_number = i + 1
            filename_input_number = pattern_input_number + self.nb_of_files
            pattern_input_port = node.get_input_port(pattern_input_number)
            name = pattern_input_port["name"]
            self.methods_mapping[pattern_input_number] = {}

            def pattern_changed(self, pattern_input_number):
                pattern = self.files_mapping[pattern_input_number][
                    "pattern_lineedit"
                ].text()
                self.files_mapping[pattern_input_number][
                    "pattern_lineedit_cursorPosition"
                ] = self.files_mapping[pattern_input_number][
                    "pattern_lineedit"
                ].cursorPosition()
                self.node.set_input(pattern_input_number, str(pattern))

            self.methods_mapping[pattern_input_number]["pattern_changed"] = (
                pattern_changed
            )

            def filename_changed(self, pattern_input_number):
                filename_input_number = pattern_input_number + self.nb_of_files
                filename = self.files_mapping[pattern_input_number][
                    "filenames_combobox"
                ].currentText()
                self.node.set_input(filename_input_number, str(filename))

            self.methods_mapping[pattern_input_number]["filename_changed"] = (
                filename_changed
            )

            @lock_notify
            def update_pattern_input(self, pattern_input_number):
                pattern = self.node.get_input(pattern_input_number)
                if self.updating or pattern is None:
                    return
                self.updating = True
                self.files_mapping[pattern_input_number]["pattern_lineedit"].setText(
                    pattern
                )
                self.updating = False

            self.methods_mapping[pattern_input_number]["update_pattern_input"] = (
                update_pattern_input
            )

            @lock_notify
            def update_filenames_combobox(self, pattern_input_number):
                package_name = self.node.get_input(0)
                pattern = self.node.get_input(pattern_input_number)
                if package_name is None or pattern is None:
                    return
                directory_path = path.path(node._pd([package_name])[0])
                if self.updating:
                    return
                self.updating = True
                self.files_mapping[pattern_input_number]["filenames_combobox"].clear()
                self.files_mapping[pattern_input_number][
                    "filename2filepath_mapping"
                ].clear()
                if directory_path.isdir():
                    paths = directory_path.glob(pattern)
                else:
                    paths = []
                filenames = [a_path.name for a_path in paths]
                self.files_mapping[pattern_input_number][
                    "filename2filepath_mapping"
                ].update(dict(list(zip(filenames, paths))))
                self.files_mapping[pattern_input_number]["filenames_combobox"].addItems(
                    list(
                        self.files_mapping[pattern_input_number][
                            "filename2filepath_mapping"
                        ].keys()
                    )
                )
                self.updating = False

            self.methods_mapping[pattern_input_number]["update_filenames_combobox"] = (
                update_filenames_combobox
            )

            @lock_notify
            def update_filepath_output(self, pattern_input_number):
                output_filepath = None
                package_name = self.node.get_input(0)
                filename_input_number = pattern_input_number + self.nb_of_files
                input_filename = self.node.get_input(filename_input_number)
                if package_name is None or input_filename is None:
                    return
                if self.updating:
                    return
                self.updating = True
                input_filename_index = self.files_mapping[pattern_input_number][
                    "filenames_combobox"
                ].findText(input_filename)
                self.files_mapping[pattern_input_number][
                    "filenames_combobox"
                ].setCurrentIndex(input_filename_index)
                output_filepath = self.files_mapping[pattern_input_number][
                    "filename2filepath_mapping"
                ].get(input_filename)
                self.node.set_output(i, output_filepath)
                self.updating = False

            self.methods_mapping[pattern_input_number]["update_filepath_output"] = (
                update_filepath_output
            )

            self.files_mapping[pattern_input_number] = {}
            self.files_mapping[pattern_input_number]["label"] = QLabel(name, self)
            self.gridlayout.addWidget(
                self.files_mapping[pattern_input_number]["label"], i + 2, 0
            )
            self.files_mapping[pattern_input_number]["pattern_lineedit"] = QLineEdit(
                self
            )
            self.files_mapping[pattern_input_number][
                "pattern_lineedit_cursorPosition"
            ] = self.files_mapping[pattern_input_number][
                "pattern_lineedit"
            ].cursorPosition()
            self.gridlayout.addWidget(
                self.files_mapping[pattern_input_number]["pattern_lineedit"],
                i + 2,
                1,
                1,
                2,
            )
            self.pattern_lineedit_signal_mapper = QSignalMapper(self)
            self.connect(
                self.files_mapping[pattern_input_number]["pattern_lineedit"],
                SIGNAL("textChanged(QString)"),
                self.pattern_lineedit_signal_mapper,
                SLOT("map()"),
            )
            self.pattern_lineedit_signal_mapper.setMapping(
                self.files_mapping[pattern_input_number]["pattern_lineedit"],
                pattern_input_number,
            )
            self.connect(
                self.pattern_lineedit_signal_mapper,
                SIGNAL("mapped(int)"),
                self.map_to_pattern_method,
            )
            self.files_mapping[pattern_input_number]["filenames_combobox"] = QComboBox(
                self
            )
            self.gridlayout.addWidget(
                self.files_mapping[pattern_input_number]["filenames_combobox"],
                i + 2,
                3,
                1,
                2,
            )
            self.filenames_combobox_signal_mapper = QSignalMapper(self)
            self.connect(
                self.files_mapping[pattern_input_number]["filenames_combobox"],
                SIGNAL("currentIndexChanged(QString)"),
                self.filenames_combobox_signal_mapper,
                SLOT("map()"),
            )
            self.filenames_combobox_signal_mapper.setMapping(
                self.files_mapping[pattern_input_number]["filenames_combobox"],
                pattern_input_number,
            )
            self.connect(
                self.filenames_combobox_signal_mapper,
                SIGNAL("mapped(int)"),
                self.map_to_filenames_method,
            )
            self.files_mapping[pattern_input_number]["filename2filepath_mapping"] = {}

        self.setWindowTitle("select_multiple_files")
        self.setGeometry(250, 200, 350, 550)
        self.updating = False

        self.notify(node, ("input_modified", "all"))
        self.notify(node, ("caption_modified", node.get_caption()))

    def map_to_pattern_method(self, pattern_input_number):
        self.methods_mapping[pattern_input_number]["pattern_changed"](
            self, pattern_input_number
        )

    def map_to_filenames_method(self, pattern_input_number):
        self.methods_mapping[pattern_input_number]["filename_changed"](
            self, pattern_input_number
        )

    def notify(self, sender, event):
        """Update the widgets according to the notification sent by the node"""

        if event[0] == "caption_modified":
            self.window().setWindowTitle(event[1])

        if event[0] == "input_modified":
            if event[1] == 0:
                self.update_package_input()
                for pattern_number, methods in self.methods_mapping.items():
                    methods["update_filenames_combobox"](self, pattern_number)
                    methods["update_filepath_output"](self, pattern_number)
                self.package_lineedit.setCursorPosition(
                    self.cursorPosition_package_lineedit
                )
            elif event[1] <= self.nb_of_files:
                self.methods_mapping[event[1]]["update_pattern_input"](self, event[1])
                self.methods_mapping[event[1]]["update_filenames_combobox"](
                    self, event[1]
                )
                self.methods_mapping[event[1]]["update_filepath_output"](self, event[1])
                self.files_mapping[event[1]]["pattern_lineedit"].setCursorPosition(
                    self.files_mapping[event[1]]["pattern_lineedit_cursorPosition"]
                )
            elif event[1] < self.nb_of_inputs:
                pattern_number = event[1] - self.nb_of_files
                self.methods_mapping[pattern_number]["update_filepath_output"](
                    self, pattern_number
                )
            elif event[1] == "all":
                self.update_package_input()
                for pattern_number, methods in self.methods_mapping.items():
                    methods["update_pattern_input"](self, pattern_number)
                    methods["update_filenames_combobox"](self, pattern_number)
                    methods["update_filepath_output"](self, pattern_number)
                    self.files_mapping[pattern_number][
                        "pattern_lineedit"
                    ].setCursorPosition(
                        self.files_mapping[pattern_number][
                            "pattern_lineedit_cursorPosition"
                        ]
                    )
                self.package_lineedit.setCursorPosition(
                    self.cursorPosition_package_lineedit
                )

    @lock_notify
    def update_package_input(self):
        package_name = self.node.get_input(0)
        if self.updating or package_name is None:
            return
        self.updating = True
        self.package_lineedit.setText(package_name)
        self.updating = False

    def package_changed(self, package_name):
        self.cursorPosition_package_lineedit = self.package_lineedit.cursorPosition()
        self.node.set_input(0, str(package_name))


class select_multiple_files(Node):
    """Permits to select files which are available in a Python package and which
    names satisfy Unix style patterns.
    """

    def __init__(self, inputs, outputs):
        Node.__init__(self, inputs, outputs)
        self._concatenated_paths = None
        from openalea.file.files import PackageDir

        self._pd = PackageDir()

    def __call__(self, inputs):
        package_name = inputs[0]
        nb_of_inputs = self.get_nb_input()
        nb_of_files = (nb_of_inputs - 1) / 2
        filenames = inputs[nb_of_files + 1 : nb_of_inputs]
        self._concatenated_paths = []
        if package_name is not None:
            directory_path = path.path(self._pd([package_name])[0])
            for filename in filenames:
                if filename is not None and filename != "" and directory_path.isdir():
                    concatenated_path = directory_path / filename
                    if not concatenated_path.isfile():
                        concatenated_path = None
                else:
                    concatenated_path = None
                self._concatenated_paths.append(concatenated_path)
            return self._concatenated_paths
        else:
            return self._concatenated_paths
