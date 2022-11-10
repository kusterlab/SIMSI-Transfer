from PyQt5 import QtWidgets
from PyQt5.QtCore import QObject, QThread, pyqtSignal, Qt  # https://realpython.com/python-pyqt-qthread/

import sys
import logging
import time
import multiprocessing
from logging.handlers import QueueListener

import simsi_transfer.main as simsi_transfer
import simsi_transfer.utils.multiprocessing_pool as pool

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def run_simsi_transfer(mq_txt_dir, raw_dir, output_dir, metafile_path, tmt_params, other_params, extra_params):
    parameters = []
    if mq_txt_dir:
        parameters.extend(['--mq_txt_folder', mq_txt_dir])
    if raw_dir:
        parameters.extend(['--raw_folder', raw_dir])
    if output_dir:
        parameters.extend(['--output_folder', output_dir])
    if metafile_path:
        parameters.extend(['--meta_input_file', metafile_path])
    if tmt_params:
        parameters.extend(tmt_params)
    if other_params:
        parameters.extend(other_params)
    if extra_params:
        parameters.extend(extra_params.split())
    try:
        simsi_transfer.main(parameters)
        # logger.info(parameters)
    except SystemExit as e:
        logger.info(f"Error while running SIMSI-Transfer, exited with error code {e}.")
    except Exception as e:
        logger.info(f"Error while running SIMSI-Transfer: {e}")


# https://stackoverflow.com/questions/28655198/best-way-to-display-logs-in-pyqt#60528393
class QTextEditLogger(logging.Handler, QObject):
    appendPlainText = pyqtSignal(str)

    def __init__(self, parent):
        # initialize logging.Handler
        super().__init__()
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        self.setFormatter(formatter)

        # initialize QObject
        QObject.__init__(self)
        self.widget = QtWidgets.QPlainTextEdit(parent)
        self.widget.setReadOnly(True)
        self.appendPlainText.connect(self.widget.appendPlainText)

    def emit(self, record):
        self.appendPlainText.emit(self.format(record))


# https://stackoverflow.com/questions/53288877/python-multiprocessing-sending-child-process-logging-to-gui-running-in-parent
class LogEmitter(QObject):
    sigLog = pyqtSignal(str)


class LogHandler(logging.Handler):
    def __init__(self):
        super().__init__()
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        self.setFormatter(formatter)
        self.emitter = LogEmitter()

    def emit(self, record):
        msg = self.format(record)
        self.emitter.sigLog.emit(msg)


class FileSelect(QtWidgets.QWidget):
    def __init__(self, file_type, file_extensions, file_hint='', folder_select=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file_type = file_type
        self.file_hint = file_hint
        self.file_extensions = file_extensions

        self.label_text = f"Select {self.file_type} file"
        if folder_select:
            self.label_text = f"Select {self.file_type} folder"

        self.file_hint_text = ""
        if len(file_hint) > 0:
            self.file_hint_text = f'<br><font color="grey">{self.file_hint}</font>'
        self.label = QtWidgets.QLabel(self.label_text + self.file_hint_text)

        self.hbox_layout = QtWidgets.QHBoxLayout()
        self.hbox_layout.setContentsMargins(0, 0, 0, 0)

        self.line_edit = QtWidgets.QLineEdit()
        self.browse_button = QtWidgets.QPushButton("Browse")
        if folder_select:
            self.browse_button.clicked.connect(self.select_dir)
        else:
            self.browse_button.clicked.connect(self.select_file)
        self.hbox_layout.addWidget(self.line_edit, stretch=1)
        self.hbox_layout.addWidget(self.browse_button)

        self.setLayout(self.hbox_layout)

    def select_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, self.label_text, '', self.file_extensions)
        self.line_edit.setText(filename)

    def select_dir(self):
        output_dir = QtWidgets.QFileDialog.getExistingDirectory(self, self.label_text, '',
                                                                QtWidgets.QFileDialog.ShowDirsOnly)
        self.line_edit.setText(output_dir)

    def get_file(self):
        return self.line_edit.text()

    def setButtonsEnabled(self, enable):
        self.browse_button.setEnabled(enable)


class TMTGroup(QtWidgets.QGroupBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.tmt_group_layout = QtWidgets.QGridLayout()
        self.setLayout(self.tmt_group_layout)

        self.ms_level_label = QtWidgets.QLabel("TMT MS level")
        # self.ms_level_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.ms_level_select = QtWidgets.QComboBox()
        self.ms_level_select.addItems(['MS2', 'MS3'])
        self.ms_level_select.setCurrentText('MS3')

        self.requantify_label = QtWidgets.QLabel("Requantify TMT?")
        self.requantify_checkbox = QtWidgets.QCheckBox()

        self.tmt_group_layout.addWidget(self.ms_level_label, 0, 0)
        self.tmt_group_layout.addWidget(self.ms_level_select, 0, 1)
        self.tmt_group_layout.addWidget(self.requantify_label, 0, 3)
        self.tmt_group_layout.addWidget(self.requantify_checkbox, 0, 4)

        for col in range(5):
            self.tmt_group_layout.setColumnStretch(col, 1)

    def get_params(self):
        returnval = ["--tmt_ms_level", str(self.ms_level_select.currentText()).lower()]
        if self.requantify_checkbox.isChecked():
            returnval.append("--tmt_requantify")
        return returnval


class ParameterGroup(QtWidgets.QGroupBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.parameter_group_layout = QtWidgets.QGridLayout()
        self.setLayout(self.parameter_group_layout)

        self.filter_decoy_label = QtWidgets.QLabel("Remove decoys")
        self.filter_decoy_checkbox = QtWidgets.QCheckBox()

        # self.plotting_columns_label = QtWidgets.QLabel("Keep plotting columns")
        # self.plotting_columns_checkbox = QtWidgets.QCheckBox()

        self.ambiguity_label = QtWidgets.QLabel("PSM ambiguity")
        self.ambiguity_select = QtWidgets.QComboBox()
        self.ambiguity_select.addItems(['majority', 'all', 'none'])
        self.ambiguity_select.setCurrentText('majority')

        self.threads_label = QtWidgets.QLabel("CPU threads")
        self.threads_spinbox = QtWidgets.QSpinBox()
        self.threads_spinbox.setValue(1)
        self.threads_spinbox.setRange(1, 100)

        self.stringency_label = QtWidgets.QLabel("MaRaCluster stringencies")
        self.stringency_line = QtWidgets.QLineEdit('10, 15, 20')

        self.parameter_group_layout.addWidget(self.filter_decoy_label, 0, 0)
        self.parameter_group_layout.addWidget(self.filter_decoy_checkbox, 0, 1)
        # self.parameter_group_layout.addWidget(self.plotting_columns_label, 0, 3)
        # self.parameter_group_layout.addWidget(self.plotting_columns_checkbox, 0, 4)
        self.parameter_group_layout.addWidget(self.stringency_label, 0, 3)
        self.parameter_group_layout.addWidget(self.stringency_line, 0, 4)

        self.parameter_group_layout.addWidget(self.ambiguity_label, 1, 0)
        self.parameter_group_layout.addWidget(self.ambiguity_select, 1, 1)

        self.parameter_group_layout.addWidget(self.threads_label, 1, 3)
        self.parameter_group_layout.addWidget(self.threads_spinbox, 1, 4)

        for col in range(5):
            self.parameter_group_layout.setColumnStretch(col, 1)

    def get_params(self):
        returnval = [
            '--stringencies', str(self.stringency_line.text()),
            '--num_threads', str(self.threads_spinbox.value()),
            '--ambiguity_decision', str(self.ambiguity_select.currentText())
        ]
        if self.filter_decoy_checkbox.isChecked():
            returnval.append('--filter_decoys')
        # if self.plotting_columns_checkbox.isChecked():
        #     returnval.append('--add_plotting_columns')
        return returnval

class MainWindow(QtWidgets.QWidget):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        layout = QtWidgets.QFormLayout()

        self.setWindowTitle("SIMSI-Transfer")


        self.tabs = QtWidgets.QTabWidget()
        self._add_single_search_input_tab()
        self._add_metafile_input_tab()

        layout.addRow(self.tabs)
        self._add_output_dir_field(layout)

        self.tmt_group = TMTGroup('TMT parameters')
        self.parameter_group = ParameterGroup('SIMSI parameters')

        layout.addRow(self.tmt_group)
        layout.addRow(self.parameter_group)

        self._add_extra_params_field(layout)
        self._add_buttons(layout)
        self._add_log_textarea(layout)

        self.setLayout(layout)

        # sets up handler that will be used by QueueListener
        # which will update the LogDialoag
        handler = LogHandler()
        handler.emitter.sigLog.connect(self.log_text_area.widget.appendPlainText)

        self.q = multiprocessing.Queue()
        self.ql = QueueListener(self.q, handler)
        self.ql.start()

        self.pool = pool.JobPool(processes=1, warningFilter="default", queue=self.q)

        self.resize(700, self.height())

    def _add_buttons(self, layout):
        self.help_button = QtWidgets.QPushButton("Help")
        self.help_button.clicked.connect(self.run_simsi_help)

        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.clicked.connect(self.run_simsi_transfer)

        layout.addRow(self.help_button, self.run_button)

    def _add_single_search_input_tab(self):
        self.singlesearch_tab = QtWidgets.QWidget()

        self.singlesearch_layout = QtWidgets.QFormLayout()
        self._add_mq_txt_dir_field(self.singlesearch_layout)
        self._add_raw_dir_field(self.singlesearch_layout)

        self.singlesearch_tab.setLayout(self.singlesearch_layout)
        self.tabs.addTab(self.singlesearch_tab, "Single search input")

    def _add_metafile_input_tab(self):
        self.metafile_tab = QtWidgets.QWidget()
        self.metafile_widget = FileSelect('metafile', 'Comma-separated file (*.txt *.csv *.tsv)')

        self.metafile_layout = QtWidgets.QFormLayout()
        self.metafile_layout.addRow(self.metafile_widget.label, self.metafile_widget)
        self.metafile_tab.setLayout(self.metafile_layout)
        self.tabs.addTab(self.metafile_tab, "Metafile input")

    def _add_mq_txt_dir_field(self, layout):
        # evidence.txt input
        self.mq_txt_dir_label = QtWidgets.QLabel("Select MaxQuant combined/txt folder")
        # self.mq_txt_dir_label.setMargin(10)

        self.mq_txt_dir_widget = QtWidgets.QWidget()

        self.mq_txt_dir_hbox_layout = QtWidgets.QHBoxLayout()
        self.mq_txt_dir_hbox_layout.setContentsMargins(0, 0, 0, 0)

        self.mq_txt_dir_line_edit = QtWidgets.QLineEdit()
        self.mq_txt_dir_browse_button = QtWidgets.QPushButton("Browse")
        self.mq_txt_dir_browse_button.clicked.connect(self.get_mq_txt_dir)
        self.mq_txt_dir_hbox_layout.addWidget(self.mq_txt_dir_line_edit, stretch=1)
        self.mq_txt_dir_hbox_layout.addWidget(self.mq_txt_dir_browse_button)

        self.mq_txt_dir_widget.setLayout(self.mq_txt_dir_hbox_layout)

        layout.addRow(self.mq_txt_dir_label, self.mq_txt_dir_widget)

    def _add_raw_dir_field(self, layout):
        # fasta file input
        self.raw_dir_label = QtWidgets.QLabel("Select folder with RAW files")
        # self.raw_dir_label.setMargin(10)

        self.raw_dir_widget = QtWidgets.QWidget()

        self.raw_dir_hbox_layout = QtWidgets.QHBoxLayout()
        self.raw_dir_hbox_layout.setContentsMargins(0, 0, 0, 0)

        self.raw_dir_line_edit = QtWidgets.QLineEdit()
        self.raw_dir_browse_button = QtWidgets.QPushButton("Browse")
        self.raw_dir_browse_button.clicked.connect(self.get_raw_dir)
        self.raw_dir_hbox_layout.addWidget(self.raw_dir_line_edit, stretch=1)
        self.raw_dir_hbox_layout.addWidget(self.raw_dir_browse_button)

        self.raw_dir_widget.setLayout(self.raw_dir_hbox_layout)

        layout.addRow(self.raw_dir_label, self.raw_dir_widget)

    def _add_output_dir_field(self, layout):
        # fasta file input
        self.output_dir_label = QtWidgets.QLabel("Select output folder")
        # self.output_dir_label.setMargin(10)

        self.output_dir_widget = QtWidgets.QWidget()

        self.output_dir_hbox_layout = QtWidgets.QHBoxLayout()
        self.output_dir_hbox_layout.setContentsMargins(0, 0, 0, 0)

        self.output_dir_line_edit = QtWidgets.QLineEdit()
        self.output_dir_browse_button = QtWidgets.QPushButton("Browse")
        self.output_dir_browse_button.clicked.connect(self.get_output_dir)
        self.output_dir_hbox_layout.addWidget(self.output_dir_line_edit, stretch=1)
        self.output_dir_hbox_layout.addWidget(self.output_dir_browse_button)

        self.output_dir_widget.setLayout(self.output_dir_hbox_layout)

        layout.addRow(self.output_dir_label, self.output_dir_widget)

    def _add_extra_params_field(self, layout):
        self.args_label = QtWidgets.QLabel("Additional parameters")

        self.args_line_edit = QtWidgets.QLineEdit()

        layout.addRow(self.args_label, self.args_line_edit)

    def _add_log_textarea(self, layout):
        self.log_text_area = QTextEditLogger(self)
        self.log_text_area.setLevel(logging.INFO)
        logger.addHandler(self.log_text_area)

        layout.addRow(self.log_text_area.widget)

    def get_mq_txt_dir(self):
        mq_txt_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select MaxQuant combined/txt folder', '',
                                                                QtWidgets.QFileDialog.ShowDirsOnly)
        self.mq_txt_dir_line_edit.setText(mq_txt_dir)

    def get_metafile_path(self):
        metafile_path = QtWidgets.QFileDialog.getOpenFileName(self, 'Select path to metafile', '')

    def get_raw_dir(self):
        raw_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select folder with RAW files', '',
                                                             QtWidgets.QFileDialog.ShowDirsOnly)
        self.raw_dir_line_edit.setText(raw_dir)

    def get_output_dir(self):
        output_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output folder', '',
                                                                QtWidgets.QFileDialog.ShowDirsOnly)
        self.output_dir_line_edit.setText(output_dir)

    def set_buttons_enabled_state(self, enable):
        self.mq_txt_dir_browse_button.setEnabled(enable)
        self.raw_dir_browse_button.setEnabled(enable)
        self.output_dir_browse_button.setEnabled(enable)
        # self.run_button.setEnabled(enable)
        # Cannot stop a QThread if it doesn't have an own event loop
        self.run_button.clicked.disconnect()
        if enable:
            self.run_button.setText("Run")
            self.run_button.clicked.connect(self.run_simsi_transfer)
        else:
            self.run_button.setText("Stop")
            self.run_button.clicked.connect(self.stop_simsi_transfer)

    def run_simsi_help(self):
        import webbrowser
        webbrowser.open('https://github.com/kusterlab/SIMSI-Transfer')
        # self.pool.applyAsync(run_simsi_transfer, ('', '', '', '', [], [], '--help'), callback=self.on_simsi_finished)

    def run_simsi_transfer(self):
        # initialize all parameters as empty values and then override?
        mq_txt_dir = ''
        raw_dir = ''
        metafile_path = False
        if self.tabs.currentIndex() == 0:
            mq_txt_dir = self.mq_txt_dir_line_edit.text()
            raw_dir = self.raw_dir_line_edit.text()
        else:
            metafile_path = self.metafile_widget.get_file()

        output_dir = self.output_dir_line_edit.text()
        tmt_params = self.tmt_group.get_params()
        other_params = self.parameter_group.get_params()

        extra_params = self.args_line_edit.text()

        self.set_buttons_enabled_state(False)
        self.pool.applyAsync(run_simsi_transfer,
                             (mq_txt_dir, raw_dir, output_dir, metafile_path, tmt_params, other_params, extra_params),
                             callback=self.on_simsi_finished)

    def on_simsi_finished(self, return_code):
        self.set_buttons_enabled_state(True)

    def stop_simsi_transfer(self):
        self.pool.stopPool()
        self.on_simsi_finished(-2)

        logger.info("SIMSI-Transfer stopped by user")

        self.pool = pool.JobPool(processes=1, warningFilter="default", queue=self.q)

    def closeEvent(self, _):
        self.stop_simsi_transfer()


if __name__ == '__main__':
    if sys.platform.startswith('win'):
        # On Windows calling this function is necessary when combined with pyinstaller: https://stackoverflow.com/questions/24944558/pyinstaller-built-windows-exe-fails-with-multiprocessing
        multiprocessing.freeze_support()
    else:
        # On Linux calling this function is necessary when combined with pyqt: https://stackoverflow.com/questions/29556291/multiprocessing-with-qt-works-in-windows-but-not-linux
        multiprocessing.set_start_method('spawn')
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()
