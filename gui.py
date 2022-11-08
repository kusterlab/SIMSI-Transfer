from PyQt5 import QtWidgets
from PyQt5.QtCore import QObject, QThread, pyqtSignal # https://realpython.com/python-pyqt-qthread/

import sys
import logging
import time
import multiprocessing
from logging.handlers import QueueListener

import simsi_transfer.main as simsi_transfer
import simsi_transfer.utils.multiprocessing_pool as pool


logger = logging.getLogger()
logger.setLevel(logging.INFO)


def run_simsi_transfer(mq_txt_dir, raw_dir, output_dir, extra_params):
    try:
        simsi_transfer.main(['--mq_txt_folder', mq_txt_dir, '--raw_folder', raw_dir, '--output_folder', output_dir] + extra_params.split())
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


class MainWindow(QtWidgets.QWidget):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        layout = QtWidgets.QFormLayout()
        
        self.setWindowTitle("SIMSI-Transfer")
        
        self._add_mq_txt_dir_field(layout)
        self._add_raw_dir_field(layout)
        self._add_output_dir_field(layout)
        self._add_extra_params_field(layout)
        self._add_run_button(layout)
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

    def _add_mq_txt_dir_field(self, layout):
        # evidence.txt input
        self.mq_txt_dir_label = QtWidgets.QLabel("Select MaxQuant combined/txt folder")
        #self.mq_txt_dir_label.setMargin(10)
        
        self.mq_txt_dir_widget = QtWidgets.QWidget()
        
        self.mq_txt_dir_hbox_layout = QtWidgets.QHBoxLayout()
        self.mq_txt_dir_hbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.mq_txt_dir_line_edit = QtWidgets.QLineEdit()
        self.mq_txt_dir_browse_button = QtWidgets.QPushButton("Browse")
        self.mq_txt_dir_browse_button.clicked.connect(self.get_mq_txt_dir)
        self.mq_txt_dir_hbox_layout.addWidget(self.mq_txt_dir_line_edit, stretch = 1)
        self.mq_txt_dir_hbox_layout.addWidget(self.mq_txt_dir_browse_button)
        
        self.mq_txt_dir_widget.setLayout(self.mq_txt_dir_hbox_layout)
        
        layout.addRow(self.mq_txt_dir_label, self.mq_txt_dir_widget)

    def _add_raw_dir_field(self, layout):
        # fasta file input
        self.raw_dir_label = QtWidgets.QLabel("Select folder with RAW files")
        #self.raw_dir_label.setMargin(10)
        
        self.raw_dir_widget = QtWidgets.QWidget()
        
        self.raw_dir_hbox_layout = QtWidgets.QHBoxLayout()
        self.raw_dir_hbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.raw_dir_line_edit = QtWidgets.QLineEdit()
        self.raw_dir_browse_button = QtWidgets.QPushButton("Browse")
        self.raw_dir_browse_button.clicked.connect(self.get_raw_dir)
        self.raw_dir_hbox_layout.addWidget(self.raw_dir_line_edit, stretch = 1)
        self.raw_dir_hbox_layout.addWidget(self.raw_dir_browse_button)
        
        self.raw_dir_widget.setLayout(self.raw_dir_hbox_layout)
        
        layout.addRow(self.raw_dir_label, self.raw_dir_widget)
    
    def _add_output_dir_field(self, layout):
        # fasta file input
        self.output_dir_label = QtWidgets.QLabel("Select output folder")
        #self.output_dir_label.setMargin(10)
        
        self.output_dir_widget = QtWidgets.QWidget()
        
        self.output_dir_hbox_layout = QtWidgets.QHBoxLayout()
        self.output_dir_hbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.output_dir_line_edit = QtWidgets.QLineEdit()
        self.output_dir_browse_button = QtWidgets.QPushButton("Browse")
        self.output_dir_browse_button.clicked.connect(self.get_output_dir)
        self.output_dir_hbox_layout.addWidget(self.output_dir_line_edit, stretch = 1)
        self.output_dir_hbox_layout.addWidget(self.output_dir_browse_button)
        
        self.output_dir_widget.setLayout(self.output_dir_hbox_layout)
        
        layout.addRow(self.output_dir_label, self.output_dir_widget)

    def _add_extra_params_field(self, layout):
        # additional parameters input, TODO: make user friendly options for each important parameter
        self.args_label = QtWidgets.QLabel("Additional parameters")
        #self.args_label.setMargin(10)
        
        self.args_line_edit = QtWidgets.QLineEdit()
        
        layout.addRow(self.args_label, self.args_line_edit)
    
    def _add_run_button(self, layout):    
        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.clicked.connect(self.run_simsi_transfer)
        self.run_button.setContentsMargins(20,100,20,100)
        
        layout.addRow(self.run_button)
    
    def _add_log_textarea(self, layout):    
        self.log_text_area = QTextEditLogger(self)
        self.log_text_area.setLevel(logging.INFO)
        logger.addHandler(self.log_text_area)
             
        layout.addRow(self.log_text_area.widget)
        
    def get_mq_txt_dir(self):
        mq_txt_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select MaxQuant combined/txt folder' , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.mq_txt_dir_line_edit.setText(mq_txt_dir)
    
    def get_raw_dir(self):
        raw_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select folder with RAW files' , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.raw_dir_line_edit.setText(raw_dir)
    
    def get_output_dir(self):
        output_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output folder' , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.output_dir_line_edit.setText(output_dir)
    
    def set_buttons_enabled_state(self, enable):
        self.mq_txt_dir_browse_button.setEnabled(enable)
        self.raw_dir_browse_button.setEnabled(enable)
        self.output_dir_browse_button.setEnabled(enable)
        #self.run_button.setEnabled(enable)
        # Cannot stop a QThread if it doesn't have an own event loop
        self.run_button.clicked.disconnect()
        if enable:
            self.run_button.setText("Run")
            self.run_button.clicked.connect(self.run_simsi_transfer)
        else:
            self.run_button.setText("Stop")
            self.run_button.clicked.connect(self.stop_simsi_transfer)
        
    def run_simsi_transfer(self):
        mq_txt_dir = self.mq_txt_dir_line_edit.text()
        raw_dir = self.raw_dir_line_edit.text()
        output_dir = self.output_dir_line_edit.text()
        extra_params = self.args_line_edit.text()
        
        self.set_buttons_enabled_state(False)
        self.pool.applyAsync(run_simsi_transfer, (mq_txt_dir, raw_dir, output_dir, extra_params), callback=self.on_simsi_finished)
    
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
