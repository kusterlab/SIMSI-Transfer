from PyQt5 import QtWidgets

import sys
import logging
import time
import datetime

import simsi_transfer.main as simsi_transfer

import numpy as np


logging.basicConfig(level=logging.INFO)


class GuiLogger(logging.Handler):
    def emit(self, record):
        self.edit.append(self.format(record))  # implementation of append_line omitted
        QtWidgets.QApplication.processEvents()


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
        self.run_button.clicked.connect(self.process)
        self.run_button.setContentsMargins(20,100,20,100)
        
        layout.addRow(self.run_button)
    
    def _add_log_textarea(self, layout):    
        self.log_text_edit = QtWidgets.QTextEdit()
        
        self.logger = GuiLogger()
        self.logger.edit = self.log_text_edit
        logging.getLogger().addHandler(self.logger)
        
        layout.addRow(self.log_text_edit)
        
    def get_mq_txt_dir(self):
        mq_txt_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select MaxQuant combined/txt folder' , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.mq_txt_dir_line_edit.setText(mq_txt_dir)
    
    def get_raw_dir(self):
        raw_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select folder with RAW files' , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.raw_dir_line_edit.setText(raw_dir)
    
    def get_output_dir(self):
        output_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output folder' , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.output_dir_line_edit.setText(output_dir)
    
    def process(self):
        start = time.time()
        
        self.log_text_edit.append('Running...')
        
        self.log_text_edit.append('Transferring PSMs with SIMSI-Transfer')
        QtWidgets.QApplication.processEvents()
        self.run_simsi_transfer()
        
        total_time = round(time.time() - start)
        total_time = str(datetime.timedelta(seconds=total_time))
        
        self.log_text_edit.append('=== Done ===')
        QtWidgets.QApplication.processEvents()
        
        logging.info("SIMSI-Transfer completed in %s", total_time)
        
    def run_simsi_transfer(self):
        mq_txt_dir = self.mq_txt_dir_line_edit.text()
        raw_dir = self.raw_dir_line_edit.text()
        output_dir = self.output_dir_line_edit.text()
        extra_params = self.args_line_edit.text()
        simsi_transfer.main(['--mq_txt_folder', mq_txt_dir, '--raw_folder', raw_dir, '--output_folder', output_dir] + extra_params.split())
        

if __name__ == '__main__':
    if sys.platform.startswith('win'):
        # On Windows calling this function is necessary when combined with pyinstaller
        multiprocessing.freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()
