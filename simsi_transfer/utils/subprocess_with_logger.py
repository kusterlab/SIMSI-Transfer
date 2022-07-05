import sys
import logging
import subprocess

logger = logging.getLogger(__name__)


def run_simple(cmd):
    logger.info(f"Running subprocess with command: {cmd}")
    process = subprocess.run(
        cmd.split(),
        shell=True,
        check=True)


def run(cmd):
    logger.info(f"Running subprocess with command: {cmd}")
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        shell=True)
    with process.stdout:
        log_subprocess_output(process.stdout)
    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Issues encountered while running subprocess with command\nReturn code: {process.returncode}\nCommand: {cmd}")


def log_subprocess_output(pipe):
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        logger.info(line.decode("utf-8").rstrip())
