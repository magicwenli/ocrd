import logging
import sys
from logging.handlers import TimedRotatingFileHandler

FORMATTER = logging.Formatter(
    "%(asctime)s - %(name)s - [%(levelname)s] - %(message)s")
LOG_FILE = "my_app.log"

DEBUG = logging.DEBUG
INFO = logging.INFO
WARNING = logging.WARNING
ERROR = logging.ERROR
CRITICAL = logging.CRITICAL


def get_console_handler():
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(FORMATTER)
    return console_handler


def get_file_handler(file):
    if file is True:
        file_handler = TimedRotatingFileHandler(LOG_FILE, when='midnight')
    else:
        file_handler = TimedRotatingFileHandler(file, when='midnight')
    file_handler.setFormatter(FORMATTER)
    return file_handler


def get_logger(logger_name, log_level=DEBUG, file_handler=False):
    logger = logging.getLogger(logger_name)

    logger.setLevel(log_level)  # better to have too much log than not enough

    logger.addHandler(get_console_handler())
    if(file_handler):
        logger.addHandler(get_file_handler(file_handler))

    # with this pattern, it's rarely necessary to propagate the error up to parent
    logger.propagate = False

    return logger
