
import sys
import logging


def setup_logging(verbosity: int = 0) -> None:
    """
    Configure logging based on verbosity level.
    Args:
        verbosity: 0=ERROR, 1=WARNING, 2=INFO, 3+=DEBUG
    """

    ERROR_CODE   = 0
    WARNING_CODE = 1
    INFO_CODE    = 2
    DEBUG_CODE   = 3

    # Cap at DEBUG for higher verbosity
    verbosity = min(verbosity, DEBUG_CODE)

    # Map verbosity to log levels
    level_map = {
        ERROR_CODE:   logging.ERROR,   # default: only errors
        WARNING_CODE: logging.WARNING, # -v: warnings and above
        INFO_CODE:    logging.INFO,    # -vv: info and above
        DEBUG_CODE:   logging.DEBUG,   # -vvv: debug and above
    }

    log_level = level_map.get(
        verbosity,
        logging.DEBUG
    )

    # Configure root logger
    datefmt = '%Y-%m-%d %H:%M:%S'
    if verbosity < DEBUG_CODE:
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt=datefmt,
            stream=sys.stderr
        )
    else:
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(filename)s:%(funcName)s - %(message)s',
            datefmt=datefmt,
            stream=sys.stderr
        )
    # end if
# end def
