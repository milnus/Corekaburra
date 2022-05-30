import sys
import logging
import os
from logging import getLogger


def exit_with_error(message, exit_status, logger):
    """
    Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.
    :param message: Message to give the user upon exit
    :param exit_status: Status returned as exit status
    :param logger: Logger for program
    :return: None
    """

    logger.error(message)
    print(f"Corekaburra ERROR: {message}, exiting", file=sys.stderr)
    sys.exit(exit_status)
