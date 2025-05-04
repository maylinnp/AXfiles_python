import argparse
from exceptions import TitrationDataMissing
import os, sys
from util import get_matching_files
from extract_data import titration_data
import logging
from solutions import *
import math
import numpy as np
from scipy.stats import linregress

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)

# Create a formatter
formatter = logging.Formatter("%(levelname)s - %(message)s")

# Set the formatter for the stream handler
stream_handler.setFormatter(formatter)

# Add the stream handler to the logger
logger.addHandler(stream_handler)


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    "-p",
    "--path",
    help="path to one file or all files in a folder",
    default=argparse.SUPPRESS,
)
args = parser.parse_args()


class TitrateAX:
    def __init__(
        self,
        path: str = None,
    ):
        logger.info("Let's measure AX!!!\n\n")
        # initialize
        if os.path.isfile(path):
            logger.info(f"Processing single file {path}")
            self.file = path
            self.path = None
        elif os.path.isdir(path):
            logger.info(f"Processing all AX files in direcyory {path}")
            self.file = None
            self.path = path

    def titrate(self):
        if self.file:
            titration_files = [self.file]
        else:
            titration_files = self._process_inputs()
        logger.info(f"Number of files slated for processing: {len(titration_files)}")

        for file in titration_files:
            self.process_titration(file)

    def process_titration(self, file: str):
        logger.debug(f"Processing this file now: {file}")
        sample, HCl_titration_data, NaOH_titration_data = titration_data(file)
        logger.debug(f"Processing this file now: {file}")
        print(NaOH_titration_data.titrant.name)

    def _process_inputs(self) -> list:
        """
        Checks inputs and makes list of files specified by
        input from path, pattern, and extension

        Raises:
            TitrationDataMissing

        Returns:
            list: valid files
        """
        titration_files = get_matching_files(self.path, "", "csv")
        # catch if only one calibration file
        if not titration_files:
            raise TitrationDataMissing("No titration data files in the provided path.")
        else:
            return titration_files


try:
    titration = TitrateAX(**vars(args))
except TypeError as e:
    logger.critical(f"Error in command line inputs: {e}")
    sys.exit(1)

titration.titrate()
