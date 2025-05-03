import argparse
from exceptions import CalibrationDataMissing
import os, sys
from util import get_matching_files
from extract_data import NaOH_calibration_data
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
    help="path to NaOH calibration data files",
    default=argparse.SUPPRESS,
)
parser.add_argument(
    "-id",
    "--titration_id",
    help="batch identifier that will uniquely identify all files that are used to calculate avg cNaOH",
    default=argparse.SUPPRESS,
)
parser.add_argument(
    "-ext",
    "--file_extension",
    help="optional file extension if not using csv",
    default="csv",
)
parser.add_argument(
    "-hcl",
    "--hcl_concentration",
    help="optional HCl concentration, if e.g., incorrect information was entered during calibration (or missing). Will be treated as mol/kg-sol",
)

args = parser.parse_args()


class CalibrateNaOH:
    def __init__(
        self,
        path: str,
        titration_id: str,
        hcl_concentration: float = None,
        file_extension: str = None,
    ):
        logger.info("Let's get calibrating!\n\n")
        # initialize
        self.path = path
        self.titration_id = titration_id
        self.hcl_concentration = hcl_concentration
        self.file_extension = file_extension

    def calibrate(self):
        # will raise exception if invalid inputs
        calibration_files = self._process_inputs()
        self.calculate_first_titr_data(calibration_files[0])

    def calculate_first_titr_data(self, titration_file):
        logger.debug(f"First file: {titration_file}")
        titrant, sample, HCl_aliquot, titration_data = NaOH_calibration_data(
            titration_file
        )
        gran_function = (sample.w0 + HCl_aliquot.weight) * np.exp(
            titration_data.emf / sample.k
        )
        # grab only good data
        F1_gran = gran_function[: np.count_nonzero(gran_function > 100)]
        F1_weight = titration_data.weight[: len(F1_gran)]
        logger.info(
            f"Gran data in valid range uses {len(F1_gran)} of the total {len(gran_function)} points from the titration"
        )
        slope, intercept, r_value, p_value, std_err = linregress(F1_weight, F1_gran)
        goodness_of_fit = r_value**2
        equivalence_weight = -slope / intercept
        NaOH_concentration = (
            equivalence_weight * HCl_aliquot.concentration * HCl_aliquot.weight
        )
        logger.info(
            f"Calculated NaOH concentration from the first titration: {NaOH_concentration} mol/kg-sol"
        )
        # calculate how far past equivalence point it was titrated, which will be subtracted from the next titration
        NaOH_excess = (
            titration_data.weight[-1] - (-intercept / slope)
        ) * NaOH_concentration
        # calculate how much of the titrant in next round will be used to neutralize the excess NaOH
        HCl_neutr_weight = NaOH_excess / HCl_aliquot.concentration
        # estimate E0 from the data for better estimates next round(s)
        # E0 = E - k*log(concentration)
        E0_est = np.mean(
            titration_data.emf[: len(F1_gran)]
            - sample.k
            * np.log(
                (
                    HCl_aliquot.weight * HCl_aliquot.concentration
                    - F1_weight * NaOH_concentration
                )
                / (sample.w0 + F1_weight)
            )
        )
        logger.info(f"Etimated E0 from the first titration is {E0_est} V")

    def _process_inputs(self) -> list:
        """
        Checks inputs and makes list of files specified by
        input from path, pattern, and extension

        Raises:
            CalibrationDataMissing

        Returns:
            list: valid files
        """
        if self.hcl_concentration:
            logger.info(
                "Ready to overwrite file HCl concentration with this: ",
                self.hcl_concentration,
            )

        calibration_files = get_matching_files(
            self.path, self.titration_id, self.file_extension
        )
        # catch if only one calibration file
        if not len(calibration_files) > 1:
            raise CalibrationDataMissing(
                "Need at least two calibration files to proceed, please check your inputs."
            )
        else:
            return calibration_files


try:
    calibration = CalibrateNaOH(**vars(args))
except TypeError as e:
    logger.critical(f"Error in command line inputs: {e}")
    sys.exit(1)

calibration.calibrate()
