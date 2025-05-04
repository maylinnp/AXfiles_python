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
        e0 = []
        NaOH_concentration = []
        # will raise exception if invalid inputs
        calibration_files = self._process_inputs()
        # process first titration
        HCl_neutr_weight, E0_est, NaOH_conc_est = self.process_titration(
            calibration_files[0], first=True
        )
        e0.append(E0_est)
        NaOH_concentration.append(NaOH_conc_est)
        # process remaining titrations
        for file in calibration_files[1:]:
            HCl_neutr_weight, E0_est, NaOH_conc_est = self.process_titration(
                file, HCl_neutr_weight
            )
            e0.append(E0_est)
            NaOH_concentration.append(NaOH_conc_est)
        # disregard first titration
        NaOH_conc_mean = np.mean(NaOH_concentration[1:])
        NaOH_conc_std_percent = np.std(NaOH_concentration[1:]) / NaOH_conc_mean * 100
        logger.info(
            f"""
              The mean NaOH concentration estimated from this titration is:\n
              {NaOH_conc_mean:.5g} mol/kg-sol, with a standard deviation of
              +/-{NaOH_conc_std_percent:.2g} %.
              Update the value in the NaOH_summary file manually if needed."""
        )

    def process_titration(
        self, titration_file, HCl_neutr_weight: float = 0, first=False
    ):
        logger.debug(f"Processing file: {titration_file}")
        if first:
            titrant, sample, HCl_aliquot, titration_data = NaOH_calibration_data(
                titration_file
            )
            self.sample = sample
            self.titrant = titrant
        else:
            _, _, HCl_aliquot, titration_data = NaOH_calibration_data(
                titration_file, sample=self.sample, titrant=self.titrant
            )
        gran_function = (self.sample.w0 + HCl_aliquot.weight) * np.exp(
            titration_data.emf / self.sample.k
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
        NaOH_conc_est = (
            equivalence_weight
            * HCl_aliquot.conc
            * (HCl_aliquot.weight - HCl_neutr_weight)
        )
        logger.info(
            f"Calculated NaOH concentration from the first titration: {NaOH_conc_est} mol/kg-sol"
        )

        # calculate how much of the titrant in next round will be used to neutralize the excess NaOH
        HCl_neutr_weight = (
            (titration_data.weight[-1] - (-intercept / slope))
            * NaOH_conc_est
            / HCl_aliquot.conc
        )
        # estimate E0 from the data for better estimates next round(s)
        # E0 = E - k*log(concentration) at each titration point
        E0_est = np.mean(
            titration_data.emf[: len(F1_gran)]
            - self.sample.k
            * np.log(
                (HCl_aliquot.weight * HCl_aliquot.conc - F1_weight * NaOH_conc_est)
                / (self.sample.w0 + F1_weight)
            )
        )
        # new sample weight is original + total titrant added + acid aliquot
        self.sample.w0 = self.sample.w0 + HCl_aliquot.weight + titration_data.weight[-1]

        logger.info(f"Etimated E0 from the first titration is {E0_est} V")
        return HCl_neutr_weight, E0_est, NaOH_conc_est

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
