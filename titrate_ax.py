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
from ax_maths import *
from functools import partial


from scipy.optimize import least_squares, root

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
        logger.info("Let's measure AX!!!\n")
        # initialize
        if os.path.isfile(path):
            logger.info(f"Processing single file {path}")
            self.file = path
            self.path = None
        elif os.path.isdir(path):
            logger.info(f"Processing all AX files in direcyory {path}")
            self.file = None
            self.path = path
        else:
            raise FileNotFoundError(f"Check the provided path: {path}")

    def titrate(self):
        if self.file:
            titration_files = [self.file]
        else:
            titration_files = self._process_inputs()
        logger.info(f"Number of files slated for processing: {len(titration_files)}")

        for file in titration_files:
            self.process_titration(file)

    def process_titration(self, file: str):
        logger.debug(f"Processing this file now: {file}.")
        sample, HCl_titration_data, NaOH_titration_data = titration_data(file)
        logger.debug(f"{file} successfully parsed.")

        # nutrients and constants already in Sample()
        # CT after degas also in sample
        # find index for good fwd titration data
        HCl_titr_good_indices = find_data_in_range(3, 3.5, HCl_titration_data.pH_est)
        # estimate E0 and AT from fwd
        HCl_mass = HCl_titration_data.weight[HCl_titr_good_indices]
        emf = HCl_titration_data.emf[HCl_titr_good_indices]
        T = np.mean(HCl_titration_data.T[HCl_titr_good_indices])
        try:
            AT_est_fwd, E0_est_fwd = estimate_AT_E0(
                HCl_mass, emf, T, sample.m0, HCl_titration_data.titrant.concentration
            )
            logger.info(
                f"Estimated total alkalinity: {AT_est_fwd*1e6:.2f} umol/kg and E0: {E0_est_fwd:.4f} V"
            )
        except:
            logger.warning("Not enough data in the forward titration")
            AT_est_fwd = None
            E0_est_fwd = None

        if AT_est_fwd:
            result = least_squares(
                fun=partial(AT_residuals, sample=sample, titration=HCl_titration_data),
                x0=[1, AT_est_fwd],
                jac=partial(AT_jacobian, sample=sample, titration=HCl_titration_data),
                method="lm",
                xtol=1e-15,
                ftol=1e-15,
                gtol=1e-15,
            )
            # The result is a bit higher than the matlab function, needs more optimization
            # TODO might be issue with my constants, check solution classes
            f_fwd, AT_fwd = result.x
            E0_fwd = E0_est_fwd - k_boltz(T) * log(f_fwd)
            print(f"f = {f_fwd:.6f}, AT = {AT_fwd*1e6:.6f}")

        ## Back titration
        if E0_fwd:
            NaOH_titration_data.recalculate_pH(E0_fwd)
        # re-estimate NaOH pH data from fwd titration
        NaOH_titr_good_indices = find_data_in_range(3, 3.5, NaOH_titration_data.pH_est)

        # re-estimate the NaOH data with E0 from bwd low pH titr
        NaOH_high_pH_data = find_data_in_range(9, 10.5, NaOH_titration_data.pH_est)

        idx2 = "AT titration range for bwd 3-3.5"
        idx3 = "KW titration range during bwd, 9-10.5"

    def fwd_titration(self, titration_data: Titration, sample):
        # This method should essentially give you the AT from data processing
        # described in Dickson et al. 2003
        pass

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
