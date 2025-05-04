## Extracts data from AX titrations
# Data is saved in one structure per file
# Created 2020/09/08 MLP based on old extract-data codes
# Updated 2020/11/22 MLP: I have edited the Labview code, so that it
# applies volume correction per every 5 mL.  I also edited it so that on
# each line, the volume listed is what was added, waited 15 s, then pH on
# the corresponding line recorded (i.e., volume and pH value correspond to
# eachother).
import csv
from exceptions import *
from typing import Union
import pandas as pd
import logging
from solutions import *
import statistics
from util import get_matching_files

logger = logging.getLogger(__name__)

logger.setLevel("DEBUG")


# class Solution(object):
#     pass


# TODO make a class that has all these methods in the generic form,
# Then child classes that have extra specific methods
# and use other custom classes to hold the final data so that they don't have anyhting else
def titration_data(filename: str, burette_id: str = "dosimat 12"):
    # assumes calibration solution type found in file name
    if "nacl" in filename.lower():
        sample = NaCl()
    elif "kcl" in filename.lower():
        sample = KCl()
    else:
        sample = SW()
    # # Open filename and extract data
    with open(filename, "r") as datafile:
        csvreader = csv.reader(datafile)
        sample_info = next(csvreader)
        sample.w0 = float(sample_info[0]) / 1000
        sample.salt = float(sample_info[1])
        sample.emf0 = float(sample_info[2])
        t0 = float(sample_info[3])
        HCl_id = sample_info[7].split("-")[0]
        NaOH_id = sample_info[6].split("-")[0]

        # check if FWD
        fwd_data = list()
        for row in csvreader:
            # check if bwd
            if "BWD" in row:
                bwd_data = list()
                break
            else:
                fwd_data.append(row)
        if "bwd_data" in locals():
            for row in csvreader:
                bwd_data.append(row)

    # flag sample as Q (questionable) if temperature very out of range
    if t0 < 15 or t0 > 30:
        sample.flag = "Q"

    # check if batch and concentration agree
    titrant.concentration = float(sample_info[4])
    # TODO problem with titrant and id assignement here, plus concnetration is just approx and will be overwritten by calibration
    HCl_aliquot.concentration = float(sample_info[5])

    datatypes = {
        "time": str,
        "emf": float,
        "t_sample": float,
        "weight": float,
        "pH_est": float,
        "volume": float,
        "t_HCl": float,
        "t_NaOH": float,
        "t_air": float,
    }

    # process and pivot
    if fwd_data:
        # make dict of data
        fwd_data_processed = dict(zip(datatypes.keys(), map(list, zip(*fwd_data))))
        # type cast
        fwd_typecast = {
            key: [datatypes[key](value) for value in values]
            for key, values in fwd_data_processed.items()
        }
        HCl_aliquot.weight = fwd_typecast["weight"][0] / 1000

    else:
        raise DataMissing(f"There is no HCl data in {filename}, unable to proceed.")
        # TODO maybe if there are enough files ahead of it, use all of them and don't count the last one, but I can't keep going bc
        # need this data for the rest to be valid

    if BWD_data:
        BWD_data_processed = dict(zip(datatypes.keys(), map(list, zip(*BWD_data))))
        # type cast
        BWD_typecast = {
            key: [datatypes[key](value) for value in values]
            for key, values in BWD_data_processed.items()
        }
        # take burette name, read in burette_density, grab formula
        volume_corrected = correct_burette_volume(burette_id, BWD_typecast["volume"])

        titration_weights = v_to_w(
            volume_corrected, BWD_typecast["t_NaOH"], "NaOH", titrant.ID
        )
    else:
        raise TitrantDataMissing(
            f"There is no NaOH data in {filename}, unable to proceed."
        )

    titration_emf = BWD_typecast["emf"]
    titration_temp = BWD_typecast["t_sample"]
    titration_data = Titration(titration_weights, titration_emf, titration_temp)

    return titrant, sample, HCl_aliquot, titration_data


def NaOH_calibration_data(
    filename: str,
    burette_id: str = "dosimat 12",
    sample: Solution = None,
    titrant: Solution = None,
) -> tuple[Solution, Solution, Solution, Titration]:
    HCl_aliquot = Solution()

    FWD_data = list()
    # # Open filename and extract data
    with open(filename, "r") as datafile:
        csvreader = csv.reader(datafile)
        sample_info = next(csvreader)
        # if first time initializing sample
        if not sample:
            # assumes calibration solution type found in file name
            if "nacl" in filename.lower():
                sample = NaCl()
            elif "kcl" in filename.lower():
                sample = KCl()

            sample.w0 = float(sample_info[0]) / 1000
            sample.salt = float(sample_info[1])
            sample.emf0 = float(sample_info[2])

        if not titrant:
            titrant = Solution()
            if sample_info[6]:
                titrant.ID = sample_info[6].split("-")[0]
                # TODO missing batch/bag number
            else:
                titrant.ID = "nan"
            # TODO how does it get the HCl concentration
            # HCl_batch = sample_info[7].split("-")[0]

        # check if FWD
        for row in csvreader:
            if "BWD" in row:
                BWD_data = list()
                break
            else:
                FWD_data.append(row)
        if "BWD_data" in locals():
            for row in csvreader:
                BWD_data.append(row)

    t0 = float(sample_info[3])
    # flag sample as Q (questionable) if temperature very out of range
    if t0 < 15 or t0 > 30:
        sample.flag = "Q"

    titrant.concentration = float(sample_info[4])
    # TODO problem with titrant and id assignement here, plus concnetration is just approx and will be overwritten by calibration
    HCl_aliquot.concentration = float(sample_info[5])

    datatypes = {
        "time": str,
        "emf": float,
        "t_sample": float,
        "weight": float,
        "pH_est": float,
        "volume": float,
        "t_HCl": float,
        "t_NaOH": float,
        "t_air": float,
    }

    # process and pivot
    if FWD_data:
        # make dict of data
        FWD_data_processed = dict(zip(datatypes.keys(), map(list, zip(*FWD_data))))
        # type cast
        FWD_typecast = {
            key: [datatypes[key](value) for value in values]
            for key, values in FWD_data_processed.items()
        }
        HCl_aliquot.weight = FWD_typecast["weight"][0] / 1000

    else:
        raise DataMissing(f"There is no HCl data in {filename}, unable to proceed.")
        # TODO maybe if there are enough files ahead of it, use all of them and don't count the last one, but I can't keep going bc
        # need this data for the rest to be valid

    if BWD_data:
        BWD_data_processed = dict(zip(datatypes.keys(), map(list, zip(*BWD_data))))
        # type cast
        BWD_typecast = {
            key: [datatypes[key](value) for value in values]
            for key, values in BWD_data_processed.items()
        }
        # take burette name, read in burette_density, grab formula
        volume_corrected = correct_burette_volume(burette_id, BWD_typecast["volume"])

        titration_weights = v_to_w(
            volume_corrected, BWD_typecast["t_NaOH"], "NaOH", titrant.ID
        )
    else:
        raise TitrantDataMissing(
            f"There is no NaOH data in {filename}, unable to proceed."
        )

    titration_emf = BWD_typecast["emf"]
    titration_temp = BWD_typecast["t_sample"]
    titration_data = Titration(titration_weights, titration_emf, titration_temp)

    return titrant, sample, HCl_aliquot, titration_data


def correct_burette_volume(
    burette_id: str, volume: Union[str, list]
) -> Union[str, list]:

    x0, x1, x2, x3, x4, x5 = get_density_coefficients("burette_density", burette_id)

    return [
        x0 + x1 * vol + x2 * (vol**2) + x3 * (vol**3) + x4 * (vol**4) + x5 * (vol**5)
        for vol in volume
    ]


def v_to_w(
    volume: list[float], temperature: list[float], solution_type: str, solution_id: str
) -> list[float]:
    # TODO if solution_id is missing, use a generic formula, add flag to all results so poorer flag if not use specific coefficients

    if "-" in solution_id:
        solution_id = solution_id.split("-")[0]

    x0, x1, x2, x3, x4, x5 = get_density_coefficients(solution_type, solution_id)
    NaOH_density = [
        x0
        + x1 * temp
        + x2 * (temp**2)
        + x3 * (temp**3)
        + x4 * (temp**4)
        + x5 * (temp**5)
        for temp in temperature
    ]
    # returns the air buoyancy corrected value
    return [rho * v / 1000 for rho, v in zip(NaOH_density, volume)]


def get_density_coefficients(
    solution_type: str, id: str
) -> tuple[float, float, float, float, float, float]:
    """
    Method to get quadratic formula coefficients from csv files

    Args:
        solution_type (str):
        id (str): identifier

    Returns:
        float: coefficients appropriate for a x0 + x1 * y + x2 * (y**2) + x3 * (y**3) + x4 * (y**4) + x5 * (y**5) equation

    """
    batch_data = get_matching_files("auxiliary_data/", solution_type, "csv")[0]

    if id is None or id == "nan":
        raise CalibrationDataMissing("Solution identifier is invalid")

    df = pd.read_csv(batch_data)
    coefficients = df[df["id"] == id]

    x0 = float(coefficients["x0"].values[0])
    x1 = float(coefficients["x1"].values[0])
    x2 = float(coefficients["x2"].values[0])
    x3 = float(coefficients["x3"].values[0])
    x4 = float(coefficients["x4"].values[0])
    x5 = float(coefficients["x5"].values[0])

    return x0, x1, x2, x3, x4, x5


def get_concentration_ionicstrength(keyword: str, id: str) -> tuple[float, float]:
    """
    Method to get intrinsic solution properties from a batch data file for give solution_type

    Args:
        solution_type (str):
        id (str): identifier

    Returns:
        tuple[floats]:
    """

    if id is None or id == "nan":
        raise CalibrationDataMissing("Solution identifier is invalid")

    batch_data = get_matching_files("auxiliary_data/", keyword, "csv")[0]

    df = pd.read_csv(batch_data)
    coefficients = df[df["id"] == id]

    concentration = float(coefficients["c"].values[0])
    ionic_strength = float(coefficients["I"].values[0])

    return concentration, ionic_strength
