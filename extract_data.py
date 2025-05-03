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

logger = logging.getLogger(__name__)

logger.setLevel("DEBUG")


# class Solution(object):
#     pass


# TODO make a class that has all these methods in the generic form,
# Then child classes that have extra specific methods
# and use other custom classes to hold the final data so that they don't have anyhting else


def NaOH_calibration_data(
    filename: str, burette: str = "dosimat 12"
) -> tuple[Solution, Solution, Solution, Titration]:
    titrant = Solution()
    HCl_aliquot = Solution()

    FWD_data = list()
    # # Open filename and extract data
    with open(filename, "r") as datafile:
        # assumes calibration solution type found in file name
        if "nacl" in filename.lower():
            sample = NaCl()
        elif "kcl" in filename.lower():
            sample = KCl()
        csvreader = csv.reader(datafile)
        sample_info = next(csvreader)
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

    sample.w0 = float(sample_info[0]) / 1000
    sample.I = float(sample_info[1])
    sample.emf0 = float(sample_info[2])
    t0 = float(sample_info[3])
    # flag sample as Q (questionable) if temperature very out of range
    if t0 < 15 or t0 > 30:
        sample.flag = "Q"

    titrant.concentration = float(sample_info[4])
    HCl_aliquot.concentration = float(sample_info[5])
    if sample_info[6]:
        titrant.ID = sample_info[6]
        NaOH_batch = titrant.ID.split("-")[0]
    else:
        titrant.ID = "nan"
    # HCl_batch = sample_info[7].split("-")[0]

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
        volume_corrected = correct_burette_volume(burette, BWD_typecast["volume"])

        titration_weights = v_to_w(volume_corrected, BWD_typecast["t_NaOH"], NaOH_batch)
    else:
        raise TitrantDataMissing(
            f"There is no NaOH data in {filename}, unable to proceed."
        )

    titration_emf = BWD_typecast["emf"]
    titration_temp = BWD_typecast["t_sample"]
    titration_data = Titration(titration_weights, titration_emf, titration_temp)

    return titrant, sample, HCl_aliquot, titration_data


def correct_burette_volume(burette: str, volume: Union[str, list]) -> Union[str, list]:

    x0, x1, x2, x3, x4, x5 = get_coefficients(
        "calibration_data/burette_density.csv", "burette_id", burette
    )

    return [
        x0 + x1 * vol + x2 * (vol**2) + x3 * (vol**3) + x4 * (vol**4) + x5 * (vol**5)
        for vol in volume
    ]


def v_to_w(volume: list[float], temperature: list[float], solution_id) -> list[float]:
    # TODO if solution_id is missing, use a generic formula, add flag to all results so poorer flag if not use specific coefficients
    # corrects mL to kg
    x0, x1, x2, x3, x4, x5 = get_coefficients(
        "calibration_data/NaOH_density.csv", "NaOH_id", solution_id
    )
    NaOH_density = [
        x0
        + x1 * temp
        + x2 * (temp**2)
        + x3 * (temp**3)
        + x4 * (temp**4)
        + x5 * (temp**5)
        for temp in temperature
    ]

    return [rho * v / 1000 for rho, v in zip(NaOH_density, volume)]


def get_coefficients(file: str, keyword: str, id: str):
    """
    Method to get quadratic formula coefficients from csv files

    Args:
        file (str): comma separated file that holds the coefficients
        keyword (str): name of column that describes the subject, for example "NaOH_ID" or "burette_id"
        id (str): identifier

    Returns:
        float: coefficients appropriate for a x0 + x1 * y + x2 * (y**2) + x3 * (y**3) + x4 * (y**4) + x5 * (y**5) equation
    """

    if id is None or id == "nan":
        raise CalibrationDataMissing("Solution identifier is invalid")

    df = pd.read_csv(file)
    coefficients = df[df[keyword] == id]

    x0 = float(coefficients["x0"].values[0])
    x1 = float(coefficients["x1"].values[0])
    x2 = float(coefficients["x2"].values[0])
    x3 = float(coefficients["x3"].values[0])
    x4 = float(coefficients["x4"].values[0])
    x5 = float(coefficients["x5"].values[0])

    return x0, x1, x2, x3, x4, x5
