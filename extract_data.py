## Extracts data from AX titrations
# Data is saved in one structure per file
# Created 2020/09/08 MLP based on old extract-data codes
# Updated 2020/11/22 MLP: I have edited the Labview code, so that it
# applies volume correction per every 5 mL.  I also edited it so that on
# each line, the volume listed is what was added, waited 15 s, then pH on
# the corresponding line recorded (i.e., volume and pH value correspond to
# eachother).
import csv
import re
import os
from exceptions import *
from typing import Union
import pandas as pd
import logging
from solutions import *
import statistics
from util import get_matching_files

logger = logging.getLogger(__name__)

logger.setLevel("DEBUG")

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


# TODO make a class that has all these methods in the generic form,
# Then child classes that have extra specific methods
# and use other custom classes to hold the final data so that they don't have anyhting else
def titration_data(
    filename: str, burette_id: str = "dosimat 12"
) -> tuple[Solution, Titration, Titration]:
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
        if sample.type.lower() == "sw":
            sample.S = float(sample_info[1])
        else:
            sample.I = float(sample_info[1])
        sample.emf0 = float(sample_info[2])
        t0 = float(sample_info[3])
        HCl_id = sample_info[7].split("-")[0]
        NaOH_id = sample_info[6].split("-")[0]
        # check for nutrient data
        match = re.match(r"\d{8} ([^-]+)-[A-Za-z]\.csv", os.path.basename(filename))
        if match:
            sample.id = match.group(1)
        else:
            sample.id = "any"

        # check if FWD
        fwd_data = list()
        for row in csvreader:
            # check if bwd
            if "BWD" in row:
                bwd_data = list()
                break
            else:
                fwd_data.append(row)
        if fwd_data:
            # get titrant data
            HCl_conc, HCl_I = get_concentration_ionicstrength("HCl", HCl_id)
            HCl_titrant = Titrant("HCl", HCl_id, HCl_conc, HCl_I)
            # maps the column data, assuming the order of types in datatypes
            fwd_data_processed = dict(zip(datatypes.keys(), map(list, zip(*fwd_data))))
            # type cast
            fwd_typecast = {
                key: [datatypes[key](value) for value in values]
                for key, values in fwd_data_processed.items()
            }
            HCl_vol_burette_corrected = correct_burette_volume(
                burette_id, fwd_typecast["volume"]
            )
            HCl_weights = v_to_w(
                HCl_vol_burette_corrected,
                fwd_typecast["t_HCl"],
                "HCl",
                HCl_id,
            )
            titration_emf = fwd_typecast["emf"]
            titration_temp = fwd_typecast["t_sample"]
            HCl_titration_data = Titration(
                HCl_weights, titration_emf, titration_temp, HCl_titrant
            )
        else:
            HCl_titration_data = None

        if "bwd_data" in locals():
            for row in csvreader:
                bwd_data.append(row)
        if bwd_data:
            NaOH_conc, NaOH_I = get_concentration_ionicstrength("NaOH", NaOH_id)
            NaOH_titrant = Titrant("NaOH", NaOH_id, NaOH_conc, NaOH_I)
            # maps the column data, assuming the order of types in datatypes
            bwd_data_processed = dict(zip(datatypes.keys(), map(list, zip(*bwd_data))))
            # type cast
            bwd_typecast = {
                key: [datatypes[key](value) for value in values]
                for key, values in bwd_data_processed.items()
            }
            NaOH_vol_burette_corrected = correct_burette_volume(
                burette_id, bwd_typecast["volume"]
            )
            NaOH_weights = v_to_w(
                NaOH_vol_burette_corrected,
                bwd_typecast["t_NaOH"],
                "NaOH",
                NaOH_id,
            )
            titration_emf = bwd_typecast["emf"]
            titration_temp = bwd_typecast["t_sample"]
            NaOH_titration_data = Titration(
                NaOH_weights, titration_emf, titration_temp, NaOH_titrant
            )
        else:
            NaOH_titration_data = None

    # flag sample as Q (questionable) if temperature very out of range
    if t0 < 15 or t0 > 30:
        sample.flag = "Q"

    return sample, HCl_titration_data, NaOH_titration_data


def NaOH_calibration_data(
    filename: str,
    burette_id: str = "dosimat 12",
    sample: Solution = None,
    titrant: Solution = None,
) -> tuple[Solution, Solution, Solution, Titration]:
    HCl_aliquot = Solution()

    fwd_data = list()
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
            sample.salt_value = float(sample_info[1])
            sample.salt_type = "ionic strength"
            sample.emf0 = float(sample_info[2])

        if not titrant:
            titrant = Solution()
            if sample_info[6]:
                titrant.id = sample_info[6].split("-")[0]
            else:
                titrant.id = "nan"

        # check if FWD
        for row in csvreader:
            if "BWD" in row:
                bwd_data = list()
                break
            else:
                fwd_data.append(row)
        if "bwd_data" in locals():
            for row in csvreader:
                bwd_data.append(row)

    t0 = float(sample_info[3])
    # flag sample as Q (questionable) if temperature very out of range
    if t0 < 15 or t0 > 30:
        sample.flag = "Q"

    HCl_aliquot.conc = float(sample_info[5])

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

    if bwd_data:
        bwd_data_processed = dict(zip(datatypes.keys(), map(list, zip(*bwd_data))))
        # type cast
        bwd_typecast = {
            key: [datatypes[key](value) for value in values]
            for key, values in bwd_data_processed.items()
        }
        # take burette name, read in burette_density, grab formula
        volume_corrected = correct_burette_volume(burette_id, bwd_typecast["volume"])

        titration_weights = v_to_w(
            volume_corrected, bwd_typecast["t_NaOH"], "NaOH", titrant.id
        )
    else:
        raise TitrantDataMissing(
            f"There is no NaOH data in {filename}, unable to proceed."
        )

    titration_emf = bwd_typecast["emf"]
    titration_temp = fwd_typecast["t_sample"]
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
