## Extracts data from AX titrations
# Data is saved in one structure per file
# Created 2020/09/08 MLP based on old extract-data codes
# Updated 2020/11/22 MLP: I have edited the Labview code, so that it
# applies volume correction per every 5 mL.  I also edited it so that on
# each line, the volume listed is what was added, waited 15 s, then pH on
# the corresponding line recorded (i.e., volume and pH value correspond to
# eachother).
import csv
from itertools import zip_longest


def NaOH_calibration_data(filename: str, burette: str = "dosimat 12") -> dict:
    FWD_data = list()
    # # Open filename and extract data
    with open(filename, "r") as datafile:
        csvreader = csv.reader(datafile)
        sample_info = next(csvreader)
        print(sample_info)
        # check if FWD
        for row in csvreader:
            if "BWD" in row:
                print("There is back titration data")
                BWD_data = list()
                break
            else:
                FWD_data.append(row)
        if "BWD_data" in locals():
            for row in csvreader:
                BWD_data.append(row)

    w0 = sample_info[0]
    I = sample_info[1]
    emf0 = sample_info[2]
    t0 = sample_info[3]
    CNaOH = sample_info[4]
    CHCl = sample_info[5]
    if sample_info[6]:
        NaOH_ID = sample_info[6]
        NaOH_batch = NaOH_ID.split("-")[0]
    else:
        NaOH_ID = None
    HCl_ID = sample_info[7]

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
        print(FWD_data_processed)

        # TODO this will always just be weight of acid, so probably doesn't need this complexity, but save and use for later
    if BWD_data:
        BWD_data_processed = dict(zip(datatypes.keys(), map(list, zip(*BWD_data))))
        # type cast
        BWD_typecast = {
            key: [datatypes[key](value) for value in values]
            for key, values in BWD_data_processed.items()
        }
        # take burette name, read in burette_density, grab formula
        volume_corrected = correct_burette_volume(burette, BWD_typecast["volume"])

        NaOH_weight = v_to_w(volume_corrected, BWD_typecast["t_NaOH"], NaOH_batch)


def correct_burette_volume(burette: str, volume: str | list) -> str | list:

    x0, x1, x2, x3, x4, x5 = get_coefficients(
        "calibration_data/burette_density.csv", "burette_id", burette
    )

    return [
        x0 + x1 * vol + x2 * (vol**2) + x3 * (vol**3) + x4 * (vol**4) + x5 * (vol**5)
        for vol in volume
    ]


def v_to_w(volume: list[float], temperature: list[float], solution_id) -> list[float]:
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
    import pandas as pd

    df = pd.read_csv(file)
    coefficients = df[df[keyword] == id]

    x0 = float(coefficients["x0"].values[0])
    x1 = float(coefficients["x1"].values[0])
    x2 = float(coefficients["x2"].values[0])
    x3 = float(coefficients["x3"].values[0])
    x4 = float(coefficients["x4"].values[0])
    x5 = float(coefficients["x5"].values[0])

    return x0, x1, x2, x3, x4, x5


# Now I have all data, then I need to convert volume and temp into weight titrants


#     # # See if there is BWD titration data present
#     Index = find(contains(Lines,'BWD'))
#     if isempty(Index)
#         Index = length(Lines)
#         A = 0
#     else
#         A = 1
#     end
#     # # Extract FWD titration data
#     clear i
#         R         = strsplit(Lines{2,1},',')
#         wHCl = str2double(R{6})/1000

#     clear R
#     # # BWD data, if it exists
#     if A == 0
#         emf = nan
#         t   = nan
#         w   = nan
#     else
#         for i = Index+1:length(Lines)
#         R         = strsplit(Lines{i,1},',')
#         emf(i-Index) = str2double(R{2})
#         t(i-Index)   = str2double(R{3})
#         V2uncorr(i-Index)   = str2double(R{6})
#         Bur2T(i-Index)= str2double(R{8})
#         end
#     end

#     end

# return [w0,I,emf0,t0,CNaOH,CHCl,wHCl,emf,t,w,NaOH_ID]
