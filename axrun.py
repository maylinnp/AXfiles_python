# script for testing AX code

import extract_data

titrant, sample, HCl_aliquot = extract_data.NaOH_calibration_data(
    "/Users/maylinnp/personal_code/AX_R/NaOHdata/20210607 NaCl-38-A.csv"
)

print("titrant", titrant.__dict__)
print("\nsample", sample.__dict__)
print("\nHCl_aliquot", HCl_aliquot.__dict__)
