from pathlib import Path
import pandas as pd


def get_matching_files(folder: str, pattern: str, extension: str) -> list[str]:
    """
    Method that gets all files in a folder with a given unifying pattern and
    specific extension. Will return a list sorted by filename

    Args:
        folder (str): folder in which to look for files
        pattern (str): identifying file pattern
        extension (str): file extension

    Returns:
        list: ordered by filename
    """
    folder_path = Path(folder)
    return sorted(
        [
            str(file)
            for file in folder_path.glob(f"*{pattern}*{extension}")
            if file.is_file()
        ]
    )


def get_system_constant(constant: str) -> float:
    """
    Method to get system constant

    Args:
        constant (str): what constant to get system value for

    Returns:
        float: value of constant

    """
    df = pd.read_csv("auxiliary_data/system_constants.csv")
    row = df[df["constant"] == constant]

    value = float(row["value"].values[0])

    return value
