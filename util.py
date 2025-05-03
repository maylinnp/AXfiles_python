from pathlib import Path


def get_matching_files(folder: str, pattern: str, extension: str) -> list:
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
