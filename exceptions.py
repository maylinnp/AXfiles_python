class DataMissing(Exception):
    pass


class FileMissing(Exception):
    pass


class CalibrationDataMissing(DataMissing):
    pass


class TitrantDataMissing(DataMissing):
    pass


class TitrationDataMissing(DataMissing):
    pass
