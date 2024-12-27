# class to hold titrant data


class titrant:

    def __init__(self, name: str, concentration: float):
        self.name = name
        self.c = concentration
        self._density = 1
        self.weight = list()

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self):
        # this should read from a datasheet the first time, and save the density equation in memory
        if self._density == 1:
            # check if density euqation in datasheet
            # return density equation
            pass

        a = 0
        b = 0
        c = 1
        d = 0

        return lambda weight: weight * a**3 + weight * b**2 + weight * c + d
