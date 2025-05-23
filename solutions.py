# base class for any solution
from statistics import mean
from math import exp, log, log10
from typing import Union
import numpy as np
from util import *
from ax_maths import k_boltz
import csv
import os

try:
    import gsw
except:

    def SW_rho():
        return 1.026

    print(
        "GSW is not installed, will use a simple formula for seawater density if necessary"
    )


# TODO the point of the properties here is that some vlaues can be (re)calculated on the fly when the nsolution changes,
# including temperature or adding stuff so that volume and concentartion and ionic strength/salinity changes
# TODO add all constants to all relevant solutions
# TODO check if I am able to calc stuff on the fly by supplying a new temperature

# TODO implement all constants for all solutions (or void options)
salinity_aliases = ["salinity", "s", "sal"]
ionic_strength_aliases = ["ionic strength", "i", "ionic", "ionic_strength"]

nutrient_path = "auxiliary_data/nutrients.csv"


class Solution:

    # everything is estimated at 20 deg. C
    def __init__(
        self,
        type: str = "solution",
        salt_value: float = 0,
        salt_type: str = "salinity",
        temp: float = 20,
        id: str = "any",
    ):
        self.type = type
        self.id = id
        # parse salt
        if salt_type.lower() in salinity_aliases:
            self.S = salt_value
        elif salt_type.lower() in ionic_strength_aliases:
            self.I = salt_value
        # parse temp
        if temp > 200:
            self.T = temp
            self.t = self.T - 273.15
        else:
            self.t = temp
            self.T = self.t + 273.15

        self.flag = "A"
        self.conc = None
        self.weight = None
        self.w0 = None
        self.emf0 = None
        self.CT_degas = 2.5e-6
        self.SiT = 1e-6
        self.PT = 0e-6

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = value
        # check for nutrients
        self._look_for_nutrients()

    @property
    def S(self):
        return self._S

    @S.setter
    def S(self, value):
        self._S = value
        self._I = 19.924 * value / (1000 - 1.005 * value)

    @property
    def I(self):
        return self._I

    @I.setter
    def I(self, value):
        self._I = value
        self._S = 1000 * value / (1.005 * value + 19.924)

    @property
    def m0(self):
        return self.w0 * ((1 - (0.0012013 / 8)) / (1 - (0.0012013 / self.rho)))

    @property
    def k(self):
        return 8.31451 * self.T / 96484.56

    @property
    def KW(self):
        return 10 ^ -14

    # K1 and K2 is assumed very similar in NaCl and KCl, only property of ionic strength
    @property
    def K1(self):
        # i is the ionic strength that is presumed to be a property of the solution using this version of K1
        # by necessity on H+(free)
        m = self.I * 1000 / (1000 - self.I * (22.99 + 35.45))
        A = 35.2911 * m**0.5 + 0.8491 * m - 0.32 * m**1.5 + 0.055 * m**2
        B = -1583.09 * m**0.5
        C = -5.4366 * m**0.5
        pK1 = (
            -402.56788
            + 11656.46 / self.T
            + 72.173 * log(self.T)
            - 0.161325 * self.T
            + 7.5526e-5 * self.T**2
        )
        return 10 ** -(A + B / self.T + C * log(self.T) + pK1)

    @property
    def K2(self):
        # i is the ionic strength that is presumed to be a property of the solution using this version of K1
        # by necessity on H+(free)
        m = self.I * 1000 / (1000 - self.I * (22.99 + 35.45))
        A = 38.2746 * m**0.5 + 1.6057 * m - 0.647 * m**1.5 + 0.113 * m**2
        B = -1738.16 * m**0.5
        C = -6.0346 * m**0.5
        pK2 = -122.4994 + 5811.18 / self.T + 20.5263 * log(self.T) - 0.0120897 * self.T
        return 10 ** -(A + B / self.T + C * log(self.T) + pK2)

    def nutrient(self):
        # TODO maybe a way to get nutrient concentration is if solution is identified
        # by some name, and if that name is found in nutrients then give value, else
        # nutrients are 0 (i.e., "any")
        pass

    @property
    def KS(self):
        pass

    @property
    def KF(self):
        pass

    @property
    def KB(self):
        pass

    @property
    def KSi(self):
        pass

    @property
    def KP1(self):
        pass

    @property
    def KP2(self):
        pass

    @property
    def KP3(self):
        pass

    @property
    def KNH4(self):
        pass

    @property
    def KNO2(self):
        pass

    def _look_for_nutrients(self):
        if not os.path.exists(nutrient_path):
            return
        with open(nutrient_path, newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            next(reader)
            for row in reader:
                if row and row[reader.fieldnames[0]] == self.id:
                    self.SiT = float(row["silicate"]) * 1e-6
                    self.PT = float(row["phosphate"]) * 1e-6
                    break


class NaCl(Solution):

    def __init__(self, concentration: float = None):
        super().__init__("NaCl")
        self.c = concentration
        self.I = self.c  # calculate from concentration

        # assume these have not been added
        self.ST = 0
        self.FT = 0
        self.BT = 0
        self.rho = 1

    @property
    def KW(self):
        p00 = 14.83
        p10 = -0.4914
        p01 = -0.0471
        p20 = 0.3917
        p11 = 0.0001724
        p02 = 0.0004381
        p30 = -0.07454
        p21 = 3.419e-06
        p12 = -1.081e-05
        p03 = -4.382e-06
        return 10 ** -(
            p00
            + p10 * self.c
            + p01 * self.t
            + p20 * self.c**2
            + p11 * self.c * self.t
            + p02 * self.t**2
            + p30 * self.c**3
            + p21 * self.c**2 * self.t
            + p12 * self.c * self.t**2
            + p03 * self.t**3
        )


class KCl(Solution):

    def __init__(self, concentration: float = None):
        super().__init__("KCl")
        self.c = concentration
        self.I = self.c  # calculate from concentration

        # assume these have not been added
        self.ST = 0
        self.FT = 0
        self.BT = 0
        self.rho = 1

    @property
    def KW(self):
        p00 = 14.86
        p10 = -1.062
        p01 = -0.0479
        p20 = 2.031
        p11 = 0.0005349
        p02 = 0.000483
        p30 = -1.095
        p21 = -0.0005724
        p12 = -9.838e-06
        p03 = -5.091e-06
        return 10 ** -(
            p00
            + p10 * self.c
            + p01 * self.t
            + p20 * self.c**2
            + p11 * self.c * self.t
            + p02 * self.t**2
            + p30 * self.c**3
            + p21 * self.c**2 * self.t
            + p12 * self.c * self.t**2
            + p03 * self.t**3
        )


class SW(Solution):

    def __init__(
        self,
    ):
        super().__init__("SW")
        self.CT_degas = 4e-6
        self.SiT = 5e-6
        self.PT = 0.5e-6

    @property
    def KW(self):
        # used with H+(tot) --> not the original equation
        return exp(
            -13847.26 / self.T
            + 148.9652
            - 23.652 * log(self.T)
            + (118.67 / self.T - 5.977 + 1.0495 * log(self.T)) * self.S**0.5
            - 0.01615 * self.S
        )

    @property
    def K1(self):
        # used with H+(tot)
        return 10 ** (
            -3633.86 / self.T
            + 61.2172
            - 9.67770 * log(self.T)
            + 0.011555 * self.S
            - 0.0001152 * self.S**2
        )

    @property
    def K2(self):
        # used with H+(tot)
        return 10 ** (
            -471.78 / self.T
            - 25.9290
            + 3.16967 * log(self.T)
            + 0.01781 * self.S
            - 0.0001122 * self.S**2
        )

    @property
    def ST(self):
        return 0.14 / 96.062 * self.S / 1.80655

    @property
    def FT(self):
        return 0.000067 / 18.998 * self.S / 1.80655

    @property
    def BT(self, BTrat: str = "uppstrom", ratio: float = 0.5):
        # ratio will favor Lee, i.e., higher ratio will weight Lee BT/S heavier
        if ratio > 1:
            raise Exception("Ratio for BT/S cannot be above 1.")
        if BTrat.lower() == "uppstrom":
            return self._BT_uppstrom()
        elif BTrat.lower() == "lee":
            return self._BT_lee()
        elif BTrat.lower() == "mix":
            return self._BT_lee() * ratio + self._BT_uppstrom() * (1 - ratio)

    def _BT_uppstrom(self):
        return 0.000232 / 10.811 * self.S / 1.80655

    def _BT_lee(self):
        return 0.0002414 / 10.811 * self.S / 1.80655

    @property
    def KS(self):
        # used with H+(free)
        return exp(
            -4276.1 / self.T
            + 141.328
            - 23.093 * log(self.T)
            + (-13856 / self.T + 324.57 - 47.986 * log(self.T)) * self.I**0.5
            + (35474 / self.T - 771.54 + 114.723 * log(self.T)) * self.I
            - 2698.0 / self.T * self.I**1.5
            + 1776.0 / self.T * self.I**2
            + log(1 - 0.001005 * self.S)
        )

    @property
    def KF(self):
        # used with H+(tot)
        return exp(874.0 / self.T - 9.68 + 0.111 * self.S**0.5)

    @property
    def KB(self):
        # used with H+(tot)
        return exp(
            (
                -8966.90
                - 2890.53 * self.S**0.5
                - 77.942 * self.S
                + 1.728 * self.S**1.5
                - 0.0996 * self.S**2
            )
            / self.T
            + 148.0248
            + 137.1942 * self.S**0.5
            + 1.62142 * self.S
            + (-24.4344 - 25.085 * self.S**0.5 - 0.2474 * self.S) * log(self.T)
            + 0.053105 * self.S**0.5 * self.T
        )

    @property
    def KSi(self):
        # from Dickson et al. 2007
        return exp(
            -8904.2 / self.T
            + 117.385
            - 19.334 * log(self.T)
            + (-458.79 / self.T + 3.5913) * self.I**0.5
            + (188.74 / self.T - 1.5998) * self.I
            + (-12.1652 / self.T + 0.07871) * self.I**2
            + log(1 - 0.001005 * self.S)
        )

    @property
    def KP1(self):
        # from Dickson et al. 2007
        return exp(
            -4576.752 / self.T
            + 115.525
            - 18.453 * log(self.T)
            + (-106.736 / self.T + 0.69171) * self.S**0.5
            + (-0.65643 / self.T - 0.01844) * self.S
        )

    @property
    def KP2(self):
        # from Dickson et al. 2007
        return exp(
            -8814.715 / self.T
            + 172.0883
            - 27.927 * log(self.T)
            + (-160.340 / self.T + 1.3566) * self.S**0.5
            + (0.37335 / self.T - 0.05778) * self.S
        )

    @property
    def KP3(self):
        # from Dickson et al. 2007
        return exp(
            -3070.75 / self.T
            - 18.141
            + (17.27039 / self.T + 2.81197) * self.S**0.5
            + (-44.99486 / self.T - 0.09984) * self.S
        )

    @property
    def rho(self):
        """In g/mL"""
        absolute_S = self.S  # gsw.SA_from_SP(self.S, 10, 32, -117)
        return gsw.density.rho(absolute_S, self.t, 0) / 1000


class Titrant(Solution):

    def __init__(
        self,
        type: str,
        id: str,
        concentration: float,
        ionic_strength: float,
    ):
        self.type = type
        self.id = id
        self.name = self.type + "-" + self.id
        self.concentration = concentration
        self.ionic_strength = ionic_strength

    @property
    def density(self):
        # this should read from a datasheet the first time, and save the density equation in memory

        a = 0
        b = 0
        c = 1
        d = 0

        return lambda weight: weight * a**3 + weight * b**2 + weight * c + d


class Titration:
    def __init__(
        self,
        weight: list[float],
        emf: list[float],
        temp: Union[float, list[float]],
        titrant: Titrant = None,
    ):
        self.weight = np.array(weight, dtype=np.float64)
        # TODO option to give back mass/air buoyancy corrected, will depend on titrant characteristics
        self.emf = np.array(emf, dtype=np.float64)
        if mean(temp) < 100:
            self.t = np.array(temp, dtype=np.float64)
            self.T = self.t + 273.15
        self.titrant = titrant
        # estimate pH from system constnats
        self.E0 = get_system_constant("E0")
        self.k = k_boltz(np.mean(self.T))
        self.pH_est = -np.log10(np.exp((self.emf - self.E0) / self.k))

    def recalculate_pH(self, new_E0):

        new_pH = self.pH_est + (self.E0 - new_E0) / self.k * log(10)
        self.E0 = new_E0
        # TODO probably need some guard here against bad E0 values
        self.pH_est = new_pH
