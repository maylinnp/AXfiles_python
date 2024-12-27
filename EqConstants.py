## SWconst calculates all acid-base constants for a given salinity and
# temperature, as well as total concentrations of conservative species.
# 2017/8/8 May-Linn Paulsen
# 2020/09/10 edited from SWconst(both regular and _NaCl) to be a function
# to be used instead, where some input determines NaCl, KCl, or SW. -MLP
# 2021/02/01 Updated to reflect what constants are on free and total H+
# scale, and able to use input as ionic strength (< 1 unit input) or
# salinity (> 1 unit input). It will take temperature and ensure it is on
# the relevant scale, and has a BT/S ratio input as wel. It includes
# functions for KW in either NaCl or KCl, and K1&K2 in chloride solution
# instead of real seawater.
from statistics import mean
from math import exp, log


def eqConstants(S_or_I, t, Sol, BTrat) -> list:
    # # for Sol, 1 = NaCl solution, 2 = real SW, 3 = KCl solution, and BTrat, 1 = Uppstrom, 2 = Liu et al.
    ## Thermodynamic
    if t < 273.15:  # check if temp in Kelvin or Celsius
        t = mean(t)
        T = t + 273.15
    else:
        T = mean(t)
        t = T - 273.15

    if S_or_I > 1:
        I = 19.924 * S_or_I / (1000 - 1.005 * S_or_I)
        S = S_or_I
    else:
        S = 1000 * S_or_I / (1.005 * S_or_I + 19.924)
        I = S_or_I

    k = 8.31451 * mean(T) / 96484.56
    ## Conservative concentrations
    if Sol == 1 or Sol == 3:
        ST = 0
        FT = 0
        BT = 0
    else:
        ST = 0.14 / 96.062 * S / 1.80655
        FT = 0.000067 / 18.998 * S / 1.80655
        if BTrat == 1:
            BT = 0.000232 / 10.811 * S / 1.80655
        elif BTrat == 2:
            BT = 0.0002414 / 10.811 * S / 1.80655

    ## Acid dissociation constants
    # Conservative species
    KB = exp(
        (-8966.90 - 2890.53 * S**0.5 - 77.942 * S + 1.728 * S**1.5 - 0.0996 * S**2) / T
        + 148.0248
        + 137.1942 * S**0.5
        + 1.62142 * S
        + (-24.4344 - 25.085 * S**0.5 - 0.2474 * S) * log(T)
        + 0.053105 * S**0.5 * T
    )  # used with H+(tot)

    KS = exp(
        -4276.1 / T
        + 141.328
        - 23.093 * log(T)
        + (-13856 / T + 324.57 - 47.986 * log(T)) * I**0.5
        + (35474 / T - 771.54 + 114.723 * log(T)) * I
        - 2698.0 / T * I**1.5
        + 1776.0 / T * I**2
        + log(1 - 0.001005 * S)
    )  # used with H+(free)

    KF = exp(874.0 / T - 9.68 + 0.111 * S**0.5)  # used with H+(tot)

    # Nutrients, I do believe all of these are on the H+(tot) scale, at least I
    # can't find a clear indication in Dickson et al. 2007 that they aren't.
    KSi = exp(
        -8904.2 / T
        + 117.385
        - 19.334 * log(T)
        + (-458.79 / T + 3.5913) * I**0.5
        + (188.74 / T - 1.5998) * I
        + (-12.1652 / T + 0.07871) * I**2
        + log(1 - 0.001005 * S)
    )

    KP1 = exp(
        -4576.752 / T
        + 115.525
        - 18.453 * log(T)
        + (-106.736 / T + 0.69171) * S**0.5
        + (-0.65643 / T - 0.01844) * S
    )
    KP2 = exp(
        -8814.715 / T
        + 172.0883
        - 27.927 * log(T)
        + (-160.340 / T + 1.3566) * S**0.5
        + (0.37335 / T - 0.05778) * S
    )
    KP3 = exp(
        -3070.75 / T
        - 18.141
        + (17.27039 / T + 2.81197) * S**0.5
        + (-44.99486 / T - 0.09984) * S
    )
    # SW vs NaCl for KW, K1, and K2
    # # KW
    if Sol == 1:
        KW = KW_in_NaCl(I, t)  # by "necessity" on H+(free)
    elif Sol == 3:
        KW = KW_in_KCl(I, t)  # by "necessity" on H+(free)
    else:
        KW = exp(
            -13847.26 / T
            + 148.9652
            - 23.652 * log(T)
            + (118.67 / T - 5.977 + 1.0495 * log(T)) * S**0.5
            - 0.01615 * S
        )  # used with H+(tot) --> not the original equation

    # # K1, K2
    if Sol == 1 or Sol == 3:
        m = I * 1000 / (1000 - I * (22.99 + 35.45))
        A1 = 35.2911 * m**0.5 + 0.8491 * m - 0.32 * m**1.5 + 0.055 * m**2
        B1 = -1583.09 * m**0.5
        C1 = -5.4366 * m**0.5
        A2 = 38.2746 * m**0.5 + 1.6057 * m - 0.647 * m**1.5 + 0.113 * m**2
        B2 = -1738.16 * m**0.5
        C2 = -6.0346 * m**0.5
        pK1 = (
            -402.56788
            + 11656.46 / T
            + 72.173 * log(T)
            - 0.161325 * T
            + 7.5526e-5 * T**2
        )
        pK2 = -122.4994 + 5811.18 / T + 20.5263 * log(T) - 0.0120897 * T
        K1 = 10 ** -(A1 + B1 / T + C1 * log(T) + pK1)  # by necessity on H+(free)
        K2 = 10 ** -(A2 + B2 / T + C2 * log(T) + pK2)  # by necessity on H+(free)
    else:
        K1 = 10 ** (
            -3633.86 / T + 61.2172 - 9.67770 * log(T) + 0.011555 * S - 0.0001152 * S**2
        )  # used with H+(tot)
        K2 = 10 ** (
            -471.78 / T - 25.9290 + 3.16967 * log(T) + 0.01781 * S - 0.0001122 * S**2
        )  # used with H+(tot)

    return [k, KW, K1, K2, ST, FT, BT, KS, KF, KB, KSi, KP1, KP2, KP3]


def KW_in_NaCl(c, t) -> float:
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
        + p10 * c
        + p01 * t
        + p20 * c**2
        + p11 * c * t
        + p02 * t**2
        + p30 * c**3
        + p21 * c**2 * t
        + p12 * c * t**2
        + p03 * t**3
    )


def KW_in_KCl(c, t) -> float:

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
        + p10 * c
        + p01 * t
        + p20 * c**2
        + p11 * c * t
        + p02 * t**2
        + p30 * c**3
        + p21 * c**2 * t
        + p12 * c * t**2
        + p03 * t**3
    )
