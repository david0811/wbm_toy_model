import numpy as np
import math


def calculate_potential_evapotranspiration(T, N, phi):
    """
    Inputs:
     - T: daily mean temperature [C]
     - N: day of year
     - phi: latitude [deg]

    Outputs:
     - daily potential evapotranspiration calculated via the Hamon method [mm]

    Notes: (e.g.) http://data.snap.uaf.edu/data/Base/AK_2km/PET/Hamon_PET_equations.pdf
    """
    # Calculate solar declination (delta)
    delta = -23.44 * math.cos(math.radians((360 / 365) * (N + 10)))

    # Calculate fractional day length (Lambda)
    Lambda = (1 / math.pi) * math.acos(
        -math.tan(math.radians(phi)) * math.tan(math.radians(delta))
    )

    # Calculate saturation vapor density (rho_sat)
    Psat = calculate_saturation_vapor_pressure(T)
    rho_sat = (2.167 * Psat) / (T + 273.15)

    # Calculate potential evapotranspiration (PET)
    PET = 330.2 * Lambda * rho_sat
    return PET


def calculate_saturation_vapor_pressure(T):
    """
    Inputs:
     - T: daily mean temperature [C]

    Outputs:
     - saturation vapor pressure via the Tetens equation [kPa]

    Notes: https://en.wikipedia.org/wiki/Tetens_equation
    """
    if T >= 0:
        Psat = 0.61078 * np.exp((17.26939 * T) / (T + 237.3))
    else:
        Psat = 0.61078 * np.exp((21.87456 * T) / (T + 265.5))
    return Psat


def calculate_canopy_evaporation(Wi, Wi_max, T):
    """
    Inputs:
     - Wi: canopy water storage [mm]
     - Wi_max: maximum canopy water storage [mm]
     - T: daily mean temperature [C]

    Outputs:
     - daily canopy evaporation [mm]
    """

    # Reference evaporation rate
    # This is different from WBM!
    # WBM calculates based on many meteo inputs
    if T < 0:
        Eow = 0.0
    else:
        Eow = 0.36

    Ec = Eow * ((Wi / Wi_max) ** 0.6666667)
    return Ec


def soil_drying_function(Ws, Wcap, alpha):
    """
    Inputs:
     - Ws: soil moisture content [mm]
     - Wcap: maximum soil moisture content [mm]
     - alpha: scaling parameter

    Outputs:
     - soil drying parameter (restricted [0-1])
    """
    g = (1 - np.exp(-alpha * Ws / Wcap)) / (1 - np.exp(-alpha))
    return g


def simulate_water_balance(
    spinup,
    Ws_init,
    Wi_init,
    Sp_init,
    Wg_init,
    P,
    T,
    Wcap_in,
    Ts,
    Tm,
    lai,
    Kc,
    alpha,
    s_c,
    r_I,
    r_FI,
    r_p,
    beta_HBV,
    gamma_S,
    gamma_G,
    phi,
):
    """
    Inputs:
     - spinup: number of years to loop though (using same input variables)
     - Ws_init: initial soil moisture content [mm]
     - Wi_init: initial canopy water storage [mm]
     - Sp_init: initial snowpack [mm]
     - Wg_init: initial groundwater [mm]
     - P: daily precipitation timeseries [mm]
     - T: daily mean temperature timeseries [C]
     - Wcap_in: maximum soil moisture content timeseries [mm]
     - Ts: snowfall threshold [C]
     - Tm: snowmelt threshold [C]
     - lai: leaf area index timeseries []
     - Kc: crop scalar factor timeseries []
     - alpha: drying function scale parameter []
     - s_c: Crop-specific threshold factor []
     - r_I: Irrigation inefficiency factor []
     - r_FI: Framer irrigation inefficiency factor [] (added by David to test some things, in WBM r_FI = 1)
     - r_p: Runoff-percolation partitioning factor []
     - beta_HBV: HBV direct recharge parameter [] (-1 to turn off)
     - gamma_S: Soil moisture surplus coefficient []
     - gamma_G: Groundwater leakage coefficient []
     - phi: Latitude [deg]

    Outputs:
     - Ws: daily soil moisture content timeseries [mm]
     - Ws_frac: daily soil moisture fraction timeseries []
     - Wg: daily groundwater timeseries [mm]
    """

    # Simulation length
    n_sim = len(P) * spinup

    # Initialize variables
    Ws = np.zeros(n_sim)  # Soil moisture
    Ws[0] = Ws_init

    Ws_frac = np.zeros(n_sim)  # Soil moisture fraction
    Ws_frac[0] = Ws_init / Wcap_in[0]

    Wi = np.zeros(n_sim)  # Canopy water storage
    Wi[0] = Wi_init

    Sp = np.zeros(n_sim)  # Snowpack
    Sp[0] = Sp_init

    Wg = np.zeros(n_sim)  # Groundwater volume
    Wg[0] = Wg_init

    Pt = 0  # Throughfall
    AET = 0  # Actual evapotranspiration
    S = 0  # Storage

    Inet = 0  # Net irrigation
    Igross = 0  # Gross irrigation
    Enb = 0  # Non-beneficial evaporation
    Rperc = 0  # Percolation

    # Can be used to look at intermediate variables
    #     S_out = []
    #     Pt_out = []
    #     Wi_out = []
    #     AET_out = []
    #     PET_out = []
    #     Ec_out = []
    #     g_out = []
    #     Ms_out = []

    for tt in range(1, n_sim):
        # tt runs from 1 to n_sim (used for output variables)
        # t runs from 0 to 365 (used for input variables since we only have 1 year)
        t = tt % 365

        # Soil moisture cap
        Wcap = Wcap_in[t]
        W_ct = s_c * Wcap  # Crop optimal soil moisture

        # Snowfall
        if T[t] < Ts:
            Ps = P[t]
            Pa = 0
            Sp[tt] = Sp[tt - 1] + Ps
        else:
            Pa = P[t]
            Sp[tt] = Sp[tt - 1]

        # Snowmelt
        if T[t] > Tm:
            Ms = 2.63 + 2.55 * T[t] + 0.0912 * T[t] * Pa
            if Ms > Sp[tt]:
                Ms = Sp[tt]
                Sp[tt] = 0
            else:
                Sp[tt] = Sp[tt] - Ms
        else:
            Ms = 0.0

        # Calculate throughfall
        Wi_max = 0.25 * lai[t]
        Ec = calculate_canopy_evaporation(Wi[tt - 1], Wi_max, T[t])

        if Wi_max < Pa + Wi[tt - 1] - Ec:
            Pt = Pa - Ec - (Wi_max - Wi[tt - 1])
        else:
            Pt = 0

        # Update canopy storage
        if Wi[tt - 1] + (Pa - Pt) - Ec <= Wi_max:
            if Wi[tt - 1] + (Pa - Pt) - Ec > 0.0:
                Wi[tt] = Wi[tt - 1] + (Pa - Pt) - Ec
            else:
                Wi[tt] = 0.0
        else:
            Wi[tt] = Wi_max

        # Calculate actual evapotranspiration
        PET = Kc[t] * calculate_potential_evapotranspiration(T[t], t, phi)
        if Pt >= PET:
            AET = PET
        elif Pt < PET:
            if Ws[tt - 1] > 0:
                AET = soil_drying_function(Ws[tt - 1], Wcap, alpha) * (PET - Pt - Ms)
            else:
                AET = 0

        # Direct groundwater recharge (HBV)
        if beta_HBV > 0.0:
            # Id = (Pt + Ms) * (Ws[tt - 1] / Wcap) ** beta_HBV
            Id = 154 * (Ws[tt - 1] / Wcap) ** beta_HBV
        else:
            Id = 0.0

        # Update soil moisture
        if Wcap < Ws[tt - 1] + (Pt - Id) - AET:
            S = Ws[tt - 1] + (Pt - Id) - AET - Wcap
        else:
            S = 0

        Ws[tt] = Ws[tt - 1] + (Pt - Id) - AET - S
        Ws_frac[tt] = Ws[tt] / Wcap

        if Ws[tt] < 0:
            Ws[tt] = 0.0

        # # Align growing season
        # if t == gs_start or t == gs_end:
        #     Ws[t] = soilMoist[t]

        # Calculate irrigation net amount
        if Ws[tt] < W_ct:
            Inet = Wcap - Ws[tt]
        else:
            Inet = 0

        # Calculate farmer net irrigation amount
        Inet = Inet * r_FI

        # Calculate gross irrigation amount
        Igross = Inet / r_I

        # Update soil moisture with irrigation
        # Ws[t] = Ws[t] + Inet

        # Calculate non-beneficial evaporation
        Enb = min(
            calculate_potential_evapotranspiration(T[t], t, phi) - AET, Igross - Inet
        )

        # Calculate percolation
        Rperc = r_p * (Igross - Inet - Enb)

        # Calculate surface runoff
        R_ro = (1 - r_p) * (Igross - Inet - Enb)

        # Update groundwater volume
        Wg[tt] = Wg[tt - 1] + gamma_S * S - gamma_G * Wg[tt - 1] - Igross + Rperc + Id
        if Wg[tt] < 0.0:
            Wg[tt] = 0.0

        # # Append all
        # S_out.append(S)
        # Wi_out.append(Wi[t])
        # Pt_out.append(Pt)
        # AET_out.append(AET)
        # PET_out.append(PET)
        # Ec_out.append(Ec)
        # g_out.append(soil_drying_function(Ws[t-1], Wcap, alpha))
        # Ms_out.append(Ms)

    # Return
    return Ws, Ws_frac, Wg
