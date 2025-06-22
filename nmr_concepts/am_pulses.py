#!/usr/bin/env python3
# Copyright: Synex Medical Inc, 2021
"""
Author: RD Majumdar
Date: 2022Sep16
Brief: Function to generate Amplitude modulated pulses
"""

import json
import os

import matplotlib.pyplot as plt
import numpy as np

# Pulse Parameters
PULSE_LENGTH = 1000  # us
NR_OF_POINTS = 256
PULSE_FUNCTION = "gaussian"  # see below for availlable funcs
SINC_TBW = 6
GAUSSIAN_TRUNCATION = 1  # in %
HERMITE_NUT_ANGLE = 90
FREQ_OFFSET = 0
INIT_PHASE = 0

SAVE_PULSE = False
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/Python Scripts/shaped_rf_tool/saved_rf_pulses"
NAME_PULSE = "HERMITE90"


def am_rf(
    func=PULSE_FUNCTION,
    pulse_length=PULSE_LENGTH,
    shape_pts=NR_OF_POINTS,
    lobes=SINC_TBW,
    trunc_lvl=GAUSSIAN_TRUNCATION,
    nut_angle=HERMITE_NUT_ANGLE,
    freq_offset=FREQ_OFFSET,
    init_ph=INIT_PHASE,
):
    """
    Generate an Amplitude Modulated (AM)
    pulse given the following parameters:
        func = Modulation function
        pulse_length = duration of RF pulse in us
        shape_pts = No. of points to use for the shape
        lobes = Number of lobes in the sinc shape
        trunc_lvl = Truncation level % for Gaussian pulse.
                    1 or 10 are standard
        nut_angle = 90 or 180. Nutation angle for Hermite pulse
        freq_offset = Offset frequency of the pulse in Hz
        init_ph = Initial phase of the pulse

        Available functions: sinc, square, gaussian, hermite, q3, 
                q5_universal_90, reburp_univ_180, uburp_univ_90

        Reference for Q3 and Q5: Emsley, L., & Bodenhausen, G. (1992).
            Optimization of shaped selective pulses for NMR using a
            quaternion description of their overall propagators. 
            J. Mag. Res. (1969), 97(1), 135–148. 
            doi:10.1016/0022-2364(92)90242-y
    """

    # Time resolution
    time_res = pulse_length / (shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    if func == "sinc":

        # Define apparent time, from -1 to +1
        t_vec = np.linspace(-1, 1, shape_pts)
        # Define the shape
        rf_form = np.sin(np.pi * lobes * 0.5 * t_vec) / t_vec

        # Normalize so that max val. is 1
        rf_form = rf_form / np.amax(rf_form)

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset

        # Calculate phase mod from freq mod
        # phs_mod = 360 * freq_offset * (time / 1e6)
        phs_mod = np.zeros(shape_pts)

        for idx in range(shape_pts):
            if rf_form[idx] < 0:
                phs_mod[idx] = 180
            else:
                phs_mod[idx] = 0

        # Add initial phase offset
        phs_mod = phs_mod + init_ph

    elif func == "gaussian":

        # Define apparent time, from -1 to +1
        t_vec = np.linspace(-1, 1, shape_pts)

        # Convert truncation level (%) to exponential scaling
        gaussian_scale = 1.0 * np.log(trunc_lvl / 100.0)

        # Define the shape
        rf_form = np.exp(gaussian_scale * t_vec * t_vec)

        # Normalize so that max val. is 1
        rf_form = rf_form / np.amax(rf_form)

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset

        # Calculate phase mod from freq mod
        # phs_mod = 360 * freq_offset * (time / 1e6)
        phs_mod = np.zeros(shape_pts)

        for idx in range(shape_pts):
            if rf_form[idx] < 0:
                phs_mod[idx] = 180
            else:
                phs_mod[idx] = 0

        # Add initial phase offset
        phs_mod = phs_mod + init_ph

    elif func == "hermite":

        # Define apparent time, from -1 to +1
        t_vec = np.linspace(-1, 1, shape_pts)

        # Define the shape
        if nut_angle == 90:
            rf_form = (1.0 - 6.003 * t_vec * t_vec) * np.exp(
                -9.0 * t_vec * t_vec
            )
        elif nut_angle == 180:
            rf_form = (1.0 - 8.604 * t_vec * t_vec) * np.exp(
                -9.0 * t_vec * t_vec
            )
        else:
            print("nut_angle can only be 90 or 180")

        # Normalize so that max val. is 1
        rf_form = rf_form / np.amax(rf_form)

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset

        # Calculate phase mod from freq mod
        # phs_mod = 360 * freq_offset * (time / 1e6)
        phs_mod = np.zeros(shape_pts)

        for idx in range(shape_pts):
            if rf_form[idx] < 0:
                phs_mod[idx] = 180
            else:
                phs_mod[idx] = 0

        # Add initial phase offset
        phs_mod = phs_mod + init_ph

    elif func == "square":

        # Define apparent time, from -1 to +1
        t_vec = np.linspace(-1, 1, shape_pts)

        # Define the shape
        rf_form = np.linspace(1, 1, shape_pts)

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset

        # Calculate phase mod from freq mod
        phs_mod = 360 * freq_offset * (time / 1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph

    elif func == "q3_universal_180":
        e = np.exp

        # paramters for q3 k=1
        t_max_1 = 0.306 * pulse_length  # tmax/tp
        w_max_1 = -4.39  # relative amp
        fwhm_1 = 0.18 * pulse_length  # fwhm/tp
        a1 = np.log(2) / (fwhm_1 / 2) ** 2
        # paramters for q3 k=2
        t_max_2 = 0.545 * pulse_length  # tmax/tp
        w_max_2 = 4.57  # relative amp
        fwhm_2 = 0.183 * pulse_length  # fwhm/tp
        a2 = np.log(2) / (fwhm_2 / 2) ** 2
        # paramters for q3 k=3
        t_max_3 = 0.804 * pulse_length  # tmax/tp
        w_max_3 = 2.6  # relative amp
        fwhm_3 = 0.245 * pulse_length  # fwhm/tp
        a3 = np.log(2) / (fwhm_3 / 2) ** 2

        g1 = w_max_1 * e(-a1 * (time - t_max_1) ** 2)
        g2 = w_max_2 * e(-a2 * (time - t_max_2) ** 2)
        g3 = w_max_3 * e(-a3 * (time - t_max_3) ** 2)

        q3 = (g1 + g2 + g3) / np.amax(g1 + g2 + g3)

        rf_form = q3

        phs_mod = np.zeros(shape_pts)

        for idx in range(shape_pts):
            if rf_form[idx] < 0:
                phs_mod[idx] = 180
            else:
                phs_mod[idx] = 0

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset

        # # Calculate phase mod from freq mod
        # phs_mod = 360 * freq_offset * (time / 1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph

    elif func == "q5_universal_90":
        e = np.exp

        # paramters for q5 k=1
        t_max_1 = 0.162 * pulse_length  # tmax/tp
        w_max_1 = -1.48  # relative amp
        fwhm_1 = 0.186 * pulse_length  # fwhm/tp
        a1 = np.log(2) / (fwhm_1 / 2) ** 2
        # paramters for q5 k=2
        t_max_2 = 0.307 * pulse_length  # tmax/tp
        w_max_2 = -4.34  # relative amp
        fwhm_2 = 0.139 * pulse_length  # fwhm/tp
        a2 = np.log(2) / (fwhm_2 / 2) ** 2
        # paramters for q5 k=3
        t_max_3 = 0.497 * pulse_length  # tmax/tp
        w_max_3 = 7.33  # relative amp
        fwhm_3 = 0.143 * pulse_length  # fwhm/tp
        a3 = np.log(2) / (fwhm_3 / 2) ** 2
        # paramters for q5 k=4
        t_max_4 = 0.525 * pulse_length  # tmax/tp
        w_max_4 = -2.3  # relative amp
        fwhm_4 = 0.29 * pulse_length  # fwhm/tp
        a4 = np.log(2) / (fwhm_4 / 2) ** 2
        # paramters for q5 k=5
        t_max_5 = 0.803 * pulse_length  # tmax/tp
        w_max_5 = 5.66  # relative amp
        fwhm_5 = 0.137 * pulse_length  # fwhm/tp
        a5 = np.log(2) / (fwhm_5 / 2) ** 2

        g1 = w_max_1 * e(-a1 * (time - t_max_1) ** 2)
        g2 = w_max_2 * e(-a2 * (time - t_max_2) ** 2)
        g3 = w_max_3 * e(-a3 * (time - t_max_3) ** 2)
        g4 = w_max_4 * e(-a4 * (time - t_max_4) ** 2)
        g5 = w_max_5 * e(-a5 * (time - t_max_5) ** 2)

        q5 = (g1 + g2 + g3 + g4 + g5) / np.amax(g1 + g2 + g3 + g4 + g5)

        rf_form = q5

        phs_mod = np.zeros(shape_pts)

        for idx in range(shape_pts):
            if rf_form[idx] < 0:
                phs_mod[idx] = 180
            else:
                phs_mod[idx] = 0

        # Frequency modulation
        freq_mod = np.zeros(shape_pts)
        freq_mod[0:shape_pts] = freq_offset

        # # Calculate phase mod from freq mod
        # phs_mod = 360 * freq_offset * (time / 1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph

    else:
        print("ERROR: Please check 'func' argument")

    # fold back phase when it exceeds 360
    phs_mod = phs_mod % 360

    return (rf_form, freq_mod, phs_mod, time)


if __name__ == "__main__":
    # Example plot
    rf, freq_mod, phs_mod, time = am_rf()
    # Convert degrees to radians for complex plot
    phs_mod_rad = np.deg2rad(phs_mod)

    # Plot the shape
    fig = plt.figure(figsize=[8, 8])
    plt_rfamp = fig.add_subplot(221)
    plt_freq_mod = fig.add_subplot(222)
    plt_phs_mod = fig.add_subplot(223)
    plt_complex = fig.add_subplot(224)

    plt_rfamp.plot(time, rf, "b")
    plt_rfamp.set_ylabel("RF amplitude")
    plt_rfamp.set_xlabel("Time (us)")

    plt_freq_mod.plot(time, freq_mod, "b")
    plt_freq_mod.set_ylabel("Frequency (Hz)")
    plt_freq_mod.set_xlabel("Time (us)")

    plt_phs_mod.plot(time, phs_mod, "b")
    plt_phs_mod.set_ylabel("Phase (deg.)")
    plt_phs_mod.set_xlabel("Time (us)")

    plt_complex.plot(time, rf * np.cos(phs_mod_rad), "b", label="B1x")
    plt_complex.plot(time, rf * np.sin(phs_mod_rad), "r", label="B1y")
    plt_complex.set_ylabel("RF amplitude")
    plt_complex.set_xlabel("Time (us)")
    plt_complex.set_ylim(-1, 1)
    plt_complex.legend()

    plt_rfamp.grid()
    plt_freq_mod.grid()
    plt_phs_mod.grid()
    plt_complex.grid()

    plt.tight_layout()
    plt.show()

    # Option to save the shape file
    if SAVE_PULSE is True:
        SAVE_PATH = os.path.join(BASE_PATH, NAME_PULSE)
        os.makedirs(SAVE_PATH, exist_ok=True)

        # save the shape file
        file_name = os.path.join(SAVE_PATH, "rf_pulse_file")
        np.savez(file_name, rf, phs_mod, freq_mod)

        # save the pulse parameters
        _global_dict = globals()
        _glob_names = list(_global_dict.keys())
        global_vars = {}

        for name in _glob_names:
            if not name.startswith("_") and isinstance(
                _global_dict[name], (str, bool, float, int, list, dict)
            ):
                global_vars[name] = _global_dict[name]

        params_file_name = os.path.join(SAVE_PATH, "rf_pulse_parameters.json")
        with open(params_file_name, "w+") as param_file:
            json.dump(global_vars, param_file, indent=4)
