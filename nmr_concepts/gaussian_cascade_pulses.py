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
PULSE_FUNCTION = "g4_excite"  # see below for availlable funcs
FREQ_OFFSET = 0
INIT_PHASE = 0

SAVE_PULSE = False
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/Python Scripts/shaped_rf_tool/saved_rf_pulses"
NAME_PULSE = "HERMITE90"


def gaussian_cascade_rf(
    func=PULSE_FUNCTION,
    pulse_length=PULSE_LENGTH,
    shape_pts=NR_OF_POINTS,
    freq_offset=FREQ_OFFSET,
    init_ph=INIT_PHASE,
):
    """
    Generate an Gaussian Cascade
    pulses given the following parameters:
        func = Modulation function
        pulse_length = duration of RF pulse in us
        freq_offset = Offset frequency of the pulse in Hz
        init_ph = Initial phase of the pulse

        Available functions: g3_invert, g4_excite

        Reference: Emsley, L., & Bodenhausen, G. (1990). 
        Gaussian pulse cascades: New analytical functions 
        for rectangular selective inversion and in-phase 
        excitation in NMR. Chemical Physics Letters, 165(6),
        469–476. doi:10.1016/0009-2614(90)87025-m 
    """

    # Time resolution
    time_res = pulse_length / (shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    if func == "g3_invert":
        e = np.exp

        # paramters for q3 k=1
        t_max_1 = 0.287 * pulse_length  # tmax/tp
        w_max_1 = -1.00  # relative amp
        fwhm_1 = 0.189 * pulse_length  # fwhm/tp
        a1 = np.log(2) / (fwhm_1 / 2) ** 2
        # paramters for q3 k=2
        t_max_2 = 0.508 * pulse_length  # tmax/tp
        w_max_2 = 1.37  # relative amp
        fwhm_2 = 0.183 * pulse_length  # fwhm/tp
        a2 = np.log(2) / (fwhm_2 / 2) ** 2
        # paramters for q3 k=3
        t_max_3 = 0.795 * pulse_length  # tmax/tp
        w_max_3 = 0.49  # relative amp
        fwhm_3 = 0.243 * pulse_length  # fwhm/tp
        a3 = np.log(2) / (fwhm_3 / 2) ** 2

        g1 = w_max_1 * e(-a1 * (time - t_max_1) ** 2)
        g2 = w_max_2 * e(-a2 * (time - t_max_2) ** 2)
        g3 = w_max_3 * e(-a3 * (time - t_max_3) ** 2)

        g3_invert = (g1 + g2 + g3) / np.amax(g1 + g2 + g3)

        rf_form = g3_invert

    elif func == "g4_excite":
        e = np.exp

        # paramters for q5 k=1
        t_max_1 = 0.177 * pulse_length  # tmax/tp
        w_max_1 = 0.62  # relative amp
        fwhm_1 = 0.172 * pulse_length  # fwhm/tp
        a1 = np.log(2) / (fwhm_1 / 2) ** 2
        # paramters for q5 k=2
        t_max_2 = 0.492 * pulse_length  # tmax/tp
        w_max_2 = 0.72  # relative amp
        fwhm_2 = 0.129 * pulse_length  # fwhm/tp
        a2 = np.log(2) / (fwhm_2 / 2) ** 2
        # paramters for q5 k=3
        t_max_3 = 0.653 * pulse_length  # tmax/tp
        w_max_3 = -0.91  # relative amp
        fwhm_3 = 0.119 * pulse_length  # fwhm/tp
        a3 = np.log(2) / (fwhm_3 / 2) ** 2
        # paramters for q5 k=4
        t_max_4 = 0.892 * pulse_length  # tmax/tp
        w_max_4 = -0.33  # relative amp
        fwhm_4 = 0.139 * pulse_length  # fwhm/tp
        a4 = np.log(2) / (fwhm_4 / 2) ** 2

        g1 = w_max_1 * e(-a1 * (time - t_max_1) ** 2)
        g2 = w_max_2 * e(-a2 * (time - t_max_2) ** 2)
        g3 = w_max_3 * e(-a3 * (time - t_max_3) ** 2)
        g4 = w_max_4 * e(-a4 * (time - t_max_4) ** 2)

        g4_excite = (g1 + g2 + g3 + g4) / np.amax(g1 + g2 + g3 + g4)

        rf_form = g4_excite

    else:
        print("ERROR: Please check 'func' argument")

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

    # fold back phase when it exceeds 360
    phs_mod = phs_mod % 360

    return (rf_form, freq_mod, phs_mod, time)


if __name__ == "__main__":
    # Example plot
    rf, freq_mod, phs_mod, time = gaussian_cascade_rf()
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
