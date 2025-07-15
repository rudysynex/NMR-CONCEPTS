#!/usr/bin/env python3
# Copyright: Synex Medical Inc, 2021
"""
Author: RD Majumdar
Date: 2022Nov01
Brief: Function to generate BURP pulses
"""

import json
import os

import matplotlib.pyplot as plt
import numpy as np

# Pulse Parameters
PULSE_LENGTH = 500  # us
NR_OF_POINTS = 1000
PULSE_FUNCTION = "qsneeze"  # see below for availlable funcs
FREQ_OFFSET = 0
INIT_PHASE = 0
SAVE_PULSE = False
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/Python Scripts/shaped_rf_tool/saved_rf_pulses"
NAME_PULSE = "iburp2_invert"


def burp_rf(
    func=PULSE_FUNCTION,
    pulse_length=PULSE_LENGTH,
    shape_pts=NR_OF_POINTS,
    freq_offset=FREQ_OFFSET,
    init_ph=INIT_PHASE,
):
    """
    Generate BURP (Band-selective Uniform Response Pure-phase)
    pulses given the following parameters:
        func = Modulation function
        pulse_length = duration of RF pulse in us
        shape_pts = No. of points to use for the shape
        freq_offset = Offset frequency of the pulse in Hz
        init_ph = Initial phase of the pulse

        Available functions: reburp_univ_180, uburp_univ_90
            iburp2_invert,eburp2_excite, qsneeze

        NOTE: These pulses are generated with the B1max (khz)
        already calculated. 

        References: 
        Geen, H., & Freeman, R. (1991). Band-selective 
        radiofrequency pulses. J. Mag. Res. (1969),
        93(1), 93–141. doi:10.1016/0022-2364(91)90034-q

        Kupce, E., & Freeman, R. (1995). Band-Selective 
        Correlation Spectroscopy. J. Mag. Res., Series A, 
        112(1), 134–137. doi:10.1006/jmra.1995.1023 
    """

    # Time resolution
    time_res = pulse_length / (shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    if func == "reburp_univ_180":

        a0 = 0.49
        an = np.array(
            [
                -1.02,
                1.11,
                -1.57,
                0.83,
                -0.42,
                0.26,
                -0.16,
                0.10,
                -0.07,
                0.04,
                -0.03,
                0.01,
                -0.02,
                0.00,
                -0.01,
            ]
        )
        bn = np.zeros(15)

    elif func == "uburp_univ_90":

        a0 = 0.27
        an = np.array(
            [
                -1.42,
                -0.37,
                -1.84,
                4.40,
                -1.19,
                0.00,
                -0.37,
                0.50,
                -0.31,
                0.18,
                -0.21,
                0.23,
                -0.12,
                0.07,
                -0.06,
                0.06,
                -0.04,
                0.03,
                -0.02,
                0.02,
            ]
        )
        bn = np.zeros(20)

    elif func == "iburp2_invert":

        a0 = 0.5
        an = np.array(
            [
                0.81,
                0.07,
                -1.25,
                -0.24,
                0.07,
                0.11,
                0.05,
                -0.02,
                -0.03,
                -0.02,
                0.00,
            ]
        )
        bn = np.array(
            [
                -0.68,
                -1.38,
                0.20,
                0.45,
                0.23,
                0.05,
                -0.04,
                -0.04,
                0.00,
                0.01,
                0.01,
            ]
        )

    elif func == "eburp2_excite":

        a0 = 0.26
        an = np.array(
            [0.91, 0.29, -1.28, -0.05, 0.04, 0.02, 0.06, 0.00, -0.02, 0.00]
        )
        bn = np.array(
            [-0.16, -1.82, 0.18, 0.42, 0.07, 0.07, -0.01, -0.04, 0.00, 0.00]
        )

    elif func == "qsneeze":

        a0 = 0.25
        an = np.array(
            [
                0.934,
                0.180,
                -1.527,
                0.003,
                0.143,
                0.050,
                0.072,
                -0.015,
                -0.040,
                -0.005,
            ]
        )
        bn = np.array(
            [
                -0.197,
                -1.772,
                0.204,
                0.619,
                0.076,
                0.039,
                -0.025,
                -0.060,
                0.005,
                0.017,
            ]
        )

    else:
        print("ERROR: Please check 'func' argument")

    two_pi = 2 * np.pi
    pulse_length = pulse_length * 1e-6  # convert to seconds
    w = two_pi / pulse_length

    amp_out = np.zeros(shape_pts)
    for step_nr in range(shape_pts):
        c = 0
        s = 0

        sh_time = step_nr / shape_pts * pulse_length
        for Fn in range(len(an)):
            c = c + an[Fn] * np.cos((Fn + 1) * sh_time * w)
            s = s + bn[Fn] * np.sin((Fn + 1) * sh_time * w)

        total = c + s + a0
        total_b1 = total * w / two_pi / 1000

        amp_out[step_nr] = total_b1

    rf_form = amp_out

    phs_mod = np.zeros(shape_pts)

    for idx in range(shape_pts):
        if rf_form[idx] < 0:
            phs_mod[idx] = 180
        else:
            phs_mod[idx] = 0

    # Frequency modulation
    freq_mod = np.zeros(shape_pts)
    freq_mod[0:shape_pts] = freq_offset

    # Add initial phase offset
    phs_mod = phs_mod + init_ph

    # fold back phase when it exceeds 360
    phs_mod = phs_mod % 360

    return (rf_form, freq_mod, phs_mod, time)


if __name__ == "__main__":
    # Example plot
    rf, freq_mod, phs_mod, time = burp_rf()
    # Convert degrees to radians for complex plot
    phs_mod_rad = np.deg2rad(phs_mod)

    # Plot the shape
    # rf = rf/np.amax(rf)
    fig = plt.figure(figsize=[8, 8])
    plt_rfamp = fig.add_subplot(221)
    plt_freq_mod = fig.add_subplot(222)
    plt_phs_mod = fig.add_subplot(223)
    plt_complex = fig.add_subplot(224)

    plt_rfamp.plot(time, rf, "b")
    plt_rfamp.set_ylabel("RF amplitude (khz)")
    plt_rfamp.set_xlabel("Time (us)")

    plt_freq_mod.plot(time, freq_mod, "b")
    plt_freq_mod.set_ylabel("Frequency (Hz)")
    plt_freq_mod.set_xlabel("Time (us)")

    plt_phs_mod.plot(time, phs_mod, "b")
    plt_phs_mod.set_ylabel("Phase (deg.)")
    plt_phs_mod.set_xlabel("Time (us)")

    plt_complex.plot(time, rf * np.cos(phs_mod_rad), "b", label="B1x")
    plt_complex.plot(time, rf * np.sin(phs_mod_rad), "r", label="B1y")
    plt_complex.set_ylabel("RF amplitude (khz)")
    plt_complex.set_xlabel("Time (us)")
    # plt_complex.set_ylim(-1, 1)
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
