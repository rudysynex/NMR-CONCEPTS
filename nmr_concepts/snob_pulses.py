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
PULSE_FUNCTION = "dsnob"  # see below for availlable funcs
FREQ_OFFSET = 0
INIT_PHASE = 0
SAVE_PULSE = False
BASE_PATH = "C:/Users/RudrakshaMajumdar/Documents/Python Scripts/shaped_rf_tool/saved_rf_pulses"
NAME_PULSE = "esnob"


def snob_rf(
    func=PULSE_FUNCTION,
    pulse_length=PULSE_LENGTH,
    shape_pts=NR_OF_POINTS,
    freq_offset=FREQ_OFFSET,
    init_ph=INIT_PHASE,
):
    """
    Generate SNOB (Selective excitatioN fOr Biochemical applications)
    pulses given the following parameters:
        func = Modulation function
        pulse_length = duration of RF pulse in us
        shape_pts = No. of points to use for the shape
        freq_offset = Offset frequency of the pulse in Hz
        init_ph = Initial phase of the pulse

        Available functions: esnob, i2snob, i3snob, rdnob, dsnob

        NOTE: These pulses are generated with the B1max (khz)
        already calculated. 

        References: 
        Kupce, E., Boyd, J., & Campbell, I. D. (1995). 
        Short selective pulses for biochemical applications. 
        Journal of Magnetic Resonance, Series B, 106(3), 300-303.
    """

    # Time resolution
    time_res = pulse_length / (shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    if func == "esnob":

        a0 = 0.7500
        an = np.array(
            [
                -0.6176,
                -0.0373,
                -0.0005,
                -0.0182,
                -0.0058,
                -0.0036,
                -0.0051,
                -0.0031,
                -0.0017,
            ]
        )
        bn = np.array(
            [
                -0.4855,
                0.1260,
                -0.0191,
                -0.0005,
                -0.0003,
                0.0017,
                -0.0013,
                0.0001,
                -0.0025,
            ]
        )

    elif func == "i2snob":

        a0 = 0.5000
        an = np.array(
            [
                -0.2687,
                -0.2972,
                0.0989,
                -0.0010 - 0.0168,
                0.0009,
                -0.0017,
                -0.0013,
                -0.0014,
            ]
        )
        bn = np.array(
            [
                -1.1461,
                0.4016,
                0.0736,
                -0.0307,
                0.0079,
                0.0062,
                0.0003,
                -0.0002,
                0.0009,
            ]
        )

    elif func == "i3snob":

        a0 = 0.5000
        an = np.array(
            [
                0.2801,
                -0.9995,
                0.1928,
                0.0967,
                -0.0480,
                -0.0148,
                0.0088 - 0.0002,
                -0.0030,
            ]
        )
        bn = np.array(
            [
                -1.1990,
                0.4893,
                0.2439,
                -0.0816,
                -0.0409,
                0.0234,
                0.0036,
                -0.0042,
                0.0001,
            ]
        )

    elif func == "rsnob":

        a0 = 0.5000
        an = np.array([-1.1472, 0.5572, -0.0829, 0.0525])
        bn = np.zeros(len(an))

    elif func == "dsnob":

        a0 = 0.5000
        an = np.array(
            [-1.0662, 0.4466, 0.1673, -0.0049, -0.0753, 0.0001, 0.0144, 0.0041]
        )
        bn = np.array(
            [-2.4513, -0.2442, 0.6025, 0.1362, -0.0521, -0.0210, 0.0070, 0.0079]
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
    rf, freq_mod, phs_mod, time = snob_rf()
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
