#!/usr/bin/env python3
# Copyright: Synex Medical Inc, 2021
"""
Author: RD Majumdar
Date: 2022Sep16
Brief: Function to generate Adiabatic Full Passage pulses
"""

import json
import os

import matplotlib.pyplot as plt
import numpy as np

# Pulse Parameters
PULSE_LENGTH = 10000  # us
NR_OF_POINTS = 256
PULSE_FUNCTION = "HSn"  # "HSn" OR "tanh/tan" OR "lorentz"
HYPSEC_ORDER =4
FREQ_OFFSET = 0
INIT_PHASE = 0
TBW = 40

SAVE_PULSE = False
BASE_PATH = (
    "/store/src/pulse-sequence/pulse_sequence/libs/pulses/saved_waveforms"
)
NAME_PULSE = "AFP_HS4_4khz_2048"


def arcsech(val):
    """Inverse hyperbolic secant"""
    return np.arccosh(1.0 / val)


def sech(arg):
    """Hyperbolic secant"""
    return 1 / np.cosh(arg)


def afp_rf(
    pulse_length=PULSE_LENGTH,
    shape_pts=NR_OF_POINTS,
    func=PULSE_FUNCTION,
    power_n=HYPSEC_ORDER,
    k=1.52,
    tbw=TBW,
    freq_offset=FREQ_OFFSET,
    init_ph=INIT_PHASE,
):
    """
    Generate Adiabatic Full Passage pulses
    given the following parameters:
        pulse_length = duration of RF pulse in us
        shape_pts = No. of points to use for the shape
        func = Modulation function, "HSn" or "tanh/tan"
                Defaults to an HS1 pulse
        power_n = Power 'n' in Rf modulation. n=1 for HS1, n=8 for HS8 etc
        k = Freq modulation factor for tanh/tan, typically 1.52
        sweep_bw = frequency sweep bandwidth in Hz
        freq_offset = Offset frequency of the pulse in Hz
        init_ph = Initial phase of the pulse

    beta = Rf modulation truncation factor,
           typically 5.3 for HSn
                     10 for tanh/tan

        Notes:
        1. Number of points for AFP pulses refers to points per AHP segment.
        2. RF and phase modulations are forced to be time-symmetric.
    """
    shape_pts = int(shape_pts / 2)
    # Time resolution
    time_res = pulse_length / (2 * shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    sweep_bw = 1e6 * (tbw / pulse_length)

    if func == "HSn":

        # Offset-independent adiabaticity (OIA) pulses based on
        # hyperbolic secant functions raised to a power n.
        # n = 1 leads to the original hyperbolic secant RF pulses.

        # Original/relevant publications:
        # 1. M. S. Silver, R. I. Joseph, D. I. Hoult,
        #     J. Magn. Reson. 59, 347-351 (1984)
        # 2. A. Tannus, M. Garwood,
        #     J. Magn. Reson. 120A, 133-137 (1996)
        # See also Chapter 5 in 'In Vivo NMR spectroscopy' book
        # by Robin de Graaf

        beta = 5.3

        # Define apparent time, from 0 to +1
        t_vec = np.linspace(0, 1, shape_pts)

        # Define the shape (first segment)
        rf_form = sech(beta * (1.0 - t_vec) ** power_n)

        # Normalize so that max val. is 1
        rf_form = rf_form / np.amax(rf_form)

        # Frequency modulation
        freq_mod = np.cumsum(rf_form * rf_form)
        # Normalize freq_mod to 0 to +1 range
        freq_mod = 1.0 - (freq_mod / np.amax(freq_mod))
        # Scale the freq_mod to full bandwidth
        freq_mod = (sweep_bw / 2.0) * freq_mod

        # Phase modulation = time itegral of frequency
        phs_mod = 360 * np.cumsum(freq_mod * (time_res / 1e6))
        # Force the 1st pont in phs_mod to zero
        phs_mod = phs_mod - phs_mod[0]

        # Create the 2nd half of the pulse from symmetry
        # and join with first half

        rf_form = np.concatenate([rf_form, np.flip(rf_form)])
        freq_mod = np.concatenate([freq_mod, -1 * np.flip(freq_mod)])
        phs_mod = np.concatenate([phs_mod, np.flip(phs_mod)])

        # Add frequency offset to freq_mod
        freq_mod = freq_mod + freq_offset

        # Add freq offset as a phase ramp to phase modulation
        phs_mod = phs_mod + 360 * freq_offset * (time / 1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph
        rf_form = rf_form - rf_form[0]

    elif func == "lorentz":
        beta = 99

        # Define apparent time, from 0 to +1
        t_vec = np.linspace(0, 1, shape_pts)

        # Define the shape (first segment)
        tau = 1.0 - t_vec
        rf_form = 1 / (1 + beta * tau ** 2)

        # Normalize so that max val. is 1
        rf_form = rf_form / np.amax(rf_form)

        # Frequency modulation
        freq_mod = (tau / (1 + beta * tau ** 2)) + (
            1 / np.sqrt(beta)
        ) * np.arctan(np.sqrt(beta) * tau)
        # Normalize freq_mod to 0 to +1 range
        # freq_mod = 1.0 - (freq_mod / np.amax(freq_mod))
        # Scale the freq_mod to full bandwidth
        freq_mod = (sweep_bw / 2.0) * freq_mod

        # Phase modulation = time itegral of frequency
        phs_mod = 360 * np.cumsum(freq_mod * (time_res / 1e6))
        # Force the 1st pont in phs_mod to zero
        phs_mod = phs_mod - phs_mod[0]

        # Create the 2nd half of the pulse from symmetry
        # and join with first half

        rf_form = np.concatenate([rf_form, np.flip(rf_form)])
        freq_mod = np.concatenate([freq_mod, -1 * np.flip(freq_mod)])
        phs_mod = np.concatenate([phs_mod, np.flip(phs_mod)])

        # Add frequency offset to freq_mod
        freq_mod = freq_mod + freq_offset

        # Add freq offset as a phase ramp to phase modulation
        phs_mod = phs_mod + 360 * freq_offset * (time / 1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph
        rf_form = rf_form - rf_form[0]

    elif func == "tanh/tan":

        # Adiabatic full passage pulses based on tanh/tan modulation that
        # closely reflect numerically optimized modulation (NOM) functions.

        # Original/relevant publications:
        # 1. M. Garwood, Y. Ke,
        #    J. Magn. Reson. 94, 511-525 (1991)
        # 2. K. Ugurbil, M. Garwood, A. Rath,
        #     J. Magn. Reson. 80, 448-469 (1988)

        # See also Chapter 5 in 'In Vivo NMR spectroscopy' book
        # by Robin de Graaf

        beta = 10

        # Define apparent time, from 0 to +1
        t_vec = np.linspace(0, 1, shape_pts)

        # Define the shape (first segment)
        rf_form = np.tanh(beta * t_vec)

        # Normalize so that max val. is 1
        rf_form = rf_form / np.amax(rf_form)

        # Frequency modulation
        freq_mod = np.tan(k * (1.0 - t_vec))
        # Normalize freq_mod to 0 to +1 range
        freq_mod = freq_mod / np.amax(freq_mod)
        # Scale the freq_mod to full bandwidth
        freq_mod = (sweep_bw / 2.0) * freq_mod

        # Phase modulation = time itegral of frequency
        phs_mod = 360 * np.cumsum(freq_mod * (time_res / 1e6))
        # Force the 1st pont in phs_mod to zero
        phs_mod = phs_mod - phs_mod[0]

        # Create the 2nd half of the pulse from symmetry
        # and join with first half

        rf_form = np.concatenate([rf_form, np.flip(rf_form)])
        freq_mod = np.concatenate([freq_mod, -1 * np.flip(freq_mod)])
        phs_mod = np.concatenate([phs_mod, np.flip(phs_mod)])

        # Add frequency offset to freq_mod
        freq_mod = freq_mod + freq_offset

        # Add freq offset as a phase ramp to phase modulation
        phs_mod = phs_mod + 360 * freq_offset * (time / 1e6)

        # Add initial phase offset
        phs_mod = phs_mod + init_ph

    else:
        print("ERROR: Please make sure func = Hsn or tanh/tan ")

    # fold back phase when it exceeds 360
    phs_mod = phs_mod % 360

    return (rf_form, freq_mod, phs_mod, time)


if __name__ == "__main__":
    # Example plot
    rf, freq_mod, phs_mod, time = afp_rf()
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
