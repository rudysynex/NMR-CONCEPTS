#!/usr/bin/env python3
# Copyright: Synex Medical Inc, 2021
"""
Author: RD Majumdar
Date: 2022Nov14
Brief: Functions to generate composite pulses
"""

import matplotlib.pyplot as plt
import numpy as np

LIST_FLIP_ANGLES = np.array([90, 180, 90])
LIST_PHASE = np.array([0, 90, 0])
DURATION_REFERENCE_90 = 100.0
NUM_POINTS_REF_90 = 50
FREQ_OFFSET = 0


def rect_rf(pulse_length, shape_pts, freq_offset, init_ph):
    """
    rectangular rf pulse
    """

    # Time resolution
    time_res = pulse_length / (shape_pts)
    # Actual time axis
    time = np.arange(0, pulse_length, time_res)

    # Define the shape
    rf_form = np.linspace(1, 1, shape_pts)

    # Frequency modulation
    freq_mod = np.zeros(shape_pts)
    freq_mod[0:shape_pts] = freq_offset

    # Calculate phase mod from freq mod
    phs_mod = 360 * freq_offset * (time / 1e6)

    # Add initial phase offset
    phs_mod = phs_mod + init_ph

    # fold back phase when it exceeds 360
    phs_mod = phs_mod % 360

    return (rf_form, phs_mod, freq_mod)


def composite_rf(
    list_flip_angle=LIST_FLIP_ANGLES,
    list_phase_angle=LIST_PHASE,
    duration_ref_90=DURATION_REFERENCE_90,
    pnts_ref_90=NUM_POINTS_REF_90,
    freq_offset=FREQ_OFFSET,
):

    """
    Generate composite pulse
    """
    list_duration = np.zeros(len(list_flip_angle))
    list_pnts = np.zeros(len(list_flip_angle))
    for i in range(len(list_flip_angle)):
        list_duration[i] = duration_ref_90 * list_flip_angle[i] / 90
        list_pnts[i] = pnts_ref_90 * list_flip_angle[i] / 90

    duration_total = np.sum(list_duration)

    am_list = []
    pm_list = []
    fm_list = []
    for i in range(len(list_flip_angle)):
        am, pm, fm = rect_rf(
            list_duration[i],
            int(list_pnts[i]),
            freq_offset,
            list_phase_angle[i],
        )
        am_list.append(am)
        pm_list.append(pm)
        fm_list.append(fm)
    am_tot = np.concatenate(np.asarray(am_list, dtype=object)).astype(float)
    pm_tot = np.concatenate(np.asarray(pm_list, dtype=object)).astype(float)
    fm_tot = np.concatenate(np.asarray(fm_list, dtype=object))
    time = np.linspace(0, duration_total, am_tot.size)
    return (am_tot, fm_tot, pm_tot, time, duration_total)


if __name__ == "__main__":
    # Example plot
    rf, _, phs_mod, time, _ = composite_rf()
    # Convert degrees to radians for complex plot
    phs_mod_rad = np.deg2rad(phs_mod)

    # Plot the shape
    fig = plt.figure()
    plt_rfamp = fig.add_subplot(311)
    plt_phs_mod = fig.add_subplot(312)
    plt_complex = fig.add_subplot(313)

    plt_rfamp.plot(time, rf, "b")
    plt_rfamp.set_ylabel("RF amplitude")
    plt_rfamp.set_xlabel("Time (us)")

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
    plt_phs_mod.grid()
    plt_complex.grid()

    plt.tight_layout()
    plt.show()
