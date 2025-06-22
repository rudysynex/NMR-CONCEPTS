import numpy as np
from qutip import *


def simulate_amp(
    pulse_length, sim_points, rf, pm, rf_min, rf_max, freq_offset, init_mag
):
    """
    A function to simulate magnetization as a function 
    of RF amplitude. In the current version, the pulse
    needs to be pre-loaded outside the function.
    """

    rf_amp = rf
    rf_dwelltime = 1e-6 * (pulse_length / rf.size)  # convert to seconds

    # Convert rf phase from degrees to radians
    rf_phase = (np.pi / 180) * pm

    # Convert freq offset from Hz to radians/sec
    rf_offset = 2 * np.pi * freq_offset

    # Initial magnetizations
    Mx_0 = init_mag[0]
    My_0 = init_mag[1]
    Mz_0 = init_mag[2]

    # x-axis step size
    rf_step = (rf_max - rf_min) / sim_points

    # rf amplitude range for x-axis
    rf_range = np.arange(rf_min, rf_max, rf_step)

    # Scale Scale RF pulse over the entire RF range and convert from Hz to rad/s
    rf_amp = (
        2 * np.pi * rf_range.reshape(-1, 1) * rf_amp.reshape(1, -1)
    )  # can also use np.outer

    # Memory allocation
    Mx = np.zeros(sim_points)
    My = np.zeros(sim_points)
    Mz = np.zeros(sim_points)

    Rx = np.identity(3)
    Ry = np.identity(3)
    Ry2 = np.identity(3)
    Rz = np.identity(3)
    Rz2 = np.identity(3)

    M = np.zeros([3, 1])

    # Start Simulation
    for rf_range_counter in range(sim_points):
        M[0, 0] = Mx_0
        M[1, 0] = My_0
        M[2, 0] = Mz_0

        for rf_pulse_counter in range(rf.size):
            term_0 = rf_amp[rf_range_counter, rf_pulse_counter] ** 2
            term_1 = rf_offset ** 2

            # B_effective
            Be = np.sqrt(term_0 + term_1) * rf_dwelltime
            alpha = np.arctan2(
                rf_offset, rf_amp[rf_range_counter, rf_pulse_counter]
            )

            # Precalculate various sin/cos terms for increased speed
            cosBe = np.cos(Be)
            sinBe = np.sin(Be)
            cosalpha = np.cos(alpha)
            sinalpha = np.sin(alpha)
            cosphi = np.cos(rf_phase[rf_pulse_counter])
            sinphi = np.sin(rf_phase[rf_pulse_counter])

            # Construct the total rotation matrix
            Rx[1, 1] = cosBe
            Rx[1, 2] = sinBe
            Rx[2, 1] = -1.0 * sinBe
            Rx[2, 2] = cosBe

            Ry[0, 0] = cosalpha
            Ry[0, 2] = -1.0 * sinalpha
            Ry[2, 0] = sinalpha
            Ry[2, 2] = cosalpha

            Ry2[0, 0] = cosalpha
            Ry2[0, 2] = sinalpha
            Ry2[2, 0] = -1.0 * sinalpha
            Ry2[2, 2] = cosalpha

            Rz[0, 0] = cosphi
            Rz[0, 1] = sinphi
            Rz[1, 0] = -1.0 * sinphi
            Rz[1, 1] = cosphi

            Rz2[0, 0] = cosphi
            Rz2[0, 1] = -1.0 * sinphi
            Rz2[1, 0] = sinphi
            Rz2[1, 1] = cosphi

            M = np.linalg.multi_dot([Rz, Ry, Rx, Ry2, Rz2, M])

        Mx[rf_range_counter] = M[0, 0]
        My[rf_range_counter] = M[1, 0]
        Mz[rf_range_counter] = M[2, 0]

    Mxy = np.sqrt(Mx ** 2 + My ** 2)

    return (Mx, My, Mz, Mxy, rf_range, rf_amp)


def simulate_freq(
    pulse_length, sim_points, rf, pm, freq_min, freq_max, rf_amplitude, init_mag
):
    """
    A function to simulate magnetization as a function 
    of Frequency Offset. In the current version, the pulse
    needs to be pre-loaded outside the function.
    """

    rf_amp = rf
    rf_dwelltime = 1e-6 * (pulse_length / rf.size)  # convert to seconds

    # Convert rf phase from degrees to radians
    rf_phase = (np.pi / 180) * pm

    # Initial magnetizations
    Mx_0 = init_mag[0]
    My_0 = init_mag[1]
    Mz_0 = init_mag[2]

    # x-axis step size
    freq_step = (freq_max - freq_min) / sim_points

    # rf frequency offset range for x-axis
    freq_range = np.arange(freq_min, freq_max, freq_step)

    # Scale Scale RF pulse over the entire RF range and convert from Hz to rad/s
    rf_amp = 2 * np.pi * rf_amplitude * rf_amp.reshape(-1, 1)

    # Memory allocation
    Mx = np.zeros(sim_points)
    My = np.zeros(sim_points)
    Mz = np.zeros(sim_points)

    Rx = np.identity(3)
    Ry = np.identity(3)
    Ry2 = np.identity(3)
    Rz = np.identity(3)
    Rz2 = np.identity(3)

    M = np.zeros([3, 1])

    # Start Simulation
    for freq_range_counter in range(sim_points):

        M[0, 0] = Mx_0
        M[1, 0] = My_0
        M[2, 0] = Mz_0

        # Convert frequency offset from Hz to rad/s
        rf_offset = 2 * np.pi * freq_range[freq_range_counter]

        for rf_pulse_counter in range(rf.size):
            term_0 = rf_amp[rf_pulse_counter] ** 2
            term_1 = rf_offset ** 2

            # B_effective
            Be = np.sqrt(term_0 + term_1) * rf_dwelltime
            alpha = np.arctan2(rf_offset, rf_amp[rf_pulse_counter])

            # Precalculate various sin/cos terms for increased speed
            cosBe = np.cos(Be)
            sinBe = np.sin(Be)
            cosalpha = np.cos(alpha)
            sinalpha = np.sin(alpha)
            cosphi = np.cos(rf_phase[rf_pulse_counter])
            sinphi = np.sin(rf_phase[rf_pulse_counter])

            # Construct the total rotation matrix
            Rx[1, 1] = cosBe
            Rx[1, 2] = sinBe
            Rx[2, 1] = -1.0 * sinBe
            Rx[2, 2] = cosBe

            Ry[0, 0] = cosalpha
            Ry[0, 2] = -1.0 * sinalpha
            Ry[2, 0] = sinalpha
            Ry[2, 2] = cosalpha

            Ry2[0, 0] = cosalpha
            Ry2[0, 2] = sinalpha
            Ry2[2, 0] = -1.0 * sinalpha
            Ry2[2, 2] = cosalpha

            Rz[0, 0] = cosphi
            Rz[0, 1] = sinphi
            Rz[1, 0] = -1.0 * sinphi
            Rz[1, 1] = cosphi

            Rz2[0, 0] = cosphi
            Rz2[0, 1] = -1.0 * sinphi
            Rz2[1, 0] = sinphi
            Rz2[1, 1] = cosphi

            M = np.linalg.multi_dot([Rz, Ry, Rx, Ry2, Rz2, M])

        Mx[freq_range_counter] = M[0, 0]
        My[freq_range_counter] = M[1, 0]
        Mz[freq_range_counter] = M[2, 0]

    Mxy = np.sqrt(Mx ** 2 + My ** 2)

    return (Mx, My, Mz, Mxy, freq_range)


def simulate_trajectory(
    pulse_length, rf, pm, rf_amplitude, freq_offset, init_mag
):
    """
    A function to simulate magnetization as a function 
    of time and on the Bloch sphere. In the current version, the pulse
    needs to be pre-loaded outside the function.
    """

    rf_amp = rf
    rf_dwelltime = 1e-6 * (pulse_length / rf.size)  # convert to seconds
    time_axis = np.arange(0, pulse_length, pulse_length / rf.size)

    # Convert rf phase from degrees to radians
    rf_phase = (np.pi / 180) * pm

    # Initial magnetizations
    Mx_0 = init_mag[0]
    My_0 = init_mag[1]
    Mz_0 = init_mag[2]

    # Scale Scale RF pulse over the entire RF range and convert from Hz to rad/s
    rf_amp = 2 * np.pi * rf_amplitude * rf_amp.reshape(-1, 1)

    # Memory allocation
    Mx = np.zeros(rf.size)
    My = np.zeros(rf.size)
    Mz = np.zeros(rf.size)

    Rx = np.identity(3)
    Ry = np.identity(3)
    Ry2 = np.identity(3)
    Rz = np.identity(3)
    Rz2 = np.identity(3)

    M = np.zeros([3, 1])

    # Start Simulation

    M[0, 0] = Mx_0
    M[1, 0] = My_0
    M[2, 0] = Mz_0

    # Convert frequency offset from Hz to rad/s
    rf_offset = 2 * np.pi * freq_offset

    for rf_pulse_counter in range(rf.size):

        term_0 = rf_amp[rf_pulse_counter] ** 2
        term_1 = rf_offset ** 2

        # B_effective
        Be = np.sqrt(term_0 + term_1) * rf_dwelltime
        alpha = np.arctan2(rf_offset, rf_amp[rf_pulse_counter])

        # Precalculate various sin/cos terms for increased speed
        cosBe = np.cos(Be)
        sinBe = np.sin(Be)
        cosalpha = np.cos(alpha)
        sinalpha = np.sin(alpha)
        cosphi = np.cos(rf_phase[rf_pulse_counter])
        sinphi = np.sin(rf_phase[rf_pulse_counter])

        # Construct the total rotation matrix
        Rx[1, 1] = cosBe
        Rx[1, 2] = sinBe
        Rx[2, 1] = -1.0 * sinBe
        Rx[2, 2] = cosBe

        Ry[0, 0] = cosalpha
        Ry[0, 2] = -1.0 * sinalpha
        Ry[2, 0] = sinalpha
        Ry[2, 2] = cosalpha

        Ry2[0, 0] = cosalpha
        Ry2[0, 2] = sinalpha
        Ry2[2, 0] = -1.0 * sinalpha
        Ry2[2, 2] = cosalpha

        Rz[0, 0] = cosphi
        Rz[0, 1] = sinphi
        Rz[1, 0] = -1.0 * sinphi
        Rz[1, 1] = cosphi

        Rz2[0, 0] = cosphi
        Rz2[0, 1] = -1.0 * sinphi
        Rz2[1, 0] = sinphi
        Rz2[1, 1] = cosphi

        M = np.linalg.multi_dot([Rz, Ry, Rx, Ry2, Rz2, M])

        Mx[rf_pulse_counter] = M[0, 0]
        My[rf_pulse_counter] = M[1, 0]
        Mz[rf_pulse_counter] = M[2, 0]

    return (Mx, My, Mz, time_axis)

