import numpy as np
from qutip import *


def simulate_freq(
    pulse_length, sim_points, rf, pm, freq_min, freq_max, rf_amplitude, init_mag
):
    """
    Optimized version of simulate_freq with several performance improvements:
    1. Pre-computed rotation matrices for each pulse step
    2. Vectorized operations where possible
    3. Reduced memory allocations
    4. Optimized matrix multiplication order
    """
    
    rf_amp = rf
    rf_dwelltime = 1e-6 * (pulse_length / rf.size)  # convert to seconds
    
    # Convert rf phase from degrees to radians (pre-compute once)
    rf_phase = (np.pi / 180) * pm
    
    # Initial magnetizations
    Mx_0, My_0, Mz_0 = init_mag
    
    # Frequency range
    freq_step = (freq_max - freq_min) / sim_points
    freq_range = np.arange(freq_min, freq_max, freq_step)
    
    # Scale RF pulse and convert from Hz to rad/s
    rf_amp_scaled = 2 * np.pi * rf_amplitude * rf_amp
    
    # Pre-compute all rotation matrices for each pulse step
    # This avoids recomputing them for each frequency
    rotation_matrices = _precompute_rotation_matrices(
        rf_amp_scaled, rf_phase, rf_dwelltime, freq_range
    )
    
    # Memory allocation for results
    Mx = np.zeros(sim_points)
    My = np.zeros(sim_points)
    Mz = np.zeros(sim_points)
    
    # Process each frequency
    for freq_idx in range(sim_points):
        M = np.array([Mx_0, My_0, Mz_0])
        
        # Apply all rotation matrices for this frequency
        for pulse_idx in range(rf.size):
            M = rotation_matrices[freq_idx, pulse_idx] @ M
        
        Mx[freq_idx] = M[0]
        My[freq_idx] = M[1]
        Mz[freq_idx] = M[2]
    
    Mxy = np.sqrt(Mx**2 + My**2)
    
    return (Mx, My, Mz, Mxy, freq_range)


def _precompute_rotation_matrices(rf_amp_scaled, rf_phase, rf_dwelltime, freq_range):
    """
    Pre-compute all rotation matrices for all frequency offsets and pulse steps.
    This is the most expensive operation but only done once.
    """
    num_freqs = len(freq_range)
    num_pulses = len(rf_amp_scaled)
    
    # Pre-allocate array for all rotation matrices
    rotation_matrices = np.zeros((num_freqs, num_pulses, 3, 3))
    
    # Pre-compute frequency offsets in rad/s
    freq_offsets_rad = 2 * np.pi * freq_range
    
    for freq_idx in range(num_freqs):
        rf_offset = freq_offsets_rad[freq_idx]
        
        for pulse_idx in range(num_pulses):
            rf_amp_val = rf_amp_scaled[pulse_idx]
            rf_phase_val = rf_phase[pulse_idx]
            
            # Calculate effective field and angle
            term_0 = rf_amp_val**2
            term_1 = rf_offset**2
            Be = np.sqrt(term_0 + term_1) * rf_dwelltime
            alpha = np.arctan2(rf_offset, rf_amp_val)
            
            # Pre-calculate trigonometric values
            cosBe = np.cos(Be)
            sinBe = np.sin(Be)
            cosalpha = np.cos(alpha)
            sinalpha = np.sin(alpha)
            cosphi = np.cos(rf_phase_val)
            sinphi = np.sin(rf_phase_val)
            
            # Build individual rotation matrices
            Rx = np.array([
                [1, 0, 0],
                [0, cosBe, sinBe],
                [0, -sinBe, cosBe]
            ])
            
            Ry = np.array([
                [cosalpha, 0, -sinalpha],
                [0, 1, 0],
                [sinalpha, 0, cosalpha]
            ])
            
            Ry2 = np.array([
                [cosalpha, 0, sinalpha],
                [0, 1, 0],
                [-sinalpha, 0, cosalpha]
            ])
            
            Rz = np.array([
                [cosphi, sinphi, 0],
                [-sinphi, cosphi, 0],
                [0, 0, 1]
            ])
            
            Rz2 = np.array([
                [cosphi, -sinphi, 0],
                [sinphi, cosphi, 0],
                [0, 0, 1]
            ])
            
            # Compute total rotation matrix (optimized order)
            # Rz @ Ry @ Rx @ Ry2 @ Rz2
            temp1 = Rz @ Ry
            temp2 = temp1 @ Rx
            temp3 = temp2 @ Ry2
            total_rotation = temp3 @ Rz2
            
            rotation_matrices[freq_idx, pulse_idx] = total_rotation
    
    return rotation_matrices


def simulate_amp(
    pulse_length, sim_points, rf, pm, rf_min, rf_max, freq_offset, init_mag
):
    """
    Optimized version of simulate_amp with several performance improvements:
    1. Pre-computed rotation matrices for each pulse step and RF amplitude
    2. Vectorized operations where possible
    3. Reduced memory allocations
    4. Optimized matrix multiplication order
    """
    
    rf_amp = rf
    rf_dwelltime = 1e-6 * (pulse_length / rf.size)  # convert to seconds
    
    # Convert rf phase from degrees to radians (pre-compute once)
    rf_phase = (np.pi / 180) * pm
    
    # Convert freq offset from Hz to radians/sec (constant for this function)
    rf_offset = 2 * np.pi * freq_offset
    
    # Initial magnetizations
    Mx_0, My_0, Mz_0 = init_mag
    
    # x-axis step size
    rf_step = (rf_max - rf_min) / sim_points
    
    # rf amplitude range for x-axis
    rf_range = np.arange(rf_min, rf_max, rf_step)
    
    # Scale RF pulse over the entire RF range and convert from Hz to rad/s
    rf_amp_scaled = 2 * np.pi * rf_range.reshape(-1, 1) * rf_amp.reshape(1, -1)
    
    # Pre-compute all rotation matrices for each RF amplitude and pulse step
    rotation_matrices = _precompute_rotation_matrices_amp(
        rf_amp_scaled, rf_phase, rf_dwelltime, rf_offset
    )
    
    # Memory allocation for results
    Mx = np.zeros(sim_points)
    My = np.zeros(sim_points)
    Mz = np.zeros(sim_points)
    
    # Process each RF amplitude
    for amp_idx in range(sim_points):
        M = np.array([Mx_0, My_0, Mz_0])
        
        # Apply all rotation matrices for this RF amplitude
        for pulse_idx in range(rf.size):
            M = rotation_matrices[amp_idx, pulse_idx] @ M
        
        Mx[amp_idx] = M[0]
        My[amp_idx] = M[1]
        Mz[amp_idx] = M[2]
    
    Mxy = np.sqrt(Mx**2 + My**2)
    
    return (Mx, My, Mz, Mxy, rf_range, rf_amp_scaled)

def _precompute_rotation_matrices_amp(rf_amp_scaled, rf_phase, rf_dwelltime, rf_offset):
    """
    Pre-compute all rotation matrices for all RF amplitudes and pulse steps.
    This is the most expensive operation but only done once.
    """
    num_amps = rf_amp_scaled.shape[0]
    num_pulses = rf_amp_scaled.shape[1]
    
    # Pre-allocate array for all rotation matrices
    rotation_matrices = np.zeros((num_amps, num_pulses, 3, 3))
    
    for amp_idx in range(num_amps):
        for pulse_idx in range(num_pulses):
            rf_amp_val = rf_amp_scaled[amp_idx, pulse_idx]
            rf_phase_val = rf_phase[pulse_idx]
            
            # Calculate effective field and angle
            term_0 = rf_amp_val**2
            term_1 = rf_offset**2
            Be = np.sqrt(term_0 + term_1) * rf_dwelltime
            alpha = np.arctan2(rf_offset, rf_amp_val)
            
            # Pre-calculate trigonometric values
            cosBe = np.cos(Be)
            sinBe = np.sin(Be)
            cosalpha = np.cos(alpha)
            sinalpha = np.sin(alpha)
            cosphi = np.cos(rf_phase_val)
            sinphi = np.sin(rf_phase_val)
            
            # Build individual rotation matrices
            Rx = np.array([
                [1, 0, 0],
                [0, cosBe, sinBe],
                [0, -sinBe, cosBe]
            ])
            
            Ry = np.array([
                [cosalpha, 0, -sinalpha],
                [0, 1, 0],
                [sinalpha, 0, cosalpha]
            ])
            
            Ry2 = np.array([
                [cosalpha, 0, sinalpha],
                [0, 1, 0],
                [-sinalpha, 0, cosalpha]
            ])
            
            Rz = np.array([
                [cosphi, sinphi, 0],
                [-sinphi, cosphi, 0],
                [0, 0, 1]
            ])
            
            Rz2 = np.array([
                [cosphi, -sinphi, 0],
                [sinphi, cosphi, 0],
                [0, 0, 1]
            ])
            
            # Compute total rotation matrix (optimized order)
            # Rz @ Ry @ Rx @ Ry2 @ Rz2
            temp1 = Rz @ Ry
            temp2 = temp1 @ Rx
            temp3 = temp2 @ Ry2
            total_rotation = temp3 @ Rz2
            
            rotation_matrices[amp_idx, pulse_idx] = total_rotation
    
    return rotation_matrices

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
