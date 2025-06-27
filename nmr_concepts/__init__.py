"""
NMR Concepts - A Python package for NMR simulation and analysis.

This package provides tools for simulating NMR experiments using Bloch equations,
including various pulse shapes and composite pulses.
"""

from .bloch_sim import (
    simulate_amp,
    simulate_freq,
    simulate_trajectory
)

__version__ = "1.0.0"
__author__ = "NMR Concepts Team"

__all__ = [
    "simulate_amp",
    "simulate_freq", 
    "simulate_trajectory"
]
