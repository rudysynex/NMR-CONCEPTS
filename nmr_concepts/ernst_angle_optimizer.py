#!/usr/bin/env python3
"""
Plot steady-state NMR signal-per-scan or sensitivity-per-unit-time vs. flip angle for given TR and T1,
and annotate the maximum (Ernst angle).
[Disclaimer: Generative AI was used to generate parts ofthis code.]

Model:
  E = exp(-TR/T1)
  E2 = exp(-TE/T2) if T2 and TE are provided, otherwise E2 = 1
  S(alpha) ∝ [ sin(alpha) * (1 - E) ] / [ 1 - E * cos(alpha) ] * E2

This is the classic Ernst steady-state formula with optional T2 decay. The maximizing flip angle is:
  alpha_E = arccos(E)

Sources: 
https://labs.dgsom.ucla.edu/file/92487/M229_Lecture4_PulseSeqGRE.pdf
https://www.mritoolbox.com/Contrastinator.html
https://aheyam.science/ernst_angle.html
"""

import math
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Global variables for TR and T1
TR = None  # Repetition time (seconds) - single value. If None, use TR_list.
TR_list = [0.1, 0.3, 0.5, 0.8, 1.0]  # List of TRs for multiple curves. If None, use TR.
T1 = 0.35  # Spin-lattice relaxation time (seconds)
T2 = 0.2  # Spin-spin relaxation time (seconds). If None, E2 = 1.
TE = 0.096  # Total Echo time (seconds). If None, E2 = 1.

signal_averaging = False  # If True, plot sensitivity-per-unit-time instead of signal-per-scan

def ernst_signal(alpha_rad: np.ndarray, TR: float, T1: float, T2: float = None, TE: float = None) -> np.ndarray:
    """Return signal-per-scan (arbitrary units) for flip angles (radians)."""
    if TR <= 0 or T1 <= 0:
        raise ValueError("TR and T1 must be positive.")
    E = math.exp(-TR / T1)
    
    # Handle T2 decay: if T2 and TE are provided, calculate E2, otherwise E2 = 1
    if T2 is not None and TE is not None and T2 > 0 and TE >= 0:
        E2 = math.exp(-TE / T2)
    else:
        E2 = 1.0
    
    # Avoid division by zero at alpha=0 when using vectorized operations
    cos_a = np.cos(alpha_rad)
    sin_a = np.sin(alpha_rad)
    if signal_averaging:
        average_factor = np.sqrt(1/TR)
        return average_factor * (sin_a * (1.0 - E)) / (1.0 - E * cos_a) * E2
    else:
        return (sin_a * (1.0 - E)) / (1.0 - E * cos_a) * E2

def ernst_angle(TR: float, T1: float) -> float:
    """Return Ernst flip angle in radians."""
    E = math.exp(-TR / T1)
    # Clamp in case of tiny numerical drift
    E = min(1.0, max(0.0, E))
    return math.acos(E)

def main():
    p = argparse.ArgumentParser(description="Plot signal per scan vs flip angle and mark the Ernst angle.")
    p.add_argument("--TR", type=float, help=f"Repetition time (seconds, default: {TR})")
    p.add_argument("--TR-list", nargs='+', type=float, help=f"List of TRs for multiple curves (default: {TR_list})")
    p.add_argument("--T1", type=float, help=f"Spin-lattice relaxation time (seconds, default: {T1})")
    p.add_argument("--T2", type=float, help=f"Spin-spin relaxation time (seconds, default: {T2})")
    p.add_argument("--TE", type=float, help=f"Total Echo time (seconds, default: {TE})")
    p.add_argument("--degmax", type=float, default=90.0, help="Max flip angle to plot (degrees, default 90)")
    p.add_argument("--outfile", type=str, default=None, help="Optional output PNG filename")
    p.add_argument("--self-test", action="store_true", help="Run a quick sanity check and exit")
    args = p.parse_args()

    if args.self_test:
        # Quick check: maximum of S(alpha) occurs at alpha ≈ arccos(E)
        test_TR, test_T1 = 0.63, 0.50
        E = math.exp(-test_TR/test_T1)
        a_grid = np.linspace(1e-6, math.pi - 1e-6, 20000)
        Sg = ernst_signal(a_grid, test_TR, test_T1)
        a_num = a_grid[np.argmax(Sg)]
        a_theory = math.acos(E)
        assert abs(a_num - a_theory) < math.radians(0.1), f"Max mismatch: {math.degrees(a_num):.3f} vs {math.degrees(a_theory):.3f}"
        print("Self-test passed.")
        return

    # Use command line arguments if provided, otherwise use global variables
    current_T1 = args.T1 if args.T1 is not None else T1
    current_T2 = args.T2 if args.T2 is not None else T2
    current_TE = args.TE if args.TE is not None else TE
    
    # Determine which TR values to use
    if args.TR_list is not None:
        # Use the provided TR list
        current_TR_list = args.TR_list
        single_TR_mode = False
    elif args.TR is not None:
        # Use single TR value
        current_TR_list = [args.TR]
        single_TR_mode = True
    else:
        # Use global TR list
        current_TR_list = TR_list
        single_TR_mode = False
    
    if current_T1 <= 0:
        raise SystemExit("Error: T1 must be positive.")
    
    if any(tr <= 0 for tr in current_TR_list):
        raise SystemExit("Error: All TR values must be positive.")
    
    # Validate T2 and TE if provided
    if current_T2 is not None and current_T2 <= 0:
        raise SystemExit("Error: T2 must be positive if provided.")
    if current_TE is not None and current_TE < 0:
        raise SystemExit("Error: TE must be non-negative if provided.")

    # Angle grid (avoid singular endpoints)
    degmax = max(10.0, min(90.0, args.degmax))
    alpha_deg = np.linspace(1e-3, degmax - 1e-3, 2000)
    alpha_rad = np.deg2rad(alpha_deg)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Colors for different TR curves
    colors = plt.cm.viridis(np.linspace(0, 1, len(current_TR_list)))
    
    for i, current_TR in enumerate(current_TR_list):
        S = ernst_signal(alpha_rad, current_TR, current_T1, current_T2, current_TE)
        
        # Find maximum and Ernst angle for this TR
        idx_max = int(np.argmax(S))
        a_max_deg = float(alpha_deg[idx_max])
        S_max = float(S[idx_max])

        a_ernst_rad = ernst_angle(current_TR, current_T1)
        a_ernst_deg = math.degrees(a_ernst_rad)

        # Plot the curve
        label = f"TR = {current_TR:.3g} s"
        ax.plot(alpha_deg, S, linewidth=2, color=colors[i], label=label)
        
        # Vertical line at Ernst angle
        ax.axvline(a_ernst_deg, linestyle="--", linewidth=1, color=colors[i], alpha=0.7)
        
        # Annotation for this TR
        ax.plot([a_max_deg], [S_max], marker="o", color=colors[i], markersize=6)
        ax.annotate(
            f"TR={current_TR:.3g}s: α≈{a_ernst_deg:.1f}°",
            xy=(a_max_deg, S_max),
            xytext=(a_max_deg + 5, S_max * 0.9),
            arrowprops=dict(arrowstyle="->", lw=1, color=colors[i]),
            va="center",
            fontsize=9,
            color=colors[i]
        )

    ax.set_xlabel("Flip angle α (degrees)")
    if signal_averaging:
        ax.set_ylabel("Sensitivity per unit time (arb. units)")
    else:
        ax.set_ylabel("Signal per scan (arb. units)")
    
    # Build title with all relevant parameters
    title_parts = []
    if single_TR_mode:
        title_parts.append(f"TR={current_TR_list[0]:.3g}s")
    else:
        title_parts.append("multiple TRs")
    
    title_parts.append(f"T1={current_T1:.3g}s")
    
    if current_T2 is not None and current_TE is not None:
        title_parts.append(f"T2={current_T2:.3g}s")
        title_parts.append(f"TE={current_TE:.3g}s")
    
    title = f"Ernst steady-state signal vs flip angle ({', '.join(title_parts)})"
    
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Optional y-limit padding
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + 0.05 * (ymax - ymin))

    plt.tight_layout()

    if args.outfile:
        fig.savefig(args.outfile, dpi=160)
        print(f"Saved: {args.outfile}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
