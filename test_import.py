#!/usr/bin/env python3
"""
Simple test script to verify nmr_concepts package import.
"""

import sys
import os

print("Testing nmr_concepts package import...")
print("=" * 40)

# Method 1: Try direct import
print("Method 1: Direct import")
try:
    from nmr_concepts.bloch_sim import simulate_freq
    print("✓ Success: Direct import worked")
except ImportError as e:
    print(f"✗ Failed: {e}")

# Method 2: Add current directory to path
print("\nMethod 2: Add current directory to path")
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from nmr_concepts.bloch_sim import simulate_freq
    print("✓ Success: Import with path modification worked")
except ImportError as e:
    print(f"✗ Failed: {e}")

# Method 3: Import from nmr_concepts directory
print("\nMethod 3: Import from nmr_concepts directory")
try:
    nmr_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nmr_concepts')
    sys.path.insert(0, nmr_path)
    from bloch_sim import simulate_freq
    print("✓ Success: Direct module import worked")
except ImportError as e:
    print(f"✗ Failed: {e}")

# Method 4: Test package import
print("\nMethod 4: Test package import")
try:
    import nmr_concepts
    print(f"✓ Success: Package imported, version: {nmr_concepts.__version__}")
    print(f"Available functions: {nmr_concepts.__all__}")
except ImportError as e:
    print(f"✗ Failed: {e}")

print("\n" + "=" * 40)
print("Import test completed!") 