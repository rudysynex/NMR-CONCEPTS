# Fixing the "ModuleNotFoundError: No module named 'nmr_concepts'" Error

## Problem

You're getting this error because Python can't find the `nmr_concepts` module. This happens when the package isn't properly installed or the Python path isn't set correctly.

## Solutions (Choose One)

### Solution 1: Install in Development Mode (Recommended)

This is the cleanest solution that makes the package available system-wide:

```bash
# Navigate to the project root directory (where setup.py is located)
cd /path/to/NMR-CONCEPTS

# Install the package in development mode
pip install -e .
```

**What this does:**
- Installs the package in "editable" mode
- Creates a link to your source code
- Makes `nmr_concepts` importable from anywhere
- Changes to your code are immediately available

### Solution 2: Add to PYTHONPATH

Add the project directory to your Python path:

```bash
# For Linux/Mac
export PYTHONPATH=$PYTHONPATH:$(pwd)

# For Windows PowerShell
$env:PYTHONPATH = "$env:PYTHONPATH;$(Get-Location)"

# For Windows Command Prompt
set PYTHONPATH=%PYTHONPATH%;%CD%
```

### Solution 3: Run from Project Root

Run your scripts from the project root directory:

```bash
# Navigate to project root
cd /path/to/NMR-CONCEPTS

# Run your script
python your_script.py
```

### Solution 4: Modify Scripts to Handle Import

Update your scripts to handle the import dynamically:

```python
import sys
import os

# Add the project root to Python path
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, project_root)

# Now import should work
from nmr_concepts.bloch_sim import simulate_freq
```

## Testing the Fix

Run the test script to verify the import works:

```bash
python test_import.py
```

You should see output like:
```
Testing nmr_concepts package import...
========================================
Method 1: Direct import
✓ Success: Direct import worked
...
```

## Using the Package

Once the import issue is fixed, you can use the package in several ways:

### Method 1: Import specific functions

```python
from nmr_concepts.bloch_sim import simulate_freq

# Use the functions
result = simulate_freq(pulse_length, sim_points, rf, pm, freq_min, freq_max, rf_amplitude, init_mag)
```

### Method 2: Import the whole package

```python
import nmr_concepts

# Use functions with package prefix
result = nmr_concepts.simulate_freq(pulse_length, sim_points, rf, pm, freq_min, freq_max, rf_amplitude, init_mag)
```

### Method 3: Import with alias

```python
import nmr_concepts as nmr

# Use with alias
result = nmr.simulate_freq(pulse_length, sim_points, rf, pm, freq_min, freq_max, rf_amplitude, init_mag)
```

## Running the Package

After fixing the import, you can use the package:

```bash
# Install dependencies first
pip install -r requirements.txt

# Run your own scripts
python your_script.py
```

## Common Issues and Solutions

### Issue: "No module named 'numpy'"
**Solution:** Install dependencies
```bash
pip install -r requirements.txt
```

### Issue: "Permission denied" when installing
**Solution:** Use user installation
```bash
pip install -e . --user
```

### Issue: Import works in terminal but not in IDE
**Solution:** 
1. Restart your IDE
2. Check IDE's Python interpreter settings
3. Make sure IDE is using the same Python environment

## Package Structure

After installation, your package structure should look like:

```
NMR-CONCEPTS/
├── nmr_concepts/
│   ├── __init__.py          # Package initialization
│   ├── bloch_sim.py         # Main simulation functions
│   ├── am_pulses.py         # AM pulse functions
│   ├── burp_pulses.py       # BURP pulse functions
│   └── ...
├── setup.py                 # Package setup
├── requirements.txt         # Dependencies
└── ...
```

## Verification

To verify everything is working:

```python
# Test basic import
import nmr_concepts
print(f"Package version: {nmr_concepts.__version__}")

# Test function import
from nmr_concepts.bloch_sim import simulate_freq
print("✓ simulate_freq imported successfully")

# Test function call (with dummy data)
import numpy as np
try:
    # Create dummy data
    pulse_length = 1000
    sim_points = 10
    rf = np.random.rand(10)
    pm = np.random.rand(10) * 360
    freq_min = -1000
    freq_max = 1000
    rf_amplitude = 100
    init_mag = [0, 0, 1]
    
    # Test function call
    result = simulate_freq(pulse_length, sim_points, rf, pm, freq_min, freq_max, rf_amplitude, init_mag)
    print("✓ simulate_freq executed successfully")
except Exception as e:
    print(f"✗ Function execution failed: {e}")
```

## Next Steps

Once the import is working:

1. **Test with your data** to verify accuracy
2. **Use the optimized functions** for better performance:
   - `simulate_freq()` - Optimized frequency sweep simulation
   - `simulate_amp()` - Optimized amplitude sweep simulation
   - `simulate_trajectory()` - Time trajectory simulation

## Support

If you're still having issues:

1. Check your Python version: `python --version`
2. Check your pip version: `pip --version`
3. Verify you're in the correct directory
4. Try creating a new virtual environment
5. Check the error messages for specific details 