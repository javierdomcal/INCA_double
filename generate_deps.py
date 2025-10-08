import os
import re
from pathlib import Path

SRC_DIR = "src"
OBJ_DIR = "obj"
DEPS_FILE = "deps.mk"

print(f"Working directory: {os.getcwd()}")
print(f"Looking for .f90 files in: {SRC_DIR}")

# Regex to match module usage
use_pattern = re.compile(r'^\s*use\s+(\w+)', re.IGNORECASE)

# System/intrinsic modules to ignore
SYSTEM_MODULES = {
    'omp_lib', 'iso_c_binding', 'iso_fortran_env', 
    'ieee_arithmetic', 'ieee_exceptions', 'ieee_features'
}

# Find all module definitions in the source directory
local_modules = set()
f90_files = list(Path(SRC_DIR).glob("*.f90"))
print(f"Found {len(f90_files)} .f90 files")

for f90_file in f90_files:
    try:
        with f90_file.open("r", encoding='utf-8', errors='ignore') as f:
            for line in f:
                # Match module definition (not module procedure)
                mod_match = re.match(r'^\s*module\s+(\w+)', line, re.IGNORECASE)
                if mod_match:
                    mod_name = mod_match.group(1).lower()
                    if mod_name != 'procedure':  # Exclude "module procedure"
                        local_modules.add(mod_name)
    except Exception as e:
        print(f"Error reading {f90_file}: {e}")

print(f"Found {len(local_modules)} local modules: {sorted(local_modules)}")

# Map to hold dependencies
deps = {}

for f90_file in f90_files:
    base = f90_file.stem
    obj_file = f"{OBJ_DIR}/{base}.o"
    deps[obj_file] = [str(f90_file)]
    
    try:
        with f90_file.open("r", encoding='utf-8', errors='ignore') as f:
            for line in f:
                match = use_pattern.match(line)
                if match:
                    mod_name = match.group(1).lower()
                    
                    # Skip system modules
                    if mod_name in SYSTEM_MODULES:
                        continue
                    
                    # Only add dependency if it's a local module
                    if mod_name in local_modules:
                        mod_file = f"{OBJ_DIR}/mod_{mod_name}.o"
                        if mod_file != obj_file and mod_file not in deps[obj_file]:
                            deps[obj_file].append(mod_file)
    except Exception as e:
        print(f"Error processing {f90_file}: {e}")

print(f"Generated dependencies for {len(deps)} object files")

# Write deps.mk
try:
    with open(DEPS_FILE, "w") as out:
        for obj, dep_list in sorted(deps.items()):
            out.write(f"{obj} : {' '.join(dep_list)}\n")
    print(f"Successfully wrote {DEPS_FILE}")
except Exception as e:
    print(f"Error writing {DEPS_FILE}: {e}")
