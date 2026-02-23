# File Formats

## Input Files

### JSON Configuration File

The main configuration file in JSON format specifies simulation parameters.

**Required Fields:**
- `sample`: String identifier for the simulation
- `mcs`: Integer number of Monte Carlo steps
- `kb`: Double Boltzmann constant (typically 1.0)
- `temperature`: Object or double specifying temperature range
- `field`: Object or double specifying magnetic field range
- `out`: String output filename (HDF5 format)

**Temperature Specification:**
```json
"temperature": 2.5  // Single temperature value
```
or
```json
"temperature": {
  "start": 0.1,
  "final": 5.0,
  "points": 50,
  "delta": 0.1,  // Alternative to points
  "cycle": false  // If true, cycle back to start
}
```

**Field Specification:**
Same format as temperature.

### Lattice File Format

The lattice file defines atoms and their interactions. Each line in the interactions section specifies an exchange coupling constant J.

**Exchange Convention:**
The exchange energy is calculated as:

```
E = -J * Σ (S_i · S_j)
```

Where:
- **J > 0**: Ferromagnetic coupling (spins align parallel to minimize energy)
- **J < 0**: Antiferromagnetic coupling (spins align antiparallel to minimize energy)

This is the standard physics convention used in statistical mechanics.

**File Structure:**
```
num_atoms num_interactions num_types
type1 type2 ...  # Atom type names
atom_index pos_x pos_y pos_z spin_norm hx hy hz type spin_model
...
interaction_index neighbor_index exchange_J
...
```

### Initial State File

Optional text file specifying initial spin configurations:

```
x1 y1 z1
x2 y2 z2
x3 y3 z3
...
```

Each line contains three floating-point numbers representing the x, y, z components of a spin vector.

### Anisotropy File

Text file specifying anisotropy terms:

**Uniaxial Anisotropy:**
```
ax ay az kan
```
Where (ax, ay, az) is the anisotropy axis unit vector and kan is the anisotropy constant.

**Cubic Anisotropy:**
```
Ax Ay Az Bx By Bz kan
```
Where A and B are orthogonal vectors defining the cubic axes.

## Output Files

### HDF5 Output

The main output file contains:

**Attributes:**
- `mcs`: Number of measurement steps (NOT total MCS - thermalization discarded in v2.3.0+)
- `seed`: Random number generator seed
- `kb`: Boltzmann constant

**Datasets:**
- `temperature`: 1D array of temperature values [n_temps]
- `field`: 1D array of field values [n_temps]
- `magnetization_x`: 2D array [n_temps, measurement_steps] (v2.3.0+)
- `magnetization_y`: 2D array [n_temps, measurement_steps] (v2.3.0+)
- `magnetization_z`: 2D array [n_temps, measurement_steps] (v2.3.0+)
- `energy`: 2D array [n_temps, measurement_steps] (v2.3.0+)
- `positions`: 2D array [n_atoms, 3] of atom positions
- `types`: 1D array [n_atoms] of atom type indices
- `finalstates`: 3D array [n_temps, n_atoms, 3] final spin configurations

**Note**: Since v2.3.0, only measurement phase (80% of MCS) is stored. Thermalization data is discarded.

### Analysis Outputs

Python analyzers produce:

1. **Mean Values File** (`*.mean`):
   ```
   # seed = 12345
   #	T	H	E	Cv	M	Mz	X
   0.1	0.0	0.0	0.0	0.0	0.95	0.0
   0.2	0.0	0.0	0.0	0.0	0.90	0.0
   ...
   ```

2. **PDF Plots**: Graphical representations of magnetization vs field/temperature

3. **XYZ Trajectories**: For visualization in molecular viewers

## Python Analysis Scripts

### vegas-analyzer-lite.py
Basic analysis producing mean magnetization and simple plots.

**Usage:**
```bash
python analyzers/vegas-analyzer-lite.py output.h5
```

**Outputs:**
- `output.mean`: Mean values text file
- `output.pdf`: Magnetization plot

### vegas-analyzer-heisenberg.py
Advanced analysis for Heisenberg models with 3D visualization.

**Usage:**
```bash
python analyzers/vegas-analyzer-heisenberg.py output.h5
```

**Outputs:**
- Multiple PDF files with different visualizations
- 3D spin configurations

### vegas-analyzer-xyz.py
Trajectory analysis for visualization in VMD, PyMOL, etc.

**Usage:**
```bash
python analyzers/vegas-analyzer-xyz.py output.h5
```

**Outputs:**
- `trajectory.xyz`: XYZ format trajectory file

## Data Processing Pipeline

1. **Simulation**: `vegas input.json` → `output.h5`
2. **Basic Analysis**: `vegas-analyzer-lite.py output.h5` → `output.mean`, `output.pdf`
3. **Advanced Analysis**: `vegas-analyzer-heisenberg.py output.h5` → detailed plots
4. **Visualization**: `vegas-analyzer-xyz.py output.h5` → `trajectory.xyz`

## Example Workflow

```bash
# 1. Run simulation
./vegas example_config.json

# 2. Basic analysis
python analyzers/vegas-analyzer-lite.py ferromagnetic_results.h5

# 3. View results
evince ferromagnetic_results.pdf
cat ferromagnetic_results.mean
```

## Troubleshooting

### Common Issues

1. **File Not Found**: Ensure all referenced files exist in the working directory
2. **Permission Denied**: Check file permissions for reading/writing
3. **Invalid Format**: Verify JSON syntax and numerical formats
4. **Memory Issues**: Reduce system size or increase swap space
5. **HDF5 Errors**: Check HDF5 library installation and version

### Error Messages

- `"A JSON file is necessary !!!"`: No input file specified
- `"The initial state file can't open or doesn't exist !!!"`: Missing initial state file
- `"Invalid number format in ..."`: Malformed numerical data
- `"Filename contains path traversal attempt"`: Security violation in filename
- `"Cannot open file: ..."`: File access error