# Warren-Cowley Short-Range Order (SRO) Parameter Calculator for High-Entropy Alloys

## Overview
This Python script calculates Warren-Cowley Short-Range Order (SRO) parameters for multi-component High-Entropy Alloys (HEAs) using the Ovito Python API. The script provides a comprehensive analysis of local atomic configurations, helping researchers understand the degree of atomic clustering or segregation in complex alloy systems.

## Features
- Calculates Warren-Cowley SRO parameters for all possible atomic pair combinations
- Supports multiple neighbor shells
- Generates a data table with SRO parameters
- Compatible with Ovito data analysis framework
- Flexible cutoff distance configuration

## Dependencies
- Ovito
- NumPy
- Python 3.x

## Installation
1. Install Ovito (https://www.ovito.org/)
2. Ensure NumPy is installed (`pip install numpy`)

## Usage
```python
# Example usage in Ovito Python scripting environment
import sys
sys.path.append('/path/to/script/directory')
from warren_cowley_sro import calculate_sro

# Load your atomic configuration
data = /* Your Ovito data object */
calculate_sro(data, min_cutoff=0.0, max_cutoff=3.0, neighbor_shell=1)
```

## Parameters
- `min_cutoff`: Minimum neighbor distance (default: 0.0)
- `max_cutoff`: Maximum neighbor distance (default: 3.0)
- `neighbor_shell`: Neighbor shell for SRO calculation (default: 1)

## Theory
Warren-Cowley SRO parameters quantify the local atomic ordering in multi-component systems. A value of 0 indicates random distribution, positive values suggest segregation, and negative values indicate clustering.

## Citing
If you use this script in your research, please cite:
- Cowley, J. M. (1950).
- https://doi.org/10.1016/j.surfcoat.2022.129136

## License
MIT

## Contributing
Contributions are welcome! Please submit pull requests or open issues on GitHub.

## Contact
Prashant Dwivedi
