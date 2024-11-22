"""
Warren-Cowley Short-Range Order (SRO) Parameter Calculator

This module provides functionality to calculate Warren-Cowley SRO parameters
for multi-component High-Entropy Alloys using the Ovito Python API.

Author: [Your Name]
Date: 2024
Version: 1.0
"""

from ovito.data import *
import numpy as np
import itertools
from typing import Optional

def calculate_sro(
    data: DataObject, 
    min_cutoff: float = 0.0, 
    max_cutoff: float = 3.0, 
    neighbor_shell: int = 1
) -> Optional[DataTable]:
    """
    Calculate Warren-Cowley Short-Range Order parameters for all atomic pair combinations.

    Args:
        data (DataObject): Ovito data object containing atomic configuration
        min_cutoff (float): Minimum neighbor distance
        max_cutoff (float): Maximum neighbor distance
        neighbor_shell (int): Neighbor shell for SRO calculation

    Returns:
        DataTable: Table containing SRO parameters for each atomic pair
    """
    # Validate input data
    if not data or not hasattr(data, 'particles'):
        raise ValueError("Invalid data object or missing particle information")

    # Prefetch particle type information
    ptypes = data.particles.particle_types
    pair_combinations = list(itertools.product(set(ptypes), set(ptypes)))
    
    # Initialize neighbor finder object
    finder = CutoffNeighborFinder(max_cutoff, data)
    
    # Create result table for SRO parameters
    table = DataTable(
        title='Warren Cowley SRO Parameters', 
        plot_mode=DataTable.PlotMode.BarChart
    )
    table.x = table.create_property('WC-SRO parameters', data=np.arange(len(pair_combinations)))
    
    y_values = []
    
    # Compute SRO for each pair combination
    for count, (type_i, type_j) in enumerate(pair_combinations):
        # Generate property name
        prop_name = _generate_property_name(
            neighbor_shell, 
            ptypes.type_by_id(type_i), 
            ptypes.type_by_id(type_j)
        )
        
        # Add type to x-axis
        table.x.types.append(ElementType(id=count, name=prop_name))
        
        # Create property to store SRO values
        prop = data.particles_.create_property(prop_name, data=np.zeros(data.particles.count, dtype=float))
        prop[:] = np.nan
        
        # Compute concentration of particle type j
        c_j = np.count_nonzero(ptypes == type_j) / data.particles.count
        
        # Compute SRO for each particle of type i
        sro_values = _compute_sro_values(finder, ptypes, type_i, type_j, min_cutoff, c_j)
        prop[:] = sro_values
        
        # Store average SRO value
        avg_sro = np.nanmean(sro_values)
        data.attributes[f"<{prop_name}>"] = avg_sro
        y_values.append(avg_sro)
    
    # Finalize table
    table.y = table.create_property('Value', data=y_values)
    data.objects.append(table)
    
    return table

def _generate_property_name(neighbor_shell: int, type_i: ElementType, type_j: ElementType) -> str:
    """
    Generate a descriptive property name for SRO parameter.
    
    Args:
        neighbor_shell (int): Neighbor shell number
        type_i (ElementType): First atomic type
        type_j (ElementType): Second atomic type
    
    Returns:
        str: Formatted property name
    """
    name_i = type_i.name if type_i.name else str(type_i.id)
    name_j = type_j.name if type_j.name else str(type_j.id)
    return f"a_{neighbor_shell}_{name_i}-{name_j}"

def _compute_sro_values(
    finder: CutoffNeighborFinder, 
    ptypes: PropertyArray, 
    type_i: int, 
    type_j: int, 
    min_cutoff: float, 
    c_j: float
) -> np.ndarray:
    """
    Compute SRO values for particles of type i.
    
    Args:
        finder (CutoffNeighborFinder): Neighbor finder object
        ptypes (PropertyArray): Particle type array
        type_i (int): First atomic type
        type_j (int): Second atomic type
        min_cutoff (float): Minimum neighbor distance
        c_j (float): Concentration of particle type j
    
    Returns:
        np.ndarray: Array of SRO values
    """
    sro_values = np.full(len(ptypes), np.nan)
    
    # Compute SRO for particles of type i
    for index in np.where(ptypes == type_i)[0]:
        neighbors = list(finder.find(index))
        valid_neighbors = [n for n in neighbors if n.distance >= min_cutoff]
        
        if not valid_neighbors:
            continue
        
        num_neighbors = len(valid_neighbors)
        num_type_j_neighbors = sum(1 for n in valid_neighbors if ptypes[n.index] == type_j)
        
        sro_values[index] = 1.0 - num_type_j_neighbors / (num_neighbors * c_j)
    
    return sro_values

# Example usage
if __name__ == "__main__":
    print("This is a module for SRO parameter calculation in HEA.")
    print("Import and use calculate_sro() function in your Ovito script.")
