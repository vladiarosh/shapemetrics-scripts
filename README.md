# shapemetrics-scripts
This repository is an addition to the original [Shapemetrics script](https://github.com/kerosuolab/shapemetrics) by the Kerosuo Lab.  
It includes several modifications and additional scripts to expand analysis capabilities using the parameter matrices derived using the original script.

## Contents

### I. `Shapemetrics_3D_corrected_and_stalk_removal.m`

This file includes a few small changes to the original 3D reconstruction script:

- Corrected the elongation formula  
- Added stalk-removal functionality
- Added variables for lower and upper bounds of cell volume in cubic microns (previously these cell size variables were in voxels only, which was unhandy to interpet during analysis)

### II. `Mapping_cells_back_subset_of_parameters.m`

A slightly modified version of a similar section in the original script.  
Instead of generating a joint heatmap using all parameters from the parameter matrix, this version allows selection of a desired subset of parameters.

### III. `Mapping_cells_back_threshold_based.m`

This script filters the z-scored parameter matrix using a given threshold value across all parameters.

### IV. `Mapping_cells_back_using_specific_range.m`

This file includes several scripts:

1. **Mapping cells back above the threshold**  
   Similar to the script in section III, but uses raw parameter values. Filtering is done prior to z-scoring.

2. **Visualize single parameters within a specific range**  
   Visualizes cells that fall within a specified range of a selected parameter (for example, cells between 200 and 250 cubic microns).

3. **Visualize a single parameter across a set of ranges**  
   Same functionality as above, but allows input of multiple ranges.

### V. `Top_and_bottom_tails.m`

This file includes scripts to map back cells that fall above the 90th percentile or below the 10th percentile for a selected parameter.  
It includes scripts for both individual sample testing and for batch processing across all sample folders.
