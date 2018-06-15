# Getting The Model

In addition to cloning this repository, it is necessary to clone the submodules"

 ``git submodule update --init``
 
The source bathymetry is uploaded to a THREDDS server at SFEI.  The script
download_bathy.py will copy that inputs-static.

# Files

derived
  Data/files derived from other inputs, but not part of a specific run.
  Currently only used for a shapefile generated from the grid.
  
inputs-static
  Files which do not change across runs, are not derived, but are instead
  referenced directly from the MDU or other parts of the model setup.
  
nudged_features.pli
  A Delft-style polyline file defining locations of sources.  This is created
  by exporting boundary condition features from Delta Shell.
  
write_grid_shp.py
  Short script which writes the shapefile for grid edges, used for loading
  grid representation into GIS.
  
runs
  Script-generated simulation setups are in subdirectories below herea
  
sfbay_freshwater
sfbay_potw
  Git submodules holding forcing data for rivers and wastewater discharges
  
lsb_dfm.py
  Main script for generating new runs.
  
lsb_v99_net.nc
  Grid, without bathymetry, in non-ugrid DFM format

lsb_v99_bathy_net.nc
  Grid, with bathymetry, in non-ugrid DFM format
  
template.mdu
  Template model definition.  Settings which needn't be set dynamically can
  be set here.  lsb_dfm.py reads this in and modifies a small subset of
  parameters dynamically.
  
# Process

## Bathymetry


# TODO

``plotting/plot_salt.py`` is old.  Move salt_plots.py into plotting, document.

Are there files in static which can be removed?
