# WBM soil moisture submodule
Python script to reproduce WBM soil moisture submodule. Equations are taken from [Grogan (2016)](https://scholars.unh.edu/dissertation/2/), [Grogan et al. (2022)](https://gmd.copernicus.org/articles/15/7287/2022/), and WBM v 23.1.0 source code.

There are some differences from full WBM:
- Here we simulate a single grid cell so there is no lateral transport of moisture
- Initial groundwater volumes can be set arbitrarily high/low. At the start/end of the growing season there is no transfer of moisture between the active soil layer and groundwater storage. This contrasts with full WBM in which some transfer occurs to ensure continuity in soilMoistFrac.
- Open water evaporation `Eow` is assumed constant here and tuned to a reasonable value based on comparing to full WBM outputs. In full WBM, it is calculated based on a suite of meteorological variables.
- Soil moisture capacity is here a timeseries input. There is no explicit representation of rootDepth or available water capacity. In comparing to full WBM the soil moisture capacity is calculated to align with full WBM outputs (i.e. computed via soilMoist/soilMoistFrac)

Agreement with full WBM is generally quite good. This can be explored in the `wbm_toy_model_eval` notebook where we compare to full WBM outputs over three central US locations (full WBM is run using Hamon PET and crop averaging). The largest differences seem to emerge at the end of the calendar year, after the growing season. 
