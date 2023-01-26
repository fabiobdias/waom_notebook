# waom_notebook
# Fabio Boeira Dias (CCRC-UNSW, previously at University of Helsinki/Finland)

Python scripts to pre- and (mostly) postprocessing of the Whole Antarctica Ocean Model (WAOM)

A bit of a background:

WAOM was developed by Ole Richter (now at AWI) and Ben Galton-Fenzi at the University of Tasmania (UTAS). The v1.0 was published (submitted) in 2022 (2021) at GMDD (https://doi.org/10.5194/gmd-15-617-2022). The model setup was further developed at the University of Helsinki/Finland by me and Petteri Uotila (with support of Ben and Ole) where modifications on the (1) surface fluxes treatment, (2) model domain, (3) surface temperature restoring. These setup changes were done aiming to reduce surface temperature bias, to increase salt input from the Tamura et al (2011, https://doi.org/10.2151/sola.2011-005) and to better resolve the Ross Gyre and circulation features on the surrounding regions.

The waom_notebook package was developed to analyse the model runs with 3 different horizontal resolution (WAOM10, WAOM4 and WAOM2 with 10, 4 and 2 km, respectively) and sensitivity experiments (WAOM4 with coarse bathymetry, WAOM4-COARSE; WAOM4 without tides, WAOM4-NOTIDES). These analyses were used to investigate the sensitivity of WMT on the Antarctica continental shelf (including ice shelf cavities). Paper currently under review (https://review.frontiersin.org/review/1027704/16/1802900/#tab/History).
