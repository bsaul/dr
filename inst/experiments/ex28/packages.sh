#!/bin/bash

scp -r /Users/bradley/Dropbox/Research/projects/active/causal/manuscripts/upstream_downstream/programs-simulation/sim_packages/* saulb@killdevil.unc.edu:/netscr/saulb/Rpackages
scp -r /Users/bradley/Dropbox/Research/projects/active/causal/manuscripts/upstream_downstream/programs-simulation/install_packages.R saulb@killdevil.unc.edu:/netscr/saulb/Rscripts

# THEN
# module add r/3.3.1
# bsub –Ip –q day R
# source('/netscr/saulb/Rscripts/install_packages.R')