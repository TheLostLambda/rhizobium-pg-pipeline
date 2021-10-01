#!/usr/bin/env bash

# 1) Run the MS1 search
# NOTE: Be sure to update the csv and ftrs filepaths in the script!
python run_pgfinder.py

# 2) Add any structures of interest to the selected_structures file

# 3) Generate graph files from the selected structures
python build_graphs.py Data/Outputs/selected_structures.csv Data/Outputs/Graphs

# 4) Run the MS2 search
# NOTE: Be sure to update the data (byspec) filepath in the script!
python run_ms2.py
