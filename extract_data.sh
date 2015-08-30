#!/usr/bin/env bash
# AUTHOR:   Shane Gordon/Mike Kuiper
# FILE:     a1_extract_all_my_data.sh
# ROLE:     This is a simple launch script to run: 
#           <main>/Scripts/Analysis_Scripts/a1_extract_all_data
#           This will extract and collate all jobs in the job directories under ../MainJob_dir 
# CREATED:  2014-10-27 22:43:42

vmd -dispdev text -e tcl_scipts/a1_extract_all_dcd_data.tcl >/dev/null
tcl_scipts/a-1_extract_all_dcd_data
