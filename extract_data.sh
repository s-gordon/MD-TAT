#!/usr/bin/env bash
# AUTHOR:   Shane Gordon/Mike Kuiper
# CREATED:  2014-10-27 22:43:42

vmd -dispdev text -e tcl_scipts/extract_all_dcd_data.tcl
./tcl_scipts/a1_extract_all_dcd_data
