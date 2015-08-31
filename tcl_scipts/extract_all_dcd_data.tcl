#!/usr/bin/env tclsh
# AUTHOR:   Shane Gordon
# FILE:     extract_all_dcd_data
# ROLE:     TODO (some explanation)
# CREATED:  2014-06-04 20:16:45
# MODIFIED: 2015-05-12 10:51:26

# BASIC USAGE 

# Call this script from VMD using:
#
#   vmd -dispdev text -e <this-script.tcl>

# SUBROUTINES ------------------------------------------------------------- {{{

proc reduced { sel fname } {
  # Write the indexes for the selection, sel, to a file fname.text and
  # also write a corresponding psf file fname.psf
  set f [open $fname.text "w"]
  puts -nonewline $f [$sel get index]
  close $f
  $sel writepsf $fname.psf
  $sel writepdb $fname.pdb
}

# }}}

mol new [ glob ../InputFiles/*psf ]
mol addfile [ glob ../InputFiles/*pdb ]

set sel [ atomselect top "protein" ]
reduced $sel no_water

exit
