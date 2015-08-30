#!/usr/bin/env tclsh
# AUTHOR:   Shane Gordon
# FILE:     a5_other_analyses.tcl
# ROLE:     TODO (some explanation)
# CREATED:  2014-06-03 22:04:49
# MODIFIED: 2015-06-17 16:53:45

# Common variables ----------------------------------------------------------- {{{

# Output directory for raw data files
set raw "raw_analysis_data"

# }}}

# ---------------------------------------------------------------------------- {{{

# read in raw data: 
puts " reading in reduced data set:"

proc dircheck { dirname } {
  if { [ file isdirectory $dirname ] } {
    file delete -force $dirname
    puts "Deleted $dirname"
  }
}

proc filecheck { filename } {
  if { [ file exists $filename ] } {
    file delete -force $filename
    puts "Deleted $filename"
  }
}

# If directory exists, blow it away and make a new one
# Else, make the directory
proc dir_make { dir } {
  if { [ file isdirectory $dir ] } {
    file delete -force $dir
    puts "Deleted $dir"
    file mkdir $dir
  } else {
    file mkdir $dir
  }
}

# load useful analysis script:  
source ../Scripts/Tcl_Scripts/analysis.tcl
source ../Scripts/Tcl_Scripts/bigdcd.tcl

# }}}

#------------------------------------------------------------------------------


foreach index [ lsort [glob no_water_*.dcd] ] {
  regexp {0.[0-9]{1,3}} $index index_no

  # Variable definitions for later
  set input "no_water_${index_no}"
  set input_psf "no_water"
  set out_dir "$raw/sim_$index_no"
  dir_make $out_dir
  set num_rmsf_windows 5
  set seltext_protein "protein"
  set seltext_backbone "protein and backbone"
  set seltext_CA "protein and name CA"
  set molid 0
  set mol [ mol new $input_psf.psf type psf waitfor all ]
  set reference_CA [ atomselect $mol "$seltext_CA" frame 0]
  set sel_protein [atomselect $mol "$seltext_protein"]
  set sel_backbone [atomselect $mol "$seltext_backbone"]
  set sel_CA [ atomselect $mol "$seltext_CA" ]
  mol addfile $input_psf.pdb first 0 last 0 waitfor all


  #------------------------------------------------------------------------------

  # RMSD scan
  # filecheck rmsd_protein_${index_no}.txt
  # bigdcd rmsdscan_bigdcd $input.dcd
  # bigdcd_wait
  # file rename "rmsd_protein.txt" "${out_dir}/rmsd_protein_${index_no}.txt"

  # Radius of gyration scan
  filecheck protein_radius_gyration_${index_no}.txt
  bigdcd rgyrscan_bigdcd $input.dcd
  bigdcd_wait
  file rename "protein_radius_gyration.txt" "${out_dir}/protein_radius_gyration_${index_no}.txt"

  # RMSF scan
  if { [file exists ./number_frames_${index_no}.txt ] } {
    set f [ open "number_frames_${index_no}.txt" r ]
    set r [ read $f ]
    close $f
    set num_frames $r
  } else {
    puts "Can't find number_frames_${index_no}.txt.\nCan't manage RMSF calculations without this info"
    exit
  }
  set rmsf_fraction_count 1
  set frame_incr [ expr $num_frames / $num_rmsf_windows ]

  while { $rmsf_fraction_count <= $num_rmsf_windows } {
    set first_frame [ expr ( $rmsf_fraction_count - 1 ) * $frame_incr ]
    set last_frame [ expr ( $first_frame + $frame_incr ) -1 ]
    mol addfile $input.dcd first $first_frame last $last_frame waitfor all
    fitframes top "protein and name CA"
    rmsfscan_range $sel_CA 1 -1 "${out_dir}/rmsf_protein_backbone_${rmsf_fraction_count}of${num_rmsf_windows}_${index_no}"
    animate delete beg 1 end -1 ;# Spare the first frame, containing ref pdb
    incr rmsf_fraction_count
  }

  # Calculates secondary structure over course of entire simulation
  # filecheck ${out_dir}/sec_structure_${index_no}.dat
  # dircheck ${out_dir}/SecondaryStructure
  # filecheck ${out_dir}/SecondaryStructure_${index_no}/betaPercent.plt
  # filecheck ${out_dir}/SecondaryStructure_${index_no}/coilPercent.plt
  # filecheck ${out_dir}/SecondaryStructure_${index_no}/helixPercent.plt
  # filecheck ${out_dir}/SecondaryStructure_${index_no}/turnPercent.plt
  # bigdcd ss_calc_bigdcd $input.dcd
  # bigdcd_wait
  # file rename -force "SecondaryStructure" "${out_dir}/SecondaryStructure_${index_no}"

  # Time-dependent SASA-scan for entire simulation
  # filecheck ${out_dir}/protein_sasa_${index_no}.txt
  # bigdcd sasa_scan_bigdcd $input.dcd
  # bigdcd_wait
  # file rename "protein_sasa.txt" "${out_dir}/protein_sasa_${index_no}.txt"

}

exit
