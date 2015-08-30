#!/usr/bin/env tclsh
# AUTHOR:   Shane Gordon
# FILE:     analysis_rmsf.tcl
# ROLE:     TODO (some explanation)
# CREATED:  2015-06-18 20:43:31
# MODIFIED: 2015-07-12 10:59:46

# BASIC USAGE 

# Call this script from VMD using:
#
#   vmd -dispdev text -e <this-script.tcl> -args <traj-file> <out-dir> <sel>
#
# Where:
#   <traj-file> is a single trajectory file (must exit-slash any spaces in the
#   name)
#   <out-dir> is the directory into which the output data file will be moved
#   <sel> is an optional atom selection to be used for alignment

# VARIABLES --------------------------------------------------------------- {{{

set input_psf "no_water"
set seltext "protein"
set seltext_protein "protein"
set seltext_backbone "$seltext_protein and backbone"
set seltext_CA "$seltext_protein and name CA"
set molid 0
set dcd [ lindex $argv 0 ]
set raw [ lindex $argv 1 ]
set align_seltext [ lrange $argv 2 end ]

# }}}

# SUBROUTINES ------------------------------------------------------------- {{{


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


proc dir_make { dir } {
  # If directory exists, blow it away and make a new one
  # Else, make the directory
  if { [ file isdirectory $dir ] } {
    puts "Directory $dir already exists!"
  } else {
    file mkdir $dir
  }
}


proc rmsfscan { sel fname } {
  # per residue rmsf to file
  set rlist [$sel get resid]
  set rmsf [measure rmsf $sel]
  set s [lindex $rlist 0]
  set n [lindex $rlist end]
  set f [open $fname.txt "w"]
  for {set i $s} {$i <= $n} {incr i} {
    puts $f "$i [lindex $rmsf [expr $i - $s]]"
  }
  close $f
}


proc rmsfscan_range { sel start_frame end_frame fname } {
  # per residue rmsf to file
  set rlist [$sel get resid]
  set rmsf [measure rmsf $sel first $start_frame last $end_frame ]
  set s [lindex $rlist 0]
  set n [lindex $rlist end]
  set f [open $fname.txt "w"]
  for {set i $s} {$i <= $n} {incr i} {
    puts $f "$i [lindex $rmsf [expr $i - $s]]"
  }
  close $f
}

# }}}

# Sequentially read in args as $1, $2, etc.
regexp {0.[0-9]{1,3}} $dcd index_no

# Variable definitions for later
set input $dcd
set out_dir $raw/sim_$index_no
set mol [ mol new $input_psf.psf type psf waitfor all ]
set reference_CA [ atomselect $mol $seltext_CA frame 0 ]
set sel_protein [ atomselect $mol $seltext_protein ]
set sel_backbone [ atomselect $mol $seltext_backbone ]
set sel_CA [ atomselect $mol $seltext_CA ]
set num_rmsf_windows 5
mol addfile $input_psf.pdb first 0 last 0 waitfor all

# RMSF scan
if { [file exists ./number_frames_$index_no.txt ] } {
  set f [ open "number_frames_${index_no}.txt" r ]
  set r [ read $f ]
  close $f
  set num_frames $r
} else {
  puts "Can't find number_frames_${index_no}.txt"
  puts "Can't manage RMSF calculations without this info"
  exit
}
set rmsf_fraction_count 1
set frame_incr [ expr $num_frames / $num_rmsf_windows ]

while { $rmsf_fraction_count <= $num_rmsf_windows } {
  set first_frame [ expr {($rmsf_fraction_count - 1) * $frame_incr} ]
  set last_frame [ expr {($first_frame + $frame_incr ) - 1} ]
  mol addfile $dcd first $first_frame last $last_frame waitfor all
  fitframes top "$align_seltext and name CA"
  rmsfscan_range $sel_CA 1 -1 $out_dir/rmsf_protein_backbone_${rmsf_fraction_count}of${num_rmsf_windows}_$index_no
  animate delete beg 1 end -1 ;# Spare the first frame, containing ref pdb
  incr rmsf_fraction_count
}
mol addfile $dcd first 0 last -1 waitfor all
fitframes top $seltext_CA
rmsfscan_range $sel_CA 1 -1 $out_dir/rmsf_all_protein_backbone_$index_no
animate delete beg 1 end -1 ;# Spare the first frame, containing ref pdb

exit
