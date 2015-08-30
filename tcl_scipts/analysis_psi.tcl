#!/usr/bin/env tclsh
# AUTHOR:   Shane Gordon
# FILE:     analysis_psi.tcl
# ROLE:     TODO (some explanation)
# CREATED:  2015-06-18 20:48:18

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

# BigDCD-v2 --------------------------------------------------------------- {{{

# Justin Gullingsrud
# updates by Axel Kohlmeyer
# vmd@ks.uiuc.edu

proc bigdcd { script type args } {
  global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame bigdcd_running

  set bigdcd_running 1
  set bigdcd_frame 0
  set bigdcd_firstframe [molinfo top get numframes]
  set bigdcd_proc $script

  # backwards "compatibility". type flag is omitted.
  if {[file exists $type]} { 
    set args [linsert $args 0 $type] 
    set type auto
  }

  uplevel #0 trace variable vmd_frame w bigdcd_callback
  foreach dcd $args {
    if { $type == "auto" } {
      mol addfile $dcd waitfor 0
    } else {
      mol addfile $dcd type $type waitfor 0
    }
}
after idle bigdcd_wait
}

proc bigdcd_callback { tracedvar mol op } {
  global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame
  set msg {}

  # If we're out of frames, we're also done 
  # AK: (can this happen at all these days???). XXX
  set thisframe $vmd_frame($mol)
  if { $thisframe < $bigdcd_firstframe } {
    puts "end of frames"
    bigdcd_done
    return
}

incr bigdcd_frame
if { [catch {uplevel #0 $bigdcd_proc $bigdcd_frame} msg] } { 
  puts stderr "bigdcd aborting at frame $bigdcd_frame\n$msg"
  bigdcd_done
  return
}
animate delete beg $thisframe end $thisframe $mol
return $msg
}

proc bigdcd_done { } {
  global bigdcd_running

  if {$bigdcd_running > 0} then {
    uplevel #0 trace vdelete vmd_frame w bigdcd_callback
    puts "bigdcd_done"
    set bigdcd_running 0
}
}

proc bigdcd_wait { } {
  global bigdcd_running bigdcd_frame
  while {$bigdcd_running > 0} {
    global bigdcd_oldframe
    set bigdcd_oldframe $bigdcd_frame
    # run global processing hooks (including loading of scheduled frames)
    display update ui
    # if we have read a new frame during then the two should be different.
    if { $bigdcd_oldframe == $bigdcd_frame } {bigdcd_done}
}
}



# }}}

# VARIABLES --------------------------------------------------------------- {{{

set input_psf "no_water"
set seltext "protein"
set seltext_protein "protein"
set seltext_backbone "protein and backbone"
set seltext_CA "protein and name CA"
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


proc psi_bigdcd { frame } {
  global sel_CA
  set f [ open "psi.txt" a ]
  set dihedral_long [ $sel_CA get {psi} ]
  foreach dh $dihedral_long { lappend dihedral [format "%.2f" $dh] }
  puts $f "$frame $dihedral"
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
mol addfile $input_psf.pdb first 0 last 0 waitfor all

# Time-dependent SASA-scan for entire simulation
filecheck $out_dir/phi.txt
bigdcd psiscan_bigdcd $input
bigdcd_wait
file rename -force psi.txt $out_dir/psi.txt

exit
