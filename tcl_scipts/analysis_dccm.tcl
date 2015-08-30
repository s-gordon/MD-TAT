#!/usr/bin/env tclsh
# AUTHOR:   Shane Gordon
# FILE:     analysis_dccm.tcl
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


proc fitframes { molid seltext } {
  # This takes a selection and fits that selection for every frame in the
  # molecule (all atoms are moved, but the fit is based on the selection).
  #
  # For example:  fitframes top "protein"
  #
  # -Jim
  set ref [atomselect $molid $seltext frame 0]
  set sel [atomselect $molid $seltext]
  set all [atomselect $molid all]
  set n [molinfo $molid get numframes]

  for { set i 1 } { $i < $n } { incr i } {
    $sel frame $i
    $all frame $i
    $all move [measure fit $sel $ref]
  }
  $ref delete
  $sel delete
  $all delete
  return
}


proc align_dcd { dcd seltext } {
  set traj temp.dcd
  set dcd_i $dcd
  mol addfile $dcd_i waitfor all
  fitframes top $align_seltext
  animate write dcd $traj waitfor all top
  animate delete all
}

# }}}

# cross_corr2 ------------------------------------------------------------- {{{

#
# This version of the cross-correlation script is capable of processing very
# large trajectories.  It uses the embedded bigdcd script to process each frame
# separately.
# Usage:
# set sel [atomselect 0 "valid selection"]
# bigcorr $sel ccmat cmlast reslist sarray file1.dcd file2.dcd ...
# bigdcd_wait
# cross_corr_finis file.name ccmat reslist
#
global ccmat cmlast reslist sarray
#
proc bigcorr { selection ccmat cmlast reslist sarray args } {
  global bigdcd_frame vmd_frame bigdcd_running
  global ocnt 
  upvar #0 $ccmat ccor
  upvar #0 $cmlast comlast
  upvar #0 $reslist rlist
  upvar #0 $sarray sel
  #
  # Initialize the arrays for the cross-correlation calculation
  #
  set ocnt 0
  #
  # get the list of residues
  #
  set rlist [lsort -integer -unique -index 0 [$selection get {residue segname resid}]]
  #
  # get the atom selections
  #
  foreach res $rlist {
    set r [lindex $res 0]
    set sel($r) [atomselect 0 "residue $r"]
    $sel($r) global
    #
    # initialize the cross-correlation matrix
    #
    foreach resprime $rlist {
      set rprime [lindex $resprime 0]
      if { $rprime > $r } { break }
      set ccor($r,$rprime) 0.0
    }
  }
  #
  # set initial values for the callback procedure 
  #
  set bigdcd_running 1
  set bigdcd_frame 0
  #
  # set a trace on each frame read and then
  # loop over the input dcd files  
  #
  uplevel #0 trace variable vmd_frame w bigcorr_frame
  foreach dcd $args {
    mol addfile $dcd type dcd waitfor 0
  }
  after idle bigdcd_wait
}

proc bigcorr_frame { tracedvar mol op } {
  global bigdcd_frame vmd_frame
  global ccmat cmlast reslist sarray
  set msg {}
  #
  # compute the cross-correlation for this frame
  #
  set thisframe $vmd_frame($mol)
  incr bigdcd_frame
  if { [catch {cross_corr_frame $bigdcd_frame ccmat cmlast reslist sarray} msg] } { 
    puts stderr "bigdcd aborting at frame $bigdcd_frame\n$msg"
    bigdcd_done
    return
  }
  #
  # delete this frame
  #
  animate delete beg 0 end 0 top
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

  proc cross_corr_frame {nframe ccmat cmlast reslist sarray} { 
  #
  # Makes a matrix of all residue center-of-mass to residue center-of-mass
  # cross correlations.
  #
  global ocnt bigdcd_frame 
  upvar #0 $ccmat ccor
  upvar #0 $cmlast comlast
  upvar #0 $reslist rlist
  upvar #0 $sarray sel
  #
  # print the frame number to show progress (this routine takes a while)
  #
  if { ![expr $bigdcd_frame%10] } {
    incr ocnt
    if { $ocnt < 10 } {    
      puts -nonewline "$bigdcd_frame.."
    } else {
      set ocnt 0
      puts "$bigdcd_frame.."
    }
  }
  #
  # get the initial center-of-mass positions
  #
  if { $bigdcd_frame == 1 } {
    foreach res $rlist {
      set r [lindex $res 0]
      set comlast($r) [measure center $sel($r) weight mass]
    }
  } else {
  #
  # Compute the center of mass for each residue
  # and the delta from the previous position.
  # 
    foreach res $rlist {
      set r [lindex $res 0]
      set com($r) [measure center $sel($r) weight mass]
      set delta($r) [vecnorm [vecsub $com($r) $comlast($r)]]
      #
      # compute the cross-correlation
      #
      foreach resprime $rlist {
        set rprime [lindex $resprime 0]
        if { $rprime > $r } { break }
        set ccor($r,$rprime) [expr $ccor($r,$rprime) + \
          [vecdot $delta($r) $delta($rprime)]] 
      }
      #
      # update the last values
      #
      set comlast($r) $com($r)
    }
  }
  } 

  proc cross_corr_finis {filename ccmat reslist} { 
  #
  # Computes the normalized matrix and prepares it for plotting.
  # Saves the output to a file. 
  # Symmetrize and normalize the correlation matrix
  #
  global bigdcd_frame 
  upvar #0 $ccmat ccor
  upvar #0 $reslist rlist
  #
  # normalize the results
  #
  set norm [expr (1./($bigdcd_frame-1.0))]
  foreach res $rlist {
    set r [lindex $res 0]
    foreach resprime $rlist {
      set rprime [lindex $resprime 0]
      if { $rprime > $r } { break }
      set ccor($r,$rprime) [expr $norm * $ccor($r,$rprime)]
      #
      # for plotting purposes, zero out negative entries above the
      # diagonal and positive entries below
      #
      if { $ccor($r,$rprime) < 0.0 } {
        set ccor($rprime,$r) $ccor($r,$rprime)
        set ccor($r,$rprime) 0.0
      } elseif { $r != $rprime } {
        set ccor($rprime,$r) 0.0
      }
    }
  }
  #
  # open the file for output
  #
  set outfile [open $filename w] 
  #  
  # loop over each residue and print the matrix 
  set rlast [lindex [lindex $rlist end] 0]
  foreach res $rlist { 
    set r [lindex $res 0]
    set ocnt 0
    foreach resprime $rlist {
      set rprime [lindex $resprime 0]
      incr ocnt
      if { $rprime != $rlast } {      
        puts -nonewline $outfile [format "%7.3f" $ccor($r,$rprime)]
      } else { 
        set ocnt 0
        puts $outfile [format "%7.3f" $ccor($r,$rprime)]
      }
    }
  }
  puts $outfile " " 
  # close the file
  flush $outfile
  close $outfile 
  return 1
  } 

  proc rd_ccmat {filename ccmat seglist reslist} {
  #
  # Read the cross-correlation matrix from a file
  #
  # Usage:  From VMD command window or Tkcon
  # set ccmat(1,1) 0
  # set seglist {}
  # set reslist {}
  # set nres [rd_ccmat $filename ccmat seglist reslist]    == Invoke this routine
  #
  upvar $ccmat c
  upvar $seglist slist
  upvar $reslist rlist
  #
  # open the input file
  #
  set infile [open $filename r]
  #
  # get the total number of residues
  #
  gets $infile line
  set nres [lindex $line 1]
  #
  # get the list of residues and segment names. These
  # were written 8 to a line.
  #
  set nline [expr $nres / 8 ]
  set nn [expr $nline * 8]
  if { $nn < $nline } {
    incr nline
  }
  set slist {}
  set rlist {}
  set ocnt 7
  for {set j 0} {$j < $nres} {incr j} {
    incr ocnt
    if { $ocnt == 8 } {
      gets $infile line
      set ocnt 0
    }  
    lappend slist [lindex $line [expr 2*$ocnt]]
    lappend rlist [lindex $line [expr 2*$ocnt + 1]]
  }
  #
  # Skip the header
  #
  gets $infile line
  #
  # Read the matrix elements
  #
  set j 0
  set k 0
  while { 1 } {
    gets $infile line
    if { [eof $infile] } {break}
    set M [llength $line]
    for {set m 0} {$m < $M} {incr m} {
      set c($j,$k) [lindex $line $m]
      incr k
      if { $k == $nres } {
        set k 0
        incr j
      }   
    }
  }
  # close the file
  close $infile
  return $nres
  }

  proc draw_corr { ccmat seglist reslist color thresh } {
  #
  upvar $ccmat c
  upvar $seglist slist
  upvar $reslist rlist
  #
  set nres [llength $rlist]
  #
  # Connect the CAs that are correlated above the threshold
  # ATOMSELECT is expensive, so only get the positions of CAs
  # that are above the threshold.
  # Mark the backbone atoms as not found
  #
  for {set j 0} {$j < $nres} {incr j} {
    set ca($j) -9999
  }  
  set negthresh [expr $thresh * -1.0]
  draw color $color
  for {set j 0} {$j < $nres} {incr j} {
    for {set k [expr $j + 1]} {$k < $nres} {incr k} {
    #
    # Note:  The correlation matrix is split into positive (upper)
    # and negative (lower) halves for plotting purposes. Summing
    # upper and lower yields the proper matrix element.
    #
    set elem [expr $c($j,$k) + $c($k,$j)]
    if { ($elem > $thresh) || ($elem < $negthresh) } {
      if { $ca($j) == -9999 } {
        set ca($j) [atomselect top "segname [lindex $slist $j] \
          and resid [lindex $rlist $j] and name CA C5'"]
        if { [llength [$ca($j) list]] == 1 } {
          set xca($j) [eval "vecscale [$ca($j) get {x y z}] 1.0"]
        }  
      }   
      if { $ca($k) == -9999 } {
        set ca($k) [atomselect top "segname [lindex $slist $k] \
          and resid [lindex $rlist $k] and name CA C5'"]
        if { [llength [$ca($k) list]] == 1 } { 
          set xca($k) [eval "vecscale [$ca($k) get {x y z}] 1.0"]
        }
      }   
      if { [llength [$ca($j) list]] == 1 && [llength [$ca($k) list]] == 1 } {
        draw line $xca($j) $xca($k) width 2     
      }  
    }
    }
  }
  }  

  proc corr_check { ccmat jclist kclist nres thresh } {
  #
  upvar $ccmat c
  upvar $jclist jlist
  upvar $kclist klist
  #
  set negthresh [expr $thresh*-1.0]
  for {set j 0} {$j < $nres} {incr j} {
    for { set k 0} {$k < $nres} {incr k} {
      if { $j != $k } {
        if { $c($j,$k) > $thresh } {
          lappend jlist $j
          lappend klist $k
        }
        if { $c($j,$k) < $negthresh } {
          lappend jlist $j
          lappend klist $k
        }  
      }  
    }
  }
}          

# }}}

# Sequentially read in args as $1, $2, etc.
regexp {0.[0-9]{1,3}} $dcd index_no

# Variable definitions for later
set input_t temp.dcd
set input $dcd
set o dccm_$index_no.dat
set out_dir $raw/sim_$index_no
set mol [ mol new $input_psf.psf type psf waitfor all ]
set reference_CA [ atomselect $mol $seltext_CA frame 0 ]
set sel_protein [ atomselect $mol $seltext_protein ]
set sel_backbone [ atomselect $mol $seltext_backbone ]
set sel_CA [ atomselect $mol $seltext_CA ]
mol addfile $input_psf.pdb first 0 last 0 waitfor all

# dccm
filecheck $out_dir/$o
set mol [ mol new $input_psf.psf type psf waitfor all ]
align_dcd $dcd $seltext
set sel [ atomselect top $seltext_protein ]
bigcorr $sel ccmat cmlast reslist sarray $input_t
bigdcd_wait
cross_corr_finis $o ccmat reslist
file rename -force $o $out_dir/$o

exit
