#!/usr/local/bin/vmd

set base [molinfo top]
mol delete all

set mof MIL-101-primitive
set resolution 202
set first 0 
set last 0
set step 1

#temp in K and press in Pa
set Temp {298.15}
set Press {100000}


source [file join [file dirname [info script]] displaying.tcl]
source [file join [file dirname [info script]] colors.tcl]
source [file join [file dirname [info script]] loading.tcl]
source [file join [file dirname [info script]] rendering.tcl]

display_multi_cell
colors

#loadFrame $mof $res $first $last $step (default: loadFrame $mof 12 0 -1 1)
#loadMolecules $mof $t $p $res $first $last $step 
#(default: loadMolecules $mof $t $p 12 0 -1 1)

loadFrame $mof $resolution $first $last $step
set count [molinfo top]
foreach t $Temp {
	foreach p $Press {
		loadMolecules $mof $t $p $resolution $first $last $step
		incr count
		#rendering_high_res multi-$mof-$t\K-$p\Pa
		#mol showrep $count 0 off				
	}
}


