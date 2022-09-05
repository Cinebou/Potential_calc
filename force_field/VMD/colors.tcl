#!/usr/local/bin/vmd

#This proc writes colors of the atoms in your simulation.
#Frame is colored by elements and written in the script
#and adsorbate molecules are colored by Name. 
#To add a color simply input in the RGB code and add the atom to
#appropriate type.


proc colors {} {
color scale colors RWB {1.0 0.0 0.0} {1.0 1.0 1.0} {0.0 0.0 1.0}
	#set scale method RWB
	set colorcmds {
		{color Element {Cr} 0}
		{color Element {C} 1}
		{color Element {O} 2}
		{color Element {H} 3}
		{color Element {F} 4}
		{color Name {C} 5}
		{color Name {O} 6}
		{color Display {Background} white}
	}

	foreach colcmd $colorcmds {
		set val [catch {eval $colcmd}]
	}	
	color change rgb 0 0.2 0.8 0.8
	color change rgb 1 1.0 1.0 1.0
	color change rgb 2 1.0 0.6 0.0
	color change rgb 3 0.2 0.2 0.2
	color change rgb 4 1.0 0.6 0.8
	color change rgb 5 0.3 0.3 0.3
	color change rgb 6 1.0 0.0 0.0

}