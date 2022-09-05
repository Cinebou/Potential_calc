#!/usr/local/bin/vmd
source [file join [file dirname [info script]] rendering.tcl]

#Use this file to load the molecules in VMD
#For importing steps of MD make changes here.

proc loadFrame {frame {res 12} {first 0} {last -1} {step 1}} {
	mol new $frame.pdb first $first last $last step $step filebonds 1 autobonds 1 waitfor all
	mol delrep 0 top
	#mol representation Licorice 0.300000 $res $res
	#mol representation VDW 0.500000 $res
	mol representation Bonds 0.1
	set frame [atomselect top "all"]
	set C_frame [atomselect top "all {type C}"]
	set O_frame [atomselect top "all {type O}"]
	mol color Element
	mol material AOChalky
	mol addrep top
	mol rename top frame
}

proc loadMolecules {frame t p {res 12} {first 0} {last -1} {step 1}} {
	set files [lsort [glob $t\_$p\.pdb]]
	set count [molinfo top]
	foreach f $files {
		incr count
		mol new $f first $first last $last step $step filebonds 1 autobonds 1 waitfor all
		mol delrep 0 top
		#atomselect top "all name C"
		mol representation VDW 0.5 $res
		mol color Name
		mol material AOChalky
		mol addrep top
		mol top $count
		set molName [concat $frame-$t\K-$p\Pa]
		mol rename top $molName
		#rendering $frame-$t\K-$p\Pa
		mol showrep $count 0 off	
	}
}

proc loadFrame_multi {frame {res 12} {first 0} {last -1} {step 1}} {
	mol new $frame.pdb first $first last $last step $step filebonds 1 autobonds 1 waitfor all
	mol delrep 0 top
	#mol representation Licorice 0.300000 $res $res
	mol representation Lines 2.000
	atomselect macro frame {all}
	mol color Element
	mol material AOChalky
	mol addrep top
	set selframe [atomselect top "all"]
	$selframe moveby {12.8345 8.2 0}
	mol showperiodic top 0 XY 
	mol numperiodic top 0 1
	mol rename top frame
	
}

proc loadMolecules_multi {frame t p {res 12} {first 0} {last -1} {step 1}} {
	set files [lsort [glob Movie_$frame\_1.1.1_$t\.000000_$p\.000000_component_CO2_0.pdb]]
	foreach f $files {
		incr count
		mol new $f first $first last $last step $step filebonds 1 autobonds 1 waitfor all
		mol delrep 0 top
		atomselect top "all name C"
		mol representation CPK 2.0 0.3 $res $res
		mol color Name
		mol material AOChalky
		mol addrep top
		set selmol [atomselect top "all"]
		$selmol moveby {12.8345 8.2 0}
		mol showperiodic top 0 XY
		mol numperiodic top 0 1
		#mol top $count
		set molName [concat $frame-$t\K-$p\Pa]
		mol rename top $molName			
	}
}