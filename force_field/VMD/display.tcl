#!/usr/local/bin/vmd

#set the display properties here. 
#for high quality rendering it is advised to keep 
#shadows and ambient occulsion on at the given values but 
#you can play around and decide the best combination for you

proc display_single_cell {} {
display 
axes location off
light 0 on
light 1 on
light 2 off
light 3 off
}

proc display_multi_cell {} {
display eyesep 0.065000
display focallength 2.000000
display height 10.000000
display distance -2.000000
display projection Orthographic
display nearclip set 0.010000
display farclip  set 10.000000
display depthcue on
display cuestart 0.500000
display cueend 10.000000
display cuestart 0.500000
display cueend 10.000000
display cuemode Exp2
display	cuedensity 0.47
display shadows on
display ambientocclusion on
display aoambient 0.80
display aodirect 0.30
axes location off
light 0 on
light 1 on
light 2 off
light 3 off
}