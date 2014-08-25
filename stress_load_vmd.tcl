set numframes [molinfo top get numframes]
set numatoms [molinfo top get numatoms]
set fp [open "stress_output.dat" r]
echo $numatoms

for {set i 0} {$i < $numframes} {incr i} {
animate goto $i
set sel [atomselect top "all"] 
$sel frame $i

set chrg {}
set chrg [concat $chrg [gets $fp]]

$sel set user $chrg

}

mol modcolor 0 0 User 
mol colupdate 0 0 1 
mol scaleminmax 0 0 0.0 $numframes 
