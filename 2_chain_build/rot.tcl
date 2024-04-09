#set sel [atomselect top "residue 0 to 108"]
set sel [atomselect top "residue 109 to 217"]
set com [measure center $sel weight mass]
set matrix [transaxis y 2]
#set traslation_matrix [transoffset {1 0 0}]
$sel moveby [vecscale -1.0 $com]
$sel move $matrix
$sel moveby $com
