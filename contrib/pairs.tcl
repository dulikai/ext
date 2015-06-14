#! tcl-vmd


atomselect macro RING "name NA1 NA2 CW1 CW2 CR1"
atomselect macro CARB "name O21 O22 C"

# select the most close cation : 
# sort it and get close n
# then select the most close anion : 
# sort it

# select near residue for grp1 in 'grp2'
# return a list of neighbor along with min-dist within xdist distance..

proc near_res { grp1 grp2 cutoff nres } {
    set sel0 [atomselect top $grp1]
    # find neighbor
    set sel1 [atomselect top "($grp2 and not $grp1) and (RING or CARB) and same residue as within $cutoff of $grp1"]
    set myids [lsort -unique [$sel1 get resid]]
    set num [llength $myids]
    puts "~ $num neighbors found"
    # then, measure the neighbors
    foreach id $myids {
        set sel [atomselect top "resid $id and (CARB or RING)"]
        set mydist [min_dist $sel0 $sel]
        set i_dist [lindex $mydist 0]
        set sm [lindex $i_dist 4]
        lappend mylist [list $id $sm]
        $sel delete
    }    
    $sel0 delete
    # the list
    set mylist [lsort -real -index 1 $mylist]
    # puts $mylist
    set ionlist [lrange $mylist 0 $nres-1]
    set myresid {}
    set aver 0.0
    if { $nres == 0 } { set nres 1 }
    foreach i_list $ionlist {
        set id [lindex $i_list 0]
        set aver [expr $aver + [lindex $i_list 1]/$nres]
        lappend myresid $id
    }
    lappend myresid $aver
    return $myresid
}


proc find_pairs { type0 type1 type2 cutoff nres } {
    # type 0: resid 1; type 1: resname Bmim; type 2 resname Gly    
    set sel [atomselect top "$type0"]
    set myids [lsort -unique [$sel get resid]]

    set pairs {}
    foreach id $myids {
        set grp1 "resid $id"
        set grp2 $type1
        # -1 because type 0 is same as type 1
        set myid1 [near_res  $grp1 $grp2 $cutoff [expr $nres-1]]
        # puts "why $id $myid1 "
        set aver1 [lindex $myid1 end] 
        set myid1 [lrange $myid1 0 end-1]
        set newgrp "resid $myid1"
        set grp2 $type2
        set myid2 [near_res  $newgrp $grp2 $cutoff $nres]  
        set aver2 [lindex $myid2 end]
        set myid2 [lrange $myid2 0 end-1]
        set aver [expr ($aver1+$aver2)/2.0]
        set newgrp "$myid2 $myid1"
       
        set one_list [lappend newgrp $id $aver]
        lappend pairs $one_list 
    }
    $sel delete
    set pairs [lsort -real -index end $pairs]
    # puts $pairs
    return $pairs
}

proc min_dist { sel1 sel2 } {
    set myndx1 [lsort -unique [$sel1 get index]]
    set myndx2 [lsort -unique [$sel2 get index]]
    set i_res [lsort -unique [$sel1 get resid]]
    set j_res [lsort -unique [$sel2 get resid]]
    set mydist {}
    foreach i_ndx $myndx1 {
        foreach j_ndx $myndx2 {
            set d [measure bond [list $i_ndx $j_ndx]]
            lappend mydist [list $i_res $j_res $i_ndx $j_ndx $d]
            # puts "$i_ndx $j_ndx"
        }
    }
    set mydist [lsort -real -index 4 $mydist]
    return $mydist    
}


proc dump_xyz { pairs type1 type2 filename } {
    # working directory
    set mydir "mydir"
    file mkdir $mydir
    set i_frame [molinfo top get frame]
    foreach one_list $pairs {
        set aver [lindex $one_list end]
        set mycenter [lindex $one_list end-1]
        set myresid [lrange $one_list 0 end-1]
        
        set sel [atomselect top "resid $myresid"]
        set res [lsort -unique [$sel get resid]]
        # puts $res
        puts ">>($i_frame) resid $myresid"
        puts "aver dist: $aver"
        set myfile "${filename}.${i_frame}.${mycenter}"
        $sel writexyz "./${mydir}/${myfile}"
        # addtional info.
        set sel [atomselect top "resid $myresid and $type1"]
        set n_cation [llength [lsort -unique [$sel get resid]]]
        set cation [lsort -unique [$sel get resid]]
        set sel [atomselect top "resid $myresid and $type2"]
        set n_anion [llength [lsort -unique [$sel get resid]]]
        set anion [lsort -unique [$sel get resid]]
        $sel delete
        
        set myfile "typenumber.${i_frame}.${mycenter}"
        set out [open "${mydir}/$myfile" w]
        puts $out "FRAME:$i_frame $n_cation $n_anion"
        puts $out "TOT: resid $cation $anion"
        puts $out "TYPEA: resid $cation"
        puts $out "TYPEB: resid $anion"
        close $out       
    }
    return
}

proc mycat { prefix } {
    file delete $prefix

    set mydir "mydir"
    set out [open $prefix w]
    fconfigure $out -translation binary
    foreach fname [glob -nocomplain -type f "./${mydir}/${prefix}.*"] {
        set in [open $fname]
        fconfigure $in -translation binary
        fcopy $in $out
        close $in
#        file delete $fname
    }
    close $out
}


proc get_center { type1 rad } {
    set sel [atomselect top $type1]
    set center [measure center $sel]
    set x [lindex $center 0]
    set y [lindex $center 1]
    set z [lindex $center 2]
    set r2 [expr $rad * $rad]
    set sel [atomselect top "resname Bmim and (x-$x)*(x-$x)+(y-$y)*(y-$y)+(z-$z)*(z-$z) < $r2"]
    set res [lsort -unique [$sel get resid]]
    set num [llength $res]
    # puts "mynum $num, $res"
    
    return $res  
}

proc traj_find { type1 type2 cutoff nres} {
    set sel [atomselect top $type1]
    set n [molinfo top get numframes]
    # set sel [atomselect top "resid 100"]
    for { set i 0 } { $i < $n } { incr i } {
        animate goto $i
        $sel update
        
        set myid [get_center $type1 10.0]
        set num [llength $myid]
        set ran [expr {int(rand()*($num-1))}]
        set rid [lindex $myid $ran]
        puts "my random number, resid: $ran $rid: $num"
        set type0 "resid $rid"
        
        puts "$i th jobs"
        set mypairs [find_pairs $type0 $type1 $type2 $cutoff $nres ]
        dump_xyz $mypairs $type1 $type2 "dump.xyz"
    }
}



# main prog
#mol new md.gro
#mol addfile trajout.xtc waitfor -1
source traj.tcl


set filename "dump.xyz"
set type1 "resname Bmim"
set type2 "resname Gly"
set nres 5
set cutoff 10.0

traj_find $type1 $type2 $cutoff $nres

mycat $filename
mycat "typenumber"
#file delete -force "mydir"

mol delete top


