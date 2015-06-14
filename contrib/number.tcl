#! tcl-vmd


atomselect macro RING "name NA1 NA2 CW1 CW2 CR1"
atomselect macro CARB "name O21 O22 C CTA N3 HN31 HN32"
atomselect macro CARB "name O21 O22 C"



# select near residue for grp1 in 'grp2'
# return a list of neighbor along with min-dist within xdist distance..
proc near_res { grp1 grp2 cutoff nres } {
    # init
    set sel0 [atomselect top $grp1]    
    puts "working for $grp1"
    set c0 [measure center $sel0]
    # find neighbor
    set sel1 [atomselect top "($grp2 and (CARB or RING)) and same residue as within $cutoff of $grp1"]
    set myids [lsort -unique [$sel1 get resid]]
    set num [llength $myids]
    puts "~ $num neighbors ($grp2) find for $grp1"
    # loop
    set mylist {}
    foreach id $myids {
        set sel [atomselect top "(CARB or RING) and resid $id"]
        set mydist [min_dist $sel0 $sel]
        set i_dist [lindex $mydist 0]
        set sm [lindex $i_dist 4]
        set c1  [measure center $sel]
        set sc [vecdist $c0 $c1]
        #puts $sm
        # puts "$id $sm"
        lappend mylist [list $id $sm]
        $sel delete
    }    
    $sel0 delete

    set mylist [lsort -real -index 1 $mylist]
    set ionlist [lrange $mylist 0 $nres-1]
    puts [llength $ionlist]
    return $ionlist
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


proc find_pairs { type1 type2 cutoff nres } {
    set sel [atomselect top "$type1"]
    set myids [lsort -unique [$sel get resid]]

    set grp2 $type2
    set pairs {}
    foreach id $myids {
        set grp1 "resid $id"
        set mydist [near_res $grp1 $grp2 $cutoff $nres]
        # set one_dist [lindex $mydist 0]
        # the first n residue selected
        set one_list [lrange $mydist 0 end]
        set aver 0.0
        foreach ilist $one_list {
            set dist [lindex $ilist end]
            set aver [expr $aver + $dist / $nres]
        }
        set one_list [linsert $one_list end $id $aver]
        # puts $one_list
        # set one_dist [linsert $one_dist 0 $id]
        lappend pairs $one_list
    }
    $sel delete
    set pairs [lsort -real -index end $pairs]
    #puts [lindex $pairs 0]
    return $pairs
}


proc dump_xyz { pairs type1 type2 filename } {
    # working directory
    set mydir "mydir"
    #file mkdir $mydir
    set i_frame [molinfo top get frame]
    foreach one_list $pairs {
        set myneighbor [lrange $one_list 0 end-2]
        set mycenter [lindex $one_list end-1]
        set aver [lindex $one_list end]
        # puts $myneighbor
        set other $mycenter
        foreach i_list $myneighbor {
            set i [lindex $i_list 0]
            set d [lindex $i_list 1]
            append other " $i"
            # puts "$i $d"
            #puts "Hello: xx $i xx $d xx"
        }
        #puts [llength $other]
        set sel [atomselect top "resid $other"]
        set res [lsort -unique [$sel get resid]]
        #puts $res
        puts "resid $other"
        puts "aver dist: $aver"
        set myfile "${filename}.${i_frame}.${mycenter}"
        $sel writexyz "${mydir}/${myfile}"
        
        # additional info.
        set sel [atomselect top "$type1"]
        set res [lsort -unique [$sel get resname]]
        set sel [atomselect top "resid $other and resname $res"]
        set n_cation [llength [lsort -unique [$sel get resid]]]
        set cation [lsort -unique [$sel get resid]]
        set sel [atomselect top "resid $other and $type2"]
        set n_anion [llength [lsort -unique [$sel get resid]]]
        set anion [lsort -unique [$sel get resid]]
        # set sel [atomselect top "resid $other"]
        # set n_atom [llength [$sel get index]]
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
        #file delete $fname
    }
    close $out
}

proc get_center { type1 rad } {
    set sel [atomselect top $type1]
    $sel update
    set center [measure center $sel]
    set x [lindex $center 0]
    set y [lindex $center 1]
    set z [lindex $center 2]
    set r2 [expr $rad * $rad]
    set sel [atomselect top "$type1 and (x-$x)*(x-$x)+(y-$y)*(y-$y)+(z-$z)*(z-$z) < $r2"]
    set res [lsort -unique [$sel get resid]]
    set num [llength $res]
    # puts "mynum $num, $res"
    
    return $res  
}


proc traj_find { type1 type2 cutoff nres filename } {
    # working directory
    set mydir "mydir"
    if { [file isdirectory $mydir] } {
        file delete -force $mydir
    }
    file mkdir $mydir

    set sel [atomselect top $type1]
    set n [molinfo top get numframes]
    # set sel [atomselect top "resid 100"]
    for { set i 0 } { $i < $n } { incr i } {
        animate goto $i
        $sel update
        
        set myid [get_center $type1 8.0]
        set num [llength $myid]
        set ran [expr {int(rand()*($num-1))}]
        set rid [lindex $myid $ran]
        puts "my random number, resid: $ran $rid: $num"
        set type0 "resid $rid"

        set typex $type0
        
        puts "$i th jobs"
        set mypairs [find_pairs $typex $type2 $cutoff $nres]
        dump_xyz $mypairs $typex $type2 $filename
    }
}



# main prog
#set datadir ..
#mol new "$datadir/md.gro"
#mol addfile "$datadir/wrap.xtc" waitfor -1
source traj.tcl

if { [file exists "config.tcl"] } {
source config.tcl
} else {
set filename "dump.xyz"
set type1 "resname Bmim"
set type2 "resname Gly"
set nres 3
set cutoff 8.0
}

traj_find $type1 $type2 $cutoff $nres $filename

mycat $filename
mycat "typenumber"
#file delete -force "mydir"

mol delete top

exit



