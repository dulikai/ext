#! tcl-vmd


atomselect macro RING "name NA1 NA2 CW1 CW2 CR1"
atomselect macro CARB "name O21 O22 C"


set grp1 "RING and resid 1"
set grp2 "CARB"

proc once { grp1 grp2 } {
    # init
    set sel0 [atomselect top $grp1]    
    puts "working for $grp1"
    set c0 [measure center $sel0]
    # find neighbor
    set sel1 [atomselect top "$grp2 and same residue as within 5.0 of $grp1"]
    set myids [lsort -unique [$sel1 get resid]]
    set num [llength $myids]
    puts "~ $num neighbors ($grp2) find for $grp1"
    # loop
    set mylist {}
    foreach id $myids {
    set sel [atomselect top "CARB and resid $id"]
    set c1  [measure center $sel]
    set s [vecdist $c0 $c1]
    #puts "$id $s"
    lappend mylist [list $id $s]
    }    
    set mylist [lsort -real -index 1 $mylist]


    return $mylist
}

proc min_dist { grp1 grp2 } {
    set sel1 [atomselect top $grp1]
    set sel2 [atomselect top $grp2]
    set myndx1 [lsort -unique [$sel1 get index]]
    set myndx2 [lsort -unique [$sel2 get index]]
    set mydist {}
    foreach i_ndx $myndx1 {
        foreach j_ndx $myndx2 {
            set d [measure bond [list $i_ndx $j_ndx]]
            lappend mydist [list $grp1 $grp2 $i_ndx $j_ndx $d]
            # puts "$i_ndx $j_ndx"
        }
    }
    set mydist [lsort -real -index 4 $mydist]
    return $mydist    
}

#once $grp1 $grp2

set grp1 "RING and resid 1"
set grp2 "CARB"
set grp1 "resid 1"
set grp2 "resid 603"
set mydist [min_dist $grp1 $grp2]
puts [lindex $mydist 0]
puts [lindex $mydist end]


