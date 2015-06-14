#! tcl-vmd

#
# need XYdist Ytype 
#
#
# this atom X list
# select elements
set sel [atomselect top "element Br and (altloc '' or altloc A)"]
set xatom_list [$sel get index]
set res_name [lsort -unique [$sel get resname]]
# $sel delete

puts "final selected residues <resname>: $res_name"


set tmplist {}
set minResAtom 5
# that atom A
foreach xatom $xatom_list {
    puts "index X $xatom"
    # this atom X
    set xsel [atomselect top "index $xatom"]
    set xres [$xsel get resname]
    # puts "xsel: $xsel"
    # that atom A
    set aatom [$xsel getbonds]
    set naatom [llength [lindex $aatom 0]]
    if { $naatom != 1 } {
        puts "too many connected atoms ($naatom) for $aatom"
    } else {
        set rsel [atomselect top "resname $xres"]
        set natom [ $rsel num]
        if { $natom < $minResAtom } {
            puts "$xres $natom"
            continue
        }
        $rsel delete
        puts "index A-X $aatom $xatom"
        set dist [measure bond [list $xatom $aatom]]
        puts $dist
        set halo [list $xres $aatom $xatom $dist]
        lappend tmplist $halo
    }
    $xsel delete
}

# puts $tmplist

# =============
# CL-- O N S
# BR
# I
#
set xy_list [list [list 3.27 3.30 3.55] [list 3.40 3.37 3.65] [list 3.50 3.53 3.78]]
set xpi_list [list [list 4.2 60.0 120.0] [list 4.3 60.0 120.0] [list 4.5 60.0 120.0]]
# Phe, Tyr, His, and Trp
# =============


set halo_list {}
if { ![info exists $Ytype] } {
    puts "Ytype is $Ytype"
} else {
    set Ytype "O"
}
if { ![info exists $XYdist] } {
    puts "XYdist is $XYdist"
} else {
    set XYdist 3.37
}


foreach halo $tmplist {
    set xres [lindex $halo 0]
    set aatom [lindex $halo 1]
    set xatom [lindex $halo 2]
    set length [lindex $halo 3]
    # that atom Y
    # there may be more than one Y atom for a A-X case..
    set ysel [atomselect top "not resname $xres and element $Ytype and within $XYdist of index $xatom"]
    # set ysel [atomselect top "not resname $xres and within $XYdist of index $xatom"]
    # puts "ysel: $ysel"
    set yatom_list [$ysel get index]
    set nyatom [$ysel num]
    if { $n_yatom > 0 } {
        foreach yatom $yatom_list {
            set ysel2 [atomselect top "index $yatom"]
            set yres [$ysel2 get resname]
            set ytype [$ysel2 get name]
            set dist [measure bond [list $xatom $yatom]]
            set angle [measure angle [list $aatom $xatom $yatom]]
            set triatoms [list $xres $yres $aatom $xatom $yatom $length $dist $angle]
            lappend halo_list $triatoms
            puts "index A-X -- Y: $aatom $xatom $yatom ($ytype)"
        } 
    } else {
        puts "NO satisfied Y atom.."
    }
}
$sel delete

puts $halo_list

