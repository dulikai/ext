set datadir ..
mol new "$datadir/md.gro"
mol addfile "$datadir/wrap.xtc" waitfor -1

