mol load xyz geom_9_7.xyz
set num [molinfo top get numframes]
mol drawframes top 0 0:$num
color Display Background white
render snapshot a.gif
