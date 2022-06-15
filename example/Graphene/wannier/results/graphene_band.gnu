set style data dots
set nokey
set xrange [0: 4.01729]
set yrange [-10.09857 : 10.68117]
set arrow from  1.47043, -10.09857 to  1.47043,  10.68117 nohead
set arrow from  2.31938, -10.09857 to  2.31938,  10.68117 nohead
set xtics ("G"  0.00000,"M"  1.47043,"K"  2.31938,"G"  4.01729)
 plot "graphene_band.dat"
