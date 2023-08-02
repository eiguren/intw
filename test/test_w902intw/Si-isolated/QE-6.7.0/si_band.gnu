set style data dots
set nokey
set xrange [0: 3.81716]
set yrange [ -6.82388 :  7.22059]
set arrow from  1.00778,  -6.82388 to  1.00778,   7.22059 nohead
set arrow from  2.17146,  -6.82388 to  2.17146,   7.22059 nohead
set arrow from  2.58288,  -6.82388 to  2.58288,   7.22059 nohead
set xtics ("L"  0.00000,"G"  1.00778,"X"  2.17146,"K"  2.58288,"G"  3.81716)
 plot "si_band.dat"
