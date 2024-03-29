 set terminal postscript "fontsize" 26
 set output "com_lda_bnd_LAPW.eps"
set data style lines
set nokey
set xrange [0: 0.77897]
set yrange [-204.08025 : 204.08025]
set arrow from  0.12525,-204.08025 to  0.12525, 204.08025 nohead
set arrow from  0.21381,-204.08025 to  0.21381, 204.08025 nohead
set arrow from  0.30238,-204.08025 to  0.30238, 204.08025 nohead
set arrow from  0.41085,-204.08025 to  0.41085, 204.08025 nohead
set arrow from  0.51932,-204.08025 to  0.51932, 204.08025 nohead
set arrow from  0.62778,-204.08025 to  0.62778, 204.08025 nohead
set arrow from  0.69041,-204.08025 to  0.69041, 204.08025 nohead
set xtics (" G "  0.00000," H "  0.12525," N "  0.21381," G "  0.30238," P "  0.41085," H "  0.51932," P "  0.62778," N "  0.69041,"   "  0.77897)
 plot "com_dft_band_LAPW.dat" using 2:3 with lines lw 3 lt 1 lc 1
 unset output
 unset terminal
