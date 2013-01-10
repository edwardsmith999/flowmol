tcol = 2
KEcol = 6
PEcol = 9

set terminal wxt size 1200,600
set datafile separator ";"
set multiplot layout 2,2 columnsfirst title "Energies vs. Time (Restarts)"
set key at screen 0.54,0.95
set key horizontal

set title "Parallel to parallel"
set xlabel "t(LJU)"
plot "./100_p.macro" u (column(tcol)):(column(KEcol)) w p pt 4 lc 3 not,\
     "./100_p.macro" u (column(tcol)):(column(PEcol)) w p pt 4 lc 1 not,\
     "./100-200_p2p.macro" u (column(tcol)):(column(KEcol)) w p pt 9 lc 3 not,\
     "./100-200_p2p.macro" u (column(tcol)):(column(PEcol)) w p pt 9 lc 1 not,\
     "./200_p.macro" u (column(tcol)):(column(KEcol)) w l lt 1 lc 5 t "KE",\
     "./200_p.macro" u (column(tcol)):(column(PEcol)) w l lt 1 lc 8 t "PE"
  
set title "Serial to serial"
set xlabel "t(LJU)"
plot "./100_s.macro" u (column(tcol)):(column(KEcol)) w p pt 4 lc 3 not,\
     "./100_s.macro" u (column(tcol)):(column(PEcol)) w p pt 4 lc 1 not,\
     "./100-200_s2s.macro" u (column(tcol)):(column(KEcol)) w p pt 9 lc 3 not,\
     "./100-200_s2s.macro" u (column(tcol)):(column(PEcol)) w p pt 9 lc 1 not,\
     "./200_s.macro" u (column(tcol)):(column(KEcol)) w l lt 1 lc 5 not,\
     "./200_s.macro" u (column(tcol)):(column(PEcol)) w l lt 1 lc 8 not
  
set title "Serial to parallel"
set xlabel "t(LJU)"
plot "./100_s.macro" u (column(tcol)):(column(KEcol)) w p pt 4 lc 3 not,\
     "./100_s.macro" u (column(tcol)):(column(PEcol)) w p pt 4 lc 1 not,\
     "./100-200_s2p.macro" u (column(tcol)):(column(KEcol)) w p pt 9 lc 3 not,\
     "./100-200_s2p.macro" u (column(tcol)):(column(PEcol)) w p pt 9 lc 1 not,\
     "./200_s.macro" u (column(tcol)):(column(KEcol)) w l lt 1 lc 5 not,\
     "./200_s.macro" u (column(tcol)):(column(PEcol)) w l lt 1 lc 8 not
  
set title "Parallel to serial"
set xlabel "t(LJU)"
plot "./100_p.macro" u (column(tcol)):(column(KEcol)) w p pt 4 lc 3 not,\
     "./100_p.macro" u (column(tcol)):(column(PEcol)) w p pt 4 lc 1 not,\
     "./100-200_p2s.macro" u (column(tcol)):(column(KEcol)) w p pt 9 lc 3 not,\
     "./100-200_p2s.macro" u (column(tcol)):(column(PEcol)) w p pt 9 lc 1 not,\
     "./200_p.macro" u (column(tcol)):(column(KEcol)) w l lt 1 lc 5 not,\
     "./200_p.macro" u (column(tcol)):(column(PEcol)) w l lt 1 lc 8 not
