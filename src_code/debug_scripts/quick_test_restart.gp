tcol = 2
KEcol = 6
PEcol = 7

set terminal wxt 0 size 1400,800
set datafile separator ";"
set multiplot layout 2,2 columnsfirst title "Kinetic Energies vs. Time (Restarts)"

set title "Parallel to parallel"
set xlabel "t(LJU)"
set ylabel "KE"
plot [2:] "./start_p.macro" u (column(tcol)):(column(KEcol)) w p pt 4 lc 3 not,\
          "./restart_p2p.macro" u (column(tcol)):(column(KEcol)) w p pt 9 lc 3 not,\
          "./full_p.macro" u (column(tcol)):(column(KEcol)) w l lt 1 lc 5 not 
  
set title "Serial to serial"
set xlabel "t(LJU)"
plot [2:] "./start_s.macro" u (column(tcol)):(column(KEcol)) w p pt 4 lc 3 not,\
          "./restart_s2s.macro" u (column(tcol)):(column(KEcol)) w p pt 9 lc 3 not,\
          "./full_s.macro" u (column(tcol)):(column(KEcol)) w l lt 1 lc 5 not

set title "Serial to parallel"
set xlabel "t(LJU)"
plot [2:] "./start_s.macro" u (column(tcol)):(column(KEcol)) w p pt 4 lc 3 not,\
          "./restart_s2p.macro" u (column(tcol)):(column(KEcol)) w p pt 9 lc 3 not,\
          "./full_s.macro" u (column(tcol)):(column(KEcol)) w l lt 1 lc 5 not
  
set title "Parallel to serial"
set xlabel "t(LJU)"
plot [2:] "./start_p.macro" u (column(tcol)):(column(KEcol)) w p pt 4 lc 3 not,\
          "./restart_p2s.macro" u (column(tcol)):(column(KEcol)) w p pt 9 lc 3 not,\
          "./full_p.macro" u (column(tcol)):(column(KEcol)) w l lt 1 lc 5 not

unset multiplot

set terminal wxt 1 size 1400,800
set datafile separator ";"
set multiplot layout 2,2 columnsfirst title "Potential Energies vs. Time (Restarts)"
  
set title "Parallel to parallel"
set xlabel "t(LJU)"
set ylabel "PE"
plot [2:] "./start_p.macro" u (column(tcol)):(column(PEcol)) w p pt 4 lc 1 not,\
          "./restart_p2p.macro" u (column(tcol)):(column(PEcol)) w p pt 9 lc 1 not,\
          "./full_p.macro" u (column(tcol)):(column(PEcol)) w l lt 1 lc 8 not 

set title "Serial to serial"
set xlabel "t(LJU)"
plot [2:] "./start_s.macro" u (column(tcol)):(column(PEcol)) w p pt 4 lc 1 not,\
          "./restart_s2s.macro" u (column(tcol)):(column(PEcol)) w p pt 9 lc 1 not,\
          "./full_s.macro" u (column(tcol)):(column(PEcol)) w l lt 1 lc 8 not

set title "Serial to parallel"
set xlabel "t(LJU)"
plot [2:] "./start_s.macro" u (column(tcol)):(column(PEcol)) w p pt 4 lc 1 not,\
          "./restart_s2p.macro" u (column(tcol)):(column(PEcol)) w p pt 9 lc 1 not,\
          "./full_s.macro" u (column(tcol)):(column(PEcol)) w l lt 1 lc 8 not

set title "Parallel to serial"
set xlabel "t(LJU)"
plot [2:] "./start_p.macro" u (column(tcol)):(column(PEcol)) w p pt 4 lc 1 not,\
          "./restart_p2s.macro" u (column(tcol)):(column(PEcol)) w p pt 9 lc 1 not,\
          "./full_p.macro" u (column(tcol)):(column(PEcol)) w l lt 1 lc 8 not

