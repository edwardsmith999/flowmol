set term wxt 1
set palette gray
splot "info_olap_md" u 3:4:6 w p ps 5 pt 5 palette notitle

set term wxt 2
set view map
set palette rgbformulae 33,13,10
splot [][][0:] "info_olap_md" u 3:4:6 w p ps 2 pt 5 palette notitle

