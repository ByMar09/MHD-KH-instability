reset

set terminal pngcairo  enhanced color font "arial,12" fontscale 1.0 size 910, 910 
set output '00_cut_Bx_test.png' 
#unset key
#set ticslevel 0
set title "Field Bx in y=0  S_l = 10^6 k = 12 (256X32) HLLC" font "Times-Roman,20"  tc rgb "gray"

#Sl = 1e6
set xrange [ -0.2 : 0.2 ] noreverse nowriteback
set xtics    -0.2,0.1,0.2

#Sl = 1e5
#set xrange [ -0.430887 : 0.430887 ] noreverse nowriteback
#set xtics    -0.4,0.2,0.4

#Sl = 1e4
#set xrange [ -0.928318 : 0.928318 ] noreverse nowriteback
#set xtics    -0.9,0.3,0.9


set xlabel "x" font "Times-Roman,20" tc rgb "gray" rotate by 00 offset 0,0,0
set ylabel "Bx" font "Times-Roman,22" tc rgb "gray" rotate by 00 offset 0,0,0
#set yrange [0.8:1.2] noreverse nowriteback
#set ytics  0.8, 0.2, 1.2 font "Times-Roman,18"
set key top right

set lmargin 12
set tmargin 4
set bmargin 4

#set pointintervalbox 3

#set style line 1 lc rgb '#0060ad' lt 1 lw  pi -1 pt 7 ps 1.0   # --- blue
#set style line 1 lt 1 lw 1 pt 3 ps 0.5


plot  'TM_1D_Bx_t_000.dat'   title 't = 00.0' w l ls 1 linecolor rgb "dark-green"  lw 1,\
      'TM_1D_Bx_t_025.dat'   title 't = 02.5' w l ls 1 linecolor rgb "red"         lw 1,\
      'TM_1D_Bx_t_050.dat'   title 't = 05.0' w l ls 1 linecolor rgb "orange"      lw 1,\
      'TM_1D_Bx_t_075.dat'   title 't = 07.5' w l ls 1 linecolor rgb "violet"      lw 1,\
      'TM_1D_Bx_t_100.dat'   title 't = 10.0' w l ls 1 linecolor rgb "light-blue"  lw 1,\

