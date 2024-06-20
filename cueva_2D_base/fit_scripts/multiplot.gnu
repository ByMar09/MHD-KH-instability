reset
#set terminal postscript eps enhanced color solid font "Times-Roman,15" 
#set output '2d_multiplot.eps' 
#set size 1.5,1.0
set terminal pngcairo  enhanced color font "arial,12" fontscale 1.0 size 1400, 1900	
set output '00_2d_multiplot_test.png' 
set pm3d map
set encoding utf8
set view map
set border 4095 front linetype -1 linewidth 1.000
#set border 10 front ls 11
#set tmargin at screen 0.95
#set bmargin at screen 0.05
#set rmargin at screen 0.05
#set lmargin 00
#unset surface
unset key
set style data pm3d
set style function pm3d
set ticslevel 0
#set pm3d implicit at st
set palette rgb 33,13,10
#set pal gray

set macros

# Margins for each row resp. column

LMARGIN = "set lmargin at screen 0.70; set rmargin at screen 0.90"

#TMARGIN = "set tmargin at screen 0.15; set bmargin at screen 0.55"
#BMARGIN = "set tmargin at screen 0.25; set bmargin at screen 0.10"
#LMARGIN = "set lmargin at screen 0.10; set rmargin at screen 0.30"
#RMARGIN = "set lmargin at screen 0.5; set rmargin at screen 0.90"


### Start multiplot (3x3 layout)
set multiplot layout 4,3 rowsfirst title "Tearing Mode t = 2.5  S_l= 10^6 k= 14 {/Symbol \153} = 1.0 {/Symbol \264} 10^4 (512X64) HLLC"  font "arial,20"

# OJO estilo de linea en: (eps lw  0.005) (png lw  0.1)

#Sl = 1e6 k=02
#set xrange [ -0.2 : 0.2 ] noreverse nowriteback
#set xtics    -0.2,0.1,0.2
#set yrange [  0.0 : 3.1416] noreverse nowriteback
#set ytics     0.0,1.0,3.0

#Sl = 1e6 k=04
#set xrange [ -0.2 : 0.2 ] noreverse nowriteback
#set xtics    -0.2,0.1,0.2
#set yrange [  0.0 : 1.5708] noreverse nowriteback
#set ytics     0.0,0.5,1.5

#Sl = 1e6 k=06
#set xrange [ -0.2 : 0.2 ] noreverse nowriteback
#set xtics    -0.2,0.1,0.2
#set yrange [  0.0 : 1.0472] noreverse nowriteback
#set ytics     0.0,0.5,1.5

#Sl = 1e6 k=08
#set xrange [ -0.2 : 0.2 ] noreverse nowriteback
#set xtics    -0.2,0.1,0.2
#set yrange [  0.0 : 0.785398] noreverse nowriteback
#set ytics     0.0,0.3,1.5

#Sl = 1e6 k=10
#set xrange [ -0.2 : 0.2 ] noreverse nowriteback
#set xtics    -0.2,0.1,0.2
#set yrange [  0.0 : 0.628319] noreverse nowriteback
#set ytics     0.0,0.2,0.6

#Sl = 1e6 k=12
set xrange [ -0.2 : 0.2 ] noreverse nowriteback
set xtics    -0.2,0.1,0.2
set yrange [  0.0 : 0.523599] noreverse nowriteback
set ytics     0.0,0.25,0.5

#Sl = 1e6 k=14
#set xrange [ -0.2 : 0.2 ] noreverse nowriteback
#set xtics    -0.2,0.1,0.2
#set yrange [  0.0 : 0.448799] noreverse nowriteback
~set ytics     0.0,0.2,0.4

#Sl = 1e6 k=16
#set xrange [ -0.2 : 0.2 ] noreverse nowriteback
#set xtics    -0.2,0.1,0.2
#set yrange [  0.0 : 0.392699] noreverse nowriteback
#set ytics     0.0,0.15,0.3





#Sl = 1e5
#set xrange [ -0.430887 : 0.430887 ] noreverse nowriteback
#set xtics    -0.4,0.2,0.4
#set yrange [  0.0 : 0.658747 ] noreverse nowriteback
#set ytics     0.0,0.3,0.6

#Sl = 1e4
#set xrange [ -0.928318 : 0.928318 ] noreverse nowriteback
#set xtics    -0.9,0.3,0.9
#set yrange [  0.0 : 0.966908 ] noreverse nowriteback
#set ytics     0.0,0.3,0.9

# set zrange [ -6e-4 : 6e-4  ] noreverse nowriteback 
# set ztics    -6.e-4, 2e-4, 6e-4
# set colorbox vertical
# set cbtics autofreq  norangelimit font ",10"
#  set cbrange [  -2.5e-3 : 2.5e-3  ]
#  set cbtics     -2.e-3, 1e-3 , 2.e-3

# estilo de linea contornos  
#set style line 1 lc rgb "#000000" lw 0.2
#set style increment userstyle
#unset clabel
#contornos
#set contour base 
# numero de contornos
#set cntrparam levels incremental -3.0, 0.5, 6.0

# Logaritmic in data
# u 1:2:(log($3))

# Log |E dot B |
#Log  {/Symbol \275} E {/Symbol \267} B {/Symbol \275} 

# --- GRAPH (1,1)

set samples 200
set isosamples 200

#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced
set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced


set label 1 'Bx' at graph 0.05,0.85 font ',14' front
splot "TM_2D_Bx_t_100.dat" u 1:2:3


# --- GRAPH (1,2)

set samples 200
set isosamples 200

#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel

set label 1 'By ' at graph 0.05,0.85 font ',14' front
splot "TM_2D_By_t_100.dat" u 1:2:3


# --- GRAPH (1,3)

set samples 200
set isosamples 200

#set cbrange [  0 : 1  ]
#set cbtics     0, 0.5 , 1

unset title
set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 


set label 1 'Bz' at graph 0.05,0.85 font ',14' front
splot "TM_2D_Bz_t_100.dat" u 1:2:3


# --- GRAPH (2,1)

#@TMARGIN 

set samples 200
set isosamples 200

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced


#set cbrange [  -2.5e-5 : 2.5e-5  ]
#set cbtics     -2.e-5, 1e-5 , 2.e-5

#set format cb "%.1t{/Symbol \264}10^{%L}" 
#set cbtics add ("0" 0) 

set label 1 'Vx' at graph 0.05,0.85 font ',14' front
splot "TM_2D_Vx_t_100.dat" u 1:2:3

# --- GRAPH (2,2)

unset format

set samples 200
set isosamples 200

#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

unset title
set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 



set label 1 'Vy' at graph 0.05,0.85 font ',14' front
splot "TM_2D_Vy_t_100.dat" u 1:2:3

# --- GRAPH (2,3)

#@LMARGIN 

set samples 200
set isosamples 200

#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 

set label 1 'Vz' at graph 0.05,0.85 font ',14' front
splot "TM_2D_Vz_t_100.dat" u 1:2:3

# --- GRAPH (3,1)

#@BMARGIN 


set samples 200
set isosamples 200

#set cbrange [  0.49 : 0.57  ]
#set cbtics     0.49, 0.02 , 0.57

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced

set label 1 'P' at graph 0.05,0.85 font ',14' front
splot "TM_2D_P_t_100.dat" u 1:2:3

# --- GRAPH (3,2)


set samples 200
set isosamples 200

#set cbrange [  1 : 1.014  ]
#set cbtics     1, 0.002 , 1.014

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 

set label 1 'Rho' at graph 0.05,0.85 font ',14' front
splot "TM_2D_Rho_t_100.dat" u 1:2:3


# --- GRAPH (3,3)


set samples 200
set isosamples 200

#set cbrange [  0.5 : 0.514  ]
#set cbtics     0.5, 0.002 , 0.514

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 

set label 1 'P/Rho' at graph 0.05,0.85 font ',14' front
splot "TM_2D_p-rho_t_100.dat" u 1:2:3

# --- GRAPH (4,1)

#@BMARGIN 


set samples 200
set isosamples 200

#set cbrange [  0 : 120  ]
#set cbtics     0, 30 , 120

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced

set label 1 'Jz' at graph 0.05,0.85 font ',14' front
splot "TM_2D_Jz_t_100.dat" u 1:2:3

# --- GRAPH (4,2)

set samples 200
set isosamples 200

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

set format cb "%.1t{/Symbol \264}10^{%L}" 
set cbtics add ("0" 0) 

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 

set label 1 'E {/Symbol \267} B ' at graph 0.05,0.85 font ',14' front
splot "TM_2D_EdotB_t_100.dat" u 1:2:3


unset multiplot
### End multiplot

