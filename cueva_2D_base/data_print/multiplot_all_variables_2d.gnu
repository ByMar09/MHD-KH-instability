reset
set terminal pngcairo  enhanced color size 1500, 2300	
set output '00_2d_test.png' 
set pm3d map
set encoding utf8
set view map
set border 4095 front linetype -1 linewidth 1.000
#unset surface
unset key
set style data pm3d
set style function pm3d
set ticslevel 0
#set pm3d implicit at st
set palette rgb 33,13,10
#set pal gray


######## README ###########


# function to margin and gap between rows and columns (Christoph)


init_margins(left, right, bottom, top, dx, dy, rows, cols) = \
  sprintf('left_margin = %f; right_margin = %f; top_margin = %f; bottom_margin = %f; ', left, right, top, bottom) . \
  sprintf('col_count = %d; row_count = %d; gap_size_x = %f; gap_size_y = %f', cols, rows, dx, dy)

get_lmargin(col) = (left_margin + (col - 1) * (gap_size_x + ((right_margin - left_margin)-(col_count - 1) * gap_size_x)/col_count))
get_rmargin(col) = (left_margin + (col - 1) * gap_size_x + col * ((right_margin - left_margin)-(col_count - 1) * gap_size_x)/col_count)
get_tmargin(row) = (top_margin  - (row - 1) * gap_size_y - (row-1) * ((top_margin - bottom_margin  - gap_size_y * row_count) / row_count))
get_bmargin(row) = (top_margin  - (row - 1) * gap_size_y -  row    * ((top_margin - bottom_margin  - gap_size_y * row_count) / row_count))
set_margins(col, row) = \
  sprintf('set lmargin at screen %f;', get_lmargin(col)) . \
  sprintf('set rmargin at screen %f;', get_rmargin(col)) . \
  sprintf('set tmargin at screen %f;', get_tmargin(row)) . \
  sprintf('set bmargin at screen %f;', get_bmargin(row))  


######## README ###########

### Start multiplot (5x3 layout)
set multiplot layout 5,3 rowsfirst 

eval(init_margins(0.05, 0.925, 0.05, 0.95, 0.1, 0.005,5, 3))


unset xlabel
#set xrange [ -6.0 : 6.0 ] noreverse nowriteback
#set xtics    -6.0, 3.0, 6.0
#set yrange [ -6.0 : 6.0 ] noreverse nowriteback
#set ytics    -6.0, 3.0, 6.0
#set zlabel "P" font "arial,10"
#set zrange [ -6.000 : 6.000  ] noreverse nowriteback 
#set cblabel ""  font "arial,10"
set colorbox vertical

set format cb "%.0e"
set format x "%g"
set format y "%g"


# --- GRAPH (1,1)

eval(set_margins(1,1))

#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced

set format x ""

set label 1 'Bx' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($31)


# --- GRAPH (1,2)

eval(set_margins(2,1))


#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

unset xlabel
unset ylabel

set format y ""

set label 1 'By ' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($32)


# --- GRAPH (1,3)


eval(set_margins(3,1))

#set cbrange [  0 : 1  ]
#set cbtics     0, 0.5 , 1

unset title
unset xlabel
unset ylabel 

set format y ""

set label 1 'Bz' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($33)


# --- GRAPH (2,1)

eval(set_margins(1,2))

unset xlabel
set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced


set format y "%g"


#set cbrange [  -2.5e-5 : 2.5e-5  ]
#set cbtics     -2.e-5, 1e-5 , 2.e-5

#set format cb "%.1t{/Symbol \264}10^{%L}" 
#set cbtics add ("0" 0) 

set label 1 'Vx' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($34)

# --- GRAPH (2,2)

eval(set_margins(2,2))

#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

unset title
unset ylabel 

set format y ""

set label 1 'Vy' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($35)

# --- GRAPH (2,3)

eval(set_margins(3,2))

#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

unset ylabel

set format y ""

set label 1 'Vz' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($36)

# --- GRAPH (3,1)


eval(set_margins(1,3))

#set cbrange [  0.49 : 0.57  ]
#set cbtics     0.49, 0.02 , 0.57

set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced

set format y "%g"

set label 1 'Ex' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($37)

# --- GRAPH (3,2)

eval(set_margins(2,3))

#set cbrange [  1 : 1.014  ]
#set cbtics     1, 0.002 , 1.014

unset ylabel

set format y ""

set label 1 'Ey' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($38)


# --- GRAPH (3,3)

eval(set_margins(3,3))


#set cbrange [  0.5 : 0.514  ]
#set cbtics     0.5, 0.002 , 0.514

unset ylabel

set format y ""

set label 1 'Ez' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($39)

# --- GRAPH (4,1)

eval(set_margins(1,4))

#set cbrange [  0 : 120  ]
#set cbtics     0, 30 , 120

set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced

set format y "%g"

set label 1 'Q' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($40)

# --- GRAPH (4,2)


eval(set_margins(2,4))

#set format cb "%.3f"

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

unset ylabel

set format y ""

set label 1 'Rho ' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($41)

# --- GRAPH (4,3)

eval(set_margins(3,4))

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

#set format x "%g"
#set xtics -0.2,0.1,0.2

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel

set format y ""
set format x ""


set label 1 'P ' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($42)

# --- GRAPH (5,1)


eval(set_margins(1,5))

#set format cb "%.0e"

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel

set format y "%g"
set format x "%g"

set label 1 'psi-->(divE) ' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($29)

# --- GRAPH (5,2)

eval(set_margins(2,5))

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel

set format y ""

set label 1 'phi-->(divB) ' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):($30)

# --- GRAPH (5,3)

eval(set_margins(3,5))

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel

set format y ""

set label 1 'V^2 ' at graph 0.05,0.85 font ',14' front
splot "file3.dat" u ($27):($26):(($34**2 + $35**2 + $36**2))


unset multiplot
### End multiplot

