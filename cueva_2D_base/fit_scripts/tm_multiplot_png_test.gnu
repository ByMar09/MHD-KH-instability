reset
set terminal pngcairo  color enhanced font "arial,12" fontscale 1.0 size 1600, 500
set output '00_profiles_tm_test.png' 

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


set termoption dashed

unset colorbox
#
# reset linetypes to base dash patterns
#
#set for [i=1:5] linetype i dt i

set linestyle 1  lt 1  lc rgb "red"           ps 0.5  lw 1 
set linestyle 2  lt 1  lc rgb "dark-orange"     ps 0.5  lw 1
set linestyle 3  lt 1  lc rgb "orange"          ps 0.5  lw 1 
set linestyle 4  lt 1  lc rgb "green"         ps 0.5  lw 1
set linestyle 5  lt 1  lc rgb "dark-green"        ps 0.5  lw 1
set linestyle 6  lt 1  lc rgb "dark-cyan"       ps 0.5  lw 1
set linestyle 7  lt 1  lc rgb "cyan"    ps 0.5  lw 1
set linestyle 8  lt 1  lc rgb "blue"         ps 0.5  lw 1
set linestyle 9  lt 1  lc rgb "dark-blue"          ps 0.5  lw 1 

set style line 020 lt 0 lc rgb "gray" lw 0.5

show style line

### Start multiplot (1x3 layout)
set multiplot layout 1,3 rowsfirst 

eval(init_margins(0.08, 0.94, 0.1, 0.95, 0.05, 0.05,1, 3))  

set grid back  

set mxtics 5
set mytics 5

set format  y '%g'
set format  x '%g'

### Statistic Bx data
stats '1d_data_cut_y_bx_t_1.dat' u 1:2 nooutput
max_bx_t_1 = STATS_max_y
stats '1d_data_cut_y_bx_t_2.dat' u 1:2 nooutput
max_bx_t_2 = STATS_max_y
stats '1d_data_cut_y_bx_t_3.dat' u 1:2 nooutput
max_bx_t_3 = STATS_max_y
stats '1d_data_cut_y_bx_t_4.dat' u 1:2 nooutput
max_bx_t_4 = STATS_max_y
stats '1d_data_cut_y_bx_t_5.dat' u 1:2 nooutput
max_bx_t_5 = STATS_max_y
stats '1d_data_cut_y_bx_t_6.dat' u 1:2 nooutput
max_bx_t_6 = STATS_max_y
stats '1d_data_cut_y_bx_t_7.dat' u 1:2 nooutput
max_bx_t_7 = STATS_max_y
stats '1d_data_cut_y_bx_t_8.dat' u 1:2 nooutput
max_bx_t_8 = STATS_max_y
stats '1d_data_cut_y_bx_t_9.dat' u 1:2 nooutput
max_bx_t_9 = STATS_max_y
### Statistic Vx data
stats '1d_data_cut_y_vx_t_1.dat' u 1:2 nooutput
max_vx_t_1 = STATS_max_y
stats '1d_data_cut_y_vx_t_2.dat' u 1:2 nooutput
max_vx_t_2 = STATS_max_y
stats '1d_data_cut_y_vx_t_3.dat' u 1:2 nooutput
max_vx_t_3 = STATS_max_y
stats '1d_data_cut_y_vx_t_4.dat' u 1:2 nooutput
max_vx_t_4 = STATS_max_y
stats '1d_data_cut_y_vx_t_5.dat' u 1:2 nooutput
max_vx_t_5 = STATS_max_y
stats '1d_data_cut_y_vx_t_6.dat' u 1:2 nooutput
max_vx_t_6 = STATS_max_y
stats '1d_data_cut_y_vx_t_7.dat' u 1:2 nooutput
max_vx_t_7 = STATS_max_y
stats '1d_data_cut_y_vx_t_8.dat' u 1:2 nooutput
max_vx_t_8 = STATS_max_y
stats '1d_data_cut_y_vx_t_9.dat' u 1:2 nooutput
max_vx_t_9 = STATS_max_y
### Statistic Vy data
stats '1d_data_cut_y_vy_t_1.dat' u 1:2 nooutput
max_vy_t_1 = STATS_max_y
stats '1d_data_cut_y_vy_t_2.dat' u 1:2 nooutput
max_vy_t_2 = STATS_max_y
stats '1d_data_cut_y_vy_t_3.dat' u 1:2 nooutput
max_vy_t_3 = STATS_max_y
stats '1d_data_cut_y_vy_t_4.dat' u 1:2 nooutput
max_vy_t_4 = STATS_max_y
stats '1d_data_cut_y_vy_t_5.dat' u 1:2 nooutput
max_vy_t_5 = STATS_max_y
stats '1d_data_cut_y_vy_t_6.dat' u 1:2 nooutput
max_vy_t_6 = STATS_max_y
stats '1d_data_cut_y_vy_t_7.dat' u 1:2 nooutput
max_vy_t_7 = STATS_max_y
stats '1d_data_cut_y_vy_t_8.dat' u 1:2 nooutput
max_vy_t_8 = STATS_max_y
stats '1d_data_cut_y_vy_t_9.dat' u 1:2 nooutput
max_vy_t_9 = STATS_max_y

# --- GRAPH (1,1)
eval(set_margins(1,1))

set ytics 0.0,0.2,1
set xtics -0.2,0.1,0.2
set xrange[-0.2:0.2]
set yrange[0:1]

set ytics add ('0.0' 0)

set xlabel 'x'  

unset key

set label 1 'B_x' at graph 0.775,0.85 front font "arial,22"

plot   '1d_data_cut_y_bx_t_1.dat'  u ($1):($2/ abs(max_bx_t_1 ))  title 't=1'    w l ls 1 ,\
       '1d_data_cut_y_bx_t_2.dat'  u ($1):($2/ abs(max_bx_t_2 ))  title 't=2'    w l ls 2 ,\
       '1d_data_cut_y_bx_t_3.dat'  u ($1):($2/ abs(max_bx_t_3 ))  title 't=3'    w l ls 3 ,\
       '1d_data_cut_y_bx_t_4.dat'  u ($1):($2/ abs(max_bx_t_4 ))  title 't=4'    w l ls 4 ,\
       '1d_data_cut_y_bx_t_5.dat'  u ($1):($2/ abs(max_bx_t_5 ))  title 't=5'    w l ls 5 ,\
       '1d_data_cut_y_bx_t_6.dat'  u ($1):($2/ abs(max_bx_t_6 ))  title 't=6'    w l ls 6 ,\
       '1d_data_cut_y_bx_t_7.dat'  u ($1):($2/ abs(max_bx_t_7 ))  title 't=7'    w l ls 7 ,\
       '1d_data_cut_y_bx_t_8.dat'  u ($1):($2/ abs(max_bx_t_8 ))  title 't=8'    w l ls 8 ,\
       '1d_data_cut_y_bx_t_9.dat'  u ($1):($2/ abs(max_bx_t_9 ))  title 't=9'    w l ls 9
       
       
# --- GRAPH (1,2)

set style data linespoints

eval(set_margins(2,1))

set ytics -1,0.5,1
#set xtics 0.4,0.2,0.8
#set xrange[0.4:0.8]
set yrange[-1:1]
#set ytics add ('0.0' 0)

set key bottom left

set label 1 'v_x' at graph 0.775,0.85 front font "arial,22"

plot   '1d_data_cut_y_vx_t_1.dat'  u ($1):($2/ abs(max_vx_t_1 ))  title 't=1'    w l ls 1 ,\
       '1d_data_cut_y_vx_t_2.dat'  u ($1):($2/ abs(max_vx_t_2 ))  title 't=2'    w l ls 2 ,\
       '1d_data_cut_y_vx_t_3.dat'  u ($1):($2/ abs(max_vx_t_3 ))  title 't=3'    w l ls 3 ,\
       '1d_data_cut_y_vx_t_4.dat'  u ($1):($2/ abs(max_vx_t_4 ))  title 't=4'    w l ls 4 ,\
       '1d_data_cut_y_vx_t_5.dat'  u ($1):($2/ abs(max_vx_t_5 ))  title 't=5'    w l ls 5 ,\
       '1d_data_cut_y_vx_t_6.dat'  u ($1):($2/ abs(max_vx_t_6 ))  title 't=6'    w l ls 6 ,\
       '1d_data_cut_y_vx_t_7.dat'  u ($1):($2/ abs(max_vx_t_7 ))  title 't=7'    w l ls 7 ,\
       '1d_data_cut_y_vx_t_8.dat'  u ($1):($2/ abs(max_vx_t_8 ))  title 't=8'    w l ls 8 ,\
       '1d_data_cut_y_vx_t_9.dat'  u ($1):($2/ abs(max_vx_t_9 ))  title 't=9'    w l ls 9 

# --- GRAPH (1,3)

set style data linespoints

eval(set_margins(3,1))

set ytics -0.4,0.2,1
#set xtics 0.4,0.1,0.8
#set xrange[0.4:0.8]
set yrange[-0.4:1]
#set ytics add ('0.0' 0)

unset key

set label 1 'v_y' at graph 0.775,0.85 front font "arial,22"

plot   '1d_data_cut_y_vy_t_1.dat'  u ($1):($2/ abs(max_vy_t_1 ))  title 't=1'    w l ls 1 ,\
       '1d_data_cut_y_vy_t_2.dat'  u ($1):($2/ abs(max_vy_t_2 ))  title 't=2'    w l ls 2 ,\
       '1d_data_cut_y_vy_t_3.dat'  u ($1):($2/ abs(max_vy_t_3 ))  title 't=3'    w l ls 3 ,\
       '1d_data_cut_y_vy_t_4.dat'  u ($1):($2/ abs(max_vy_t_4 ))  title 't=4'    w l ls 4 ,\
       '1d_data_cut_y_vy_t_5.dat'  u ($1):($2/ abs(max_vy_t_5 ))  title 't=5'    w l ls 5 ,\
       '1d_data_cut_y_vy_t_6.dat'  u ($1):($2/ abs(max_vy_t_6 ))  title 't=6'    w l ls 6 ,\
       '1d_data_cut_y_vy_t_7.dat'  u ($1):($2/ abs(max_vy_t_7 ))  title 't=7'    w l ls 7 ,\
       '1d_data_cut_y_vy_t_8.dat'  u ($1):($2/ abs(max_vy_t_8 ))  title 't=8'    w l ls 8 ,\
       '1d_data_cut_y_vy_t_9.dat'  u ($1):($2/ abs(max_vy_t_9 ))  title 't=9'    w l ls 9 

unset multiplot
### End multiplot

