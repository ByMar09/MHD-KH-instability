reset
set terminal epslatex  color size 11, 7  header  "\\newcommand{\\lr}[0]{\\Large}\n\\newcommand{\\ft}[0]{\\footnotesize}\n\\newcommand{\\hg}[0]{\\Huge}" 
set output '00_rr_compare_only_mirk.tex' 


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


set pointsize 0.7

set termoption dashed

# line styles


set style line 2  pt 2  ps 0.7 lt 1 lc rgb '#E41A1C' lw 0.5 # red
set style line 3  pt 3  ps 0.7 lt 1 lc rgb '#377EB8' lw 0.5 # blue
set style line 4  pt 4  ps 0.7 lt 1 lc rgb '#005A32' lw 0.5 # dark yellow-green  # '#4DAF4A' lw 0.5 # green
set style line 5  pt 5  ps 0.7 lt 1 lc rgb '#984EA3' lw 0.5 # purple
set style line 6  pt 6  ps 0.7 lt 1 lc rgb '#FF7F00' lw 0.5 # orange
set style line 7  pt 7  ps 0.7 lt 1 lc rgb '#66A61E' lw 0.5 # dark lime green    #'#FFFF33' lw 0.5 # yellow
set style line 8  pt 8  ps 0.7 lt 1 lc rgb '#A65628' lw 0.5 # brown
set style line 9  pt 9  ps 0.7 lt 1 lc rgb '#F781BF' lw 0.5 # pink
set style line 10 pt 10 ps 0.7 lt 1 lc rgb '#666666' lw 0.5 # dark gray
set style line 11 pt 13 ps 0.7 lt 1 lc rgb '#BDB76B' lw 0.5 # gold


# set style line 2  pt 2  ps 0.7 lt 1 lc rgb '#1B9E77' lw 0.5 # dark teal
# set style line 3  pt 3  ps 0.7 lt 1 lc rgb '#D95F02' lw 0.5 # dark orange
# set style line 4  pt 4  ps 0.7 lt 1 lc rgb '#7570B3' lw 0.5 # dark lilac
# set style line 5  pt 5  ps 0.7 lt 1 lc rgb '#E7298A' lw 0.5 # dark magenta
# set style line 6  pt 6  ps 0.7 lt 1 lc rgb '#66A61E' lw 0.5 # dark lime green
# set style line 7  pt 7  ps 0.7 lt 1 lc rgb '#E6AB02' lw 0.5 # dark banana
# set style line 8  pt 8  ps 0.7 lt 1 lc rgb '#A6761D' lw 0.5 # dark tan
# set style line 9  pt 9  ps 0.7 lt 1 lc rgb '#666666' lw 0.5 # dark gray
# set style line 10 pt 10 ps 0.7 lt 1 lc rgb '#D73027' lw 0.5 # red
# set style line 11 pt 11 ps 0.7 lt 1 lc rgb '#542788' lw 0.5 # dark purple


# set style line 2  pt 2  ps 0.7 lt 1 lc rgb '#1B9E77' lw 0.5 # dark teal
# set style line 3  pt 3  ps 0.7 lt 1 lc rgb '#D95F02' lw 0.5 # dark orange
# set style line 4  pt 4  ps 0.7 lt 1 lc rgb '#7570B3' lw 0.5 # dark lilac
# set style line 5  pt 5  ps 0.7 lt 1 lc rgb '#E7298A' lw 0.5 # dark magenta
# set style line 6  pt 6  ps 0.7 lt 1 lc rgb '#66A61E' lw 0.5 # dark lime green
# set style line 7  pt 7  ps 0.7 lt 1 lc rgb '#E6AB02' lw 0.5 # dark banana
# set style line 8  pt 8  ps 0.7 lt 1 lc rgb '#A6761D' lw 0.5 # dark tan
# set style line 9  pt 9  ps 0.7 lt 1 lc rgb '#666666' lw 0.5 # dark gray
# set style line 10 pt 10 ps 0.7 lt 1 lc rgb '#D73027' lw 0.5 # red
# set style line 11 pt 11 ps 0.7 lt 1 lc rgb '#542788' lw 0.5 # dark purple




# set linestyle 1  lt 22 pt 11  ps 0.7 lc rgb "#9e0142"   lw 0.5 
# set linestyle 2  lt 22 pt 22  ps 0.7 lc rgb "#d53e4f"   lw 0.5 
# set linestyle 3  lt 33 pt 33  ps 0.7 lc rgb "#f46d43"   lw 0.5
# set linestyle 4  lt 44 pt 44  ps 0.7 lc rgb "#fdae61"   lw 0.5
# set linestyle 5  lt 55 pt 55  ps 0.7 lc rgb "#fee08b"   lw 0.5
# set linestyle 6  lt 66 pt 66  ps 0.7 lc rgb "#e6f598"   lw 0.5
# set linestyle 7  lt 77 pt 77  ps 0.7 lc rgb "#abdda4"   lw 0.5
# set linestyle 8  lt 88 pt 88  ps 0.7 lc rgb "#66c2a5"   lw 0.5
# set linestyle 9  lt 99 pt 133  ps 0.7 lc rgb "#006837"   lw 0.5
# set linestyle 10  lt 99 pt 100  ps 0.7 lc rgb "cyan"     lw 0.5
# set linestyle 11  lt 99 pt 111  ps 0.7 lc rgb "#3288bd"  lw 0.5
# set linestyle 12  lt 99 pt 122  ps 0.7 lc rgb "#5e4fa2"  lw 0.5
# set linestyle 13  lt 99 pt 99   ps 0.7 lc rgb "#762a83"  lw 0.5

### Start multiplot (2x3 layout)
set multiplot layout 2,3 rowsfirst 

eval(init_margins(0.08, 0.94, 0.1, 0.95, 0.07, 0.025,2, 3))


set grid back 

set mxtics 5
set mytics 5

# --- GRAPH (1,1)

eval(set_margins(1,1))
set xrange[-0.5:0.5]
set xtics -0.6,0.2,0.6
set yrange[0:1.6]
set ytics 0.0,0.4,1.60

set format x ""
set format y '\lr %g'

unset key

set label 1 '\lr $\sigma = 10^6$' at graph 0.05,0.1 font 'Helvetica,18' front
set label 2 '\lr $p_g$'           at graph 0.05,0.9 font 'Helvetica,18' front

plot   'mirk/1st/sig_1e6/1d_data_cut_x_p.dat'         title 'MIRK1-OLD'      w p  ls 2 ,\
       'mirk/2nd/sig_1e6/1d_data_cut_x_p.dat'         title 'MIRK2-OLD'      w p  ls 3 ,\
       'imex/SSP2_LUM/sig_1e6/1d_data_cut_x_p.dat'    title 'SSP2(332)-LUM'  w p  ls 8 ,\
       'mirk/1st_new/sig_1e6/1d_data_cut_x_p.dat'     title 'MIRK1-NEW'      w p  ls 4 ,\
       'mirk/2nd_new/sig_1e6/1d_data_cut_x_p.dat'     title 'MIRK2-NEW'      w p  ls 6 ,\




# --- GRAPH (1,2)

set style data linespoints

eval(set_margins(2,1))

set yrange[0:1.6]
set ytics 0.0,0.4,1.60

set format x ""
set format y '\lr %g'

unset key

set label 1 '\lr $\sigma = 10^3$' at graph 0.05,0.1 font 'Helvetica,18' front


plot   'mirk/1st/sig_1e3/1d_data_cut_x_p.dat'         title 'MIRK1-OLD'      w p  ls 2 ,\
       'mirk/2nd/sig_1e3/1d_data_cut_x_p.dat'         title 'MIRK2-OLD'      w p  ls 3 ,\
       'imex/SSP2_LUM/sig_1e3/1d_data_cut_x_p.dat'    title 'SSP2(332)-LUM'  w p  ls 8 ,\
       'mirk/1st_new/sig_1e3/1d_data_cut_x_p.dat'     title 'MIRK1-NEW'      w p  ls 4 ,\
       'mirk/2nd_new/sig_1e3/1d_data_cut_x_p.dat'     title 'MIRK2-NEW'      w p  ls 6 ,\



# --- GRAPH (1,3)

set style data linespoints

eval(set_margins(3,1))


set yrange[0:2]
set ytics 0.0,0.4,2

set format x ""
#set format y ""

unset key

set label 1 '\lr $\sigma = 10^1$' at graph 0.05,0.1 font 'Helvetica,18' front


plot   'mirk/1st/sig_1e1/1d_data_cut_x_p.dat'         title 'MIRK1-OLD'      w p  ls 2 ,\
       'mirk/2nd/sig_1e1/1d_data_cut_x_p.dat'         title 'MIRK2-OLD'      w p  ls 3 ,\
       'imex/SSP2_LUM/sig_1e1/1d_data_cut_x_p.dat'    title 'SSP2(332)-LUM'  w p  ls 8 ,\
       'mirk/1st_new/sig_1e1/1d_data_cut_x_p.dat'     title 'MIRK1-NEW'      w p  ls 4 ,\
       'mirk/2nd_new/sig_1e1/1d_data_cut_x_p.dat'     title 'MIRK2-NEW'      w p  ls 6 ,\



# --- GRAPH (2,1)


set style data linespoints

eval(set_margins(1,2))

set xtics -0.6,0.2,0.6
set xrange[-0.5:0.5]
set yrange[-0.31:0.31]
set ytics  -0.3,0.15,0.3

set format x '\lr %g'
set format y '\lr %g'

unset key

set xlabel '\hg $y$'

set label 1 '' at graph 0.1,0.85 font 'Helvetica,18' front
set label 2 '\lr $E_z$'           at graph 0.05,0.9 font 'Helvetica,18' front

plot   1/0 w p ls 2 ps 1.5 title 'MIRK1-OLD',  'mirk/1st/sig_1e6/1d_data_cut_x_ez.dat'      w p  ls 2 notitle,\
       1/0 w p ls 3 ps 1.5 title 'MIRK2-OLD',  'mirk/2nd/sig_1e6/1d_data_cut_x_ez.dat'      w p  ls 3 notitle,\
       1/0 w p ls 8 ps 1.5 title 'SSP2(332)-LUM',    'imex/SSP2_LUM/sig_1e6/1d_data_cut_x_ez.dat' w p  ls 8 notitle,\
       1/0 w p ls 4 ps 1.5 title 'MIRK1-NEW',  'mirk/1st_new/sig_1e6/1d_data_cut_x_ez.dat'  w p  ls 4 notitle,\
       1/0 w p ls 6 ps 1.5 title 'MIRK1-NEW',  'mirk/2nd_new/sig_1e6/1d_data_cut_x_ez.dat'  w p  ls 6 notitle


# --- GRAPH (2,2)

set style data linespoints

eval(set_margins(2,2))

set yrange[-0.31:0.31]
set ytics  -0.3,0.15,0.3

set key bottom right font "arial,8"

set xlabel '\hg $y$'

set format x '\lr %g'
set format y '\lr %g'

unset key

set label 1 '' at graph 0.1,0.85 font 'Helvetica,18' front

plot   1/0 w p ls 2 ps 1.5 title 'MIRK1-OLD',  'mirk/1st/sig_1e3/1d_data_cut_x_ez.dat'      w p  ls 2 notitle,\
       1/0 w p ls 3 ps 1.5 title 'MIRK2-OLD',  'mirk/2nd/sig_1e3/1d_data_cut_x_ez.dat'      w p  ls 3 notitle,\
       1/0 w p ls 8 ps 1.5 title 'SSP2(332)-LUM',    'imex/SSP2_LUM/sig_1e3/1d_data_cut_x_ez.dat' w p  ls 8 notitle,\
       1/0 w p ls 4 ps 1.5 title 'MIRK1-NEW',  'mirk/1st_new/sig_1e3/1d_data_cut_x_ez.dat'  w p  ls 4 notitle,\
       1/0 w p ls 6 ps 1.5 title 'MIRK1-NEW',  'mirk/2nd_new/sig_1e3/1d_data_cut_x_ez.dat'  w p  ls 6 notitle


# --- GRAPH (2,3)

set style data linespoints

eval(set_margins(3,2))

set yrange[-0.15:0.15]
set ytics  -0.15,0.075,0.15

set key bottom right
set key spacing 1.3


set xlabel '\hg $y$'

set format x '\lr %g'
set format y '\lr %g'

set label 1 '' at graph 0.1,0.85 font 'Helvetica,18' front


plot   1/0 w p ls 2 ps 1.5 title 'MIRK1-OLD',  'mirk/1st/sig_1e1/1d_data_cut_x_ez.dat'      w p  ls 2 notitle,\
       1/0 w p ls 3 ps 1.5 title 'MIRK2-OLD',  'mirk/2nd/sig_1e1/1d_data_cut_x_ez.dat'      w p  ls 3 notitle,\
       1/0 w p ls 8 ps 1.5 title 'SSP2(332)-LUM',    'imex/SSP2_LUM/sig_1e1/1d_data_cut_x_ez.dat' w p  ls 8 notitle,\
       1/0 w p ls 4 ps 1.5 title 'MIRK1-NEW',  'mirk/1st_new/sig_1e1/1d_data_cut_x_ez.dat'  w p  ls 4 notitle,\
       1/0 w p ls 6 ps 1.5 title 'MIRK2-NEW',  'mirk/2nd_new/sig_1e1/1d_data_cut_x_ez.dat'  w p  ls 6 notitle


unset multiplot
### End multiplot

