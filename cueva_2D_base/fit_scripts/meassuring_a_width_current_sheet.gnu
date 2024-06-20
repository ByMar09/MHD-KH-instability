reset

B0 = 1.0


stats "1d_data_cut_y_by.dat" u ($1 < 0 ? $1 : 1/0):( $2 <= B0 * tanh(1.0) && $2 >= - B0 * tanh(1.0) ? $2 : 1/0)

x1 = STATS_min_x
x2 = STATS_max_x

print "============================================="
print "============================================="
print " "
print "current sheet width a:"
print  (x2-x1)/2.
print " "
print "============================================="
print "============================================="

set arrow from x1,graph(0,0) to x1,graph(1,1) nohead lw 1.5 lc rgb "orange"
set arrow from x2,graph(0,0) to x2,graph(1,1) nohead lw 1.5 lc rgb "orange"

set grid back

set mxtics 10
set mytics 10

set key bottom left

plot "1d_data_cut_y_by.dat" u ($1):($2) title "by" w lp lt 9  lc rgb "dark-green" lw 1.2
