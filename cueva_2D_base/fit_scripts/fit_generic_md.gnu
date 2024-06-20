reset
set fit errorvariables

Lx = 1.
kappa = 2. * pi / Lx
eta_p = 1.e-5

# Interval time to made the fit adjust

t1 = 20.0
t2 = 25.0


f(x) = a1 *x + b1
fit [t1:t2] f(x) "data_fit.dat" u ($2):(log($3)) via a1, b1


# Vertical lines, point times  where we made the adjust

set arrow from t1,graph(0,0) to t1,graph(1,1) nohead lw 1.5 lc rgb "orange"
set arrow from t2,graph(0,0) to t2,graph(1,1) nohead lw 1.5 lc rgb "orange"


eta_t   = - 0.5 * a1/kappa**2
eta_err = - 0.5 * a1_err/kappa**2
eta_n   = eta_t - eta_p


set grid back

set mxtics 5
set mytics 5

set title  ' fit'

set xlabel 'time'  font "Times-Roman,16"

#set xrange [0:200]

set key bottom left spacing 1.0

set term x11 0

plot    'data_fit.dat' u ($2):(log($3))  title '' w l lt 163  lc rgb "green" lw 1,\
        f(x) title 'fit' lt 164 linecolor rgb "dark-green" lw 1

# Now we are going to calculate 0.5 * \partial_t (\ln(\abs(bx^2)))


r(x) = 0.5 * a1

x0=NaN

y0=NaN

set yrange [0:0.1]

set key at graph 0.25,0.05

# To open another terminal named 1

set term x11 1

plot "data_fit.dat" u (dx=$2-x0,x0=$2,$2-dx/2):(dy=(log(abs($3))/2.0)-y0,y0=(log(abs($3))/2.0),dy/dx) w l lc rgb "red" lw 1.2 notitle,\
     r(x) title "growth rate" w l lc rgb "green" lw 1.2


print "============================================="
print "============================================="
print " "
print " t1 =  ", t1
print " t2 =  ", t2
print " "
print " "
print "Box length  Lx =", Lx
print "Wave number kx =", kappa
print " "
print " "
print " a1 = ", - 0.5 * a1/kappa**2, "  +/-  ", - 0.5 * a1_err/kappa**2
print " "
print " "
print " physical resistivity    eta_p = ", eta_p
print " total resistivity       eta_t = ", eta_t, "  +/-  ", eta_err
print " numerical resistivity   eta*  = ", eta_n, "  +/-  ", eta_err
print " "
print "============================================="
print "============================================="
print "========= Have a good day! =================="
print "============================================="
