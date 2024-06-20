reset
set fit errorvariables

#set terminal pngcairo  color enhanced font "Times New Roman, 16"  fontscale 1.0 size 1800, 1800
#set output '00_fit_generic.png' 

Lx = 1.0e0
kappa = 2.0 * pi / Lx
eta_p = 1.0e-9


t1 = 4.3
t2 = 17.0


f(x) = a1 *x + b1
fit [t1:t2] f(x) "data_fit.dat" u ($2):(log($3)) via a1, b1



set arrow from t1,graph(0,0) to t1,graph(1,1) nohead lw 1.5 lc rgb "orange"
set arrow from t2,graph(0,0) to t2,graph(1,1) nohead lw 1.5 lc rgb "orange"


set grid back 

set mxtics 5
set mytics 5

set title  ' '

set xlabel 'time'  font "Times-Roman,16"

set key bottom left spacing 1.0 


###################### growth exponent ##########################

growth_exp1     = - 0.5 * a1
growth_exp1_err =   0.5 * a1_err

eta_plus_vis_1_n       = 2. * growth_exp1    /kappa**2 - eta_p
eta_plus_vis_1_n_err   = 2. * growth_exp1_err/kappa**2 


print "============================================="
print "============================================="
print "========= Have a good day! =================="
print "============================================="
print " physical resistivity    eta_p = ", eta_p
print "============================================="
print " "
print " Alfvén damping rate D_{A}            =", growth_exp1,       "  +/-", growth_exp1_err
print " numerical resistivity plus viscosity =", eta_plus_vis_1_n,  "  +/-", eta_plus_vis_1_n_err
print " fit parameter b                      =", b1,                "  +/-", b1_err
print " fit parameter a                      =", a1,                "  +/-", a1_err
print " "
print "============================================="
print "============================================="
print "========= Have a day! =================="
print "============================================="


set term x11 0

plot	'data_fit.dat' u ($2):(log($3))  title '' w l lt 163  lc rgb "green" lw 1,\
	f(x) title 'fit' lt 164 linecolor rgb "dark-green" lw 1

# Now we are going to calculate 0.5 * \partial_t (\ln(\abs(bx^2)))


r(x) = 0.5 * a1

x0=NaN

y0=NaN

#set yrange [0:0.1]

set key at graph 0.25,0.05

# To open another terminal named 1

set term x11 1

plot "data_fit.dat" u (dx=$2-x0,x0=$2,$2-dx/2):(dy=(log(abs($3))/2.0)-y0,y0=(log(abs($3))/2.0),dy/dx) w l lc rgb "red" lw 1.2 notitle,\
     r(x) title "growth rate" w l lc rgb "green" lw 1.2