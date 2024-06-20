reset
set fit errorvariables



rho_up    =  5e-1
B0        =  sqrt(5.0)
P0        =  B0**2/10.0
sigma_m   =  B0**2/rho_up
beta_m    =  2. * P0 / B0**2
VA_rec_tearing = 1. / sqrt(1./sigma_m + 2. * beta_m + 1.)

S_a   = 5e3
a_tm  = 0.05
La    = S_a**(1./3.) * a_tm

kappa = 1./a_tm
Lx    = 80. * a_tm
fac   = 0.1 * Lx  / (512.)
eta_p = (a_tm * VA_rec_tearing)/S_a

# Interval time to made the fit adjust

t1 = 20.0
t2 = 25.0


f(x) = a1 *x + b1
fit [t1:t2] f(x) "delta_rec_t_2D.dat" u ($2):(log($3)) via a1, b1

g(x) = a2 *x + b2
fit [t1:t2] g(x) "TM_2D_bx_integral_vol.dat" u ($2):(log($3)) via a2, b2

h(x) = a2a *x + b2a
fit [t1:t2] h(x) "TM_2D_bx_max_val.dat" u ($2):(log($3)) via a2a, b2a

i1(x) = a3 * x + b3
fit [t1:t2] i1(x) "TM_2D_ex_integral_vol.dat" u ($2):(log($3)) via a3, b3

i(x) = a4 *x + b4
fit [t1:t2] i(x) "TM_2D_ey_integral_vol.dat" u ($2):(log($3)) via a4, b4

j(x) = a5 *x + b5
fit [t1:t2] j(x) "TM_2D_ez_integral_vol.dat" u ($2):(log($3)) via a5, b5


# Vertical lines, point times  where we made the adjust

set arrow from t1,graph(0,0) to t1,graph(1,1) nohead lw 1.5 lc rgb "orange"
set arrow from t2,graph(0,0) to t2,graph(1,1) nohead lw 1.5 lc rgb "orange"


grow_bx_exp       =  0.5 * a2
grow_bx_err       =  0.5 * a2_err 

grow_max_bx_exp   =  a2a
grow_max_bx_err   =  a2a_err 

grow_ex_exp       =  a3
grow_ex_err       =  a3_err 

grow_ey_exp       =  a4
grow_ey_err       =  a4_err 

grow_ez_exp       =  a5
grow_ez_err       =  a5_err 


eta_t   = a2/kappa**2
eta_err = a2_err/kappa**2
eta_n   = eta_t - eta_p



set grid back

set mxtics 5
set mytics 5

set title  ' fit'

set xlabel 'time'  font "Times-Roman,16"

set key bottom left spacing 1.0

# To open a terminal named 0

set term x11 0

plot "TM_2D_bx_integral_vol.dat" u ($2):(log($3))  title '' w l lt 163  lc rgb "green" lw 1.2,\
     g(x) title 'fit' lt 164 linecolor rgb "dark-green" lw 1.2


# Now we are going to calculate 0.5 * \partial_t (\ln(\abs(bx^2)))


r(x) = 0.5 * a2

x0=NaN

y0=NaN

set yrange [0:1]

set key at graph 0.25,0.05

# To open another terminal named 1

set term x11 1

plot "TM_2D_bx_integral_vol.dat" u (dx=$2-x0,x0=$2,$2-dx/2):(dy=(log(abs($3))/2.0)-y0,y0=(log(abs($3))/2.0),dy/dx) w l lw 1.2 notitle,\
     r(x) title "growth rate" w l lc rgb "green" lw 1.2


print "============================================="
print "============================================="
print " "
print " t1 =  ", t1
print " t2 =  ", t2
print " "
print " "
print " grow exponent \gamma_bx     = ", grow_bx_exp, "  +/-  ", grow_bx_err
print " "
print " "
print " grow exponent max(val) \gamma_bx     = ", grow_max_bx_exp, "  +/-  ", grow_max_bx_err
print " "
print " "
print " grow exponent \gamma_ex     = ", grow_ex_exp, "  +/-  ", grow_ex_err
print " "
print " "
print " grow exponent \gamma_ey     = ", grow_ey_exp, "  +/-  ", grow_ey_err
print " "
print " "
print " grow exponent \gamma_ez     = ", grow_ez_exp, "  +/-  ", grow_ez_err
print " "
#print "============================================="
#print " physical resistivity    eta_p = ", eta_p
#print " total resistivity       eta_t = ", eta_t, "  +/-  ", eta_err
#print " numerical resistivity   eta*  = ", eta_n, "  +/-  ", eta_err
print "============================================="
print "============================================="
print "========= Have a good day! =================="
print "============================================="
print "============================================="
