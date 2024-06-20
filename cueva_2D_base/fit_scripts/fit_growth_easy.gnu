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

t1 = 10.0
t2 = 20.0

t3 = 20.0
t4 = 30.0

t5 = 20.0
t6 = 40.0


g(x) = a2 *x + b2
fit [t1:t2] g(x) "TM_2D_bx_integral_vol.dat" u ($2):(log($3)) via a2, b2

g1(x) = a21 *x + b21
fit [t3:t4] g1(x) "TM_2D_bx_integral_vol.dat" u ($2):(log($3)) via a21, b21

g2(x) = a22 *x + b22
fit [t5:t6] g2(x) "TM_2D_bx_integral_vol.dat" u ($2):(log($3)) via a22, b22





h(x) = a2a *x + b2a
fit [t1:t2] h(x) "TM_2D_bx_max_val.dat" u ($2):(log($3)) via a2a, b2a


# Vertical lines, point times  where we made the adjust

set arrow from t1,graph(0,0) to t1,graph(1,1) nohead lw 1.5 lc rgb "orange"
set arrow from t2,graph(0,0) to t2,graph(1,1) nohead lw 1.5 lc rgb "orange"


grow_bx_exp       =  0.5 * a2
grow_bx_err       =  0.5 * a2_err 

grow_bx1_exp       =  0.5 * a21
grow_bx1_err       =  0.5 * a21_err 

grow_bx2_exp       =  0.5 * a22
grow_bx2_err       =  0.5 * a22_err 

grow_max_bx_exp   =  a2a
grow_max_bx_err   =  a2a_err 


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

s(x) = 0.5 * a21

q(x) = 0.5 * a22

x0=NaN

y0=NaN

set yrange [0:1]

set key top right

# To open another terminal named 1

set term x11 1

plot "TM_2D_bx_integral_vol.dat" u (dx=$2-x0,x0=$2,$2-dx/2):(dy=(log(abs($3))/2.0)-y0,y0=(log(abs($3))/2.0),dy/dx) w l lw 1.2 notitle,\
     r(x) title "growth rate " w l lc rgb "green" lw 1.2,\
     s(x) title "growth rate 1" w l lc rgb "cyan" lw 1.2,\
     q(x) title "growth rate 2" w l lc rgb "blue" lw 1.2


print "============================================="
print "============================================="
print " "
print " grow rate gamma ", "            ", "        t1", "         ","  t2"
print   grow_bx_exp, "                  ",  t1, "           ",  t2
print   grow_bx1_exp, "                  ",  t3, "           ",  t4
print   grow_bx2_exp, "                  ",  t5, "           ",  t6
print " "
print " "
print " grow exponent \gamma_bx     = ", grow_bx_exp, "  +/-  ", grow_bx_err
print " "
print " "
print " grow exponent max(val) \gamma_bx     = ", grow_max_bx_exp, "  +/-  ", grow_max_bx_err
print " "
print " "
print " "
print "============================================="
print "========= Have a good day! =================="
print "============================================="
