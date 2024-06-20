reset
set fit errorvariables



rho_up    =  1.
B0        =  1.
P0        =  B0**2/2
sigma_m   =  B0**2/rho_up
beta_m    =  2. * P0 / B0**2
VA_rec_tearing = 0.5  

#1. / sqrt(1./sigma_m + 2. * beta_m + 1.)

S_l   = 1e6
a_tm  = 0.05
Lr    = S_l**(1./3.) * a_tm

kappa = 1./a_tm
Lx    = 80. * a_tm
fac   = 0.1 * Lx  / (512.)
eta_p = 1e-9 #(Lr * VA_rec_tearing)/S_l

t1 = 20.0
t2 = 22.0

# define the functions depending on the current number
fstr(N) = sprintf("f%d(x) = a%d*x + b%d", N, N, N)

grow_bx_exp(N) = sprintf("0.5 * a%d", N)
grow_bx_err(N) = sprintf("0.5 * a%d_err", N)


# The fitting string for a specific file and the related function
fitstr(N) = sprintf("fit [2.0 * N : 2.0 * N + 2.0] f%d(x) 'TM_2D_bx_integral_vol.dat' u ($2):(log($3)) via a%d,b%d", N, N, N, N)

N = 10

# Do all the fits
do for [i=1:N-2] {

    eval(fstr(i))
    eval(fitstr(i))
    eval(grow_bx_exp(i))
    eval(grow_bx_err(i))	


}


#f(x) = a1 *x + b1
#fit [t1:t2] f(x) "delta_rec_t_2D.dat" u ($2):(log($3)) via a1, b1

#g(x) = a2 *x + b2
#fit [t1:t2] g(x) "TM_2D_bx_integral_vol.dat" u ($2):(log($3)) via a2, b2

#h(x) = a2a *x + b2a
#fit [t1:t2] h(x) "TM_2D_bx_max_val.dat" u ($2):(log($3)) via a2a, b2a

#i1(x) = a3 * x + b3
#fit [t1:t2] i1(x) "TM_2D_ex_integral_vol.dat" u ($2):(log($3)) via a3, b3

#i(x) = a4 *x + b4
#fit [t1:t2] i(x) "TM_2D_ey_integral_vol.dat" u ($2):(log($3)) via a4, b4

#j(x) = a5 *x + b5
#fit [t1:t2] j(x) "TM_2D_ez_integral_vol.dat" u ($2):(log($3)) via a5, b5



#grow_bx_exp       =  0.5 * a2
#grow_bx_err       =  0.5 * a2_err 

#grow_max_bx_exp   =  a2a
#grow_max_bx_err   =  a2a_err 

#grow_ex_exp       =  a3
#grow_ex_err       =  a3_err 

#grow_ey_exp       =  a4
#grow_ey_err       =  a4_err 

#grow_ez_exp       =  a5
#grow_ez_err       =  a5_err 


#eta_t   = a1/kappa**2
#eta_err = a1_err/kappa**2
#eta_n   = eta_t - eta_p



set grid back

set mxtics 5
set mytics 5

set title  ' fit'

set xlabel 'time'  font "Times-Roman,16"

set key bottom left spacing 1.0

#set yrange [-16.0:0.0]


# construct the complete plotting string
plotstr = "plot "
do for [i=1:n] {
    plotstr = plotstr . sprintf("f%d(x), 'file%d.dat'%s ", i, i, (i == n) ? "" : ", ")
}

eval(plotstr)

#plot "TM_2D_bx_integral_vol.dat" u ($2):(log($3))  title '' w l lt 163  lc rgb "green" lw 1,\
#     g(x) title 'fit' lt 164 linecolor rgb "dark-green" lw 1




print "============================================="
print "============================================="
print " "
do for [i=1:N] {
print  "grow_bx_exp(i)", grow_bx_exp(i)
print " "
}
print "============================================="
print "============================================="
print "========= Have a good day! =================="
print "============================================="
