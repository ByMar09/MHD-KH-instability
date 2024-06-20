==============================================
CUEVA CODE
==============================================

Cueva, from the acronyms, Computer Unit for EnVieroments in Astrophysics, solve numerically the relativistic non-ideal magnetohydrodynamic (RRMHD) system of equations.

Time integration is carried out with the Implicit-Explicit Runge-Kutta methods (IMEX-RK) as recommended by Palenzuela et al. (2009) and the novel Minimally Implicit Runge-Kutta methods (MIRK) created by Aloy and Cordero-Carrión (2016). The formulation of the RRMHD equations is based on finite-volume methods and numerical fluxes are calculated according to the prescription of the Local Lax-Friedrichs (LLF) approximate Riemann solver introduced by Alic et al. (2007) or Harten-Lax-Van Leer (HLL) Harten et al. (1983) approximate Riemann solver, as described in LeVeque (2002) or HLLC approximate Riemann solver for RRMHD equations presented in Miranda-Aranguren et al. (2018). The code employs different slope limiters for reconstruction techniques, imposing a total-variation-diminishing (TVD) property on the methods. The slopes of picewise-linear reconstruction type used are, the minmod slope limiter (MinMod) van Leer (1979), the monotonised central-difference limiter (MCL) van Leer (1977) and the superbee limiter (SBL) Roe (1985), all of second space acurracy order. Also we implement the fifth (MP5), seventh (MP7) and nineth (MP9) order, monotonicity-preserving schemes Suresh and Huynh (1997).


Cueva is free for everyone to use, with no restrictions. If this code is useful for your research, we invite you to cite the relevant paper (Miranda-Aranguren, S, et. al. “An HLLC Riemann Solver for Resistive Relativistic Magnetohydrodynamics.” Monthly Notices of the Royal Astronomical Society 476.3 (2018): 3837–3860) and hope that you can also contribute to the project.


==============================================
HOW TO COMPILE AND RUN
==============================================

Prerequisities: you need to have installed the programs *fortran and *OMP

-------------------
For UNIX machines
-------------------

Compile: you must to use the "make" command in your shell prompt as,
```
make
```

Run: to execute cueva code

1) first we need to establish the number of threads as
```
export OMP_NUM_THREADS=2
```
2) then to run in background we use the "nohup" command as
```
nohup time ./cueva.out  >& cueva.log &
```
3) to see the print exit in log file, we use
```
tail -f cueva.log
```

Note: for clean job directory

1) to clean $OBJETS files, use
```
make clean
```
2) to clean $OBJETS plus data files,  use
```
make veryclean
```

==============================================
AUTHOR
==============================================

        Sergio Miranda-Aranguren 
        Burjassot-Valencia-España (2018)

Under License GPLv3

