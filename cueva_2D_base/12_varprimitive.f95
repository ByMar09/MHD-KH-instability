
  !     *******************************************************************
  !     Subroutine that calculates the primitive final variables
  !     *******************************************************************

  subroutine varprimitivecardano

    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

       call boundary_conserved
       call boundary_electric
       call boundary_flux_conserved


!$OMP   PARALLEL PRIVATE(TOL,n_elec,CHA_ELEC_1,CHA_ELEC_2,CHA_ELEC_3, Vx_elec_0, Vy_elec_0, Vz_elec_0,W_drift),                   &
!$OMP & PRIVATE(Ex_elec_0, Ey_elec_0, Ez_elec_0, Bx_elec, By_elec, Bz_elec, D_elec, tau_elec, enthpy_elec, B2_drift, v_drift),    &
!$OMP & PRIVATE(Sx_elec, Sy_elec, Sz_elec, rho_elec, p_elec, V2c, W_elec, E2B2, gamma_1, C_1, C_2, a_4, a_3, a_2, a_1, a_0),      &
!$OMP & PRIVATE(b_3, b_2, b_1, b_0, d_2, d_1, d_0, e, f, g, Q_1, R_1, T_1, z_1a, s_1a, t_1a, z_1b, s_1b, t_1b, z_1c, s_1c, t_1c), &
!$OMP & PRIVATE(det, s, t, z_1, xx1, AA_1, BB_1, Q_2, y_1, y_2, Wnro, nrw, CHA1, Wnr, pnro, nr, epsilonnr, Cs2, CHA, pnr)


!___________________________________________________________________________________________________________________________________
       !Headers for control files

!$OMP SINGLE

     if (h == 1) write(301,*) "iteration    Delt     contador   i    j     RK-step     x(i)     y(j)   rho*       rho^n"
     if (h == 1) write(301,*) "========================================================================================"

     if (h == 1) write(302,*) "iteration    Delt     contador   i    j     RK-step     x(i)     y(j)   p*          p^n "
     if (h == 1) write(302,*) "========================================================================================"

     if (h == 1) write(303,*) "iteration    Delt     contador   i    j     RK-step     x(i)     y(j)   v*^2  v_drift^2 "
     if (h == 1) write(303,*) "========================================================================================"

!$OMP END SINGLE NOWAIT
!___________________________________________________________________________________________________________________________________


  
!$OMP DO     

    do j=0+fix,jmax
       do i=0,imax

     TOL   = 1.d-10

     n_elec     = 1
     CHA_ELEC_1 = 1.d0
     CHA_ELEC_2 = 1.d0
     CHA_ELEC_3 = 1.d0


     if (l_swap == 2) then

        Vx_elec_0   = Vxint (i,j,1,lmax) 
        Vy_elec_0   = Vyint (i,j,1,lmax)
        Vz_elec_0   = Vzint (i,j,1,lmax)

        Ex_elec_0   = Ex    (i,j,1) 
        Ey_elec_0   = Ey    (i,j,1)
        Ez_elec_0   = Ez    (i,j,1)

        Bx_elec     = Bx    (i,j,1)
        By_elec     = By    (i,j,1)
        Bz_elec     = Bz    (i,j,1)

        D_elec      = D     (i,j,1)
        tau_elec    = tau   (i,j,1)
        enthpy_elec = enthpy(i,j,1)

        Sx_elec     = Sx    (i,j,1)
        Sy_elec     = Sy    (i,j,1)
        Sz_elec     = Sz    (i,j,1)

        rho_elec    = rhoint(i,j,1,lmax)
        p_elec      = pint  (i,j,1,lmax)


      else if (l_swap == 1) then

        Vx_elec_0   = Vx    (i,j,1) 
        Vy_elec_0   = Vy    (i,j,1)
        Vz_elec_0   = Vz    (i,j,1)

        Ex_elec_0   = Ex    (i,j,1) 
        Ey_elec_0   = Ey    (i,j,1)
        Ez_elec_0   = Ez    (i,j,1)

        Bx_elec     = Bx    (i,j,1)
        By_elec     = By    (i,j,1)
        Bz_elec     = Bz    (i,j,1)

        D_elec      = D     (i,j,1)
        tau_elec    = tau   (i,j,1)
        enthpy_elec = enthpy(i,j,1)

        Sx_elec     = Sx    (i,j,1)
        Sy_elec     = Sy    (i,j,1)
        Sz_elec     = Sz    (i,j,1)

        rho_elec    = rho   (i,j,1)
        p_elec      = p     (i,j,1)


     else if (l_swap == 0) then

        Vx_elec_0   = Vxint    (i,j,1,l) 
        Vy_elec_0   = Vyint    (i,j,1,l)
        Vz_elec_0   = Vzint    (i,j,1,l)

        Ex_elec_0   = Exint    (i,j,1,l) 
        Ey_elec_0   = Eyint    (i,j,1,l)
        Ez_elec_0   = Ezint    (i,j,1,l)

        Bx_elec     = Bxint    (i,j,1,l)
        By_elec     = Byint    (i,j,1,l)
        Bz_elec     = Bzint    (i,j,1,l)

        D_elec      = DDint    (i,j,1,l)
        tau_elec    = tauint   (i,j,1,l)
        enthpy_elec = enthpyint(i,j,1,l)

        Sx_elec     = Sxint    (i,j,1,l)
        Sy_elec     = Syint    (i,j,1,l)
        Sz_elec     = Szint    (i,j,1,l)

        rho_elec    = rhoint   (i,j,1,l)
        p_elec      = pint     (i,j,1,l)

     else

        write(*,*) "this parameter l_swap is not correctly"
        write(*,*) "l_swap=", l_swap
        stop

     end if



     if (SLC == 1) sigma = sigma_loc(i,j,1)


     E2B2= 0.5d0 * (Ex_elec_0**2 + Ey_elec_0**2 + Ez_elec_0**2 + &
                    Bx_elec**2   + By_elec**2   + Bz_elec**2   )
    
     gamma_1 = (gamma - 1.d0) / gamma

    !Drift velocity approximation
    
    if (VDRIFT == 1 .and. p_elec .le. (gamma - 1.d0) * E2B2) then !

        
       B2_drift    = Bx_elec**2 + By_elec**2 + Bz_elec**2


       Vx_elec     = (Ey_elec_0 * Bz_elec - Ez_elec_0 * By_elec) / E2B2
       Vy_elec     = (Ez_elec_0 * Bx_elec - Ex_elec_0 * Bz_elec) / E2B2
       Vz_elec     = (Ex_elec_0 * By_elec - Ey_elec_0 * Bx_elec) / E2B2

       
       v_drift     = sqrt( Vx_elec**2 + Vy_elec**2 + Vz_elec**2 )

       
          if (v_drift**2 .ge. 1.d0 ) then

!!$             write(*,*) "******************************"
!!$             write(*,*) "Super luminal velocity"
!!$             write(*,*) "Drift velocity approximation do not work INITIAL ¡¡¡"
!!$             write(*,*) " (h,i,j,v_drift**2) ---> ", h, i, j, v_drift**2 
!!$             write(*,*) " "
!!$             write(*,*) "******************************"
!!$             stop

             eps1 = 1.d-6
             
             print*, "Drift velocity approximation do not work in subroutine 12.varprimitive.f95"
             print*, "Before fixed, (v_drift, h, i ,j, l) --->", v_drift, h, i, j, l
             v_drift = max(-1.d0 + eps1 , min( 1.d0 - eps1, v_drift ))
             print*, "After fixed, (v_drift, h, i ,j, l) --->", v_drift, h, i, j, l
             
          end if

       V2c         = Vx_elec**2 + Vy_elec**2 + Vz_elec**2

       W_drift     = 1.d0/sqrt(1.d0 - v_drift**2)

       rho_elec    =  D_elec / W_drift

!!$       enthpy_elec = (tau_elec  - E2B2 - gamma_1 * rho_elec) / &
!!$                     (W_drift**2 - gamma_1 ) 
!!$
!!$       pnro        = gamma_1 * (enthpy_elec - rho_elec) !
!!$       
       p_elec      = 1.d-2 * E2B2 ! max(1.d-2 * pnro, p_elec) !

       Vx_elec_0   = Vx_elec
       Vy_elec_0   = Vy_elec
       Vz_elec_0   = Vz_elec

       W_elec = W_drift

       Wnro   = W_elec

       contador_v = contador_v + 1
         
       write(303,*) h, h * Delt, contador_v, i, j, l, 0.5d0 * Lx - i * Delx, j * Dely, V2c, v_drift**2

!!$       print*, "Passing through Drift Velocity approximation"



          
       else !else DRIFT velocity approximation (Normal precedure for less magnetized plasmas)


     do  while (max(CHA_ELEC_1,CHA_ELEC_2,CHA_ELEC_3) > TOL .and. n_elec < nrmax)

       CHA   = 1.d0

       !     Lorentz Factor

       V2c= Vx_elec_0**2 + Vy_elec_0**2 + Vz_elec_0**2

        W_elec= 1.d0/sqrt(1.d0 - V2c)    

    if (V2c > 1) then

       write(*,*) "Super luminal velocity, Cardano initial step"
       write(*,*) " (V^2, h , i, j, l) ---> ", V2c , h, i, j, l
       write(*,*) "FIXED IT"
       
       V2c  = max(-1.d0 + eps1 , min( 1.d0 - eps1, V2c ))
    
       W_elec= 1.d0/sqrt(1.d0 - V2c)    
       
    end if



!         quartic equation for Lorentz fator 
!  
!      a_4 W^4 + a_3 W^3 + 2_2 W^2 + a_1 W + a_0 = 0

    E2B2= 0.5d0 * (Ex_elec_0**2 + Ey_elec_0**2 + Ez_elec_0**2 + &
                   Bx_elec**2   + By_elec**2   + Bz_elec**2   )

    gamma_1 = (gamma - 1.d0) / gamma

    C_1 =  (Sx_elec-Ey_elec_0*Bz_elec+Ez_elec_0*By_elec)**2 &
         + (Sy_elec-Ez_elec_0*Bx_elec+Ex_elec_0*Bz_elec)**2 &
         + (Sz_elec-Ex_elec_0*By_elec+Ey_elec_0*Bx_elec)**2


    C_2 =   tau_elec - E2B2


    a_4 = C_1 - C_2**2
    a_3 = 2.d0 * gamma_1 * C_2 * D_elec
    a_2 = C_2**2 - 2.d0 * gamma_1 * C_1 - gamma_1**2 * D_elec**2
    a_1 = - 2.d0 * gamma_1 * C_2 * D_elec
    a_0 = gamma_1**2 * (C_1 + D_elec**2)



! Ferrari-Cardano derivation of the quartic formula
!
!  given the quartic equation
!
!   W^4 + b_3 W^3 + b_2 W^2 + b_1 W + b_0 = 0

    b_3 = a_3 / a_4
    b_2 = a_2 / a_4
    b_1 = a_1 / a_4
    b_0 = a_0 / a_4

 

! Apply the Tchirnhaus transformation W ---> y - b_3/4 to obtain
!
!   y^4 + d_2 y^2 + d_1 y + d_0 = 0


    d_2 = b_2 - (3.d0 * b_3**2)/8.d0
    d_1 = b_1 - ( b_3 * b_2)/2.d0 + (b_3**3) / 8.d0
    d_0 = b_0 - ( b_3 * b_1)/4.d0 + (b_3**2 * b_2)/16.d0 - (3.d0 * b_3**4)/256.d0



    ! conditional for low or high velocities

    if ( low_vel == 1 .and. V2c .le. 1.d-2 ) then 


       V2c  = -(a_4 + a_3 + a_2 + a_1 + a_0) / &
            (2.d0 * a_4 + 1.5d0 * a_3 + a_2 + 0.5d0 * a_1)

       if (abs(V2c) .le. 1.d-10 ) V2c = 0.d0


       if ((2.d0 * a_4 + 1.5d0 * a_3 + a_2 + 0.5d0 * a_1) .eq. 0.d0 )  then

          write(*,*) "Huston we have a problem, V^2 indetermined cte/0 in (12_varprimitive.f95)", V2c, i, j, l, h
          stop

       end if
       
!       write(*,*) "passing throught condicional low velocities (12_varprimitive.f95)", V2c, i, j, l, h

       if ( V2c .lt. 0.d0) then

          write(*,*) "Huston we have a problem, V^2 < 0 in condicional low velocities (12_varprimitive.f95)", V2c, i, j, l, h
          stop

       end if

       Wnro = 1.d0/sqrt(1.d0 - V2c)  
       

    else
    

! For unespecified value z
!
!  ( y^2 +d_2 + z)^2 = (d_2 + 2 z) y^2 - d_1 y + (z^2 + 2 d_2 z + d_2^2 - d_0)  Equation(1)
!
! The goal is now to choose a value for z which makes the right hand side a perfet square
! that is, the discriminant in the quadratic polynomial in y must be zero
!
!        -8 z^3 - 20 d_2 z^2 + (8 d_0 - 16 d_2^2) z + (d_1^2 + 4 d_2 d_0 - 4 d_2^3) = 0
!          
!           z^3 + e z^2 + f z + g = 0


    e = ( 20.d0 * d_2) / 8.d0
    f = - (d_0 - 2.d0 * d_2**2)
    g = - ( d_1**2 + 4.d0 * d_2 * d_0 - 4.d0 * d_2**3)/ 8.d0



! Choosing the real root for z

    Q_1 = (e**2 - 3.d0 * f) / 9.d0
    R_1 = (2.d0 * e**3 - 9.d0 * e * f + 27.d0 * g) / 54.d0

    if ( R_1**2 .lt. Q_1**3) then

       T_1 = acos( R_1 / abs(sqrt(Q_1**3)))

       z_1a = - 2.d0 * abs(sqrt(Q_1)) * cos((T_1 )/3.d0) - e/3.d0 
       s_1a =   d_2     + 2.d0 * z_1a  
       t_1a =   z_1a**2 + 2.d0 * d_2  * z_1a + d_2**2 - d_0    
  
       z_1b = - 2.d0 * abs(sqrt(Q_1)) * cos((T_1 + 2.d0 * pi )/3.d0) - e/3.d0
       s_1b =   d_2     + 2.d0 * z_1b
       t_1b =   z_1b**2 + 2.d0 * d_2  * z_1b + d_2**2 - d_0

       z_1c = - 2.d0 * abs(sqrt(Q_1)) * cos((T_1 - 2.d0 * pi )/3.d0) - e/3.d0
       s_1c =   d_2     + 2.d0 * z_1c
       t_1c =   z_1c**2 + 2.d0 * d_2  * z_1c + d_2**2 - d_0

       if      ( s_1a .ge. 0.d0 .and.  t_1a .ge. 0.d0) then
        
          s   = abs(sqrt( s_1a ))
          t   = abs(sqrt( t_1a ))
          det = s**2 - 4.d0 * (d_2 + z_1a - t)
          z_1 = z_1a

          if ( det  .lt. 0.d0) then

                write(*,*) " det < 0 para z_1a (h,i,j,l, det) --->", h,i,j,l,det
                print*, "D --->", D_elec
                print*, "tau --->", tau_elec
                print*, "Sx --->", Sx_elec
                print*, "Sy --->", Sy_elec
                print*, "Sz --->", Sz_elec
                stop
!                det = abs(det)

          end if

       else if ( s_1b .ge. 0.d0 .and.  t_1b .ge. 0.d0) then

          s   = abs(sqrt( s_1b ))
          t   = abs(sqrt( t_1b ))
          det = s**2 - 4.d0 * (d_2 + z_1b - t)
          z_1 = z_1b

          if ( det  .lt. 0.d0) then

                write(*,*) " det < 0 para z_1b  (h,i,j,l, det) --->", h,i,j,l,det
                print*, "D --->", D_elec
                print*, "tau --->", tau_elec
                print*, "Sx --->", Sx_elec
                print*, "Sy --->", Sy_elec
                print*, "Sz --->", Sz_elec
                stop
!                det = abs(det)

             end if
          
       else if ( s_1c .ge. 0.d0 .and.  t_1c .ge. 0.d0) then
        
          s   = abs(sqrt( s_1c ))
          t   = abs(sqrt( t_1c ))
          det = s**2 - 4.d0 * (d_2 + z_1c - t)
          z_1 = z_1c

             if ( det  .lt. 0.d0) then

                write(*,*) " det < 0 para z_1c (h,i,j,l, det) --->", h,i,j,l,det
                print*, "D --->", D_elec
                print*, "tau --->", tau_elec
                print*, "Sx --->", Sx_elec
                print*, "Sy --->", Sy_elec
                print*, "Sz --->", Sz_elec
                stop
!                det = abs(det)

             end if
       else
          write(*,*) "Que hacemos todo cardano IMAGINARIO FINAL", i,j,l,h
          stop
       end if

    else

       xx1  =    abs(R_1)      + abs (sqrt(R_1**2 - Q_1**3))
       AA_1 = - sign(1.d0,R_1) * sign(abs (xx1**THIRD), xx1)

       if ( AA_1 .ne. 0.d0) then

          BB_1 = Q_1 / AA_1

       else

          BB_1 = 0.d0

       end if

       z_1 = ( AA_1 + BB_1) - e/3.d0

       s   = abs(sqrt( d_2    + 2.d0 * z_1 ))
       t   = abs(sqrt( z_1**2 + 2.d0 * d_2 * z_1 + d_2**2 - d_0 ))
       det = s**2 - 4.d0 * (d_2 + z_1 - t)

    end if


! We may write equation (1) as
!
!     ( y^2 + d_2 + z_1)^2 = ( sy + t)^2  Equation (2)

! And taking the squared root of both side in Eq (2) and solving the resulting quadratic equation
!
! y^2 - s y + (d_2 + z_1 - t) = 0

     Q_2 = - 0.5d0 * ( - s + sign(1.d0,-s) * abs(sqrt( det )) )

! then the two roots are  y_1 = Q_2 and y_2 = (d_2 + z_1 - t) / Q_2  !Numerical Recipes 3er edition page 227

! Since we use de Tchirnhaus transformation W ---> y - b_3/4 we obtain

     y_1 = Q_2 - b_3/4.d0
     y_2 = (d_2 + z_1 - t) / Q_2 - b_3/4.d0


     if ( y_1 .ge. 1.d0 .and. y_2 .ge. 1.d0) then

!!$        write(*,*) "Huston we have a problem y_1 > 1 and y_2 > 1 FINAL", i,j,l,h
!!$        stop

        Wnro = min(y_1,y_2)

     end if

     if ( y_1 .le. 1.d0 .and. y_2 .le. 1.d0) then

        write(*,*) "Huston we have a problem y_1< 1 and y_2 < 1 (h,i,j,l,y1,y2) --->", h,i,j,l,y_1,y_2
        stop

     end if


     if ( y_1 .ge. 1.d0 .and. y_2 .le. 1.d0) then

        Wnro = y_1 

     end if

    if ( y_1 .le. 1.d0 .and. y_2 .ge. 1.d0) then

        Wnro  = y_2

     end if


   end if ! end conditional for low or high velocities
  

!******************************************************
     !      Newton Raphson for Lorentz factor
!******************************************************

     nrw           = 1
     CHA1          = 1.d0

     do  while (CHA1 > TOL .and. nrw < nrmax)

        Wnr  = Wnro - (Wnro**4        + b_3  * Wnro**3 + b_2 * Wnro**2 + b_1 * Wnro + b_0) / &
                      (4.d0 * Wnro**3 + 3.d0 * b_3 * Wnro**2 + 2.d0 * b_2 * Wnro    + b_1 )

        !     Convergence test  CHA= |p(n)-p(n+1)|/(0.5*(p(n)+p(n+1))) !Toro E.F (4.45)

        CHA1 = abs(Wnr-Wnro) / (0.5d0 * (Wnr+ Wnro))

        Wnro = Wnr

        nrw  = nrw + 1 

     end do     ! for nrw
!******************************************************

    W_elec = Wnro

! Once the Lorentz factor is know, we can recover the primitive variables as:

     rho_elec    =  D_elec / W_elec


     if (rho_elec.lt.0.d0) then

!        print*, "Huston we have negative density! (rho,D,i,j,l) ---> ", rho_elec,D_elec,i,j,l

        ! Truco chungo para evitar rho<0

        contador_r = contador_r + 1

       write(301,*) h, h * Delt, contador_r, i, j, l, 0.5d0 * Lx - i * Delx, j * Dely, rho_elec, rho(i,j,1)
        
       rho_elec    = max(rho(i,j,1) * 1.d-1, rho_elec)

     end if

     enthpy_elec = (tau_elec  - E2B2 - gamma_1 * rho_elec) / &
                   (W_elec**2 - gamma_1 ) 

     p_elec      = gamma_1 * (enthpy_elec - rho_elec)

!******************************************************
     !      Newton Raphson for Pressure
!******************************************************
     pnro        =  p_elec 
     nr          =  1

      !     Specific internal energy

     epsilonnr   = (tau_elec - E2B2 - D_elec * W_elec + pnro * (1.d0 - W_elec**2)) / &
                   (D_elec   * W_elec)

       !     Local speed of the fluid

     Cs2         = (gamma * (gamma - 1.d0) * epsilonnr) /  &
                   (1.d0  +  gamma         * epsilonnr)          
     
     CHA = 1.d0

       do  while (CHA .ge. TOL .and. nr .le. nrmax)

          pnr =  pnro  - ( (gamma-1.d0) * rho_elec * epsilonnr - pnro ) / ( V2c * Cs2 - 1.d0 )

        !     Convergence test  CHA= |p(n)-p(n+1)|/(0.5*(p(n)+p(n+1))) !Toro E.F (4.45)

          CHA = abs(pnr-pnro) / (0.5d0* (pnr+pnro))

          pnro =  pnr

!       write(*,*) "NR Pressure", pnr, nr, i

          nr = nr+1

       end do     ! for nr
!******************************************************



    if( pnro .lt. 0.d0) then

!       print*, "Huston we have negative pressure! (p,rho,W,i,j,l) --->", pnro,rho_elec,W_elec,i,j,l

       ! Truco chungo para evitar p < 0
       !ADVICE: This fix-up must be done over conserved variables and not primitives

       contador_p = contador_p + 1

       write(302,*) h, h * Delt, contador_p, i, j, l, 0.5d0 * Lx - i * Delx, j * Dely, pnro, p(i,j,1)
        
       pnro    = max(p(i,j,1) * 1.d-4, pnro)

!       print*, "valor de fix-up para presion asignado pnro =", pnro
       
    end if


    p_elec      =  pnro

 
    tau_elec    =  E2B2 + enthpy_elec * W_elec**2 - p_elec

    Vx_elec     = (Sx_elec  - Ey_elec_0*Bz_elec + Ez_elec_0*By_elec) /       &
                  (tau_elec - E2B2 + p_elec)
    Vy_elec     = (Sy_elec  - Ez_elec_0*Bx_elec + Ex_elec_0*Bz_elec) /       &
                  (tau_elec - E2B2 + p_elec)
    Vz_elec     = (Sz_elec  - Ex_elec_0*By_elec + Ey_elec_0*Bx_elec) /       &
                  (tau_elec - E2B2 + p_elec)

    V2c = Vx_elec**2 + Vy_elec**2 + Vz_elec**2

    W_elec= 1.d0/sqrt(1.d0-V2c)

          
    if (V2c .ge. 1.d0 ) then
       
       ErotB    = Ey_elec_0 * Bz_elec - Ez_elec_0 * By_elec + &
                  Ez_elec_0 * Bx_elec - Ex_elec_0 * Bz_elec + &
                  Ex_elec_0 * By_elec - Ey_elec_0 * Bx_elec

       B2_drift = Bx_elec**2 + By_elec**2 + Bz_elec**2

       if (B2_drift / rho_elec .ge. 10.d0) then

          v_drift = ErotB / B2_drift
          
          W_drift = 1.d0/sqrt(1.d0 - v_drift**2)

          rho_elec    =  D_elec / W_drift

          enthpy_elec = (tau_elec  - E2B2 - gamma_1 * rho_elec) / &
                        (W_drift**2 - gamma_1 ) 

          p_elec      = gamma_1 * (enthpy_elec - rho_elec)

          tau_elec    =  E2B2 + enthpy_elec * W_elec**2 - p_elec

          Vx_elec     = (Sx_elec  - Ey_elec_0*Bz_elec + Ez_elec_0*By_elec) /       &
               (tau_elec - E2B2 + p_elec)
          Vy_elec     = (Sy_elec  - Ez_elec_0*Bx_elec + Ex_elec_0*Bz_elec) /       &
               (tau_elec - E2B2 + p_elec)
          Vz_elec     = (Sz_elec  - Ex_elec_0*By_elec + Ey_elec_0*Bx_elec) /       &
               (tau_elec - E2B2 + p_elec)

          V2c    = Vx_elec**2 + Vy_elec**2 + Vz_elec**2

          W_elec = 1.d0/sqrt(1.d0-V2c)

          contador_v = contador_v + 1
         
          write(303,*) h, h * Delt, contador_v, i, j, l, 0.5d0 * Lx - i * Delx, j * Dely, V2c, v_drift**2

!          print*, "Passing through Drift Velocity approximation"


          if (V2c .ge. 1.d0 ) then

             write(*,*) "******************************"
             write(*,*) "Super luminal velocity"
             write(*,*) "Drift velocity approximation do not work ¡¡¡"
             write(*,*) " (h,i,v2) ---> ", h, i, V2c
             write(*,*) "The value of magnetization is:"
             write(*,*) " sigma_m = |B^2| / rho = ", B2_drift / rho_elec
             write(*,*) " "
             write(*,*) "******************************"
             stop

          end if
          
       else

          write(*,*) "******************************"
          write(*,*) "Super luminal velocity, Cardano final step"
          write(*,*) " (h,i,v2) ---> ", h, i, V2c 
          write(*,*) "It is not correct to do the Drift Velocity approximation"
          write(*,*) "The value of magnetization is to low"
          write(*,*) " sigma_m = |B^2| / rho = ", B2_drift / rho_elec
          write(*,*) "******************************"
          stop

       end if

       
    end if



 
    if (MIRK == 0 .and. REC_ELEC == 1 ) then 

       if (l_swap == 2) l = lmax


    !     Auxiliar intermediate value of E (solution of the implicit part)

    a(l)    = Abt(l,l) * sigma * Delt

    matx(l) = W_elec**2 * a(l) + W_elec * a(l)**2      &
            + W_elec    + a(l)

    vrotB_x = (Bz_elec*Vy_elec-By_elec*Vz_elec)
    vrotB_y = (Bx_elec*Vz_elec-Bz_elec*Vx_elec)
    vrotB_z = (By_elec*Vx_elec-Bx_elec*Vy_elec)

    Ex_elec = (1./matx(l)) * (                                               &
         (a(l) + W_elec    + a(l) * W_elec**2 * Vx_elec**2 )* Ex_elec_0 +    &
          a(l) * W_elec**2 * Vx_elec          * Vy_elec     * Ey_elec_0 +    &
          a(l) * W_elec**2 * Vx_elec          * Vz_elec     * Ez_elec_0 -    &
          Abt(l,l) * sigma * W_elec           * Delt        * (              &
         (a(l) + W_elec    + a(l) * W_elec**2 * Vx_elec**2 )* vrotB_x   +    &
          a(l) * W_elec**2 * Vx_elec          * Vy_elec     * vrotB_y   +    &
          a(l) * W_elec**2 * Vx_elec          * Vz_elec     * vrotB_z    ))

    Ey_elec = (1./matx(l)) * (                                               &
          a(l) * W_elec**2 * Vx_elec          * Vy_elec     * Ex_elec_0 +    &
         (a(l) + W_elec    + a(l) * W_elec**2 * Vy_elec**2 )* Ey_elec_0 +    &
          a(l) * W_elec**2 * Vy_elec          * Vz_elec     * Ez_elec_0 -    &
          Abt(l,l) * sigma * W_elec           * Delt        * (              &
          a(l) * W_elec**2 * Vx_elec          * Vy_elec     * vrotB_x   +    &
         (a(l) + W_elec    + a(l) * W_elec**2 * Vy_elec**2 )* vrotB_y   +    &
          a(l) * W_elec**2 * Vy_elec          * Vz_elec     * vrotB_z    ))

    Ez_elec = (1./matx(l)) * (                                               &
          a(l) * W_elec**2 * Vx_elec          * Vz_elec     * Ex_elec_0 +    &
          a(l) * W_elec**2 * Vy_elec          * Vz_elec     * Ey_elec_0 +    &
         (a(l) + W_elec    + a(l) * W_elec**2 * Vz_elec**2 )* Ez_elec_0 -    &
          Abt(l,l) * sigma * W_elec           * Delt        * (              &
          a(l) * W_elec**2 * Vx_elec          * Vz_elec     * vrotB_x   +    &
          a(l) * W_elec**2 * Vy_elec          * Vz_elec     * vrotB_y   +    &
         (a(l) + W_elec    + a(l) * W_elec**2 * Vz_elec**2 )* vrotB_z    ))

    CHA_ELEC_1 = abs(Vx_elec - Vx_elec_0) / (0.5d0* (Vx_elec + Vx_elec_0))
    CHA_ELEC_2 = abs(Vy_elec - Vy_elec_0) / (0.5d0* (Vy_elec + Vy_elec_0))
    CHA_ELEC_3 = abs(Vz_elec - Vz_elec_0) / (0.5d0* (Vz_elec + Vz_elec_0))


    Ex_elec_0 = Ex_elec
    Ey_elec_0 = Ey_elec
    Ez_elec_0 = Ez_elec

   else


     CHA_ELEC_1 = 1.d-11
     CHA_ELEC_2 = 1.d-11
     CHA_ELEC_3 = 1.d-11
  
   end if


111 Vx_elec_0 = Vx_elec
    Vy_elec_0 = Vy_elec
    Vz_elec_0 = Vz_elec

    n_elec = n_elec + 1

   end do     ! end do while for n_elec



  end if ! end conditional DRIFT velocity



   if (l_swap == 1 .or. l_swap == 2) then 

      Vx (i,j,1) = Vx_elec_0
      Vy (i,j,1) = Vy_elec_0 
      Vz (i,j,1) = Vz_elec_0 

      Ex (i,j,1) = Ex_elec_0 
      Ey (i,j,1) = Ey_elec_0
      Ez (i,j,1) = Ez_elec_0 

      rho(i,j,1) = rho_elec  
      p  (i,j,1) = p_elec

      W  (i,j,1) = W_elec

   else if (l_swap == 0) then

      Vxint (i,j,1,l)  = Vx_elec_0   
      Vyint (i,j,1,l)  = Vy_elec_0  
      Vzint (i,j,1,l)  = Vz_elec_0  

      Exint (i,j,1,l)  = Ex_elec_0   
      Eyint (i,j,1,l)  = Ey_elec_0  
      Ezint (i,j,1,l)  = Ez_elec_0  

      rhoint(i,j,1,l)  = rho_elec  
      pint  (i,j,1,l)  = p_elec    

      Wint  (i,j,1,l)  = W_elec 

     else

        write(*,*) "this l_swap parameter is not correctly"
        write(*,*) "l_swap=", l_swap
        stop

     end if

!          print*, i,l, Vx_elec_0, "after cardano"


         end do ! for i
     end do ! for j

!$OMP END DO
!$OMP END PARALLEL
     
  end subroutine varprimitivecardano
