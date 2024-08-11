
  !     *******************************************************************
  !     Subroutine that calculates the primitive final variables
  !     *******************************************************************

  subroutine varprimitive_nr

    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

       call boundary_conserved
       call boundary_electric
       call boundary_flux_conserved


!$OMP   PARALLEL PRIVATE(TOL,n_elec,CHA_ELEC_1,CHA_ELEC_2,CHA_ELEC_3, Vx_elec_0, Vy_elec_0, Vz_elec_0),                           &
!$OMP & PRIVATE(Ex_elec_0, Ey_elec_0, Ez_elec_0, Bx_elec, By_elec, Bz_elec, D_elec, tau_elec, enthpy_elec, W_elec_2),             &
!$OMP & PRIVATE(Sx_elec, Sy_elec, Sz_elec, rho_elec, p_elec, V2c, W_elec, E2B2, gamma_1, C_1, C_2, a_4, a_3, a_2, a_1, a_0),      &
!$OMP & PRIVATE(b_3, b_2, b_1, b_0, d_2, d_1, d_0, e, f, g, Q_1, R_1, T_1, z_1a, s_1a, t_1a, z_1b, s_1b, t_1b, z_1c, s_1c, t_1c), &
!$OMP & PRIVATE(det, s, t, z_1, xx1, AA_1, BB_1, Q_2, y_1, y_2, Wnro, nrw, CHA1, Wnr, pnro, nr, epsilonnr, Cs2, CHA, pnr)
    
!$OMP DO     

    do j=0+fix,jmax
       do i=0,imax

     TOL   = 1.d-10

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


    Vx_elec     = (Sx_elec  - Ey_elec_0*Bz_elec + Ez_elec_0*By_elec) /       &
                  (tau_elec - E2B2 + p_elec)
    Vy_elec     = (Sy_elec  - Ez_elec_0*Bx_elec + Ex_elec_0*Bz_elec) /       &
                  (tau_elec - E2B2 + p_elec)
    Vz_elec     = (Sz_elec  - Ex_elec_0*By_elec + Ey_elec_0*Bx_elec) /       &
                  (tau_elec - E2B2 + p_elec)

    V2c         = Vx_elec**2 + Vy_elec**2 + Vz_elec**2

    W_elec      = 1.d0/sqrt(1.d0-V2c)
    W_elec_2    = 1.d0/    (1.d0-V2c)


! Once the Lorentz factor is know, we can recover the primitive variables as:

     rho_elec    =  D_elec / W_elec

     if (rho_elec.lt.0.d0) then
        print*, "DENSIDAD negativa ¡¡¡ ", rho_elec, i
     end if

     !     Specific internal energy

     epsilonnr   = (tau_elec - E2B2 - D_elec * W_elec + p_elec * (1.d0 - W_elec_2)) / &
                   (D_elec   * W_elec)

     !     Local speed of the fluid

     Cs2         = (gamma * (gamma - 1.d0) * epsilonnr) /  &
                   (1.d0  +  gamma         * epsilonnr)

     !******************************************************
     !      Newton Raphson for Pressure
     !******************************************************
     
     CHA   = 1.d0
     nr    =  1

     pnro  =  p_elec   !seed value for pressure
     
       do  while (CHA .ge. TOL .and. nr .le. nrmax)

          pnr =  pnro  - ( (gamma-1.d0) * rho_elec * epsilonnr - pnro ) / ( V2c * Cs2 - 1.d0 )

        !     Convergence test  CHA= |p(n)-p(n+1)|/(0.5*(p(n)+p(n+1))) !Toro E.F (4.45)

          CHA = abs(pnr-pnro) / (0.5d0* (pnr+pnro))

          pnro =  pnr

!       write(*,*) "NR Pressure", pnr, nr, i

          nr = nr+1

       end do     ! for nr
!******************************************************

!!$    if( pnro .lt. 0.d0) then
!!$       print*, "presion negativa ¡¡¡", pnro,rho_elec,W_elec,i,j,l,l_swap
!!$    end if


       p_elec      =  pnro


!     enthpy_elec = (tau_elec  - E2B2 - gamma_1 * rho_elec) / &
!                   (W_elec_2 - gamma_1 )        

!    tau_elec    =  E2B2 + enthpy_elec * W_elec**2 - p_elec

    Vx_elec     = (Sx_elec  - Ey_elec_0*Bz_elec + Ez_elec_0*By_elec) /       &
                  (tau_elec - E2B2 + p_elec)
    Vy_elec     = (Sy_elec  - Ez_elec_0*Bx_elec + Ex_elec_0*Bz_elec) /       &
                  (tau_elec - E2B2 + p_elec)
    Vz_elec     = (Sz_elec  - Ex_elec_0*By_elec + Ey_elec_0*Bx_elec) /       &
                  (tau_elec - E2B2 + p_elec)

    V2c         = Vx_elec**2 + Vy_elec**2 + Vz_elec**2

    W_elec      = 1.d0/sqrt(1.d0-V2c)
    W_elec_2    = 1.d0/    (1.d0-V2c)

    rho_elec    =  D_elec / W_elec

    if (V2c .ge. 1.d0 ) then
       write(*,*) "******************************"
       write(*,*) "Super luminal velocity, Newton-Raphson final step"
       write(*,*) " V^2 ---> ", h, i, V2c 
       write(*,*) "Super luminal velocity Newton-Raphson  final step"
       write(*,*) "******************************"
       stop
    end if

    CHA_ELEC_1 = 1.d0
    CHA_ELEC_2 = 1.d0
    CHA_ELEC_3 = 1.d0

    n_elec     = 1
    
    do  while (max(CHA_ELEC_1,CHA_ELEC_2,CHA_ELEC_3) > TOL .and. n_elec < nrmax)
 
       if (MIRK == 0 .and. REC_ELEC == 1 ) then

          if (l_swap == 2) l = lmax

    !     Auxiliar intermediate value of E (solution of the implicit part)

    a(l)    = Abt(l,l) * sigma * Delt

    matx(l) = W_elec_2 * a(l) + W_elec * a(l)**2      &
            + W_elec    + a(l)

    vrotB_x = (Bz_elec*Vy_elec-By_elec*Vz_elec)
    vrotB_y = (Bx_elec*Vz_elec-Bz_elec*Vx_elec)
    vrotB_z = (By_elec*Vx_elec-Bx_elec*Vy_elec)

    Ex_elec = (1./matx(l)) * (                                               &
         (a(l) + W_elec    + a(l) * W_elec_2 * Vx_elec**2 )* Ex_elec_0 +    &
          a(l) * W_elec_2  * Vx_elec         * Vy_elec     * Ey_elec_0 +    &
          a(l) * W_elec_2  * Vx_elec         * Vz_elec     * Ez_elec_0 -    &
          Abt(l,l) * sigma * W_elec          * Delt        * (              &
         (a(l) + W_elec    + a(l) * W_elec_2 * Vx_elec**2 )* vrotB_x   +    &
          a(l) * W_elec_2  * Vx_elec         * Vy_elec     * vrotB_y   +    &
          a(l) * W_elec_2  * Vx_elec         * Vz_elec     * vrotB_z    ))

    Ey_elec = (1./matx(l)) * (                                               &
          a(l) * W_elec_2  * Vx_elec         * Vy_elec     * Ex_elec_0 +    &
         (a(l) + W_elec    + a(l) * W_elec_2 * Vy_elec**2 )* Ey_elec_0 +    &
          a(l) * W_elec_2  * Vy_elec         * Vz_elec     * Ez_elec_0 -    &
          Abt(l,l) * sigma * W_elec          * Delt        * (              &
          a(l) * W_elec_2  * Vx_elec         * Vy_elec     * vrotB_x   +    &
         (a(l) + W_elec    + a(l) * W_elec_2 * Vy_elec**2 )* vrotB_y   +    &
          a(l) * W_elec_2  * Vy_elec         * Vz_elec     * vrotB_z    ))

    Ez_elec = (1./matx(l)) * (                                               &
          a(l) * W_elec_2  * Vx_elec         * Vz_elec     * Ex_elec_0 +    &
          a(l) * W_elec_2  * Vy_elec         * Vz_elec     * Ey_elec_0 +    &
         (a(l) + W_elec    + a(l) * W_elec_2 * Vz_elec**2 )* Ez_elec_0 -    &
          Abt(l,l) * sigma * W_elec          * Delt        * (              &
          a(l) * W_elec_2  * Vx_elec         * Vz_elec     * vrotB_x   +    &
          a(l) * W_elec_2  * Vy_elec         * Vz_elec     * vrotB_y   +    &
         (a(l) + W_elec    + a(l) * W_elec_2 * Vz_elec**2 )* vrotB_z    ))

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

      n_elec = n_elec + 1

    end do     ! end do while for n_elec


111 Vx_elec_0 = Vx_elec
    Vy_elec_0 = Vy_elec
    Vz_elec_0 = Vz_elec


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


         end do ! for i
     end do ! for j

!$OMP END DO
!$OMP END PARALLEL
     
   end subroutine varprimitive_nr
