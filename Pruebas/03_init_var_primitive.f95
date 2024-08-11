
!     ********************
!     Inicializar variables primitivas
!     ********************


subroutine init_var_primitive

  use scalar
  use parameters
  use threevectors
  use funciones

  implicit none


  !     Let''s be paranoid: init every 3D variable to zero

  if (DIM == 1) then


     do i=-6,imax+6

        psi(i,1,1) = 0.d0
        phi(i,1,1) = 0.d0
        Bx(i,1,1)  = 0.d0
        By(i,1,1)  = 0.d0
        Bz(i,1,1)  = 0.d0
        Ex(i,1,1)  = 0.d0
        Ey(i,1,1)  = 0.d0
        Ez(i,1,1)  = 0.d0
        Vx(i,1,1)  = 0.d0
        Vy(i,1,1)  = 0.d0
        Vz(i,1,1)  = 0.d0
        q(i,1,1)   = 0.d0
        p(i,1,1)   = 0.d0
        rho(i,1,1) = 0.d0

     end do   !i


  else if (DIM ==2) then



     do i=-6,imax+6
        do j=-6,jmax+6

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  = 0.d0
           By(i,j,1)  = 0.d0
           Bz(i,j,1)  = 0.d0
           Ex(i,j,1)  = 0.d0
           Ey(i,j,1)  = 0.d0
           Ez(i,j,1)  = 0.d0
           Vx(i,j,1)  = 0.d0
           Vy(i,j,1)  = 0.d0
           Vz(i,j,1)  = 0.d0
           q(i,j,1)   = 0.d0
           p(i,j,1)   = 0.d0
           rho(i,j,1) = 0.d0

        end do  !j
     end do   !i

  else if (DIM ==3) then



     do i=-6,imax+6
        do j=-6,jmax+6
           do k=-6,kmax+6

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              Bx(i,j,1)  = 0.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Ex(i,j,1)  = 0.d0
              Ey(i,j,1)  = 0.d0
              Ez(i,j,1)  = 0.d0
              Vx(i,j,1)  = 0.d0
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0
              q(i,j,1)   = 0.d0
              p(i,j,1)   = 0.d0
              rho(i,j,1) = 0.d0

           end do !k
        end do  !j
     end do   !i


  end if


  if (DIM == 1 .and. TEST == 1) then

     !  SCS TEST

     do i=0,imax


        posx = - 1.5d0 + i * Delx
        x    =   0.5d0 * sqrt(sigma) * posx

        psi(i,1,1) = 0.d0
        phi(i,1,1) = 0.d0
        Bx(i,1,1)  = 0.d0
        By(i,1,1)  = errorf(x) !0.d0 !
        Bz(i,1,1)  = 0.d0        !errorf(x,i) !
        Ex(i,1,1)  = 0.d0
        Ey(i,1,1)  = 0.d0
        Ez(i,1,1)  = 0.d0
        Vx(i,1,1)  = 0.d0
        Vy(i,1,1)  = 0.d0
        Vz(i,1,1)  = 0.d0
        q(i,1,1)   = 0.d0
        p(i,1,1)   = 50.d0
        rho(i,1,1) = 1.d0

        !_________________________________________________________________________________________________________________________
        write (26,130)  -1.5d0 + i * Delx, By(i,1,1)
        !_________________________________________________________________________________________________________________________


     end do   !i


  else if (DIM == 1 .and. TEST == 2) then

     ! CPAW TEST

     do i=0,imax


        posx = - 0.5d0 * Lx + i * Delx

        psi(i,1,1) = 0.d0
        phi(i,1,1) = 0.d0
        q(i,1,1)   = 0.d0
        p(i,1,1)   = P_up
        rho(i,1,1) = rho_up

        epsi_A     = P_up    / ((gamma - 1.d0  ) * rho_up )
        enthalpy_A = rho_up *  ( 1.d0  + epsi_A) + P_up


        VA = sqrt( 2.d0 * B0**2 /      ( ( enthalpy_A + B0**2 * (1.d0 + etha_A**2)) * &
             ( 1.d0 + sqrt  (1.d0 -( 2.d0 * etha_A * B0**2 / (enthalpy_A + B0**2 * (1.d0 + etha_A**2)))**2))))
!!$
!!$            VA =9.81d-002
!!$           VA = 0.1d0

        Bx(i,1,1)  = B0
        By(i,1,1)  = etha_A * B0 * cos(k_x * posx)
        Bz(i,1,1)  = etha_A * B0 * sin(k_x * posx)
        Vx(i,1,1)  = 0.d0
        Vy(i,1,1)  = - etha_A * VA * cos(k_x * posx)
        Vz(i,1,1)  = - etha_A * VA * sin(k_x * posx)
        Ex(i,1,1)  =   0.d0           !-(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))   !
        Ey(i,1,1)  =   VA * etha_A * B0 * sin(k_x * posx) !-(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))   !
        Ez(i,1,1)  = - VA * etha_A * B0 * cos(k_x * posx) !-(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))   !


        vrotB_x = (Bz(i,1,1)*Vy(i,1,1)-By(i,1,1)*Vz(i,1,1))
        vrotB_y = (Bx(i,1,1)*Vz(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
        vrotB_z = (By(i,1,1)*Vx(i,1,1)-Bx(i,1,1)*Vy(i,1,1))

        !           print*, i, vrotB_x



        !_________________________________________________________________________________________________________________________

        write (20,130)  - 0.5d0 * Longx + i * Delx, By(i,1,1)
        write (21,130)  - 0.5d0 * Longx + i * Delx, Bz(i,1,1)
        write (22,130)  - 0.5d0 * Longx + i * Delx, Vy(i,1,1)
        write (23,130)  - 0.5d0 * Longx + i * Delx, Vz(i,1,1)
        !_________________________________________________________________________________________________________________________


        !       print*, 'ALFVEN VELOCITY', VA


     end do   !i


  else if (DIM == 1 .and. TEST == 3) then

     ! RST_Palenzuela


     do i=0,imax

!!$          ! RST1-B0 OK (ST1-B0 hllc paper)
!!$          
!!$           if (i < 0.5 * imax) then
!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 0.d0
!!$           By(i,1,1)  = 0.5d0
!!$           Bz(i,1,1)  = 0.d0
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0
!!$           Vx(i,1,1)  = 0.d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.d0
!!$           rho(i,1,1) = 1.d0
!!$        else
!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 0.d0
!!$           By(i,1,1)  =-0.5d0
!!$           Bz(i,1,1)  = 0.d0
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0
!!$           Vx(i,1,1)  = 0.d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 0.1d0
!!$           rho(i,1,1) = 0.125d0
!!$
!!$        end if


        ! RST1 OK (RST1 Anton et. al. 2010 and ST1 hllc paper) sigma=1d4

        if (i < 0.5 * imax) then

           psi(i,1,1) = 0.d0
           phi(i,1,1) = 0.d0
           Bx(i,1,1)  = 0.5d0
           By(i,1,1)  = 1.0d0
           Bz(i,1,1)  = 0.d0
           Ex(i,1,1)  = 0.d0
           Ey(i,1,1)  = 0.d0
           Ez(i,1,1)  = 0.d0
           Vx(i,1,1)  = 0.d0
           Vy(i,1,1)  = 0.d0
           Vz(i,1,1)  = 0.d0
           q(i,1,1)   = 0.d0
           p(i,1,1)   = 1.d0
           rho(i,1,1) = 1.d0

        else

           psi(i,1,1) = 0.d0
           phi(i,1,1) = 0.d0
           Bx(i,1,1)  = 0.5d0
           By(i,1,1)  =-1.0d0
           Bz(i,1,1)  = 0.d0
           Ex(i,1,1)  = 0.d0
           Ey(i,1,1)  = 0.d0
           Ez(i,1,1)  = 0.d0
           Vx(i,1,1)  = 0.d0
           Vy(i,1,1)  = 0.d0
           Vz(i,1,1)  = 0.d0
           q(i,1,1)   = 0.d0
           p(i,1,1)   = 0.1d0
           rho(i,1,1) = 0.125d0

        end if

!!$        ! RST2   (RST Bucciantini and Del Zanna 2013 and ST2 hllc paper)  ok sigma = 1d4 (hllc)
!!$        
!!$        if (i < 0.5 * imax) then
!!$           psi(i,1,1) =  0.00d0
!!$           phi(i,1,1) =  0.00d0
!!$           Bx(i,1,1)  =  2.00d0
!!$           By(i,1,1)  =  0.30d0
!!$           Bz(i,1,1)  =  0.30d0
!!$           Vx(i,1,1)  =  0.40d0
!!$           Vy(i,1,1)  =  0.30d0
!!$           Vz(i,1,1)  =  0.20d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   =  0.00d0
!!$           p(i,1,1)   =  0.95d0
!!$           rho(i,1,1) =  1.08d0
!!$        else
!!$           psi(i,1,1) =  0.00d0
!!$           phi(i,1,1) =  0.00d0
!!$           Bx(i,1,1)  =  2.00d0
!!$           By(i,1,1)  = -0.70d0
!!$           Bz(i,1,1)  =  0.50d0
!!$           Vx(i,1,1)  = -0.45d0
!!$           Vy(i,1,1)  = -0.20d0
!!$           Vz(i,1,1)  =  0.20d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   =  0.00d0
!!$           p(i,1,1)   =  1.00d0
!!$           rho(i,1,1) =  1.00d0
!!$
!!$        end if



!!$           ! RST3  OK  (RST3 Anton et. al. 2010 and ST3 hllc paper)
!!$
!!$          if (i < 0.5 * imax) then
!!$             
!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 10.d0
!!$           By(i,1,1)  = 7.d0
!!$           Bz(i,1,1)  = 7.d0
!!$           Vx(i,1,1)  = 0.999d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = 0.d0
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.d-1
!!$           rho(i,1,1) = 1.d0
!!$           
!!$        else
!!$           
!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 10.d0
!!$           By(i,1,1)  = -7.d0
!!$           Bz(i,1,1)  = -7.d0
!!$           Vx(i,1,1)  = -0.999d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = 0.d0
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.d-1
!!$           rho(i,1,1) = 1.d0
!!$
!!$        end if



!!$        ! RST4 OK (RST4 Anton et. al. 2010 and ST4 hllc paper) sigma=1d4
!!$          
!!$        if (i < 0.5 * imax) then
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 1.0d0
!!$           By(i,1,1)  = 6.0d0
!!$           Bz(i,1,1)  = 2.0d0
!!$           Vx(i,1,1)  = 0.0d0
!!$           Vy(i,1,1)  = 0.3d0
!!$           Vz(i,1,1)  = 0.4d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.0d0
!!$           p(i,1,1)   = 5.0d0
!!$           rho(i,1,1) = 1.0d0
!!$        else
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 1.0d0
!!$           By(i,1,1)  = 5.0d0
!!$           Bz(i,1,1)  = 2.0d0
!!$           Vx(i,1,1)  = 0.0d0
!!$           Vy(i,1,1)  = 0.0d0
!!$           Vz(i,1,1)  = 0.0d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 5.3d0
!!$           rho(i,1,1) = 0.9d0
!!$
!!$        end if


!!$          ! RST5   (RST2 Mignone & Bodo 2006 and ST5 hllc paper)  ok sigma = 1d4 (hllc)
!!$
!!$          if (i < 0.5 * imax) then
!!$             
!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 5.d0
!!$           By(i,1,1)  = 6.d0
!!$           Bz(i,1,1)  = 6.d0
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0
!!$           Vx(i,1,1)  = 0.d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 3.d1
!!$           rho(i,1,1) = 1.d0
!!$           
!!$        else
!!$           
!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 5.d0
!!$           By(i,1,1)  = 7.d-1
!!$           Bz(i,1,1)  = 7.d-1
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0
!!$           Vx(i,1,1)  = 0.d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.d0
!!$           rho(i,1,1) = 1.d0
!!$
!!$        end if


!!$        ! Contact Wave (CW Anton et. al. 2010 and CW1 hllc paper) ok sig=1d4
!!$          
!!$        if (i < 0.5 * imax) then
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 5.0d0
!!$           By(i,1,1)  = 1.0d0
!!$           Bz(i,1,1)  = 0.5d0
!!$           Vx(i,1,1)  = 0.0d0
!!$           Vy(i,1,1)  = 0.7d0
!!$           Vz(i,1,1)  = 0.2d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.0d0
!!$           p(i,1,1)   = 1.0d0
!!$           rho(i,1,1) = 1.0d1
!!$        else
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 5.0d0
!!$           By(i,1,1)  = 1.0d0
!!$           Bz(i,1,1)  = 0.5d0
!!$           Vx(i,1,1)  = 0.0d0
!!$           Vy(i,1,1)  = 0.7d0
!!$           Vz(i,1,1)  = 0.2d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.0d0
!!$           rho(i,1,1) = 1.0d0
!!$
!!$        end if


!!$          !Contact Wave 2 (CW V. Honkkila & P. Janhunen and CW2 hllc paper) ok sig=1d4
!!$          
!!$        if (i < 0.5 * imax) then
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 1.0d0
!!$           By(i,1,1)  = 1.0d0
!!$           Bz(i,1,1)  = 0.0d0
!!$           Vx(i,1,1)  = 0.2d0
!!$           Vy(i,1,1)  = 0.0d0
!!$           Vz(i,1,1)  = 0.0d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.0d0
!!$           p(i,1,1)   = 1.0d0
!!$           rho(i,1,1) = 1.0d0
!!$        else
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 1.0d0
!!$           By(i,1,1)  = 1.0d0
!!$           Bz(i,1,1)  = 0.0d0
!!$           Vx(i,1,1)  = 0.2d0
!!$           Vy(i,1,1)  = 0.0d0
!!$           Vz(i,1,1)  = 0.0d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.0d0
!!$           rho(i,1,1) = 0.125d0
!!$
!!$        end if



!!$          !Rotational Wave (RW Anton et. al. 2010 and RW hllc paper) ok ssp2-lum
!!$          
!!$        if (i < 0.5 * imax) then
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 2.4d0
!!$           By(i,1,1)  = 1.0d0
!!$           Bz(i,1,1)  =-1.6d0
!!$           Vx(i,1,1)  = 0.4d0
!!$           Vy(i,1,1)  =-0.3d0
!!$           Vz(i,1,1)  = 0.5d0
        !           Ex(i,1,1)  = 0.d0
        !           Ey(i,1,1)  = 0.d0
        !           Ez(i,1,1)  = 0.d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.0d0
!!$           p(i,1,1)   = 1.0d0
!!$           rho(i,1,1) = 1.0d0
!!$        else
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 2.4d0
!!$           By(i,1,1)  =-0.1d0
!!$           Bz(i,1,1)  =-2.178213d0
!!$           Vx(i,1,1)  = 0.377347d0
!!$           Vy(i,1,1)  =-0.482389d0
!!$           Vz(i,1,1)  = 0.424190d0
        !           Ex(i,1,1)  = 0.d0
        !           Ey(i,1,1)  = 0.d0
        !           Ez(i,1,1)  = 0.d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.0d0
!!$           rho(i,1,1) = 1.0d0
!!$
!!$        end if


!!$           ! RST7 NO   (RST3 Mignone & Bodo 2006) NO      
!!$
!!$          if (i < 0.5 * imax) then
!!$             
!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 10.d0
!!$           By(i,1,1)  = 7.d0
!!$           Bz(i,1,1)  = 7.d0
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0
!!$           Vx(i,1,1)  = 0.d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.d3
!!$           rho(i,1,1) = 1.d0
!!$           
!!$        else
!!$           
!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 10.d0
!!$           By(i,1,1)  = 7.d-1
!!$           Bz(i,1,1)  = 7.d-1
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0
!!$           Vx(i,1,1)  = 0.d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.d-1
!!$           rho(i,1,1) = 1.d0
!!$
!!$        end if


     end do   !i

  else if (DIM ==1 .and. TEST == 4) then

     !RST_BZ13 (RST2 Anton et. al. 2010)

     do i=0,imax


        if (i < 0.5 * imax) then
           psi(i,1,1) =  0.00d0
           phi(i,1,1) =  0.00d0
           Bx(i,1,1)  =  2.00d0
           By(i,1,1)  =  0.30d0
           Bz(i,1,1)  =  0.30d0
           Vx(i,1,1)  =  0.40d0
           Vy(i,1,1)  =  0.30d0
           Vz(i,1,1)  =  0.20d0
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0
           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
           q(i,1,1)   =  0.00d0
           p(i,1,1)   =  0.95d0
           rho(i,1,1) =  1.08d0
        else
           psi(i,1,1) =  0.00d0
           phi(i,1,1) =  0.00d0
           Bx(i,1,1)  =  2.00d0
           By(i,1,1)  = -0.70d0
           Bz(i,1,1)  =  0.50d0
           Vx(i,1,1)  = -0.45d0
           Vy(i,1,1)  = -0.20d0
           Vz(i,1,1)  =  0.20d0
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0
           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
           q(i,1,1)   =  0.00d0
           p(i,1,1)   =  1.00d0
           rho(i,1,1) =  1.00d0

        end if


!!$!          Komissarov collision test x € [-2,2] t=1.22
!!$
!!$        if (i < 0.5 * imax) then
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 10.0d0
!!$           By(i,1,1)  = 10.0d0
!!$           Bz(i,1,1)  = 0.0d0
!!$           Vx(i,1,1)  = 5.0d0/sqrt(26.d0)
!!$           Vy(i,1,1)  = 0.0d0
!!$           Vz(i,1,1)  = 0.0d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.0d0
!!$           p(i,1,1)   = 1.0d0
!!$           rho(i,1,1) = 1.0d0
!!$        else
!!$           psi(i,1,1) = 0.0d0
!!$           phi(i,1,1) = 0.0d0
!!$           Bx(i,1,1)  = 10.0d0
!!$           By(i,1,1)  =-10.0d0
!!$           Bz(i,1,1)  = 0.0d0
!!$           Vx(i,1,1)  =- 5.0d0/sqrt(26.d0)
!!$           Vy(i,1,1)  = 0.0d0
!!$           Vz(i,1,1)  = 0.0d0
!!$           Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$           Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$           Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$           q(i,1,1)   = 0.d0
!!$           p(i,1,1)   = 1.0d0
!!$           rho(i,1,1) = 1.0d0
!!$
!!$        end if

     end do   !i

  else if (DIM ==1 .and. TEST == 5) then


     !////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     !RST_AlMi B orthogonal to V

!!$           Ws     = 2.d1
!!$           Vs     = sqrt(1.d0 - 1.d0/Ws**2)
!!$           psi_a  = 0.d0
!!$           phi_a  = 0.d0
!!$
!!$           Vx_a   = 0.55d0
!!$           Vz_a   = 0.25d0
!!$           W_a    = sqrt(1.d0/(1.d0 - Vx_a**2 - Vz_a**2))
!!$
!!$           rho_a  = 1.0d0
!!$           p_a    = 2.5d0
!!$           h_a    = 1.0d0 + (gamma/(gamma -1.d0))*(p_a/rho_a)
!!$
!!$           By_a   =  0.1d0
!!$           Ex_a   =  0.1d0 !By_a * Vz_a !0.1d0 !
!!$           Ey_a   =  0.1d0
!!$           Ez_a   =  0.1d0 !-By_a * Vx_a !0.1d0 !
!!$
!!$           b_a    = By_a/W_a
!!$           hh_a   = h_a + (By_a**2 - (Ex_a**2 + Ey_a**2 + Ez_a**2))/rho_a
!!$           ph_a   = p_a + 0.5d0 * (By_a**2 - (Ex_a**2 + Ey_a**2 + Ez_a**2))
!!$           q_a    = 0.0d0
!!$
!!$           J_inv  = Ws * rho_a * W_a * (Vs - Vx_a) ! J_inv = Ws * D * (Vs-Vx)
!!$
!!$           Jx_a   = sigma * W_a * ( Ex_a - Vz_a*By_a - (Ex_a*Vx_a + Ez_a*Vz_a)* Vx_a) + q_a * Vx_a
!!$           Jy_a   = sigma * W_a *   Ey_a
!!$           Jz_a   = sigma * W_a * ( Ez_a + Vx_a*By_a - (Ex_a*Vx_a + Ez_a*Vz_a)* Vz_a) + q_a * Vz_a
!!$
!!$           p_b    =  4.5d0
!!$           psi_b  =  0.0d0
!!$           phi_b  =  0.0d0
!!$!___________________________ R-H condition with out source  __________________ ___________________________________
!!$
!!$           By_b   = By_a !0.7d0 !By_a 
!!$           Ez_b   = Ez_a !+ Vs* (By_a - By_b)
!!$           Ex_b   = Ex_a 
!!$           Ey_b   = Ey_a 
!!$           q_b    =  q_a
!!$
!!$
!!$!___________________________ Using Jeffrey definition conservation law (3.1.5) ___________________________________
!!$
!!$           By_b   = By_a + Jz_a * Delx / (Ws * (Vs**2-1.d0))
!!$           Ez_b   = Ez_a - (By_b - By_a) / Vs - Jx_a * Delx / (Ws * Vs)
!!$           Ex_b   = Ex_a - (psi_b - psi_a) - 0.5d0 * (q_b - q_a) * Delx / Ws & 
!!$                         - 0.5d0 * kappa * (psi_b - psi_a) * Delx/ (Ws * (Vs-1.d0) ) &
!!$                         + (q_a - kappa * psi_a) * Delx / (Ws * (Vs - 1.d0)) - Jx_a * Delx / (Ws * (Vs-1.d0) )
!!$           Ey_b   = Ey_a - Jy_a * Delx / (Ws * Vs)
!!$
!!$!__________________________ Using LeVeque source as delta function (17.59) ________________________________________
!!$           By_b   = By_a +  sigma * W_a   * (Ez_a + Vx_a*By_a - (Ex_a*Vx_a + Ez_a*Vz_a)* Vz_a)/(Ws*(Vs**2 - 1.d0)) &
!!$                         -  q_a   * Vz_a  / (Ws * (Vs**2 - 1.d0))                          
!!$           Ez_b   = Ez_a - (By_a - By_b) - sigma * W_a *( Ez_a + Vx_a*By_a-(Ex_a*Vx_a+Ez_a*Vz_a)*Vz_a)/(Ws*(Vs - 1.d0)) &
!!$                         + q_a * Vz_a / (Ws*(Vs - 1.d0))
!!$           Ex_b   = Ex_a + (q_a   * (1.d0 - Vx_a)- sigma * W_a *( Ex_a - Vz_a*By_a-(Ex_a*Vx_a+Ez_a*Vz_a)*Vx_a))/(Ws*(Vs - 1.d0))
!!$           Ey_b   = Ey_a - sigma  * W_a * Ey_a /(Ws*(Vs     - 1.d0))
!!$!___________________________________________________________________________________________________________________
!!$
!!$           TOL1   = 1d-10
!!$           ERR_BY = 1.d0
!!$
!!$     do  while ( ERR_BY > TOL1 )
!!$
!!$
!!$           ph_b   = p_b + 0.5d0 * (By_b**2 - (Ex_b**2 + Ey_b**2 + Ez_b**2))
!!$
!!$ !           Vx_b   = (h_a * W_a * Vx_a + Ws *((ph_b-ph_a) + (Ez_b**2-Ez_a**2))/J_inv  &
!!$ !                 + (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * (By_a*Ez_a - By_b*Ez_b))    &
!!$ !                 / (h_a * W_a + (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * (ph_b-ph_a) + Ws * (By_a*Ez_a - By_b * Ez_b) /J_inv )                 
!!$
!!$ !          Vz_b   = (h_a * W_a * Vz_a + (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * (By_a*Ex_a - By_b*Ex_b) &
!!$ !                 + Ws * ( Ex_a*Ez_a-Ex_b*Ez_b)/J_inv)                                               &
!!$ !                 / (h_a * W_a + (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * (ph_b-ph_a) + Ws * (By_a*Ez_a - By_b * Ez_b)/J_inv)              
!!$
!!$           Vx_b   = (hh_a * W_a * Vx_a + Ws *( (ph_b-ph_a) + (Ez_b**2-Ez_a**2) - By_a*Ez_a*Vx_a + By_b*Ez_b*Vs)/J_inv &
!!$                  - (By_a**2 * Vx_a + By_a * Ez_a)/(rho_a  * W_a) )                                                   &
!!$                  / (hh_a * W_a - (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * (ph_a-ph_b)                                    &
!!$                  +  Ws   * (By_a*Ez_a-By_b*Ez_b)/J_inv + Ws * (By_a**2 * Vx_a - By_b**2 * Vs)/J_inv )
!!$
!!$           Vz_b   = (hh_a * W_a * Vz_a + Ws*(By_a*Ex_a*Vx_a-By_b*Ex_b*Vs)/J_inv + Ws * (Ex_a*Ez_a-Ex_b*Ez_b)/J_inv    &
!!$                  + (By_a * Ex_a-By_a**2 * Vz_a)/(rho_a*W_a))                                                         &
!!$                  / (hh_a * W_a -(Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * (ph_a-ph_b)                                     &
!!$                  +  Ws   *(By_a*Ez_a-By_b*Ez_b)/J_inv + Ws * (By_a**2*Vx_a-By_b**2*Vs)/J_inv)
!!$
!!$           if (Vx_b .gt. 1.d0 .or. Vz_b .gt. 1.d0) then
!!$
!!$              write(*,*) "Vx_b and Vz_b = ", Vx_b, Vz_b
!!$              stop
!!$
!!$           else 
!!$
!!$              continue
!!$
!!$           end if
!!$
!!$           W_b    = sqrt(1.d0/(1.d0 - Vx_b**2 - Vz_b**2))
!!$
!!$           rho_b  =  J_inv/(Ws*W_b*(Vs - Vx_b))
!!$
!!$           h_b    =  1.d0 + (gamma/(gamma -1.d0))*(p_b/rho_b)
!!$
!!$           Jx_b   = sigma * W_b * ( Ex_b - Vz_b*By_b - (Ex_b*Vx_b + Ez_b*Vz_b)* Vx_b) + q_b * Vx_b
!!$           Jy_b   = sigma * W_b *   Ey_b
!!$           Jz_b   = sigma * W_b * ( Ez_b + Vx_b*By_b - (Ex_b*Vx_b + Ez_b*Vz_b)* Vz_b) + q_b * Vz_b
!!$
!!$ !___________________________ Using Jeffrey definition conservation law (3.1.5) ___________________________________
!!$
!!$          q_b    = q_a - (Jx_a - Jx_b)/Vs !rho_b * W_b * ( q_a/(rho_a * W_a) + Ws * ( Jx_a - Jx_b - (q_a * Vx_a - q_b * Vx_b) ) / (J_inv) ) !
!!$
!!$
!!$!____________________________ Using LeVeque source as delta function (17.59) ________________________________________
!!$
!!$
!!$           By_b1   = By_a + (Ez_a -  Ez_b)  + (Jz_a  - Jz_b ) / (Ws *(1.d0 + Vs))
!!$
!!$           Ez_b1   = Ez_a + (By_a -  By_b1) * Vs   !Ez_a + (By_a -  By_b1 ) + (Jz_a  - Jz_b ) / (Ws *(1.d0 + Vs)) !Ez_a + (By_a -  By_b1)  / Vs + (Jz_a  - Jz_b ) / (Ws * Vs) !
!!$           Ex_b1   = Ex_a + (psi_a - psi_b) - ( (q_a - kappa  * psi_a) - (q_b - kappa * psi_b) )  / (Ws *(1.d0 + Vs) ) &
!!$                                            + (Jx_a  - Jx_b)  / (Ws *(1.d0 + Vs) )
!!$           Ey_b1   = Ey_a                   + (Jy_a  - Jy_b ) / (Ws *(1.d0 + Vs))  
!!$
!!$           ERR_BY  = abs(By_b - By_b1) / (0.5d0* (By_b + By_b1))
!!$
!!$           By_b   = By_b1 
!!$           Ez_b   = Ez_b1 
!!$           Ex_b   = Ex_b1 
!!$           Ey_b   = Ey_b1
!!$           Vx_b   = Vx_b
!!$           Vz_b   = Vz_b
!!$           rho_b  = rho_b
!!$
!!$        end do


     !////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     !RST_AlMi General case

     Ws     = 2.d1
     Vs     = sqrt(1.d0 - 1.d0/Ws**2)
     psi_a  = 0.0d0
     phi_a  = 0.0d0

     Vx_a   = 0.50d0
     Vy_a   = 0.250d0
     Vz_a   = 0.50d0
     W_a    = sqrt(1.d0/(1.d0 - Vx_a**2 - Vy_a**2- Vz_a**2))

     rho_a  = 0.1d0
     p_a    = 1.0d0
     h_a    = 1.0d0 + (gamma/(gamma -1.d0))*(p_a/rho_a)

     Bx_a   =  0.0001d0
     By_a   =  0.0001d0
     Bz_a   =  0.0001d0

     Ex_a   =  0.00163135d0 ! By_a * Vz_a - Bz_a * Vy_a !
     Ey_a   =  0.00113389d0 ! Bz_a * Vx_a - Bx_a * Vz_a ! 
     Ez_a   =  0.00158135d0  ! Bx_a * Vy_a - By_a * Vx_a !  

     Ba2    = (Bx_a**2 + By_a**2 + Bz_a**2)
     Ea2    = (Ex_a**2 + Ey_a**2 + Ez_a**2)

     b_a    = By_a/W_a
     hh_a   = h_a + ( Ba2 - Ea2)/rho_a
     ph_a   = p_a + 0.5d0 * ( Ba2 - Ea2)
     q_a    = 0.0d0

     J_inv  = Ws * rho_a * W_a * (Vs - Vx_a) ! J_inv = Ws * D * (Vs-Vx)

     EdotV_a  = (Ex_a*Vx_a + Ey_a*Vy_a + Ez_a*Vz_a)

     Jx_a   = sigma * W_a * ( Ex_a + (Bz_a * Vy_a - By_a * Vz_a) - EdotV_a * Vx_a) + q_a * Vx_a
     Jy_a   = sigma * W_a * ( Ey_a + (Bx_a * Vz_a - Bz_a * Vx_a) - EdotV_a * Vy_a) + q_a * Vy_a
     Jz_a   = sigma * W_a * ( Ez_a + (By_a * Vx_a - Bx_a * Vy_a) - EdotV_a * Vz_a) + q_a * Vz_a

     print*, "garbage"
     print*, sigma, ( (Bz_a * Vy_a - By_a * Vz_a) - EdotV_a * Vx_a)
     print*, "garbage"

     p_b    =  2.5d0
     psi_b  =  psi_a
     phi_b  =  phi_a
     !___________________________ R-H condition with out source  __________________ ___________________________________


     Bx_b   =  Bx_a
     By_b   =  By_a
     Bz_b   =  Bz_a
     Ex_b   =  Ex_a
     Ey_b   =  Ey_a
     Ez_b   =  Ez_a



     TOL1    = 1d-10
     ERR_B0 = 1.d0
     ERR_B1 = 1.d0
     ERR_B2 = 1.d0

     do  while ( ERR_B0 > TOL1 .and. ERR_B1 > TOL1 .and. ERR_B2 > TOL1 )


        Bb2    = (Bx_b**2 + By_b**2 + Bz_b**2)
        Eb2    = (Ex_b**2 + Ey_b**2 + Ez_b**2)

        ph_b   = p_b + 0.5d0 * (Bb2 - Eb2)


        Vx_b   = (hh_a * W_a * Vx_a - W_a**2         * ( Ba2 - Ea2 ) * Vx_a / (rho_a * W_a)                                 &
             -  Ws   *((ph_a-ph_b)- (Bx_a**2-Bx_b**2) + (Ey_a**2-Ey_b**2) + (Ez_a**2-Ez_b**2) ) /J_inv                    &
             + (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a))  * ( (Bz_a * Ey_a - By_a * Ez_a) - (Bz_b * Ey_b - By_b * Ez_b) ))       &
             / (hh_a * W_a        - W_a**2         * ( Ba2 - Ea2 )       / (rho_a * W_a)                                  &
             - Ws             * ( (Bz_a * Ey_a - By_a * Ez_a) - (Bz_b * Ey_b - By_b * Ez_b) ) /J_inv &
             - (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * ( (ph_a-ph_b) - (Ba2 - Bb2) ) )

        Vy_b   = (hh_a * W_a * Vy_a - W_a**2         * ( Ba2 - Ea2 ) * Vy_a / (rho_a * W_a)                                 &
             +  Ws   * ((Ex_a*Ey_a+ Bx_a*By_a)     - ( Ex_b * Ey_b + Bx_b * By_b)  ) /J_inv                               &
             + (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a))  * ( (Bx_a * Ez_a - Bz_a * Ex_a) - (Bx_b * Ez_b - Bz_b * Ex_b) ))       &
             / (hh_a * W_a        - W_a**2         * ( Ba2 - Ea2 )       / (rho_a * W_a)                                  &
             - Ws             * ( (Bz_a * Ey_a - By_a * Ez_a) - (Bz_b * Ey_b - By_b * Ez_b) ) /J_inv &
             - (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * ( (ph_a-ph_b) - (Ba2 - Bb2) ) )
        Vz_b   = (hh_a * W_a * Vz_a - W_a**2         * ( Ba2 - Ea2 ) * Vz_a / (rho_a * W_a)                                 &
             +  Ws   * ((Ex_a*Ez_a+ Bx_a*Bz_a)     - (Ex_b * Ez_b + Bx_b * Bz_b)   ) /J_inv                               &
             + (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a))  * ( (By_a * Ex_a - Bx_a * Ey_a) - (By_b * Ex_b - Bx_b * Ey_b) ))       & 
             / (hh_a * W_a        - W_a**2         * ( Ba2 - Ea2 )       / (rho_a * W_a)                                  &
             - Ws             * ( (Bz_a * Ey_a - By_a * Ez_a) - (Bz_b * Ey_b - By_b * Ez_b) ) /J_inv &
             - (Vx_a*Ws/J_inv + 1.d0/(rho_a*W_a)) * ( (ph_a-ph_b) - (Ba2 - Bb2) ) )

        Vb2    = sqrt(Vx_b**2 + Vy_b**2 +  Vz_b**2)

        if (Vb2 .gt. 1.d0) then

           write(*,*) "stop superluminal velocity, Vb2 = ", Vb2
           stop

        else 

           continue

        end if

        W_b    = sqrt(1.d0/(1.d0 - Vx_b**2 - Vy_b**2 - Vz_b**2))

        rho_b  =  J_inv/(Ws*W_b*(Vs - Vx_b))

        h_b    =  1.d0 + (gamma/(gamma -1.d0))*(p_b/rho_b)


        EdotV_b  = (Ex_b*Vx_b + Ey_b*Vy_b + Ez_b*Vz_b)

        Jx_b   = sigma * W_b * ( Ex_b + (Bz_b * Vy_b - By_b * Vz_b) - EdotV_b * Vx_b) + q_b * Vx_b
        Jy_b   = sigma * W_b * ( Ey_b + (Bx_b * Vz_b - Bz_b * Vx_b) - EdotV_b * Vy_b) + q_b * Vy_b
        Jz_b   = sigma * W_b * ( Ez_b + (By_b * Vx_b - Bx_b * Vy_b) - EdotV_b * Vz_b) + q_b * Vz_b

        !___________________________ Using Jeffrey definition conservation law (3.1.5) ___________________________________

        q_b    = q_a - (Jx_a - Jx_b)/Vs 


        !____________________________ Using LeVeque source as delta function (17.59) ________________________________________


        Bx_b1   = Bx_a + (phi_a - phi_b)  + kappa * (phi_a - phi_b) / (Ws *(1.d0 + Vs))
        By_b1   = By_a + (Ez_a -  Ez_b)  + (Jz_a  - Jz_b ) / (Ws *(1.d0 + Vs)) ! By_a + Vs * (Ez_a - Ez_b) + (Jz_a  - Jz_b )/Ws !
        Bz_b1   = Bz_a - (Ey_a -  Ey_b)  - (Jy_a  - Jy_b ) / (Ws *(1.d0 + Vs)) ! Bz_a - Vs * (Ey_a - Ey_b) - (Jy_a  - Jy_b )/Ws !

        Ex_b1   = Ex_a + (psi_a - psi_b) - ( (q_a - kappa  * psi_a) - (q_b - kappa * psi_b) )  / (Ws * (1.d0 + Vs) ) &
             + (Jx_a  - Jx_b)  / (Ws *(1.d0 + Vs) )
        Ey_b1   = Ey_a - (Bz_a -  Bz_b1) * Vs   
        Ez_b1   = Ez_a + (By_a -  By_b1) * Vs   

!!$           Ex_b   =  By_b * Vz_b - Bz_b * Vy_b !0.1d0 
!!$           Ey_b   =  Bz_b * Vx_b - Bx_b * Vz_b !0.1d0
!!$           Ez_b   =  Bx_b * Vy_b - By_b * Vx_b !0.1d0 



        ERR_B0  = abs(By_b - By_b1) / (0.5d0* (By_b + By_b1))
        ERR_B1  = abs(By_b - By_b1) / (0.5d0* (By_b + By_b1))
        ERR_B2  = abs(By_b - By_b1) / (0.5d0* (By_b + By_b1))

        Bx_b   = Bx_b1 
        By_b   = By_b1 
        Bz_b   = Bz_b1 
        Ez_b   = Ez_b1 
        Ex_b   = Ex_b1 
        Ey_b   = Ey_b1
        Vx_b   = Vx_b
        Vz_b   = Vz_b
        rho_b  = rho_b

     end do


     !****************************** ANALITICAL SOLUTION ***********************************************************************

     do i=0,imax

        pos_shock = 0.5d0 + Vs * Delt * hmax

        if (i < floor(pos_shock/Delx)) then

           psi(i,1,1) = psi_a
           phi(i,1,1) = phi_a
           Bx(i,1,1)  = Bx_a
           By(i,1,1)  = By_a
           Bz(i,1,1)  = Bz_a
           Vx(i,1,1)  = Vx_a
           Vy(i,1,1)  = Vy_a
           Vz(i,1,1)  = Vz_a
           Ex(i,1,1)  = Ex_a
           Ey(i,1,1)  = Ey_a
           Ez(i,1,1)  = Ez_a
           q(i,1,1)   = q_a
           p(i,1,1)   = p_a
           rho(i,1,1) = rho_a
        else
           psi(i,1,1) = psi_b
           phi(i,1,1) = phi_b
           Bx(i,1,1)  = Bx_b
           By(i,1,1)  = By_b
           Bz(i,1,1)  = Bz_b
           Vx(i,1,1)  = Vx_b
           Vy(i,1,1)  = Vy_b
           Vz(i,1,1)  = Vz_b
           Ex(i,1,1)  = Ex_b
           Ey(i,1,1)  = Ey_b
           Ez(i,1,1)  = Ez_b
           q(i,1,1)   = q_b
           p(i,1,1)   = p_b
           rho(i,1,1) = rho_b

        end if


        !_________________________________________________________________________________________________________________________

        write (40,130)  i * Delx, Bx (i,1,1) 
        write (41,130)  i * Delx, By (i,1,1) 
        write (42,130)  i * Delx, Bz (i,1,1) 
        write (43,130)  i * Delx, Vx (i,1,1) 
        write (44,130)  i * Delx, Vy (i,1,1) 
        write (45,130)  i * Delx, Vz (i,1,1) 
        write (46,130)  i * Delx, Ex (i,1,1) 
        write (47,130)  i * Delx, Ey (i,1,1) 
        write (48,130)  i * Delx, Ez (i,1,1) 
        write (49,130)  i * Delx, p  (i,1,1) 
        write (50,130)  i * Delx, rho(i,1,1)
        write (51,130)  i * Delx, q  (i,1,1)


        W(i,1,1)       = 1.d0/sqrt(1.d0-(Vx (i,1,1)**2 + Vz (i,1,1)**2))

        ph(i,1,1)      = p(i,1,1) + 0.5d0 * (By(i,1,1)**2 - (Ex(i,1,1)**2 +Ey(i,1,1)**2 +Ez(i,1,1)**2) )

        epsiln(i,1,1) = p(i,1,1)/((gamma-1.d0)*rho(i,1,1))
        enthpy(i,1,1)  = rho(i,1,1) * (1.d0+epsiln(i,1,1)) + p(i,1,1) 

        hhat(i,1,1)    = ( enthpy(i,1,1) + By (i,1,1)**2/W(i,1,1)**2 ) / rho(i,1,1)

        write (52,130)  i * Delx, ph  (i,1,1) 
        write (53,130)  i * Delx, hhat(i,1,1) 
        !_________________________________________________________________________________________________________________________



     end do   !i

     !************************************************************************************************************************


     do i=0,imax


        if (i < 0.5 * imax) then
           psi(i,1,1) = psi_a
           phi(i,1,1) = phi_a
           Bx(i,1,1)  = Bx_a
           By(i,1,1)  = By_a
           Bz(i,1,1)  = Bz_a
           Vx(i,1,1)  = Vx_a
           Vy(i,1,1)  = Vy_a
           Vz(i,1,1)  = Vz_a
           Ex(i,1,1)  = Ex_a
           Ey(i,1,1)  = Ey_a
           Ez(i,1,1)  = Ez_a
           q(i,1,1)   = q_a
           p(i,1,1)   = p_a
           rho(i,1,1) = rho_a
        else
           psi(i,1,1) = psi_b
           phi(i,1,1) = phi_b
           Bx(i,1,1)  = Bx_b
           By(i,1,1)  = By_b
           Bz(i,1,1)  = Bz_b
           Vx(i,1,1)  = Vx_b
           Vy(i,1,1)  = Vy_b
           Vz(i,1,1)  = Vz_b
           Ex(i,1,1)  = Ex_b
           Ey(i,1,1)  = Ey_b
           Ez(i,1,1)  = Ez_b
           q(i,1,1)   = q_b
           p(i,1,1)   = p_b
           rho(i,1,1) = rho_b

        end if
        !_________________________________________________________________________________________________________________________

        write (20,130)  i * Delx, By(i,1,1)
        !_________________________________________________________________________________________________________________________


     end do   !i

  else if (DIM == 1 .and. TEST == 6) then

     ! Reconnection 1D test
     ! Aloy and Mimica (2010)

     do i=0,imax


        posx   = - 0.5d0 * Ly + i * Delx
        posy   =   0.0d0 

        r      =   sqrt(1 -(posy**2/Lr**2))

        if (SLC == 0)  sigma = sigma_0

        psi(i,1,1) = 0.d0
        phi(i,1,1) = 0.d0
        Bx(i,1,1)  =  (B0*posy*delta_rec/Lr**2) * (1.d0+log(cosh(posx/delta_rec))/r) !& 
        By(i,1,1)  =   B0 * r * tanh(posx/delta_rec) 
        Bz(i,1,1)  = 0.d0
        Vx(i,1,1)  = 0.d0
        Vy(i,1,1)  = 0.d0
        Vz(i,1,1)  = 0.d0
        Ex(i,1,1)  = -(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))
        Ey(i,1,1)  = -(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
        Ez(i,1,1)  = -(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
        q(i,1,1)   = 0.d0
        rho(i,1,1) = rho_up
        p  (i,1,1) = P0 * r**2 / (cosh(posx/delta_rec))**2 + P_up


     end do   !i

  else if (DIM == 1 .and. TEST == 7) then

     ! TM V1 Tomek (2016)

     do i=0,imax


        posx   =   - 0.5d0 * Longx + i * Delx
        posy   =   0.0d0 

        if (SLC == 0)  sigma = sigma_0

!!$           psi(i,1,1) = 0.d0
!!$           phi(i,1,1) = 0.d0
!!$           Bx(i,1,1)  = 0.d0
!!$           By(i,1,1)  = B0 * tanh(posx/delta_rec) 
!!$           Bz(i,1,1)  = B0 / cosh(posx/delta_rec) 
!!$           Vx(i,1,1)  = 0.d0
!!$           Vy(i,1,1)  = 0.d0
!!$           Vz(i,1,1)  = 0.d0
!!$           Ex(i,1,1)  = 0.d0
!!$           Ey(i,1,1)  = 0.d0
!!$           Ez(i,1,1)  = 0.d0 
!!$           q(i,1,1)   = 0.d0
!!$           rho(i,1,1) = rho_up
!!$           p  (i,1,1) = P0


        psi(i,1,1) = 0.d0
        phi(i,1,1) = 0.d0


        Bx(i,1,1)  = epsilon_tm * B0 / cosh(posx / a_tm)
        By(i,1,1)  = B0 * tanh(posx / a_tm) !+ epsilon_tm * B0 * &
        !sin(k_tm * posy) * tanh(posx / a_tm) / (cosh(posx / a_tm) * kap_tm)
        Bz(i,1,1)  = B0 / cosh(posx / a_tm)

        Vx(i,1,1)  = 0.d0
        Vy(i,1,1)  = 0.d0
        Vz(i,1,1)  = 0.d0

        Ex(i,1,1)  = 0.d0
        Ey(i,1,1)  = 0.d0
        Ez(i,1,1)  = 0.d0

        q(i,1,1)   = 0.d0
        rho(i,1,1) = rho_up
        p  (i,1,1) = P0 



        !_________________________________________________________________________________________________________________________

        write (20,130)  - 0.5d0 * Longx + i * Delx, Bx(i,1,1)
        write (21,130)  - 0.5d0 * Longx + i * Delx, By(i,1,1)
        write (22,130)  - 0.5d0 * Longx + i * Delx, Bz(i,1,1)
        !_________________________________________________________________________________________________________________________


     end do   !i


  else if (DIM == 1 .and. TEST == 8) then

     ! 1D MD

     sum_b20   = 0.d0

     do i=0,imax


        posx   =   - 0.5d0 * Longx + i * Delx

        if (SLC == 0)  sigma = sigma_0

        psi(i,1,1) = 0.d0
        phi(i,1,1) = 0.d0
        Bx(i,1,1)  = 0.d0
        By(i,1,1)  = B0 * sin(k_x * posx)
        Bz(i,1,1)  = B0 * cos(k_x * posx)
        Vx(i,1,1)  = V0
        Vy(i,1,1)  = 0.d0
        Vz(i,1,1)  = 0.d0
        Ex(i,1,1)  = 0.d0
        Ey(i,1,1)  = 0.d0
        Ez(i,1,1)  = 0.d0
        q(i,1,1)   = 0.d0
        rho(i,1,1) = rho_up
        p  (i,1,1) = P_up


        sum_b20   = sum_b20 + (By(i,1,1)**2 + Bz(i,1,1)**2) * Delx 
        !_________________________________________________________________________________________________________________________

        write (20,130)  - 0.5d0 * Longx + i * Delx, Vx(i,1,1)
        write (21,130)  - 0.5d0 * Longx + i * Delx, By(i,1,1)
        write (22,130)  - 0.5d0 * Longx + i * Delx, Bz(i,1,1)
        !_________________________________________________________________________________________________________________________


     end do   !i

  else if (DIM == 1 .and. TEST == 9) then

     ! 1D SL

     do i=0,imax

        posx   =   - 0.5d0 * Longx + i * Delx

        if (SLC == 0)  sigma = sigma_0

        psi(i,1,1) = 0.d0
        phi(i,1,1) = 0.d0
        Bx(i,1,1)  = 0.d0
        By(i,1,1)  = 0.d0
        Bz(i,1,1)  = 0.d0
        Vx(i,1,1)  = 0.d0
        Vy(i,1,1)  = V0 * sin(k_x * posx) 
        Vz(i,1,1)  = 0.d0
        Ex(i,1,1)  = 0.d0
        Ey(i,1,1)  = 0.d0
        Ez(i,1,1)  = 0.d0 
        q(i,1,1)   = 0.d0
        rho(i,1,1) = rho_up
        p  (i,1,1) = P_up

        !_________________________________________________________________________________________________________________________

        write (20,130)  - 0.5d0 * Longx + i * Delx, Vx(i,1,1)
        write (21,130)  - 0.5d0 * Longx + i * Delx, Vy(i,1,1)
        write (22,130)  - 0.5d0 * Longx + i * Delx, Vz(i,1,1)
        !_________________________________________________________________________________________________________________________



     end do   !i

  else if (DIM == 1 .and. TEST == 10) then

     ! CPAW-2 TEST

     do i=0,imax


        posx = - 0.5d0 + i * Delx

        psi(i,1,1) = 0.d0
        phi(i,1,1) = 0.d0
        q(i,1,1)   = 0.d0
        p(i,1,1)   = P_up
        rho(i,1,1) = rho_up


        epsi_A     = P_up    / ((gamma-1.d0)    * rho_up )
        enthalpy_A = rho_up * ( 1.d0 + epsi_A) + P_up


        VA = sqrt(2.d0 * B0**2 /    ( ( enthalpy_A + B0**2 * (1.d0 + etha_A**2)) * &
             ( 1.d0 + sqrt(1.d0-(2.d0 * etha_A * B0**2 / (enthalpy_A + B0**2 * (1.d0 + etha_A**2)))**2))))

        !VA = B0 / sqrt(rho_up) ! classical expression


        Bx(i,1,1)  =   B0
        By(i,1,1)  =   B0 * etha_A * cos(k_x * posx)
        Bz(i,1,1)  =   B0 * etha_A * sin(k_x * posx)
        Vx(i,1,1)  =   0.d0
        Vy(i,1,1)  = - VA * etha_A * cos(k_x * posx)
        Vz(i,1,1)  = - VA * etha_A * sin(k_x * posx)
        Ex(i,1,1)  =   0.d0           !-(Vy(i,1,1)*Bz(i,1,1)-By(i,1,1)*Vz(i,1,1))   !
        Ey(i,1,1)  =   VA * Bz(i,1,1) !-(Vz(i,1,1)*Bx(i,1,1)-Bz(i,1,1)*Vx(i,1,1))   !
        Ez(i,1,1)  = - VA * By(i,1,1) !-(Vx(i,1,1)*By(i,1,1)-Bx(i,1,1)*Vy(i,1,1))   !

        !_________________________________________________________________________________________________________________________

        write (20,130)  - 0.5d0 * Longx + i * Delx, By(i,1,1)
        write (21,130)  - 0.5d0 * Longx + i * Delx, Bz(i,1,1)
        write (22,130)  - 0.5d0 * Longx + i * Delx, Vy(i,1,1)
        write (23,130)  - 0.5d0 * Longx + i * Delx, Vz(i,1,1)
        !_________________________________________________________________________________________________________________________


     end do   !i



  else if (DIM == 2 .and. TEST == 20) then

     ! CE TEST

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx = - 6.d0 + i * Delx
           posy = - 6.d0 + j * Dely

           r   =   sqrt(posx**2 + posy**2)

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  = 0.05d0 ! 0.1d0 !1.d0 ! 
           By(i,j,1)  = 0.d0
           Bz(i,j,1)  = 0.d0
           Ex(i,j,1)  = 0.d0
           Ey(i,j,1)  = 0.d0
           Ez(i,j,1)  = 0.d0
           Vx(i,j,1)  = 0.d0
           Vy(i,j,1)  = 0.d0
           Vz(i,j,1)  = 0.d0
           q(i,j,1)   = 0.d0


           if ( r .le. 0.8d0) then

              rho(i,j,1)   = 1.d-2
              p  (i,j,1)   = 1.d0

           else if ( r .gt. 0.8d0 .and. r .le. 1.d0) then

              rho(i,j,1)   = 1.d-2  * exp(- beta  * (r-0.8d0))
              p  (i,j,1)   = 1.d0   * exp(- alpha * (r-0.8d0))

           else
              rho(i,j,1)   =  1.d-3 !Palenzuela \\! 1.d-4 ! Komissarov \\!
              p  (i,j,1)   =  1.d-3 !Palenzuela \\! 3.d-5 ! Komissarov \\!1.d-5 ! Komissarov modificado \\!

           end if

           !        end do !k
        end do  !j
     end do   !i

  else if (DIM == 2 .and. TEST == 21) then

     ! Rotor TEST

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx  = - 0.5d0 + i * Delx
           posy  = - 0.5d0 + j * Dely

           r     =   sqrt(posx**2 + posy**2)

           !---------------------------------------------
           ! OMEGA CONDITIONS works fine with
           !
           ! Grid 100X100 :All schemes \Omega 8.5
           !
           ! Grid 200X200 :
           !
           ! SSP2(2,2,2)   omega = 7.0
           ! SSP2(3,3,2)   omega = 7.5
           ! SSP2H(3,3,2)  omega = 6.5
           ! DP2-A(2,4,2)  omega = 6.0
           !
           ! Grid 400X400 : All schemes \Omega 6.0; its possible to run with \Omega 8.5 if the rotation is counterclockwise and CFL=0.1
           ! Grid 400X400 : All schemes \Omega 8.5  witout EGLM CFL =0.1

           omega = 8.5d0

           !---------------------------------------------

           if ( r .le. 0.079d0) then

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              Bx(i,j,1)  = 1.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Vx(i,j,1)  =  - posy  * omega  !original  (-)
              Vy(i,j,1)  =    posx  * omega  !original  (+)
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
              q(i,j,1)   = 0.d0
              rho(i,j,1) = 1.d1
              p  (i,j,1) = 1.d0

           else if ( r .gt. 0.079d0 .and. r .le. 0.1d0) then

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              Bx(i,j,1)  = 1.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Vx(i,j,1)  =  - posy  * omega  !* exp(- beta  * (r-0.079d0))  !original  (-)
              Vy(i,j,1)  =    posx  * omega  !* exp(- beta  * (r-0.079d0))  !original  (+)
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
              q(i,j,1)   = 0.d0
              rho(i,j,1) = 1.d1 !* exp(- alpha * (r-0.079d0))
              p  (i,j,1) = 1.d0


           else

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              Bx(i,j,1)  = 1.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Vx(i,j,1)  = 0.d0
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1)) !0.d0
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1)) !0.d0
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1)) !0.d0
              q(i,j,1)   = 0.d0
              rho(i,j,1) = 1.d0
              p  (i,j,1) = 1.d0

           end if

           !        end do !k
        end do  !j
     end do   !i

  else if (DIM == 2 .and. TEST == 22) then

     ! Cylindrical Star

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx  = - 3.0d0 + i * Delx
           posy  = - 3.0d0 + j * Dely

           r     =   sqrt(posx**2 + posy**2)

           omega = 0.1d0

!!$           if ( r .le. 0.7d0) then

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           rho(i,j,1) = 1.d0 * exp(-(r/0.7d0)**2)
           p  (i,j,1) = rho(i,j,1)**gamma
           Bx(i,j,1)  = 0.d0
           By(i,j,1)  = 0.d0
           Bz(i,j,1)  = 2.d0* rho(i,j,1) * 0.05d0 * (1.d0 - (r/0.7d0)**2)
           Vx(i,j,1)  = rho(i,j,1) * omega * posy
           Vy(i,j,1)  =-rho(i,j,1) * omega * posx
           Vz(i,j,1)  = 0.d0
           Ex(i,j,1)  =-(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
           Ey(i,j,1)  =-(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
           Ez(i,j,1)  =-(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))

           if (i == 0 .or. j== 0) then
              q(i,j,1)   = 0.d0
           else
              q(i,j,1)   = (Ex(i,j,1) - Ex(i-1,j,1)) / Delx + (Ey(i,j,1) - Ey(i,j-1,1)) / Dely &
                   + (Ez(i,j,1) - Ez(i,j  ,1)) / Delz
           end if

!!$           else
!!$
!!$              psi(i,j,1) = 0.d0
!!$              phi(i,j,1) = 0.d0
!!$              rho(i,j,1) = 0.01d0
!!$              p  (i,j,1) = rho(i,j,1)**gamma
!!$              Bx(i,j,1)  = 0.d0
!!$              By(i,j,1)  = 0.d0
!!$              Bz(i,j,1)  = 2.d0* rho(i,j,1) * 0.05d0 * (1.d0 - (r/0.7d0)**2)
!!$              Vx(i,j,1)  = 0.d0
!!$              Vy(i,j,1)  = 0.d0
!!$              Vz(i,j,1)  = 0.d0
!!$              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
!!$              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
!!$              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
!!$              q(i,j,1)   = 0.d0
!!$
!!$           end if

           !           if (j == 0.5 * jmax) print*, Vx(i,j,1), i


           !        end do !k
        end do  !j
     end do   !i

  else if (DIM == 2 .and. TEST == 23) then

     ! Slab Jet TEST (SJ)

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx  = i * Delx
           posy  = - 14.d0 + j * Dely

           r     =   sqrt(posx**2 + posy**2)

           if ( i .le. 1 .and. posy .gt. -1.d0 .and. posy .le. 1.d0) then

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              Bx(i,j,1)  = 1.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Vx(i,j,1)  = 0.998752d0
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
              q(i,j,1)   = 0.d0
              rho(i,j,1) = 1.d-1
              p  (i,j,1) = 1.d-2

           else

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              Bx(i,j,1)  = 1.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Vx(i,j,1)  = 0.d0
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
              q(i,j,1)   = 0.d0
              rho(i,j,1) = 1.d1
              p  (i,j,1) = 1.d-2

           end if

           !        end do !k
        end do  !j
     end do   !i

  else if (DIM == 2 .and. TEST == 24) then

     ! Relativistic Shock Cloud Interaction (CSI), Mignone and Bodo (2006)

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx = i * Delx
           posy = j * Dely



           if ( posx .le. 0.6d0) then

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              Bx(i,j,1)  = 0.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  =-2.12971d0 
              Vx(i,j,1)  = 0.d0
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
              q(i,j,1)   = 0.d0
              rho(i,j,1) = 42.5942d0
              p  (i,j,1) = 127.9483d0

           else if ( posx .gt. 0.6d0) then

              r   =   sqrt((posx-0.8d0)**2 + (posy-0.5)**2)

              if ( r .le. 0.15d0) then

                 rho(i,j,1) = 10.d0

              else


                 psi(i,j,1) = 0.d0
                 phi(i,j,1) = 0.d0
                 Bx(i,j,1)  = 0.d0
                 By(i,j,1)  = 0.d0
                 Bz(i,j,1)  = 0.5d0
                 Vx(i,j,1)  =-0.994987
                 Vy(i,j,1)  = 0.d0
                 Vz(i,j,1)  = 0.d0
                 Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
                 Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
                 Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
                 q(i,j,1)   = 0.d0
                 rho(i,j,1) = 1.d0
                 p  (i,j,1) = 1.d-3 



              end if

           end if

           if (j == 0) then
              write(32,*) ''
           end if

           write (32,131)  i * Delx, j * Dely, rho(i,j,1)


           !        end do !k
        end do  !j
     end do   !i

     !************************************************************

  else if (DIM == 2 .and. TEST == 25) then

     ! RMI TEST

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx = i * Delx
           posy = j * Dely

           Ws     = 1.5d0
           Vs     = sqrt(1.d0-1.d0/Ws**2)
           rho_a  = 1.d0
           W_a    = 1.d0
           p_a    = 0.1d0
           p_b    = 1.35d0
           Vx_a   = 0.d0
           J_inv  = Ws * rho_a * W_a * (Vs - Vx_a) ! J_inv = Ws * D * (Vs-Vx)

           eps_a  = p_a/((gamma-1.d0)*rho_a)
           h_a    = ( 1.d0 + eps_a ) + p_a / rho_a

           Vx_b   = (h_a * W_a * Vx_a + Ws * (p_b - p_a)/J_inv)/(h_a * W_a + (p_b - p_a)*(Ws * Vx_a/J_inv + 1/(rho_a * W_a)))

           W_b    = sqrt(1.d0 / (1.d0 - Vx_b**2))

           rho_b  = J_inv/(Ws * W_b * (Vs - Vx_b))

           frontera = x_0 + a_rmi * sin(pi/2 + 2*pi * posy / lambda)



           if ( posx .lt. 1.d0) then

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              q(i,j,1)   = 0.d0
              rho(i,j,1) = rho_b
              p  (i,j,1) = p_b
              Bx(i,j,1)  = 0.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0 ! Its possible introduce a magnetic field ussing the invariant f= Ws B (Vs-Vx)
              Vx(i,j,1)  = Vx_b
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0 ! this component of velocity must be zero for use the Riemman problem as Romero et. al (2005)
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))


           else if ( posx .eq. 1.d0) then

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              q(i,j,1)   = 0.d0
              rho(i,j,1) = rho_a
              p  (i,j,1) = p_a
              Bx(i,j,1)  = 0.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Vx(i,j,1)  = Vs
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))

           else if ( posx .le. frontera  ) then

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              q(i,j,1)   = 0.d0
              rho(i,j,1) = rho_a
              p  (i,j,1) = p_a
              Bx(i,j,1)  = 0.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Vx(i,j,1)  = 0.d0
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))


           else

              psi(i,j,1) = 0.d0
              phi(i,j,1) = 0.d0
              q(i,j,1)   = 0.d0
              rho(i,j,1) = 35.d0
              p  (i,j,1) = 0.1d0
              Bx(i,j,1)  = 0.d0
              By(i,j,1)  = 0.d0
              Bz(i,j,1)  = 0.d0
              Vx(i,j,1)  = 0.d0
              Vy(i,j,1)  = 0.d0
              Vz(i,j,1)  = 0.d0
              Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
              Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
              Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))


           end if

           if (j == 0) then
              write(32,*) ''
           end if

           write (32,131)  i * Delx, j * Dely, rho(i,j,1)


           !        end do !k
        end do  !j
     end do   !i

  else if (DIM == 2 .and. TEST == 26) then

     !Magnetic Diffusion (MD 2D)


     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx   = - 0.5d0 * Lx + i * Delx
           posy  =  - 0.5d0 * Lx + j * Dely


           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx (i,j,1) =- B0 * sin(k_x * posx + k_y * posy) / sqrt(2.d0)
           By (i,j,1) =  B0 * sin(k_x * posx + k_y * posy) / sqrt(2.d0)
           Bz (i,j,1) =  B0 * cos(k_x * posx + k_y * posy)
           Vx (i,j,1) = 0.d0
           Vy (i,j,1) = 0.d0
           Vz (i,j,1) = 0.d0
           Ex (i,j,1) = 0.d0
           Ey (i,j,1) = 0.d0
           Ez (i,j,1) = 0.d0
           q  (i,j,1) = 0.d0
           rho(i,j,1) = rho_up
           p  (i,j,1) = P_up

           if (j == jmed) then
              !_________________________________________________________________________________________________________________________

              write (120,130)  - 0.5d0 * Lx + i * Delx, Bx(i,jmed,1)
              write (121,130)  - 0.5d0 * Lx + i * Delx, By(i,jmed,1)
              write (122,130)  - 0.5d0 * Lx + i * Delx, Bz(i,jmed,1)
              !_________________________________________________________________________________________________________________________

           end if

           !        end do !k
        end do  !j
     end do   !i

     !*****************************************************************************

  else if (DIM == 2 .and. TEST == 27) then

     !Shear Layer (SL 2D)


     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx   = - 0.5d0 * Lx + i * Delx
           posy  =  - 0.5d0 * Lx + j * Dely


           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx (i,j,1) = 0.d0
           By (i,j,1) = 0.d0
           Bz (i,j,1) = 0.d0
           Vx (i,j,1) = V0 * sin(k_x * posx + k_y * posy) / sqrt(2.d0)
           Vy (i,j,1) =-V0 * sin(k_x * posx + k_y * posy) / sqrt(2.d0)
           Vz (i,j,1) = 0.d0
           Ex (i,j,1) = 0.d0
           Ey (i,j,1) = 0.d0
           Ez (i,j,1) = 0.d0
           q  (i,j,1) = 0.d0
           rho(i,j,1) = rho_up
           p  (i,j,1) = P_up

           if (j == jmed) then
              !_________________________________________________________________________________________________________________________

              write (120,130)  - 0.5d0 * Lx + i * Delx, Vx(i,jmed,1)
              write (121,130)  - 0.5d0 * Lx + i * Delx, Vy(i,jmed,1)
              !_________________________________________________________________________________________________________________________

           end if

           !        end do !k
        end do  !j
     end do   !i


     !*****************************************************************************

  else if (DIM == 2 .and. TEST == 30) then

     ! Reconecction TEST Aloy & Mimica (2009)

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx   = - 0.5d0 * Lx + i * Delx
           posy   = - 0.5d0 * Ly + j * Dely

           r      =   sqrt(1.d0 - (posy**2/Ly**2))

           if (SLC == 1) sigma_loc(i,j,1) = sigma_0 * sigma_1 * (cosh(4.d0 * sqrt(posx**2+posy**2)/delta_rec))**2 / &
                (sigma_1 * (cosh(4.d0 * sqrt(posx**2+posy**2) /delta_rec))**2 + (sigma_0 - sigma_1) )

           if (SLC == 1)  sigma = sigma_loc(i,j,1) 

           if (SLC == 0)  sigma = sigma_0

           ! Aloy and Mimica (2010)

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  =  (B0*posy*delta_rec/Ly**2) * (1.d0+(log(abs(cosh(posx/delta_rec)))/r)) !& 
           By(i,j,1)  =   B0 * r * tanh(posx/delta_rec) 
           Bz(i,j,1)  = 0.d0
           Vx(i,j,1)  = 0.d0
           Vy(i,j,1)  = 0.d0
           Vz(i,j,1)  = 0.d0
           Ex(i,j,1)  = 0.d0 !-(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
           Ey(i,j,1)  = 0.d0 !-(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
           Ez(i,j,1)  = 0.d0 !-(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
           q(i,j,1)   = 0.d0
           rho(i,j,1) = rho_up
           p  (i,j,1) = P0 * r**2 / (cosh(posx/delta_rec))**2 + P_up 


           ! Aloy and Mimica (2010) with current density Jz


!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0
!!$           Bx(i,j,1)  =   B0 * r * tanh(posy/delta_rec) 
!!$           By(i,j,1)  =  (B0 * posx * delta_rec/Lr**2) * (1.d0+log(cosh(posy/delta_rec))/r) 
!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0 !-(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
!!$           Ey(i,j,1)  = 0.d0 !-(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
!!$           Ez(i,j,1)  = B0*delta_rec*(1.d0+log(cosh(posy/delta_rec))/r+posx**2*log(cosh(posy/delta_rec))/(Lr**2*r**3)) &
!!$                      / ( sigma * Lr**2)                                                                  &
!!$                      - B0 * r / (sigma * delta_rec * (cosh(posy/delta_rec))**2)  !-(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = rho_up
!!$
!!$           p  (i,j,1) = P0 * ( -( r * tanh(posy/delta_rec) )**2  &
!!$                      - ( (posx * delta_rec/Lr**2) * (1.d0+log(cosh(posy/delta_rec))/r) )**2 ) + P_up


!!$           p  (i,j,1) = P0 * ( -( r * tanh(posy/delta_rec) )**2  &
!!$                      - ( (posx * delta_rec/Lr**2) * (1.d0+log(cosh(posy/delta_rec))/r) )**2 ) + P_up

!!$           p  (i,j,1) = 2.d0  * P0 * ( - ( r * tanh(posy/delta_rec))**2                                                             &
!!$                                   - ( (posx * delta_rec/Lr**2) * (1.d0+log(cosh(posy/delta_rec))/r) )**2                           &
!!$                                   - (  delta_rec*(1.d0+log(cosh(posy/delta_rec))/r+posx**2*log(cosh(posy/delta_rec))/(Lr**2*r**3)) &
!!$                                   / ( sigma * Lr**2) )**2                                                                          &
!!$                                   + ( r / (sigma * delta_rec * (cosh(posy/delta_rec))**2) )**2 )  

!!$           p  (i,j,1) = P0 * (delta_rec**3 * log(cosh(posy/delta_rec)) * (2.d0 * r**3 + log(cosh(posy/delta_rec)))         &
!!$                           +   Lr**2 * r**4 * delta_rec /  (cosh(posy/delta_rec))**2 ) / (2.d0 * Lr**2 * r**2 * delta_rec) &
!!$                      - P0 * ( Lr**2 * delta_rec**2 * (log(cosh(posy/delta_rec)))**2/(2.d0 * Lr**2 * r**2)                 &
!!$                           +  0.25d0 * posx**2 * (delta_rec**2 * (1.d0 + cosh( 2.d0 * posy/delta_rec))                     &
!!$                           -  2.d0   * Lr**2 * log(cosh(posy/delta_rec))) / (cosh(posy/delta_rec))**2                      &
!!$                           + (2.d0   * Lr**2 * r**2 + 3.d0 * posx**2 * delta_rec**2 * log(cosh(posy/delta_rec))            &
!!$                           +  3.d0   * posx**2 * delta_rec**2 * cosh( 2.d0 * posy/delta_rec) * log(cosh(posy/delta_rec)) ) &
!!$                           / (6.d0   * (cosh(posy/delta_rec))**2 ) ) / Lr**4 + P_up

!!$           p  (i,j,1) = - P0 * ( posx**2 * delta_rec**3 * log(cosh(posy/delta_rec)) * (2.d0 * r + log(cosh(posy/delta_rec)))   &
!!$                        -  Lr**4 * r**4 * delta_rec /(cosh(posy/delta_rec))**2 ) / (Lr**4 * r**2 * delta_rec)                  &
!!$                        -  P0 * ( (delta_rec**2 * log(cosh(posy/delta_rec))                                                    &
!!$                        * (2.d0 * posx**2 * r + Lr**2 * log(cosh(posy/delta_rec)) ) )/r**2                                     &
!!$                        - 0.5d0 * posx**2 * ((Lr**2 - delta_rec**2) * cosh( posy/delta_rec)                                    &
!!$                        - (Lr**2 + delta_rec**2)/(cosh(posy/delta_rec))**2) )  / Lr**4  + P_up                 


           !           p  (i,j,1) = P0 * r**2 / (cosh(posy/delta_rec))**2 + P_up
!!$           p  (i,j,1) = 0.5d0 * P0 * delta_rec**2 * log(cosh(posy/delta_rec)) * (2.d0 * r**3 + log(cosh(posy/delta_rec)) ) &
!!$                      / (Lr**2 * r**2)                                                                                        &
!!$                      - P0 * ( THIRD * r**3 + 0.5d0 * r**2 * log(cosh(posy/delta_rec))) / (cosh(posy/delta_rec))**2  &
!!$                      - P0 * ( posx * delta_rec * (1.d0 + log(cosh(posy/delta_rec))/r)/Lr**2 + r * tanh(posy/delta_rec)) + P_up
!!$                     - B0 * delta_rec * (1.d0 + ( posx**2/(Lr**2 * r**3) + 1.d0/r) * log(cosh(posy/delta_rec)))              &
!!$                      + B0 * r / (delta_rec * sigma * (cosh(posy/delta_rec))**2)                                              &

           ! Takahashi and Ohsuga (2013)


!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0
!!$           Bx(i,j,1)  = 2.d0 * epsilon_tm * B0 * posy * exp(-(posx**2 + posy**2)/delta_rec**2)/delta_rec
!!$           By(i,j,1)  = B0 * tanh(posx/delta_rec) &
!!$                      - 2.d0 * epsilon_tm * B0 * posx * exp(-(posx**2 + posy**2)/delta_rec**2)/delta_rec
!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = 0.d0
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = P0 / (cosh(posx/delta_rec))**2 + rho_up
!!$           p  (i,j,1) = P0 / (cosh(posx/delta_rec))**2 + P_up !P_up


           ! Zenitani et al. (2010)

!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0
!!$           Bx(i,j,1)  = 0.03d0 * B0 * posy * exp(-( posx**2 + posy**2 )/4.d0) 
!!$           By(i,j,1)  =  B0 * tanh(posx) &
!!$                       - 0.03d0 * B0 * posx * exp(-( posx**2 + posy**2 )/4.d0) 
!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = B0 / (sigma * (cosh(posx))**2)
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = P0 / (cosh(posx))**2 + rho_up
!!$           p  (i,j,1) = P0 / (cosh(posx))**2 + P_up


!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0
!!$           Bx(i,j,1)  = B0 * tanh(2.d0 * posy/delta_rec) &
!!$                           + 2.d0 * 0.03d0 * B0 * posy * exp(-( posx**2 + posy**2 )/delta_rec**2)/delta_rec ! 
!!$           By(i,j,1)  =    - 2.d0 * 0.03d0 * B0 * posx * exp(-( posx**2 + posy**2 )/delta_rec**2)/delta_rec ! 
!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = - 2.d0 * B0 / (delta_rec * (cosh(2.d0 * posy/delta_rec))**2)
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = P0 * (tanh(2.d0 * posy/delta_rec))**2 / theta + rho_up
!!$           p  (i,j,1) =-P0 * (tanh(2.d0 * posy/delta_rec))**2 + P_up



!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0
!!$           Bx(i,j,1)  = B0 * tanh(posy/delta_rec) - 2.d0 * B0 * posy * exp(-( posx**2 + posy**2 )/delta_rec)/delta_rec ! in Aloy & Mimica * r 
!!$           By(i,j,1)  = 2.d0 * B0 * posx * exp(-( posx**2 + posy**2 )/delta_rec)/delta_rec ! (B0*posx*delta_rec/Lr**2) * (1+log(cosh(posy/delta_rec))/r) + ... or B0 * r / cosh(posy/delta_rec) +  ...
!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
!!$           Ey(i,j,1)  = -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
!!$           Ez(i,j,1)  = -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = - B0 / (sigma * delta_rec**2 * (cosh(posy/delta_rec))**2) !& 
!!$!                        - 4.d0 * B0 * (posx**2 + posy**2)  * exp(-(posx**2+posy**2)/delta_rec)/ (sigma * delta_rec**2)
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = rho_up
!!$           p  (i,j,1) = P0     / (cosh(posy/delta_rec))**2 + P_up ! in Aloy & Mimica * r**2 


           ! Aloy and Zenitani  mixed


!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0

!!$           Bx(i,j,1)  = B0 * r * tanh(posy/delta_rec) - 2.d0 * B0 * posy * exp(-(posx**2 + posy**2 ))
!!$           By(i,j,1)  = 2.d0 * B0 * posx * exp(-( posx**2 + posy**2 ))

!!$           Bx(i,j,1)  = B0 * r * tanh(posy/delta_rec) - 2.d0 * B0 * r**2 * posy * exp(- (posx**2 + posy**2 ))
!!$           By(i,j,1)  = (B0*posx*delta_rec/Lr**2) * (1.d0+log(cosh(posy/delta_rec))/r) &
!!$                        + 2.d0 * B0 * posx * (1.d0 / Lr**2 + r**2) * exp(-( posx**2 + posy**2 ))

!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = B0 * r / (sigma * delta_rec * cosh(posy/delta_rec)**2)!& 
!!$                        !+ 0.015d0 * B0 * (posx**2 + posy**2)  * exp(-0.25d0*(posx**2+posy**2)) / sigma_loc(i,j,1)
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = rho_up !P0 / (cosh(posy/delta_rec))**2 + rho_up
!!$           p  (i,j,1) = P0 * r**2 / (cosh(posy/delta_rec))**2 + P_up !P0 / (cosh(posy/delta_rec))**2 + rho_up !P_up

           ! E. Cazzola et al. (2015)


!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0

!!$           Bx(i,j,1)  = B0 * (tanh(posy/delta_rec) + 0.5d0) - B0 * cos(posx ) &
!!$                           *  exp(-(posx**2 + posy**2 )) * (sin(posy) + 2.d0 * posy * cos(posy))
!!$           By(i,j,1)  = B0 *  cos(posy ) * exp(-(posx**2 + posy**2 )) * (sin(posx) + 2.d0 * posx * cos(posx))

!!$           Bx(i,j,1)  = B0   * (tanh(posy/delta_rec) + 0.5d0) - 2.d0 * B0 * cos(2.d0 * pi * posx /Lr) &
!!$                             *  exp(-(posx**2 + posy**2 )/delta_rec**2)                &
!!$                             * ( pi * sin(2.d0*pi*posy/Lr)/Lr + posy * cos(2.d0*pi*posy/Lr)/delta_rec**2)
!!$           By(i,j,1)  = 2.d0 * B0 * cos(2.d0 * pi * posy /Lr) &
!!$                             *  exp(-(posx**2 + posy**2 )/delta_rec**2)                &
!!$                             * ( pi * sin(2.d0*pi*posx/Lr)/Lr + posx * cos(2.d0*pi*posx/Lr)/delta_rec**2)

!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = B0 / (delta_rec * sigma * cosh(posy/delta_rec)**2)
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = rho_up * ( 1.d0 - THIRD * tanh(posy/delta_rec) - THIRD * (tanh(posy/delta_rec))**2)
!!$           p  (i,j,1) = P0     * ( 1.d0 - THIRD * tanh(posy/delta_rec) - THIRD * (tanh(posy/delta_rec))**2) + P_up !P0 / (cosh(posy/delta_rec))**2 + rho_up !P_up



           ! Miranda's attempt (2015)

!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0
!!$           Bx(i,j,1)  = 0.5d0 * B0 * posx**2 * tanh(posy/delta_rec) / (Lr**2 * cosh(posy/delta_rec)) &
!!$                              - 2.d0 * B0 * posy * exp(-( posx**2 + posy**2 )/delta_rec)/delta_rec 
!!$           By(i,j,1)  = B0 * posx    * delta_rec            / (Lr**2 * cosh(posy/delta_rec))   &
!!$                              + 2.d0 * B0 * posx * exp(-( posx**2 + posy**2 )/delta_rec)/delta_rec 
!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = 0.5d0 * B0 * ( 2.d0 * delta_rec / cosh(posy/delta_rec) &
!!$                      + posx**2 * (tanh(posy/delta_rec)-1.d0/cosh(posy/delta_rec))/(cosh(posy/delta_rec))**2)/(sigma*Lr**2)
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = rho_up
!!$           p  (i,j,1) = P0     / (cosh(posy/delta_rec))**2 + P_up 

           ! Miranda's attempt II (2015)


!!$           psi(i,j,1) = 0.d0
!!$           phi(i,j,1) = 0.d0
!!$           Bx(i,j,1)  = - THIRD * B0 * r**3 * tanh(posy/delta_rec) / cosh(posy/delta_rec) ! &
!!$!                              - 2.d0 * B0 * posy * exp(-( posx**2 + posy**2 )/delta_rec)/delta_rec 
!!$           By(i,j,1)  = B0 * posx    * delta_rec  * r / (Lr**2 * cosh(posy/delta_rec))   !&
!!$!                              + 2.d0 * B0 * posx * exp(-( posx**2 + posy**2 )/delta_rec)/delta_rec 
!!$           Bz(i,j,1)  = 0.d0
!!$           Vx(i,j,1)  = 0.d0
!!$           Vy(i,j,1)  = 0.d0
!!$           Vz(i,j,1)  = 0.d0
!!$           Ex(i,j,1)  = 0.d0
!!$           Ey(i,j,1)  = 0.d0
!!$           Ez(i,j,1)  = B0 * ( r * delta_rec - posx**2/r) / ( Lr**2 * sigma * cosh(posy/delta_rec)) &
!!$                        - (B0 * r**3) * ((tanh(posy/delta_rec))**2 - (1.d0  / cosh(posy/delta_rec))**2 )         &
!!$                        / (3.d0 * sigma)
!!$           q(i,j,1)   = 0.d0
!!$           rho(i,j,1) = rho_up
!!$           p  (i,j,1) = P0 * r**2 / (cosh(posy/delta_rec))**2 + P_up 



           !        end do !k
        end do  !j
     end do   !i

     !*****************************************************************************

  else if (DIM == 2 .and. TEST == 31) then

     ! Reconecction TEST Zenitani et al. (2010)

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax

           posx   = - 0.5d0 * Lx + i * Delx
           posy   = - 0.5d0 * Ly + j * Dely

           r      =   sqrt(1 -(posx**2/Lr**2))

           if (SLC == 1) sigma_loc(i,j,1) = sigma_0 * sigma_1 / (sigma_1 + (sigma_0 - sigma_1) &
                * (1.d0/(cosh(sqrt(posx**2+posy**2)))**2))

           if (SLC == 1)  sigma = sigma_loc(i,j,1) 

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  = B0 * tanh(posy) - 0.03d0 * B0 * posy * exp(-0.25d0 * ( posx**2 + posy**2 ))
           By(i,j,1)  = 0.03d0 * B0 * posx * exp(-0.25d0 * ( posx**2 + posy**2 ))
           Bz(i,j,1)  = 0.d0
           Vx(i,j,1)  = 0.d0
           Vy(i,j,1)  = 0.d0
           Vz(i,j,1)  = 0.d0
           Ex(i,j,1)  = 0.d0
           Ey(i,j,1)  = 0.d0
           Ez(i,j,1)  = - B0 / (sigma * cosh(posy)**2) & 
                - 0.015d0 * B0 * (posx**2 + posy**2)  * exp(-0.25d0*(posx**2+posy**2)) / sigma
           q(i,j,1)   = 0.d0
           rho(i,j,1) = P0 / (cosh(posy))**2 + rho_up
           p  (i,j,1) = P0 / (cosh(posy))**2 + rho_up !P_up

           !        end do !k
        end do  !j
     end do   !i

     !*****************************************************************************

  else if (DIM == 2 .and. TEST == 32) then

     ! TM 2D Del Zanna et al. (2016)

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx   = - 0.5d0 * Lx + i * Delx
           posy   =   j * Dely

           if (SLC == 0)  sigma = sigma_0

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0


           if ( posy .le. 0.5d0 * Ly) then

              Bx(i,j,1)  = epsilon_tm * B0 * 0.5d0 * &
                   (cos(k_tm * posy) - cos(k_tm * (posy + 0.5d0 * Ly))) / cosh(posx / a_tm)

              By(i,j,1)  = B0 * tanh(posx / a_tm) + epsilon_tm * B0 * 0.5d0    * &
                   (sin(k_tm * posy) - sin (k_tm * (posy + 0.5d0 * Ly))) * &
                   tanh(posx / a_tm) / (cosh(posx / a_tm) * kap_tm)

           else if ( posy .gt. 0.5d0 * Ly) then

              Bx(i,j,1)  = epsilon_tm * B0 * 0.5d0 * &
                   (cos(k_tm * posy) - cos(k_tm * (posy - 0.5d0 * Ly))) / cosh(posx / a_tm)

              By(i,j,1)  = B0 * tanh(posx / a_tm) + epsilon_tm * B0 * 0.5d0    * &
                   (sin(k_tm * posy) - sin (k_tm * (posy - 0.5d0 * Ly))) * &
                   tanh(posx / a_tm) / (cosh(posx / a_tm) * kap_tm)

           end if

!!$           Bx(i,j,1)  = epsilon_tm * B0 * cos(k_tm * posy) / cosh(posx / a_tm)
!!$           By(i,j,1)  = B0 * tanh(posx / a_tm) + epsilon_tm * B0 * &
!!$                        sin(k_tm * posy) * tanh(posx / a_tm) / (cosh(posx / a_tm) * kap_tm)

           Bz(i,j,1)  = B0 / cosh( posx / a_tm)

           Vx(i,j,1)  = 0.d0
           Vy(i,j,1)  = 0.d0
           Vz(i,j,1)  = 0.d0

           Ex(i,j,1)  = 0.d0
           Ey(i,j,1)  = 0.d0
           Ez(i,j,1)  = 0.d0

           q(i,j,1)   = 0.d0
           rho(i,j,1) = rho_up
           p  (i,j,1) = P0 


           !        end do !k
        end do  !j
     end do   !i

     !*****************************************************************************

  else if (DIM == 2 .and. TEST == 33) then

     ! Resistive Double Tearing Mode (RDTM)
     ! ( Hubert Baty  (2017))

     do i=0,imax
        do j=0,jmax
           !        do k=0,kmax


           posx   = - 0.5d0 * Lx + i * Delx
           posy   =   j * Dely

           if (SLC == 0)  sigma = sigma_0

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0


           if ( posx .ge. 0.d0 ) then

              Bx(i,j,1)  = epsilon_tm * B0 * cos(k_tm * posy) / cosh( (posx - l_tm) / a_tm)

              By(i,j,1)  = B0 * ( 1.d0 + tanh( (posx - l_tm) / a_tm) - tanh( (posx + l_tm) / a_tm) ) + &
                   B0 * epsilon_tm * sin(k_tm * posy) * &
                   tanh( (posx - l_tm) / a_tm) / (cosh( (posx - l_tm) / a_tm) * kap_tm)

           else if ( posx .lt. 0.d0 ) then

              Bx(i,j,1)  = epsilon_tm * B0 * cos(k_tm * posy) / cosh( (posx + l_tm) / a_tm)

              By(i,j,1)  = B0 * ( 1.d0 + tanh( (posx - l_tm) / a_tm) - tanh( (posx + l_tm) / a_tm) ) + &
                   B0 * epsilon_tm * sin(k_tm * posy) * &
                   tanh( (posx + l_tm) / a_tm) / (cosh( (posx + l_tm) / a_tm) * kap_tm)

           end if

           !           IDTM: FF Sergio setup
           !           Bz(i,j,1)  = B0 * ( 1.d0/ cosh( (posx - l_tm) / a_tm) - 1.d0/ cosh( (posx + l_tm) / a_tm) )


           !          IDTM: FF Tomek setup (a)
           Bz(i,j,1)  = B0 * sqrt(1.d0/( cosh( (l_tm - posx)/a_tm) )**2 + 1.d0/( cosh( (l_tm + posx) / a_tm) )**2      + &
                exp(- 2.d0 * l_tm / a_tm) * 1.d0/( sinh( l_tm/a_tm) ) * 1.d0/( cosh( l_tm/a_tm) )**2 * sinh(posx/a_tm) * &
                ( 1.d0/ cosh( (l_tm - posx) / a_tm) - 1.d0/ cosh( (l_tm + posx) / a_tm) ) )


           !          IDTM: FF Tomek setup (b)
           !           Bz(i,j,1)  = B0 * sqrt( 1.d0 - ( 1.d0 + tanh( (posx - l_tm) / a_tm) - tanh( (posx + l_tm) / a_tm) )**2)

           !           Bz(i,j,1)  = 0.d0

           Vx(i,j,1)  = 0.d0
           Vy(i,j,1)  = 0.d0
           Vz(i,j,1)  = 0.d0

           Ex(i,j,1)  = 0.d0
           Ey(i,j,1)  = 0.d0
           Ez(i,j,1)  = 0.d0

           q(i,j,1)   = 0.d0
           rho(i,j,1) = rho_up
           p  (i,j,1) = P_up !- 0.5d0 * &
           !( B0 * ( 1.d0 + tanh( (posx - l_tm) / a_tm) - tanh( (posx + l_tm) / a_tm) ) )**2


           !        end do !k
        end do  !j
     end do   !i

  else if (DIM == 2 .and. TEST == 33) then

     ! KHI (Miranda 2022)

     do i=0,imax
        do j=0,jmax

           posx   = - 0.5d0 * Lx + i * Delx
           posy   = - 1.0d0 * Ly + j * Dely

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  = 0.d0
           By(i,j,1)  = sqrt(0.02d0)
           Bz(i,j,1)  = sqrt(2.0d0)
           
           Vx(i,j,1)  = V0   * vsh * sin(k_kh * posy) * exp(-(posx/alpha_kh)**2)

           Vy(i,j,1)  = vsh  * tanh(posx/a_kh)

                            
           if ( Vy(i,j,1) .ge. 0.0d0 ) then
              
              rho(i,j,1) = rho_up !* tanh((posx - 0.5d0)/a_kh)
              
           else  if ( Vy(i,j,1) .lt. 0.0d0) then
              
              rho(i,j,1) =  1d-3 * rho_up !* tanh((posx + 0.5d0)/a_kh)
           
           end if
           
           Vz(i,j,1)  = 0.d0
           q(i,j,1)   = 0.d0
           p(i,j,1)   = P0 
           Ex(i,j,1)  = 0.d0! -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
           Ey(i,j,1)  = 0.d0! -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
           Ez(i,j,1)  = 0.d0! -(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))

        end do  !j
     end do   !i
  else if (DIM == 2 .and. TEST == 33) then

     ! KHI (Simulación Diamagnética)

     

     do i=0,imax
        do j=0,jmax

           posx   = - 0.5d0 * Lx + i * Delx
           posy   = - 1.0d0 * Ly + j * Dely

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  = 0.d0
           By(i,j,1)  = sqrt(0.02d0)
           Bz(i,j,1)  = sqrt(2.0d0)

           Vx(i,j,1)  = V0   * vsh * sin(k_kh * posy) * exp(-(posx/alpha_kh)**2)

           Vy(i,j,1)  = vsh  * tanh(posx/a_kh)

           if ( Vy(i,j,1) .ge. 0.0d0 ) then
              rho(i,j,1) = rho_up
           else if ( Vy(i,j,1) .lt. 0.0d0) then
              rho(i,j,1) = 1d-3 * rho_up
           end if

           Vz(i,j,1)  = 0.d0
           q(i,j,1)   = 0.d0
           p(i,j,1)   = P0
           Ex(i,j,1)  = 0.d0
           Ey(i,j,1)  = 0.d0
           Ez(i,j,1)  = 0.d0

        end do  !j
     end do   !i
  else if (DIM == 2 .and. TEST == 34) then

     ! KHI (Simulación Paramagnética)

     

     do i=0,imax
        do j=0,jmax

           posx   = - 0.5d0 * Lx + i * Delx
           posy   = - 1.0d0 * Ly + j * Dely

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  = 0.d0
           By(i,j,1)  = sqrt(0.02d0)
           Bz(i,j,1)  = sqrt(2.0d0)

           Vx(i,j,1)  = V0   * vsh * sin(k_kh * posy) * exp(-(posx/alpha_kh)**2)

           Vy(i,j,1)  = vsh  * tanh(posx/a_kh)

           if ( Vy(i,j,1) .ge. 0.0d0 ) then
              rho(i,j,1) = rho_up
           else if ( Vy(i,j,1) .lt. 0.0d0) then
              rho(i,j,1) = 1d-3 * rho_up
           end if

           Vz(i,j,1)  = 0.d0
           q(i,j,1)   = 0.d0
           p(i,j,1)   = P0
           Ex(i,j,1)  = 0.d0
           Ey(i,j,1)  = 0.d0
           Ez(i,j,1)  = 0.d0

        end do  !j
     end do   !i
  else if (DIM == 2 .and. TEST == 37) then

     ! KHI (Simulación sin Susceptibilidad Magnética)

     

     do i=0,imax
        do j=0,jmax

           posx   = - 0.5d0 * Lx + i * Delx
           posy   = - 1.0d0 * Ly + j * Dely

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  = 0.d0
           By(i,j,1)  = sqrt(0.02d0)
           Bz(i,j,1)  = sqrt(2.0d0)

           Vx(i,j,1)  = V0   * vsh * sin(k_kh * posy) * exp(-(posx/alpha_kh)**2)

           Vy(i,j,1)  = vsh  * tanh(posx/a_kh)

           if ( Vy(i,j,1) .ge. 0.0d0 ) then
              rho(i,j,1) = rho_up
           else if ( Vy(i,j,1) .lt. 0.0d0) then
              rho(i,j,1) = 1d-3 * rho_up
           end if

           Vz(i,j,1)  = 0.d0
           q(i,j,1)   = 0.d0
           p(i,j,1)   = P0
           Ex(i,j,1)  = 0.d0
           Ey(i,j,1)  = 0.d0
           Ez(i,j,1)  = 0.d0

        end do  !j
     end do   !i

  else if (DIM == 2 .and. TEST == 38) then

     ! DUAL KHI (Muzino 2013)

     do i=0,imax
        do j=0,jmax

           posx   = - 0.5d0 * Lx + i * Delx
           posy   = - 0.5d0 * Ly + j * Dely

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0
           Bx(i,j,1)  = 0.d0
           By(i,j,1)  = sqrt(0.02d0)
           Bz(i,j,1)  = sqrt(2.0d0)
           
           if ( posx .ge. 0.0d0 ) then
              
              Vx(i,j,1)  = V0   * vsh * sin(k_kh * posy) * exp(-((posx-0.5d0)/alpha_kh)**2)
              
              Vy(i,j,1)  = vsh  * tanh((posx - 0.5d0)/a_kh)

                            
           else  if ( posx .lt. 0.0d0 ) then
              
              Vx(i,j,1)  = - V0   * vsh * sin(k_kh * posy) *  exp(-((posx+0.5d0)/alpha_kh)**2)
              
              Vy(i,j,1)  = - vsh  * tanh((posx + 0.5d0)/a_kh)
           
           end if
           
           if ( Vy(i,j,1) .ge. 0.0d0 ) then
              
              rho(i,j,1) = rho_up !* tanh((posx - 0.5d0)/a_kh)
              
           else  if ( Vy(i,j,1) .lt. 0.0d0) then
              
              rho(i,j,1) =  1d-2 * rho_up !* tanh((posx + 0.5d0)/a_kh)
           
           end if
           
           Vz(i,j,1)  = 0.d0
           q(i,j,1)   = 0.d0
           p(i,j,1)   = P0 
           Ex(i,j,1)  = 0.d0! -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
           Ey(i,j,1)  = 0.d0! -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
           Ez(i,j,1)  = 0.d0!-(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))

        end do  !j
     end do   !i

  else if (DIM == 2 .and. TEST == 39) then

     ! Orszag–Tang vortex (J. Núñez-de la Rosa, and C.-D. Munz 2016)

     do i=0,imax
        do j=0,jmax

           posx   = i * Delx
           posy   = j * Dely

           psi(i,j,1) = 0.d0
           phi(i,j,1) = 0.d0

           Bx(i,j,1)  = - sin(2.d0 * pi * posy)
           By(i,j,1)  =   sin(4.d0 * pi * posx)
           Bz(i,j,1)  = 0.d0
           
           Vx(i,j,1)  = - V0 * sin(2.d0 * pi * posy)
           Vy(i,j,1)  =   V0 * sin(2.d0 * pi * posx)
           Vz(i,j,1)  = 0.d0

              
           rho(i,j,1) = gamma**2
           p(i,j,1)   = gamma
           q(i,j,1)   = 0.d0

           Ex(i,j,1)  = 0.d0! -(Vy(i,j,1)*Bz(i,j,1)-By(i,j,1)*Vz(i,j,1))
           Ey(i,j,1)  = 0.d0! -(Vz(i,j,1)*Bx(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
           Ez(i,j,1)  = 0.d0!-(Vx(i,j,1)*By(i,j,1)-Bx(i,j,1)*Vy(i,j,1))

        end do  !j
     end do   !i


     
     !*****************************************************************************



  else


     write(*,*) "STOP sunroutine int_var_primitive"
     write(*,*) "This TEST is not implemented yet"
     stop

  end if


  !_____________________________ Open files initial values ____________________



  open (unit = 101, file = "init_psi.dat")
  open (unit = 102, file = "init_phi.dat")
  open (unit = 103, file = "init_vx.dat")
  open (unit = 104, file = "init_vy.dat")
  open (unit = 105, file = "init_vz.dat")
  open (unit = 106, file = "init_bx.dat")
  open (unit = 107, file = "init_by.dat")
  open (unit = 108, file = "init_bz.dat")
  open (unit = 109, file = "init_ex.dat")
  open (unit = 110, file = "init_ey.dat")
  open (unit = 111, file = "init_ez.dat")
  open (unit = 112, file = "init_p.dat")
  open (unit = 113, file = "init_q.dat")
  open (unit = 114, file = "init_rho.dat")




  if (DIM == 1) then

     do i=0,imax

        write (101,*)  - 0.5d0 * Longx + i * Delx,  psi(i,1,1)   
        write (102,*)  - 0.5d0 * Longx + i * Delx,  phi(i,1,1)   
        write (103,*)  - 0.5d0 * Longx + i * Delx,  Vx (i,1,1)   
        write (104,*)  - 0.5d0 * Longx + i * Delx,  Vy (i,1,1)   
        write (105,*)  - 0.5d0 * Longx + i * Delx,  Vz (i,1,1)
        write (106,*)  - 0.5d0 * Longx + i * Delx,  Bx (i,1,1)
        write (107,*)  - 0.5d0 * Longx + i * Delx,  By (i,1,1)
        write (108,*)  - 0.5d0 * Longx + i * Delx,  Bz (i,1,1)
        write (109,*)  - 0.5d0 * Longx + i * Delx,  Ex (i,1,1)
        write (110,*)  - 0.5d0 * Longx + i * Delx,  Ey (i,1,1)
        write (111,*)  - 0.5d0 * Longx + i * Delx,  Ez (i,1,1)
        write (112,*)  - 0.5d0 * Longx + i * Delx,  p  (i,1,1)
        write (113,*)  - 0.5d0 * Longx + i * Delx,  q  (i,1,1)
        write (114,*)  - 0.5d0 * Longx + i * Delx,  rho(i,1,1)

     end do ! for i


  else if (DIM == 2) then

     do i=0,imax
        do j=0,jmax


           if (j == 0) then
              write (101 ,*)  ''
              write (102,*)  ''
              write (103,*)  ''
              write (104,*)  ''
              write (105,*)  ''
              write (106,*)  ''
              write (107,*)  ''
              write (108,*)  ''
              write (109,*)  ''
              write (110,*)  ''
              write (111,*) ''
              write (112,*) ''
              write (113,*) ''
              write (114,*) ''
           end if



           write (101,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  psi(i,j,1)   
           write (102,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  phi(i,j,1)   
           write (103,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Vx (i,j,1)   
           write (104,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Vy (i,j,1)   
           write (105,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Vz (i,j,1)
           write (106,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Bx (i,j,1)
           write (107,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  By (i,j,1)
           write (108,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Bz (i,j,1)
           write (109,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Ex (i,j,1)
           write (110,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Ey (i,j,1)
           write (111,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Ez (i,j,1)
           write (112,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  p  (i,j,1)
           write (113,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  q  (i,j,1)
           write (114,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  rho(i,j,1)


        end do
     end do

  end if


  close (101) 
  close (102) 
  close (103) 
  close (104) 
  close (105) 
  close (106) 
  close (107) 
  close (108) 
  close (109) 
  close (110) 
  close (111) 
  close (112) 
  close (113) 
  close (114) 






  !---------FORMAT----------
130 format (E21.16,E21.16)
131 format (E21.16,E21.16,E21.16)
123 format (E21.16,E21.16,E21.16,E21.16,E21.16,E21.16,E21.16,E21.16,E21.16,E21.16,E21.16,E21.16,E21.16,E21.16)    
  !---------FORMAT----------

end subroutine init_var_primitive
