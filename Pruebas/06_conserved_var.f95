  !     *******************************************************************
  !     Subroutine which Starts Conserved Variables 
  !     *******************************************************************

  subroutine startconvar

    use scalar
    use parameters
    use threevectors

    implicit none

!$OMP PARALLEL
    
    if (DIM == 1) then

!$OMP DO 

     do i=0,imax

              !     Velocity module

              V2(i,1,1)= Vx(i,1,1)**2 +Vy(i,1,1)**2 +Vz(i,1,1)**2 

              !     Lorentz factor

              W(i,1,1) = 1.d0/sqrt(1.d0-V2(i,1,1))


              !     Electric field module

              E2(i,1,1)=Ex(i,1,1)**2 + Ey(i,1,1)**2 + Ez(i,1,1)**2

              !     Magnetic field module 

              B2(i,1,1)=Bx(i,1,1)**2 + By(i,1,1)**2 + Bz(i,1,1)**2

              !     Specific internal energy ! Toro E. F (1.48)

              epsiln(i,1,1) = p(i,1,1)/((gamma-1.d0)*rho(i,1,1))

              !     Entalphy 

              enthpy(i,1,1) = rho(i,1,1) * (1.d0+epsiln(i,1,1)) + p(i,1,1) 

              ! Valores iniciales de los campos conservados derivados de los valores iniciales de las variables primitivas 

              !     Conserved mass density (D = \rho W) 

              D(i,1,1) = rho(i,1,1) * W(i,1,1)

              !     Conserved Energy density \tau = 1/2 (E^2+B^2)+ h W^2- P

              tau(i,1,1) = 0.5 * (E2(i,1,1) + B2(i,1,1)) &
                               +  enthpy(i,1,1) * W(i,1,1)**2 - p(i,1,1)

              !     Conserved Momentum S = E X B + h W^2 V 

              Sx(i,1,1) = (Bz(i,1,1)*Ey(i,1,1) -       &
                   By(i,1,1)*Ez(i,1,1))+enthpy(i,1,1)* &
                   W(i,1,1)**2 * Vx(i,1,1)
              Sy(i,1,1) = (Bx(i,1,1)*Ez(i,1,1) -       &
                   Bz(i,1,1)*Ex(i,1,1))+enthpy(i,1,1)* &
                   W(i,1,1)**2 * Vy(i,1,1)

              Sz(i,1,1) = (By(i,1,1)*Ex(i,1,1) -       &
                   Bx(i,1,1)*Ey(i,1,1))+enthpy(i,1,1)* &
                   W(i,1,1)**2 * Vz(i,1,1)

              
              if (TEST == 3 .and. SLC == 1) sigma_loc(i,1,1) = sigma_0 * D(i,1,1)**gamma_slc

     end do

!$OMP END DO

    else if (DIM == 2) then

!$OMP DO 

     do j=0,jmax
        do i=0,imax
!           do k=1,kmax


              !     Velocity module

              V2(i,j,1)= Vx(i,j,1)**2 +Vy(i,j,1)**2 +Vz(i,j,1)**2 

              !     Lorentz factor

              W(i,j,1) = 1.d0/sqrt(1.d0-V2(i,j,1))


              !     Electric field module

              E2(i,j,1)=Ex(i,j,1)**2 + Ey(i,j,1)**2 + Ez(i,j,1)**2

              !     Magnetic field module 

              B2(i,j,1)=Bx(i,j,1)**2 + By(i,j,1)**2 + Bz(i,j,1)**2

              !     Specific internal energy ! Toro E. F (1.48)

              epsiln(i,j,1) = p(i,j,1)/((gamma-1.d0)*rho(i,j,1))

              !     Entalphy 

              enthpy(i,j,1) = rho(i,j,1) * (1.d0+epsiln(i,j,1)) + p(i,j,1) 

              ! Valores iniciales de los campos conservados derivados de los valores iniciales de las variables primitivas 

              !     Conserved mass density (D = \rho W) 

              D(i,j,1) = rho(i,j,1) * W(i,j,1)

              !     Conserved Energy density \tau = 1/2 (E^2+B^2)+ h W^2- P

              tau(i,j,1) = 0.5d0 * (E2(i,j,1) + B2(i,j,1)) &
                                 +  enthpy(i,j,1) * W(i,j,1)**2 - p(i,j,1)

              !     Conserved Momentum S = E X B + h W^2 V 

              Sx(i,j,1) = (Bz(i,j,1)*Ey(i,j,1) -       &
                   By(i,j,1)*Ez(i,j,1))+enthpy(i,j,1)* &
                   W(i,j,1)**2 * Vx(i,j,1)
              Sy(i,j,1) = (Bx(i,j,1)*Ez(i,j,1) -       &
                   Bz(i,j,1)*Ex(i,j,1))+enthpy(i,j,1)* &
                   W(i,j,1)**2 * Vy(i,j,1)

              Sz(i,j,1) = (By(i,j,1)*Ex(i,j,1) -       &
                   Bx(i,j,1)*Ey(i,j,1))+enthpy(i,j,1)* &
                   W(i,j,1)**2 * Vz(i,j,1)


!           end do
        end do
     end do

!$OMP END DO

   else

       write(*,*)  "STOP: subroutine startconvar "
       write(*,*)  " this dimension is not implemented yet"
       stop

    end if

!$OMP END PARALLEL
    
  end subroutine startconvar
