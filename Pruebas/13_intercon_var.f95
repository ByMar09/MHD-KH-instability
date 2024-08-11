
  !     *******************************************************************
  !     Subroutine which Starts Intermediate Conserved Variables 
  !     *******************************************************************

  subroutine interconvar

    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

!$OMP PARALLEL PRIVATE(V2int)
!$OMP DO 

    do j=0+fix,jmax
       do i=0,imax

              !     Velocity module

              V2int = Vxint(i,j,1,l)**2 +Vyint(i,j,1,l)**2 +Vzint(i,j,1,l)**2 

              !     Lorentz factor

              Wint(i,j,1,l)  = 1.d0/sqrt(1-V2int)


              !     Electric field module

              E2int(i,j,1,l) = Exint(i,j,1,l)**2 + Eyint(i,j,1,l)**2 + Ezint(i,j,1,l)**2

              !     Magnetic field module 

              B2int(i,j,1,l) = Bxint(i,j,1,l)**2 + Byint(i,j,1,l)**2 +Bzint(i,j,1,l)**2

              !     Specific internal energy ! Toro E. F (1.48)

              epsilonint(i,j,1,l) = pint(i,j,1,l)/((gamma-1.d0)*rhoint(i,j,1,l))

              !     Entalphy 

              enthpyint(i,j,1,l) = rhoint(i,j,1,l) * (1.d0+epsilonint(i,j,1,l)) + pint(i,j,1,l) 

              ! Valores iniciales de los campos conservados derivados de los valores iniciales de las variables primitivas 

              !     Conserved mass density (D = \rho W) 

              DDint(i,j,1,l) = rhoint(i,j,1,l) * Wint(i,j,1,l)

              !     Conserved Energy density \tau = 1/2 (E^2+B^2)+ h W^2- P

              tauint(i,j,1,l) = 0.5d0 * (E2int(i,j,1,l) + B2int(i,j,1,l) ) + enthpyint(i,j,1,l) * Wint(i,j,1,l)**2 - pint(i,j,1,l)

              !     Conserved Momentum S = E X B + h W^2 V 

              Sxint(i,j,1,l) = (Bzint(i,j,1,l)*Eyint(i,j,1,l) - Byint(i,j,1,l)*Ezint(i,j,1,l)) + &
                                enthpyint(i,j,1,l)*Wint(i,j,1,l)**2 * Vxint(i,j,1,l)

              Syint(i,j,1,l) = (Bxint(i,j,1,l)*Ezint(i,j,1,l) - Bzint(i,j,1,l)*Exint(i,j,1,l)) + &
                                enthpyint(i,j,1,l)*Wint(i,j,1,l)**2 * Vyint(i,j,1,l)

              Szint(i,j,1,l) = (Byint(i,j,1,l)*Exint(i,j,1,l) - Bxint(i,j,1,l)*Eyint(i,j,1,l)) + &
                                enthpyint(i,j,1,l)*Wint(i,j,1,l)**2 * Vzint(i,j,1,l)


        end do !i
     end do !j

!$OMP END DO
!$OMP END PARALLEL

     call boundary_conserved
     call boundary_electric

   end subroutine interconvar
