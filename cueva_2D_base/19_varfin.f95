  !     *******************************************************************
  !     Subroutine that calculates the value of the final variables
  !     *******************************************************************

  subroutine varfin

    use scalar
    use parameters
    use threevectors

    implicit none

!$OMP PARALLEL
!$OMP DO 

    do j=0+fix,jmax
       do i=0,imax


! _________________________________________________________          

!  For MIRK 1st we have to keep the values of B in time "n"  when include c2 parameter         

    Bxn(i,j,1) = Bx (i,j,1)  
    Byn(i,j,1) = By (i,j,1)  
    Bzn(i,j,1) = Bz (i,j,1)

    Exn(i,j,1) = Ex (i,j,1)  
    Eyn(i,j,1) = Ey (i,j,1)  
    Ezn(i,j,1) = Ez (i,j,1)  

! _________________________________________________________
         

    psi(i,j,1) = psi(i,j,1)  + Delt * psifinsum(i,j,1)
    phi(i,j,1) = phi(i,j,1)  + Delt * phifinsum(i,j,1)

    Bx (i,j,1) = Bx (i,j,1)  + Delt * Bxfinsum (i,j,1)
    By (i,j,1) = By (i,j,1)  + Delt * Byfinsum (i,j,1)
    Bz (i,j,1) = Bz (i,j,1)  + Delt * Bzfinsum (i,j,1)

    q  (i,j,1) = q  (i,j,1)  + Delt * qfinsum  (i,j,1)

    D  (i,j,1) = D  (i,j,1)  + Delt * Dfinsum  (i,j,1)

    tau(i,j,1) = tau(i,j,1)  + Delt * taufinsum(i,j,1)

    Sx (i,j,1) = Sx (i,j,1)  + Delt * Sxfinsum (i,j,1)
    Sy (i,j,1) = Sy (i,j,1)  + Delt * Syfinsum (i,j,1)
    Sz (i,j,1) = Sz (i,j,1)  + Delt * Szfinsum (i,j,1)


      end do ! for i
    end do ! for j

!$OMP END DO
!$OMP END PARALLEL

  end subroutine varfin
