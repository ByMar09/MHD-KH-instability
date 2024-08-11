  !     *******************************************************************
  !     Subroutine that calculates the value of intermediate auxiliary variables
  !     *******************************************************************

  subroutine  varint

    use scalar
    use parameters
    use threevectors
    use fourvectors
    implicit none

    call imex_schemes

!$OMP PARALLEL
!$OMP DO 

    do j=0+fix,jmax
       do i=0,imax

    !     Explicit intermediate values of augmented fields

    psiint(i,j,1,l) = psi(i,j,1) + Delt * psiintsum(i,j,1) 

    phiint(i,j,1,l) = phi(i,j,1) + Delt * phiintsum(i,j,1) 

    !     Explicit intermediate values of magnetic field

    Bxint(i,j,1,l)  =  Bx(i,j,1) + Delt * Bxintsum(i,j,1)

    Byint(i,j,1,l)  =  By(i,j,1) + Delt * Byintsum(i,j,1)

    Bzint(i,j,1,l)  =  Bz(i,j,1) + Delt * Bzintsum(i,j,1)


    !     Explicit intermediate value of electric field

    Exast(i,j,1,l)  =  Ex(i,j,1) + Delt * Exastsum(i,j,1)

    Eyast(i,j,1,l)  =  Ey(i,j,1) + Delt * Eyastsum(i,j,1)

    Ezast(i,j,1,l)  =  Ez(i,j,1) + Delt * Ezastsum(i,j,1)

!        print*, i, l, Ex(i,j,1), Ey(i,j,1), Ez(i,j,1), Exast(i,j,1,l), Eyast(i,j,1,l), Ezast(i,j,1,l)
  

    !     Explicit Intermediate value of charge density 

    qint(i,j,1,l) =  q(i,j,1) + Delt * qintsum(i,j,1)

    !     Mass Conservation Explicit intermediate value of D

    DDint(i,j,1,l) =  D(i,j,1) + Delt * Dintsum(i,j,1)

    !     Energy Conservation Explicit intermediate value of tau

    tauint(i,j,1,l) = tau(i,j,1) + Delt * tauintsum(i,j,1)

    !    Explicit intermediate value of S

    Sxint(i,j,1,l) = Sx(i,j,1) + Delt * Sxintsum(i,j,1)
    Syint(i,j,1,l) = Sy(i,j,1) + Delt * Syintsum(i,j,1)
    Szint(i,j,1,l) = Sz(i,j,1) + Delt * Szintsum(i,j,1)


        end do !i
     end do !j

!$OMP END DO
!$OMP END PARALLEL

     call boundary_conserved


  end subroutine varint
