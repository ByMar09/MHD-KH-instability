
 subroutine  swap2

    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

!$OMP PARALLEL
!$OMP DO 

    do j=0+fix,jmax
       do i=0,imax


    rho(i,j,1)     =  rhoint(i,j,1,lmax)  

    p  (i,j,1)     =  pint  (i,j,1,lmax)

    Vx (i,j,1)     =  Vxint (i,j,1,lmax)
    Vy (i,j,1)     =  Vyint (i,j,1,lmax)
    Vz (i,j,1)     =  Vzint (i,j,1,lmax)

    W  (i,j,1)     =  Wint  (i,j,1,lmax)

    enthpy(i,j,1)  =  enthpyint(i,j,1,lmax)

    epsiln(i,j,1) =  epsilonint(i,j,1,lmax)

    
        end do ! for i
     end do ! for j

!$OMP END DO
!$OMP END PARALLEL

  end subroutine swap2
