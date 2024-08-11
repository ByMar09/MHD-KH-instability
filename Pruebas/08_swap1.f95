
  subroutine  swap1

    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

!$OMP PARALLEL
!$OMP DO 

    do j=0+fix,jmax
       do i=0,imax



    rhoint(i,j,1,1)     = rho(i,j,1) 

    pint  (i,j,1,1)     = p  (i,j,1) 

    Vxint (i,j,1,1)     = Vx (i,j,1) 
    Vyint (i,j,1,1)     = Vy (i,j,1) 
    Vzint (i,j,1,1)     = Vz (i,j,1)

    Wint  (i,j,1,1)     = W  (i,j,1)

    enthpyint(i,j,1,1)  = enthpy(i,j,1) 

    epsilonint(i,j,1,1) = epsiln(i,j,1) 

        end do ! for i
     end do ! for j

!$OMP END DO
!$OMP END PARALLEL

  end subroutine swap1
