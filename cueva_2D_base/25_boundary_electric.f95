  subroutine boundary_electric

    use scalar
    use parameters
    use threevectors 
    use fourvectors

    implicit none

!$OMP PARALLEL PRIVATE(posy)

    if (DIM == 1 .and. BOUND == 1 .or. BOUND == 3 .or. BOUND == 4) then

       ! Copy boundary 1D

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

!!$            Exint(0,1,1,l)  = Exint(1,1,1,l)  
!!$            Eyint(0,1,1,l)  = Eyint(1,1,1,l)  
!!$            Ezint(0,1,1,l)  = Ezint(1,1,1,l)  

            Exint(-1,1,1,l)  = Exint(0,1,1,l)  
            Eyint(-1,1,1,l)  = Eyint(0,1,1,l)  
            Ezint(-1,1,1,l)  = Ezint(0,1,1,l)  
 
            Exint(-2,1,1,l)  = Exint(-1,1,1,l)  
            Eyint(-2,1,1,l)  = Eyint(-1,1,1,l)  
            Ezint(-2,1,1,l)  = Ezint(-1,1,1,l)  

            Exint(-3,1,1,l)  = Exint(-2,1,1,l)  
            Eyint(-3,1,1,l)  = Eyint(-2,1,1,l)  
            Ezint(-3,1,1,l)  = Ezint(-2,1,1,l)  

            Exint(-4,1,1,l)  = Exint(-3,1,1,l)  
            Eyint(-4,1,1,l)  = Eyint(-3,1,1,l)  
            Ezint(-4,1,1,l)  = Ezint(-3,1,1,l)

            Exint(-5,1,1,l)  = Exint(-4,1,1,l)  
            Eyint(-5,1,1,l)  = Eyint(-4,1,1,l)  
            Ezint(-5,1,1,l)  = Ezint(-4,1,1,l)

            Exint(-6,1,1,l)  = Exint(-5,1,1,l)  
            Eyint(-6,1,1,l)  = Eyint(-5,1,1,l)  
            Ezint(-6,1,1,l)  = Ezint(-5,1,1,l)



            ! Right Boundary
            !---------------

!!$            Exint(imax  ,1,1,l)  = Exint(imax-1,1,1,l)  
!!$            Eyint(imax  ,1,1,l)  = Eyint(imax-1,1,1,l)  
!!$            Ezint(imax  ,1,1,l)  = Ezint(imax-1,1,1,l)  

            Exint(imax+1,1,1,l)  = Exint(imax  ,1,1,l)  
            Eyint(imax+1,1,1,l)  = Eyint(imax  ,1,1,l)  
            Ezint(imax+1,1,1,l)  = Ezint(imax  ,1,1,l)  
  
            Exint(imax+2,1,1,l)  = Exint(imax+1,1,1,l)  
            Eyint(imax+2,1,1,l)  = Eyint(imax+1,1,1,l)  
            Ezint(imax+2,1,1,l)  = Ezint(imax+1,1,1,l)

            Exint(imax+3,1,1,l)  = Exint(imax+2,1,1,l)  
            Eyint(imax+3,1,1,l)  = Eyint(imax+2,1,1,l)  
            Ezint(imax+3,1,1,l)  = Ezint(imax+2,1,1,l)

            Exint(imax+4,1,1,l)  = Exint(imax+3,1,1,l)  
            Eyint(imax+4,1,1,l)  = Eyint(imax+3,1,1,l)  
            Ezint(imax+4,1,1,l)  = Ezint(imax+3,1,1,l)

            Exint(imax+5,1,1,l)  = Exint(imax+4,1,1,l)  
            Eyint(imax+5,1,1,l)  = Eyint(imax+4,1,1,l)  
            Ezint(imax+5,1,1,l)  = Ezint(imax+4,1,1,l)

            Exint(imax+6,1,1,l)  = Exint(imax+5,1,1,l)  
            Eyint(imax+6,1,1,l)  = Eyint(imax+5,1,1,l)  
            Ezint(imax+6,1,1,l)  = Ezint(imax+5,1,1,l)


    else if (DIM == 1 .and. BOUND == 2) then

       ! Periodic boundary 1D

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

!!$            Exint( 0,1,1,l)  = Exint(imax-1,1,1,l)  
!!$            Eyint( 0,1,1,l)  = Eyint(imax-1,1,1,l)  
!!$            Ezint( 0,1,1,l)  = Ezint(imax-1,1,1,l)  

            Exint(-1,1,1,l)  = Exint(imax-1,1,1,l)  
            Eyint(-1,1,1,l)  = Eyint(imax-1,1,1,l)  
            Ezint(-1,1,1,l)  = Ezint(imax-1,1,1,l)  
 
            Exint(-2,1,1,l)  = Exint(imax-2,1,1,l)  
            Eyint(-2,1,1,l)  = Eyint(imax-2,1,1,l)  
            Ezint(-2,1,1,l)  = Ezint(imax-2,1,1,l)  

            Exint(-3,1,1,l)  = Exint(imax-3,1,1,l)  
            Eyint(-3,1,1,l)  = Eyint(imax-3,1,1,l)  
            Ezint(-3,1,1,l)  = Ezint(imax-3,1,1,l)  

            Exint(-4,1,1,l)  = Exint(imax-4,1,1,l)  
            Eyint(-4,1,1,l)  = Eyint(imax-4,1,1,l)  
            Ezint(-4,1,1,l)  = Ezint(imax-4,1,1,l)

            Exint(-5,1,1,l)  = Exint(imax-5,1,1,l)  
            Eyint(-5,1,1,l)  = Eyint(imax-5,1,1,l)  
            Ezint(-5,1,1,l)  = Ezint(imax-5,1,1,l)

            Exint(-6,1,1,l)  = Exint(imax-6,1,1,l)  
            Eyint(-6,1,1,l)  = Eyint(imax-6,1,1,l)  
            Ezint(-6,1,1,l)  = Ezint(imax-6,1,1,l)  


            ! Right Boundary
            !---------------

!!$            Exint(imax  ,1,1,l)  = Exint(1,1,1,l)  
!!$            Eyint(imax  ,1,1,l)  = Eyint(1,1,1,l)  
!!$            Ezint(imax  ,1,1,l)  = Ezint(1,1,1,l)  

            Exint(imax+1,1,1,l)  = Exint(1,1,1,l)  
            Eyint(imax+1,1,1,l)  = Eyint(1,1,1,l)  
            Ezint(imax+1,1,1,l)  = Ezint(1,1,1,l)  
  
            Exint(imax+2,1,1,l)  = Exint(2,1,1,l)  
            Eyint(imax+2,1,1,l)  = Eyint(2,1,1,l)  
            Ezint(imax+2,1,1,l)  = Ezint(2,1,1,l) 

            Exint(imax+3,1,1,l)  = Exint(3,1,1,l)  
            Eyint(imax+3,1,1,l)  = Eyint(3,1,1,l)  
            Ezint(imax+3,1,1,l)  = Ezint(3,1,1,l) 

            Exint(imax+4,1,1,l)  = Exint(4,1,1,l)  
            Eyint(imax+4,1,1,l)  = Eyint(4,1,1,l)  
            Ezint(imax+4,1,1,l)  = Ezint(4,1,1,l)

            Exint(imax+5,1,1,l)  = Exint(5,1,1,l)  
            Eyint(imax+5,1,1,l)  = Eyint(5,1,1,l)  
            Ezint(imax+5,1,1,l)  = Ezint(5,1,1,l)

            Exint(imax+6,1,1,l)  = Exint(6,1,1,l)  
            Eyint(imax+6,1,1,l)  = Eyint(6,1,1,l)  
            Ezint(imax+6,1,1,l)  = Ezint(6,1,1,l) 


     else if (DIM == 1 .and. BOUND == 8) then

       ! Perfect Conductor and anti-simetric reflection TM_III

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

!!$            Exint(0,1,1,l)  = Exint(     1,1,1,l)  
!!$            Eyint(0,1,1,l)  = Eyint(imax-1,1,1,l)  
!!$            Ezint(0,1,1,l)  = Ezint(imax-1,1,1,l)  

            Exint(-1,1,1,l)  = 2.d0 * Exint(0,1,1,l) - Exint(1,1,1,l)
            Eyint(-1,1,1,l)  = 2.d0 * Eyint(0,1,1,l) - Eyint(1,1,1,l)
            Ezint(-1,1,1,l)  = 2.d0 * Ezint(0,1,1,l) - Ezint(1,1,1,l)
 
            Exint(-2,1,1,l)  = 2.d0 * Exint(0,1,1,l) - Exint(2,1,1,l)
            Eyint(-2,1,1,l)  = 2.d0 * Eyint(0,1,1,l) - Eyint(2,1,1,l)
            Ezint(-2,1,1,l)  = 2.d0 * Ezint(0,1,1,l) - Ezint(2,1,1,l)



            ! Right Boundary
            !---------------

!!$            Exint(imax,1,1,l)  = Exint(imax-1,1,1,l)  
!!$            Eyint(imax,1,1,l)  = Eyint(1,1,1,l)  
!!$            Ezint(imax,1,1,l)  = Ezint(1,1,1,l)  

            Exint(imax+1,1,1,l)  = 2.d0 * Exint(imax,1,1,l) - Exint(imax-1,1,1,l)
            Eyint(imax+1,1,1,l)  = 2.d0 * Eyint(imax,1,1,l) - Eyint(imax-1,1,1,l)
            Ezint(imax+1,1,1,l)  = 2.d0 * Ezint(imax,1,1,l) - Ezint(imax-1,1,1,l)
  
            Exint(imax+2,1,1,l)  = 2.d0 * Exint(imax,1,1,l) - Exint(imax-2,1,1,l)
            Eyint(imax+2,1,1,l)  = 2.d0 * Eyint(imax,1,1,l) - Eyint(imax-2,1,1,l)
            Ezint(imax+2,1,1,l)  = 2.d0 * Ezint(imax,1,1,l) - Ezint(imax-2,1,1,l)


    else if (DIM == 2 .and. BOUND == 1 .or. BOUND == 3) then

       ! Copy boundary 2D

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO
       
      do j=-6,jmax+6

            !-----------------
            ! Left Boundary
            !---------------

!!$            Exint(0,j,1,l)  = Exint( 1,j,1,l) 
!!$            Eyint(0,j,1,l)  = Eyint( 1,j,1,l) 
!!$            Ezint(0,j,1,l)  = Ezint( 1,j,1,l) 
 
            Exint(-1,j,1,l)  = Exint( 0,j,1,l) 
            Eyint(-1,j,1,l)  = Eyint( 0,j,1,l) 
            Ezint(-1,j,1,l)  = Ezint( 0,j,1,l) 

            Exint(-2,j,1,l)  = Exint(-1,j,1,l) 
            Eyint(-2,j,1,l)  = Eyint(-1,j,1,l) 
            Ezint(-2,j,1,l)  = Ezint(-1,j,1,l) 

            Exint(-3,j,1,l)  = Exint(-2,j,1,l) 
            Eyint(-3,j,1,l)  = Eyint(-2,j,1,l) 
            Ezint(-3,j,1,l)  = Ezint(-2,j,1,l) 

            Exint(-4,j,1,l)  = Exint(-3,j,1,l) 
            Eyint(-4,j,1,l)  = Eyint(-3,j,1,l) 
            Ezint(-4,j,1,l)  = Ezint(-3,j,1,l)

            Exint(-5,j,1,l)  = Exint(-4,j,1,l) 
            Eyint(-5,j,1,l)  = Eyint(-4,j,1,l) 
            Ezint(-5,j,1,l)  = Ezint(-4,j,1,l)

            Exint(-6,j,1,l)  = Exint(-5,j,1,l) 
            Eyint(-6,j,1,l)  = Eyint(-5,j,1,l) 
            Ezint(-6,j,1,l)  = Ezint(-5,j,1,l) 

            !-----------------
            ! Right Boundary
            !---------------

!!$            Exint(imax,j,1,l)  = Exint(imax-1,j,1,l)  
!!$            Eyint(imax,j,1,l)  = Eyint(imax-1,j,1,l)  
!!$            Ezint(imax,j,1,l)  = Ezint(imax-1,j,1,l)  

            Exint(imax+1,j,1,l)  = Exint(imax,j,1,l)  
            Eyint(imax+1,j,1,l)  = Eyint(imax,j,1,l)  
            Ezint(imax+1,j,1,l)  = Ezint(imax,j,1,l)  

            Exint(imax+2,j,1,l)  = Exint(imax+1,j,1,l)  
            Eyint(imax+2,j,1,l)  = Eyint(imax+1,j,1,l)  
            Ezint(imax+2,j,1,l)  = Ezint(imax+1,j,1,l)  

            Exint(imax+3,j,1,l)  = Exint(imax+2,j,1,l)  
            Eyint(imax+3,j,1,l)  = Eyint(imax+2,j,1,l)  
            Ezint(imax+3,j,1,l)  = Ezint(imax+2,j,1,l)  

            Exint(imax+4,j,1,l)  = Exint(imax+3,j,1,l)  
            Eyint(imax+4,j,1,l)  = Eyint(imax+3,j,1,l)  
            Ezint(imax+4,j,1,l)  = Ezint(imax+3,j,1,l)

            Exint(imax+5,j,1,l)  = Exint(imax+4,j,1,l)  
            Eyint(imax+5,j,1,l)  = Eyint(imax+4,j,1,l)  
            Ezint(imax+5,j,1,l)  = Ezint(imax+4,j,1,l)

            Exint(imax+6,j,1,l)  = Exint(imax+5,j,1,l)  
            Eyint(imax+6,j,1,l)  = Eyint(imax+5,j,1,l)  
            Ezint(imax+6,j,1,l)  = Ezint(imax+5,j,1,l)  

         end do

!$OMP END DO         

  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO
         
         do i=-6,imax+6

            ! Left Boundary
            !---------------

!!$            Exint (i, 0,1,l) = Exint(i, 1,1,l)  
!!$            Eyint (i, 0,1,l) = Eyint(i, 1,1,l)  
!!$            Ezint (i, 0,1,l) = Ezint(i, 1,1,l)  

            Exint (i,-1,1,l) = Exint(i, 0,1,l)  
            Eyint (i,-1,1,l) = Eyint(i, 0,1,l)  
            Ezint (i,-1,1,l) = Ezint(i, 0,1,l)  

            Exint (i,-2,1,l) = Exint(i,-1,1,l)  
            Eyint (i,-2,1,l) = Eyint(i,-1,1,l)  
            Ezint (i,-2,1,l) = Ezint(i,-1,1,l)  

            Exint (i,-3,1,l) = Exint(i,-2,1,l)  
            Eyint (i,-3,1,l) = Eyint(i,-2,1,l)  
            Ezint (i,-3,1,l) = Ezint(i,-2,1,l)  

            Exint (i,-4,1,l) = Exint(i,-3,1,l)  
            Eyint (i,-4,1,l) = Eyint(i,-3,1,l)  
            Ezint (i,-4,1,l) = Ezint(i,-3,1,l)

            Exint (i,-5,1,l) = Exint(i,-4,1,l)  
            Eyint (i,-5,1,l) = Eyint(i,-4,1,l)  
            Ezint (i,-5,1,l) = Ezint(i,-4,1,l)

            Exint (i,-6,1,l) = Exint(i,-5,1,l)  
            Eyint (i,-6,1,l) = Eyint(i,-5,1,l)  
            Ezint (i,-6,1,l) = Ezint(i,-5,1,l)  
            

            ! Right Boundary
            !---------------

!!$            Exint (i,jmax,1,l) = Exint(i,jmax-1,1,l)  
!!$            Eyint (i,jmax,1,l) = Eyint(i,jmax-1,1,l)  
!!$            Ezint (i,jmax,1,l) = Ezint(i,jmax-1,1,l)  

            Exint (i,jmax+1,1,l) = Exint(i,jmax,1,l)  
            Eyint (i,jmax+1,1,l) = Eyint(i,jmax,1,l)  
            Ezint (i,jmax+1,1,l) = Ezint(i,jmax,1,l)  

            Exint (i,jmax+2,1,l) = Exint(i,jmax+1,1,l)  
            Eyint (i,jmax+2,1,l) = Eyint(i,jmax+1,1,l)  
            Ezint (i,jmax+2,1,l) = Ezint(i,jmax+1,1,l)  

            Exint (i,jmax+3,1,l) = Exint(i,jmax+2,1,l)  
            Eyint (i,jmax+3,1,l) = Eyint(i,jmax+2,1,l)  
            Ezint (i,jmax+3,1,l) = Ezint(i,jmax+2,1,l)  

            Exint (i,jmax+4,1,l) = Exint(i,jmax+3,1,l)  
            Eyint (i,jmax+4,1,l) = Eyint(i,jmax+3,1,l)  
            Ezint (i,jmax+4,1,l) = Ezint(i,jmax+3,1,l)

            Exint (i,jmax+5,1,l) = Exint(i,jmax+4,1,l)  
            Eyint (i,jmax+5,1,l) = Eyint(i,jmax+4,1,l)  
            Ezint (i,jmax+5,1,l) = Ezint(i,jmax+4,1,l)

            Exint (i,jmax+6,1,l) = Exint(i,jmax+5,1,l)  
            Eyint (i,jmax+6,1,l) = Eyint(i,jmax+5,1,l)  
            Ezint (i,jmax+6,1,l) = Ezint(i,jmax+5,1,l)  

         end do

!$OMP END DO
         

    else if (DIM == 2 .and. BOUND == 2) then

       ! Periodic boundary 2D


  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
      do j=-6,jmax+6

            !-----------------
            ! Left Boundary
            !---------------

!!$            Exint( 0,j,1,l)  = Exint(imax-1,j,1,l) 
!!$            Eyint( 0,j,1,l)  = Eyint(imax-1,j,1,l) 
!!$            Ezint( 0,j,1,l)  = Ezint(imax-1,j,1,l) 
 
            Exint(-1,j,1,l)  = Exint(imax-1,j,1,l) 
            Eyint(-1,j,1,l)  = Eyint(imax-1,j,1,l) 
            Ezint(-1,j,1,l)  = Ezint(imax-1,j,1,l) 

            Exint(-2,j,1,l)  = Exint(imax-2,j,1,l) 
            Eyint(-2,j,1,l)  = Eyint(imax-2,j,1,l) 
            Ezint(-2,j,1,l)  = Ezint(imax-2,j,1,l) 

            Exint(-3,j,1,l)  = Exint(imax-3,j,1,l) 
            Eyint(-3,j,1,l)  = Eyint(imax-3,j,1,l) 
            Ezint(-3,j,1,l)  = Ezint(imax-3,j,1,l) 

            Exint(-4,j,1,l)  = Exint(imax-4,j,1,l) 
            Eyint(-4,j,1,l)  = Eyint(imax-4,j,1,l) 
            Ezint(-4,j,1,l)  = Ezint(imax-4,j,1,l)

            Exint(-5,j,1,l)  = Exint(imax-5,j,1,l) 
            Eyint(-5,j,1,l)  = Eyint(imax-5,j,1,l) 
            Ezint(-5,j,1,l)  = Ezint(imax-5,j,1,l)
            
            Exint(-6,j,1,l)  = Exint(imax-6,j,1,l) 
            Eyint(-6,j,1,l)  = Eyint(imax-6,j,1,l) 
            Ezint(-6,j,1,l)  = Ezint(imax-6,j,1,l)
            


            !-----------------
            ! Right Boundary
            !---------------

!!$            Exint(imax  ,j,1,l)  = Exint(1,j,1,l)  
!!$            Eyint(imax  ,j,1,l)  = Eyint(1,j,1,l)  
!!$            Ezint(imax  ,j,1,l)  = Ezint(1,j,1,l) 

            Exint(imax+1,j,1,l)  = Exint(1,j,1,l)  
            Eyint(imax+1,j,1,l)  = Eyint(1,j,1,l)  
            Ezint(imax+1,j,1,l)  = Ezint(1,j,1,l)  

            Exint(imax+2,j,1,l)  = Exint(2,j,1,l)  
            Eyint(imax+2,j,1,l)  = Eyint(2,j,1,l)  
            Ezint(imax+2,j,1,l)  = Ezint(2,j,1,l)  

            Exint(imax+3,j,1,l)  = Exint(3,j,1,l)  
            Eyint(imax+3,j,1,l)  = Eyint(3,j,1,l)  
            Ezint(imax+3,j,1,l)  = Ezint(3,j,1,l)  

            Exint(imax+4,j,1,l)  = Exint(4,j,1,l)  
            Eyint(imax+4,j,1,l)  = Eyint(4,j,1,l)  
            Ezint(imax+4,j,1,l)  = Ezint(4,j,1,l)

            Exint(imax+5,j,1,l)  = Exint(5,j,1,l)  
            Eyint(imax+5,j,1,l)  = Eyint(5,j,1,l)  
            Ezint(imax+5,j,1,l)  = Ezint(5,j,1,l)

            Exint(imax+6,j,1,l)  = Exint(6,j,1,l)  
            Eyint(imax+6,j,1,l)  = Eyint(6,j,1,l)  
            Ezint(imax+6,j,1,l)  = Ezint(6,j,1,l)  

         end do

!$OMP END DO
         
  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
         
         do i=-6,imax+6

            !-----------------
            ! Left Boundary
            !---------------

!!$            Exint (i, 0,1,l) = Exint(i,jmax-1,1,l)  
!!$            Eyint (i, 0,1,l) = Eyint(i,jmax-1,1,l)  
!!$            Ezint (i, 0,1,l) = Ezint(i,jmax-1,1,l)  

            Exint (i,-1,1,l) = Exint(i,jmax-1,1,l)  
            Eyint (i,-1,1,l) = Eyint(i,jmax-1,1,l)  
            Ezint (i,-1,1,l) = Ezint(i,jmax-1,1,l)  

            Exint (i,-2,1,l) = Exint(i,jmax-2,1,l)  
            Eyint (i,-2,1,l) = Eyint(i,jmax-2,1,l)  
            Ezint (i,-2,1,l) = Ezint(i,jmax-2,1,l)  

            Exint (i,-3,1,l) = Exint(i,jmax-3,1,l)  
            Eyint (i,-3,1,l) = Eyint(i,jmax-3,1,l)  
            Ezint (i,-3,1,l) = Ezint(i,jmax-3,1,l)  

            Exint (i,-4,1,l) = Exint(i,jmax-4,1,l)  
            Eyint (i,-4,1,l) = Eyint(i,jmax-4,1,l)  
            Ezint (i,-4,1,l) = Ezint(i,jmax-4,1,l)

            Exint (i,-5,1,l) = Exint(i,jmax-5,1,l)  
            Eyint (i,-5,1,l) = Eyint(i,jmax-5,1,l)  
            Ezint (i,-5,1,l) = Ezint(i,jmax-5,1,l)

            Exint (i,-6,1,l) = Exint(i,jmax-6,1,l)  
            Eyint (i,-6,1,l) = Eyint(i,jmax-6,1,l)  
            Ezint (i,-6,1,l) = Ezint(i,jmax-6,1,l)  
            
            !-----------------
            ! Right Boundary
            !---------------

!!$            Exint (i,jmax  ,1,l) = Exint(i,1,1,l)  
!!$            Eyint (i,jmax  ,1,l) = Eyint(i,1,1,l)  
!!$            Ezint (i,jmax  ,1,l) = Ezint(i,1,1,l)  

            Exint (i,jmax+1,1,l) = Exint(i,1,1,l)  
            Eyint (i,jmax+1,1,l) = Eyint(i,1,1,l)  
            Ezint (i,jmax+1,1,l) = Ezint(i,1,1,l)  

            Exint (i,jmax+2,1,l) = Exint(i,2,1,l)  
            Eyint (i,jmax+2,1,l) = Eyint(i,2,1,l)  
            Ezint (i,jmax+2,1,l) = Ezint(i,2,1,l)  

            Exint (i,jmax+3,1,l) = Exint(i,3,1,l)  
            Eyint (i,jmax+3,1,l) = Eyint(i,3,1,l)  
            Ezint (i,jmax+3,1,l) = Ezint(i,3,1,l)  

            Exint (i,jmax+4,1,l) = Exint(i,4,1,l)  
            Eyint (i,jmax+4,1,l) = Eyint(i,4,1,l)  
            Ezint (i,jmax+4,1,l) = Ezint(i,4,1,l)

            Exint (i,jmax+5,1,l) = Exint(i,5,1,l)  
            Eyint (i,jmax+5,1,l) = Eyint(i,5,1,l)  
            Ezint (i,jmax+5,1,l) = Ezint(i,5,1,l)

            Exint (i,jmax+6,1,l) = Exint(i,6,1,l)  
            Eyint (i,jmax+6,1,l) = Eyint(i,6,1,l)  
            Ezint (i,jmax+6,1,l) = Ezint(i,6,1,l)

         end do

!$OMP END DO

         
    else if (DIM == 2 .and. BOUND == 4) then

       ! Slab Jet Boundary

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
      do j=-6,jmax+6

           posy  = - 14.d0 + j * Dely

           !-----------------           
           ! Left Boundary
           ! -----------------

         if (posy .gt. -1.d0 .and. posy .le. 1.d0) then

              if (REC_PRIM == 0) then

            Exint(0,j,1,l)  = -(Vy(0,j,1)*Bzint(0,j,1,l)-Byint(0,j,1,l)*Vz(0,j,1)) 
            Eyint(0,j,1,l)  = -(Vz(0,j,1)*Bxint(0,j,1,l)-Bzint(0,j,1,l)*Vx(0,j,1)) 
            Ezint(0,j,1,l)  = -(Vx(0,j,1)*Byint(0,j,1,l)-Bxint(0,j,1,l)*Vy(0,j,1)) 

            Exint(-1,j,1,l)  = -(Vy(-1,j,1)*Bzint(-1,j,1,l)-Byint(-1,j,1,l)*Vz(-1,j,1)) 
            Eyint(-1,j,1,l)  = -(Vz(-1,j,1)*Bxint(-1,j,1,l)-Bzint(-1,j,1,l)*Vx(-1,j,1)) 
            Ezint(-1,j,1,l)  = -(Vx(-1,j,1)*Byint(-1,j,1,l)-Bxint(-1,j,1,l)*Vy(-1,j,1)) 

            Exint(-2,j,1,l)  = -(Vy(-2,j,1)*Bzint(-2,j,1,l)-Byint(-2,j,1,l)*Vz(-2,j,1)) 
            Eyint(-2,j,1,l)  = -(Vz(-2,j,1)*Bxint(-2,j,1,l)-Bzint(-2,j,1,l)*Vx(-2,j,1)) 
            Ezint(-2,j,1,l)  = -(Vx(-2,j,1)*Byint(-2,j,1,l)-Bxint(-2,j,1,l)*Vy(-2,j,1)) 

               else if (REC_PRIM == 1) then

            Exint(0,j,1,l)  = -(Vyint(0,j,1,l)*Bzint(0,j,1,l)-Byint(0,j,1,l)*Vzint(0,j,1,l)) 
            Eyint(0,j,1,l)  = -(Vzint(0,j,1,l)*Bxint(0,j,1,l)-Bzint(0,j,1,l)*Vxint(0,j,1,l)) 
            Ezint(0,j,1,l)  = -(Vxint(0,j,1,l)*Byint(0,j,1,l)-Bxint(0,j,1,l)*Vyint(0,j,1,l)) 

            Exint(-1,j,1,l)  = -(Vyint(-1,j,1,l)*Bzint(-1,j,1,l)-Byint(-1,j,1,l)*Vzint(-1,j,1,l)) 
            Eyint(-1,j,1,l)  = -(Vzint(-1,j,1,l)*Bxint(-1,j,1,l)-Bzint(-1,j,1,l)*Vxint(-1,j,1,l)) 
            Ezint(-1,j,1,l)  = -(Vxint(-1,j,1,l)*Byint(-1,j,1,l)-Bxint(-1,j,1,l)*Vyint(-1,j,1,l)) 

            Exint(-2,j,1,l)  = -(Vyint(-2,j,1,l)*Bzint(-2,j,1,l)-Byint(-2,j,1,l)*Vzint(-2,j,1,l)) 
            Eyint(-2,j,1,l)  = -(Vzint(-2,j,1,l)*Bxint(-2,j,1,l)-Bzint(-2,j,1,l)*Vxint(-2,j,1,l)) 
            Ezint(-2,j,1,l)  = -(Vxint(-2,j,1,l)*Byint(-2,j,1,l)-Bxint(-2,j,1,l)*Vyint(-2,j,1,l)) 

                 else 

            write(*,*) "STOP: subroutine boundary_electric"
            write(*,*) "Revise REC_PRIM parameter"
            stop

                 end if

          else

            Exint(0,j,1,l)  =  Exint(1,j,1,l)
            Eyint(0,j,1,l)  = -Eyint(1,j,1,l)
            Ezint(0,j,1,l)  = -Ezint(1,j,1,l)

            Exint(-1,j,1,l)  = Exint(0,j,1,l)
            Eyint(-1,j,1,l)  = Eyint(0,j,1,l)
            Ezint(-1,j,1,l)  = Ezint(0,j,1,l)

            Exint(-2,j,1,l)  = Exint(-1,j,1,l)
            Eyint(-2,j,1,l)  = Eyint(-1,j,1,l)
            Ezint(-2,j,1,l)  = Ezint(-1,j,1,l)

         end if
 

            !-----------------
            ! Right Boundary
            !---------------

            Exint(imax,j,1,l)  = Exint(imax-1,j,1,l)  
            Eyint(imax,j,1,l)  = Eyint(imax-1,j,1,l)  
            Ezint(imax,j,1,l)  = Ezint(imax-1,j,1,l)  

            Exint(imax+1,j,1,l)  = Exint(imax,j,1,l)  
            Eyint(imax+1,j,1,l)  = Eyint(imax,j,1,l)  
            Ezint(imax+1,j,1,l)  = Ezint(imax,j,1,l)  

            Exint(imax+2,j,1,l)  = Exint(imax+1,j,1,l)  
            Eyint(imax+2,j,1,l)  = Eyint(imax+1,j,1,l)  
            Ezint(imax+2,j,1,l)  = Ezint(imax+1,j,1,l)  

         end do

!$OMP END DO
         
  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
         
         do i=-6,imax+6

            ! Left Boundary
            !---------------

            Exint (i,0,1,l) = Exint(i,1,1,l)  
            Eyint (i,0,1,l) = Eyint(i,1,1,l)  
            Ezint (i,0,1,l) = Ezint(i,1,1,l)  

            Exint (i,-1,1,l) = Exint(i,0,1,l)  
            Eyint (i,-1,1,l) = Eyint(i,0,1,l)  
            Ezint (i,-1,1,l) = Ezint(i,0,1,l)  

            Exint (i,-2,1,l) = Exint(i,-1,1,l)  
            Eyint (i,-2,1,l) = Eyint(i,-1,1,l)  
            Ezint (i,-2,1,l) = Ezint(i,-1,1,l)  
            

            ! Right Boundary
            !---------------

            Exint (i,jmax,1,l) = Exint(i,jmax-1,1,l)  
            Eyint (i,jmax,1,l) = Eyint(i,jmax-1,1,l)  
            Ezint (i,jmax,1,l) = Ezint(i,jmax-1,1,l)  

            Exint (i,jmax+1,1,l) = Exint(i,jmax,1,l)  
            Eyint (i,jmax+1,1,l) = Eyint(i,jmax,1,l)  
            Ezint (i,jmax+1,1,l) = Ezint(i,jmax,1,l)  

            Exint (i,jmax+2,1,l) = Exint(i,jmax,1,l)  
            Eyint (i,jmax+2,1,l) = Eyint(i,jmax,1,l)  
            Ezint (i,jmax+2,1,l) = Ezint(i,jmax,1,l)  


         end do

!$OMP END DO
         
    else if (DIM == 2 .and. BOUND == 5) then

       ! Cloud Shock Interaction

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
      do j=-6,jmax+6

            !-----------------
            ! Left Boundary
            !---------------

            Exint(0,j,1,l)  = Exint(1,j,1,l) 
            Eyint(0,j,1,l)  = Eyint(1,j,1,l) 
            Ezint(0,j,1,l)  = Ezint(1,j,1,l) 
 
            Exint(-1,j,1,l)  = Exint(0,j,1,l) 
            Eyint(-1,j,1,l)  = Eyint(0,j,1,l) 
            Ezint(-1,j,1,l)  = Ezint(0,j,1,l) 

            Exint(-2,j,1,l)  = Exint(-1,j,1,l) 
            Eyint(-2,j,1,l)  = Eyint(-1,j,1,l) 
            Ezint(-2,j,1,l)  = Ezint(-1,j,1,l) 

            !-----------------
            ! Right Boundary
            !---------------

            if (REC_PRIM == 0) then

            Exint(imax,j,1,l)  = -(Vy(imax,j,1)*Bzint(imax,j,1,l)-Byint(imax,j,1,l)*Vz(imax,j,1)) !Exint(imax-1,j,1,l)  !
            Eyint(imax,j,1,l)  = -(Vz(imax,j,1)*Bxint(imax,j,1,l)-Bzint(imax,j,1,l)*Vx(imax,j,1)) !Eyint(imax-1,j,1,l)  !
            Ezint(imax,j,1,l)  = -(Vx(imax,j,1)*Byint(imax,j,1,l)-Bxint(imax,j,1,l)*Vy(imax,j,1)) !Ezint(imax-1,j,1,l)  !

            Exint(imax+1,j,1,l)  = Exint(imax,j,1,l)  !-(Vy(imax+1,j,1)*Bzint(imax+1,j,1,l)-Byint(imax+1,j,1,l)*Vz(imax+1,j,1)) 
            Eyint(imax+1,j,1,l)  = Eyint(imax,j,1,l)  !-(Vz(imax+1,j,1)*Bxint(imax+1,j,1,l)-Bzint(imax+1,j,1,l)*Vx(imax+1,j,1)) 
            Ezint(imax+1,j,1,l)  = Ezint(imax,j,1,l)  !-(Vx(imax+1,j,1)*Byint(imax+1,j,1,l)-Bxint(imax+1,j,1,l)*Vy(imax+1,j,1)) 

            Exint(imax+2,j,1,l)  = Exint(imax+1,j,1,l)  !-(Vy(imax+2,j,1)*Bzint(imax+2,j,1,l)-Byint(imax+2,j,1,l)*Vz(imax+2,j,1)) 
            Eyint(imax+2,j,1,l)  = Eyint(imax+1,j,1,l)  !-(Vz(imax+2,j,1)*Bxint(imax+2,j,1,l)-Bzint(imax+2,j,1,l)*Vx(imax+2,j,1)) 
            Ezint(imax+2,j,1,l)  = Ezint(imax+1,j,1,l)  !-(Vx(imax+2,j,1)*Byint(imax+2,j,1,l)-Bxint(imax+2,j,1,l)*Vy(imax+2,j,1)) 

               else if (REC_PRIM == 1) then

            Exint(imax,j,1,l)  = -(Vyint(imax,j,1,l)*Bzint(imax,j,1,l)-Byint(imax,j,1,l)*Vzint(imax,j,1,l)) !Exint(imax-1,j,1,l)  !
            Eyint(imax,j,1,l)  = -(Vzint(imax,j,1,l)*Bxint(imax,j,1,l)-Bzint(imax,j,1,l)*Vxint(imax,j,1,l)) !Eyint(imax-1,j,1,l)  !
            Ezint(imax,j,1,l)  = -(Vxint(imax,j,1,l)*Byint(imax,j,1,l)-Bxint(imax,j,1,l)*Vyint(imax,j,1,l)) !Ezint(imax-1,j,1,l)  !

            Exint(imax+1,j,1,l)  = Exint(imax,j,1,l)  !-(Vyint(imax+1,j,1,l)*Bzint(imax+1,j,1,l)-Byint(imax+1,j,1,l)*Vzint(imax+1,j,1,l)) 
            Eyint(imax+1,j,1,l)  = Eyint(imax,j,1,l)  !-(Vzint(imax+1,j,1,l)*Bxint(imax+1,j,1,l)-Bzint(imax+1,j,1,l)*Vxint(imax+1,j,1,l)) 
            Ezint(imax+1,j,1,l)  = Ezint(imax,j,1,l)  !-(Vxint(imax+1,j,1,l)*Byint(imax+1,j,1,l)-Bxint(imax+1,j,1,l)*Vyint(imax+1,j,1,l)) 

            Exint(imax+2,j,1,l)  = Exint(imax+1,j,1,l)  !-(Vyint(imax+2,j,1,l)*Bzint(imax+2,j,1,l)-Byint(imax+2,j,1,l)*Vzint(imax+2,j,1,l)) 
            Eyint(imax+2,j,1,l)  = Eyint(imax+1,j,1,l)  !-(Vzint(imax+2,j,1,l)*Bxint(imax+2,j,1,l)-Bzint(imax+2,j,1,l)*Vxint(imax+2,j,1,l)) 
            Ezint(imax+2,j,1,l)  = Ezint(imax+1,j,1,l)  !-(Vxint(imax+2,j,1,l)*Byint(imax+2,j,1,l)-Bxint(imax+2,j,1,l)*Vyint(imax+2,j,1,l)) 

                 else 

            write(*,*) "STOP: subroutine boundary_electric"
            write(*,*) "Revise REC_PRIM parameter"
            stop

                 end if


              end do

!$OMP END DO              


  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
              
         do i=-6,imax+6

            ! Left Boundary
            !---------------

            Exint (i,0,1,l) = Exint(i,1,1,l)  
            Eyint (i,0,1,l) = Eyint(i,1,1,l)  
            Ezint (i,0,1,l) = Ezint(i,1,1,l)  

            Exint (i,-1,1,l) = Exint(i,0,1,l)  
            Eyint (i,-1,1,l) = Eyint(i,0,1,l)  
            Ezint (i,-1,1,l) = Ezint(i,0,1,l)  

            Exint (i,-2,1,l) = Exint(i,-1,1,l)  
            Eyint (i,-2,1,l) = Eyint(i,-1,1,l)  
            Ezint (i,-2,1,l) = Ezint(i,-1,1,l)  
            

            ! Right Boundary
            !---------------

            Exint (i,jmax,1,l) = Exint(i,jmax-1,1,l)  
            Eyint (i,jmax,1,l) = Eyint(i,jmax-1,1,l)  
            Ezint (i,jmax,1,l) = Ezint(i,jmax-1,1,l)  

            Exint (i,jmax+1,1,l) = Exint(i,jmax,1,l)  
            Eyint (i,jmax+1,1,l) = Eyint(i,jmax,1,l)  
            Ezint (i,jmax+1,1,l) = Ezint(i,jmax,1,l)  

            Exint (i,jmax+2,1,l) = Exint(i,jmax+1,1,l)  
            Eyint (i,jmax+2,1,l) = Eyint(i,jmax+1,1,l)  
            Ezint (i,jmax+2,1,l) = Ezint(i,jmax+1,1,l)  

         end do

!$OMP END DO
         
    else if (DIM == 2 .and. BOUND == 6) then

       ! Tearing-Mode
       ! Richtmyer-Meshkov instability Boundary RM

       ! Periodic Boundary y-direction

  !-----------------
  ! Y Direction
  !-----------------

       ! Left Boundary
       ! -----------------

!$OMP DO  ORDERED           

       do i=-stencil,imax+stencil
          do j = 1, stencil

!$OMP ORDERED                

            Exint (i,0-j,1,l) = Exint (i,jmax-j,1,l)  
            Eyint (i,0-j,1,l) = Eyint (i,jmax-j,1,l)  
            Ezint (i,0-j,1,l) = Ezint (i,jmax-j,1,l)

!$OMP END ORDERED
             
          end do
       end do

!$OMP END DO
       
       ! Right  Boundary
       ! -----------------
       
!$OMP DO  ORDERED 

       do i=-stencil,imax+stencil
          do j = 1, stencil

!$OMP ORDERED
             
            Exint (i,jmax+j,1,l) = Exint (i,0+j,1,l)  
            Eyint (i,jmax+j,1,l) = Eyint (i,0+j,1,l)  
            Ezint (i,jmax+j,1,l) = Ezint (i,0+j,1,l)  

!$OMP END ORDERED
            
          end do
       end do

!$OMP END DO
         
       ! Copy Boundary x-direction

  !-----------------
  ! X Direction
  !-----------------

           ! Left  Boundary
           ! -----------------

!$OMP DO  ORDERED 
       
       do i = 1,stencil
          do j=-stencil,jmax+stencil

!$OMP ORDERED             

             Exint (0-i,j,1,l) = Exint ( 0,j,1,l)
             Eyint (0-i,j,1,l) = Eyint ( 0,j,1,l) 
             Ezint (0-i,j,1,l) = Ezint ( 0,j,1,l)

!$OMP END ORDERED             

        end do
     end do
     
!$OMP END DO

 
           ! Right  Boundary
           ! -----------------

!$OMP DO  ORDERED 
        
        do i = 1,stencil
           do j=-stencil,jmax+stencil

!$OMP ORDERED                  

            Exint (imax+i,j,1,l) = Exint (imax  ,j,1,l) 
            Eyint (imax+i,j,1,l) = Eyint (imax  ,j,1,l) 
            Ezint (imax+i,j,1,l) = Ezint (imax  ,j,1,l)

!$OMP END ORDERED               
            
        end do
     end do

!$OMP END DO


!****************************************************************


    else if (DIM == 2 .and. BOUND == 7) then

       ! Copy boundary 2D

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
      do j=-6,jmax+6

            !-----------------
            ! Left Boundary
            !---------------

            Exint(0,j,1,l)  = - Exint(1,j,1,l) 
            Eyint(0,j,1,l)  = - Eyint(1,j,1,l) 
            Ezint(0,j,1,l)  = - Ezint(1,j,1,l) 
 
            Exint(-1,j,1,l)  = Exint(0,j,1,l) 
            Eyint(-1,j,1,l)  = Eyint(0,j,1,l) 
            Ezint(-1,j,1,l)  = Ezint(0,j,1,l) 

            Exint(-2,j,1,l)  = Exint(-1,j,1,l) 
            Eyint(-2,j,1,l)  = Eyint(-1,j,1,l) 
            Ezint(-2,j,1,l)  = Ezint(-1,j,1,l) 

            !-----------------
            ! Right Boundary
            !---------------

            Exint(imax,j,1,l)    = Exint(imax-1,j,1,l)  
            Eyint(imax,j,1,l)    = Eyint(imax-1,j,1,l)  
            Ezint(imax,j,1,l)    = Ezint(imax-1,j,1,l)  

            Exint(imax+1,j,1,l)  = Exint(imax  ,j,1,l)  
            Eyint(imax+1,j,1,l)  = Eyint(imax  ,j,1,l)  
            Ezint(imax+1,j,1,l)  = Ezint(imax  ,j,1,l)  

            Exint(imax+2,j,1,l)  = Exint(imax+1,j,1,l)  
            Eyint(imax+2,j,1,l)  = Eyint(imax+1,j,1,l)  
            Ezint(imax+2,j,1,l)  = Ezint(imax+1,j,1,l)  

         end do

!$OMP END DO
         
  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
         
         do i=-6,imax+6


            ! Left Boundary
            !---------------

            Exint (i,0,1,l) = Exint(i,1,1,l)  
            Eyint (i,0,1,l) = Eyint(i,1,1,l)  
            Ezint (i,0,1,l) = Ezint(i,1,1,l)  

            Exint (i,-1,1,l) = Exint(i,0,1,l)  
            Eyint (i,-1,1,l) = Eyint(i,0,1,l)  
            Ezint (i,-1,1,l) = Ezint(i,0,1,l)  

            Exint (i,-2,1,l) = Exint(i,-1,1,l)  
            Eyint (i,-2,1,l) = Eyint(i,-1,1,l)  
            Ezint (i,-2,1,l) = Ezint(i,-1,1,l)  
            

            ! Right Boundary
            !---------------

            Exint (i,jmax,1,l)   = Exint(i,jmax-1,1,l)  
            Eyint (i,jmax,1,l)   = Eyint(i,jmax-1,1,l)  
            Ezint (i,jmax,1,l)   = Ezint(i,jmax-1,1,l)  

            Exint (i,jmax+1,1,l) = Exint(i,jmax  ,1,l)  
            Eyint (i,jmax+1,1,l) = Eyint(i,jmax  ,1,l)  
            Ezint (i,jmax+1,1,l) = Ezint(i,jmax  ,1,l)  

            Exint (i,jmax+2,1,l) = Exint(i,jmax+1,1,l)  
            Eyint (i,jmax+2,1,l) = Eyint(i,jmax+1,1,l)  
            Ezint (i,jmax+2,1,l) = Ezint(i,jmax+1,1,l)  

         end do

!$OMP END DO
         
     else

       write(*,*) "STOP: subroutine boundary_electric"
       write(*,*) "This combination of Dimension and Boundary is not correctly"
       write(*,*) "DIM =", DIM, "BOUND =", BOUND
       stop

    end if
   
!$OMP END PARALLEL


  end subroutine boundary_electric
