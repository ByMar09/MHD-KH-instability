   subroutine source_boundary

   use parameters
   use fourvectors

!$OMP PARALLEL PRIVATE(posy)
   
    if (DIM == 1 .and. BOUND == 1) then

      ! copy boundary 1D

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

!!$            Exastsour(0,1,1,l) = Exastsour(1,1,1,l) 
!!$            Eyastsour(0,1,1,l) = Eyastsour(1,1,1,l) 
!!$            Ezastsour(0,1,1,l) = Ezastsour(1,1,1,l)  


            Exastsour(-1,1,1,l) = Exastsour(0,1,1,l) 
            Eyastsour(-1,1,1,l) = Eyastsour(0,1,1,l) 
            Ezastsour(-1,1,1,l) = Ezastsour(0,1,1,l)  

            Exastsour(-2,1,1,l) = Exastsour(-1,1,1,l) 
            Eyastsour(-2,1,1,l) = Eyastsour(-1,1,1,l) 
            Ezastsour(-2,1,1,l) = Ezastsour(-1,1,1,l)

            Exastsour(-3,1,1,l) = Exastsour(-2,1,1,l) 
            Eyastsour(-3,1,1,l) = Eyastsour(-2,1,1,l) 
            Ezastsour(-3,1,1,l) = Ezastsour(-2,1,1,l)

            Exastsour(-4,1,1,l) = Exastsour(-3,1,1,l) 
            Eyastsour(-4,1,1,l) = Eyastsour(-3,1,1,l) 
            Ezastsour(-4,1,1,l) = Ezastsour(-3,1,1,l)

            Exastsour(-5,1,1,l) = Exastsour(-4,1,1,l) 
            Eyastsour(-5,1,1,l) = Eyastsour(-4,1,1,l) 
            Ezastsour(-5,1,1,l) = Ezastsour(-4,1,1,l)

            Exastsour(-6,1,1,l) = Exastsour(-5,1,1,l) 
            Eyastsour(-6,1,1,l) = Eyastsour(-5,1,1,l) 
            Ezastsour(-6,1,1,l) = Ezastsour(-5,1,1,l)

            

            ! Right Boundary
            !---------------

!!$            Exastsour(imax,1,1,l) = Exastsour(imax-1,1,1,l) 
!!$            Eyastsour(imax,1,1,l) = Eyastsour(imax-1,1,1,l) 
!!$            Ezastsour(imax,1,1,l) = Ezastsour(imax-1,1,1,l)  

            Exastsour(imax+1,1,1,l) = Exastsour(imax  ,1,1,l) 
            Eyastsour(imax+1,1,1,l) = Eyastsour(imax  ,1,1,l) 
            Ezastsour(imax+1,1,1,l) = Ezastsour(imax  ,1,1,l)  

            Exastsour(imax+2,1,1,l) = Exastsour(imax+1,1,1,l) 
            Exastsour(imax+2,1,1,l) = Exastsour(imax+1,1,1,l) 
            Exastsour(imax+2,1,1,l) = Exastsour(imax+1,1,1,l)

            Exastsour(imax+3,1,1,l) = Exastsour(imax+2,1,1,l) 
            Exastsour(imax+3,1,1,l) = Exastsour(imax+2,1,1,l) 
            Exastsour(imax+3,1,1,l) = Exastsour(imax+2,1,1,l)

            Exastsour(imax+4,1,1,l) = Exastsour(imax+3,1,1,l) 
            Exastsour(imax+4,1,1,l) = Exastsour(imax+3,1,1,l) 
            Exastsour(imax+4,1,1,l) = Exastsour(imax+3,1,1,l)

            Exastsour(imax+5,1,1,l) = Exastsour(imax+4,1,1,l) 
            Exastsour(imax+5,1,1,l) = Exastsour(imax+4,1,1,l) 
            Exastsour(imax+5,1,1,l) = Exastsour(imax+4,1,1,l)

            Exastsour(imax+6,1,1,l) = Exastsour(imax+5,1,1,l) 
            Exastsour(imax+6,1,1,l) = Exastsour(imax+5,1,1,l) 
            Exastsour(imax+6,1,1,l) = Exastsour(imax+5,1,1,l)  


    else if (DIM == 1 .and. BOUND == 2) then

      ! Periodic boundary 1D

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

!!$            Exastsour(0,1,1,l) = Exastsour(imax-1,1,1,l) 
!!$            Eyastsour(0,1,1,l) = Eyastsour(imax-1,1,1,l) 
!!$            Ezastsour(0,1,1,l) = Ezastsour(imax-1,1,1,l)  


            Exastsour(-1,1,1,l) = Exastsour(imax-1,1,1,l) 
            Eyastsour(-1,1,1,l) = Eyastsour(imax-1,1,1,l) 
            Ezastsour(-1,1,1,l) = Ezastsour(imax-1,1,1,l)  

            Exastsour(-2,1,1,l) = Exastsour(imax-2,1,1,l) 
            Eyastsour(-2,1,1,l) = Eyastsour(imax-2,1,1,l) 
            Ezastsour(-2,1,1,l) = Ezastsour(imax-2,1,1,l)

            Exastsour(-3,1,1,l) = Exastsour(imax-3,1,1,l) 
            Eyastsour(-3,1,1,l) = Eyastsour(imax-3,1,1,l) 
            Ezastsour(-3,1,1,l) = Ezastsour(imax-3,1,1,l)

            Exastsour(-4,1,1,l) = Exastsour(imax-4,1,1,l) 
            Eyastsour(-4,1,1,l) = Eyastsour(imax-4,1,1,l) 
            Ezastsour(-4,1,1,l) = Ezastsour(imax-4,1,1,l)

            Exastsour(-5,1,1,l) = Exastsour(imax-5,1,1,l) 
            Eyastsour(-5,1,1,l) = Eyastsour(imax-5,1,1,l) 
            Ezastsour(-5,1,1,l) = Ezastsour(imax-5,1,1,l)

            Exastsour(-6,1,1,l) = Exastsour(imax-6,1,1,l) 
            Eyastsour(-6,1,1,l) = Eyastsour(imax-6,1,1,l) 
            Ezastsour(-6,1,1,l) = Ezastsour(imax-6,1,1,l)

            ! Right Boundary
            !---------------

!!$            Exastsour(imax,1,1,l) = Exastsour(1,1,1,l) 
!!$            Eyastsour(imax,1,1,l) = Eyastsour(1,1,1,l) 
!!$            Ezastsour(imax,1,1,l) = Ezastsour(1,1,1,l)  

            Exastsour(imax+1,1,1,l) = Exastsour(1,1,1,l) 
            Eyastsour(imax+1,1,1,l) = Eyastsour(1,1,1,l) 
            Ezastsour(imax+1,1,1,l) = Ezastsour(1,1,1,l)  

            Exastsour(imax+2,1,1,l) = Exastsour(2,1,1,l) 
            Exastsour(imax+2,1,1,l) = Exastsour(2,1,1,l) 
            Exastsour(imax+2,1,1,l) = Exastsour(2,1,1,l)

            Exastsour(imax+3,1,1,l) = Exastsour(3,1,1,l) 
            Exastsour(imax+3,1,1,l) = Exastsour(3,1,1,l) 
            Exastsour(imax+3,1,1,l) = Exastsour(3,1,1,l)

            Exastsour(imax+4,1,1,l) = Exastsour(4,1,1,l) 
            Exastsour(imax+4,1,1,l) = Exastsour(4,1,1,l) 
            Exastsour(imax+4,1,1,l) = Exastsour(4,1,1,l)

            Exastsour(imax+5,1,1,l) = Exastsour(5,1,1,l) 
            Exastsour(imax+5,1,1,l) = Exastsour(5,1,1,l) 
            Exastsour(imax+5,1,1,l) = Exastsour(5,1,1,l)

            Exastsour(imax+6,1,1,l) = Exastsour(6,1,1,l) 
            Exastsour(imax+6,1,1,l) = Exastsour(6,1,1,l) 
            Exastsour(imax+6,1,1,l) = Exastsour(6,1,1,l)

    else if (DIM == 2 .and. BOUND == 1 .or. BOUND == 3) then

       ! Copy boundary 2D

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
     do j=-6,jmax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------

!!$            Exastsour( 0,j,1,l) = Exastsour( 1,j,1,l)
!!$            Eyastsour( 0,j,1,l) = Eyastsour( 1,j,1,l)
!!$            Ezastsour( 0,j,1,l) = Ezastsour( 1,j,1,l) 

            Exastsour(-1,j,1,l) = Exastsour( 0,j,1,l)
            Eyastsour(-1,j,1,l) = Eyastsour( 0,j,1,l)
            Ezastsour(-1,j,1,l) = Ezastsour( 0,j,1,l) 

            Exastsour(-2,j,1,l) = Exastsour(-1,j,1,l)
            Eyastsour(-2,j,1,l) = Eyastsour(-1,j,1,l)
            Ezastsour(-2,j,1,l) = Ezastsour(-1,j,1,l)

            Exastsour(-3,j,1,l) = Exastsour(-2,j,1,l)
            Eyastsour(-3,j,1,l) = Eyastsour(-2,j,1,l)
            Ezastsour(-3,j,1,l) = Ezastsour(-2,j,1,l)

            Exastsour(-4,j,1,l) = Exastsour(-3,j,1,l)
            Eyastsour(-4,j,1,l) = Eyastsour(-3,j,1,l)
            Ezastsour(-4,j,1,l) = Ezastsour(-3,j,1,l)

            Exastsour(-5,j,1,l) = Exastsour(-4,j,1,l)
            Eyastsour(-5,j,1,l) = Eyastsour(-4,j,1,l)
            Ezastsour(-5,j,1,l) = Ezastsour(-4,j,1,l)

            Exastsour(-6,j,1,l) = Exastsour(-5,j,1,l)
            Eyastsour(-6,j,1,l) = Eyastsour(-5,j,1,l)
            Ezastsour(-6,j,1,l) = Ezastsour(-5,j,1,l)


            ! Right Boundary
            !---------------

!!$            Exastsour(imax  ,j,1,l) = Exastsour(imax-1,j,1,l) 
!!$            Eyastsour(imax  ,j,1,l) = Eyastsour(imax-1,j,1,l) 
!!$            Ezastsour(imax  ,j,1,l) = Ezastsour(imax-1,j,1,l)  

            Exastsour(imax+1,j,1,l) = Exastsour(imax  ,j,1,l) 
            Eyastsour(imax+1,j,1,l) = Eyastsour(imax  ,j,1,l) 
            Ezastsour(imax+1,j,1,l) = Ezastsour(imax  ,j,1,l)  

            Exastsour(imax+2,j,1,l) = Exastsour(imax+1,j,1,l)  
            Eyastsour(imax+2,j,1,l) = Eyastsour(imax+1,j,1,l)  
            Ezastsour(imax+2,j,1,l) = Ezastsour(imax+1,j,1,l)
            
            Exastsour(imax+3,j,1,l) = Exastsour(imax+2,j,1,l)  
            Eyastsour(imax+3,j,1,l) = Eyastsour(imax+2,j,1,l)  
            Ezastsour(imax+3,j,1,l) = Ezastsour(imax+2,j,1,l)

            Exastsour(imax+4,j,1,l) = Exastsour(imax+3,j,1,l)  
            Eyastsour(imax+4,j,1,l) = Eyastsour(imax+3,j,1,l)  
            Ezastsour(imax+4,j,1,l) = Ezastsour(imax+3,j,1,l)

            Exastsour(imax+5,j,1,l) = Exastsour(imax+4,j,1,l)  
            Eyastsour(imax+5,j,1,l) = Eyastsour(imax+4,j,1,l)  
            Ezastsour(imax+5,j,1,l) = Ezastsour(imax+4,j,1,l)

            Exastsour(imax+6,j,1,l) = Exastsour(imax+5,j,1,l)  
            Eyastsour(imax+6,j,1,l) = Eyastsour(imax+5,j,1,l)  
            Ezastsour(imax+6,j,1,l) = Ezastsour(imax+5,j,1,l)  
            

!         end do
      end do

!$OMP END DO
      
  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
      
      do i=-6,imax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------

!!$            Exastsour(i, 0,1,l) = Exastsour(i, 1,1,l) 
!!$            Eyastsour(i, 0,1,l) = Eyastsour(i, 1,1,l) 
!!$            Ezastsour(i, 0,1,l) = Ezastsour(i, 1,1,l)  

            Exastsour(i,-1,1,l) = Exastsour(i, 0,1,l) 
            Eyastsour(i,-1,1,l) = Eyastsour(i, 0,1,l) 
            Ezastsour(i,-1,1,l) = Ezastsour(i, 0,1,l)  

            Exastsour(i,-2,1,l) = Exastsour(i,-1,1,l) 
            Eyastsour(i,-2,1,l) = Eyastsour(i,-1,1,l) 
            Ezastsour(i,-2,1,l) = Ezastsour(i,-1,1,l)

            Exastsour(i,-3,1,l) = Exastsour(i,-2,1,l) 
            Eyastsour(i,-3,1,l) = Eyastsour(i,-2,1,l) 
            Ezastsour(i,-3,1,l) = Ezastsour(i,-2,1,l)

            Exastsour(i,-4,1,l) = Exastsour(i,-3,1,l) 
            Eyastsour(i,-4,1,l) = Eyastsour(i,-3,1,l) 
            Ezastsour(i,-4,1,l) = Ezastsour(i,-3,1,l)

            Exastsour(i,-5,1,l) = Exastsour(i,-4,1,l) 
            Eyastsour(i,-5,1,l) = Eyastsour(i,-4,1,l) 
            Ezastsour(i,-5,1,l) = Ezastsour(i,-4,1,l)

            Exastsour(i,-6,1,l) = Exastsour(i,-5,1,l) 
            Eyastsour(i,-6,1,l) = Eyastsour(i,-5,1,l) 
            Ezastsour(i,-6,1,l) = Ezastsour(i,-5,1,l)  


            ! Right Boundary
            !---------------

!!$            Exastsour(i,jmax  ,1,l) = Exastsour(i,jmax-1,1,l) 
!!$            Eyastsour(i,jmax  ,1,l) = Eyastsour(i,jmax-1,1,l) 
!!$            Ezastsour(i,jmax  ,1,l) = Ezastsour(i,jmax-1,1,l)  

            Exastsour(i,jmax+1,1,l) = Exastsour(i,jmax  ,1,l) 
            Eyastsour(i,jmax+1,1,l) = Eyastsour(i,jmax  ,1,l) 
            Ezastsour(i,jmax+1,1,l) = Ezastsour(i,jmax  ,1,l)  

            Exastsour(i,jmax+2,1,l) = Exastsour(i,jmax+1,1,l) 
            Eyastsour(i,jmax+2,1,l) = Eyastsour(i,jmax+1,1,l) 
            Ezastsour(i,jmax+2,1,l) = Ezastsour(i,jmax+1,1,l)

            Exastsour(i,jmax+3,1,l) = Exastsour(i,jmax+2,1,l) 
            Eyastsour(i,jmax+3,1,l) = Eyastsour(i,jmax+2,1,l) 
            Ezastsour(i,jmax+3,1,l) = Ezastsour(i,jmax+2,1,l)

            Exastsour(i,jmax+4,1,l) = Exastsour(i,jmax+3,1,l) 
            Eyastsour(i,jmax+4,1,l) = Eyastsour(i,jmax+3,1,l) 
            Ezastsour(i,jmax+4,1,l) = Ezastsour(i,jmax+3,1,l)

            Exastsour(i,jmax+5,1,l) = Exastsour(i,jmax+4,1,l) 
            Eyastsour(i,jmax+5,1,l) = Eyastsour(i,jmax+4,1,l) 
            Ezastsour(i,jmax+5,1,l) = Ezastsour(i,jmax+4,1,l)

            Exastsour(i,jmax+6,1,l) = Exastsour(i,jmax+5,1,l) 
            Eyastsour(i,jmax+6,1,l) = Eyastsour(i,jmax+5,1,l) 
            Ezastsour(i,jmax+6,1,l) = Ezastsour(i,jmax+5,1,l)

!         end do
      end do

!$OMP END DO
      

    else if (DIM == 2 .and. BOUND == 2) then

       ! Periodic boundary 2D

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
     do j=-6,jmax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------

!!$            Exastsour( 0,j,1,l) = Exastsour(imax-1,j,1,l)
!!$            Eyastsour( 0,j,1,l) = Eyastsour(imax-1,j,1,l)
!!$            Ezastsour( 0,j,1,l) = Ezastsour(imax-1,j,1,l) 

            Exastsour(-1,j,1,l) = Exastsour(imax-1,j,1,l)
            Eyastsour(-1,j,1,l) = Eyastsour(imax-1,j,1,l)
            Ezastsour(-1,j,1,l) = Ezastsour(imax-1,j,1,l) 

            Exastsour(-2,j,1,l) = Exastsour(imax-2,j,1,l)
            Eyastsour(-2,j,1,l) = Eyastsour(imax-2,j,1,l)
            Ezastsour(-2,j,1,l) = Ezastsour(imax-2,j,1,l)

            Exastsour(-3,j,1,l) = Exastsour(imax-3,j,1,l)
            Eyastsour(-3,j,1,l) = Eyastsour(imax-3,j,1,l)
            Ezastsour(-3,j,1,l) = Ezastsour(imax-3,j,1,l)

            Exastsour(-4,j,1,l) = Exastsour(imax-4,j,1,l)
            Eyastsour(-4,j,1,l) = Eyastsour(imax-4,j,1,l)
            Ezastsour(-4,j,1,l) = Ezastsour(imax-4,j,1,l)

            Exastsour(-5,j,1,l) = Exastsour(imax-5,j,1,l)
            Eyastsour(-5,j,1,l) = Eyastsour(imax-5,j,1,l)
            Ezastsour(-5,j,1,l) = Ezastsour(imax-5,j,1,l)

            Exastsour(-6,j,1,l) = Exastsour(imax-6,j,1,l)
            Eyastsour(-6,j,1,l) = Eyastsour(imax-6,j,1,l)
            Ezastsour(-6,j,1,l) = Ezastsour(imax-6,j,1,l) 


            ! Right Boundary
            !---------------

!!$            Exastsour(imax  ,j,1,l) = Exastsour(1,j,1,l)  
!!$            Eyastsour(imax  ,j,1,l) = Eyastsour(1,j,1,l)  
!!$            Ezastsour(imax  ,j,1,l) = Ezastsour(1,j,1,l)  

            Exastsour(imax+1,j,1,l) = Exastsour(1,j,1,l) 
            Eyastsour(imax+1,j,1,l) = Eyastsour(1,j,1,l) 
            Ezastsour(imax+1,j,1,l) = Ezastsour(1,j,1,l)  
 
            Exastsour(imax+2,j,1,l) = Exastsour(2,j,1,l)  
            Eyastsour(imax+2,j,1,l) = Eyastsour(2,j,1,l)  
            Ezastsour(imax+2,j,1,l) = Ezastsour(2,j,1,l)

            Exastsour(imax+3,j,1,l) = Exastsour(3,j,1,l)  
            Eyastsour(imax+3,j,1,l) = Eyastsour(3,j,1,l)  
            Ezastsour(imax+3,j,1,l) = Ezastsour(3,j,1,l)

            Exastsour(imax+4,j,1,l) = Exastsour(4,j,1,l)  
            Eyastsour(imax+4,j,1,l) = Eyastsour(4,j,1,l)  
            Ezastsour(imax+4,j,1,l) = Ezastsour(4,j,1,l)

            Exastsour(imax+5,j,1,l) = Exastsour(5,j,1,l)  
            Eyastsour(imax+5,j,1,l) = Eyastsour(5,j,1,l)  
            Ezastsour(imax+5,j,1,l) = Ezastsour(5,j,1,l)

            Exastsour(imax+6,j,1,l) = Exastsour(6,j,1,l)  
            Eyastsour(imax+6,j,1,l) = Eyastsour(6,j,1,l)  
            Ezastsour(imax+6,j,1,l) = Ezastsour(6,j,1,l) 

!         end do
      end do

!$OMP END DO
      
  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
      
      do i=-6,imax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------

!!$            Exastsour(i, 0,1,l) = Exastsour(i,jmax-1,1,l) 
!!$            Eyastsour(i, 0,1,l) = Eyastsour(i,jmax-1,1,l) 
!!$            Ezastsour(i, 0,1,l) = Ezastsour(i,jmax-1,1,l)  

            Exastsour(i,-1,1,l) = Exastsour(i,jmax-1,1,l) 
            Eyastsour(i,-1,1,l) = Eyastsour(i,jmax-1,1,l) 
            Ezastsour(i,-1,1,l) = Ezastsour(i,jmax-1,1,l)  

            Exastsour(i,-2,1,l) = Exastsour(i,jmax-2,1,l) 
            Eyastsour(i,-2,1,l) = Eyastsour(i,jmax-2,1,l) 
            Ezastsour(i,-2,1,l) = Ezastsour(i,jmax-2,1,l)

            Exastsour(i,-3,1,l) = Exastsour(i,jmax-3,1,l) 
            Eyastsour(i,-3,1,l) = Eyastsour(i,jmax-3,1,l) 
            Ezastsour(i,-3,1,l) = Ezastsour(i,jmax-3,1,l)

            Exastsour(i,-4,1,l) = Exastsour(i,jmax-4,1,l) 
            Eyastsour(i,-4,1,l) = Eyastsour(i,jmax-4,1,l) 
            Ezastsour(i,-4,1,l) = Ezastsour(i,jmax-4,1,l)

            Exastsour(i,-5,1,l) = Exastsour(i,jmax-5,1,l) 
            Eyastsour(i,-5,1,l) = Eyastsour(i,jmax-5,1,l) 
            Ezastsour(i,-5,1,l) = Ezastsour(i,jmax-5,1,l)

            Exastsour(i,-6,1,l) = Exastsour(i,jmax-6,1,l) 
            Eyastsour(i,-6,1,l) = Eyastsour(i,jmax-6,1,l) 
            Ezastsour(i,-6,1,l) = Ezastsour(i,jmax-6,1,l) 

            ! Right Boundary
            !---------------

!!$            Exastsour(i,jmax  ,1,l) = Exastsour(i,1,1,l) 
!!$            Eyastsour(i,jmax  ,1,l) = Eyastsour(i,1,1,l) 
!!$            Ezastsour(i,jmax  ,1,l) = Ezastsour(i,1,1,l)  

            Exastsour(i,jmax+1,1,l) = Exastsour(i,1,1,l) 
            Eyastsour(i,jmax+1,1,l) = Eyastsour(i,1,1,l) 
            Ezastsour(i,jmax+1,1,l) = Ezastsour(i,1,1,l)  
 
            Exastsour(i,jmax+2,1,l) = Exastsour(i,2,1,l) 
            Eyastsour(i,jmax+2,1,l) = Eyastsour(i,2,1,l) 
            Ezastsour(i,jmax+2,1,l) = Ezastsour(i,2,1,l)

            Exastsour(i,jmax+3,1,l) = Exastsour(i,3,1,l) 
            Eyastsour(i,jmax+3,1,l) = Eyastsour(i,3,1,l) 
            Ezastsour(i,jmax+3,1,l) = Ezastsour(i,3,1,l) 

            Exastsour(i,jmax+4,1,l) = Exastsour(i,4,1,l) 
            Eyastsour(i,jmax+4,1,l) = Eyastsour(i,4,1,l) 
            Ezastsour(i,jmax+4,1,l) = Ezastsour(i,4,1,l)

            Exastsour(i,jmax+5,1,l) = Exastsour(i,5,1,l) 
            Eyastsour(i,jmax+5,1,l) = Eyastsour(i,5,1,l) 
            Ezastsour(i,jmax+5,1,l) = Ezastsour(i,5,1,l)

            Exastsour(i,jmax+6,1,l) = Exastsour(i,6,1,l) 
            Eyastsour(i,jmax+6,1,l) = Eyastsour(i,6,1,l) 
            Ezastsour(i,jmax+6,1,l) = Ezastsour(i,6,1,l) 
            
            
!         end do
      end do

!$OMP END DO
      
    else if (DIM == 2 .and. BOUND == 4) then

       ! Slab Jet

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
     do j=-6,jmax+6
!      do k=1,kmax

           posy  = - 14.d0 + j * Dely
           
           ! Left Boundary
           ! -----------------

           if (posy .gt. -1.d0 .and. posy .le. 1.d0) then


            Exastsour( 0,j,1,l) = Exastsour( 1,j,1,l)
            Eyastsour( 0,j,1,l) = Eyastsour( 1,j,1,l)
            Ezastsour( 0,j,1,l) = Ezastsour( 1,j,1,l) 

            Exastsour(-1,j,1,l) = Exastsour( 0,j,1,l)
            Eyastsour(-1,j,1,l) = Eyastsour( 0,j,1,l)
            Ezastsour(-1,j,1,l) = Ezastsour( 0,j,1,l) 

            Exastsour(-2,j,1,l) = Exastsour(-1,j,1,l)
            Eyastsour(-2,j,1,l) = Eyastsour(-1,j,1,l)
            Ezastsour(-2,j,1,l) = Ezastsour(-1,j,1,l) 


           else 


            Exastsour( 0,j,1,l) = Exastsour( 1,j,1,l)
            Eyastsour( 0,j,1,l) = Eyastsour( 1,j,1,l)
            Ezastsour( 0,j,1,l) = Ezastsour( 1,j,1,l) 

            Exastsour(-1,j,1,l) = Exastsour( 0,j,1,l)
            Eyastsour(-1,j,1,l) = Eyastsour( 0,j,1,l)
            Ezastsour(-1,j,1,l) = Ezastsour( 0,j,1,l) 

            Exastsour(-2,j,1,l) = Exastsour(-1,j,1,l)
            Eyastsour(-2,j,1,l) = Eyastsour(-1,j,1,l)
            Ezastsour(-2,j,1,l) = Ezastsour(-1,j,1,l) 


            end if


         ! Right Boundary
         !---------------


            Exastsour(imax  ,j,1,l) = Exastsour(imax-1,j,1,l) 
            Eyastsour(imax  ,j,1,l) = Eyastsour(imax-1,j,1,l) 
            Ezastsour(imax  ,j,1,l) = Ezastsour(imax-1,j,1,l)  

            Exastsour(imax+1,j,1,l) = Exastsour(imax,j,1,l) 
            Eyastsour(imax+1,j,1,l) = Eyastsour(imax,j,1,l) 
            Ezastsour(imax+1,j,1,l) = Ezastsour(imax,j,1,l)  

            Exastsour(imax+2,j,1,l) = Exastsour(imax,j,1,l)  
            Eyastsour(imax+2,j,1,l) = Eyastsour(imax,j,1,l)  
            Ezastsour(imax+2,j,1,l) = Ezastsour(imax,j,1,l)  

!         end do
      end do

!$OMP END DO
      
  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
      
      do i=-6,imax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------

            Exastsour(i, 0,1,l) = Exastsour(i, 1,1,l) 
            Eyastsour(i, 0,1,l) = Eyastsour(i, 1,1,l) 
            Ezastsour(i, 0,1,l) = Ezastsour(i, 1,1,l)  

            Exastsour(i,-1,1,l) = Exastsour(i, 0,1,l) 
            Eyastsour(i,-1,1,l) = Eyastsour(i, 0,1,l) 
            Ezastsour(i,-1,1,l) = Ezastsour(i, 0,1,l)  

            Exastsour(i,-2,1,l) = Exastsour(i,-1,1,l) 
            Eyastsour(i,-2,1,l) = Eyastsour(i,-1,1,l) 
            Ezastsour(i,-2,1,l) = Ezastsour(i,-1,1,l)  


            ! Right Boundary
            !---------------

            Exastsour(i,jmax  ,1,l) = Exastsour(i,jmax-1,1,l) 
            Eyastsour(i,jmax  ,1,l) = Eyastsour(i,jmax-1,1,l) 
            Ezastsour(i,jmax  ,1,l) = Ezastsour(i,jmax-1,1,l)  

            Exastsour(i,jmax+1,1,l) = Exastsour(i,jmax  ,1,l) 
            Eyastsour(i,jmax+1,1,l) = Eyastsour(i,jmax  ,1,l) 
            Ezastsour(i,jmax+1,1,l) = Ezastsour(i,jmax  ,1,l)  

            Exastsour(i,jmax+2,1,l) = Exastsour(i,jmax+1,1,l) 
            Eyastsour(i,jmax+2,1,l) = Eyastsour(i,jmax+1,1,l) 
            Ezastsour(i,jmax+2,1,l) = Ezastsour(i,jmax+1,1,l)  

!         end do
      end do

!$OMP END DO
      
    else if (DIM == 2 .and. BOUND == 5) then

       ! Cloud Shock Interction

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
     do j=-6,jmax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------

            Exastsour( 0,j,1,l) = Exastsour( 1,j,1,l)
            Eyastsour( 0,j,1,l) = Eyastsour( 1,j,1,l)
            Ezastsour( 0,j,1,l) = Ezastsour( 1,j,1,l) 

            Exastsour(-1,j,1,l) = Exastsour( 0,j,1,l)
            Eyastsour(-1,j,1,l) = Eyastsour( 0,j,1,l)
            Ezastsour(-1,j,1,l) = Ezastsour( 0,j,1,l) 

            Exastsour(-2,j,1,l) = Exastsour(-1,j,1,l)
            Eyastsour(-2,j,1,l) = Eyastsour(-1,j,1,l)
            Ezastsour(-2,j,1,l) = Ezastsour(-1,j,1,l) 


            ! Right Boundary
            !---------------


            Exastsour(imax  ,j,1,l) = Exastsour(imax-1,j,1,l) 
            Eyastsour(imax  ,j,1,l) = Eyastsour(imax-1,j,1,l) 
            Ezastsour(imax  ,j,1,l) = Ezastsour(imax-1,j,1,l)  

            Exastsour(imax+1,j,1,l) = Exastsour(imax  ,j,1,l) 
            Eyastsour(imax+1,j,1,l) = Eyastsour(imax  ,j,1,l) 
            Ezastsour(imax+1,j,1,l) = Ezastsour(imax  ,j,1,l)  

            Exastsour(imax+2,j,1,l) = Exastsour(imax+1,j,1,l)  
            Eyastsour(imax+2,j,1,l) = Eyastsour(imax+1,j,1,l)  
            Ezastsour(imax+2,j,1,l) = Ezastsour(imax+1,j,1,l)  

!         end do
      end do

!$OMP END DO
      
  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
      
      do i=-6,imax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------


            Exastsour(i, 0,1,l) = Exastsour(i, 1,1,l) 
            Eyastsour(i, 0,1,l) = Eyastsour(i, 1,1,l) 
            Ezastsour(i, 0,1,l) = Ezastsour(i, 1,1,l)  

            Exastsour(i,-1,1,l) = Exastsour(i, 0,1,l) 
            Eyastsour(i,-1,1,l) = Eyastsour(i, 0,1,l) 
            Ezastsour(i,-1,1,l) = Ezastsour(i, 0,1,l)  

            Exastsour(i,-2,1,l) = Exastsour(i,-1,1,l) 
            Eyastsour(i,-2,1,l) = Eyastsour(i,-1,1,l) 
            Ezastsour(i,-2,1,l) = Ezastsour(i,-1,1,l)  

            ! Right Boundary
            !---------------

            Exastsour(i,jmax  ,1,l) = Exastsour(i,jmax-1,1,l) 
            Eyastsour(i,jmax  ,1,l) = Eyastsour(i,jmax-1,1,l) 
            Ezastsour(i,jmax  ,1,l) = Ezastsour(i,jmax-1,1,l)  

            Exastsour(i,jmax+1,1,l) = Exastsour(i,jmax  ,1,l) 
            Eyastsour(i,jmax+1,1,l) = Eyastsour(i,jmax  ,1,l) 
            Ezastsour(i,jmax+1,1,l) = Ezastsour(i,jmax  ,1,l)  

            Exastsour(i,jmax+2,1,l) = Exastsour(i,jmax+1,1,l) 
            Eyastsour(i,jmax+2,1,l) = Eyastsour(i,jmax+1,1,l) 
            Ezastsour(i,jmax+2,1,l) = Ezastsour(i,jmax+1,1,l) 

!         end do
      end do

!$OMP END DO
      
    else if (DIM == 2 .and. BOUND == 6) then

       ! Tearing Mode
       ! Richtmyer-Meshkov instability Boundary RM

       ! Copy boundary x-direction

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
    do j=-6,jmax+6
!       do k=1,kmax

           ! Left Boundary
           ! -----------------

!!$           Exastsour( 0,j,1,l) = Exastsour(1,j,1,l) 
!!$           Eyastsour( 0,j,1,l) = Eyastsour(1,j,1,l) 
!!$           Ezastsour( 0,j,1,l) = Ezastsour(1,j,1,l) 

           Exastsour(-1,j,1,l) = Exastsour(0,j,1,l) 
           Eyastsour(-1,j,1,l) = Eyastsour(0,j,1,l) 
           Ezastsour(-1,j,1,l) = Ezastsour(0,j,1,l) 

           Exastsour(-2,j,1,l) = Exastsour(-1,j,1,l) 
           Eyastsour(-2,j,1,l) = Eyastsour(-1,j,1,l) 
           Ezastsour(-2,j,1,l) = Ezastsour(-1,j,1,l)

           Exastsour(-3,j,1,l) = Exastsour(-2,j,1,l) 
           Eyastsour(-3,j,1,l) = Eyastsour(-2,j,1,l) 
           Ezastsour(-3,j,1,l) = Ezastsour(-2,j,1,l)

           Exastsour(-4,j,1,l) = Exastsour(-3,j,1,l) 
           Eyastsour(-4,j,1,l) = Eyastsour(-3,j,1,l) 
           Ezastsour(-4,j,1,l) = Ezastsour(-3,j,1,l)

           Exastsour(-5,j,1,l) = Exastsour(-4,j,1,l) 
           Eyastsour(-5,j,1,l) = Eyastsour(-4,j,1,l) 
           Ezastsour(-5,j,1,l) = Ezastsour(-4,j,1,l)

           Exastsour(-6,j,1,l) = Exastsour(-5,j,1,l) 
           Eyastsour(-6,j,1,l) = Eyastsour(-5,j,1,l) 
           Ezastsour(-6,j,1,l) = Ezastsour(-5,j,1,l) 

           ! Right  Boundary
           ! -----------------

!!$           Exastsour(imax  ,j,1,l) = Exastsour(imax-1,j,1,l) 
!!$           Eyastsour(imax  ,j,1,l) = Eyastsour(imax-1,j,1,l) 
!!$           Ezastsour(imax  ,j,1,l) = Ezastsour(imax-1,j,1,l) 

           Exastsour(imax+1,j,1,l) = Exastsour(imax ,j,1,l) 
           Eyastsour(imax+1,j,1,l) = Eyastsour(imax ,j,1,l) 
           Ezastsour(imax+1,j,1,l) = Ezastsour(imax ,j,1,l) 

           Exastsour(imax+2,j,1,l) = Exastsour(imax+1,j,1,l)
           Eyastsour(imax+2,j,1,l) = Eyastsour(imax+1,j,1,l)
           Ezastsour(imax+2,j,1,l) = Ezastsour(imax+1,j,1,l)

           Exastsour(imax+3,j,1,l) = Exastsour(imax+2,j,1,l)
           Eyastsour(imax+3,j,1,l) = Eyastsour(imax+2,j,1,l)
           Ezastsour(imax+3,j,1,l) = Ezastsour(imax+2,j,1,l)

           Exastsour(imax+4,j,1,l) = Exastsour(imax+3,j,1,l)
           Eyastsour(imax+4,j,1,l) = Eyastsour(imax+3,j,1,l)
           Ezastsour(imax+4,j,1,l) = Ezastsour(imax+3,j,1,l)

           Exastsour(imax+5,j,1,l) = Exastsour(imax+4,j,1,l)
           Eyastsour(imax+5,j,1,l) = Eyastsour(imax+4,j,1,l)
           Ezastsour(imax+5,j,1,l) = Ezastsour(imax+4,j,1,l)

           Exastsour(imax+6,j,1,l) = Exastsour(imax+5,j,1,l)
           Eyastsour(imax+6,j,1,l) = Eyastsour(imax+5,j,1,l)
           Ezastsour(imax+6,j,1,l) = Ezastsour(imax+5,j,1,l)

           
!        end do
     end do

!$OMP END DO
     
     ! Periodic Boundary y-direction

  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
     
    do i=-6,imax+6
!       do k=1,kmax

           ! Left Boundary
           ! -----------------

!!$           Exastsour(i, 0,1,l) = Exastsour(i,jmax-1,1,l)
!!$           Eyastsour(i, 0,1,l) = Eyastsour(i,jmax-1,1,l)
!!$           Ezastsour(i, 0,1,l) = Ezastsour(i,jmax-1,1,l)

           Exastsour(i,-1,1,l) = Exastsour(i,jmax-1,1,l)
           Eyastsour(i,-1,1,l) = Eyastsour(i,jmax-1,1,l)
           Ezastsour(i,-1,1,l) = Ezastsour(i,jmax-1,1,l)

           Exastsour(i,-2,1,l) = Exastsour(i,jmax-2,1,l)
           Eyastsour(i,-2,1,l) = Eyastsour(i,jmax-2,1,l)
           Ezastsour(i,-2,1,l) = Ezastsour(i,jmax-2,1,l)

           Exastsour(i,-3,1,l) = Exastsour(i,jmax-3,1,l)
           Eyastsour(i,-3,1,l) = Eyastsour(i,jmax-3,1,l)
           Ezastsour(i,-3,1,l) = Ezastsour(i,jmax-3,1,l)

           Exastsour(i,-4,1,l) = Exastsour(i,jmax-4,1,l)
           Eyastsour(i,-4,1,l) = Eyastsour(i,jmax-4,1,l)
           Ezastsour(i,-4,1,l) = Ezastsour(i,jmax-4,1,l)

           Exastsour(i,-5,1,l) = Exastsour(i,jmax-5,1,l)
           Eyastsour(i,-5,1,l) = Eyastsour(i,jmax-5,1,l)
           Ezastsour(i,-5,1,l) = Ezastsour(i,jmax-5,1,l)

           Exastsour(i,-6,1,l) = Exastsour(i,jmax-6,1,l)
           Eyastsour(i,-6,1,l) = Eyastsour(i,jmax-6,1,l)
           Ezastsour(i,-6,1,l) = Ezastsour(i,jmax-6,1,l)


           ! Right  Boundary
           ! -----------------

!!$           Exastsour(i,jmax  ,1,l) = Exastsour(i,1,1,l) 
!!$           Eyastsour(i,jmax  ,1,l) = Eyastsour(i,1,1,l) 
!!$           Ezastsour(i,jmax  ,1,l) = Ezastsour(i,1,1,l) 

           Exastsour(i,jmax+1,1,l) = Exastsour(i,1,1,l) 
           Eyastsour(i,jmax+1,1,l) = Eyastsour(i,1,1,l) 
           Ezastsour(i,jmax+1,1,l) = Ezastsour(i,1,1,l) 

           Exastsour(i,jmax+2,1,l) = Exastsour(i,2,1,l) 
           Eyastsour(i,jmax+2,1,l) = Eyastsour(i,2,1,l) 
           Ezastsour(i,jmax+2,1,l) = Ezastsour(i,2,1,l)

           Exastsour(i,jmax+3,1,l) = Exastsour(i,3,1,l) 
           Eyastsour(i,jmax+3,1,l) = Eyastsour(i,3,1,l) 
           Ezastsour(i,jmax+3,1,l) = Ezastsour(i,3,1,l)

           Exastsour(i,jmax+4,1,l) = Exastsour(i,4,1,l) 
           Eyastsour(i,jmax+4,1,l) = Eyastsour(i,4,1,l) 
           Ezastsour(i,jmax+4,1,l) = Ezastsour(i,4,1,l)

           Exastsour(i,jmax+5,1,l) = Exastsour(i,5,1,l) 
           Eyastsour(i,jmax+5,1,l) = Eyastsour(i,5,1,l) 
           Ezastsour(i,jmax+5,1,l) = Ezastsour(i,5,1,l)

           Exastsour(i,jmax+6,1,l) = Exastsour(i,6,1,l) 
           Eyastsour(i,jmax+6,1,l) = Eyastsour(i,6,1,l) 
           Ezastsour(i,jmax+6,1,l) = Ezastsour(i,6,1,l) 

!        end do
     end do

!$OMP END DO
     
!**********************************************************************

    else if (DIM == 2 .and. BOUND == 7) then

       ! Copy boundary 2D

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 
       
     do j=-6,jmax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------

            Exastsour( 0,j,1,l) = - Exastsour( 1,j,1,l)
            Eyastsour( 0,j,1,l) = - Eyastsour( 1,j,1,l)
            Ezastsour( 0,j,1,l) = - Ezastsour( 1,j,1,l) 

            Exastsour(-1,j,1,l) = Exastsour( 0,j,1,l)
            Eyastsour(-1,j,1,l) = Eyastsour( 0,j,1,l)
            Ezastsour(-1,j,1,l) = Ezastsour( 0,j,1,l) 

            Exastsour(-2,j,1,l) = Exastsour(-1,j,1,l)
            Eyastsour(-2,j,1,l) = Eyastsour(-1,j,1,l)
            Ezastsour(-2,j,1,l) = Ezastsour(-1,j,1,l) 


            ! Right Boundary
            !---------------

            Exastsour(imax  ,j,1,l) = Exastsour(imax-1,j,1,l) 
            Eyastsour(imax  ,j,1,l) = Eyastsour(imax-1,j,1,l) 
            Ezastsour(imax  ,j,1,l) = Ezastsour(imax-1,j,1,l)  

            Exastsour(imax+1,j,1,l) = Exastsour(imax  ,j,1,l) 
            Eyastsour(imax+1,j,1,l) = Eyastsour(imax  ,j,1,l) 
            Ezastsour(imax+1,j,1,l) = Ezastsour(imax  ,j,1,l)  

            Exastsour(imax+2,j,1,l) = Exastsour(imax+1,j,1,l)  
            Eyastsour(imax+2,j,1,l) = Eyastsour(imax+1,j,1,l)  
            Ezastsour(imax+2,j,1,l) = Ezastsour(imax+1,j,1,l)  

!         end do
      end do

!$OMP END DO
      
  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 
      
      do i=-6,imax+6
!      do k=1,kmax

            ! Left Boundary
            !---------------

            Exastsour(i, 0,1,l) = Exastsour(i, 1,1,l) 
            Eyastsour(i, 0,1,l) = Eyastsour(i, 1,1,l) 
            Ezastsour(i, 0,1,l) = Ezastsour(i, 1,1,l)  

            Exastsour(i,-1,1,l) = Exastsour(i, 0,1,l) 
            Eyastsour(i,-1,1,l) = Eyastsour(i, 0,1,l) 
            Ezastsour(i,-1,1,l) = Ezastsour(i, 0,1,l)  

            Exastsour(i,-2,1,l) = Exastsour(i,-1,1,l) 
            Eyastsour(i,-2,1,l) = Eyastsour(i,-1,1,l) 
            Ezastsour(i,-2,1,l) = Ezastsour(i,-1,1,l)  


            ! Right Boundary
            !---------------

            Exastsour(i,jmax  ,1,l) = Exastsour(i,jmax-1,1,l) 
            Eyastsour(i,jmax  ,1,l) = Eyastsour(i,jmax-1,1,l) 
            Ezastsour(i,jmax  ,1,l) = Ezastsour(i,jmax-1,1,l)  

            Exastsour(i,jmax+1,1,l) = Exastsour(i,jmax  ,1,l) 
            Eyastsour(i,jmax+1,1,l) = Eyastsour(i,jmax  ,1,l) 
            Ezastsour(i,jmax+1,1,l) = Ezastsour(i,jmax  ,1,l)  

            Exastsour(i,jmax+2,1,l) = Exastsour(i,jmax+1,1,l) 
            Eyastsour(i,jmax+2,1,l) = Eyastsour(i,jmax+1,1,l) 
            Ezastsour(i,jmax+2,1,l) = Ezastsour(i,jmax+1,1,l) 

!         end do
      end do

!$OMP END DO
      
    else

       write(*,*) "STOP: subroutine source_boundary"
       write(*,*) "This combination of Dimension and Boundary is not correctly"
       write(*,*) "DIM =", DIM, "BOUND =", BOUND
       stop

    end if

!$OMP END PARALLEL
    
   end subroutine source_boundary
