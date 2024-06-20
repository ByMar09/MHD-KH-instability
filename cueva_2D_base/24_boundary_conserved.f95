
  subroutine boundary_conserved

    use scalar
    use parameters
    use threevectors 
    use fourvectors

    implicit none

!$OMP PARALLEL
 
    if (DIM == 1 .and. BOUND == 1) then

       ! Copy boundary 1D

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

!!$            psiint( 0,1,1,l) = psiint(1,1,1,l) 
!!$            phiint( 0,1,1,l) = phiint(1,1,1,l) 
!!$            Bxint ( 0,1,1,l) = Bxint (1,1,1,l)  
!!$            Byint ( 0,1,1,l) = Byint (1,1,1,l)  
!!$            Bzint ( 0,1,1,l) = Bzint (1,1,1,l)  
!!$            qint  ( 0,1,1,l) = qint  (1,1,1,l)   
!!$            DDint ( 0,1,1,l) = DDint (1,1,1,l)  
!!$            tauint( 0,1,1,l) = tauint(1,1,1,l) 
!!$            Sxint ( 0,1,1,l) = Sxint (1,1,1,l)  
!!$            Syint ( 0,1,1,l) = Syint (1,1,1,l)  
!!$            Szint ( 0,1,1,l) = Szint (1,1,1,l)  


            psiint(-1,1,1,l) = psiint(0,1,1,l) 
            phiint(-1,1,1,l) = phiint(0,1,1,l) 
            Bxint (-1,1,1,l) = Bxint (0,1,1,l)  
            Byint (-1,1,1,l) = Byint (0,1,1,l)  
            Bzint (-1,1,1,l) = Bzint (0,1,1,l)  
            qint  (-1,1,1,l) = qint  (0,1,1,l)   
            DDint (-1,1,1,l) = DDint (0,1,1,l)  
            tauint(-1,1,1,l) = tauint(0,1,1,l) 
            Sxint (-1,1,1,l) = Sxint (0,1,1,l)  
            Syint (-1,1,1,l) = Syint (0,1,1,l)  
            Szint (-1,1,1,l) = Szint (0,1,1,l)  


            psiint(-2,1,1,l) = psiint(-1,1,1,l) 
            phiint(-2,1,1,l) = phiint(-1,1,1,l) 
            Bxint (-2,1,1,l) = Bxint (-1,1,1,l)  
            Byint (-2,1,1,l) = Byint (-1,1,1,l)  
            Bzint (-2,1,1,l) = Bzint (-1,1,1,l)  
            qint  (-2,1,1,l) = qint  (-1,1,1,l)   
            DDint (-2,1,1,l) = DDint (-1,1,1,l)  
            tauint(-2,1,1,l) = tauint(-1,1,1,l) 
            Sxint (-2,1,1,l) = Sxint (-1,1,1,l)  
            Syint (-2,1,1,l) = Syint (-1,1,1,l)  
            Szint (-2,1,1,l) = Szint (-1,1,1,l)


            psiint(-3,1,1,l) = psiint(-2,1,1,l) 
            phiint(-3,1,1,l) = phiint(-2,1,1,l) 
            Bxint (-3,1,1,l) = Bxint (-2,1,1,l)  
            Byint (-3,1,1,l) = Byint (-2,1,1,l)  
            Bzint (-3,1,1,l) = Bzint (-2,1,1,l)  
            qint  (-3,1,1,l) = qint  (-2,1,1,l)   
            DDint (-3,1,1,l) = DDint (-2,1,1,l)  
            tauint(-3,1,1,l) = tauint(-2,1,1,l) 
            Sxint (-3,1,1,l) = Sxint (-2,1,1,l)  
            Syint (-3,1,1,l) = Syint (-2,1,1,l)  
            Szint (-3,1,1,l) = Szint (-2,1,1,l) 


            psiint(-4,1,1,l) = psiint(-3,1,1,l) 
            phiint(-4,1,1,l) = phiint(-3,1,1,l) 
            Bxint (-4,1,1,l) = Bxint (-3,1,1,l)  
            Byint (-4,1,1,l) = Byint (-3,1,1,l)  
            Bzint (-4,1,1,l) = Bzint (-3,1,1,l)  
            qint  (-4,1,1,l) = qint  (-3,1,1,l)   
            DDint (-4,1,1,l) = DDint (-3,1,1,l)  
            tauint(-4,1,1,l) = tauint(-3,1,1,l) 
            Sxint (-4,1,1,l) = Sxint (-3,1,1,l)  
            Syint (-4,1,1,l) = Syint (-3,1,1,l)  
            Szint (-4,1,1,l) = Szint (-3,1,1,l)

            
            psiint(-5,1,1,l) = psiint(-4,1,1,l) 
            phiint(-5,1,1,l) = phiint(-4,1,1,l) 
            Bxint (-5,1,1,l) = Bxint (-4,1,1,l)  
            Byint (-5,1,1,l) = Byint (-4,1,1,l)  
            Bzint (-5,1,1,l) = Bzint (-4,1,1,l)  
            qint  (-5,1,1,l) = qint  (-4,1,1,l)   
            DDint (-5,1,1,l) = DDint (-4,1,1,l)  
            tauint(-5,1,1,l) = tauint(-4,1,1,l) 
            Sxint (-5,1,1,l) = Sxint (-4,1,1,l)  
            Syint (-5,1,1,l) = Syint (-4,1,1,l)  
            Szint (-5,1,1,l) = Szint (-4,1,1,l)

            psiint(-6,1,1,l) = psiint(-5,1,1,l) 
            phiint(-6,1,1,l) = phiint(-5,1,1,l) 
            Bxint (-6,1,1,l) = Bxint (-5,1,1,l)  
            Byint (-6,1,1,l) = Byint (-5,1,1,l)  
            Bzint (-6,1,1,l) = Bzint (-5,1,1,l)  
            qint  (-6,1,1,l) = qint  (-5,1,1,l)   
            DDint (-6,1,1,l) = DDint (-5,1,1,l)  
            tauint(-6,1,1,l) = tauint(-5,1,1,l) 
            Sxint (-6,1,1,l) = Sxint (-5,1,1,l)  
            Syint (-6,1,1,l) = Syint (-5,1,1,l)  
            Szint (-6,1,1,l) = Szint (-5,1,1,l)  



            ! Right Boundary
            !---------------

!!$            psiint(imax  ,1,1,l) = psiint(imax-1,1,1,l) 
!!$            phiint(imax  ,1,1,l) = phiint(imax-1,1,1,l) 
!!$            Bxint (imax  ,1,1,l) = Bxint (imax-1,1,1,l)  
!!$            Byint (imax  ,1,1,l) = Byint (imax-1,1,1,l)  
!!$            Bzint (imax  ,1,1,l) = Bzint (imax-1,1,1,l)  
!!$            qint  (imax  ,1,1,l) = qint  (imax-1,1,1,l)   
!!$            DDint (imax  ,1,1,l) = DDint (imax-1,1,1,l)  
!!$            tauint(imax  ,1,1,l) = tauint(imax-1,1,1,l) 
!!$            Sxint (imax  ,1,1,l) = Sxint (imax-1,1,1,l)  
!!$            Syint (imax  ,1,1,l) = Syint (imax-1,1,1,l)  
!!$            Szint (imax  ,1,1,l) = Szint (imax-1,1,1,l) 


            psiint(imax+1,1,1,l) = psiint(imax  ,1,1,l) 
            phiint(imax+1,1,1,l) = phiint(imax  ,1,1,l) 
            Bxint (imax+1,1,1,l) = Bxint (imax  ,1,1,l)  
            Byint (imax+1,1,1,l) = Byint (imax  ,1,1,l)  
            Bzint (imax+1,1,1,l) = Bzint (imax  ,1,1,l)  
            qint  (imax+1,1,1,l) = qint  (imax  ,1,1,l)   
            DDint (imax+1,1,1,l) = DDint (imax  ,1,1,l)  
            tauint(imax+1,1,1,l) = tauint(imax  ,1,1,l) 
            Sxint (imax+1,1,1,l) = Sxint (imax  ,1,1,l)  
            Syint (imax+1,1,1,l) = Syint (imax  ,1,1,l)  
            Szint (imax+1,1,1,l) = Szint (imax  ,1,1,l) 


            psiint(imax+2,1,1,l) = psiint(imax+1,1,1,l) 
            phiint(imax+2,1,1,l) = phiint(imax+1,1,1,l) 
            Bxint (imax+2,1,1,l) = Bxint (imax+1,1,1,l)  
            Byint (imax+2,1,1,l) = Byint (imax+1,1,1,l)  
            Bzint (imax+2,1,1,l) = Bzint (imax+1,1,1,l)  
            qint  (imax+2,1,1,l) = qint  (imax+1,1,1,l)   
            DDint (imax+2,1,1,l) = DDint (imax+1,1,1,l)  
            tauint(imax+2,1,1,l) = tauint(imax+1,1,1,l) 
            Sxint (imax+2,1,1,l) = Sxint (imax+1,1,1,l)  
            Syint (imax+2,1,1,l) = Syint (imax+1,1,1,l)  
            Szint (imax+2,1,1,l) = Szint (imax+1,1,1,l) 


            psiint(imax+3,1,1,l) = psiint(imax+2,1,1,l) 
            phiint(imax+3,1,1,l) = phiint(imax+2,1,1,l) 
            Bxint (imax+3,1,1,l) = Bxint (imax+2,1,1,l)  
            Byint (imax+3,1,1,l) = Byint (imax+2,1,1,l)  
            Bzint (imax+3,1,1,l) = Bzint (imax+2,1,1,l)  
            qint  (imax+3,1,1,l) = qint  (imax+2,1,1,l)   
            DDint (imax+3,1,1,l) = DDint (imax+2,1,1,l)  
            tauint(imax+3,1,1,l) = tauint(imax+2,1,1,l) 
            Sxint (imax+3,1,1,l) = Sxint (imax+2,1,1,l)  
            Syint (imax+3,1,1,l) = Syint (imax+2,1,1,l)  
            Szint (imax+3,1,1,l) = Szint (imax+2,1,1,l) 


            psiint(imax+4,1,1,l) = psiint(imax+3,1,1,l) 
            phiint(imax+4,1,1,l) = phiint(imax+3,1,1,l) 
            Bxint (imax+4,1,1,l) = Bxint (imax+3,1,1,l)  
            Byint (imax+4,1,1,l) = Byint (imax+3,1,1,l)  
            Bzint (imax+4,1,1,l) = Bzint (imax+3,1,1,l)  
            qint  (imax+4,1,1,l) = qint  (imax+3,1,1,l)   
            DDint (imax+4,1,1,l) = DDint (imax+3,1,1,l)  
            tauint(imax+4,1,1,l) = tauint(imax+3,1,1,l) 
            Sxint (imax+4,1,1,l) = Sxint (imax+3,1,1,l)  
            Syint (imax+4,1,1,l) = Syint (imax+3,1,1,l)  
            Szint (imax+4,1,1,l) = Szint (imax+3,1,1,l)


            psiint(imax+5,1,1,l) = psiint(imax+4,1,1,l) 
            phiint(imax+5,1,1,l) = phiint(imax+4,1,1,l) 
            Bxint (imax+5,1,1,l) = Bxint (imax+4,1,1,l)  
            Byint (imax+5,1,1,l) = Byint (imax+4,1,1,l)  
            Bzint (imax+5,1,1,l) = Bzint (imax+4,1,1,l)  
            qint  (imax+5,1,1,l) = qint  (imax+4,1,1,l)   
            DDint (imax+5,1,1,l) = DDint (imax+4,1,1,l)  
            tauint(imax+5,1,1,l) = tauint(imax+4,1,1,l) 
            Sxint (imax+5,1,1,l) = Sxint (imax+4,1,1,l)  
            Syint (imax+5,1,1,l) = Syint (imax+4,1,1,l)  
            Szint (imax+5,1,1,l) = Szint (imax+4,1,1,l)

            
            psiint(imax+6,1,1,l) = psiint(imax+5,1,1,l) 
            phiint(imax+6,1,1,l) = phiint(imax+5,1,1,l) 
            Bxint (imax+6,1,1,l) = Bxint (imax+5,1,1,l)  
            Byint (imax+6,1,1,l) = Byint (imax+5,1,1,l)  
            Bzint (imax+6,1,1,l) = Bzint (imax+5,1,1,l)  
            qint  (imax+6,1,1,l) = qint  (imax+5,1,1,l)   
            DDint (imax+6,1,1,l) = DDint (imax+5,1,1,l)  
            tauint(imax+6,1,1,l) = tauint(imax+5,1,1,l) 
            Sxint (imax+6,1,1,l) = Sxint (imax+5,1,1,l)  
            Syint (imax+6,1,1,l) = Syint (imax+5,1,1,l)  
            Szint (imax+6,1,1,l) = Szint (imax+5,1,1,l) 



    else if (DIM == 1 .and. BOUND == 2) then

       ! Periodic boundary 1D

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

!!$            psiint( 0,1,1,l) = psiint(imax-1,1,1,l) 
!!$            phiint( 0,1,1,l) = phiint(imax-1,1,1,l) 
!!$            Bxint ( 0,1,1,l) = Bxint (imax-1,1,1,l)  
!!$            Byint ( 0,1,1,l) = Byint (imax-1,1,1,l)  
!!$            Bzint ( 0,1,1,l) = Bzint (imax-1,1,1,l)  
!!$            qint  ( 0,1,1,l) = qint  (imax-1,1,1,l)   
!!$            DDint ( 0,1,1,l) = DDint (imax-1,1,1,l)  
!!$            tauint( 0,1,1,l) = tauint(imax-1,1,1,l) 
!!$            Sxint ( 0,1,1,l) = Sxint (imax-1,1,1,l)  
!!$            Syint ( 0,1,1,l) = Syint (imax-1,1,1,l)  
!!$            Szint ( 0,1,1,l) = Szint (imax-1,1,1,l) 


            psiint(-1,1,1,l) = psiint(imax-1,1,1,l) 
            phiint(-1,1,1,l) = phiint(imax-1,1,1,l) 
            Bxint (-1,1,1,l) = Bxint (imax-1,1,1,l)  
            Byint (-1,1,1,l) = Byint (imax-1,1,1,l)  
            Bzint (-1,1,1,l) = Bzint (imax-1,1,1,l)  
            qint  (-1,1,1,l) = qint  (imax-1,1,1,l)   
            DDint (-1,1,1,l) = DDint (imax-1,1,1,l)  
            tauint(-1,1,1,l) = tauint(imax-1,1,1,l) 
            Sxint (-1,1,1,l) = Sxint (imax-1,1,1,l)  
            Syint (-1,1,1,l) = Syint (imax-1,1,1,l)  
            Szint (-1,1,1,l) = Szint (imax-1,1,1,l)

            Exast (-1,j,1,l)  = Exast (imax-1,1,1,l)
            Eyast (-1,j,1,l)  = Eyast (imax-1,1,1,l)
            Ezast (-1,j,1,l)  = Ezast (imax-1,1,1,l)
            

            psiint(-2,1,1,l) = psiint(imax-2,1,1,l) 
            phiint(-2,1,1,l) = phiint(imax-2,1,1,l) 
            Bxint (-2,1,1,l) = Bxint (imax-2,1,1,l)  
            Byint (-2,1,1,l) = Byint (imax-2,1,1,l)  
            Bzint (-2,1,1,l) = Bzint (imax-2,1,1,l)  
            qint  (-2,1,1,l) = qint  (imax-2,1,1,l)   
            DDint (-2,1,1,l) = DDint (imax-2,1,1,l)  
            tauint(-2,1,1,l) = tauint(imax-2,1,1,l) 
            Sxint (-2,1,1,l) = Sxint (imax-2,1,1,l)  
            Syint (-2,1,1,l) = Syint (imax-2,1,1,l)  
            Szint (-2,1,1,l) = Szint (imax-2,1,1,l)

            Exast (-2,j,1,l)  = Exast (imax-2,1,1,l)
            Eyast (-2,j,1,l)  = Eyast (imax-2,1,1,l)
            Ezast (-2,j,1,l)  = Ezast (imax-2,1,1,l)

            psiint(-3,1,1,l) = psiint(imax-3,1,1,l) 
            phiint(-3,1,1,l) = phiint(imax-3,1,1,l) 
            Bxint (-3,1,1,l) = Bxint (imax-3,1,1,l)  
            Byint (-3,1,1,l) = Byint (imax-3,1,1,l)  
            Bzint (-3,1,1,l) = Bzint (imax-3,1,1,l)  
            qint  (-3,1,1,l) = qint  (imax-3,1,1,l)   
            DDint (-3,1,1,l) = DDint (imax-3,1,1,l)  
            tauint(-3,1,1,l) = tauint(imax-3,1,1,l) 
            Sxint (-3,1,1,l) = Sxint (imax-3,1,1,l)  
            Syint (-3,1,1,l) = Syint (imax-3,1,1,l)  
            Szint (-3,1,1,l) = Szint (imax-3,1,1,l)

            Exast (-3,j,1,l)  = Exast (imax-3,1,1,l)
            Eyast (-3,j,1,l)  = Eyast (imax-3,1,1,l)
            Ezast (-3,j,1,l)  = Ezast (imax-3,1,1,l)


            psiint(-4,1,1,l) = psiint(imax-4,1,1,l) 
            phiint(-4,1,1,l) = phiint(imax-4,1,1,l) 
            Bxint (-4,1,1,l) = Bxint (imax-4,1,1,l)  
            Byint (-4,1,1,l) = Byint (imax-4,1,1,l)  
            Bzint (-4,1,1,l) = Bzint (imax-4,1,1,l)  
            qint  (-4,1,1,l) = qint  (imax-4,1,1,l)   
            DDint (-4,1,1,l) = DDint (imax-4,1,1,l)  
            tauint(-4,1,1,l) = tauint(imax-4,1,1,l) 
            Sxint (-4,1,1,l) = Sxint (imax-4,1,1,l)  
            Syint (-4,1,1,l) = Syint (imax-4,1,1,l)  
            Szint (-4,1,1,l) = Szint (imax-4,1,1,l)

            Exast (-4,j,1,l)  = Exast (imax-4,1,1,l)
            Eyast (-4,j,1,l)  = Eyast (imax-4,1,1,l)
            Ezast (-4,j,1,l)  = Ezast (imax-4,1,1,l)

            
            psiint(-5,1,1,l) = psiint(imax-5,1,1,l) 
            phiint(-5,1,1,l) = phiint(imax-5,1,1,l) 
            Bxint (-5,1,1,l) = Bxint (imax-5,1,1,l)  
            Byint (-5,1,1,l) = Byint (imax-5,1,1,l)  
            Bzint (-5,1,1,l) = Bzint (imax-5,1,1,l)  
            qint  (-5,1,1,l) = qint  (imax-5,1,1,l)   
            DDint (-5,1,1,l) = DDint (imax-5,1,1,l)  
            tauint(-5,1,1,l) = tauint(imax-5,1,1,l) 
            Sxint (-5,1,1,l) = Sxint (imax-5,1,1,l)  
            Syint (-5,1,1,l) = Syint (imax-5,1,1,l)  
            Szint (-5,1,1,l) = Szint (imax-5,1,1,l)

            Exast (-5,j,1,l) = Exast (imax-5,1,1,l)
            Eyast (-5,j,1,l) = Eyast (imax-5,1,1,l)
            Ezast (-5,j,1,l) = Ezast (imax-5,1,1,l)

            
            psiint(-6,1,1,l) = psiint(imax-6,1,1,l) 
            phiint(-6,1,1,l) = phiint(imax-6,1,1,l) 
            Bxint (-6,1,1,l) = Bxint (imax-6,1,1,l)  
            Byint (-6,1,1,l) = Byint (imax-6,1,1,l)  
            Bzint (-6,1,1,l) = Bzint (imax-6,1,1,l)  
            qint  (-6,1,1,l) = qint  (imax-6,1,1,l)   
            DDint (-6,1,1,l) = DDint (imax-6,1,1,l)  
            tauint(-6,1,1,l) = tauint(imax-6,1,1,l) 
            Sxint (-6,1,1,l) = Sxint (imax-6,1,1,l)  
            Syint (-6,1,1,l) = Syint (imax-6,1,1,l)  
            Szint (-6,1,1,l) = Szint (imax-6,1,1,l)

            
            Exast (-6,j,1,l) = Exast (imax-6,1,1,l)
            Eyast (-6,j,1,l) = Eyast (imax-6,1,1,l)
            Ezast (-6,j,1,l) = Ezast (imax-6,1,1,l)

            
            ! Right Boundary
            !---------------

!!$            psiint(imax  ,1,1,l) = psiint(1,1,1,l) 
!!$            phiint(imax  ,1,1,l) = phiint(1,1,1,l) 
!!$            Bxint (imax  ,1,1,l) = Bxint (1,1,1,l)  
!!$            Byint (imax  ,1,1,l) = Byint (1,1,1,l)  
!!$            Bzint (imax  ,1,1,l) = Bzint (1,1,1,l)  
!!$            qint  (imax  ,1,1,l) = qint  (1,1,1,l)   
!!$            DDint (imax  ,1,1,l) = DDint (1,1,1,l)  
!!$            tauint(imax  ,1,1,l) = tauint(1,1,1,l) 
!!$            Sxint (imax  ,1,1,l) = Sxint (1,1,1,l)  
!!$            Syint (imax  ,1,1,l) = Syint (1,1,1,l)  
!!$            Szint (imax  ,1,1,l) = Szint (1,1,1,l) 


            psiint(imax+1,1,1,l) = psiint(1,1,1,l) 
            phiint(imax+1,1,1,l) = phiint(1,1,1,l) 
            Bxint (imax+1,1,1,l) = Bxint (1,1,1,l)  
            Byint (imax+1,1,1,l) = Byint (1,1,1,l)  
            Bzint (imax+1,1,1,l) = Bzint (1,1,1,l)  
            qint  (imax+1,1,1,l) = qint  (1,1,1,l)   
            DDint (imax+1,1,1,l) = DDint (1,1,1,l)  
            tauint(imax+1,1,1,l) = tauint(1,1,1,l) 
            Sxint (imax+1,1,1,l) = Sxint (1,1,1,l)  
            Syint (imax+1,1,1,l) = Syint (1,1,1,l)  
            Szint (imax+1,1,1,l) = Szint (1,1,1,l)

            Exast (imax+1,j,1,l) = Exast (1,1,1,l)
            Eyast (imax+1,j,1,l) = Eyast (1,1,1,l)
            Ezast (imax+1,j,1,l) = Ezast (1,1,1,l)

            psiint(imax+2,1,1,l) = psiint(2,1,1,l) 
            phiint(imax+2,1,1,l) = phiint(2,1,1,l) 
            Bxint (imax+2,1,1,l) = Bxint (2,1,1,l)  
            Byint (imax+2,1,1,l) = Byint (2,1,1,l)  
            Bzint (imax+2,1,1,l) = Bzint (2,1,1,l)  
            qint  (imax+2,1,1,l) = qint  (2,1,1,l)   
            DDint (imax+2,1,1,l) = DDint (2,1,1,l)  
            tauint(imax+2,1,1,l) = tauint(2,1,1,l) 
            Sxint (imax+2,1,1,l) = Sxint (2,1,1,l)  
            Syint (imax+2,1,1,l) = Syint (2,1,1,l)  
            Szint (imax+2,1,1,l) = Szint (2,1,1,l)

            Exast (imax+2,j,1,l) = Exast (2,1,1,l)
            Eyast (imax+2,j,1,l) = Eyast (2,1,1,l)
            Ezast (imax+2,j,1,l) = Ezast (2,1,1,l)


            psiint(imax+3,1,1,l) = psiint(3,1,1,l) 
            phiint(imax+3,1,1,l) = phiint(3,1,1,l) 
            Bxint (imax+3,1,1,l) = Bxint (3,1,1,l)  
            Byint (imax+3,1,1,l) = Byint (3,1,1,l)  
            Bzint (imax+3,1,1,l) = Bzint (3,1,1,l)  
            qint  (imax+3,1,1,l) = qint  (3,1,1,l)   
            DDint (imax+3,1,1,l) = DDint (3,1,1,l)  
            tauint(imax+3,1,1,l) = tauint(3,1,1,l) 
            Sxint (imax+3,1,1,l) = Sxint (3,1,1,l)  
            Syint (imax+3,1,1,l) = Syint (3,1,1,l)  
            Szint (imax+3,1,1,l) = Szint (3,1,1,l)

            Exast (imax+3,j,1,l) = Exast (3,1,1,l)
            Eyast (imax+3,j,1,l) = Eyast (3,1,1,l)
            Ezast (imax+3,j,1,l) = Ezast (3,1,1,l)


            psiint(imax+4,1,1,l) = psiint(4,1,1,l) 
            phiint(imax+4,1,1,l) = phiint(4,1,1,l) 
            Bxint (imax+4,1,1,l) = Bxint (4,1,1,l)  
            Byint (imax+4,1,1,l) = Byint (4,1,1,l)  
            Bzint (imax+4,1,1,l) = Bzint (4,1,1,l)  
            qint  (imax+4,1,1,l) = qint  (4,1,1,l)   
            DDint (imax+4,1,1,l) = DDint (4,1,1,l)  
            tauint(imax+4,1,1,l) = tauint(4,1,1,l) 
            Sxint (imax+4,1,1,l) = Sxint (4,1,1,l)  
            Syint (imax+4,1,1,l) = Syint (4,1,1,l)  
            Szint (imax+4,1,1,l) = Szint (4,1,1,l)

            Exast (imax+4,j,1,l) = Exast (4,1,1,l)
            Eyast (imax+4,j,1,l) = Eyast (4,1,1,l)
            Ezast (imax+4,j,1,l) = Ezast (4,1,1,l)

            
            psiint(imax+5,1,1,l) = psiint(5,1,1,l) 
            phiint(imax+5,1,1,l) = phiint(5,1,1,l) 
            Bxint (imax+5,1,1,l) = Bxint (5,1,1,l)  
            Byint (imax+5,1,1,l) = Byint (5,1,1,l)  
            Bzint (imax+5,1,1,l) = Bzint (5,1,1,l)  
            qint  (imax+5,1,1,l) = qint  (5,1,1,l)   
            DDint (imax+5,1,1,l) = DDint (5,1,1,l)  
            tauint(imax+5,1,1,l) = tauint(5,1,1,l) 
            Sxint (imax+5,1,1,l) = Sxint (5,1,1,l)  
            Syint (imax+5,1,1,l) = Syint (5,1,1,l)  
            Szint (imax+5,1,1,l) = Szint (5,1,1,l)

            Exast (imax+5,j,1,l) = Exast (5,1,1,l)
            Eyast (imax+5,j,1,l) = Eyast (5,1,1,l)
            Ezast (imax+5,j,1,l) = Ezast (5,1,1,l)

            psiint(imax+6,1,1,l) = psiint(6,1,1,l) 
            phiint(imax+6,1,1,l) = phiint(6,1,1,l) 
            Bxint (imax+6,1,1,l) = Bxint (6,1,1,l)  
            Byint (imax+6,1,1,l) = Byint (6,1,1,l)  
            Bzint (imax+6,1,1,l) = Bzint (6,1,1,l)  
            qint  (imax+6,1,1,l) = qint  (6,1,1,l)   
            DDint (imax+6,1,1,l) = DDint (6,1,1,l)  
            tauint(imax+6,1,1,l) = tauint(6,1,1,l) 
            Sxint (imax+6,1,1,l) = Sxint (6,1,1,l)  
            Syint (imax+6,1,1,l) = Syint (6,1,1,l)  
            Szint (imax+6,1,1,l) = Szint (6,1,1,l)

            Exast (imax+6,j,1,l) = Exast (6,1,1,l)
            Eyast (imax+6,j,1,l) = Eyast (6,1,1,l)
            Ezast (imax+6,j,1,l) = Ezast (6,1,1,l)

    else if (DIM == 1 .and. BOUND == 3) then


       ! Copy boundary 1D

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

            psiint(0,1,1,l) = psi(0,1,1) !psiint(1,1,1,l) 
            phiint(0,1,1,l) = phi(0,1,1) !phiint(1,1,1,l) 
            Bxint(0,1,1,l)  = Bx (0,1,1) !Bxint(1,1,1,l)  
            Byint(0,1,1,l)  = By (0,1,1) !Byint(1,1,1,l)  
            Bzint(0,1,1,l)  = Bz (0,1,1) !Bzint(1,1,1,l)  
            qint(0,1,1,l)   = q  (0,1,1) !qint(1,1,1,l)   
            DDint(0,1,1,l)  = D  (0,1,1) !DDint(1,1,1,l)  
            tauint(0,1,1,l) = tau(0,1,1) !tauint(1,1,1,l) 
            Sxint(0,1,1,l)  = Sx (0,1,1) !Sxint(1,1,1,l)  
            Syint(0,1,1,l)  = Sy (0,1,1) !Syint(1,1,1,l)  
            Szint(0,1,1,l)  = Sz (0,1,1) !Szint(1,1,1,l)  


            psiint(-1,1,1,l) = psiint(0,1,1,l) 
            phiint(-1,1,1,l) = phiint(0,1,1,l) 
            Bxint(-1,1,1,l)  = Bxint(0,1,1,l)  
            Byint(-1,1,1,l)  = Byint(0,1,1,l)  
            Bzint(-1,1,1,l)  = Bzint(0,1,1,l)  
            qint(-1,1,1,l)   = qint(0,1,1,l)   
            DDint(-1,1,1,l)  = DDint(0,1,1,l)  
            tauint(-1,1,1,l) = tauint(0,1,1,l) 
            Sxint(-1,1,1,l)  = Sxint(0,1,1,l)  
            Syint(-1,1,1,l)  = Syint(0,1,1,l)  
            Szint(-1,1,1,l)  = Szint(0,1,1,l)  

            psiint(-2,1,1,l) = psiint(-1,1,1,l) 
            phiint(-2,1,1,l) = phiint(-1,1,1,l) 
            Bxint(-2,1,1,l)  = Bxint(-1,1,1,l)  
            Byint(-2,1,1,l)  = Byint(-1,1,1,l)  
            Bzint(-2,1,1,l)  = Bzint(-1,1,1,l)  
            qint(-2,1,1,l)   = qint(-1,1,1,l)   
            DDint(-2,1,1,l)  = DDint(-1,1,1,l)  
            tauint(-2,1,1,l) = tauint(-1,1,1,l) 
            Sxint(-2,1,1,l)  = Sxint(-1,1,1,l)  
            Syint(-2,1,1,l)  = Syint(-1,1,1,l)  
            Szint(-2,1,1,l)  = Szint(-1,1,1,l)  

            ! Right Boundary
            !---------------

            psiint(imax,1,1,l) = psi(imax,1,1) !psiint(imax-1,1,1,l) 
            phiint(imax,1,1,l) = phi(imax,1,1) !phiint(imax-1,1,1,l) 
            Bxint(imax,1,1,l)  = Bx (imax,1,1) !Bxint(imax-1,1,1,l)  
            Byint(imax,1,1,l)  = By (imax,1,1) !Byint(imax-1,1,1,l)  
            Bzint(imax,1,1,l)  = Bz (imax,1,1) !Bzint(imax-1,1,1,l)  
            qint(imax,1,1,l)   = q  (imax,1,1) !qint(imax-1,1,1,l)   
            DDint(imax,1,1,l)  = D  (imax,1,1) !DDint(imax-1,1,1,l)  
            tauint(imax,1,1,l) = tau(imax,1,1) !tauint(imax-1,1,1,l) 
            Sxint(imax,1,1,l)  = Sx (imax,1,1) !Sxint(imax-1,1,1,l)  
            Syint(imax,1,1,l)  = Sy (imax,1,1) !Syint(imax-1,1,1,l)  
            Szint(imax,1,1,l)  = Sz (imax,1,1) !Szint(imax-1,1,1,l) 

            psiint(imax+1,1,1,l) = psiint(imax,1,1,l) 
            phiint(imax+1,1,1,l) = phiint(imax,1,1,l) 
            Bxint(imax+1,1,1,l)  = Bxint(imax,1,1,l)  
            Byint(imax+1,1,1,l)  = Byint(imax,1,1,l)  
            Bzint(imax+1,1,1,l)  = Bzint(imax,1,1,l)  
            qint(imax+1,1,1,l)   = qint(imax,1,1,l)   
            DDint(imax+1,1,1,l)  = DDint(imax,1,1,l)  
            tauint(imax+1,1,1,l) = tauint(imax,1,1,l) 
            Sxint(imax+1,1,1,l)  = Sxint(imax,1,1,l)  
            Syint(imax+1,1,1,l)  = Syint(imax,1,1,l)  
            Szint(imax+1,1,1,l)  = Szint(imax,1,1,l)  

            psiint(imax+2,1,1,l) = psiint(imax,1,1,l) 
            phiint(imax+2,1,1,l) = phiint(imax,1,1,l) 
            Bxint(imax+2,1,1,l)  = Bxint(imax,1,1,l)  
            Byint(imax+2,1,1,l)  = Byint(imax,1,1,l)  
            Bzint(imax+2,1,1,l)  = Bzint(imax,1,1,l)  
            qint(imax+2,1,1,l)   = qint(imax,1,1,l)   
            DDint(imax+2,1,1,l)  = DDint(imax,1,1,l)  
            tauint(imax+2,1,1,l) = tauint(imax,1,1,l) 
            Sxint(imax+2,1,1,l)  = Sxint(imax,1,1,l)  
            Syint(imax+2,1,1,l)  = Syint(imax,1,1,l)  
            Szint(imax+2,1,1,l)  = Szint(imax,1,1,l) 


    else if (DIM == 1 .and. BOUND == 4) then


      ! Copy boundary 1D + perfect conductor

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

            psiint(0,1,1,l) = psiint(1,1,1,l) 
            phiint(0,1,1,l) = phiint(1,1,1,l) 
            Bxint(0,1,1,l)  = 0.d0 !Bxint(1,1,1,l)  
            Byint(0,1,1,l)  = Byint(1,1,1,l)  
            Bzint(0,1,1,l)  = Bzint(1,1,1,l)  
            qint(0,1,1,l)   = qint(1,1,1,l)   
            DDint(0,1,1,l)  = DDint(1,1,1,l)  
            tauint(0,1,1,l) = tauint(1,1,1,l) 
            Sxint(0,1,1,l)  = Sxint(1,1,1,l)  
            Syint(0,1,1,l)  = Syint(1,1,1,l)  
            Szint(0,1,1,l)  = Szint(1,1,1,l)  


            psiint(-1,1,1,l) = psiint(0,1,1,l) 
            phiint(-1,1,1,l) = phiint(0,1,1,l) 
            Bxint(-1,1,1,l)  = Bxint(0,1,1,l)  
            Byint(-1,1,1,l)  = Byint(0,1,1,l)  
            Bzint(-1,1,1,l)  = Bzint(0,1,1,l)  
            qint(-1,1,1,l)   = qint(0,1,1,l)   
            DDint(-1,1,1,l)  = DDint(0,1,1,l)  
            tauint(-1,1,1,l) = tauint(0,1,1,l) 
            Sxint(-1,1,1,l)  = Sxint(0,1,1,l)  
            Syint(-1,1,1,l)  = Syint(0,1,1,l)  
            Szint(-1,1,1,l)  = Szint(0,1,1,l)  

            psiint(-2,1,1,l) = psiint(-1,1,1,l) 
            phiint(-2,1,1,l) = phiint(-1,1,1,l) 
            Bxint(-2,1,1,l)  = Bxint(-1,1,1,l)  
            Byint(-2,1,1,l)  = Byint(-1,1,1,l)  
            Bzint(-2,1,1,l)  = Bzint(-1,1,1,l)  
            qint(-2,1,1,l)   = qint(-1,1,1,l)   
            DDint(-2,1,1,l)  = DDint(-1,1,1,l)  
            tauint(-2,1,1,l) = tauint(-1,1,1,l) 
            Sxint(-2,1,1,l)  = Sxint(-1,1,1,l)  
            Syint(-2,1,1,l)  = Syint(-1,1,1,l)  
            Szint(-2,1,1,l)  = Szint(-1,1,1,l)  

            ! Right Boundary
            !---------------

            psiint(imax,1,1,l) = psiint(imax-1,1,1,l) 
            phiint(imax,1,1,l) = phiint(imax-1,1,1,l) 
            Bxint(imax,1,1,l)  = 0.d0 !Bxint(imax-1,1,1,l)  
            Byint(imax,1,1,l)  = Byint(imax-1,1,1,l)  
            Bzint(imax,1,1,l)  = Bzint(imax-1,1,1,l)  
            qint(imax,1,1,l)   = qint(imax-1,1,1,l)   
            DDint(imax,1,1,l)  = DDint(imax-1,1,1,l)  
            tauint(imax,1,1,l) = tauint(imax-1,1,1,l) 
            Sxint(imax,1,1,l)  = Sxint(imax-1,1,1,l)  
            Syint(imax,1,1,l)  = Syint(imax-1,1,1,l)  
            Szint(imax,1,1,l)  = Szint(imax-1,1,1,l) 

            psiint(imax+1,1,1,l) = psiint(imax,1,1,l) 
            phiint(imax+1,1,1,l) = phiint(imax,1,1,l) 
            Bxint(imax+1,1,1,l)  = Bxint(imax,1,1,l)  
            Byint(imax+1,1,1,l)  = Byint(imax,1,1,l)  
            Bzint(imax+1,1,1,l)  = Bzint(imax,1,1,l)  
            qint(imax+1,1,1,l)   = qint(imax,1,1,l)   
            DDint(imax+1,1,1,l)  = DDint(imax,1,1,l)  
            tauint(imax+1,1,1,l) = tauint(imax,1,1,l) 
            Sxint(imax+1,1,1,l)  = Sxint(imax,1,1,l)  
            Syint(imax+1,1,1,l)  = Syint(imax,1,1,l)  
            Szint(imax+1,1,1,l)  = Szint(imax,1,1,l)  

            psiint(imax+2,1,1,l) = psiint(imax,1,1,l) 
            phiint(imax+2,1,1,l) = phiint(imax,1,1,l) 
            Bxint(imax+2,1,1,l)  = Bxint(imax,1,1,l)  
            Byint(imax+2,1,1,l)  = Byint(imax,1,1,l)  
            Bzint(imax+2,1,1,l)  = Bzint(imax,1,1,l)  
            qint(imax+2,1,1,l)   = qint(imax,1,1,l)   
            DDint(imax+2,1,1,l)  = DDint(imax,1,1,l)  
            tauint(imax+2,1,1,l) = tauint(imax,1,1,l) 
            Sxint(imax+2,1,1,l)  = Sxint(imax,1,1,l)  
            Syint(imax+2,1,1,l)  = Syint(imax,1,1,l)  
            Szint(imax+2,1,1,l)  = Szint(imax,1,1,l) 


       ! Perfect Conductor and anti-simetric reflection TM_III

  !-----------------
  ! X Direction
  !-----------------
            ! Left Boundary
            !---------------

!!$            psiint(0,1,1,l) = psiint(imax-1,1,1,l) 
!!$            phiint(0,1,1,l) = phiint(imax-1,1,1,l) 
!!$            Bxint(0,1,1,l)  = 0.d0
!!$            Byint(0,1,1,l)  = Byint(imax-1,1,1,l)  
!!$            Bzint(0,1,1,l)  = Bzint(imax-1,1,1,l)  
!!$            qint(0,1,1,l)   = qint(imax-1,1,1,l)   
!!$            DDint(0,1,1,l)  = DDint(imax-1,1,1,l)  
!!$            tauint(0,1,1,l) = tauint(imax-1,1,1,l) 
!!$            Sxint(0,1,1,l)  = Sxint(     1,1,1,l)  
!!$            Syint(0,1,1,l)  = Syint(imax-1,1,1,l)  
!!$            Szint(0,1,1,l)  = Szint(imax-1,1,1,l)  


!!$            psiint(-1,1,1,l) = 2.d0 * psiint(0,1,1,l) - psiint(1,1,1,l) 
!!$            phiint(-1,1,1,l) = 2.d0 * phiint(0,1,1,l) - phiint(1,1,1,l) 
!!$            Bxint(-1,1,1,l)  = - Bxint(1,1,1,l)
!!$            Byint(-1,1,1,l)  = 2.d0 * Byint (0,1,1,l) - Byint (1,1,1,l) 
!!$            Bzint(-1,1,1,l)  = 2.d0 * Bzint (0,1,1,l) - Bzint (1,1,1,l) 
!!$            qint(-1,1,1,l)   = 2.d0 * qint  (0,1,1,l) - qint  (1,1,1,l)  
!!$            DDint(-1,1,1,l)  = 2.d0 * DDint (0,1,1,l) - DDint (1,1,1,l)
!!$            tauint(-1,1,1,l) = 2.d0 * tauint(0,1,1,l) - tauint(1,1,1,l)
!!$            Sxint(-1,1,1,l)  = 2.d0 * Sxint (0,1,1,l) - Sxint (1,1,1,l)
!!$            Syint(-1,1,1,l)  = 2.d0 * Syint (0,1,1,l) - Syint (1,1,1,l) 
!!$            Szint(-1,1,1,l)  = 2.d0 * Szint (0,1,1,l) - Szint (1,1,1,l)  
!!$
!!$            psiint(-2,1,1,l) = 2.d0 * psiint(0,1,1,l) - psiint(2,1,1,l) 
!!$            phiint(-2,1,1,l) = 2.d0 * phiint(0,1,1,l) - phiint(2,1,1,l) 
!!$            Bxint(-2,1,1,l)  = - Bxint(2,1,1,l)
!!$            Byint(-2,1,1,l)  = 2.d0 * Byint (0,1,1,l) - Byint (2,1,1,l) 
!!$            Bzint(-2,1,1,l)  = 2.d0 * Bzint (0,1,1,l) - Bzint (2,1,1,l) 
!!$            qint(-2,1,1,l)   = 2.d0 * qint  (0,1,1,l) - qint  (2,1,1,l)  
!!$            DDint(-2,1,1,l)  = 2.d0 * DDint (0,1,1,l) - DDint (2,1,1,l)
!!$            tauint(-2,1,1,l) = 2.d0 * tauint(0,1,1,l) - tauint(2,1,1,l)
!!$            Sxint(-2,1,1,l)  = 2.d0 * Sxint (0,1,1,l) - Sxint (2,1,1,l)
!!$            Syint(-2,1,1,l)  = 2.d0 * Syint (0,1,1,l) - Syint (2,1,1,l) 
!!$            Szint(-2,1,1,l)  = 2.d0 * Szint (0,1,1,l) - Szint (2,1,1,l) 

            ! Right Boundary
            !---------------

!!$            psiint(imax,1,1,l) = psiint(1,1,1,l) 
!!$            phiint(imax,1,1,l) = phiint(1,1,1,l) 
!!$            Bxint(imax,1,1,l)  = 0.d0
!!$            Byint(imax,1,1,l)  = Byint(1,1,1,l)  
!!$            Bzint(imax,1,1,l)  = Bzint(1,1,1,l)  
!!$            qint(imax,1,1,l)   = qint(1,1,1,l)   
!!$            DDint(imax,1,1,l)  = DDint(1,1,1,l)  
!!$            tauint(imax,1,1,l) = tauint(1,1,1,l) 
!!$            Sxint(imax,1,1,l)  = Sxint(imax-1,1,1,l)  
!!$            Syint(imax,1,1,l)  = Syint(1,1,1,l)  
!!$            Szint(imax,1,1,l)  = Szint(1,1,1,l) 


!!$            psiint(imax+1,1,1,l) = 2.d0 * psiint(imax,1,1,l) - psiint(imax-1,1,1,l) 
!!$            phiint(imax+1,1,1,l) = 2.d0 * phiint(imax,1,1,l) - phiint(imax-1,1,1,l) 
!!$            Bxint(imax+1,1,1,l)  = - Bxint(imax-1,1,1,l)
!!$            Byint(imax+1,1,1,l)  = 2.d0 * Byint (imax,1,1,l) - Byint (imax-1,1,1,l) 
!!$            Bzint(imax+1,1,1,l)  = 2.d0 * Bzint (imax,1,1,l) - Bzint (imax-1,1,1,l) 
!!$            qint(imax+1,1,1,l)   = 2.d0 * qint  (imax,1,1,l) - qint  (imax-1,1,1,l)  
!!$            DDint(imax+1,1,1,l)  = 2.d0 * DDint (imax,1,1,l) - DDint (imax-1,1,1,l)
!!$            tauint(imax+1,1,1,l) = 2.d0 * tauint(imax,1,1,l) - tauint(imax-1,1,1,l)
!!$            Sxint(imax+1,1,1,l)  = 2.d0 * Sxint (imax,1,1,l) - Sxint (imax-1,1,1,l)
!!$            Syint(imax+1,1,1,l)  = 2.d0 * Syint (imax,1,1,l) - Syint (imax-1,1,1,l) 
!!$            Szint(imax+1,1,1,l)  = 2.d0 * Szint (imax,1,1,l) - Szint (imax-1,1,1,l)  
!!$
!!$            psiint(imax+2,1,1,l) = 2.d0 * psiint(imax,1,1,l) - psiint(imax-2,1,1,l) 
!!$            phiint(imax+2,1,1,l) = 2.d0 * phiint(imax,1,1,l) - phiint(imax-2,1,1,l) 
!!$            Bxint(imax+2,1,1,l)  = - Bxint(imax-2,1,1,l)
!!$            Byint(imax+2,1,1,l)  = 2.d0 * Byint (imax,1,1,l) - Byint (imax-2,1,1,l) 
!!$            Bzint(imax+2,1,1,l)  = 2.d0 * Bzint (imax,1,1,l) - Bzint (imax-2,1,1,l) 
!!$            qint(imax+2,1,1,l)   = 2.d0 * qint  (imax,1,1,l) - qint  (imax-2,1,1,l)  
!!$            DDint(imax+2,1,1,l)  = 2.d0 * DDint (imax,1,1,l) - DDint (imax-2,1,1,l)
!!$            tauint(imax+2,1,1,l) = 2.d0 * tauint(imax,1,1,l) - tauint(imax-2,1,1,l)
!!$            Sxint(imax+2,1,1,l)  = 2.d0 * Sxint (imax,1,1,l) - Sxint (imax-2,1,1,l)
!!$            Syint(imax+2,1,1,l)  = 2.d0 * Syint (imax,1,1,l) - Syint (imax-2,1,1,l) 
!!$            Szint(imax+2,1,1,l)  = 2.d0 * Szint (imax,1,1,l) - Szint (imax-2,1,1,l) 


    else if (DIM == 2 .and. BOUND == 1) then

       ! Copy boundary 2D

  !-----------------
  ! X Direction
  !-----------------

 
!$OMP DO
       
     do j=-6,jmax+6
!      do k=1,kmax

 

            ! Left Boundary
            !---------------

!!$            psiint( 0,j,1,l) = psiint( 1,j,1,l)
!!$            phiint( 0,j,1,l) = phiint( 1,j,1,l)
!!$            Bxint ( 0,j,1,l) = Bxint ( 1,j,1,l) 
!!$            Byint ( 0,j,1,l) = Byint ( 1,j,1,l) 
!!$            Bzint ( 0,j,1,l) = Bzint ( 1,j,1,l) 
!!$            qint  ( 0,j,1,l) = qint  ( 1,j,1,l) 
!!$            DDint ( 0,j,1,l) = DDint ( 1,j,1,l)
!!$            tauint( 0,j,1,l) = tauint( 1,j,1,l)
!!$            Sxint ( 0,j,1,l) = Sxint ( 1,j,1,l) 
!!$            Syint ( 0,j,1,l) = Syint ( 1,j,1,l) 
!!$            Szint ( 0,j,1,l) = Szint ( 1,j,1,l) 


            psiint(-1,j,1,l) = psiint( 0,j,1,l)
            phiint(-1,j,1,l) = phiint( 0,j,1,l)
            Bxint (-1,j,1,l) = Bxint ( 0,j,1,l) 
            Byint (-1,j,1,l) = Byint ( 0,j,1,l) 
            Bzint (-1,j,1,l) = Bzint ( 0,j,1,l) 
            qint  (-1,j,1,l) = qint  ( 0,j,1,l) 
            DDint (-1,j,1,l) = DDint ( 0,j,1,l)
            tauint(-1,j,1,l) = tauint( 0,j,1,l)
            Sxint (-1,j,1,l) = Sxint ( 0,j,1,l) 
            Syint (-1,j,1,l) = Syint ( 0,j,1,l) 
            Szint (-1,j,1,l) = Szint ( 0,j,1,l) 


            psiint(-2,j,1,l) = psiint(-1,j,1,l)
            phiint(-2,j,1,l) = phiint(-1,j,1,l)
            Bxint (-2,j,1,l) = Bxint (-1,j,1,l) 
            Byint (-2,j,1,l) = Byint (-1,j,1,l) 
            Bzint (-2,j,1,l) = Bzint (-1,j,1,l) 
            qint  (-2,j,1,l) = qint  (-1,j,1,l) 
            DDint (-2,j,1,l) = DDint (-1,j,1,l)
            tauint(-2,j,1,l) = tauint(-1,j,1,l)
            Sxint (-2,j,1,l) = Sxint (-1,j,1,l) 
            Syint (-2,j,1,l) = Syint (-1,j,1,l) 
            Szint (-2,j,1,l) = Szint (-1,j,1,l) 


            psiint(-3,j,1,l) = psiint(-2,j,1,l)
            phiint(-3,j,1,l) = phiint(-2,j,1,l)
            Bxint (-3,j,1,l) = Bxint (-2,j,1,l) 
            Byint (-3,j,1,l) = Byint (-2,j,1,l) 
            Bzint (-3,j,1,l) = Bzint (-2,j,1,l) 
            qint  (-3,j,1,l) = qint  (-2,j,1,l) 
            DDint (-3,j,1,l) = DDint (-2,j,1,l)
            tauint(-3,j,1,l) = tauint(-2,j,1,l)
            Sxint (-3,j,1,l) = Sxint (-2,j,1,l) 
            Syint (-3,j,1,l) = Syint (-2,j,1,l) 
            Szint (-3,j,1,l) = Szint (-2,j,1,l) 


            psiint(-4,j,1,l) = psiint(-3,j,1,l)
            phiint(-4,j,1,l) = phiint(-3,j,1,l)
            Bxint (-4,j,1,l) = Bxint (-3,j,1,l) 
            Byint (-4,j,1,l) = Byint (-3,j,1,l) 
            Bzint (-4,j,1,l) = Bzint (-3,j,1,l) 
            qint  (-4,j,1,l) = qint  (-3,j,1,l) 
            DDint (-4,j,1,l) = DDint (-3,j,1,l)
            tauint(-4,j,1,l) = tauint(-3,j,1,l)
            Sxint (-4,j,1,l) = Sxint (-3,j,1,l) 
            Syint (-4,j,1,l) = Syint (-3,j,1,l) 
            Szint (-4,j,1,l) = Szint (-3,j,1,l)

            psiint(-5,j,1,l) = psiint(-4,j,1,l)
            phiint(-5,j,1,l) = phiint(-4,j,1,l)
            Bxint (-5,j,1,l) = Bxint (-4,j,1,l) 
            Byint (-5,j,1,l) = Byint (-4,j,1,l) 
            Bzint (-5,j,1,l) = Bzint (-4,j,1,l) 
            qint  (-5,j,1,l) = qint  (-4,j,1,l) 
            DDint (-5,j,1,l) = DDint (-4,j,1,l)
            tauint(-5,j,1,l) = tauint(-4,j,1,l)
            Sxint (-5,j,1,l) = Sxint (-4,j,1,l) 
            Syint (-5,j,1,l) = Syint (-4,j,1,l) 
            Szint (-5,j,1,l) = Szint (-4,j,1,l)

            psiint(-6,j,1,l) = psiint(-5,j,1,l)
            phiint(-6,j,1,l) = phiint(-5,j,1,l)
            Bxint (-6,j,1,l) = Bxint (-5,j,1,l) 
            Byint (-6,j,1,l) = Byint (-5,j,1,l) 
            Bzint (-6,j,1,l) = Bzint (-5,j,1,l) 
            qint  (-6,j,1,l) = qint  (-5,j,1,l) 
            DDint (-6,j,1,l) = DDint (-5,j,1,l)
            tauint(-6,j,1,l) = tauint(-5,j,1,l)
            Sxint (-6,j,1,l) = Sxint (-5,j,1,l) 
            Syint (-6,j,1,l) = Syint (-5,j,1,l) 
            Szint (-6,j,1,l) = Szint (-5,j,1,l) 



            ! Right Boundary
            !---------------

!!$            psiint(imax  ,j,1,l) = psiint(imax-1,j,1,l) 
!!$            phiint(imax  ,j,1,l) = phiint(imax-1,j,1,l) 
!!$            Bxint (imax  ,j,1,l) = Bxint (imax-1,j,1,l)  
!!$            Byint (imax  ,j,1,l) = Byint (imax-1,j,1,l)  
!!$            Bzint (imax  ,j,1,l) = Bzint (imax-1,j,1,l)  
!!$            qint  (imax  ,j,1,l) = qint  (imax-1,j,1,l)  
!!$            DDint (imax  ,j,1,l) = DDint (imax-1,j,1,l) 
!!$            tauint(imax  ,j,1,l) = tauint(imax-1,j,1,l)
!!$            Sxint (imax  ,j,1,l) = Sxint (imax-1,j,1,l) 
!!$            Syint (imax  ,j,1,l) = Syint (imax-1,j,1,l) 
!!$            Szint (imax  ,j,1,l) = Szint (imax-1,j,1,l) 


            psiint(imax+1,j,1,l) = psiint(imax  ,j,1,l) 
            phiint(imax+1,j,1,l) = phiint(imax  ,j,1,l) 
            Bxint (imax+1,j,1,l) = Bxint (imax  ,j,1,l)  
            Byint (imax+1,j,1,l) = Byint (imax  ,j,1,l)  
            Bzint (imax+1,j,1,l) = Bzint (imax  ,j,1,l)  
            qint  (imax+1,j,1,l) = qint  (imax  ,j,1,l)  
            DDint (imax+1,j,1,l) = DDint (imax  ,j,1,l) 
            tauint(imax+1,j,1,l) = tauint(imax  ,j,1,l)
            Sxint (imax+1,j,1,l) = Sxint (imax  ,j,1,l) 
            Syint (imax+1,j,1,l) = Syint (imax  ,j,1,l) 
            Szint (imax+1,j,1,l) = Szint (imax  ,j,1,l) 


            psiint(imax+2,j,1,l) = psiint(imax+1,j,1,l)  
            phiint(imax+2,j,1,l) = phiint(imax+1,j,1,l)  
            Bxint (imax+2,j,1,l) = Bxint (imax+1,j,1,l)  
            Byint (imax+2,j,1,l) = Byint (imax+1,j,1,l)  
            Bzint (imax+2,j,1,l) = Bzint (imax+1,j,1,l)  
            qint  (imax+2,j,1,l) = qint  (imax+1,j,1,l)  
            DDint (imax+2,j,1,l) = DDint (imax+1,j,1,l) 
            tauint(imax+2,j,1,l) = tauint(imax+1,j,1,l)
            Sxint (imax+2,j,1,l) = Sxint (imax+1,j,1,l) 
            Syint (imax+2,j,1,l) = Syint (imax+1,j,1,l) 
            Szint (imax+2,j,1,l) = Szint (imax+1,j,1,l) 


            psiint(imax+3,j,1,l) = psiint(imax+2,j,1,l)  
            phiint(imax+3,j,1,l) = phiint(imax+2,j,1,l)  
            Bxint (imax+3,j,1,l) = Bxint (imax+2,j,1,l)  
            Byint (imax+3,j,1,l) = Byint (imax+2,j,1,l)  
            Bzint (imax+3,j,1,l) = Bzint (imax+2,j,1,l)  
            qint  (imax+3,j,1,l) = qint  (imax+2,j,1,l)  
            DDint (imax+3,j,1,l) = DDint (imax+2,j,1,l) 
            tauint(imax+3,j,1,l) = tauint(imax+2,j,1,l)
            Sxint (imax+3,j,1,l) = Sxint (imax+2,j,1,l) 
            Syint (imax+3,j,1,l) = Syint (imax+2,j,1,l) 
            Szint (imax+3,j,1,l) = Szint (imax+2,j,1,l) 


            psiint(imax+4,j,1,l) = psiint(imax+3,j,1,l)  
            phiint(imax+4,j,1,l) = phiint(imax+3,j,1,l)  
            Bxint (imax+4,j,1,l) = Bxint (imax+3,j,1,l)  
            Byint (imax+4,j,1,l) = Byint (imax+3,j,1,l)  
            Bzint (imax+4,j,1,l) = Bzint (imax+3,j,1,l)  
            qint  (imax+4,j,1,l) = qint  (imax+3,j,1,l)  
            DDint (imax+4,j,1,l) = DDint (imax+3,j,1,l) 
            tauint(imax+4,j,1,l) = tauint(imax+3,j,1,l)
            Sxint (imax+4,j,1,l) = Sxint (imax+3,j,1,l) 
            Syint (imax+4,j,1,l) = Syint (imax+3,j,1,l) 
            Szint (imax+4,j,1,l) = Szint (imax+3,j,1,l)

            psiint(imax+5,j,1,l) = psiint(imax+4,j,1,l)  
            phiint(imax+5,j,1,l) = phiint(imax+4,j,1,l)  
            Bxint (imax+5,j,1,l) = Bxint (imax+4,j,1,l)  
            Byint (imax+5,j,1,l) = Byint (imax+4,j,1,l)  
            Bzint (imax+5,j,1,l) = Bzint (imax+4,j,1,l)  
            qint  (imax+5,j,1,l) = qint  (imax+4,j,1,l)  
            DDint (imax+5,j,1,l) = DDint (imax+4,j,1,l) 
            tauint(imax+5,j,1,l) = tauint(imax+4,j,1,l)
            Sxint (imax+5,j,1,l) = Sxint (imax+4,j,1,l) 
            Syint (imax+5,j,1,l) = Syint (imax+4,j,1,l) 
            Szint (imax+5,j,1,l) = Szint (imax+4,j,1,l)

            psiint(imax+6,j,1,l) = psiint(imax+5,j,1,l)  
            phiint(imax+6,j,1,l) = phiint(imax+5,j,1,l)  
            Bxint (imax+6,j,1,l) = Bxint (imax+5,j,1,l)  
            Byint (imax+6,j,1,l) = Byint (imax+5,j,1,l)  
            Bzint (imax+6,j,1,l) = Bzint (imax+5,j,1,l)  
            qint  (imax+6,j,1,l) = qint  (imax+5,j,1,l)  
            DDint (imax+6,j,1,l) = DDint (imax+5,j,1,l) 
            tauint(imax+6,j,1,l) = tauint(imax+5,j,1,l)
            Sxint (imax+6,j,1,l) = Sxint (imax+5,j,1,l) 
            Syint (imax+6,j,1,l) = Syint (imax+5,j,1,l) 
            Szint (imax+6,j,1,l) = Szint (imax+5,j,1,l)  


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

!!$            psiint(i, 0,1,l) = psiint(i, 1,1,l) 
!!$            phiint(i, 0,1,l) = phiint(i, 1,1,l) 
!!$            Bxint (i, 0,1,l) = Bxint (i, 1,1,l)  
!!$            Byint (i, 0,1,l) = Byint (i, 1,1,l)  
!!$            Bzint (i, 0,1,l) = Bzint (i, 1,1,l)  
!!$            qint  (i, 0,1,l) = qint  (i, 1,1,l)   
!!$            DDint (i, 0,1,l) = DDint (i, 1,1,l)  
!!$            tauint(i, 0,1,l) = tauint(i, 1,1,l) 
!!$            Sxint (i, 0,1,l) = Sxint (i, 1,1,l)  
!!$            Syint (i, 0,1,l) = Syint (i, 1,1,l)  
!!$            Szint (i, 0,1,l) = Szint (i, 1,1,l)  


            psiint(i,-1,1,l) = psiint(i, 0,1,l) 
            phiint(i,-1,1,l) = phiint(i, 0,1,l) 
            Bxint (i,-1,1,l) = Bxint (i, 0,1,l)  
            Byint (i,-1,1,l) = Byint (i, 0,1,l)  
            Bzint (i,-1,1,l) = Bzint (i, 0,1,l)  
            qint  (i,-1,1,l) = qint  (i, 0,1,l)   
            DDint (i,-1,1,l) = DDint (i, 0,1,l)  
            tauint(i,-1,1,l) = tauint(i, 0,1,l) 
            Sxint (i,-1,1,l) = Sxint (i, 0,1,l)  
            Syint (i,-1,1,l) = Syint (i, 0,1,l)  
            Szint (i,-1,1,l) = Szint (i, 0,1,l)  


            psiint(i,-2,1,l) = psiint(i,-1,1,l) 
            phiint(i,-2,1,l) = phiint(i,-1,1,l) 
            Bxint (i,-2,1,l) = Bxint (i,-1,1,l)  
            Byint (i,-2,1,l) = Byint (i,-1,1,l)  
            Bzint (i,-2,1,l) = Bzint (i,-1,1,l)  
            qint  (i,-2,1,l) = qint  (i,-1,1,l)   
            DDint (i,-2,1,l) = DDint (i,-1,1,l)  
            tauint(i,-2,1,l) = tauint(i,-1,1,l) 
            Sxint (i,-2,1,l) = Sxint (i,-1,1,l)  
            Syint (i,-2,1,l) = Syint (i,-1,1,l)  
            Szint (i,-2,1,l) = Szint (i,-1,1,l)  


            psiint(i,-3,1,l) = psiint(i,-2,1,l) 
            phiint(i,-3,1,l) = phiint(i,-2,1,l) 
            Bxint (i,-3,1,l) = Bxint (i,-2,1,l)  
            Byint (i,-3,1,l) = Byint (i,-2,1,l)  
            Bzint (i,-3,1,l) = Bzint (i,-2,1,l)  
            qint  (i,-3,1,l) = qint  (i,-2,1,l)   
            DDint (i,-3,1,l) = DDint (i,-2,1,l)  
            tauint(i,-3,1,l) = tauint(i,-2,1,l) 
            Sxint (i,-3,1,l) = Sxint (i,-2,1,l)  
            Syint (i,-3,1,l) = Syint (i,-2,1,l)  
            Szint (i,-3,1,l) = Szint (i,-2,1,l)  


            psiint(i,-4,1,l) = psiint(i,-3,1,l) 
            phiint(i,-4,1,l) = phiint(i,-3,1,l) 
            Bxint (i,-4,1,l) = Bxint (i,-3,1,l)  
            Byint (i,-4,1,l) = Byint (i,-3,1,l)  
            Bzint (i,-4,1,l) = Bzint (i,-3,1,l)  
            qint  (i,-4,1,l) = qint  (i,-3,1,l)   
            DDint (i,-4,1,l) = DDint (i,-3,1,l)  
            tauint(i,-4,1,l) = tauint(i,-3,1,l) 
            Sxint (i,-4,1,l) = Sxint (i,-3,1,l)  
            Syint (i,-4,1,l) = Syint (i,-3,1,l)  
            Szint (i,-4,1,l) = Szint (i,-3,1,l)

            
            psiint(i,-5,1,l) = psiint(i,-4,1,l) 
            phiint(i,-5,1,l) = phiint(i,-4,1,l) 
            Bxint (i,-5,1,l) = Bxint (i,-4,1,l)  
            Byint (i,-5,1,l) = Byint (i,-4,1,l)  
            Bzint (i,-5,1,l) = Bzint (i,-4,1,l)  
            qint  (i,-5,1,l) = qint  (i,-4,1,l)   
            DDint (i,-5,1,l) = DDint (i,-4,1,l)  
            tauint(i,-5,1,l) = tauint(i,-4,1,l) 
            Sxint (i,-5,1,l) = Sxint (i,-4,1,l)  
            Syint (i,-5,1,l) = Syint (i,-4,1,l)  
            Szint (i,-5,1,l) = Szint (i,-4,1,l)

            
            psiint(i,-6,1,l) = psiint(i,-5,1,l) 
            phiint(i,-6,1,l) = phiint(i,-5,1,l) 
            Bxint (i,-6,1,l) = Bxint (i,-5,1,l)  
            Byint (i,-6,1,l) = Byint (i,-5,1,l)  
            Bzint (i,-6,1,l) = Bzint (i,-5,1,l)  
            qint  (i,-6,1,l) = qint  (i,-5,1,l)   
            DDint (i,-6,1,l) = DDint (i,-5,1,l)  
            tauint(i,-6,1,l) = tauint(i,-5,1,l) 
            Sxint (i,-6,1,l) = Sxint (i,-5,1,l)  
            Syint (i,-6,1,l) = Syint (i,-5,1,l)  
            Szint (i,-6,1,l) = Szint (i,-5,1,l)
            
            ! Right Boundary
            !---------------

!!$            psiint(i,jmax  ,1,l) = psiint(i,jmax-1,1,l) 
!!$            phiint(i,jmax  ,1,l) = phiint(i,jmax-1,1,l) 
!!$            Bxint (i,jmax  ,1,l) = Bxint (i,jmax-1,1,l)  
!!$            Byint (i,jmax  ,1,l) = Byint (i,jmax-1,1,l)  
!!$            Bzint (i,jmax  ,1,l) = Bzint (i,jmax-1,1,l)  
!!$            qint  (i,jmax  ,1,l) = qint  (i,jmax-1,1,l)   
!!$            DDint (i,jmax  ,1,l) = DDint (i,jmax-1,1,l)  
!!$            tauint(i,jmax  ,1,l) = tauint(i,jmax-1,1,l) 
!!$            Sxint (i,jmax  ,1,l) = Sxint (i,jmax-1,1,l)  
!!$            Syint (i,jmax  ,1,l) = Syint (i,jmax-1,1,l)  
!!$            Szint (i,jmax  ,1,l) = Szint (i,jmax-1,1,l)  

            psiint(i,jmax+1,1,l) = psiint(i,jmax  ,1,l) 
            phiint(i,jmax+1,1,l) = phiint(i,jmax  ,1,l) 
            Bxint (i,jmax+1,1,l) = Bxint (i,jmax  ,1,l)  
            Byint (i,jmax+1,1,l) = Byint (i,jmax  ,1,l)  
            Bzint (i,jmax+1,1,l) = Bzint (i,jmax  ,1,l)  
            qint  (i,jmax+1,1,l) = qint  (i,jmax  ,1,l)   
            DDint (i,jmax+1,1,l) = DDint (i,jmax  ,1,l)  
            tauint(i,jmax+1,1,l) = tauint(i,jmax  ,1,l) 
            Sxint (i,jmax+1,1,l) = Sxint (i,jmax  ,1,l)  
            Syint (i,jmax+1,1,l) = Syint (i,jmax  ,1,l)  
            Szint (i,jmax+1,1,l) = Szint (i,jmax  ,1,l)  

            psiint(i,jmax+2,1,l) = psiint(i,jmax+1,1,l) 
            phiint(i,jmax+2,1,l) = phiint(i,jmax+1,1,l) 
            Bxint (i,jmax+2,1,l) = Bxint (i,jmax+1,1,l)  
            Byint (i,jmax+2,1,l) = Byint (i,jmax+1,1,l)  
            Bzint (i,jmax+2,1,l) = Bzint (i,jmax+1,1,l)  
            qint  (i,jmax+2,1,l) = qint  (i,jmax+1,1,l)   
            DDint (i,jmax+2,1,l) = DDint (i,jmax+1,1,l)  
            tauint(i,jmax+2,1,l) = tauint(i,jmax+1,1,l) 
            Sxint (i,jmax+2,1,l) = Sxint (i,jmax+1,1,l)  
            Syint (i,jmax+2,1,l) = Syint (i,jmax+1,1,l)  
            Szint (i,jmax+2,1,l) = Szint (i,jmax+1,1,l) 

            psiint(i,jmax+3,1,l) = psiint(i,jmax+2,1,l) 
            phiint(i,jmax+3,1,l) = phiint(i,jmax+2,1,l) 
            Bxint (i,jmax+3,1,l) = Bxint (i,jmax+2,1,l)  
            Byint (i,jmax+3,1,l) = Byint (i,jmax+2,1,l)  
            Bzint (i,jmax+3,1,l) = Bzint (i,jmax+2,1,l)  
            qint  (i,jmax+3,1,l) = qint  (i,jmax+2,1,l)   
            DDint (i,jmax+3,1,l) = DDint (i,jmax+2,1,l)  
            tauint(i,jmax+3,1,l) = tauint(i,jmax+2,1,l) 
            Sxint (i,jmax+3,1,l) = Sxint (i,jmax+2,1,l)  
            Syint (i,jmax+3,1,l) = Syint (i,jmax+2,1,l)  
            Szint (i,jmax+3,1,l) = Szint (i,jmax+2,1,l) 

            psiint(i,jmax+4,1,l) = psiint(i,jmax+3,1,l) 
            phiint(i,jmax+4,1,l) = phiint(i,jmax+3,1,l) 
            Bxint (i,jmax+4,1,l) = Bxint (i,jmax+3,1,l)  
            Byint (i,jmax+4,1,l) = Byint (i,jmax+3,1,l)  
            Bzint (i,jmax+4,1,l) = Bzint (i,jmax+3,1,l)  
            qint  (i,jmax+4,1,l) = qint  (i,jmax+3,1,l)   
            DDint (i,jmax+4,1,l) = DDint (i,jmax+3,1,l)  
            tauint(i,jmax+4,1,l) = tauint(i,jmax+3,1,l) 
            Sxint (i,jmax+4,1,l) = Sxint (i,jmax+3,1,l)  
            Syint (i,jmax+4,1,l) = Syint (i,jmax+3,1,l)  
            Szint (i,jmax+4,1,l) = Szint (i,jmax+3,1,l)

            psiint(i,jmax+5,1,l) = psiint(i,jmax+4,1,l) 
            phiint(i,jmax+5,1,l) = phiint(i,jmax+4,1,l) 
            Bxint (i,jmax+5,1,l) = Bxint (i,jmax+4,1,l)  
            Byint (i,jmax+5,1,l) = Byint (i,jmax+4,1,l)  
            Bzint (i,jmax+5,1,l) = Bzint (i,jmax+4,1,l)  
            qint  (i,jmax+5,1,l) = qint  (i,jmax+4,1,l)   
            DDint (i,jmax+5,1,l) = DDint (i,jmax+4,1,l)  
            tauint(i,jmax+5,1,l) = tauint(i,jmax+4,1,l) 
            Sxint (i,jmax+5,1,l) = Sxint (i,jmax+4,1,l)  
            Syint (i,jmax+5,1,l) = Syint (i,jmax+4,1,l)  
            Szint (i,jmax+5,1,l) = Szint (i,jmax+4,1,l)

            psiint(i,jmax+6,1,l) = psiint(i,jmax+5,1,l) 
            phiint(i,jmax+6,1,l) = phiint(i,jmax+5,1,l) 
            Bxint (i,jmax+6,1,l) = Bxint (i,jmax+5,1,l)  
            Byint (i,jmax+6,1,l) = Byint (i,jmax+5,1,l)  
            Bzint (i,jmax+6,1,l) = Bzint (i,jmax+5,1,l)  
            qint  (i,jmax+6,1,l) = qint  (i,jmax+5,1,l)   
            DDint (i,jmax+6,1,l) = DDint (i,jmax+5,1,l)  
            tauint(i,jmax+6,1,l) = tauint(i,jmax+5,1,l) 
            Sxint (i,jmax+6,1,l) = Sxint (i,jmax+5,1,l)  
            Syint (i,jmax+6,1,l) = Syint (i,jmax+5,1,l)  
            Szint (i,jmax+6,1,l) = Szint (i,jmax+5,1,l) 


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

!!$            psiint( 0,j,1,l) = psiint(imax-1,j,1,l)
!!$            phiint( 0,j,1,l) = phiint(imax-1,j,1,l)
!!$            Bxint ( 0,j,1,l) = Bxint (imax-1,j,1,l) 
!!$            Byint ( 0,j,1,l) = Byint (imax-1,j,1,l) 
!!$            Bzint ( 0,j,1,l) = Bzint (imax-1,j,1,l) 
!!$            qint  ( 0,j,1,l) = qint  (imax-1,j,1,l) 
!!$            DDint ( 0,j,1,l) = DDint (imax-1,j,1,l)
!!$            tauint( 0,j,1,l) = tauint(imax-1,j,1,l)
!!$            Sxint ( 0,j,1,l) = Sxint (imax-1,j,1,l) 
!!$            Syint ( 0,j,1,l) = Syint (imax-1,j,1,l) 
!!$            Szint ( 0,j,1,l) = Szint (imax-1,j,1,l) 

            psiint(-1,j,1,l) = psiint(imax-1,j,1,l)
            phiint(-1,j,1,l) = phiint(imax-1,j,1,l)
            Bxint (-1,j,1,l) = Bxint (imax-1,j,1,l) 
            Byint (-1,j,1,l) = Byint (imax-1,j,1,l) 
            Bzint (-1,j,1,l) = Bzint (imax-1,j,1,l) 
            qint  (-1,j,1,l) = qint  (imax-1,j,1,l) 
            DDint (-1,j,1,l) = DDint (imax-1,j,1,l)
            tauint(-1,j,1,l) = tauint(imax-1,j,1,l)
            Sxint (-1,j,1,l) = Sxint (imax-1,j,1,l) 
            Syint (-1,j,1,l) = Syint (imax-1,j,1,l) 
            Szint (-1,j,1,l) = Szint (imax-1,j,1,l) 

            psiint(-2,j,1,l) = psiint(imax-2,j,1,l)
            phiint(-2,j,1,l) = phiint(imax-2,j,1,l)
            Bxint (-2,j,1,l) = Bxint (imax-2,j,1,l) 
            Byint (-2,j,1,l) = Byint (imax-2,j,1,l) 
            Bzint (-2,j,1,l) = Bzint (imax-2,j,1,l) 
            qint  (-2,j,1,l) = qint  (imax-2,j,1,l) 
            DDint (-2,j,1,l) = DDint (imax-2,j,1,l)
            tauint(-2,j,1,l) = tauint(imax-2,j,1,l)
            Sxint (-2,j,1,l) = Sxint (imax-2,j,1,l) 
            Syint (-2,j,1,l) = Syint (imax-2,j,1,l) 
            Szint (-2,j,1,l) = Szint (imax-2,j,1,l) 

            psiint(-3,j,1,l) = psiint(imax-3,j,1,l)
            phiint(-3,j,1,l) = phiint(imax-3,j,1,l)
            Bxint (-3,j,1,l) = Bxint (imax-3,j,1,l) 
            Byint (-3,j,1,l) = Byint (imax-3,j,1,l) 
            Bzint (-3,j,1,l) = Bzint (imax-3,j,1,l) 
            qint  (-3,j,1,l) = qint  (imax-3,j,1,l) 
            DDint (-3,j,1,l) = DDint (imax-3,j,1,l)
            tauint(-3,j,1,l) = tauint(imax-3,j,1,l)
            Sxint (-3,j,1,l) = Sxint (imax-3,j,1,l) 
            Syint (-3,j,1,l) = Syint (imax-3,j,1,l) 
            Szint (-3,j,1,l) = Szint (imax-3,j,1,l) 

            psiint(-4,j,1,l) = psiint(imax-4,j,1,l)
            phiint(-4,j,1,l) = phiint(imax-4,j,1,l)
            Bxint (-4,j,1,l) = Bxint (imax-4,j,1,l) 
            Byint (-4,j,1,l) = Byint (imax-4,j,1,l) 
            Bzint (-4,j,1,l) = Bzint (imax-4,j,1,l) 
            qint  (-4,j,1,l) = qint  (imax-4,j,1,l) 
            DDint (-4,j,1,l) = DDint (imax-4,j,1,l)
            tauint(-4,j,1,l) = tauint(imax-4,j,1,l)
            Sxint (-4,j,1,l) = Sxint (imax-4,j,1,l) 
            Syint (-4,j,1,l) = Syint (imax-4,j,1,l) 
            Szint (-4,j,1,l) = Szint (imax-4,j,1,l)

            psiint(-5,j,1,l) = psiint(imax-5,j,1,l)
            phiint(-5,j,1,l) = phiint(imax-5,j,1,l)
            Bxint (-5,j,1,l) = Bxint (imax-5,j,1,l) 
            Byint (-5,j,1,l) = Byint (imax-5,j,1,l) 
            Bzint (-5,j,1,l) = Bzint (imax-5,j,1,l) 
            qint  (-5,j,1,l) = qint  (imax-5,j,1,l) 
            DDint (-5,j,1,l) = DDint (imax-5,j,1,l)
            tauint(-5,j,1,l) = tauint(imax-5,j,1,l)
            Sxint (-5,j,1,l) = Sxint (imax-5,j,1,l) 
            Syint (-5,j,1,l) = Syint (imax-5,j,1,l) 
            Szint (-5,j,1,l) = Szint (imax-5,j,1,l)

            psiint(-6,j,1,l) = psiint(imax-6,j,1,l)
            phiint(-6,j,1,l) = phiint(imax-6,j,1,l)
            Bxint (-6,j,1,l) = Bxint (imax-6,j,1,l) 
            Byint (-6,j,1,l) = Byint (imax-6,j,1,l) 
            Bzint (-6,j,1,l) = Bzint (imax-6,j,1,l) 
            qint  (-6,j,1,l) = qint  (imax-6,j,1,l) 
            DDint (-6,j,1,l) = DDint (imax-6,j,1,l)
            tauint(-6,j,1,l) = tauint(imax-6,j,1,l)
            Sxint (-6,j,1,l) = Sxint (imax-6,j,1,l) 
            Syint (-6,j,1,l) = Syint (imax-6,j,1,l) 
            Szint (-6,j,1,l) = Szint (imax-6,j,1,l) 

            ! Right Boundary
            !---------------

!!$            psiint(imax  ,j,1,l) = psiint(1,j,1,l)  
!!$            phiint(imax  ,j,1,l) = phiint(1,j,1,l)  
!!$            Bxint (imax  ,j,1,l) = Bxint (1,j,1,l)  
!!$            Byint (imax  ,j,1,l) = Byint (1,j,1,l)  
!!$            Bzint (imax  ,j,1,l) = Bzint (1,j,1,l)  
!!$            qint  (imax  ,j,1,l) = qint  (1,j,1,l)  
!!$            DDint (imax  ,j,1,l) = DDint (1,j,1,l) 
!!$            tauint(imax  ,j,1,l) = tauint(1,j,1,l)
!!$            Sxint (imax  ,j,1,l) = Sxint (1,j,1,l) 
!!$            Syint (imax  ,j,1,l) = Syint (1,j,1,l) 
!!$            Szint (imax  ,j,1,l) = Szint (1,j,1,l) 

            psiint(imax+1,j,1,l) = psiint(1,j,1,l) 
            phiint(imax+1,j,1,l) = phiint(1,j,1,l) 
            Bxint (imax+1,j,1,l) = Bxint (1,j,1,l)  
            Byint (imax+1,j,1,l) = Byint (1,j,1,l)  
            Bzint (imax+1,j,1,l) = Bzint (1,j,1,l)  
            qint  (imax+1,j,1,l) = qint  (1,j,1,l)  
            DDint (imax+1,j,1,l) = DDint (1,j,1,l) 
            tauint(imax+1,j,1,l) = tauint(1,j,1,l)
            Sxint (imax+1,j,1,l) = Sxint (1,j,1,l) 
            Syint (imax+1,j,1,l) = Syint (1,j,1,l) 
            Szint (imax+1,j,1,l) = Szint (1,j,1,l) 

            psiint(imax+2,j,1,l) = psiint(2,j,1,l)  
            phiint(imax+2,j,1,l) = phiint(2,j,1,l)  
            Bxint (imax+2,j,1,l) = Bxint (2,j,1,l)  
            Byint (imax+2,j,1,l) = Byint (2,j,1,l)  
            Bzint (imax+2,j,1,l) = Bzint (2,j,1,l)  
            qint  (imax+2,j,1,l) = qint  (2,j,1,l)  
            DDint (imax+2,j,1,l) = DDint (2,j,1,l) 
            tauint(imax+2,j,1,l) = tauint(2,j,1,l)
            Sxint (imax+2,j,1,l) = Sxint (2,j,1,l) 
            Syint (imax+2,j,1,l) = Syint (2,j,1,l) 
            Szint (imax+2,j,1,l) = Szint (2,j,1,l) 

            psiint(imax+3,j,1,l) = psiint(3,j,1,l)  
            phiint(imax+3,j,1,l) = phiint(3,j,1,l)  
            Bxint (imax+3,j,1,l) = Bxint (3,j,1,l)  
            Byint (imax+3,j,1,l) = Byint (3,j,1,l)  
            Bzint (imax+3,j,1,l) = Bzint (3,j,1,l)  
            qint  (imax+3,j,1,l) = qint  (3,j,1,l)  
            DDint (imax+3,j,1,l) = DDint (3,j,1,l) 
            tauint(imax+3,j,1,l) = tauint(3,j,1,l)
            Sxint (imax+3,j,1,l) = Sxint (3,j,1,l) 
            Syint (imax+3,j,1,l) = Syint (3,j,1,l) 
            Szint (imax+3,j,1,l) = Szint (3,j,1,l) 

            psiint(imax+4,j,1,l) = psiint(4,j,1,l)  
            phiint(imax+4,j,1,l) = phiint(4,j,1,l)  
            Bxint (imax+4,j,1,l) = Bxint (4,j,1,l)  
            Byint (imax+4,j,1,l) = Byint (4,j,1,l)  
            Bzint (imax+4,j,1,l) = Bzint (4,j,1,l)  
            qint  (imax+4,j,1,l) = qint  (4,j,1,l)  
            DDint (imax+4,j,1,l) = DDint (4,j,1,l) 
            tauint(imax+4,j,1,l) = tauint(4,j,1,l)
            Sxint (imax+4,j,1,l) = Sxint (4,j,1,l) 
            Syint (imax+4,j,1,l) = Syint (4,j,1,l) 
            Szint (imax+4,j,1,l) = Szint (4,j,1,l)

            psiint(imax+5,j,1,l) = psiint(5,j,1,l)  
            phiint(imax+5,j,1,l) = phiint(5,j,1,l)  
            Bxint (imax+5,j,1,l) = Bxint (5,j,1,l)  
            Byint (imax+5,j,1,l) = Byint (5,j,1,l)  
            Bzint (imax+5,j,1,l) = Bzint (5,j,1,l)  
            qint  (imax+5,j,1,l) = qint  (5,j,1,l)  
            DDint (imax+5,j,1,l) = DDint (5,j,1,l) 
            tauint(imax+5,j,1,l) = tauint(5,j,1,l)
            Sxint (imax+5,j,1,l) = Sxint (5,j,1,l) 
            Syint (imax+5,j,1,l) = Syint (5,j,1,l) 
            Szint (imax+5,j,1,l) = Szint (5,j,1,l)

            psiint(imax+6,j,1,l) = psiint(6,j,1,l)  
            phiint(imax+6,j,1,l) = phiint(6,j,1,l)  
            Bxint (imax+6,j,1,l) = Bxint (6,j,1,l)  
            Byint (imax+6,j,1,l) = Byint (6,j,1,l)  
            Bzint (imax+6,j,1,l) = Bzint (6,j,1,l)  
            qint  (imax+6,j,1,l) = qint  (6,j,1,l)  
            DDint (imax+6,j,1,l) = DDint (6,j,1,l) 
            tauint(imax+6,j,1,l) = tauint(6,j,1,l)
            Sxint (imax+6,j,1,l) = Sxint (6,j,1,l) 
            Syint (imax+6,j,1,l) = Syint (6,j,1,l) 
            Szint (imax+6,j,1,l) = Szint (6,j,1,l) 

 
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

!!$            psiint(i, 0,1,l) = psiint(i,jmax-1,1,l) 
!!$            phiint(i, 0,1,l) = phiint(i,jmax-1,1,l) 
!!$            Bxint (i, 0,1,l) = Bxint (i,jmax-1,1,l)  
!!$            Byint (i, 0,1,l) = Byint (i,jmax-1,1,l)  
!!$            Bzint (i, 0,1,l) = Bzint (i,jmax-1,1,l)  
!!$            qint  (i, 0,1,l) = qint  (i,jmax-1,1,l)   
!!$            DDint (i, 0,1,l) = DDint (i,jmax-1,1,l)  
!!$            tauint(i, 0,1,l) = tauint(i,jmax-1,1,l) 
!!$            Sxint (i, 0,1,l) = Sxint (i,jmax-1,1,l)  
!!$            Syint (i, 0,1,l) = Syint (i,jmax-1,1,l)  
!!$            Szint (i, 0,1,l) = Szint (i,jmax-1,1,l)  


            psiint(i,-1,1,l) = psiint(i,jmax-1,1,l) 
            phiint(i,-1,1,l) = phiint(i,jmax-1,1,l) 
            Bxint (i,-1,1,l) = Bxint (i,jmax-1,1,l)  
            Byint (i,-1,1,l) = Byint (i,jmax-1,1,l)  
            Bzint (i,-1,1,l) = Bzint (i,jmax-1,1,l)  
            qint  (i,-1,1,l) = qint  (i,jmax-1,1,l)   
            DDint (i,-1,1,l) = DDint (i,jmax-1,1,l)  
            tauint(i,-1,1,l) = tauint(i,jmax-1,1,l) 
            Sxint (i,-1,1,l) = Sxint (i,jmax-1,1,l)  
            Syint (i,-1,1,l) = Syint (i,jmax-1,1,l)  
            Szint (i,-1,1,l) = Szint (i,jmax-1,1,l)  


            psiint(i,-2,1,l) = psiint(i,jmax-2,1,l) 
            phiint(i,-2,1,l) = phiint(i,jmax-2,1,l) 
            Bxint (i,-2,1,l) = Bxint (i,jmax-2,1,l)  
            Byint (i,-2,1,l) = Byint (i,jmax-2,1,l)  
            Bzint (i,-2,1,l) = Bzint (i,jmax-2,1,l)  
            qint  (i,-2,1,l) = qint  (i,jmax-2,1,l)   
            DDint (i,-2,1,l) = DDint (i,jmax-2,1,l)  
            tauint(i,-2,1,l) = tauint(i,jmax-2,1,l) 
            Sxint (i,-2,1,l) = Sxint (i,jmax-2,1,l)  
            Syint (i,-2,1,l) = Syint (i,jmax-2,1,l)  
            Szint (i,-2,1,l) = Szint (i,jmax-2,1,l)  

            psiint(i,-3,1,l) = psiint(i,jmax-3,1,l) 
            phiint(i,-3,1,l) = phiint(i,jmax-3,1,l) 
            Bxint (i,-3,1,l) = Bxint (i,jmax-3,1,l)  
            Byint (i,-3,1,l) = Byint (i,jmax-3,1,l)  
            Bzint (i,-3,1,l) = Bzint (i,jmax-3,1,l)  
            qint  (i,-3,1,l) = qint  (i,jmax-3,1,l)   
            DDint (i,-3,1,l) = DDint (i,jmax-3,1,l)  
            tauint(i,-3,1,l) = tauint(i,jmax-3,1,l) 
            Sxint (i,-3,1,l) = Sxint (i,jmax-3,1,l)  
            Syint (i,-3,1,l) = Syint (i,jmax-3,1,l)  
            Szint (i,-3,1,l) = Szint (i,jmax-3,1,l)  

            psiint(i,-4,1,l) = psiint(i,jmax-4,1,l) 
            phiint(i,-4,1,l) = phiint(i,jmax-4,1,l) 
            Bxint (i,-4,1,l) = Bxint (i,jmax-4,1,l)  
            Byint (i,-4,1,l) = Byint (i,jmax-4,1,l)  
            Bzint (i,-4,1,l) = Bzint (i,jmax-4,1,l)  
            qint  (i,-4,1,l) = qint  (i,jmax-4,1,l)   
            DDint (i,-4,1,l) = DDint (i,jmax-4,1,l)  
            tauint(i,-4,1,l) = tauint(i,jmax-4,1,l) 
            Sxint (i,-4,1,l) = Sxint (i,jmax-4,1,l)  
            Syint (i,-4,1,l) = Syint (i,jmax-4,1,l)  
            Szint (i,-4,1,l) = Szint (i,jmax-4,1,l)
            
            psiint(i,-5,1,l) = psiint(i,jmax-5,1,l) 
            phiint(i,-5,1,l) = phiint(i,jmax-5,1,l) 
            Bxint (i,-5,1,l) = Bxint (i,jmax-5,1,l)  
            Byint (i,-5,1,l) = Byint (i,jmax-5,1,l)  
            Bzint (i,-5,1,l) = Bzint (i,jmax-5,1,l)  
            qint  (i,-5,1,l) = qint  (i,jmax-5,1,l)   
            DDint (i,-5,1,l) = DDint (i,jmax-5,1,l)  
            tauint(i,-5,1,l) = tauint(i,jmax-5,1,l) 
            Sxint (i,-5,1,l) = Sxint (i,jmax-5,1,l)  
            Syint (i,-5,1,l) = Syint (i,jmax-5,1,l)  
            Szint (i,-5,1,l) = Szint (i,jmax-5,1,l)

            psiint(i,-6,1,l) = psiint(i,jmax-6,1,l) 
            phiint(i,-6,1,l) = phiint(i,jmax-6,1,l) 
            Bxint (i,-6,1,l) = Bxint (i,jmax-6,1,l)  
            Byint (i,-6,1,l) = Byint (i,jmax-6,1,l)  
            Bzint (i,-6,1,l) = Bzint (i,jmax-6,1,l)  
            qint  (i,-6,1,l) = qint  (i,jmax-6,1,l)   
            DDint (i,-6,1,l) = DDint (i,jmax-6,1,l)  
            tauint(i,-6,1,l) = tauint(i,jmax-6,1,l) 
            Sxint (i,-6,1,l) = Sxint (i,jmax-6,1,l)  
            Syint (i,-6,1,l) = Syint (i,jmax-6,1,l)  
            Szint (i,-6,1,l) = Szint (i,jmax-6,1,l)  

            ! Right Boundary
            !---------------

!!$            psiint(i,jmax  ,1,l) = psiint(i,1,1,l) 
!!$            phiint(i,jmax  ,1,l) = phiint(i,1,1,l) 
!!$            Bxint (i,jmax  ,1,l) = Bxint (i,1,1,l)  
!!$            Byint (i,jmax  ,1,l) = Byint (i,1,1,l)  
!!$            Bzint (i,jmax  ,1,l) = Bzint (i,1,1,l)  
!!$            qint  (i,jmax  ,1,l) = qint  (i,1,1,l)   
!!$            DDint (i,jmax  ,1,l) = DDint (i,1,1,l)  
!!$            tauint(i,jmax  ,1,l) = tauint(i,1,1,l) 
!!$            Sxint (i,jmax  ,1,l) = Sxint (i,1,1,l)  
!!$            Syint (i,jmax  ,1,l) = Syint (i,1,1,l)  
!!$            Szint (i,jmax  ,1,l) = Szint (i,1,1,l)  

            psiint(i,jmax+1,1,l) = psiint(i,1,1,l) 
            phiint(i,jmax+1,1,l) = phiint(i,1,1,l) 
            Bxint (i,jmax+1,1,l) = Bxint (i,1,1,l)  
            Byint (i,jmax+1,1,l) = Byint (i,1,1,l)  
            Bzint (i,jmax+1,1,l) = Bzint (i,1,1,l)  
            qint  (i,jmax+1,1,l) = qint  (i,1,1,l)   
            DDint (i,jmax+1,1,l) = DDint (i,1,1,l)  
            tauint(i,jmax+1,1,l) = tauint(i,1,1,l) 
            Sxint (i,jmax+1,1,l) = Sxint (i,1,1,l)  
            Syint (i,jmax+1,1,l) = Syint (i,1,1,l)  
            Szint (i,jmax+1,1,l) = Szint (i,1,1,l)  

 
            psiint(i,jmax+2,1,l) = psiint(i,2,1,l) 
            phiint(i,jmax+2,1,l) = phiint(i,2,1,l) 
            Bxint (i,jmax+2,1,l) = Bxint (i,2,1,l)  
            Byint (i,jmax+2,1,l) = Byint (i,2,1,l)  
            Bzint (i,jmax+2,1,l) = Bzint (i,2,1,l)  
            qint  (i,jmax+2,1,l) = qint  (i,2,1,l)   
            DDint (i,jmax+2,1,l) = DDint (i,2,1,l)  
            tauint(i,jmax+2,1,l) = tauint(i,2,1,l) 
            Sxint (i,jmax+2,1,l) = Sxint (i,2,1,l)  
            Syint (i,jmax+2,1,l) = Syint (i,2,1,l)  
            Szint (i,jmax+2,1,l) = Szint (i,2,1,l) 

            psiint(i,jmax+3,1,l) = psiint(i,3,1,l) 
            phiint(i,jmax+3,1,l) = phiint(i,3,1,l) 
            Bxint (i,jmax+3,1,l) = Bxint (i,3,1,l)  
            Byint (i,jmax+3,1,l) = Byint (i,3,1,l)  
            Bzint (i,jmax+3,1,l) = Bzint (i,3,1,l)  
            qint  (i,jmax+3,1,l) = qint  (i,3,1,l)   
            DDint (i,jmax+3,1,l) = DDint (i,3,1,l)  
            tauint(i,jmax+3,1,l) = tauint(i,3,1,l) 
            Sxint (i,jmax+3,1,l) = Sxint (i,3,1,l)  
            Syint (i,jmax+3,1,l) = Syint (i,3,1,l)  
            Szint (i,jmax+3,1,l) = Szint (i,3,1,l) 

            psiint(i,jmax+4,1,l) = psiint(i,4,1,l) 
            phiint(i,jmax+4,1,l) = phiint(i,4,1,l) 
            Bxint (i,jmax+4,1,l) = Bxint (i,4,1,l)  
            Byint (i,jmax+4,1,l) = Byint (i,4,1,l)  
            Bzint (i,jmax+4,1,l) = Bzint (i,4,1,l)  
            qint  (i,jmax+4,1,l) = qint  (i,4,1,l)   
            DDint (i,jmax+4,1,l) = DDint (i,4,1,l)  
            tauint(i,jmax+4,1,l) = tauint(i,4,1,l) 
            Sxint (i,jmax+4,1,l) = Sxint (i,4,1,l)  
            Syint (i,jmax+4,1,l) = Syint (i,4,1,l)  
            Szint (i,jmax+4,1,l) = Szint (i,4,1,l)

            psiint(i,jmax+5,1,l) = psiint(i,5,1,l) 
            phiint(i,jmax+5,1,l) = phiint(i,5,1,l) 
            Bxint (i,jmax+5,1,l) = Bxint (i,5,1,l)  
            Byint (i,jmax+5,1,l) = Byint (i,5,1,l)  
            Bzint (i,jmax+5,1,l) = Bzint (i,5,1,l)  
            qint  (i,jmax+5,1,l) = qint  (i,5,1,l)   
            DDint (i,jmax+5,1,l) = DDint (i,5,1,l)  
            tauint(i,jmax+5,1,l) = tauint(i,5,1,l) 
            Sxint (i,jmax+5,1,l) = Sxint (i,5,1,l)  
            Syint (i,jmax+5,1,l) = Syint (i,5,1,l)  
            Szint (i,jmax+5,1,l) = Szint (i,5,1,l)

            psiint(i,jmax+6,1,l) = psiint(i,6,1,l) 
            phiint(i,jmax+6,1,l) = phiint(i,6,1,l) 
            Bxint (i,jmax+6,1,l) = Bxint (i,6,1,l)  
            Byint (i,jmax+6,1,l) = Byint (i,6,1,l)  
            Bzint (i,jmax+6,1,l) = Bzint (i,6,1,l)  
            qint  (i,jmax+6,1,l) = qint  (i,6,1,l)   
            DDint (i,jmax+6,1,l) = DDint (i,6,1,l)  
            tauint(i,jmax+6,1,l) = tauint(i,6,1,l) 
            Sxint (i,jmax+6,1,l) = Sxint (i,6,1,l)  
            Syint (i,jmax+6,1,l) = Syint (i,6,1,l)  
            Szint (i,jmax+6,1,l) = Szint (i,6,1,l) 


  
!         end do
      end do

 !$OMP END DO


    else if (DIM == 2 .and. BOUND == 3) then

       ! Reconnection boundary 

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO  

     do j=-6,jmax+6
!      do k=1,kmax

 

        ! Left Boundary
        !---------------

            psiint( 0,j,1,l) = psi(0,j,1)  !psiint( 1,j,1,l)
            phiint( 0,j,1,l) = phi(0,j,1)  !phiint( 1,j,1,l)
            Bxint ( 0,j,1,l) = Bx (0,j,1)  !Bxint ( 1,j,1,l) 
            Byint ( 0,j,1,l) = By (0,j,1)  !Byint ( 1,j,1,l) 
            Bzint ( 0,j,1,l) = Bz (0,j,1)  !Bzint ( 1,j,1,l) 
            qint  ( 0,j,1,l) = q  (0,j,1)  !qint  ( 1,j,1,l) 
            DDint ( 0,j,1,l) = D  (0,j,1)  !DDint ( 1,j,1,l)
            tauint( 0,j,1,l) = tau(0,j,1)  !tauint( 1,j,1,l)
            Sxint ( 0,j,1,l) = Sx (0,j,1)  !Sxint ( 1,j,1,l) 
            Syint ( 0,j,1,l) = Sy (0,j,1)  !Syint ( 1,j,1,l) 
            Szint ( 0,j,1,l) = Sz (0,j,1)  !Szint ( 1,j,1,l) 


            psiint(-1,j,1,l) = psiint( 0,j,1,l)
            phiint(-1,j,1,l) = phiint( 0,j,1,l)
            Bxint (-1,j,1,l) = Bxint ( 0,j,1,l) 
            Byint (-1,j,1,l) = Byint ( 0,j,1,l) 
            Bzint (-1,j,1,l) = Bzint ( 0,j,1,l) 
            qint  (-1,j,1,l) = qint  ( 0,j,1,l) 
            DDint (-1,j,1,l) = DDint ( 0,j,1,l)
            tauint(-1,j,1,l) = tauint( 0,j,1,l)
            Sxint (-1,j,1,l) = Sxint ( 0,j,1,l) 
            Syint (-1,j,1,l) = Syint ( 0,j,1,l) 
            Szint (-1,j,1,l) = Szint ( 0,j,1,l) 


            psiint(-2,j,1,l) = psiint(-1,j,1,l)
            phiint(-2,j,1,l) = phiint(-1,j,1,l)
            Bxint (-2,j,1,l) = Bxint (-1,j,1,l) 
            Byint (-2,j,1,l) = Byint (-1,j,1,l) 
            Bzint (-2,j,1,l) = Bzint (-1,j,1,l) 
            qint  (-2,j,1,l) = qint  (-1,j,1,l) 
            DDint (-2,j,1,l) = DDint (-1,j,1,l)
            tauint(-2,j,1,l) = tauint(-1,j,1,l)
            Sxint (-2,j,1,l) = Sxint (-1,j,1,l) 
            Syint (-2,j,1,l) = Syint (-1,j,1,l) 
            Szint (-2,j,1,l) = Szint (-1,j,1,l) 


            psiint(-3,j,1,l) = psiint(-2,j,1,l)
            phiint(-3,j,1,l) = phiint(-2,j,1,l)
            Bxint (-3,j,1,l) = Bxint (-2,j,1,l) 
            Byint (-3,j,1,l) = Byint (-2,j,1,l) 
            Bzint (-3,j,1,l) = Bzint (-2,j,1,l) 
            qint  (-3,j,1,l) = qint  (-2,j,1,l) 
            DDint (-3,j,1,l) = DDint (-2,j,1,l)
            tauint(-3,j,1,l) = tauint(-2,j,1,l)
            Sxint (-3,j,1,l) = Sxint (-2,j,1,l) 
            Syint (-3,j,1,l) = Syint (-2,j,1,l) 
            Szint (-3,j,1,l) = Szint (-2,j,1,l) 


            psiint(-4,j,1,l) = psiint(-3,j,1,l)
            phiint(-4,j,1,l) = phiint(-3,j,1,l)
            Bxint (-4,j,1,l) = Bxint (-3,j,1,l) 
            Byint (-4,j,1,l) = Byint (-3,j,1,l) 
            Bzint (-4,j,1,l) = Bzint (-3,j,1,l) 
            qint  (-4,j,1,l) = qint  (-3,j,1,l) 
            DDint (-4,j,1,l) = DDint (-3,j,1,l)
            tauint(-4,j,1,l) = tauint(-3,j,1,l)
            Sxint (-4,j,1,l) = Sxint (-3,j,1,l) 
            Syint (-4,j,1,l) = Syint (-3,j,1,l) 
            Szint (-4,j,1,l) = Szint (-3,j,1,l)


            psiint(-5,j,1,l) = psiint(-4,j,1,l)
            phiint(-5,j,1,l) = phiint(-4,j,1,l)
            Bxint (-5,j,1,l) = Bxint (-4,j,1,l) 
            Byint (-5,j,1,l) = Byint (-4,j,1,l) 
            Bzint (-5,j,1,l) = Bzint (-4,j,1,l) 
            qint  (-5,j,1,l) = qint  (-4,j,1,l) 
            DDint (-5,j,1,l) = DDint (-4,j,1,l)
            tauint(-5,j,1,l) = tauint(-4,j,1,l)
            Sxint (-5,j,1,l) = Sxint (-4,j,1,l) 
            Syint (-5,j,1,l) = Syint (-4,j,1,l) 
            Szint (-5,j,1,l) = Szint (-4,j,1,l)


            psiint(-6,j,1,l) = psiint(-5,j,1,l)
            phiint(-6,j,1,l) = phiint(-5,j,1,l)
            Bxint (-6,j,1,l) = Bxint (-5,j,1,l) 
            Byint (-6,j,1,l) = Byint (-5,j,1,l) 
            Bzint (-6,j,1,l) = Bzint (-5,j,1,l) 
            qint  (-6,j,1,l) = qint  (-5,j,1,l) 
            DDint (-6,j,1,l) = DDint (-5,j,1,l)
            tauint(-6,j,1,l) = tauint(-5,j,1,l)
            Sxint (-6,j,1,l) = Sxint (-5,j,1,l) 
            Syint (-6,j,1,l) = Syint (-5,j,1,l) 
            Szint (-6,j,1,l) = Szint (-5,j,1,l) 

            ! Right Boundary
            !---------------
            

            psiint(imax  ,j,1,l) = psi(imax,j  ,1) !psiint(imax-1,j,1,l) 
            phiint(imax  ,j,1,l) = phi(imax,j  ,1) !phiint(imax-1,j,1,l) 
            Bxint (imax  ,j,1,l) = Bx (imax,j  ,1) !Bxint (imax-1,j,1,l)  
            Byint (imax  ,j,1,l) = By (imax,j  ,1) !Byint (imax-1,j,1,l)  
            Bzint (imax  ,j,1,l) = Bz (imax,j  ,1) !Bzint (imax-1,j,1,l)  
            qint  (imax  ,j,1,l) = q  (imax,j  ,1) !qint  (imax-1,j,1,l)  
            DDint (imax  ,j,1,l) = D  (imax,j  ,1) !DDint (imax-1,j,1,l) 
            tauint(imax  ,j,1,l) = tau(imax,j  ,1) !tauint(imax-1,j,1,l)
            Sxint (imax  ,j,1,l) = Sx (imax,j  ,1) !Sxint (imax-1,j,1,l) 
            Syint (imax  ,j,1,l) = Sy (imax,j  ,1) !Syint (imax-1,j,1,l) 
            Szint (imax  ,j,1,l) = Sz (imax,j  ,1) !Szint (imax-1,j,1,l) 


            psiint(imax+1,j,1,l) = psiint(imax  ,j,1,l) 
            phiint(imax+1,j,1,l) = phiint(imax  ,j,1,l) 
            Bxint (imax+1,j,1,l) = Bxint (imax  ,j,1,l)  
            Byint (imax+1,j,1,l) = Byint (imax  ,j,1,l)  
            Bzint (imax+1,j,1,l) = Bzint (imax  ,j,1,l)  
            qint  (imax+1,j,1,l) = qint  (imax  ,j,1,l)  
            DDint (imax+1,j,1,l) = DDint (imax  ,j,1,l) 
            tauint(imax+1,j,1,l) = tauint(imax  ,j,1,l)
            Sxint (imax+1,j,1,l) = Sxint (imax  ,j,1,l) 
            Syint (imax+1,j,1,l) = Syint (imax  ,j,1,l) 
            Szint (imax+1,j,1,l) = Szint (imax  ,j,1,l) 

            psiint(imax+2,j,1,l) = psiint(imax+1,j,1,l)  
            phiint(imax+2,j,1,l) = phiint(imax+1,j,1,l)  
            Bxint (imax+2,j,1,l) = Bxint (imax+1,j,1,l)  
            Byint (imax+2,j,1,l) = Byint (imax+1,j,1,l)  
            Bzint (imax+2,j,1,l) = Bzint (imax+1,j,1,l)  
            qint  (imax+2,j,1,l) = qint  (imax+1,j,1,l)  
            DDint (imax+2,j,1,l) = DDint (imax+1,j,1,l) 
            tauint(imax+2,j,1,l) = tauint(imax+1,j,1,l)
            Sxint (imax+2,j,1,l) = Sxint (imax+1,j,1,l) 
            Syint (imax+2,j,1,l) = Syint (imax+1,j,1,l) 
            Szint (imax+2,j,1,l) = Szint (imax+1,j,1,l) 

            psiint(imax+3,j,1,l) = psiint(imax+2,j,1,l)  
            phiint(imax+3,j,1,l) = phiint(imax+2,j,1,l)  
            Bxint (imax+3,j,1,l) = Bxint (imax+2,j,1,l)  
            Byint (imax+3,j,1,l) = Byint (imax+2,j,1,l)  
            Bzint (imax+3,j,1,l) = Bzint (imax+2,j,1,l)  
            qint  (imax+3,j,1,l) = qint  (imax+2,j,1,l)  
            DDint (imax+3,j,1,l) = DDint (imax+2,j,1,l) 
            tauint(imax+3,j,1,l) = tauint(imax+2,j,1,l)
            Sxint (imax+3,j,1,l) = Sxint (imax+2,j,1,l) 
            Syint (imax+3,j,1,l) = Syint (imax+2,j,1,l) 
            Szint (imax+3,j,1,l) = Szint (imax+2,j,1,l) 

            psiint(imax+4,j,1,l) = psiint(imax+3,j,1,l)  
            phiint(imax+4,j,1,l) = phiint(imax+3,j,1,l)  
            Bxint (imax+4,j,1,l) = Bxint (imax+3,j,1,l)  
            Byint (imax+4,j,1,l) = Byint (imax+3,j,1,l)  
            Bzint (imax+4,j,1,l) = Bzint (imax+3,j,1,l)  
            qint  (imax+4,j,1,l) = qint  (imax+3,j,1,l)  
            DDint (imax+4,j,1,l) = DDint (imax+3,j,1,l) 
            tauint(imax+4,j,1,l) = tauint(imax+3,j,1,l)
            Sxint (imax+4,j,1,l) = Sxint (imax+3,j,1,l) 
            Syint (imax+4,j,1,l) = Syint (imax+3,j,1,l) 
            Szint (imax+4,j,1,l) = Szint (imax+3,j,1,l)

            psiint(imax+5,j,1,l) = psiint(imax+4,j,1,l)  
            phiint(imax+5,j,1,l) = phiint(imax+4,j,1,l)  
            Bxint (imax+5,j,1,l) = Bxint (imax+4,j,1,l)  
            Byint (imax+5,j,1,l) = Byint (imax+4,j,1,l)  
            Bzint (imax+5,j,1,l) = Bzint (imax+4,j,1,l)  
            qint  (imax+5,j,1,l) = qint  (imax+4,j,1,l)  
            DDint (imax+5,j,1,l) = DDint (imax+4,j,1,l) 
            tauint(imax+5,j,1,l) = tauint(imax+4,j,1,l)
            Sxint (imax+5,j,1,l) = Sxint (imax+4,j,1,l) 
            Syint (imax+5,j,1,l) = Syint (imax+4,j,1,l) 
            Szint (imax+5,j,1,l) = Szint (imax+4,j,1,l)

            psiint(imax+6,j,1,l) = psiint(imax+5,j,1,l)  
            phiint(imax+6,j,1,l) = phiint(imax+5,j,1,l)  
            Bxint (imax+6,j,1,l) = Bxint (imax+5,j,1,l)  
            Byint (imax+6,j,1,l) = Byint (imax+5,j,1,l)  
            Bzint (imax+6,j,1,l) = Bzint (imax+5,j,1,l)  
            qint  (imax+6,j,1,l) = qint  (imax+5,j,1,l)  
            DDint (imax+6,j,1,l) = DDint (imax+5,j,1,l) 
            tauint(imax+6,j,1,l) = tauint(imax+5,j,1,l)
            Sxint (imax+6,j,1,l) = Sxint (imax+5,j,1,l) 
            Syint (imax+6,j,1,l) = Syint (imax+5,j,1,l) 
            Szint (imax+6,j,1,l) = Szint (imax+5,j,1,l) 

 
!         end do
      end do

!$OMP END DO 

  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO  

      do i=-6,imax+6
!      do k=1,kmax

            psiint(i,0,1,l) = psiint(i, 1,1,l) 
            phiint(i,0,1,l) = phiint(i, 1,1,l) 
            Bxint (i,0,1,l) = Bxint (i, 1,1,l)  
            Byint (i,0,1,l) = Byint (i, 1,1,l)  
            Bzint (i,0,1,l) = Bzint (i, 1,1,l)  
            qint  (i,0,1,l) = qint  (i, 1,1,l)   
            DDint (i,0,1,l) = DDint (i, 1,1,l)  
            tauint(i,0,1,l) = tauint(i, 1,1,l) 
            Sxint (i,0,1,l) = Sxint (i, 1,1,l)  
            Syint (i,0,1,l) = Syint (i, 1,1,l)  
            Szint (i,0,1,l) = Szint (i, 1,1,l)  
 
            psiint(i,-1,1,l) = psiint(i, 0,1,l) 
            phiint(i,-1,1,l) = phiint(i, 0,1,l) 
            Bxint (i,-1,1,l) = Bxint (i, 0,1,l)  
            Byint (i,-1,1,l) = Byint (i, 0,1,l)  
            Bzint (i,-1,1,l) = Bzint (i, 0,1,l)  
            qint  (i,-1,1,l) = qint  (i, 0,1,l)   
            DDint (i,-1,1,l) = DDint (i, 0,1,l)  
            tauint(i,-1,1,l) = tauint(i, 0,1,l) 
            Sxint (i,-1,1,l) = Sxint (i, 0,1,l)  
            Syint (i,-1,1,l) = Syint (i, 0,1,l)  
            Szint (i,-1,1,l) = Szint (i, 0,1,l)  

            psiint(i,-2,1,l) = psiint(i,-1,1,l) 
            phiint(i,-2,1,l) = phiint(i,-1,1,l) 
            Bxint (i,-2,1,l) = Bxint (i,-1,1,l)  
            Byint (i,-2,1,l) = Byint (i,-1,1,l)  
            Bzint (i,-2,1,l) = Bzint (i,-1,1,l)  
            qint  (i,-2,1,l) = qint  (i,-1,1,l)   
            DDint (i,-2,1,l) = DDint (i,-1,1,l)  
            tauint(i,-2,1,l) = tauint(i,-1,1,l) 
            Sxint (i,-2,1,l) = Sxint (i,-1,1,l)  
            Syint (i,-2,1,l) = Syint (i,-1,1,l)  
            Szint (i,-2,1,l) = Szint (i,-1,1,l)  

            psiint(i,-3,1,l) = psiint(i,-2,1,l) 
            phiint(i,-3,1,l) = phiint(i,-2,1,l) 
            Bxint (i,-3,1,l) = Bxint (i,-2,1,l)  
            Byint (i,-3,1,l) = Byint (i,-2,1,l)  
            Bzint (i,-3,1,l) = Bzint (i,-2,1,l)  
            qint  (i,-3,1,l) = qint  (i,-2,1,l)   
            DDint (i,-3,1,l) = DDint (i,-2,1,l)  
            tauint(i,-3,1,l) = tauint(i,-2,1,l) 
            Sxint (i,-3,1,l) = Sxint (i,-2,1,l)  
            Syint (i,-3,1,l) = Syint (i,-2,1,l)  
            Szint (i,-3,1,l) = Szint (i,-2,1,l)  

            psiint(i,-4,1,l) = psiint(i,-3,1,l) 
            phiint(i,-4,1,l) = phiint(i,-3,1,l) 
            Bxint (i,-4,1,l) = Bxint (i,-3,1,l)  
            Byint (i,-4,1,l) = Byint (i,-3,1,l)  
            Bzint (i,-4,1,l) = Bzint (i,-3,1,l)  
            qint  (i,-4,1,l) = qint  (i,-3,1,l)   
            DDint (i,-4,1,l) = DDint (i,-3,1,l)  
            tauint(i,-4,1,l) = tauint(i,-3,1,l) 
            Sxint (i,-4,1,l) = Sxint (i,-3,1,l)  
            Syint (i,-4,1,l) = Syint (i,-3,1,l)  
            Szint (i,-4,1,l) = Szint (i,-3,1,l)

            psiint(i,-5,1,l) = psiint(i,-4,1,l) 
            phiint(i,-5,1,l) = phiint(i,-4,1,l) 
            Bxint (i,-5,1,l) = Bxint (i,-4,1,l)  
            Byint (i,-5,1,l) = Byint (i,-4,1,l)  
            Bzint (i,-5,1,l) = Bzint (i,-4,1,l)  
            qint  (i,-5,1,l) = qint  (i,-4,1,l)   
            DDint (i,-5,1,l) = DDint (i,-4,1,l)  
            tauint(i,-5,1,l) = tauint(i,-4,1,l) 
            Sxint (i,-5,1,l) = Sxint (i,-4,1,l)  
            Syint (i,-5,1,l) = Syint (i,-4,1,l)  
            Szint (i,-5,1,l) = Szint (i,-4,1,l)

            
            psiint(i,-6,1,l) = psiint(i,-5,1,l) 
            phiint(i,-6,1,l) = phiint(i,-5,1,l) 
            Bxint (i,-6,1,l) = Bxint (i,-5,1,l)  
            Byint (i,-6,1,l) = Byint (i,-5,1,l)  
            Bzint (i,-6,1,l) = Bzint (i,-5,1,l)  
            qint  (i,-6,1,l) = qint  (i,-5,1,l)   
            DDint (i,-6,1,l) = DDint (i,-5,1,l)  
            tauint(i,-6,1,l) = tauint(i,-5,1,l) 
            Sxint (i,-6,1,l) = Sxint (i,-5,1,l)  
            Syint (i,-6,1,l) = Syint (i,-5,1,l)  
            Szint (i,-6,1,l) = Szint (i,-5,1,l)  
 
            ! Right Boundary
            !---------------
            
            psiint(i,jmax  ,1,l) = psiint(i,jmax-1,1,l) 
            phiint(i,jmax  ,1,l) = phiint(i,jmax-1,1,l) 
            Bxint (i,jmax  ,1,l) = Bxint (i,jmax-1,1,l)  
            Byint (i,jmax  ,1,l) = Byint (i,jmax-1,1,l)  
            Bzint (i,jmax  ,1,l) = Bzint (i,jmax-1,1,l)  
            qint  (i,jmax  ,1,l) = qint  (i,jmax-1,1,l)   
            DDint (i,jmax  ,1,l) = DDint (i,jmax-1,1,l)  
            tauint(i,jmax  ,1,l) = tauint(i,jmax-1,1,l) 
            Sxint (i,jmax  ,1,l) = Sxint (i,jmax-1,1,l)  
            Syint (i,jmax  ,1,l) = Syint (i,jmax-1,1,l)  
            Szint (i,jmax  ,1,l) = Szint (i,jmax-1,1,l)  

            psiint(i,jmax+1,1,l) = psiint(i,jmax  ,1,l) 
            phiint(i,jmax+1,1,l) = phiint(i,jmax  ,1,l) 
            Bxint (i,jmax+1,1,l) = Bxint (i,jmax  ,1,l)  
            Byint (i,jmax+1,1,l) = Byint (i,jmax  ,1,l)  
            Bzint (i,jmax+1,1,l) = Bzint (i,jmax  ,1,l)  
            qint  (i,jmax+1,1,l) = qint  (i,jmax  ,1,l)   
            DDint (i,jmax+1,1,l) = DDint (i,jmax  ,1,l)  
            tauint(i,jmax+1,1,l) = tauint(i,jmax  ,1,l) 
            Sxint (i,jmax+1,1,l) = Sxint (i,jmax  ,1,l)  
            Syint (i,jmax+1,1,l) = Syint (i,jmax  ,1,l)  
            Szint (i,jmax+1,1,l) = Szint (i,jmax  ,1,l)  

            psiint(i,jmax+2,1,l) = psiint(i,jmax+1,1,l) 
            phiint(i,jmax+2,1,l) = phiint(i,jmax+1,1,l) 
            Bxint (i,jmax+2,1,l) = Bxint (i,jmax+1,1,l)  
            Byint (i,jmax+2,1,l) = Byint (i,jmax+1,1,l)  
            Bzint (i,jmax+2,1,l) = Bzint (i,jmax+1,1,l)  
            qint  (i,jmax+2,1,l) = qint  (i,jmax+1,1,l)   
            DDint (i,jmax+2,1,l) = DDint (i,jmax+1,1,l)  
            tauint(i,jmax+2,1,l) = tauint(i,jmax+1,1,l) 
            Sxint (i,jmax+2,1,l) = Sxint (i,jmax+1,1,l)  
            Syint (i,jmax+2,1,l) = Syint (i,jmax+1,1,l)  
            Szint (i,jmax+2,1,l) = Szint (i,jmax+1,1,l) 

            psiint(i,jmax+3,1,l) = psiint(i,jmax+2,1,l) 
            phiint(i,jmax+3,1,l) = phiint(i,jmax+2,1,l) 
            Bxint (i,jmax+3,1,l) = Bxint (i,jmax+2,1,l)  
            Byint (i,jmax+3,1,l) = Byint (i,jmax+2,1,l)  
            Bzint (i,jmax+3,1,l) = Bzint (i,jmax+2,1,l)  
            qint  (i,jmax+3,1,l) = qint  (i,jmax+2,1,l)   
            DDint (i,jmax+3,1,l) = DDint (i,jmax+2,1,l)  
            tauint(i,jmax+3,1,l) = tauint(i,jmax+2,1,l) 
            Sxint (i,jmax+3,1,l) = Sxint (i,jmax+2,1,l)  
            Syint (i,jmax+3,1,l) = Syint (i,jmax+2,1,l)  
            Szint (i,jmax+3,1,l) = Szint (i,jmax+2,1,l) 

            psiint(i,jmax+4,1,l) = psiint(i,jmax+3,1,l) 
            phiint(i,jmax+4,1,l) = phiint(i,jmax+3,1,l) 
            Bxint (i,jmax+4,1,l) = Bxint (i,jmax+3,1,l)  
            Byint (i,jmax+4,1,l) = Byint (i,jmax+3,1,l)  
            Bzint (i,jmax+4,1,l) = Bzint (i,jmax+3,1,l)  
            qint  (i,jmax+4,1,l) = qint  (i,jmax+3,1,l)   
            DDint (i,jmax+4,1,l) = DDint (i,jmax+3,1,l)  
            tauint(i,jmax+4,1,l) = tauint(i,jmax+3,1,l) 
            Sxint (i,jmax+4,1,l) = Sxint (i,jmax+3,1,l)  
            Syint (i,jmax+4,1,l) = Syint (i,jmax+3,1,l)  
            Szint (i,jmax+4,1,l) = Szint (i,jmax+3,1,l)

            psiint(i,jmax+5,1,l) = psiint(i,jmax+4,1,l) 
            phiint(i,jmax+5,1,l) = phiint(i,jmax+4,1,l) 
            Bxint (i,jmax+5,1,l) = Bxint (i,jmax+4,1,l)  
            Byint (i,jmax+5,1,l) = Byint (i,jmax+4,1,l)  
            Bzint (i,jmax+5,1,l) = Bzint (i,jmax+4,1,l)  
            qint  (i,jmax+5,1,l) = qint  (i,jmax+4,1,l)   
            DDint (i,jmax+5,1,l) = DDint (i,jmax+4,1,l)  
            tauint(i,jmax+5,1,l) = tauint(i,jmax+4,1,l) 
            Sxint (i,jmax+5,1,l) = Sxint (i,jmax+4,1,l)  
            Syint (i,jmax+5,1,l) = Syint (i,jmax+4,1,l)  
            Szint (i,jmax+5,1,l) = Szint (i,jmax+4,1,l)

            psiint(i,jmax+6,1,l) = psiint(i,jmax+5,1,l) 
            phiint(i,jmax+6,1,l) = phiint(i,jmax+5,1,l) 
            Bxint (i,jmax+6,1,l) = Bxint (i,jmax+5,1,l)  
            Byint (i,jmax+6,1,l) = Byint (i,jmax+5,1,l)  
            Bzint (i,jmax+6,1,l) = Bzint (i,jmax+5,1,l)  
            qint  (i,jmax+6,1,l) = qint  (i,jmax+5,1,l)   
            DDint (i,jmax+6,1,l) = DDint (i,jmax+5,1,l)  
            tauint(i,jmax+6,1,l) = tauint(i,jmax+5,1,l) 
            Sxint (i,jmax+6,1,l) = Sxint (i,jmax+5,1,l)  
            Syint (i,jmax+6,1,l) = Syint (i,jmax+5,1,l)  
            Szint (i,jmax+6,1,l) = Szint (i,jmax+5,1,l) 

 
 
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

            psiint(0,j,1,l) = 0.d0
            phiint(0,j,1,l) = 0.d0
            Bxint(0,j,1,l)  = 1.d0
            Byint(0,j,1,l)  = 0.d0
            Bzint(0,j,1,l)  = 0.d0
            qint(0,j,1,l)   = 0.d0
            DDint(0,j,1,l)  = DDint(1,j,1,l) !D(0,j,1)
            tauint(0,j,1,l) = tauint(1,j,1,l)!tau(0,j,1)
            Sxint(0,j,1,l)  = Sxint(1,j,1,l) !Sx(0,j,1) 
            Syint(0,j,1,l)  = Syint(1,j,1,l) !Sy(0,j,1) 
            Szint(0,j,1,l)  = Szint(1,j,1,l) !Sz(0,j,1) 


            psiint(-1,j,1,l) = 0.d0
            phiint(-1,j,1,l) = 0.d0
            Bxint(-1,j,1,l)  = 1.d0
            Byint(-1,j,1,l)  = 0.d0
            Bzint(-1,j,1,l)  = 0.d0
            qint(-1,j,1,l)   = 0.d0
            DDint(-1,j,1,l)  = DDint(1,j,1,l) !D(-1,j,1)
            tauint(-1,j,1,l) = tauint(1,j,1,l)!tau(-1,j,1)
            Sxint(-1,j,1,l)  = Sxint(1,j,1,l) !Sx(-1,j,1) 
            Syint(-1,j,1,l)  = Syint(1,j,1,l) !Sy(-1,j,1) 
            Szint(-1,j,1,l)  = Szint(1,j,1,l) !Sz(-1,j,1) 


            psiint(-2,j,1,l) = 0.d0
            phiint(-2,j,1,l) = 0.d0
            Bxint(-2,j,1,l)  = 1.d0
            Byint(-2,j,1,l)  = 0.d0
            Bzint(-2,j,1,l)  = 0.d0
            qint(-2,j,1,l)   = 0.d0
            DDint(-2,j,1,l)  = DDint(1,j,1,l) !D(-2,j,1)
            tauint(-2,j,1,l) = tauint(1,j,1,l)!tau(-2,j,1)
            Sxint(-2,j,1,l)  = Sxint(1,j,1,l) !Sx(-2,j,1) 
            Syint(-2,j,1,l)  = Syint(1,j,1,l) !Sy(-2,j,1) 
            Szint(-2,j,1,l)  = Szint(1,j,1,l) !Sz(-2,j,1) 

!           print*, DDint(-2,j,1,l) , i, j, h

           else 

            psiint(0,j,1,l) =  psiint(1,j,1,l)
            phiint(0,j,1,l) =  phiint(1,j,1,l)
            Bxint(0,j,1,l)  =  Bxint(1,j,1,l) 
            Byint(0,j,1,l)  = -Byint(1,j,1,l) 
            Bzint(0,j,1,l)  = -Bzint(1,j,1,l) 
            qint(0,j,1,l)   =  qint(1,j,1,l) 
            DDint(0,j,1,l)  =  DDint(1,j,1,l)
            tauint(0,j,1,l) =  tauint(1,j,1,l)
            Sxint(0,j,1,l)  =  Sxint(1,j,1,l) 
            Syint(0,j,1,l)  = -Syint(1,j,1,l) 
            Szint(0,j,1,l)  = -Szint(1,j,1,l) 


            psiint(-1,j,1,l) = psiint(0,j,1,l)
            phiint(-1,j,1,l) = phiint(0,j,1,l)
            Bxint(-1,j,1,l)  = Bxint(0,j,1,l) 
            Byint(-1,j,1,l)  = Byint(0,j,1,l) 
            Bzint(-1,j,1,l)  = Bzint(0,j,1,l) 
            qint(-1,j,1,l)   = qint(0,j,1,l) 
            DDint(-1,j,1,l)  = DDint(0,j,1,l)
            tauint(-1,j,1,l) = tauint(0,j,1,l)
            Sxint(-1,j,1,l)  = Sxint(0,j,1,l) 
            Syint(-1,j,1,l)  = Syint(0,j,1,l) 
            Szint(-1,j,1,l)  = Szint(0,j,1,l) 

            psiint(-2,j,1,l) = psiint(-1,j,1,l)
            phiint(-2,j,1,l) = phiint(-1,j,1,l)
            Bxint(-2,j,1,l)  = Bxint(-1,j,1,l) 
            Byint(-2,j,1,l)  = Byint(-1,j,1,l) 
            Bzint(-2,j,1,l)  = Bzint(-1,j,1,l) 
            qint(-2,j,1,l)   = qint(-1,j,1,l) 
            DDint(-2,j,1,l)  = DDint(-1,j,1,l)
            tauint(-2,j,1,l) = tauint(-1,j,1,l)
            Sxint(-2,j,1,l)  = Sxint(-1,j,1,l) 
            Syint(-2,j,1,l)  = Syint(-1,j,1,l) 
            Szint(-2,j,1,l)  = Szint(-1,j,1,l)

            end if


         ! Right Boundary
         !---------------

            psiint(imax,j,1,l) = psiint(imax-1,j,1,l) 
            phiint(imax,j,1,l) = phiint(imax-1,j,1,l) 
            Bxint(imax,j,1,l)  = Bxint(imax-1,j,1,l)  
            Byint(imax,j,1,l)  = Byint(imax-1,j,1,l)  
            Bzint(imax,j,1,l)  = Bzint(imax-1,j,1,l)  
            qint(imax,j,1,l)   = qint(imax-1,j,1,l)  
            DDint(imax,j,1,l)  = DDint(imax-1,j,1,l) 
            tauint(imax,j,1,l) = tauint(imax-1,j,1,l)
            Sxint(imax,j,1,l)  = Sxint(imax-1,j,1,l) 
            Syint(imax,j,1,l)  = Syint(imax-1,j,1,l) 
            Szint(imax,j,1,l)  = Szint(imax-1,j,1,l) 

            psiint(imax+1,j,1,l) = psiint(imax,j,1,l) 
            phiint(imax+1,j,1,l) = phiint(imax,j,1,l) 
            Bxint(imax+1,j,1,l)  = Bxint(imax,j,1,l)  
            Byint(imax+1,j,1,l)  = Byint(imax,j,1,l)  
            Bzint(imax+1,j,1,l)  = Bzint(imax,j,1,l)  
            qint(imax+1,j,1,l)   = qint(imax,j,1,l)  
            DDint(imax+1,j,1,l)  = DDint(imax,j,1,l) 
            tauint(imax+1,j,1,l) = tauint(imax,j,1,l)
            Sxint(imax+1,j,1,l)  = Sxint(imax,j,1,l) 
            Syint(imax+1,j,1,l)  = Syint(imax,j,1,l) 
            Szint(imax+1,j,1,l)  = Szint(imax,j,1,l) 

            psiint(imax+2,j,1,l) = psiint(imax+1,j,1,l)  
            phiint(imax+2,j,1,l) = phiint(imax+1,j,1,l)  
            Bxint(imax+2,j,1,l)  = Bxint(imax+1,j,1,l)  
            Byint(imax+2,j,1,l)  = Byint(imax+1,j,1,l)  
            Bzint(imax+2,j,1,l)  = Bzint(imax+1,j,1,l)  
            qint(imax+2,j,1,l)   = qint(imax+1,j,1,l)  
            DDint(imax+2,j,1,l)  = DDint(imax+1,j,1,l) 
            tauint(imax+2,j,1,l) = tauint(imax+1,j,1,l)
            Sxint(imax+2,j,1,l)  = Sxint(imax+1,j,1,l) 
            Syint(imax+2,j,1,l)  = Syint(imax+1,j,1,l) 
            Szint(imax+2,j,1,l)  = Szint(imax+1,j,1,l) 

 
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

            psiint(i,0,1,l) = psiint(i,1,1,l) 
            phiint(i,0,1,l) = phiint(i,1,1,l) 
            Bxint (i,0,1,l) = Bxint (i,1,1,l)  
            Byint (i,0,1,l) = Byint (i,1,1,l)  
            Bzint (i,0,1,l) = Bzint (i,1,1,l)  
            qint  (i,0,1,l) = qint  (i,1,1,l)   
            DDint (i,0,1,l) = DDint (i,1,1,l)  
            tauint(i,0,1,l) = tauint(i,1,1,l) 
            Sxint (i,0,1,l) = Sxint (i,1,1,l)  
            Syint (i,0,1,l) = Syint (i,1,1,l)  
            Szint (i,0,1,l) = Szint (i,1,1,l)  


            psiint(i,-1,1,l) = psiint(i,0,1,l) 
            phiint(i,-1,1,l) = phiint(i,0,1,l) 
            Bxint (i,-1,1,l) = Bxint (i,0,1,l)  
            Byint (i,-1,1,l) = Byint (i,0,1,l)  
            Bzint (i,-1,1,l) = Bzint (i,0,1,l)  
            qint  (i,-1,1,l) = qint  (i,0,1,l)   
            DDint (i,-1,1,l) = DDint (i,0,1,l)  
            tauint(i,-1,1,l) = tauint(i,0,1,l) 
            Sxint (i,-1,1,l) = Sxint (i,0,1,l)  
            Syint (i,-1,1,l) = Syint (i,0,1,l)  
            Szint (i,-1,1,l) = Szint (i,0,1,l)  


            psiint(i,-2,1,l) = psiint(i,-1,1,l) 
            phiint(i,-2,1,l) = phiint(i,-1,1,l) 
            Bxint (i,-2,1,l) = Bxint (i,-1,1,l)  
            Byint (i,-2,1,l) = Byint (i,-1,1,l)  
            Bzint (i,-2,1,l) = Bzint (i,-1,1,l)  
            qint  (i,-2,1,l) = qint  (i,-1,1,l)   
            DDint (i,-2,1,l) = DDint (i,-1,1,l)  
            tauint(i,-2,1,l) = tauint(i,-1,1,l) 
            Sxint (i,-2,1,l) = Sxint (i,-1,1,l)  
            Syint (i,-2,1,l) = Syint (i,-1,1,l)  
            Szint (i,-2,1,l) = Szint (i,-1,1,l)  


            ! Right Boundary
            !---------------

            psiint(i,jmax,1,l) = psiint(i,jmax-1,1,l) 
            phiint(i,jmax,1,l) = phiint(i,jmax-1,1,l) 
            Bxint (i,jmax,1,l) = Bxint(i,jmax-1,1,l)  
            Byint (i,jmax,1,l) = Byint(i,jmax-1,1,l)  
            Bzint (i,jmax,1,l) = Bzint(i,jmax-1,1,l)  
            qint  (i,jmax,1,l) = qint(i,jmax-1,1,l)   
            DDint (i,jmax,1,l) = DDint(i,jmax-1,1,l)  
            tauint(i,jmax,1,l) = tauint(i,jmax-1,1,l) 
            Sxint (i,jmax,1,l) = Sxint(i,jmax-1,1,l)  
            Syint (i,jmax,1,l) = Syint(i,jmax-1,1,l)  
            Szint (i,jmax,1,l) = Szint(i,jmax-1,1,l)  


            psiint(i,jmax+1,1,l) = psiint(i,jmax,1,l) 
            phiint(i,jmax+1,1,l) = phiint(i,jmax,1,l) 
            Bxint (i,jmax+1,1,l) = Bxint(i,jmax,1,l)  
            Byint (i,jmax+1,1,l) = Byint(i,jmax,1,l)  
            Bzint (i,jmax+1,1,l) = Bzint(i,jmax,1,l)  
            qint  (i,jmax+1,1,l) = qint(i,jmax,1,l)   
            DDint (i,jmax+1,1,l) = DDint(i,jmax,1,l)  
            tauint(i,jmax+1,1,l) = tauint(i,jmax,1,l) 
            Sxint (i,jmax+1,1,l) = Sxint(i,jmax,1,l)  
            Syint (i,jmax+1,1,l) = Syint(i,jmax,1,l)  
            Szint (i,jmax+1,1,l) = Szint(i,jmax,1,l)  


            psiint(i,jmax+2,1,l) = psiint(i,jmax+1,1,l) 
            phiint(i,jmax+2,1,l) = phiint(i,jmax+1,1,l) 
            Bxint (i,jmax+2,1,l) = Bxint(i,jmax+1,1,l)  
            Byint (i,jmax+2,1,l) = Byint(i,jmax+1,1,l)  
            Bzint (i,jmax+2,1,l) = Bzint(i,jmax+1,1,l)  
            qint  (i,jmax+2,1,l) = qint(i,jmax+1,1,l)   
            DDint (i,jmax+2,1,l) = DDint(i,jmax+1,1,l)  
            tauint(i,jmax+2,1,l) = tauint(i,jmax+1,1,l) 
            Sxint (i,jmax+2,1,l) = Sxint(i,jmax+1,1,l)  
            Syint (i,jmax+2,1,l) = Syint(i,jmax+1,1,l)  
            Szint (i,jmax+2,1,l) = Szint(i,jmax+1,1,l)  

 
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

            psiint(0,j,1,l) = psiint(1,j,1,l)
            phiint(0,j,1,l) = phiint(1,j,1,l)
            Bxint(0,j,1,l)  = Bxint(1,j,1,l)
            Byint(0,j,1,l)  = Byint(1,j,1,l) 
            Bzint(0,j,1,l)  = Bzint(1,j,1,l)
            qint(0,j,1,l)   = qint(1,j,1,l) 
            DDint(0,j,1,l)  = DDint(1,j,1,l)
            tauint(0,j,1,l) = tauint(1,j,1,l)
            Sxint(0,j,1,l)  = Sxint(1,j,1,l) 
            Syint(0,j,1,l)  = Syint(1,j,1,l)
            Szint(0,j,1,l)  = Szint(1,j,1,l) 

            psiint(-1,j,1,l) = psiint(0,j,1,l)
            phiint(-1,j,1,l) = phiint(0,j,1,l)
            Bxint(-1,j,1,l)  = Bxint(0,j,1,l) 
            Byint(-1,j,1,l)  = Byint(0,j,1,l) 
            Bzint(-1,j,1,l)  = Bzint(0,j,1,l) 
            qint(-1,j,1,l)   = qint(0,j,1,l) 
            DDint(-1,j,1,l)  = DDint(0,j,1,l)
            tauint(-1,j,1,l) = tauint(0,j,1,l)
            Sxint(-1,j,1,l)  = Sxint(0,j,1,l) 
            Syint(-1,j,1,l)  = Syint(0,j,1,l) 
            Szint(-1,j,1,l)  = Szint(0,j,1,l) 

            psiint(-2,j,1,l) = psiint(-1,j,1,l)
            phiint(-2,j,1,l) = phiint(-1,j,1,l)
            Bxint(-2,j,1,l)  = Bxint(-1,j,1,l) 
            Byint(-2,j,1,l)  = Byint(-1,j,1,l) 
            Bzint(-2,j,1,l)  = Bzint(-1,j,1,l) 
            qint(-2,j,1,l)   = qint(-1,j,1,l) 
            DDint(-2,j,1,l)  = DDint(-1,j,1,l)
            tauint(-2,j,1,l) = tauint(-1,j,1,l)
            Sxint(-2,j,1,l)  = Sxint(-1,j,1,l) 
            Syint(-2,j,1,l)  = Syint(-1,j,1,l) 
            Szint(-2,j,1,l)  = Szint(-1,j,1,l) 


            ! Right Boundary
            !---------------

            psiint(imax,j,1,l) = psi(imax,j,1) !psiint(imax-1,j,1,l)  
            phiint(imax,j,1,l) = phi(imax,j,1) !phiint(imax-1,j,1,l)
            Bxint(imax,j,1,l)  = Bx(imax,j,1)  !Bxint(imax-1,j,1,l)
            Byint(imax,j,1,l)  = By(imax,j,1)  !Byint(imax-1,j,1,l)
            Bzint(imax,j,1,l)  = Bz(imax,j,1)  !Bzint(imax-1,j,1,l)
            qint(imax,j,1,l)   = q(imax,j,1)   !qint(imax-1,j,1,l)
            DDint(imax,j,1,l)  = D(imax,j,1)   !DDint(imax-1,j,1,l)
            tauint(imax,j,1,l) = tau(imax,j,1) !tauint(imax-1,j,1,l)
            Sxint(imax,j,1,l)  = Sx(imax,j,1)  !Sxint(imax-1,j,1,l)
            Syint(imax,j,1,l)  = Sy(imax,j,1)  !Syint(imax-1,j,1,l)
            Szint(imax,j,1,l)  = Sz(imax,j,1)  !Szint(imax-1,j,1,l)

            psiint(imax+1,j,1,l) = psiint(imax,j,1,l) 
            phiint(imax+1,j,1,l) = phiint(imax,j,1,l) 
            Bxint(imax+1,j,1,l)  = Bxint(imax,j,1,l)  
            Byint(imax+1,j,1,l)  = Byint(imax,j,1,l)  
            Bzint(imax+1,j,1,l)  = Bzint(imax,j,1,l)  
            qint(imax+1,j,1,l)   = qint(imax,j,1,l)  
            DDint(imax+1,j,1,l)  = DDint(imax,j,1,l) 
            tauint(imax+1,j,1,l) = tauint(imax,j,1,l)
            Sxint(imax+1,j,1,l)  = Sxint(imax,j,1,l) 
            Syint(imax+1,j,1,l)  = Syint(imax,j,1,l) 
            Szint(imax+1,j,1,l)  = Szint(imax,j,1,l) 

            psiint(imax+2,j,1,l) = psiint(imax,j,1,l)  
            phiint(imax+2,j,1,l) = phiint(imax,j,1,l)  
            Bxint(imax+2,j,1,l)  = Bxint(imax,j,1,l)  
            Byint(imax+2,j,1,l)  = Byint(imax,j,1,l)  
            Bzint(imax+2,j,1,l)  = Bzint(imax,j,1,l)  
            qint(imax+2,j,1,l)   = qint(imax,j,1,l)  
            DDint(imax+2,j,1,l)  = DDint(imax,j,1,l) 
            tauint(imax+2,j,1,l) = tauint(imax,j,1,l)
            Sxint(imax+2,j,1,l)  = Sxint(imax,j,1,l) 
            Syint(imax+2,j,1,l)  = Syint(imax,j,1,l) 
            Szint(imax+2,j,1,l)  = Szint(imax,j,1,l) 

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

            psiint(i,0,1,l) = psiint(i,1,1,l) 
            phiint(i,0,1,l) = phiint(i,1,1,l) 
            Bxint (i,0,1,l) = Bxint(i,1,1,l)  
            Byint (i,0,1,l) = Byint(i,1,1,l)  
            Bzint (i,0,1,l) = Bzint(i,1,1,l)  
            qint  (i,0,1,l) = qint(i,1,1,l)   
            DDint (i,0,1,l) = DDint(i,1,1,l)  
            tauint(i,0,1,l) = tauint(i,1,1,l) 
            Sxint (i,0,1,l) = Sxint(i,1,1,l)  
            Syint (i,0,1,l) = Syint(i,1,1,l)  
            Szint (i,0,1,l) = Szint(i,1,1,l)  

            psiint(i,-1,1,l) = psiint(i,0,1,l) 
            phiint(i,-1,1,l) = phiint(i,0,1,l) 
            Bxint (i,-1,1,l) = Bxint(i,0,1,l)  
            Byint (i,-1,1,l) = Byint(i,0,1,l)  
            Bzint (i,-1,1,l) = Bzint(i,0,1,l)  
            qint  (i,-1,1,l) = qint(i,0,1,l)   
            DDint (i,-1,1,l) = DDint(i,0,1,l)  
            tauint(i,-1,1,l) = tauint(i,0,1,l) 
            Sxint (i,-1,1,l) = Sxint(i,0,1,l)  
            Syint (i,-1,1,l) = Syint(i,0,1,l)  
            Szint (i,-1,1,l) = Szint(i,0,1,l)  

            psiint(i,-2,1,l) = psiint(i,-1,1,l) 
            phiint(i,-2,1,l) = phiint(i,-1,1,l) 
            Bxint (i,-2,1,l) = Bxint(i,-1,1,l)  
            Byint (i,-2,1,l) = Byint(i,-1,1,l)  
            Bzint (i,-2,1,l) = Bzint(i,-1,1,l)  
            qint  (i,-2,1,l) = qint(i,-1,1,l)   
            DDint (i,-2,1,l) = DDint(i,-1,1,l)  
            tauint(i,-2,1,l) = tauint(i,-1,1,l) 
            Sxint (i,-2,1,l) = Sxint(i,-1,1,l)  
            Syint (i,-2,1,l) = Syint(i,-1,1,l)  
            Szint (i,-2,1,l) = Szint(i,-1,1,l)  


            ! Right Boundary
            !---------------

            psiint(i,jmax,1,l) = psiint(i,jmax-1,1,l) 
            phiint(i,jmax,1,l) = phiint(i,jmax-1,1,l) 
            Bxint (i,jmax,1,l) = Bxint(i,jmax-1,1,l)  
            Byint (i,jmax,1,l) = Byint(i,jmax-1,1,l)  
            Bzint (i,jmax,1,l) = Bzint(i,jmax-1,1,l)  
            qint  (i,jmax,1,l) = qint(i,jmax-1,1,l)   
            DDint (i,jmax,1,l) = DDint(i,jmax-1,1,l)  
            tauint(i,jmax,1,l) = tauint(i,jmax-1,1,l) 
            Sxint (i,jmax,1,l) = Sxint(i,jmax-1,1,l)  
            Syint (i,jmax,1,l) = Syint(i,jmax-1,1,l)  
            Szint (i,jmax,1,l) = Szint(i,jmax-1,1,l)  

            psiint(i,jmax+1,1,l) = psiint(i,jmax,1,l) 
            phiint(i,jmax+1,1,l) = phiint(i,jmax,1,l) 
            Bxint (i,jmax+1,1,l) = Bxint(i,jmax,1,l)  
            Byint (i,jmax+1,1,l) = Byint(i,jmax,1,l)  
            Bzint (i,jmax+1,1,l) = Bzint(i,jmax,1,l)  
            qint  (i,jmax+1,1,l) = qint(i,jmax,1,l)   
            DDint (i,jmax+1,1,l) = DDint(i,jmax,1,l)  
            tauint(i,jmax+1,1,l) = tauint(i,jmax,1,l) 
            Sxint (i,jmax+1,1,l) = Sxint(i,jmax,1,l)  
            Syint (i,jmax+1,1,l) = Syint(i,jmax,1,l)  
            Szint (i,jmax+1,1,l) = Szint(i,jmax,1,l)  


            psiint(i,jmax+2,1,l) = psiint(i,jmax+1,1,l) 
            phiint(i,jmax+2,1,l) = phiint(i,jmax+1,1,l) 
            Bxint (i,jmax+2,1,l) = Bxint(i,jmax+1,1,l)  
            Byint (i,jmax+2,1,l) = Byint(i,jmax+1,1,l)  
            Bzint (i,jmax+2,1,l) = Bzint(i,jmax+1,1,l)  
            qint  (i,jmax+2,1,l) = qint(i,jmax+1,1,l)   
            DDint (i,jmax+2,1,l) = DDint(i,jmax+1,1,l)  
            tauint(i,jmax+2,1,l) = tauint(i,jmax+1,1,l) 
            Sxint (i,jmax+2,1,l) = Sxint(i,jmax+1,1,l)  
            Syint (i,jmax+2,1,l) = Syint(i,jmax+1,1,l)  
            Szint (i,jmax+2,1,l) = Szint(i,jmax+1,1,l)  

 !         end do
      end do

 !$OMP END DO

    else if (DIM == 2 .and. BOUND == 6) then

       ! Tearing-Mode
       ! Richtmyer-Meshkov instability Boundary RM

!Tomek     

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

            psiint(i,0-j,1,l) = psiint(i,jmax-j,1,l) 
            phiint(i,0-j,1,l) = phiint(i,jmax-j,1,l) 
            Bxint (i,0-j,1,l) = Bxint (i,jmax-j,1,l)  
            Byint (i,0-j,1,l) = Byint (i,jmax-j,1,l)  
            Bzint (i,0-j,1,l) = Bzint (i,jmax-j,1,l)  
            qint  (i,0-j,1,l) = qint  (i,jmax-j,1,l)   
            DDint (i,0-j,1,l) = DDint (i,jmax-j,1,l)  
            tauint(i,0-j,1,l) = tauint(i,jmax-j,1,l) 
            Sxint (i,0-j,1,l) = Sxint (i,jmax-j,1,l)  
            Syint (i,0-j,1,l) = Syint (i,jmax-j,1,l)  
            Szint (i,0-j,1,l) = Szint (i,jmax-j,1,l)

!$OMP END ORDERED
             
          end do
       end do

!$OMP END DO

       
          ! Right  Boundary
          ! -----------------

!$OMP DO   ORDERED        

       do i=-stencil,imax+stencil
          do j = 1, stencil

!$OMP ORDERED  
             
            psiint(i,jmax+j,1,l) = psiint(i,0+j,1,l) 
            phiint(i,jmax+j,1,l) = phiint(i,0+j,1,l) 
            Bxint (i,jmax+j,1,l) = Bxint (i,0+j,1,l)  
            Byint (i,jmax+j,1,l) = Byint (i,0+j,1,l)  
            Bzint (i,jmax+j,1,l) = Bzint (i,0+j,1,l)  
            qint  (i,jmax+j,1,l) = qint  (i,0+j,1,l)   
            DDint (i,jmax+j,1,l) = DDint (i,0+j,1,l)  
            tauint(i,jmax+j,1,l) = tauint(i,0+j,1,l) 
            Sxint (i,jmax+j,1,l) = Sxint (i,0+j,1,l)  
            Syint (i,jmax+j,1,l) = Syint (i,0+j,1,l)  
            Szint (i,jmax+j,1,l) = Szint (i,0+j,1,l)

!$OMP END ORDERED             
             
          end do
       end do

!$OMP END DO       

       ! Copy boundary x-direction

  !-----------------
  ! X Direction
  !-----------------

           ! Left  Boundary
           ! -----------------

!$OMP DO   ORDERED  
       
       do i = 1,stencil
          do j=-stencil,jmax+stencil

!$OMP ORDERED               

            phiint(0-i,j,1,l) = phiint( 0,j,1,l)
            psiint(0-i,j,1,l) = psiint( 0,j,1,l)

            Bxint (0-i,j,1,l) = Bxint ( 0,j,1,l) 
            Byint (0-i,j,1,l) = Byint ( 0,j,1,l) 
            Bzint (0-i,j,1,l) = Bzint ( 0,j,1,l)
            
            qint  (0-i,j,1,l) = qint  ( 0,j,1,l) 
            DDint (0-i,j,1,l) = DDint ( 0,j,1,l)
            tauint(0-i,j,1,l) = tauint( 0,j,1,l)
            
            Sxint (0-i,j,1,l) = Sxint ( 0,j,1,l) 
            Syint (0-i,j,1,l) = Syint ( 0,j,1,l) 
            Szint (0-i,j,1,l) = Szint ( 0,j,1,l)

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

            phiint(imax+i,j,1,l) = phiint(imax  ,j,1,l)
            psiint(imax+i,j,1,l) = psiint(imax  ,j,1,l)    

            Bxint (imax+i,j,1,l) = Bxint (imax  ,j,1,l)
            Byint (imax+i,j,1,l) = Byint (imax  ,j,1,l) 
            Bzint (imax+i,j,1,l) = Bzint (imax  ,j,1,l)
            
            qint  (imax+i,j,1,l) = qint  (imax  ,j,1,l) 
            DDint (imax+i,j,1,l) = DDint (imax  ,j,1,l)
            tauint(imax+i,j,1,l) = tauint(imax  ,j,1,l)
            
            Sxint (imax+i,j,1,l) = Sxint (imax  ,j,1,l) 
            Syint (imax+i,j,1,l) = Syint (imax  ,j,1,l) 
            Szint (imax+i,j,1,l) = Szint (imax  ,j,1,l)

!$OMP END ORDERED               

         end do
     end do

!$OMP END DO 

!Tomek     
       
       
!*****************************************************************


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

            psiint(0,j,1,l) = - psiint(1,j,1,l)
            phiint(0,j,1,l) = - phiint(1,j,1,l)
            Bxint(0,j,1,l)  = - Bxint (1,j,1,l) 
            Byint(0,j,1,l)  = - Byint (1,j,1,l) 
            Bzint(0,j,1,l)  = - Bzint (1,j,1,l) 
            qint(0,j,1,l)   = - qint  (1,j,1,l) 
            DDint(0,j,1,l)  = - DDint (1,j,1,l)
            tauint(0,j,1,l) = - tauint(1,j,1,l)
            Sxint(0,j,1,l)  = - Sxint (1,j,1,l) 
            Syint(0,j,1,l)  = - Syint (1,j,1,l) 
            Szint(0,j,1,l)  = - Szint (1,j,1,l) 

            psiint(-1,j,1,l) = psiint(0,j,1,l)
            phiint(-1,j,1,l) = phiint(0,j,1,l)
            Bxint(-1,j,1,l)  = Bxint (0,j,1,l) 
            Byint(-1,j,1,l)  = Byint (0,j,1,l) 
            Bzint(-1,j,1,l)  = Bzint (0,j,1,l) 
            qint(-1,j,1,l)   = qint  (0,j,1,l) 
            DDint(-1,j,1,l)  = DDint (0,j,1,l)
            tauint(-1,j,1,l) = tauint(0,j,1,l)
            Sxint(-1,j,1,l)  = Sxint (0,j,1,l) 
            Syint(-1,j,1,l)  = Syint (0,j,1,l) 
            Szint(-1,j,1,l)  = Szint (0,j,1,l) 

            psiint(-2,j,1,l) = psiint(-1,j,1,l)
            phiint(-2,j,1,l) = phiint(-1,j,1,l)
            Bxint(-2,j,1,l)  = Bxint (-1,j,1,l) 
            Byint(-2,j,1,l)  = Byint (-1,j,1,l) 
            Bzint(-2,j,1,l)  = Bzint (-1,j,1,l) 
            qint(-2,j,1,l)   = qint  (-1,j,1,l) 
            DDint(-2,j,1,l)  = DDint (-1,j,1,l)
            tauint(-2,j,1,l) = tauint(-1,j,1,l)
            Sxint(-2,j,1,l)  = Sxint (-1,j,1,l) 
            Syint(-2,j,1,l)  = Syint (-1,j,1,l) 
            Szint(-2,j,1,l)  = Szint (-1,j,1,l) 


            ! Right Boundary
            !---------------

            psiint(imax,j,1,l) = psiint(imax-1,j,1,l) 
            phiint(imax,j,1,l) = phiint(imax-1,j,1,l) 
            Bxint(imax,j,1,l)  = Bxint (imax-1,j,1,l)  
            Byint(imax,j,1,l)  = Byint (imax-1,j,1,l)  
            Bzint(imax,j,1,l)  = Bzint (imax-1,j,1,l)  
            qint(imax,j,1,l)   = qint  (imax-1,j,1,l)  
            DDint(imax,j,1,l)  = DDint (imax-1,j,1,l) 
            tauint(imax,j,1,l) = tauint(imax-1,j,1,l)
            Sxint(imax,j,1,l)  = Sxint (imax-1,j,1,l) 
            Syint(imax,j,1,l)  = Syint (imax-1,j,1,l) 
            Szint(imax,j,1,l)  = Szint (imax-1,j,1,l) 

            psiint(imax+1,j,1,l) = psiint(imax,j,1,l) 
            phiint(imax+1,j,1,l) = phiint(imax,j,1,l) 
            Bxint(imax+1,j,1,l)  = Bxint (imax,j,1,l)  
            Byint(imax+1,j,1,l)  = Byint (imax,j,1,l)  
            Bzint(imax+1,j,1,l)  = Bzint (imax,j,1,l)  
            qint(imax+1,j,1,l)   = qint  (imax,j,1,l)  
            DDint(imax+1,j,1,l)  = DDint (imax,j,1,l) 
            tauint(imax+1,j,1,l) = tauint(imax,j,1,l)
            Sxint(imax+1,j,1,l)  = Sxint (imax,j,1,l) 
            Syint(imax+1,j,1,l)  = Syint (imax,j,1,l) 
            Szint(imax+1,j,1,l)  = Szint (imax,j,1,l) 

            psiint(imax+2,j,1,l) = psiint(imax+1,j,1,l)  
            phiint(imax+2,j,1,l) = phiint(imax+1,j,1,l)  
            Bxint(imax+2,j,1,l)  = Bxint (imax+1,j,1,l)  
            Byint(imax+2,j,1,l)  = Byint (imax+1,j,1,l)  
            Bzint(imax+2,j,1,l)  = Bzint (imax+1,j,1,l)  
            qint(imax+2,j,1,l)   = qint  (imax+1,j,1,l)  
            DDint(imax+2,j,1,l)  = DDint (imax+1,j,1,l) 
            tauint(imax+2,j,1,l) = tauint(imax+1,j,1,l)
            Sxint(imax+2,j,1,l)  = Sxint (imax+1,j,1,l) 
            Syint(imax+2,j,1,l)  = Syint (imax+1,j,1,l) 
            Szint(imax+2,j,1,l)  = Szint (imax+1,j,1,l) 

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

            psiint(i,0,1,l) = psiint(i,1,1,l) 
            phiint(i,0,1,l) = phiint(i,1,1,l) 
            Bxint (i,0,1,l) = Bxint (i,1,1,l)  
            Byint (i,0,1,l) = Byint (i,1,1,l)  
            Bzint (i,0,1,l) = Bzint (i,1,1,l)  
            qint  (i,0,1,l) = qint  (i,1,1,l)   
            DDint (i,0,1,l) = DDint (i,1,1,l)  
            tauint(i,0,1,l) = tauint(i,1,1,l) 
            Sxint (i,0,1,l) = Sxint (i,1,1,l)  
            Syint (i,0,1,l) = Syint (i,1,1,l)  
            Szint (i,0,1,l) = Szint (i,1,1,l)  

            psiint(i,-1,1,l) = psiint(i,0,1,l) 
            phiint(i,-1,1,l) = phiint(i,0,1,l) 
            Bxint (i,-1,1,l) = Bxint (i,0,1,l)  
            Byint (i,-1,1,l) = Byint (i,0,1,l)  
            Bzint (i,-1,1,l) = Bzint (i,0,1,l)  
            qint  (i,-1,1,l) = qint  (i,0,1,l)   
            DDint (i,-1,1,l) = DDint (i,0,1,l)  
            tauint(i,-1,1,l) = tauint(i,0,1,l) 
            Sxint (i,-1,1,l) = Sxint (i,0,1,l)  
            Syint (i,-1,1,l) = Syint (i,0,1,l)  
            Szint (i,-1,1,l) = Szint (i,0,1,l)  

            psiint(i,-2,1,l) = psiint(i,-1,1,l) 
            phiint(i,-2,1,l) = phiint(i,-1,1,l) 
            Bxint (i,-2,1,l) = Bxint (i,-1,1,l)  
            Byint (i,-2,1,l) = Byint (i,-1,1,l)  
            Bzint (i,-2,1,l) = Bzint (i,-1,1,l)  
            qint  (i,-2,1,l) = qint  (i,-1,1,l)   
            DDint (i,-2,1,l) = DDint (i,-1,1,l)  
            tauint(i,-2,1,l) = tauint(i,-1,1,l) 
            Sxint (i,-2,1,l) = Sxint (i,-1,1,l)  
            Syint (i,-2,1,l) = Syint (i,-1,1,l)  
            Szint (i,-2,1,l) = Szint (i,-1,1,l)  


            ! Right Boundary
            !---------------

            psiint(i,jmax,1,l) = psiint(i,jmax-1,1,l) 
            phiint(i,jmax,1,l) = phiint(i,jmax-1,1,l) 
            Bxint (i,jmax,1,l) = Bxint (i,jmax-1,1,l)  
            Byint (i,jmax,1,l) = Byint (i,jmax-1,1,l)  
            Bzint (i,jmax,1,l) = Bzint (i,jmax-1,1,l)  
            qint  (i,jmax,1,l) = qint  (i,jmax-1,1,l)   
            DDint (i,jmax,1,l) = DDint (i,jmax-1,1,l)  
            tauint(i,jmax,1,l) = tauint(i,jmax-1,1,l) 
            Sxint (i,jmax,1,l) = Sxint (i,jmax-1,1,l)  
            Syint (i,jmax,1,l) = Syint (i,jmax-1,1,l)  
            Szint (i,jmax,1,l) = Szint (i,jmax-1,1,l)  

            psiint(i,jmax+1,1,l) = psiint(i,jmax,1,l) 
            phiint(i,jmax+1,1,l) = phiint(i,jmax,1,l) 
            Bxint (i,jmax+1,1,l) = Bxint (i,jmax,1,l)  
            Byint (i,jmax+1,1,l) = Byint (i,jmax,1,l)  
            Bzint (i,jmax+1,1,l) = Bzint (i,jmax,1,l)  
            qint  (i,jmax+1,1,l) = qint  (i,jmax,1,l)   
            DDint (i,jmax+1,1,l) = DDint (i,jmax,1,l)  
            tauint(i,jmax+1,1,l) = tauint(i,jmax,1,l) 
            Sxint (i,jmax+1,1,l) = Sxint (i,jmax,1,l)  
            Syint (i,jmax+1,1,l) = Syint (i,jmax,1,l)  
            Szint (i,jmax+1,1,l) = Szint (i,jmax,1,l)  


            psiint(i,jmax+2,1,l) = psiint(i,jmax+1,1,l) 
            phiint(i,jmax+2,1,l) = phiint(i,jmax+1,1,l) 
            Bxint (i,jmax+2,1,l) = Bxint (i,jmax+1,1,l)  
            Byint (i,jmax+2,1,l) = Byint (i,jmax+1,1,l)  
            Bzint (i,jmax+2,1,l) = Bzint (i,jmax+1,1,l)  
            qint  (i,jmax+2,1,l) = qint  (i,jmax+1,1,l)   
            DDint (i,jmax+2,1,l) = DDint (i,jmax+1,1,l)  
            tauint(i,jmax+2,1,l) = tauint(i,jmax+1,1,l) 
            Sxint (i,jmax+2,1,l) = Sxint (i,jmax+1,1,l)  
            Syint (i,jmax+2,1,l) = Syint (i,jmax+1,1,l)  
            Szint (i,jmax+2,1,l) = Szint (i,jmax+1,1,l) 


!         end do
      end do
 
!$OMP END DO

    else

       write(*,*) "STOP: subroutine boundary_conserved"
       write(*,*) "This combination of Dimension and Boundary is not correctly"
       write(*,*) "DIM =", DIM, "BOUND =", BOUND
       stop

    end if
   

!$OMP END PARALLEL


  end subroutine boundary_conserved
