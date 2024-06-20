  subroutine boundary

    use scalar
    use parameters
    use threevectors 

    implicit none

 !$OMP   PARALLEL PRIVATE(posx,posy,r)
 

    if (DIM == 1 .and. BOUND == 1 ) then

       ! Copy boundary 1D


  !-----------------
  ! X Direction
  !-----------------

           ! Left Boundary
           ! -----------------

!!$           psi( 0,1,1) = psi( 1,1,1)
!!$           phi( 0,1,1) = phi( 1,1,1) 
!!$           Bx ( 0,1,1) = Bx ( 1,1,1) 
!!$           By ( 0,1,1) = By ( 1,1,1) 
!!$           Bz ( 0,1,1) = Bz ( 1,1,1)
!!$           Ex ( 0,1,1) = Ex ( 1,1,1)  
!!$           Ey ( 0,1,1) = Ey ( 1,1,1)  
!!$           Ez ( 0,1,1) = Ez ( 1,1,1)  
!!$           Vx ( 0,1,1) = Vx ( 1,1,1)  
!!$           Vy ( 0,1,1) = Vy ( 1,1,1)  
!!$           Vz ( 0,1,1) = Vz ( 1,1,1)  
!!$           q  ( 0,1,1) = q  ( 1,1,1)   
!!$           p  ( 0,1,1) = p  ( 1,1,1)   
!!$           rho( 0,1,1) = rho( 1,1,1) 

           psi(-1,1,1) = psi( 0,1,1)
           phi(-1,1,1) = phi( 0,1,1) 
           Bx (-1,1,1) = Bx ( 0,1,1) 
           By (-1,1,1) = By ( 0,1,1) 
           Bz (-1,1,1) = Bz ( 0,1,1)
           Ex (-1,1,1) = Ex ( 0,1,1)  
           Ey (-1,1,1) = Ey ( 0,1,1)  
           Ez (-1,1,1) = Ez ( 0,1,1)  
           Vx (-1,1,1) = Vx ( 0,1,1)  
           Vy (-1,1,1) = Vy ( 0,1,1)  
           Vz (-1,1,1) = Vz ( 0,1,1)  
           q  (-1,1,1) = q  ( 0,1,1)   
           p  (-1,1,1) = p  ( 0,1,1)   
           rho(-1,1,1) = rho( 0,1,1) 

           psi(-2,1,1) = psi(-1,1,1)
           phi(-2,1,1) = phi(-1,1,1) 
           Bx (-2,1,1) = Bx (-1,1,1) 
           By (-2,1,1) = By (-1,1,1) 
           Bz (-2,1,1) = Bz (-1,1,1)
           Ex (-2,1,1) = Ex (-1,1,1)  
           Ey (-2,1,1) = Ey (-1,1,1)  
           Ez (-2,1,1) = Ez (-1,1,1)  
           Vx (-2,1,1) = Vx (-1,1,1)  
           Vy (-2,1,1) = Vy (-1,1,1)  
           Vz (-2,1,1) = Vz (-1,1,1)  
           q  (-2,1,1) = q  (-1,1,1)   
           p  (-2,1,1) = p  (-1,1,1)   
           rho(-2,1,1) = rho(-1,1,1) 

           psi(-3,1,1) = psi(-2,1,1)
           phi(-3,1,1) = phi(-2,1,1) 
           Bx (-3,1,1) = Bx (-2,1,1) 
           By (-3,1,1) = By (-2,1,1) 
           Bz (-3,1,1) = Bz (-2,1,1)
           Ex (-3,1,1) = Ex (-2,1,1)  
           Ey (-3,1,1) = Ey (-2,1,1)  
           Ez (-3,1,1) = Ez (-2,1,1)  
           Vx (-3,1,1) = Vx (-2,1,1)  
           Vy (-3,1,1) = Vy (-2,1,1)  
           Vz (-3,1,1) = Vz (-2,1,1)  
           q  (-3,1,1) = q  (-2,1,1)   
           p  (-3,1,1) = p  (-2,1,1)   
           rho(-3,1,1) = rho(-2,1,1) 

           psi(-4,1,1) = psi(-3,1,1)
           phi(-4,1,1) = phi(-3,1,1) 
           Bx (-4,1,1) = Bx (-3,1,1) 
           By (-4,1,1) = By (-3,1,1) 
           Bz (-4,1,1) = Bz (-3,1,1)
           Ex (-4,1,1) = Ex (-3,1,1)  
           Ey (-4,1,1) = Ey (-3,1,1)  
           Ez (-4,1,1) = Ez (-3,1,1)  
           Vx (-4,1,1) = Vx (-3,1,1)  
           Vy (-4,1,1) = Vy (-3,1,1)  
           Vz (-4,1,1) = Vz (-3,1,1)  
           q  (-4,1,1) = q  (-3,1,1)   
           p  (-4,1,1) = p  (-3,1,1)   
           rho(-4,1,1) = rho(-3,1,1)
           
           psi(-5,1,1) = psi(-4,1,1)
           phi(-5,1,1) = phi(-4,1,1) 
           Bx (-5,1,1) = Bx (-4,1,1) 
           By (-5,1,1) = By (-4,1,1) 
           Bz (-5,1,1) = Bz (-4,1,1)
           Ex (-5,1,1) = Ex (-4,1,1)  
           Ey (-5,1,1) = Ey (-4,1,1)  
           Ez (-5,1,1) = Ez (-4,1,1)  
           Vx (-5,1,1) = Vx (-4,1,1)  
           Vy (-5,1,1) = Vy (-4,1,1)  
           Vz (-5,1,1) = Vz (-4,1,1)  
           q  (-5,1,1) = q  (-4,1,1)   
           p  (-5,1,1) = p  (-4,1,1)   
           rho(-5,1,1) = rho(-4,1,1)

           
           psi(-6,1,1) = psi(-5,1,1)
           phi(-6,1,1) = phi(-5,1,1) 
           Bx (-6,1,1) = Bx (-5,1,1) 
           By (-6,1,1) = By (-5,1,1) 
           Bz (-6,1,1) = Bz (-5,1,1)
           Ex (-6,1,1) = Ex (-5,1,1)  
           Ey (-6,1,1) = Ey (-5,1,1)  
           Ez (-6,1,1) = Ez (-5,1,1)  
           Vx (-6,1,1) = Vx (-5,1,1)  
           Vy (-6,1,1) = Vy (-5,1,1)  
           Vz (-6,1,1) = Vz (-5,1,1)  
           q  (-6,1,1) = q  (-5,1,1)   
           p  (-6,1,1) = p  (-5,1,1)   
           rho(-6,1,1) = rho(-5,1,1)

           
           ! Right  Boundary
           ! -----------------

!!$           psi(imax  ,1,1) = psi(imax-1,1,1) 
!!$           phi(imax  ,1,1) = phi(imax-1,1,1) 
!!$           Bx (imax  ,1,1) = Bx (imax-1,1,1)  
!!$           By (imax  ,1,1) = By (imax-1,1,1)  
!!$           Bz (imax  ,1,1) = Bz (imax-1,1,1)  
!!$           Ex (imax  ,1,1) = Ex (imax-1,1,1)  
!!$           Ey (imax  ,1,1) = Ey (imax-1,1,1)  
!!$           Ez (imax  ,1,1) = Ez (imax-1,1,1)  
!!$           Vx (imax  ,1,1) = Vx (imax-1,1,1)  
!!$           Vy (imax  ,1,1) = Vy (imax-1,1,1)  
!!$           Vz (imax  ,1,1) = Vz (imax-1,1,1)  
!!$           q  (imax  ,1,1) = q  (imax-1,1,1)   
!!$           p  (imax  ,1,1) = p  (imax-1,1,1)   
!!$           rho(imax  ,1,1) = rho(imax-1,1,1) 

           psi(imax+1,1,1) = psi(imax  ,1,1) 
           phi(imax+1,1,1) = phi(imax  ,1,1) 
           Bx (imax+1,1,1) = Bx (imax  ,1,1)  
           By (imax+1,1,1) = By (imax  ,1,1)  
           Bz (imax+1,1,1) = Bz (imax  ,1,1)  
           Ex (imax+1,1,1) = Ex (imax  ,1,1)  
           Ey (imax+1,1,1) = Ey (imax  ,1,1)  
           Ez (imax+1,1,1) = Ez (imax  ,1,1)  
           Vx (imax+1,1,1) = Vx (imax  ,1,1)  
           Vy (imax+1,1,1) = Vy (imax  ,1,1)  
           Vz (imax+1,1,1) = Vz (imax  ,1,1)  
           q  (imax+1,1,1) = q  (imax  ,1,1)   
           p  (imax+1,1,1) = p  (imax  ,1,1)   
           rho(imax+1,1,1) = rho(imax  ,1,1) 

           psi(imax+2,1,1) = psi(imax+1,1,1) 
           phi(imax+2,1,1) = phi(imax+1,1,1) 
           Bx (imax+2,1,1) = Bx (imax+1,1,1)  
           By (imax+2,1,1) = By (imax+1,1,1)  
           Bz (imax+2,1,1) = Bz (imax+1,1,1)  
           Ex (imax+2,1,1) = Ex (imax+1,1,1)  
           Ey (imax+2,1,1) = Ey (imax+1,1,1)  
           Ez (imax+2,1,1) = Ez (imax+1,1,1)  
           Vx (imax+2,1,1) = Vx (imax+1,1,1)  
           Vy (imax+2,1,1) = Vy (imax+1,1,1)  
           Vz (imax+2,1,1) = Vz (imax+1,1,1)  
           q  (imax+2,1,1) = q  (imax+1,1,1)   
           p  (imax+2,1,1) = p  (imax+1,1,1)   
           rho(imax+2,1,1) = rho(imax+1,1,1) 

           psi(imax+3,1,1) = psi(imax+2,1,1) 
           phi(imax+3,1,1) = phi(imax+2,1,1) 
           Bx (imax+3,1,1) = Bx (imax+2,1,1)  
           By (imax+3,1,1) = By (imax+2,1,1)  
           Bz (imax+3,1,1) = Bz (imax+2,1,1)  
           Ex (imax+3,1,1) = Ex (imax+2,1,1)  
           Ey (imax+3,1,1) = Ey (imax+2,1,1)  
           Ez (imax+3,1,1) = Ez (imax+2,1,1)  
           Vx (imax+3,1,1) = Vx (imax+2,1,1)  
           Vy (imax+3,1,1) = Vy (imax+2,1,1)  
           Vz (imax+3,1,1) = Vz (imax+2,1,1)  
           q  (imax+3,1,1) = q  (imax+2,1,1)   
           p  (imax+3,1,1) = p  (imax+2,1,1)   
           rho(imax+3,1,1) = rho(imax+2,1,1) 

           psi(imax+4,1,1) = psi(imax+3,1,1) 
           phi(imax+4,1,1) = phi(imax+3,1,1) 
           Bx (imax+4,1,1) = Bx (imax+3,1,1)  
           By (imax+4,1,1) = By (imax+3,1,1)  
           Bz (imax+4,1,1) = Bz (imax+3,1,1)  
           Ex (imax+4,1,1) = Ex (imax+3,1,1)  
           Ey (imax+4,1,1) = Ey (imax+3,1,1)  
           Ez (imax+4,1,1) = Ez (imax+3,1,1)  
           Vx (imax+4,1,1) = Vx (imax+3,1,1)  
           Vy (imax+4,1,1) = Vy (imax+3,1,1)  
           Vz (imax+4,1,1) = Vz (imax+3,1,1)  
           q  (imax+4,1,1) = q  (imax+3,1,1)   
           p  (imax+4,1,1) = p  (imax+3,1,1)   
           rho(imax+4,1,1) = rho(imax+3,1,1)

           psi(imax+5,1,1) = psi(imax+4,1,1) 
           phi(imax+5,1,1) = phi(imax+4,1,1) 
           Bx (imax+5,1,1) = Bx (imax+4,1,1)  
           By (imax+5,1,1) = By (imax+4,1,1)  
           Bz (imax+5,1,1) = Bz (imax+4,1,1)  
           Ex (imax+5,1,1) = Ex (imax+4,1,1)  
           Ey (imax+5,1,1) = Ey (imax+4,1,1)  
           Ez (imax+5,1,1) = Ez (imax+4,1,1)  
           Vx (imax+5,1,1) = Vx (imax+4,1,1)  
           Vy (imax+5,1,1) = Vy (imax+4,1,1)  
           Vz (imax+5,1,1) = Vz (imax+4,1,1)  
           q  (imax+5,1,1) = q  (imax+4,1,1)   
           p  (imax+5,1,1) = p  (imax+4,1,1)   
           rho(imax+5,1,1) = rho(imax+4,1,1)

           psi(imax+6,1,1) = psi(imax+5,1,1) 
           phi(imax+6,1,1) = phi(imax+5,1,1) 
           Bx (imax+6,1,1) = Bx (imax+5,1,1)  
           By (imax+6,1,1) = By (imax+5,1,1)  
           Bz (imax+6,1,1) = Bz (imax+5,1,1)  
           Ex (imax+6,1,1) = Ex (imax+5,1,1)  
           Ey (imax+6,1,1) = Ey (imax+5,1,1)  
           Ez (imax+6,1,1) = Ez (imax+5,1,1)  
           Vx (imax+6,1,1) = Vx (imax+5,1,1)  
           Vy (imax+6,1,1) = Vy (imax+5,1,1)  
           Vz (imax+6,1,1) = Vz (imax+5,1,1)  
           q  (imax+6,1,1) = q  (imax+5,1,1)   
           p  (imax+6,1,1) = p  (imax+5,1,1)   
           rho(imax+6,1,1) = rho(imax+5,1,1) 
          

    else if (DIM == 1 .and. BOUND == 2) then

       ! Periodic boundary

  !-----------------
  ! X Direction
  !-----------------

           ! Left Boundary
           ! -----------------
!!$           psi(0,1,1) = psi(imax-1,1,1)
!!$           phi(0,1,1) = phi(imax-1,1,1) 
!!$           Bx(0,1,1)  = Bx (imax-1,1,1) 
!!$           By(0,1,1)  = By (imax-1,1,1) 
!!$           Bz(0,1,1)  = Bz (imax-1,1,1)
!!$           Ex(0,1,1)  = Ex (imax-1,1,1)  
!!$           Ey(0,1,1)  = Ey (imax-1,1,1)  
!!$           Ez(0,1,1)  = Ez (imax-1,1,1)  
!!$           Vx(0,1,1)  = Vx (imax-1,1,1)  
!!$           Vy(0,1,1)  = Vy (imax-1,1,1)  
!!$           Vz(0,1,1)  = Vz (imax-1,1,1)  
!!$           q(0,1,1)   = q  (imax-1,1,1)   
!!$           p(0,1,1)   = p  (imax-1,1,1)   
!!$           rho(0,1,1) = rho(imax-1,1,1) 

           psi(-1,1,1) = psi(imax-1,1,1)
           phi(-1,1,1) = phi(imax-1,1,1) 
           Bx(-1,1,1)  = Bx (imax-1,1,1) 
           By(-1,1,1)  = By (imax-1,1,1) 
           Bz(-1,1,1)  = Bz (imax-1,1,1)
           Ex(-1,1,1)  = Ex (imax-1,1,1)  
           Ey(-1,1,1)  = Ey (imax-1,1,1)  
           Ez(-1,1,1)  = Ez (imax-1,1,1)  
           Vx(-1,1,1)  = Vx (imax-1,1,1)  
           Vy(-1,1,1)  = Vy (imax-1,1,1)  
           Vz(-1,1,1)  = Vz (imax-1,1,1)  
           q(-1,1,1)   = q  (imax-1,1,1)   
           p(-1,1,1)   = p  (imax-1,1,1)   
           rho(-1,1,1) = rho(imax-1,1,1) 

           psi(-2,1,1) = psi(imax-2,1,1)
           phi(-2,1,1) = phi(imax-2,1,1) 
           Bx(-2,1,1)  = Bx (imax-2,1,1) 
           By(-2,1,1)  = By (imax-2,1,1) 
           Bz(-2,1,1)  = Bz (imax-2,1,1)
           Ex(-2,1,1)  = Ex (imax-2,1,1)  
           Ey(-2,1,1)  = Ey (imax-2,1,1)  
           Ez(-2,1,1)  = Ez (imax-2,1,1)  
           Vx(-2,1,1)  = Vx (imax-2,1,1)  
           Vy(-2,1,1)  = Vy (imax-2,1,1)  
           Vz(-2,1,1)  = Vz (imax-2,1,1)  
           q(-2,1,1)   = q  (imax-2,1,1)   
           p(-2,1,1)   = p  (imax-2,1,1)   
           rho(-2,1,1) = rho(imax-2,1,1) 

           psi(-3,1,1) = psi(imax-3,1,1)
           phi(-3,1,1) = phi(imax-3,1,1) 
           Bx(-3,1,1)  = Bx (imax-3,1,1) 
           By(-3,1,1)  = By (imax-3,1,1) 
           Bz(-3,1,1)  = Bz (imax-3,1,1)
           Ex(-3,1,1)  = Ex (imax-3,1,1)  
           Ey(-3,1,1)  = Ey (imax-3,1,1)  
           Ez(-3,1,1)  = Ez (imax-3,1,1)  
           Vx(-3,1,1)  = Vx (imax-3,1,1)  
           Vy(-3,1,1)  = Vy (imax-3,1,1)  
           Vz(-3,1,1)  = Vz (imax-3,1,1)  
           q(-3,1,1)   = q  (imax-3,1,1)   
           p(-3,1,1)   = p  (imax-3,1,1)   
           rho(-3,1,1) = rho(imax-3,1,1) 

           psi(-4,1,1) = psi(imax-4,1,1)
           phi(-4,1,1) = phi(imax-4,1,1) 
           Bx(-4,1,1)  = Bx (imax-4,1,1) 
           By(-4,1,1)  = By (imax-4,1,1) 
           Bz(-4,1,1)  = Bz (imax-4,1,1)
           Ex(-4,1,1)  = Ex (imax-4,1,1)  
           Ey(-4,1,1)  = Ey (imax-4,1,1)  
           Ez(-4,1,1)  = Ez (imax-4,1,1)  
           Vx(-4,1,1)  = Vx (imax-4,1,1)  
           Vy(-4,1,1)  = Vy (imax-4,1,1)  
           Vz(-4,1,1)  = Vz (imax-4,1,1)  
           q(-4,1,1)   = q  (imax-4,1,1)   
           p(-4,1,1)   = p  (imax-4,1,1)   
           rho(-4,1,1) = rho(imax-4,1,1)

           psi(-5,1,1) = psi(imax-5,1,1)
           phi(-5,1,1) = phi(imax-5,1,1) 
           Bx(-5,1,1)  = Bx (imax-5,1,1) 
           By(-5,1,1)  = By (imax-5,1,1) 
           Bz(-5,1,1)  = Bz (imax-5,1,1)
           Ex(-5,1,1)  = Ex (imax-5,1,1)  
           Ey(-5,1,1)  = Ey (imax-5,1,1)  
           Ez(-5,1,1)  = Ez (imax-5,1,1)  
           Vx(-5,1,1)  = Vx (imax-5,1,1)  
           Vy(-5,1,1)  = Vy (imax-5,1,1)  
           Vz(-5,1,1)  = Vz (imax-5,1,1)  
           q(-5,1,1)   = q  (imax-5,1,1)   
           p(-5,1,1)   = p  (imax-5,1,1)   
           rho(-5,1,1) = rho(imax-5,1,1)

           psi(-6,1,1) = psi(imax-6,1,1)
           phi(-6,1,1) = phi(imax-6,1,1) 
           Bx(-6,1,1)  = Bx (imax-6,1,1) 
           By(-6,1,1)  = By (imax-6,1,1) 
           Bz(-6,1,1)  = Bz (imax-6,1,1)
           Ex(-6,1,1)  = Ex (imax-6,1,1)  
           Ey(-6,1,1)  = Ey (imax-6,1,1)  
           Ez(-6,1,1)  = Ez (imax-6,1,1)  
           Vx(-6,1,1)  = Vx (imax-6,1,1)  
           Vy(-6,1,1)  = Vy (imax-6,1,1)  
           Vz(-6,1,1)  = Vz (imax-6,1,1)  
           q(-6,1,1)   = q  (imax-6,1,1)   
           p(-6,1,1)   = p  (imax-6,1,1)   
           rho(-6,1,1) = rho(imax-6,1,1) 


           ! Right  Boundary
           ! -----------------
!!$           psi(imax,1,1)  = psi(1,1,1) 
!!$           phi(imax,1,1)  = phi(1,1,1) 
!!$           Bx (imax,1,1)  = Bx (1,1,1)  
!!$           By (imax,1,1)  = By (1,1,1)  
!!$           Bz (imax,1,1)  = Bz (1,1,1)  
!!$           Ex (imax,1,1)  = Ex (1,1,1)  
!!$           Ey (imax,1,1)  = Ey (1,1,1)  
!!$           Ez (imax,1,1)  = Ez (1,1,1)  
!!$           Vx (imax,1,1)  = Vx (1,1,1)  
!!$           Vy (imax,1,1)  = Vy (1,1,1)  
!!$           Vz (imax,1,1)  = Vz (1,1,1)  
!!$           q  (imax,1,1)  = q  (1,1,1)   
!!$           p  (imax,1,1)  = p  (1,1,1)   
!!$           rho(imax,1,1)  = rho(1,1,1) 

           psi(imax+1,1,1)  = psi(1,1,1) 
           phi(imax+1,1,1)  = phi(1,1,1) 
           Bx (imax+1,1,1)  = Bx (1,1,1)  
           By (imax+1,1,1)  = By (1,1,1)  
           Bz (imax+1,1,1)  = Bz (1,1,1)  
           Ex (imax+1,1,1)  = Ex (1,1,1)  
           Ey (imax+1,1,1)  = Ey (1,1,1)  
           Ez (imax+1,1,1)  = Ez (1,1,1)  
           Vx (imax+1,1,1)  = Vx (1,1,1)  
           Vy (imax+1,1,1)  = Vy (1,1,1)  
           Vz (imax+1,1,1)  = Vz (1,1,1)  
           q  (imax+1,1,1)  = q  (1,1,1)   
           p  (imax+1,1,1)  = p  (1,1,1)   
           rho(imax+1,1,1)  = rho(1,1,1) 

           psi(imax+2,1,1)  = psi(2,1,1) 
           phi(imax+2,1,1)  = phi(2,1,1) 
           Bx (imax+2,1,1)  = Bx (2,1,1)  
           By (imax+2,1,1)  = By (2,1,1)  
           Bz (imax+2,1,1)  = Bz (2,1,1)  
           Ex (imax+2,1,1)  = Ex (2,1,1)  
           Ey (imax+2,1,1)  = Ey (2,1,1)  
           Ez (imax+2,1,1)  = Ez (2,1,1)  
           Vx (imax+2,1,1)  = Vx (2,1,1)  
           Vy (imax+2,1,1)  = Vy (2,1,1)  
           Vz (imax+2,1,1)  = Vz (2,1,1)  
           q  (imax+2,1,1)  = q  (2,1,1)   
           p  (imax+2,1,1)  = p  (2,1,1)   
           rho(imax+2,1,1)  = rho(2,1,1) 

           psi(imax+3,1,1)  = psi(3,1,1) 
           phi(imax+3,1,1)  = phi(3,1,1) 
           Bx (imax+3,1,1)  = Bx (3,1,1)  
           By (imax+3,1,1)  = By (3,1,1)  
           Bz (imax+3,1,1)  = Bz (3,1,1)  
           Ex (imax+3,1,1)  = Ex (3,1,1)  
           Ey (imax+3,1,1)  = Ey (3,1,1)  
           Ez (imax+3,1,1)  = Ez (3,1,1)  
           Vx (imax+3,1,1)  = Vx (3,1,1)  
           Vy (imax+3,1,1)  = Vy (3,1,1)  
           Vz (imax+3,1,1)  = Vz (3,1,1)  
           q  (imax+3,1,1)  = q  (3,1,1)   
           p  (imax+3,1,1)  = p  (3,1,1)   
           rho(imax+3,1,1)  = rho(3,1,1) 

           psi(imax+4,1,1)  = psi(4,1,1) 
           phi(imax+4,1,1)  = phi(4,1,1) 
           Bx (imax+4,1,1)  = Bx (4,1,1)  
           By (imax+4,1,1)  = By (4,1,1)  
           Bz (imax+4,1,1)  = Bz (4,1,1)  
           Ex (imax+4,1,1)  = Ex (4,1,1)  
           Ey (imax+4,1,1)  = Ey (4,1,1)  
           Ez (imax+4,1,1)  = Ez (4,1,1)  
           Vx (imax+4,1,1)  = Vx (4,1,1)  
           Vy (imax+4,1,1)  = Vy (4,1,1)  
           Vz (imax+4,1,1)  = Vz (4,1,1)  
           q  (imax+4,1,1)  = q  (4,1,1)   
           p  (imax+4,1,1)  = p  (4,1,1)   
           rho(imax+4,1,1)  = rho(4,1,1)

           psi(imax+5,1,1)  = psi(5,1,1) 
           phi(imax+5,1,1)  = phi(5,1,1) 
           Bx (imax+5,1,1)  = Bx (5,1,1)  
           By (imax+5,1,1)  = By (5,1,1)  
           Bz (imax+5,1,1)  = Bz (5,1,1)  
           Ex (imax+5,1,1)  = Ex (5,1,1)  
           Ey (imax+5,1,1)  = Ey (5,1,1)  
           Ez (imax+5,1,1)  = Ez (5,1,1)  
           Vx (imax+5,1,1)  = Vx (5,1,1)  
           Vy (imax+5,1,1)  = Vy (5,1,1)  
           Vz (imax+5,1,1)  = Vz (5,1,1)  
           q  (imax+5,1,1)  = q  (5,1,1)   
           p  (imax+5,1,1)  = p  (5,1,1)   
           rho(imax+5,1,1)  = rho(5,1,1)

           psi(imax+6,1,1)  = psi(6,1,1) 
           phi(imax+6,1,1)  = phi(6,1,1) 
           Bx (imax+6,1,1)  = Bx (6,1,1)  
           By (imax+6,1,1)  = By (6,1,1)  
           Bz (imax+6,1,1)  = Bz (6,1,1)  
           Ex (imax+6,1,1)  = Ex (6,1,1)  
           Ey (imax+6,1,1)  = Ey (6,1,1)  
           Ez (imax+6,1,1)  = Ez (6,1,1)  
           Vx (imax+6,1,1)  = Vx (6,1,1)  
           Vy (imax+6,1,1)  = Vy (6,1,1)  
           Vz (imax+6,1,1)  = Vz (6,1,1)  
           q  (imax+6,1,1)  = q  (6,1,1)   
           p  (imax+6,1,1)  = p  (6,1,1)   
           rho(imax+6,1,1)  = rho(6,1,1) 

    else if (DIM == 1 .and. BOUND == 3) then

       ! Reconnection boundary Aloy & Mimica

  !-----------------
  ! x Direction
  !-----------------

           posx   = - 0.5d0 * Ly 
           posy   =   0.0d0 

           r      =   sqrt(1 -(posy**2/Lr**2))


           ! Left Boundary
           ! -----------------


           psi(0,1,1) = 0.d0 !psi(1, 1,1)
           phi(0,1,1) = 0.d0 !phi(1, 1,1) 
           Bx (0,1,1) = (B0*posy*delta_rec/Lr**2) * (1.d0+log(cosh(posx/delta_rec))/r) 
           By (0,1,1) = B0 * r * tanh(posx/delta_rec) 
           Bz (0,1,1) = 0.d0 !Bz (1, 1,1)
           Vx (0,1,1) = 0.d0 !Vx (1, 1,1)  
           Vy (0,1,1) = Vy (1, 1,1)  
           Vz (0,1,1) = 0.d0 !Vz (1, 1,1)
           Ex (0,1,1) = Ex (1, 1,1)  
           Ey (0,1,1) = Ey (1, 1,1)  
           Ez (0,1,1) = Ez (1, 1,1)
           q  (0,1,1) = 0.d0   !q  (1, 1,1)   
           rho(0,1,1) = rho_up !rho(1, 1,1)
           p  (0,1,1) = P0 * r**2 / (cosh(posx/delta_rec))**2 + P_up


           psi(-1,1,1) = psi(0,1,1)
           phi(-1,1,1) = phi(0,1,1) 
           Bx (-1,1,1) = Bx (0,1,1) 
           By (-1,1,1) = By (0,1,1) 
           Bz (-1,1,1) = Bz (0,1,1)
           Vx (-1,1,1) = Vx (0,1,1)  
           Vy (-1,1,1) = Vy (0,1,1)  
           Vz (-1,1,1) = Vz (0,1,1)
           Ex (-1,1,1) = Ex (0,1,1) 
           Ey (-1,1,1) = Ey (0,1,1) 
           Ez (-1,1,1) = Ez (0,1,1) 
           q  (-1,1,1) = q  (0,1,1)   
           p  (-1,1,1) = p  (0,1,1) 
           rho(-1,1,1) = rho(0,1,1)

           psi(-2,1,1) = psi(-1,1,1)
           phi(-2,1,1) = phi(-1,1,1) 
           Bx (-2,1,1) = Bx (-1,1,1)
           By (-2,1,1) = By (-1,1,1)
           Bz (-2,1,1) = Bz (-1,1,1)
           Vx (-2,1,1) = Vx (-1,1,1)  
           Vy (-2,1,1) = Vy (-1,1,1)  
           Vz (-2,1,1) = Vz (-1,1,1) 
           Ex (-2,1,1) = Ex (-1,1,1)
           Ey (-2,1,1) = Ey (-1,1,1)
           Ez (-2,1,1) = Ez (-1,1,1)
           q  (-2,1,1) = q  (-1,1,1)   
           p  (-2,1,1) = p  (-1,1,1)
           rho(-2,1,1) = rho(-1,1,1)

           ! Right  Boundary
           ! -----------------

           posx   = - 0.5d0 * Ly + imax * Delx
           posy   =   0.0d0 

           r      =   sqrt(1 -(posy**2/Lr**2))


           psi(imax,1,1) = 0.d0 !psi(imax-1,1,1) 
           phi(imax,1,1) = 0.d0 !phi(imax-1,1,1) 
           Bx (imax,1,1) = (B0*posy*delta_rec/Lr**2) * (1.d0+log(cosh(posx/delta_rec))/r) 
           By (imax,1,1) =  B0 * r * tanh(posx/delta_rec) 
           Bz (imax,1,1) = 0.d0 !Bz (imax-1,1,1)  
           Vx (imax,1,1) = 0.d0 !Vx (imax-1,1,1)  
           Vy (imax,1,1) = Vy (imax-1,1,1)  
           Vz (imax,1,1) = 0.d0 !Vz (imax-1,1,1) 
           Ex (imax,1,1) = Ex (imax-1,1,1)
           Ey (imax,1,1) = Ey (imax-1,1,1)
           Ez (imax,1,1) = Ez (imax-1,1,1)
           q  (imax,1,1) = 0.d0   !q  (imax-1,1,1)   
           rho(imax,1,1) = rho_up !rho(imax-1,1,1)
           p  (imax,1,1) = P0 * r**2 / (cosh(posx/delta_rec))**2 + P_up


           psi(imax+1,1,1) = psi(imax,1,1) 
           phi(imax+1,1,1) = phi(imax,1,1) 
           Bx (imax+1,1,1) = Bx (imax,1,1)
           By (imax+1,1,1) = By (imax,1,1)
           Bz (imax+1,1,1) = Bz (imax,1,1)  
           Vx (imax+1,1,1) = Vx (imax,1,1)  
           Vy (imax+1,1,1) = Vy (imax,1,1)  
           Vz (imax+1,1,1) = Vz (imax,1,1) 
           Ex (imax+1,1,1) = Ex (imax,1,1)
           Ey (imax+1,1,1) = Ey (imax,1,1)
           Ez (imax+1,1,1) = Ez (imax,1,1)
           q  (imax+1,1,1) = q  (imax,1,1)   
           p  (imax+1,1,1) = p  (imax,1,1)
           rho(imax+1,1,1) = rho(imax,1,1)


           psi(imax+2,1,1) = psi(imax+1,1,1) 
           phi(imax+2,1,1) = phi(imax+1,1,1) 
           Bx (imax+2,1,1) = Bx (imax+1,1,1)
           By (imax+2,1,1) = By (imax+1,1,1)
           Bz (imax+2,1,1) = Bz (imax+1,1,1)  
           Vx (imax+2,1,1) = Vx (imax+1,1,1)  
           Vy (imax+2,1,1) = Vy (imax+1,1,1)  
           Vz (imax+2,1,1) = Vz (imax+1,1,1) 
           Ex (imax+2,1,1) = Ex (imax+1,1,1)
           Ey (imax+2,1,1) = Ey (imax+1,1,1)
           Ez (imax+2,1,1) = Ez (imax+1,1,1)
           q  (imax+2,1,1) = q  (imax+1,1,1)   
           p  (imax+2,1,1) = p  (imax+1,1,1)
           rho(imax+2,1,1) = rho(imax+1,1,1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    else if (DIM == 1 .and. BOUND == 4) then


       ! Copy boundary 1D + perfect conductor


  !-----------------
  ! X Direction
  !-----------------

           ! Left Boundary
           ! -----------------

           psi(0,1,1) = psi(1,1,1)
           phi(0,1,1) = phi(1,1,1) 
           Bx(0,1,1)  = 0.d0 !Bx(1,1,1) 
           By(0,1,1)  = By(1,1,1) 
           Bz(0,1,1)  = Bz(1,1,1)
           Ex(0,1,1)  = Ex(1,1,1)  
           Ey(0,1,1)  = Ey(1,1,1)  
           Ez(0,1,1)  = Ez(1,1,1)  
           Vx(0,1,1)  = 0.d0 !Vx(1,1,1)  
           Vy(0,1,1)  = Vy(1,1,1)  
           Vz(0,1,1)  = Vz(1,1,1)  
           q(0,1,1)   = q(1,1,1)   
           p(0,1,1)   = p(1,1,1)   
           rho(0,1,1) = rho(1,1,1) 

           psi(-1,1,1) = psi(0,1,1)
           phi(-1,1,1) = phi(0,1,1) 
           Bx(-1,1,1)  = Bx(0,1,1) 
           By(-1,1,1)  = By(0,1,1) 
           Bz(-1,1,1)  = Bz(0,1,1)
           Ex(-1,1,1)  = Ex(0,1,1)  
           Ey(-1,1,1)  = Ey(0,1,1)  
           Ez(-1,1,1)  = Ez(0,1,1)  
           Vx(-1,1,1)  = Vx(0,1,1)  
           Vy(-1,1,1)  = Vy(0,1,1)  
           Vz(-1,1,1)  = Vz(0,1,1)  
           q(-1,1,1)   = q(0,1,1)   
           p(-1,1,1)   = p(0,1,1)   
           rho(-1,1,1) = rho(0,1,1) 

           psi(-2,1,1) = psi(0,1,1)
           phi(-2,1,1) = phi(0,1,1) 
           Bx(-2,1,1)  = Bx(0,1,1) 
           By(-2,1,1)  = By(0,1,1) 
           Bz(-2,1,1)  = Bz(0,1,1)
           Ex(-2,1,1)  = Ex(0,1,1)  
           Ey(-2,1,1)  = Ey(0,1,1)  
           Ez(-2,1,1)  = Ez(0,1,1)  
           Vx(-2,1,1)  = Vx(0,1,1)  
           Vy(-2,1,1)  = Vy(0,1,1)  
           Vz(-2,1,1)  = Vz(0,1,1)  
           q(-2,1,1)   = q(0,1,1)   
           p(-2,1,1)   = p(0,1,1)   
           rho(-2,1,1) = rho(0,1,1) 

           psi(-3,1,1) = psi(0,1,1)
           phi(-3,1,1) = phi(0,1,1) 
           Bx(-3,1,1)  = Bx(0,1,1) 
           By(-3,1,1)  = By(0,1,1) 
           Bz(-3,1,1)  = Bz(0,1,1)
           Ex(-3,1,1)  = Ex(0,1,1)  
           Ey(-3,1,1)  = Ey(0,1,1)  
           Ez(-3,1,1)  = Ez(0,1,1)  
           Vx(-3,1,1)  = Vx(0,1,1)  
           Vy(-3,1,1)  = Vy(0,1,1)  
           Vz(-3,1,1)  = Vz(0,1,1)  
           q(-3,1,1)   = q(0,1,1)   
           p(-3,1,1)   = p(0,1,1)   
           rho(-3,1,1) = rho(0,1,1) 

           psi(-4,1,1) = psi(0,1,1)
           phi(-4,1,1) = phi(0,1,1) 
           Bx(-4,1,1)  = Bx(0,1,1) 
           By(-4,1,1)  = By(0,1,1) 
           Bz(-4,1,1)  = Bz(0,1,1)
           Ex(-4,1,1)  = Ex(0,1,1)  
           Ey(-4,1,1)  = Ey(0,1,1)  
           Ez(-4,1,1)  = Ez(0,1,1)  
           Vx(-4,1,1)  = Vx(0,1,1)  
           Vy(-4,1,1)  = Vy(0,1,1)  
           Vz(-4,1,1)  = Vz(0,1,1)  
           q(-4,1,1)   = q(0,1,1)   
           p(-4,1,1)   = p(0,1,1)   
           rho(-4,1,1) = rho(0,1,1) 

           ! Right  Boundary
           ! -----------------

           psi(imax,1,1) = psi(imax-1,1,1) 
           phi(imax,1,1) = phi(imax-1,1,1) 
           Bx(imax,1,1)  = 0.d0 !Bx(imax-1,1,1)  
           By(imax,1,1)  = By(imax-1,1,1)  
           Bz(imax,1,1)  = Bz(imax-1,1,1)  
           Ex(imax,1,1)  = Ex(imax-1,1,1)  
           Ey(imax,1,1)  = Ey(imax-1,1,1)  
           Ez(imax,1,1)  = Ez(imax-1,1,1)  
           Vx(imax,1,1)  = 0.d0 !Vx(imax-1,1,1)  
           Vy(imax,1,1)  = Vy(imax-1,1,1)  
           Vz(imax,1,1)  = Vz(imax-1,1,1)  
           q(imax,1,1)   = q(imax-1,1,1)   
           p(imax,1,1)   = p(imax-1,1,1)   
           rho(imax,1,1) = rho(imax-1,1,1) 

           psi(imax+1,1,1) = psi(imax,1,1) 
           phi(imax+1,1,1) = phi(imax,1,1) 
           Bx(imax+1,1,1)  = Bx(imax,1,1)  
           By(imax+1,1,1)  = By(imax,1,1)  
           Bz(imax+1,1,1)  = Bz(imax,1,1)  
           Ex(imax+1,1,1)  = Ex(imax,1,1)  
           Ey(imax+1,1,1)  = Ey(imax,1,1)  
           Ez(imax+1,1,1)  = Ez(imax,1,1)  
           Vx(imax+1,1,1)  = Vx(imax,1,1)  
           Vy(imax+1,1,1)  = Vy(imax,1,1)  
           Vz(imax+1,1,1)  = Vz(imax,1,1)  
           q(imax+1,1,1)   = q(imax,1,1)   
           p(imax+1,1,1)   = p(imax,1,1)   
           rho(imax+1,1,1) = rho(imax,1,1) 

           psi(imax+2,1,1) = psi(imax,1,1) 
           phi(imax+2,1,1) = phi(imax,1,1) 
           Bx(imax+2,1,1)  = Bx(imax,1,1)  
           By(imax+2,1,1)  = By(imax,1,1)  
           Bz(imax+2,1,1)  = Bz(imax,1,1)  
           Ex(imax+2,1,1)  = Ex(imax,1,1)  
           Ey(imax+2,1,1)  = Ey(imax,1,1)  
           Ez(imax+2,1,1)  = Ez(imax,1,1)  
           Vx(imax+2,1,1)  = Vx(imax,1,1)  
           Vy(imax+2,1,1)  = Vy(imax,1,1)  
           Vz(imax+2,1,1)  = Vz(imax,1,1)  
           q(imax+2,1,1)   = q(imax,1,1)   
           p(imax+2,1,1)   = p(imax,1,1)   
           rho(imax+2,1,1) = rho(imax,1,1) 

           psi(imax+3,1,1) = psi(imax,1,1) 
           phi(imax+3,1,1) = phi(imax,1,1) 
           Bx(imax+3,1,1)  = Bx(imax,1,1)  
           By(imax+3,1,1)  = By(imax,1,1)  
           Bz(imax+3,1,1)  = Bz(imax,1,1)  
           Ex(imax+3,1,1)  = Ex(imax,1,1)  
           Ey(imax+3,1,1)  = Ey(imax,1,1)  
           Ez(imax+3,1,1)  = Ez(imax,1,1)  
           Vx(imax+3,1,1)  = Vx(imax,1,1)  
           Vy(imax+3,1,1)  = Vy(imax,1,1)  
           Vz(imax+3,1,1)  = Vz(imax,1,1)  
           q(imax+3,1,1)   = q(imax,1,1)   
           p(imax+3,1,1)   = p(imax,1,1)   
           rho(imax+3,1,1) = rho(imax,1,1) 

           psi(imax+4,1,1) = psi(imax,1,1) 
           phi(imax+4,1,1) = phi(imax,1,1) 
           Bx(imax+4,1,1)  = Bx(imax,1,1)  
           By(imax+4,1,1)  = By(imax,1,1)  
           Bz(imax+4,1,1)  = Bz(imax,1,1)  
           Ex(imax+4,1,1)  = Ex(imax,1,1)  
           Ey(imax+4,1,1)  = Ey(imax,1,1)  
           Ez(imax+4,1,1)  = Ez(imax,1,1)  
           Vx(imax+4,1,1)  = Vx(imax,1,1)  
           Vy(imax+4,1,1)  = Vy(imax,1,1)  
           Vz(imax+4,1,1)  = Vz(imax,1,1)  
           q(imax+4,1,1)   = q(imax,1,1)   
           p(imax+4,1,1)   = p(imax,1,1)   
           rho(imax+4,1,1) = rho(imax,1,1) 





       ! Perfect Conductor and anti-simetric reflection TM_III

  !-----------------
  ! X Direction
  !-----------------

           ! Left Boundary
           ! -----------------
!!$           psi(0,1,1) = 2.d0 * psi(1,1,1) - psi(2,1,1)
!!$           phi(0,1,1) = 2.d0 * phi(1,1,1) - phi(2,1,1)
!!$           Bx(0,1,1)  =-Bx(2,1,1)
!!$           By(0,1,1)  = 2.d0 * By (1,1,1) - By (2,1,1) 
!!$           Bz(0,1,1)  = 2.d0 * Bz (1,1,1) - Bz (2,1,1) 
!!$           Ex(0,1,1)  = 2.d0 * Ex (1,1,1) - Ex (2,1,1)  
!!$           Ey(0,1,1)  = 2.d0 * Ey (1,1,1) - Ey (2,1,1)  
!!$           Ez(0,1,1)  = 2.d0 * Ez (1,1,1) - E2 (2,1,1)  
!!$           Vx(0,1,1)  = 0.d0
!!$           Vy(0,1,1)  = 2.d0 * Vy (1,1,1) - Vy (2,1,1)  
!!$           Vz(0,1,1)  = 2.d0 * Vz (1,1,1) - Vz (2,1,1)  
!!$           q(0,1,1)   = 2.d0 * q  (1,1,1) - q  (2,1,1)      
!!$           p(0,1,1)   = 2.d0 * p  (1,1,1) - p  (2,1,1)   
!!$           rho(0,1,1) = 2.d0 * rho(1,1,1) - rho(2,1,1) 

!!$           psi(-1,1,1) = 2.d0 * psi(0,1,1) - psi(1,1,1)
!!$           phi(-1,1,1) = 2.d0 * phi(0,1,1) - phi(1,1,1)
!!$           Bx(-1,1,1)  = - Bx(1,1,1)
!!$           By(-1,1,1)  = 2.d0 * By (0,1,1) - By (1,1,1) 
!!$           Bz(-1,1,1)  = 2.d0 * Bz (0,1,1) - Bz (1,1,1) 
!!$           Ex(-1,1,1)  = 2.d0 * Ex (0,1,1) - Ex (1,1,1)  
!!$           Ey(-1,1,1)  = 2.d0 * Ey (0,1,1) - Ey (1,1,1) 
!!$           Ez(-1,1,1)  = 2.d0 * Ez (0,1,1) - E2 (1,1,1)  
!!$           Vx(-1,1,1)  = - vx(1,1,1)
!!$           Vy(-1,1,1)  = 2.d0 * Vy (0,1,1) - Vy (1,1,1)
!!$           Vz(-1,1,1)  = 2.d0 * Vz (0,1,1) - Vz (1,1,1)
!!$           q(-1,1,1)   = 2.d0 * q  (0,1,1) - q  (1,1,1)
!!$           p(-1,1,1)   = 2.d0 * p  (0,1,1) - p  (1,1,1)
!!$           rho(-1,1,1) = 2.d0 * rho(0,1,1) - rho(1,1,1)
!!$
!!$           psi(-2,1,1) = 2.d0 * psi(0,1,1) - psi(2,1,1)
!!$           phi(-2,1,1) = 2.d0 * phi(0,1,1) - phi(2,1,1)
!!$           Bx(-2,1,1)  = - Bx(2,1,1)
!!$           By(-2,1,1)  = 2.d0 * By (0,1,1) - By (2,1,1)
!!$           Bz(-2,1,1)  = 2.d0 * Bz (0,1,1) - Bz (2,1,1)
!!$           Ex(-2,1,1)  = 2.d0 * Ex (0,1,1) - Ex (2,1,1)
!!$           Ey(-2,1,1)  = 2.d0 * Ey (0,1,1) - Ey (2,1,1)
!!$           Ez(-2,1,1)  = 2.d0 * Ez (0,1,1) - E2 (2,1,1)
!!$           Vx(-2,1,1)  = - Vx(2,1,1)
!!$           Vy(-2,1,1)  = 2.d0 * Vy (0,1,1) - Vy (2,1,1)
!!$           Vz(-2,1,1)  = 2.d0 * Vz (0,1,1) - Vz (2,1,1)
!!$           q(-2,1,1)   = 2.d0 * q  (0,1,1) - q  (2,1,1)
!!$           p(-2,1,1)   = 2.d0 * p  (0,1,1) - p  (2,1,1)
!!$           rho(-2,1,1) = 2.d0 * rho(0,1,1) - rho(2,1,1)
!!$
!!$           psi(-3,1,1) = 2.d0 * psi(0,1,1) - psi(3,1,1)
!!$           phi(-3,1,1) = 2.d0 * phi(0,1,1) - phi(3,1,1)
!!$           Bx(-3,1,1)  = - Bx(3,1,1)
!!$           By(-3,1,1)  = 2.d0 * By (0,1,1) - By (3,1,1)
!!$           Bz(-3,1,1)  = 2.d0 * Bz (0,1,1) - Bz (3,1,1)
!!$           Ex(-3,1,1)  = 2.d0 * Ex (0,1,1) - Ex (3,1,1)
!!$           Ey(-3,1,1)  = 2.d0 * Ey (0,1,1) - Ey (3,1,1)
!!$           Ez(-3,1,1)  = 2.d0 * Ez (0,1,1) - E2 (3,1,1)
!!$           Vx(-3,1,1)  = - Vx(3,1,1)
!!$           Vy(-3,1,1)  = 2.d0 * Vy (0,1,1) - Vy (3,1,1)
!!$           Vz(-3,1,1)  = 2.d0 * Vz (0,1,1) - Vz (3,1,1)
!!$           q(-3,1,1)   = 2.d0 * q  (0,1,1) - q  (3,1,1)
!!$           p(-3,1,1)   = 2.d0 * p  (0,1,1) - p  (3,1,1)
!!$           rho(-3,1,1) = 2.d0 * rho(0,1,1) - rho(3,1,1)
!!$
!!$           psi(-4,1,1) = 2.d0 * psi(0,1,1) - psi(4,1,1)
!!$           phi(-4,1,1) = 2.d0 * phi(0,1,1) - phi(4,1,1)
!!$           Bx(-4,1,1)  = - Bx(4,1,1)
!!$           By(-4,1,1)  = 2.d0 * By (0,1,1) - By (4,1,1)
!!$           Bz(-4,1,1)  = 2.d0 * Bz (0,1,1) - Bz (4,1,1)
!!$           Ex(-4,1,1)  = 2.d0 * Ex (0,1,1) - Ex (4,1,1)
!!$           Ey(-4,1,1)  = 2.d0 * Ey (0,1,1) - Ey (4,1,1)
!!$           Ez(-4,1,1)  = 2.d0 * Ez (0,1,1) - E2 (4,1,1)
!!$           Vx(-4,1,1)  = - Vx(4,1,1)
!!$           Vy(-4,1,1)  = 2.d0 * Vy (0,1,1) - Vy (4,1,1)
!!$           Vz(-4,1,1)  = 2.d0 * Vz (0,1,1) - Vz (4,1,1)
!!$           q(-4,1,1)   = 2.d0 * q  (0,1,1) - q  (4,1,1)
!!$           p(-4,1,1)   = 2.d0 * p  (0,1,1) - p  (4,1,1)
!!$           rho(-4,1,1) = 2.d0 * rho(0,1,1) - rho(4,1,1)


           ! Right  Boundary
           ! -----------------
!!$           psi(imax,1,1)  = 2.d0 * psi(imax-1,1,1) - psi(imax-2,1,1)
!!$           phi(imax,1,1)  = 2.d0 * phi(imax-1,1,1) - phi(imax-2,1,1)
!!$           Bx (imax,1,1)  = 0.d0
!!$           By (imax,1,1)  = 2.d0 * By (imax-1,1,1) - By (imax-2,1,1)  
!!$           Bz (imax,1,1)  = 2.d0 * Bz (imax-1,1,1) - Bz (imax-2,1,1)
!!$           Ex (imax,1,1)  = 2.d0 * Ex (imax-1,1,1) - Ex (imax-2,1,1)
!!$           Ey (imax,1,1)  = 2.d0 * Ey (imax-1,1,1) - Ey (imax-2,1,1)
!!$           Ez (imax,1,1)  = 2.d0 * Ez (imax-1,1,1) - Ez (imax-2,1,1)
!!$           Vx (imax,1,1)  = 0.d0
!!$           Vy (imax,1,1)  = 2.d0 * Vy (imax-1,1,1) - Vy (imax-2,1,1)
!!$           Vz (imax,1,1)  = 2.d0 * Vz (imax-1,1,1) - Vz (imax-2,1,1)
!!$           q  (imax,1,1)  = 2.d0 * q  (imax-1,1,1) - q  (imax-2,1,1)
!!$           p  (imax,1,1)  = 2.d0 * p  (imax-1,1,1) - p  (imax-2,1,1)
!!$           rho(imax,1,1)  = 2.d0 * rho(imax-1,1,1) - rho(imax-2,1,1)

!!$           psi(imax+1,1,1)  = 2.d0 * psi(imax,1,1) - psi(imax-1,1,1)
!!$           phi(imax+1,1,1)  = 2.d0 * phi(imax,1,1) - phi(imax-1,1,1)
!!$           Bx (imax+1,1,1)  = - Bx (imax-1,1,1)
!!$           By (imax+1,1,1)  = 2.d0 * By (imax,1,1) - By (imax-1,1,1)
!!$           Bz (imax+1,1,1)  = 2.d0 * Bz (imax,1,1) - Bz (imax-1,1,1)
!!$           Ex (imax+1,1,1)  = 2.d0 * Ex (imax,1,1) - Ex (imax-1,1,1)
!!$           Ey (imax+1,1,1)  = 2.d0 * Ey (imax,1,1) - Ey (imax-1,1,1)
!!$           Ez (imax+1,1,1)  = 2.d0 * Ez (imax,1,1) - Ez (imax-1,1,1)
!!$           Vx (imax+1,1,1)  = - Vx (imax-1,1,1)
!!$           Vy (imax+1,1,1)  = 2.d0 * Vy (imax,1,1) - Vy (imax-1,1,1)
!!$           Vz (imax+1,1,1)  = 2.d0 * Vz (imax,1,1) - Vz (imax-1,1,1)
!!$           q  (imax+1,1,1)  = 2.d0 * q  (imax,1,1) - q  (imax-1,1,1)
!!$           p  (imax+1,1,1)  = 2.d0 * p  (imax,1,1) - p  (imax-1,1,1)
!!$           rho(imax+1,1,1)  = 2.d0 * rho(imax,1,1) - rho(imax-1,1,1)
!!$
!!$           psi(imax+2,1,1)  = 2.d0 * psi(imax,1,1) - psi(imax-2,1,1)
!!$           phi(imax+2,1,1)  = 2.d0 * phi(imax,1,1) - phi(imax-2,1,1)
!!$           Bx (imax+2,1,1)  = - Bx (imax-2,1,1)
!!$           By (imax+2,1,1)  = 2.d0 * By (imax,1,1) - By (imax-2,1,1)
!!$           Bz (imax+2,1,1)  = 2.d0 * Bz (imax,1,1) - Bz (imax-2,1,1)
!!$           Ex (imax+2,1,1)  = 2.d0 * Ex (imax,1,1) - Ex (imax-2,1,1)
!!$           Ey (imax+2,1,1)  = 2.d0 * Ey (imax,1,1) - Ey (imax-2,1,1)
!!$           Ez (imax+2,1,1)  = 2.d0 * Ez (imax,1,1) - Ez (imax-2,1,1)
!!$           Vx (imax+2,1,1)  = - Vx (imax-2,1,1)
!!$           Vy (imax+2,1,1)  = 2.d0 * Vy (imax,1,1) - Vy (imax-2,1,1)
!!$           Vz (imax+2,1,1)  = 2.d0 * Vz (imax,1,1) - Vz (imax-2,1,1)
!!$           q  (imax+2,1,1)  = 2.d0 * q  (imax,1,1) - q  (imax-2,1,1)
!!$           p  (imax+2,1,1)  = 2.d0 * p  (imax,1,1) - p  (imax-2,1,1)
!!$           rho(imax+2,1,1)  = 2.d0 * rho(imax,1,1) - rho(imax-2,1,1)
!!$
!!$           psi(imax+3,1,1)  = 2.d0 * psi(imax,1,1) - psi(imax-3,1,1)
!!$           phi(imax+3,1,1)  = 2.d0 * phi(imax,1,1) - phi(imax-3,1,1)
!!$           Bx (imax+3,1,1)  = - Bx (imax-3,1,1)
!!$           By (imax+3,1,1)  = 2.d0 * By (imax,1,1) - By (imax-3,1,1)
!!$           Bz (imax+3,1,1)  = 2.d0 * Bz (imax,1,1) - Bz (imax-3,1,1)
!!$           Ex (imax+3,1,1)  = 2.d0 * Ex (imax,1,1) - Ex (imax-3,1,1)
!!$           Ey (imax+3,1,1)  = 2.d0 * Ey (imax,1,1) - Ey (imax-3,1,1)
!!$           Ez (imax+3,1,1)  = 2.d0 * Ez (imax,1,1) - Ez (imax-3,1,1)
!!$           Vx (imax+3,1,1)  = - Vx (imax-3,1,1)
!!$           Vy (imax+3,1,1)  = 2.d0 * Vy (imax,1,1) - Vy (imax-3,1,1)
!!$           Vz (imax+3,1,1)  = 2.d0 * Vz (imax,1,1) - Vz (imax-3,1,1)
!!$           q  (imax+3,1,1)  = 2.d0 * q  (imax,1,1) - q  (imax-3,1,1)
!!$           p  (imax+3,1,1)  = 2.d0 * p  (imax,1,1) - p  (imax-3,1,1)
!!$           rho(imax+3,1,1)  = 2.d0 * rho(imax,1,1) - rho(imax-3,1,1)
!!$
!!$           psi(imax+4,1,1)  = 2.d0 * psi(imax,1,1) - psi(imax-4,1,1)
!!$           phi(imax+4,1,1)  = 2.d0 * phi(imax,1,1) - phi(imax-4,1,1)
!!$           Bx (imax+4,1,1)  = - Bx (imax-4,1,1)
!!$           By (imax+4,1,1)  = 2.d0 * By (imax,1,1) - By (imax-4,1,1)
!!$           Bz (imax+4,1,1)  = 2.d0 * Bz (imax,1,1) - Bz (imax-4,1,1)
!!$           Ex (imax+4,1,1)  = 2.d0 * Ex (imax,1,1) - Ex (imax-4,1,1)
!!$           Ey (imax+4,1,1)  = 2.d0 * Ey (imax,1,1) - Ey (imax-4,1,1)
!!$           Ez (imax+4,1,1)  = 2.d0 * Ez (imax,1,1) - Ez (imax-4,1,1)
!!$           Vx (imax+4,1,1)  = - Vx (imax-4,1,1)
!!$           Vy (imax+4,1,1)  = 2.d0 * Vy (imax,1,1) - Vy (imax-4,1,1)
!!$           Vz (imax+4,1,1)  = 2.d0 * Vz (imax,1,1) - Vz (imax-4,1,1)
!!$           q  (imax+4,1,1)  = 2.d0 * q  (imax,1,1) - q  (imax-4,1,1)
!!$           p  (imax+4,1,1)  = 2.d0 * p  (imax,1,1) - p  (imax-4,1,1)
!!$           rho(imax+4,1,1)  = 2.d0 * rho(imax,1,1) - rho(imax-4,1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    else if (DIM == 2 .and. BOUND == 1) then

       ! Copy boundary 2D

  !-----------------
  ! X Direction
  !-----------------

 !$OMP DO 

    do j=-6,jmax+6
!       do k=1,kmax

           ! Left Boundary
           ! -----------------

!!$           psi( 0,j,1) = psi(1,j,1) 
!!$           phi( 0,j,1) = phi(1,j,1) 
!!$           Bx ( 0,j,1) = Bx (1,j,1)  
!!$           By ( 0,j,1) = By (1,j,1)  
!!$           Bz ( 0,j,1) = Bz (1,j,1)  
!!$           Ex ( 0,j,1) = Ex (1,j,1)  
!!$           Ey ( 0,j,1) = Ey (1,j,1)  
!!$           Ez ( 0,j,1) = Ez (1,j,1)  
!!$           Vx ( 0,j,1) = Vx (1,j,1)  
!!$           Vy ( 0,j,1) = Vy (1,j,1)  
!!$           Vz ( 0,j,1) = Vz (1,j,1)  
!!$           q  ( 0,j,1) = q  (1,j,1)   
!!$           p  ( 0,j,1) = p  (1,j,1)   
!!$           rho( 0,j,1) = rho(1,j,1)  

           psi(-1,j,1) = psi(0,j,1) 
           phi(-1,j,1) = phi(0,j,1) 
           Bx (-1,j,1) = Bx (0,j,1)  
           By (-1,j,1) = By (0,j,1)  
           Bz (-1,j,1) = Bz (0,j,1)  
           Ex (-1,j,1) = Ex (0,j,1)  
           Ey (-1,j,1) = Ey (0,j,1)  
           Ez (-1,j,1) = Ez (0,j,1)  
           Vx (-1,j,1) = Vx (0,j,1)  
           Vy (-1,j,1) = Vy (0,j,1)  
           Vz (-1,j,1) = Vz (0,j,1)  
           q  (-1,j,1) = q  (0,j,1)   
           p  (-1,j,1) = p  (0,j,1)   
           rho(-1,j,1) = rho(0,j,1)  


           psi(-2,j,1) = psi(-1,j,1) 
           phi(-2,j,1) = phi(-1,j,1) 
           Bx (-2,j,1) = Bx (-1,j,1)  
           By (-2,j,1) = By (-1,j,1)  
           Bz (-2,j,1) = Bz (-1,j,1)  
           Ex (-2,j,1) = Ex (-1,j,1)  
           Ey (-2,j,1) = Ey (-1,j,1)  
           Ez (-2,j,1) = Ez (-1,j,1)  
           Vx (-2,j,1) = Vx (-1,j,1)  
           Vy (-2,j,1) = Vy (-1,j,1)  
           Vz (-2,j,1) = Vz (-1,j,1)  
           q  (-2,j,1) = q  (-1,j,1)   
           p  (-2,j,1) = p  (-1,j,1)
           rho(-2,j,1) = rho(-1,j,1)

           psi(-3,j,1) = psi(-2,j,1) 
           phi(-3,j,1) = phi(-2,j,1) 
           Bx (-3,j,1) = Bx (-2,j,1)  
           By (-3,j,1) = By (-2,j,1)  
           Bz (-3,j,1) = Bz (-2,j,1)  
           Ex (-3,j,1) = Ex (-2,j,1)  
           Ey (-3,j,1) = Ey (-2,j,1)  
           Ez (-3,j,1) = Ez (-2,j,1)  
           Vx (-3,j,1) = Vx (-2,j,1)  
           Vy (-3,j,1) = Vy (-2,j,1)  
           Vz (-3,j,1) = Vz (-2,j,1)  
           q  (-3,j,1) = q  (-2,j,1)   
           p  (-3,j,1) = p  (-2,j,1)
           rho(-3,j,1) = rho(-2,j,1)

           psi(-4,j,1) = psi(-3,j,1) 
           phi(-4,j,1) = phi(-3,j,1) 
           Bx (-4,j,1) = Bx (-3,j,1)  
           By (-4,j,1) = By (-3,j,1)  
           Bz (-4,j,1) = Bz (-3,j,1)  
           Ex (-4,j,1) = Ex (-3,j,1)  
           Ey (-4,j,1) = Ey (-3,j,1)  
           Ez (-4,j,1) = Ez (-3,j,1)  
           Vx (-4,j,1) = Vx (-3,j,1)  
           Vy (-4,j,1) = Vy (-3,j,1)  
           Vz (-4,j,1) = Vz (-3,j,1)  
           q  (-4,j,1) = q  (-3,j,1)   
           p  (-4,j,1) = p  (-3,j,1)
           rho(-4,j,1) = rho(-3,j,1)

           psi(-5,j,1) = psi(-4,j,1) 
           phi(-5,j,1) = phi(-4,j,1) 
           Bx (-5,j,1) = Bx (-4,j,1)  
           By (-5,j,1) = By (-4,j,1)  
           Bz (-5,j,1) = Bz (-4,j,1)  
           Ex (-5,j,1) = Ex (-4,j,1)  
           Ey (-5,j,1) = Ey (-4,j,1)  
           Ez (-5,j,1) = Ez (-4,j,1)  
           Vx (-5,j,1) = Vx (-4,j,1)  
           Vy (-5,j,1) = Vy (-4,j,1)  
           Vz (-5,j,1) = Vz (-4,j,1)  
           q  (-5,j,1) = q  (-4,j,1)   
           p  (-5,j,1) = p  (-4,j,1)
           rho(-5,j,1) = rho(-4,j,1)

           psi(-6,j,1) = psi(-5,j,1) 
           phi(-6,j,1) = phi(-5,j,1) 
           Bx (-6,j,1) = Bx (-5,j,1)  
           By (-6,j,1) = By (-5,j,1)  
           Bz (-6,j,1) = Bz (-5,j,1)  
           Ex (-6,j,1) = Ex (-5,j,1)  
           Ey (-6,j,1) = Ey (-5,j,1)  
           Ez (-6,j,1) = Ez (-5,j,1)  
           Vx (-6,j,1) = Vx (-5,j,1)  
           Vy (-6,j,1) = Vy (-5,j,1)  
           Vz (-6,j,1) = Vz (-5,j,1)  
           q  (-6,j,1) = q  (-5,j,1)   
           p  (-6,j,1) = p  (-5,j,1)
           rho(-6,j,1) = rho(-5,j,1)



           ! Right  Boundary
           ! -----------------
!!$           psi(imax  ,j,1) = psi(imax-1,j,1) 
!!$           phi(imax  ,j,1) = phi(imax-1,j,1) 
!!$           Bx (imax  ,j,1) = Bx (imax-1,j,1) 
!!$           By (imax  ,j,1) = By (imax-1,j,1) 
!!$           Bz (imax  ,j,1) = Bz (imax-1,j,1) 
!!$           Ex (imax  ,j,1) = Ex (imax-1,j,1) 
!!$           Ey (imax  ,j,1) = Ey (imax-1,j,1) 
!!$           Ez (imax  ,j,1) = Ez (imax-1,j,1) 
!!$           Vx (imax  ,j,1) = Vx (imax-1,j,1) 
!!$           Vy (imax  ,j,1) = Vy (imax-1,j,1) 
!!$           Vz (imax  ,j,1) = Vz (imax-1,j,1) 
!!$           q  (imax  ,j,1) = q  (imax-1,j,1) 
!!$           p  (imax  ,j,1) = p  (imax-1,j,1) 
!!$           rho(imax  ,j,1) = rho(imax-1,j,1) 

           psi(imax+1,j,1) = psi(imax  ,j,1) 
           phi(imax+1,j,1) = phi(imax  ,j,1) 
           Bx (imax+1,j,1) = Bx (imax  ,j,1) 
           By (imax+1,j,1) = By (imax  ,j,1) 
           Bz (imax+1,j,1) = Bz (imax  ,j,1) 
           Ex (imax+1,j,1) = Ex (imax  ,j,1) 
           Ey (imax+1,j,1) = Ey (imax  ,j,1) 
           Ez (imax+1,j,1) = Ez (imax  ,j,1) 
           Vx (imax+1,j,1) = Vx (imax  ,j,1) 
           Vy (imax+1,j,1) = Vy (imax  ,j,1) 
           Vz (imax+1,j,1) = Vz (imax  ,j,1) 
           q  (imax+1,j,1) = q  (imax  ,j,1) 
           p  (imax+1,j,1) = p  (imax  ,j,1) 
           rho(imax+1,j,1) = rho(imax  ,j,1) 

           psi(imax+2,j,1) = psi(imax+1,j,1)
           phi(imax+2,j,1) = phi(imax+1,j,1)
           Bx (imax+2,j,1) = Bx (imax+1,j,1) 
           By (imax+2,j,1) = By (imax+1,j,1) 
           Bz (imax+2,j,1) = Bz (imax+1,j,1) 
           Ex (imax+2,j,1) = Ex (imax+1,j,1) 
           Ey (imax+2,j,1) = Ey (imax+1,j,1) 
           Ez (imax+2,j,1) = Ez (imax+1,j,1) 
           Vx (imax+2,j,1) = Vx (imax+1,j,1) 
           Vy (imax+2,j,1) = Vy (imax+1,j,1) 
           Vz (imax+2,j,1) = Vz (imax+1,j,1) 
           q  (imax+2,j,1) = q  (imax+1,j,1) 
           p  (imax+2,j,1) = p  (imax+1,j,1) 
           rho(imax+2,j,1) = rho(imax+1,j,1) 

           psi(imax+3,j,1) = psi(imax+2,j,1)
           phi(imax+3,j,1) = phi(imax+2,j,1)
           Bx (imax+3,j,1) = Bx (imax+2,j,1) 
           By (imax+3,j,1) = By (imax+2,j,1) 
           Bz (imax+3,j,1) = Bz (imax+2,j,1) 
           Ex (imax+3,j,1) = Ex (imax+2,j,1) 
           Ey (imax+3,j,1) = Ey (imax+2,j,1) 
           Ez (imax+3,j,1) = Ez (imax+2,j,1) 
           Vx (imax+3,j,1) = Vx (imax+2,j,1) 
           Vy (imax+3,j,1) = Vy (imax+2,j,1) 
           Vz (imax+3,j,1) = Vz (imax+2,j,1) 
           q  (imax+3,j,1) = q  (imax+2,j,1) 
           p  (imax+3,j,1) = p  (imax+2,j,1) 
           rho(imax+3,j,1) = rho(imax+2,j,1) 

           psi(imax+4,j,1) = psi(imax+3,j,1)
           phi(imax+4,j,1) = phi(imax+3,j,1)
           Bx (imax+4,j,1) = Bx (imax+3,j,1) 
           By (imax+4,j,1) = By (imax+3,j,1) 
           Bz (imax+4,j,1) = Bz (imax+3,j,1) 
           Ex (imax+4,j,1) = Ex (imax+3,j,1) 
           Ey (imax+4,j,1) = Ey (imax+3,j,1) 
           Ez (imax+4,j,1) = Ez (imax+3,j,1) 
           Vx (imax+4,j,1) = Vx (imax+3,j,1) 
           Vy (imax+4,j,1) = Vy (imax+3,j,1) 
           Vz (imax+4,j,1) = Vz (imax+3,j,1) 
           q  (imax+4,j,1) = q  (imax+3,j,1) 
           p  (imax+4,j,1) = p  (imax+3,j,1) 
           rho(imax+4,j,1) = rho(imax+3,j,1)

           psi(imax+5,j,1) = psi(imax+4,j,1)
           phi(imax+5,j,1) = phi(imax+4,j,1)
           Bx (imax+5,j,1) = Bx (imax+4,j,1) 
           By (imax+5,j,1) = By (imax+4,j,1) 
           Bz (imax+5,j,1) = Bz (imax+4,j,1) 
           Ex (imax+5,j,1) = Ex (imax+4,j,1) 
           Ey (imax+5,j,1) = Ey (imax+4,j,1) 
           Ez (imax+5,j,1) = Ez (imax+4,j,1) 
           Vx (imax+5,j,1) = Vx (imax+4,j,1) 
           Vy (imax+5,j,1) = Vy (imax+4,j,1) 
           Vz (imax+5,j,1) = Vz (imax+4,j,1) 
           q  (imax+5,j,1) = q  (imax+4,j,1) 
           p  (imax+5,j,1) = p  (imax+4,j,1) 
           rho(imax+5,j,1) = rho(imax+4,j,1)

           psi(imax+6,j,1) = psi(imax+5,j,1)
           phi(imax+6,j,1) = phi(imax+5,j,1)
           Bx (imax+6,j,1) = Bx (imax+5,j,1) 
           By (imax+6,j,1) = By (imax+5,j,1) 
           Bz (imax+6,j,1) = Bz (imax+5,j,1) 
           Ex (imax+6,j,1) = Ex (imax+5,j,1) 
           Ey (imax+6,j,1) = Ey (imax+5,j,1) 
           Ez (imax+6,j,1) = Ez (imax+5,j,1) 
           Vx (imax+6,j,1) = Vx (imax+5,j,1) 
           Vy (imax+6,j,1) = Vy (imax+5,j,1) 
           Vz (imax+6,j,1) = Vz (imax+5,j,1) 
           q  (imax+6,j,1) = q  (imax+5,j,1) 
           p  (imax+6,j,1) = p  (imax+5,j,1) 
           rho(imax+6,j,1) = rho(imax+5,j,1)


!        end do
     end do
 
!$OMP END DO

  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 

    do i=-6,imax+6
!       do k=1,kmax

           ! Left Boundary
           ! -----------------

!!$           psi(i, 0,1) = psi(i, 1,1)
!!$           phi(i, 0,1) = phi(i, 1,1) 
!!$           Bx (i, 0,1) = Bx (i, 1,1)
!!$           By (i, 0,1) = By (i, 1,1)
!!$           Bz (i, 0,1) = Bz (i, 1,1)
!!$           Ex (i, 0,1) = Ex (i, 1,1)  
!!$           Ey (i, 0,1) = Ey (i, 1,1)  
!!$           Ez (i, 0,1) = Ez (i, 1,1)  
!!$           Vx (i, 0,1) = Vx (i, 1,1)  
!!$           Vy (i, 0,1) = Vy (i, 1,1)  
!!$           Vz (i, 0,1) = Vz (i, 1,1)  
!!$           q  (i, 0,1) = q  (i, 1,1)   
!!$           p  (i, 0,1) = p  (i, 1,1)
!!$           rho(i, 0,1) = rho(i, 1,1)

           psi(i,-1,1) = psi(i, 0,1)
           phi(i,-1,1) = phi(i, 0,1) 
           Bx (i,-1,1) = Bx (i, 0,1)
           By (i,-1,1) = By (i, 0,1)
           Bz (i,-1,1) = Bz (i, 0,1)
           Ex (i,-1,1) = Ex (i, 0,1)  
           Ey (i,-1,1) = Ey (i, 0,1)  
           Ez (i,-1,1) = Ez (i, 0,1)  
           Vx (i,-1,1) = Vx (i, 0,1)  
           Vy (i,-1,1) = Vy (i, 0,1)  
           Vz (i,-1,1) = Vz (i, 0,1)  
           q  (i,-1,1) = q  (i, 0,1)   
           p  (i,-1,1) = p  (i, 0,1)
           rho(i,-1,1) = rho(i, 0,1)

           psi(i,-2,1) = psi(i,-1,1)
           phi(i,-2,1) = phi(i,-1,1) 
           Bx (i,-2,1) = Bx (i,-1,1)
           By (i,-2,1) = By (i,-1,1)
           Bz (i,-2,1) = Bz (i,-1,1)
           Ex (i,-2,1) = Ex (i,-1,1)  
           Ey (i,-2,1) = Ey (i,-1,1)  
           Ez (i,-2,1) = Ez (i,-1,1)  
           Vx (i,-2,1) = Vx (i,-1,1)  
           Vy (i,-2,1) = Vy (i,-1,1)  
           Vz (i,-2,1) = Vz (i,-1,1)  
           q  (i,-2,1) = q  (i,-1,1)   
           p  (i,-2,1) = p  (i,-1,1)
           rho(i,-2,1) = rho(i,-1,1)

           psi(i,-3,1) = psi(i,-2,1)
           phi(i,-3,1) = phi(i,-2,1) 
           Bx (i,-3,1) = Bx (i,-2,1)
           By (i,-3,1) = By (i,-2,1)
           Bz (i,-3,1) = Bz (i,-2,1)
           Ex (i,-3,1) = Ex (i,-2,1)  
           Ey (i,-3,1) = Ey (i,-2,1)  
           Ez (i,-3,1) = Ez (i,-2,1)  
           Vx (i,-3,1) = Vx (i,-2,1)  
           Vy (i,-3,1) = Vy (i,-2,1)  
           Vz (i,-3,1) = Vz (i,-2,1)  
           q  (i,-3,1) = q  (i,-2,1)   
           p  (i,-3,1) = p  (i,-2,1)
           rho(i,-3,1) = rho(i,-2,1)

           psi(i,-4,1) = psi(i,-3,1)
           phi(i,-4,1) = phi(i,-3,1) 
           Bx (i,-4,1) = Bx (i,-3,1)
           By (i,-4,1) = By (i,-3,1)
           Bz (i,-4,1) = Bz (i,-3,1)
           Ex (i,-4,1) = Ex (i,-3,1)  
           Ey (i,-4,1) = Ey (i,-3,1)  
           Ez (i,-4,1) = Ez (i,-3,1)  
           Vx (i,-4,1) = Vx (i,-3,1)  
           Vy (i,-4,1) = Vy (i,-3,1)  
           Vz (i,-4,1) = Vz (i,-3,1)  
           q  (i,-4,1) = q  (i,-3,1)   
           p  (i,-4,1) = p  (i,-3,1)
           rho(i,-4,1) = rho(i,-3,1)

           psi(i,-5,1) = psi(i,-4,1)
           phi(i,-5,1) = phi(i,-4,1) 
           Bx (i,-5,1) = Bx (i,-4,1)
           By (i,-5,1) = By (i,-4,1)
           Bz (i,-5,1) = Bz (i,-4,1)
           Ex (i,-5,1) = Ex (i,-4,1)  
           Ey (i,-5,1) = Ey (i,-4,1)  
           Ez (i,-5,1) = Ez (i,-4,1)  
           Vx (i,-5,1) = Vx (i,-4,1)  
           Vy (i,-5,1) = Vy (i,-4,1)  
           Vz (i,-5,1) = Vz (i,-4,1)  
           q  (i,-5,1) = q  (i,-4,1)   
           p  (i,-5,1) = p  (i,-4,1)
           rho(i,-5,1) = rho(i,-4,1)

           psi(i,-6,1) = psi(i,-5,1)
           phi(i,-6,1) = phi(i,-5,1) 
           Bx (i,-6,1) = Bx (i,-5,1)
           By (i,-6,1) = By (i,-5,1)
           Bz (i,-6,1) = Bz (i,-5,1)
           Ex (i,-6,1) = Ex (i,-5,1)  
           Ey (i,-6,1) = Ey (i,-5,1)  
           Ez (i,-6,1) = Ez (i,-5,1)  
           Vx (i,-6,1) = Vx (i,-5,1)  
           Vy (i,-6,1) = Vy (i,-5,1)  
           Vz (i,-6,1) = Vz (i,-5,1)  
           q  (i,-6,1) = q  (i,-5,1)   
           p  (i,-6,1) = p  (i,-5,1)
           rho(i,-6,1) = rho(i,-5,1)

           ! Right  Boundary
           ! -----------------

!!$           psi(i,jmax  ,1) = psi(i,jmax-1,1) 
!!$           phi(i,jmax  ,1) = phi(i,jmax-1,1) 
!!$           Bx (i,jmax  ,1) = Bx (i,jmax-1,1) 
!!$           By (i,jmax  ,1) = By (i,jmax-1,1) 
!!$           Bz (i,jmax  ,1) = Bz (i,jmax-1,1)  
!!$           Ex (i,jmax  ,1) = Ex (i,jmax-1,1)  
!!$           Ey (i,jmax  ,1) = Ey (i,jmax-1,1)  
!!$           Ez (i,jmax  ,1) = Ez (i,jmax-1,1)  
!!$           Vx (i,jmax  ,1) = Vx (i,jmax-1,1)  
!!$           Vy (i,jmax  ,1) = Vy (i,jmax-1,1)  
!!$           Vz (i,jmax  ,1) = Vz (i,jmax-1,1)  
!!$           q  (i,jmax  ,1) = q  (i,jmax-1,1)   
!!$           p  (i,jmax  ,1) = p  (i,jmax-1,1)   
!!$           rho(i,jmax  ,1) = rho(i,jmax-1,1)   

           psi(i,jmax+1,1) = psi(i,jmax  ,1) 
           phi(i,jmax+1,1) = phi(i,jmax  ,1) 
           Bx (i,jmax+1,1) = Bx (i,jmax  ,1) 
           By (i,jmax+1,1) = By (i,jmax  ,1) 
           Bz (i,jmax+1,1) = Bz (i,jmax  ,1)  
           Ex (i,jmax+1,1) = Ex (i,jmax  ,1)  
           Ey (i,jmax+1,1) = Ey (i,jmax  ,1)  
           Ez (i,jmax+1,1) = Ez (i,jmax  ,1)  
           Vx (i,jmax+1,1) = Vx (i,jmax  ,1)  
           Vy (i,jmax+1,1) = Vy (i,jmax  ,1)  
           Vz (i,jmax+1,1) = Vz (i,jmax  ,1)  
           q  (i,jmax+1,1) = q  (i,jmax  ,1)   
           p  (i,jmax+1,1) = p  (i,jmax  ,1)   
           rho(i,jmax+1,1) = rho(i,jmax  ,1)   

           psi(i,jmax+2,1) = psi(i,jmax+1,1) 
           phi(i,jmax+2,1) = phi(i,jmax+1,1) 
           Bx (i,jmax+2,1) = Bx (i,jmax+1,1) 
           By (i,jmax+2,1) = By (i,jmax+1,1) 
           Bz (i,jmax+2,1) = Bz (i,jmax+1,1)  
           Ex (i,jmax+2,1) = Ex (i,jmax+1,1)  
           Ey (i,jmax+2,1) = Ey (i,jmax+1,1)  
           Ez (i,jmax+2,1) = Ez (i,jmax+1,1)  
           Vx (i,jmax+2,1) = Vx (i,jmax+1,1)  
           Vy (i,jmax+2,1) = Vy (i,jmax+1,1)  
           Vz (i,jmax+2,1) = Vz (i,jmax+1,1)  
           q  (i,jmax+2,1) = q  (i,jmax+1,1)   
           p  (i,jmax+2,1) = p  (i,jmax+1,1)   
           rho(i,jmax+2,1) = rho(i,jmax+1,1)   

           psi(i,jmax+3,1) = psi(i,jmax+2,1) 
           phi(i,jmax+3,1) = phi(i,jmax+2,1) 
           Bx (i,jmax+3,1) = Bx (i,jmax+2,1) 
           By (i,jmax+3,1) = By (i,jmax+2,1) 
           Bz (i,jmax+3,1) = Bz (i,jmax+2,1)  
           Ex (i,jmax+3,1) = Ex (i,jmax+2,1)  
           Ey (i,jmax+3,1) = Ey (i,jmax+2,1)  
           Ez (i,jmax+3,1) = Ez (i,jmax+2,1)  
           Vx (i,jmax+3,1) = Vx (i,jmax+2,1)  
           Vy (i,jmax+3,1) = Vy (i,jmax+2,1)  
           Vz (i,jmax+3,1) = Vz (i,jmax+2,1)  
           q  (i,jmax+3,1) = q  (i,jmax+2,1)   
           p  (i,jmax+3,1) = p  (i,jmax+2,1)   
           rho(i,jmax+3,1) = rho(i,jmax+2,1)   

           psi(i,jmax+4,1) = psi(i,jmax+3,1) 
           phi(i,jmax+4,1) = phi(i,jmax+3,1) 
           Bx (i,jmax+4,1) = Bx (i,jmax+3,1) 
           By (i,jmax+4,1) = By (i,jmax+3,1) 
           Bz (i,jmax+4,1) = Bz (i,jmax+3,1)  
           Ex (i,jmax+4,1) = Ex (i,jmax+3,1)  
           Ey (i,jmax+4,1) = Ey (i,jmax+3,1)  
           Ez (i,jmax+4,1) = Ez (i,jmax+3,1)  
           Vx (i,jmax+4,1) = Vx (i,jmax+3,1)  
           Vy (i,jmax+4,1) = Vy (i,jmax+3,1)  
           Vz (i,jmax+4,1) = Vz (i,jmax+3,1)  
           q  (i,jmax+4,1) = q  (i,jmax+3,1)   
           p  (i,jmax+4,1) = p  (i,jmax+3,1)   
           rho(i,jmax+4,1) = rho(i,jmax+3,1)

           psi(i,jmax+5,1) = psi(i,jmax+4,1) 
           phi(i,jmax+5,1) = phi(i,jmax+4,1) 
           Bx (i,jmax+5,1) = Bx (i,jmax+4,1) 
           By (i,jmax+5,1) = By (i,jmax+4,1) 
           Bz (i,jmax+5,1) = Bz (i,jmax+4,1)  
           Ex (i,jmax+5,1) = Ex (i,jmax+4,1)  
           Ey (i,jmax+5,1) = Ey (i,jmax+4,1)  
           Ez (i,jmax+5,1) = Ez (i,jmax+4,1)  
           Vx (i,jmax+5,1) = Vx (i,jmax+4,1)  
           Vy (i,jmax+5,1) = Vy (i,jmax+4,1)  
           Vz (i,jmax+5,1) = Vz (i,jmax+4,1)  
           q  (i,jmax+5,1) = q  (i,jmax+4,1)   
           p  (i,jmax+5,1) = p  (i,jmax+4,1)   
           rho(i,jmax+5,1) = rho(i,jmax+4,1)

           psi(i,jmax+6,1) = psi(i,jmax+5,1) 
           phi(i,jmax+6,1) = phi(i,jmax+5,1) 
           Bx (i,jmax+6,1) = Bx (i,jmax+5,1) 
           By (i,jmax+6,1) = By (i,jmax+5,1) 
           Bz (i,jmax+6,1) = Bz (i,jmax+5,1)  
           Ex (i,jmax+6,1) = Ex (i,jmax+5,1)  
           Ey (i,jmax+6,1) = Ey (i,jmax+5,1)  
           Ez (i,jmax+6,1) = Ez (i,jmax+5,1)  
           Vx (i,jmax+6,1) = Vx (i,jmax+5,1)  
           Vy (i,jmax+6,1) = Vy (i,jmax+5,1)  
           Vz (i,jmax+6,1) = Vz (i,jmax+5,1)  
           q  (i,jmax+6,1) = q  (i,jmax+5,1)   
           p  (i,jmax+6,1) = p  (i,jmax+5,1)   
           rho(i,jmax+6,1) = rho(i,jmax+5,1)

!        end do
     end do

!$OMP END DO
 

    else if (DIM == 2 .and. BOUND == 2) then

       ! Periodic boundary 2D

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO 

    do j=-6,jmax+6
!       do k=1,kmax


           ! Left Boundary
           ! -----------------

!!$           psi( 0,j,1) = psi(imax-1,j,1) 
!!$           phi( 0,j,1) = phi(imax-1,j,1) 
!!$           Bx ( 0,j,1) = Bx (imax-1,j,1)  
!!$           By ( 0,j,1) = By (imax-1,j,1)  
!!$           Bz ( 0,j,1) = Bz (imax-1,j,1)  
!!$           Ex ( 0,j,1) = Ex (imax-1,j,1)  
!!$           Ey ( 0,j,1) = Ey (imax-1,j,1)  
!!$           Ez ( 0,j,1) = Ez (imax-1,j,1)  
!!$           Vx ( 0,j,1) = Vx (imax-1,j,1)  
!!$           Vy ( 0,j,1) = Vy (imax-1,j,1)  
!!$           Vz ( 0,j,1) = Vz (imax-1,j,1)  
!!$           q  ( 0,j,1) = q  (imax-1,j,1)   
!!$           p  ( 0,j,1) = p  (imax-1,j,1)   
!!$           rho( 0,j,1) = rho(imax-1,j,1)  

           psi(-1,j,1) = psi(imax-1,j,1) 
           phi(-1,j,1) = phi(imax-1,j,1) 
           Bx (-1,j,1) = Bx (imax-1,j,1)  
           By (-1,j,1) = By (imax-1,j,1)  
           Bz (-1,j,1) = Bz (imax-1,j,1)  
           Ex (-1,j,1) = Ex (imax-1,j,1)  
           Ey (-1,j,1) = Ey (imax-1,j,1)  
           Ez (-1,j,1) = Ez (imax-1,j,1)  
           Vx (-1,j,1) = Vx (imax-1,j,1)  
           Vy (-1,j,1) = Vy (imax-1,j,1)  
           Vz (-1,j,1) = Vz (imax-1,j,1)  
           q  (-1,j,1) = q  (imax-1,j,1)   
           p  (-1,j,1) = p  (imax-1,j,1)   
           rho(-1,j,1) = rho(imax-1,j,1)  


           psi(-2,j,1) = psi(imax-2,j,1) 
           phi(-2,j,1) = phi(imax-2,j,1) 
           Bx (-2,j,1) = Bx (imax-2,j,1)  
           By (-2,j,1) = By (imax-2,j,1)  
           Bz (-2,j,1) = Bz (imax-2,j,1)  
           Ex (-2,j,1) = Ex (imax-2,j,1)  
           Ey (-2,j,1) = Ey (imax-2,j,1)  
           Ez (-2,j,1) = Ez (imax-2,j,1)  
           Vx (-2,j,1) = Vx (imax-2,j,1)  
           Vy (-2,j,1) = Vy (imax-2,j,1)  
           Vz (-2,j,1) = Vz (imax-2,j,1)  
           q  (-2,j,1) = q  (imax-2,j,1)   
           p (-2,j,1)  = p  (imax-2,j,1)
           rho(-2,j,1) = rho(imax-2,j,1)


           psi(-3,j,1) = psi(imax-3,j,1) 
           phi(-3,j,1) = phi(imax-3,j,1) 
           Bx (-3,j,1) = Bx (imax-3,j,1)  
           By (-3,j,1) = By (imax-3,j,1)  
           Bz (-3,j,1) = Bz (imax-3,j,1)  
           Ex (-3,j,1) = Ex (imax-3,j,1)  
           Ey (-3,j,1) = Ey (imax-3,j,1)  
           Ez (-3,j,1) = Ez (imax-3,j,1)  
           Vx (-3,j,1) = Vx (imax-3,j,1)  
           Vy (-3,j,1) = Vy (imax-3,j,1)  
           Vz (-3,j,1) = Vz (imax-3,j,1)  
           q  (-3,j,1) = q  (imax-3,j,1)   
           p  (-3,j,1) = p  (imax-3,j,1)
           rho(-3,j,1) = rho(imax-3,j,1)

           psi(-4,j,1) = psi(imax-4,j,1) 
           phi(-4,j,1) = phi(imax-4,j,1) 
           Bx (-4,j,1) = Bx (imax-4,j,1)  
           By (-4,j,1) = By (imax-4,j,1)  
           Bz (-4,j,1) = Bz (imax-4,j,1)  
           Ex (-4,j,1) = Ex (imax-4,j,1)  
           Ey (-4,j,1) = Ey (imax-4,j,1)  
           Ez (-4,j,1) = Ez (imax-4,j,1)  
           Vx (-4,j,1) = Vx (imax-4,j,1)  
           Vy (-4,j,1) = Vy (imax-4,j,1)  
           Vz (-4,j,1) = Vz (imax-4,j,1)  
           q  (-4,j,1) = q  (imax-4,j,1)   
           p  (-4,j,1) = p  (imax-4,j,1)
           rho(-4,j,1) = rho(imax-4,j,1)

           psi(-5,j,1) = psi(imax-5,j,1) 
           phi(-5,j,1) = phi(imax-5,j,1) 
           Bx (-5,j,1) = Bx (imax-5,j,1)  
           By (-5,j,1) = By (imax-5,j,1)  
           Bz (-5,j,1) = Bz (imax-5,j,1)  
           Ex (-5,j,1) = Ex (imax-5,j,1)  
           Ey (-5,j,1) = Ey (imax-5,j,1)  
           Ez (-5,j,1) = Ez (imax-5,j,1)  
           Vx (-5,j,1) = Vx (imax-5,j,1)  
           Vy (-5,j,1) = Vy (imax-5,j,1)  
           Vz (-5,j,1) = Vz (imax-5,j,1)  
           q  (-5,j,1) = q  (imax-5,j,1)   
           p  (-5,j,1) = p  (imax-5,j,1)
           rho(-5,j,1) = rho(imax-5,j,1)

           psi(-6,j,1) = psi(imax-6,j,1) 
           phi(-6,j,1) = phi(imax-6,j,1) 
           Bx (-6,j,1) = Bx (imax-6,j,1)  
           By (-6,j,1) = By (imax-6,j,1)  
           Bz (-6,j,1) = Bz (imax-6,j,1)  
           Ex (-6,j,1) = Ex (imax-6,j,1)  
           Ey (-6,j,1) = Ey (imax-6,j,1)  
           Ez (-6,j,1) = Ez (imax-6,j,1)  
           Vx (-6,j,1) = Vx (imax-6,j,1)  
           Vy (-6,j,1) = Vy (imax-6,j,1)  
           Vz (-6,j,1) = Vz (imax-6,j,1)  
           q  (-6,j,1) = q  (imax-6,j,1)   
           p  (-6,j,1) = p  (imax-6,j,1)
           rho(-6,j,1) = rho(imax-6,j,1)


           ! Right  Boundary
           ! -----------------

!!$           psi(imax  ,j,1) = psi(1,j,1) 
!!$           phi(imax  ,j,1) = phi(1,j,1) 
!!$           Bx (imax  ,j,1) = Bx (1,j,1) 
!!$           By (imax  ,j,1) = By (1,j,1) 
!!$           Bz (imax  ,j,1) = Bz (1,j,1) 
!!$           Ex (imax  ,j,1) = Ex (1,j,1) 
!!$           Ey (imax  ,j,1) = Ey (1,j,1) 
!!$           Ez (imax  ,j,1) = Ez (1,j,1) 
!!$           Vx (imax  ,j,1) = Vx (1,j,1) 
!!$           Vy (imax  ,j,1) = Vy (1,j,1) 
!!$           Vz (imax  ,j,1) = Vz (1,j,1) 
!!$           q  (imax  ,j,1) = q  (1,j,1) 
!!$           p  (imax  ,j,1) = p  (1,j,1) 
!!$           rho(imax  ,j,1) = rho(1,j,1)

           psi(imax+1,j,1) = psi(1,j,1) 
           phi(imax+1,j,1) = phi(1,j,1) 
           Bx (imax+1,j,1) = Bx (1,j,1) 
           By (imax+1,j,1) = By (1,j,1) 
           Bz (imax+1,j,1) = Bz (1,j,1) 
           Ex (imax+1,j,1) = Ex (1,j,1) 
           Ey (imax+1,j,1) = Ey (1,j,1) 
           Ez (imax+1,j,1) = Ez (1,j,1) 
           Vx (imax+1,j,1) = Vx (1,j,1) 
           Vy (imax+1,j,1) = Vy (1,j,1) 
           Vz (imax+1,j,1) = Vz (1,j,1) 
           q  (imax+1,j,1) = q  (1,j,1) 
           p  (imax+1,j,1) = p  (1,j,1) 
           rho(imax+1,j,1) = rho(1,j,1) 

           psi(imax+2,j,1) = psi(2,j,1)
           phi(imax+2,j,1) = phi(2,j,1)
           Bx (imax+2,j,1) = Bx (2,j,1) 
           By (imax+2,j,1) = By (2,j,1) 
           Bz (imax+2,j,1) = Bz (2,j,1) 
           Ex (imax+2,j,1) = Ex (2,j,1) 
           Ey (imax+2,j,1) = Ey (2,j,1) 
           Ez (imax+2,j,1) = Ez (2,j,1) 
           Vx (imax+2,j,1) = Vx (2,j,1) 
           Vy (imax+2,j,1) = Vy (2,j,1) 
           Vz (imax+2,j,1) = Vz (2,j,1) 
           q  (imax+2,j,1) = q  (2,j,1) 
           p  (imax+2,j,1) = p  (2,j,1) 
           rho(imax+2,j,1) = rho(2,j,1) 

           psi(imax+3,j,1) = psi(3,j,1)
           phi(imax+3,j,1) = phi(3,j,1)
           Bx (imax+3,j,1) = Bx (3,j,1) 
           By (imax+3,j,1) = By (3,j,1) 
           Bz (imax+3,j,1) = Bz (3,j,1) 
           Ex (imax+3,j,1) = Ex (3,j,1) 
           Ey (imax+3,j,1) = Ey (3,j,1) 
           Ez (imax+3,j,1) = Ez (3,j,1) 
           Vx (imax+3,j,1) = Vx (3,j,1) 
           Vy (imax+3,j,1) = Vy (3,j,1) 
           Vz (imax+3,j,1) = Vz (3,j,1) 
           q  (imax+3,j,1) = q  (3,j,1) 
           p  (imax+3,j,1) = p  (3,j,1) 
           rho(imax+3,j,1) = rho(3,j,1) 

           psi(imax+4,j,1) = psi(4,j,1)
           phi(imax+4,j,1) = phi(4,j,1)
           Bx (imax+4,j,1) = Bx (4,j,1) 
           By (imax+4,j,1) = By (4,j,1) 
           Bz (imax+4,j,1) = Bz (4,j,1) 
           Ex (imax+4,j,1) = Ex (4,j,1) 
           Ey (imax+4,j,1) = Ey (4,j,1) 
           Ez (imax+4,j,1) = Ez (4,j,1) 
           Vx (imax+4,j,1) = Vx (4,j,1) 
           Vy (imax+4,j,1) = Vy (4,j,1) 
           Vz (imax+4,j,1) = Vz (4,j,1) 
           q  (imax+4,j,1) = q  (4,j,1) 
           p  (imax+4,j,1) = p  (4,j,1) 
           rho(imax+4,j,1) = rho(4,j,1)

           psi(imax+5,j,1) = psi(5,j,1)
           phi(imax+5,j,1) = phi(5,j,1)
           Bx (imax+5,j,1) = Bx (5,j,1) 
           By (imax+5,j,1) = By (5,j,1) 
           Bz (imax+5,j,1) = Bz (5,j,1) 
           Ex (imax+5,j,1) = Ex (5,j,1) 
           Ey (imax+5,j,1) = Ey (5,j,1) 
           Ez (imax+5,j,1) = Ez (5,j,1) 
           Vx (imax+5,j,1) = Vx (5,j,1) 
           Vy (imax+5,j,1) = Vy (5,j,1) 
           Vz (imax+5,j,1) = Vz (5,j,1) 
           q  (imax+5,j,1) = q  (5,j,1) 
           p  (imax+5,j,1) = p  (5,j,1) 
           rho(imax+5,j,1) = rho(5,j,1)

           psi(imax+6,j,1) = psi(6,j,1)
           phi(imax+6,j,1) = phi(6,j,1)
           Bx (imax+6,j,1) = Bx (6,j,1) 
           By (imax+6,j,1) = By (6,j,1) 
           Bz (imax+6,j,1) = Bz (6,j,1) 
           Ex (imax+6,j,1) = Ex (6,j,1) 
           Ey (imax+6,j,1) = Ey (6,j,1) 
           Ez (imax+6,j,1) = Ez (6,j,1) 
           Vx (imax+6,j,1) = Vx (6,j,1) 
           Vy (imax+6,j,1) = Vy (6,j,1) 
           Vz (imax+6,j,1) = Vz (6,j,1) 
           q  (imax+6,j,1) = q  (6,j,1) 
           p  (imax+6,j,1) = p  (6,j,1) 
           rho(imax+6,j,1) = rho(6,j,1) 

!        end do
     end do

!$OMP END DO 


  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO  

    do i=-6,imax+6
!       do k=1,kmax

 


           ! Left Boundary
           ! -----------------

!!$           psi(i, 0,1) = psi(i,jmax-1,1)
!!$           phi(i, 0,1) = phi(i,jmax-1,1) 
!!$           Bx (i, 0,1) = Bx (i,jmax-1,1)
!!$           By (i, 0,1) = By (i,jmax-1,1)
!!$           Bz (i, 0,1) = Bz (i,jmax-1,1)
!!$           Ex (i, 0,1) = Ex (i,jmax-1,1)  
!!$           Ey (i, 0,1) = Ey (i,jmax-1,1)  
!!$           Ez (i, 0,1) = Ez (i,jmax-1,1)  
!!$           Vx (i, 0,1) = Vx (i,jmax-1,1)  
!!$           Vy (i, 0,1) = Vy (i,jmax-1,1)  
!!$           Vz (i, 0,1) = Vz (i,jmax-1,1)  
!!$           q  (i, 0,1) = q  (i,jmax-1,1)   
!!$           p  (i, 0,1) = p  (i,jmax-1,1)
!!$           rho(i, 0,1) = rho(i,jmax-1,1)

           psi(i,-1,1) = psi(i,jmax-1,1)
           phi(i,-1,1) = phi(i,jmax-1,1) 
           Bx (i,-1,1) = Bx (i,jmax-1,1)
           By (i,-1,1) = By (i,jmax-1,1)
           Bz (i,-1,1) = Bz (i,jmax-1,1)
           Ex (i,-1,1) = Ex (i,jmax-1,1)  
           Ey (i,-1,1) = Ey (i,jmax-1,1)  
           Ez (i,-1,1) = Ez (i,jmax-1,1)  
           Vx (i,-1,1) = Vx (i,jmax-1,1)  
           Vy (i,-1,1) = Vy (i,jmax-1,1)  
           Vz (i,-1,1) = Vz (i,jmax-1,1)  
           q  (i,-1,1) = q  (i,jmax-1,1)   
           p  (i,-1,1) = p  (i,jmax-1,1)
           rho(i,-1,1) = rho(i,jmax-1,1)

           psi(i,-2,1) = psi(i,jmax-2,1)
           phi(i,-2,1) = phi(i,jmax-2,1) 
           Bx (i,-2,1) = Bx (i,jmax-2,1)
           By (i,-2,1) = By (i,jmax-2,1)
           Bz (i,-2,1) = Bz (i,jmax-2,1)
           Ex (i,-2,1) = Ex (i,jmax-2,1)  
           Ey (i,-2,1) = Ey (i,jmax-2,1)  
           Ez (i,-2,1) = Ez (i,jmax-2,1)  
           Vx (i,-2,1) = Vx (i,jmax-2,1)  
           Vy (i,-2,1) = Vy (i,jmax-2,1)  
           Vz (i,-2,1) = Vz (i,jmax-2,1)  
           q  (i,-2,1) = q  (i,jmax-2,1)   
           p  (i,-2,1) = p  (i,jmax-2,1)
           rho(i,-2,1) = rho(i,jmax-2,1)

           psi(i,-3,1) = psi(i,jmax-3,1)
           phi(i,-3,1) = phi(i,jmax-3,1) 
           Bx (i,-3,1) = Bx (i,jmax-3,1)
           By (i,-3,1) = By (i,jmax-3,1)
           Bz (i,-3,1) = Bz (i,jmax-3,1)
           Ex (i,-3,1) = Ex (i,jmax-3,1)  
           Ey (i,-3,1) = Ey (i,jmax-3,1)  
           Ez (i,-3,1) = Ez (i,jmax-3,1)  
           Vx (i,-3,1) = Vx (i,jmax-3,1)  
           Vy (i,-3,1) = Vy (i,jmax-3,1)  
           Vz (i,-3,1) = Vz (i,jmax-3,1)  
           q  (i,-3,1) = q  (i,jmax-3,1)   
           p  (i,-3,1) = p  (i,jmax-3,1)
           rho(i,-3,1) = rho(i,jmax-3,1)

           psi(i,-4,1) = psi(i,jmax-4,1)
           phi(i,-4,1) = phi(i,jmax-4,1) 
           Bx (i,-4,1) = Bx (i,jmax-4,1)
           By (i,-4,1) = By (i,jmax-4,1)
           Bz (i,-4,1) = Bz (i,jmax-4,1)
           Ex (i,-4,1) = Ex (i,jmax-4,1)  
           Ey (i,-4,1) = Ey (i,jmax-4,1)  
           Ez (i,-4,1) = Ez (i,jmax-4,1)  
           Vx (i,-4,1) = Vx (i,jmax-4,1)  
           Vy (i,-4,1) = Vy (i,jmax-4,1)  
           Vz (i,-4,1) = Vz (i,jmax-4,1)  
           q  (i,-4,1) = q  (i,jmax-4,1)   
           p  (i,-4,1) = p  (i,jmax-4,1)
           rho(i,-4,1) = rho(i,jmax-4,1)

           psi(i,-5,1) = psi(i,jmax-5,1)
           phi(i,-5,1) = phi(i,jmax-5,1) 
           Bx (i,-5,1) = Bx (i,jmax-5,1)
           By (i,-5,1) = By (i,jmax-5,1)
           Bz (i,-5,1) = Bz (i,jmax-5,1)
           Ex (i,-5,1) = Ex (i,jmax-5,1)  
           Ey (i,-5,1) = Ey (i,jmax-5,1)  
           Ez (i,-5,1) = Ez (i,jmax-5,1)  
           Vx (i,-5,1) = Vx (i,jmax-5,1)  
           Vy (i,-5,1) = Vy (i,jmax-5,1)  
           Vz (i,-5,1) = Vz (i,jmax-5,1)  
           q  (i,-5,1) = q  (i,jmax-5,1)   
           p  (i,-5,1) = p  (i,jmax-5,1)
           rho(i,-5,1) = rho(i,jmax-5,1)

           psi(i,-6,1) = psi(i,jmax-6,1)
           phi(i,-6,1) = phi(i,jmax-6,1) 
           Bx (i,-6,1) = Bx (i,jmax-6,1)
           By (i,-6,1) = By (i,jmax-6,1)
           Bz (i,-6,1) = Bz (i,jmax-6,1)
           Ex (i,-6,1) = Ex (i,jmax-6,1)  
           Ey (i,-6,1) = Ey (i,jmax-6,1)  
           Ez (i,-6,1) = Ez (i,jmax-6,1)  
           Vx (i,-6,1) = Vx (i,jmax-6,1)  
           Vy (i,-6,1) = Vy (i,jmax-6,1)  
           Vz (i,-6,1) = Vz (i,jmax-6,1)  
           q  (i,-6,1) = q  (i,jmax-6,1)   
           p  (i,-6,1) = p  (i,jmax-6,1)
           rho(i,-6,1) = rho(i,jmax-6,1)


           ! Right  Boundary
           ! -----------------

!!$           psi(i,jmax  ,1) = psi(i,1,1) 
!!$           phi(i,jmax  ,1) = phi(i,1,1) 
!!$           Bx (i,jmax  ,1) = Bx (i,1,1) 
!!$           By (i,jmax  ,1) = By (i,1,1) 
!!$           Bz (i,jmax  ,1) = Bz (i,1,1)  
!!$           Ex (i,jmax  ,1) = Ex (i,1,1)  
!!$           Ey (i,jmax  ,1) = Ey (i,1,1)  
!!$           Ez (i,jmax  ,1) = Ez (i,1,1)  
!!$           Vx (i,jmax  ,1) = Vx (i,1,1)  
!!$           Vy (i,jmax  ,1) = Vy (i,1,1)  
!!$           Vz (i,jmax  ,1) = Vz (i,1,1)  
!!$           q  (i,jmax  ,1) = q  (i,1,1)   
!!$           p  (i,jmax  ,1) = p  (i,1,1)   
!!$           rho(i,jmax  ,1) = rho(i,1,1)   

           psi(i,jmax+1,1) = psi(i,1,1) 
           phi(i,jmax+1,1) = phi(i,1,1) 
           Bx (i,jmax+1,1) = Bx (i,1,1) 
           By (i,jmax+1,1) = By (i,1,1) 
           Bz (i,jmax+1,1) = Bz (i,1,1)  
           Ex (i,jmax+1,1) = Ex (i,1,1)  
           Ey (i,jmax+1,1) = Ey (i,1,1)  
           Ez (i,jmax+1,1) = Ez (i,1,1)  
           Vx (i,jmax+1,1) = Vx (i,1,1)  
           Vy (i,jmax+1,1) = Vy (i,1,1)  
           Vz (i,jmax+1,1) = Vz (i,1,1)  
           q  (i,jmax+1,1) = q  (i,1,1)   
           p  (i,jmax+1,1) = p  (i,1,1)   
           rho(i,jmax+1,1) = rho(i,1,1)   

           psi(i,jmax+2,1) = psi(i,2,1) 
           phi(i,jmax+2,1) = phi(i,2,1) 
           Bx (i,jmax+2,1) = Bx (i,2,1) 
           By (i,jmax+2,1) = By (i,2,1) 
           Bz (i,jmax+2,1) = Bz (i,2,1)  
           Ex (i,jmax+2,1) = Ex (i,2,1)  
           Ey (i,jmax+2,1) = Ey (i,2,1)  
           Ez (i,jmax+2,1) = Ez (i,2,1)  
           Vx (i,jmax+2,1) = Vx (i,2,1)  
           Vy (i,jmax+2,1) = Vy (i,2,1)  
           Vz (i,jmax+2,1) = Vz (i,2,1)  
           q  (i,jmax+2,1) = q  (i,2,1)   
           p  (i,jmax+2,1) = p  (i,2,1)   
           rho(i,jmax+2,1) = rho(i,2,1)   

           psi(i,jmax+3,1) = psi(i,3,1) 
           phi(i,jmax+3,1) = phi(i,3,1) 
           Bx (i,jmax+3,1) = Bx (i,3,1) 
           By (i,jmax+3,1) = By (i,3,1) 
           Bz (i,jmax+3,1) = Bz (i,3,1)  
           Ex (i,jmax+3,1) = Ex (i,3,1)  
           Ey (i,jmax+3,1) = Ey (i,3,1)  
           Ez (i,jmax+3,1) = Ez (i,3,1)  
           Vx (i,jmax+3,1) = Vx (i,3,1)  
           Vy (i,jmax+3,1) = Vy (i,3,1)  
           Vz (i,jmax+3,1) = Vz (i,3,1)  
           q  (i,jmax+3,1) = q  (i,3,1)   
           p  (i,jmax+3,1) = p  (i,3,1)   
           rho(i,jmax+3,1) = rho(i,3,1)   

           psi(i,jmax+4,1) = psi(i,4,1) 
           phi(i,jmax+4,1) = phi(i,4,1) 
           Bx (i,jmax+4,1) = Bx (i,4,1) 
           By (i,jmax+4,1) = By (i,4,1) 
           Bz (i,jmax+4,1) = Bz (i,4,1)  
           Ex (i,jmax+4,1) = Ex (i,4,1)  
           Ey (i,jmax+4,1) = Ey (i,4,1)  
           Ez (i,jmax+4,1) = Ez (i,4,1)  
           Vx (i,jmax+4,1) = Vx (i,4,1)  
           Vy (i,jmax+4,1) = Vy (i,4,1)  
           Vz (i,jmax+4,1) = Vz (i,4,1)  
           q  (i,jmax+4,1) = q  (i,4,1)   
           p  (i,jmax+4,1) = p  (i,4,1)   
           rho(i,jmax+4,1) = rho(i,4,1)

           psi(i,jmax+5,1) = psi(i,5,1) 
           phi(i,jmax+5,1) = phi(i,5,1) 
           Bx (i,jmax+5,1) = Bx (i,5,1) 
           By (i,jmax+5,1) = By (i,5,1) 
           Bz (i,jmax+5,1) = Bz (i,5,1)  
           Ex (i,jmax+5,1) = Ex (i,5,1)  
           Ey (i,jmax+5,1) = Ey (i,5,1)  
           Ez (i,jmax+5,1) = Ez (i,5,1)  
           Vx (i,jmax+5,1) = Vx (i,5,1)  
           Vy (i,jmax+5,1) = Vy (i,5,1)  
           Vz (i,jmax+5,1) = Vz (i,5,1)  
           q  (i,jmax+5,1) = q  (i,5,1)   
           p  (i,jmax+5,1) = p  (i,5,1)   
           rho(i,jmax+5,1) = rho(i,5,1)

           psi(i,jmax+6,1) = psi(i,6,1) 
           phi(i,jmax+6,1) = phi(i,6,1) 
           Bx (i,jmax+6,1) = Bx (i,6,1) 
           By (i,jmax+6,1) = By (i,6,1) 
           Bz (i,jmax+6,1) = Bz (i,6,1)  
           Ex (i,jmax+6,1) = Ex (i,6,1)  
           Ey (i,jmax+6,1) = Ey (i,6,1)  
           Ez (i,jmax+6,1) = Ez (i,6,1)  
           Vx (i,jmax+6,1) = Vx (i,6,1)  
           Vy (i,jmax+6,1) = Vy (i,6,1)  
           Vz (i,jmax+6,1) = Vz (i,6,1)  
           q  (i,jmax+6,1) = q  (i,6,1)   
           p  (i,jmax+6,1) = p  (i,6,1)   
           rho(i,jmax+6,1) = rho(i,6,1) 

 !        end do
     end do

!$OMP END DO


    else if (DIM == 2 .and. BOUND == 3) then

       ! Reconnection boundary TAKAHASHI

       ! Copy boundary 2D

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO  
       
    do j=-6,jmax+6
       !       do k=1,kmax

       posx   = - 0.5d0 * Lx 
       posy   = - 0.5d0 * Ly + j * Dely

       r      =  sqrt(1.d0 - (posy**2/Ly**2))
       

           ! Left Boundary
           ! -----------------

           psi( 0,j,1) = 0.d0 !psi(1,j,1) 
           phi( 0,j,1) = 0.d0 !phi(1,j,1) 
           Bx ( 0,j,1) = (B0*posy*delta_rec/Ly**2) * (1.d0+(log(abs(cosh(posx/delta_rec)))/r))
           By ( 0,j,1) = B0 * r * tanh(posx/delta_rec) 
           Bz ( 0,j,1) = 0.d0 !Bz (1,j,1)  
           Ex ( 0,j,1) = Ex (1,j,1)  
           Ey ( 0,j,1) = Ey (1,j,1)  
           Ez ( 0,j,1) = Ez (1,j,1)  
           Vx ( 0,j,1) = Vx (1,j,1)  
           Vy ( 0,j,1) = Vy (1,j,1)  
           Vz ( 0,j,1) = Vz (1,j,1)  
           q  ( 0,j,1) = 0.d0 !q  (1,j,1)   
           p  ( 0,j,1) = P0 * r**2 / (cosh(posx/delta_rec))**2 + P_up 
           rho( 0,j,1) = rho_up

       
           psi(-1,j,1) = psi(0,j,1) 
           phi(-1,j,1) = phi(0,j,1) 
           Bx (-1,j,1) = Bx (0,j,1)  
           By (-1,j,1) = By (0,j,1)  
           Bz (-1,j,1) = Bz (0,j,1)  
           Ex (-1,j,1) = Ex (0,j,1)  
           Ey (-1,j,1) = Ey (0,j,1)  
           Ez (-1,j,1) = Ez (0,j,1)  
           Vx (-1,j,1) = Vx (0,j,1)  
           Vy (-1,j,1) = Vy (0,j,1)  
           Vz (-1,j,1) = Vz (0,j,1)  
           q  (-1,j,1) = q  (0,j,1)   
           p  (-1,j,1) = p  (0,j,1)   
           rho(-1,j,1) = rho(0,j,1)  


           psi(-2,j,1) = psi(-1,j,1) 
           phi(-2,j,1) = phi(-1,j,1) 
           Bx (-2,j,1) = Bx (-1,j,1)  
           By (-2,j,1) = By (-1,j,1)  
           Bz (-2,j,1) = Bz (-1,j,1)  
           Ex (-2,j,1) = Ex (-1,j,1)  
           Ey (-2,j,1) = Ey (-1,j,1)  
           Ez (-2,j,1) = Ez (-1,j,1)  
           Vx (-2,j,1) = Vx (-1,j,1)  
           Vy (-2,j,1) = Vy (-1,j,1)  
           Vz (-2,j,1) = Vz (-1,j,1)  
           q  (-2,j,1) = q  (-1,j,1)   
           p  (-2,j,1) = p  (-1,j,1)
           rho(-2,j,1) = rho(-1,j,1)

           psi(-3,j,1) = psi(-2,j,1) 
           phi(-3,j,1) = phi(-2,j,1) 
           Bx (-3,j,1) = Bx (-2,j,1)  
           By (-3,j,1) = By (-2,j,1)  
           Bz (-3,j,1) = Bz (-2,j,1)  
           Ex (-3,j,1) = Ex (-2,j,1)  
           Ey (-3,j,1) = Ey (-2,j,1)  
           Ez (-3,j,1) = Ez (-2,j,1)  
           Vx (-3,j,1) = Vx (-2,j,1)  
           Vy (-3,j,1) = Vy (-2,j,1)  
           Vz (-3,j,1) = Vz (-2,j,1)  
           q  (-3,j,1) = q  (-2,j,1)   
           p  (-3,j,1) = p  (-2,j,1)
           rho(-3,j,1) = rho(-2,j,1)

           psi(-4,j,1) = psi(-3,j,1) 
           phi(-4,j,1) = phi(-3,j,1) 
           Bx (-4,j,1) = Bx (-3,j,1)  
           By (-4,j,1) = By (-3,j,1)  
           Bz (-4,j,1) = Bz (-3,j,1)  
           Ex (-4,j,1) = Ex (-3,j,1)  
           Ey (-4,j,1) = Ey (-3,j,1)  
           Ez (-4,j,1) = Ez (-3,j,1)  
           Vx (-4,j,1) = Vx (-3,j,1)  
           Vy (-4,j,1) = Vy (-3,j,1)  
           Vz (-4,j,1) = Vz (-3,j,1)  
           q  (-4,j,1) = q  (-3,j,1)   
           p  (-4,j,1) = p  (-3,j,1)
           rho(-4,j,1) = rho(-3,j,1)

           psi(-5,j,1) = psi(-4,j,1) 
           phi(-5,j,1) = phi(-4,j,1) 
           Bx (-5,j,1) = Bx (-4,j,1)  
           By (-5,j,1) = By (-4,j,1)  
           Bz (-5,j,1) = Bz (-4,j,1)  
           Ex (-5,j,1) = Ex (-4,j,1)  
           Ey (-5,j,1) = Ey (-4,j,1)  
           Ez (-5,j,1) = Ez (-4,j,1)  
           Vx (-5,j,1) = Vx (-4,j,1)  
           Vy (-5,j,1) = Vy (-4,j,1)  
           Vz (-5,j,1) = Vz (-4,j,1)  
           q  (-5,j,1) = q  (-4,j,1)   
           p  (-5,j,1) = p  (-4,j,1)
           rho(-5,j,1) = rho(-4,j,1)

           psi(-6,j,1) = psi(-5,j,1) 
           phi(-6,j,1) = phi(-5,j,1) 
           Bx (-6,j,1) = Bx (-5,j,1)  
           By (-6,j,1) = By (-5,j,1)  
           Bz (-6,j,1) = Bz (-5,j,1)  
           Ex (-6,j,1) = Ex (-5,j,1)  
           Ey (-6,j,1) = Ey (-5,j,1)  
           Ez (-6,j,1) = Ez (-5,j,1)  
           Vx (-6,j,1) = Vx (-5,j,1)  
           Vy (-6,j,1) = Vy (-5,j,1)  
           Vz (-6,j,1) = Vz (-5,j,1)  
           q  (-6,j,1) = q  (-5,j,1)   
           p  (-6,j,1) = p  (-5,j,1)
           rho(-6,j,1) = rho(-5,j,1)

           ! Right  Boundary
           ! -----------------


           posx   = - 0.5d0 * Lx + imax * Delx
           posy   = - 0.5d0 * Ly + j    * Dely

           r      =  sqrt(1.d0 - (posy**2/Ly**2))


           psi(imax  ,j,1) = 0.d0 !psi(imax-1,j,1) 
           phi(imax  ,j,1) = 0.d0 !phi(imax-1,j,1) 
           Bx (imax  ,j,1) = (B0*posy*delta_rec/Ly**2) * (1.d0+(log(abs(cosh(posx/delta_rec)))/r))
           By (imax  ,j,1) = B0 * r * tanh(posx/delta_rec) 
           Bz (imax  ,j,1) = 0.d0 !Bz (imax-1,j,1) 
           Ex (imax  ,j,1) = Ex (imax-1,j,1) 
           Ey (imax  ,j,1) = Ey (imax-1,j,1) 
           Ez (imax  ,j,1) = Ez (imax-1,j,1)
           Vx (imax  ,j,1) = Vx (imax-1,j,1) 
           Vy (imax  ,j,1) = Vy (imax-1,j,1) 
           Vz (imax  ,j,1) = Vz (imax-1,j,1) 
           q  (imax  ,j,1) = 0.d0 !q  (imax-1,j,1) 
           p  (imax  ,j,1) = P0 * r**2 / (cosh(posx/delta_rec))**2 + P_up 
           rho(imax  ,j,1) = rho_up

           psi(imax+1,j,1) = psi(imax  ,j,1) 
           phi(imax+1,j,1) = phi(imax  ,j,1) 
           Bx (imax+1,j,1) = Bx (imax  ,j,1) 
           By (imax+1,j,1) = By (imax  ,j,1) 
           Bz (imax+1,j,1) = Bz (imax  ,j,1) 
           Ex (imax+1,j,1) = Ex (imax  ,j,1) 
           Ey (imax+1,j,1) = Ey (imax  ,j,1) 
           Ez (imax+1,j,1) = Ez (imax  ,j,1) 
           Vx (imax+1,j,1) = Vx (imax  ,j,1) 
           Vy (imax+1,j,1) = Vy (imax  ,j,1) 
           Vz (imax+1,j,1) = Vz (imax  ,j,1) 
           q  (imax+1,j,1) = q  (imax  ,j,1) 
           p  (imax+1,j,1) = p  (imax  ,j,1) 
           rho(imax+1,j,1) = rho(imax  ,j,1) 

           psi(imax+2,j,1) = psi(imax+1,j,1)
           phi(imax+2,j,1) = phi(imax+1,j,1)
           Bx (imax+2,j,1) = Bx (imax+1,j,1) 
           By (imax+2,j,1) = By (imax+1,j,1) 
           Bz (imax+2,j,1) = Bz (imax+1,j,1) 
           Ex (imax+2,j,1) = Ex (imax+1,j,1) 
           Ey (imax+2,j,1) = Ey (imax+1,j,1) 
           Ez (imax+2,j,1) = Ez (imax+1,j,1) 
           Vx (imax+2,j,1) = Vx (imax+1,j,1) 
           Vy (imax+2,j,1) = Vy (imax+1,j,1) 
           Vz (imax+2,j,1) = Vz (imax+1,j,1) 
           q  (imax+2,j,1) = q  (imax+1,j,1) 
           p  (imax+2,j,1) = p  (imax+1,j,1) 
           rho(imax+2,j,1) = rho(imax+1,j,1) 

           psi(imax+3,j,1) = psi(imax+2,j,1)
           phi(imax+3,j,1) = phi(imax+2,j,1)
           Bx (imax+3,j,1) = Bx (imax+2,j,1) 
           By (imax+3,j,1) = By (imax+2,j,1) 
           Bz (imax+3,j,1) = Bz (imax+2,j,1) 
           Ex (imax+3,j,1) = Ex (imax+2,j,1) 
           Ey (imax+3,j,1) = Ey (imax+2,j,1) 
           Ez (imax+3,j,1) = Ez (imax+2,j,1) 
           Vx (imax+3,j,1) = Vx (imax+2,j,1) 
           Vy (imax+3,j,1) = Vy (imax+2,j,1) 
           Vz (imax+3,j,1) = Vz (imax+2,j,1) 
           q  (imax+3,j,1) = q  (imax+2,j,1) 
           p  (imax+3,j,1) = p  (imax+2,j,1) 
           rho(imax+3,j,1) = rho(imax+2,j,1) 

           psi(imax+4,j,1) = psi(imax+3,j,1)
           phi(imax+4,j,1) = phi(imax+3,j,1)
           Bx (imax+4,j,1) = Bx (imax+3,j,1) 
           By (imax+4,j,1) = By (imax+3,j,1) 
           Bz (imax+4,j,1) = Bz (imax+3,j,1) 
           Ex (imax+4,j,1) = Ex (imax+3,j,1) 
           Ey (imax+4,j,1) = Ey (imax+3,j,1) 
           Ez (imax+4,j,1) = Ez (imax+3,j,1) 
           Vx (imax+4,j,1) = Vx (imax+3,j,1) 
           Vy (imax+4,j,1) = Vy (imax+3,j,1) 
           Vz (imax+4,j,1) = Vz (imax+3,j,1) 
           q  (imax+4,j,1) = q  (imax+3,j,1) 
           p  (imax+4,j,1) = p  (imax+3,j,1) 
           rho(imax+4,j,1) = rho(imax+3,j,1)

           psi(imax+5,j,1) = psi(imax+4,j,1)
           phi(imax+5,j,1) = phi(imax+4,j,1)
           Bx (imax+5,j,1) = Bx (imax+4,j,1) 
           By (imax+5,j,1) = By (imax+4,j,1) 
           Bz (imax+5,j,1) = Bz (imax+4,j,1) 
           Ex (imax+5,j,1) = Ex (imax+4,j,1) 
           Ey (imax+5,j,1) = Ey (imax+4,j,1) 
           Ez (imax+5,j,1) = Ez (imax+4,j,1) 
           Vx (imax+5,j,1) = Vx (imax+4,j,1) 
           Vy (imax+5,j,1) = Vy (imax+4,j,1) 
           Vz (imax+5,j,1) = Vz (imax+4,j,1) 
           q  (imax+5,j,1) = q  (imax+4,j,1) 
           p  (imax+5,j,1) = p  (imax+4,j,1) 
           rho(imax+5,j,1) = rho(imax+4,j,1)

           psi(imax+6,j,1) = psi(imax+5,j,1)
           phi(imax+6,j,1) = phi(imax+5,j,1)
           Bx (imax+6,j,1) = Bx (imax+5,j,1) 
           By (imax+6,j,1) = By (imax+5,j,1) 
           Bz (imax+6,j,1) = Bz (imax+5,j,1) 
           Ex (imax+6,j,1) = Ex (imax+5,j,1) 
           Ey (imax+6,j,1) = Ey (imax+5,j,1) 
           Ez (imax+6,j,1) = Ez (imax+5,j,1) 
           Vx (imax+6,j,1) = Vx (imax+5,j,1) 
           Vy (imax+6,j,1) = Vy (imax+5,j,1) 
           Vz (imax+6,j,1) = Vz (imax+5,j,1) 
           q  (imax+6,j,1) = q  (imax+5,j,1) 
           p  (imax+6,j,1) = p  (imax+5,j,1) 
           rho(imax+6,j,1) = rho(imax+5,j,1) 


!        end do
     end do
 
!$OMP END DO


  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 

    do i=-6,imax+6
       !       do k=1,kmax

           psi(i,0,1) = psi(i, 1,1)
           phi(i,0,1) = phi(i, 1,1) 
           Bx (i,0,1) = Bx (i, 1,1)
           By (i,0,1) = By (i, 1,1)
           Bz (i,0,1) = Bz (i, 1,1)
           Ex (i,0,1) = Ex (i, 1,1)  
           Ey (i,0,1) = Ey (i, 1,1)  
           Ez (i,0,1) = Ez (i, 1,1)  
           Vx (i,0,1) = Vx (i, 1,1)  
           Vy (i,0,1) = Vy (i, 1,1)  
           Vz (i,0,1) = Vz (i, 1,1)  
           q  (i,0,1) = q  (i, 1,1)   
           p  (i,0,1) = p  (i, 1,1)
           rho(i,0,1) = rho(i, 1,1)

           psi(i,-1,1) = psi(i, 0,1)
           phi(i,-1,1) = phi(i, 0,1) 
           Bx (i,-1,1) = Bx (i, 0,1)
           By (i,-1,1) = By (i, 0,1)
           Bz (i,-1,1) = Bz (i, 0,1)
           Ex (i,-1,1) = Ex (i, 0,1)  
           Ey (i,-1,1) = Ey (i, 0,1)  
           Ez (i,-1,1) = Ez (i, 0,1)  
           Vx (i,-1,1) = Vx (i, 0,1)  
           Vy (i,-1,1) = Vy (i, 0,1)  
           Vz (i,-1,1) = Vz (i, 0,1)  
           q  (i,-1,1) = q  (i, 0,1)   
           p  (i,-1,1) = p  (i, 0,1)
           rho(i,-1,1) = rho(i, 0,1)

           psi(i,-2,1) = psi(i,-1,1)
           phi(i,-2,1) = phi(i,-1,1) 
           Bx (i,-2,1) = Bx (i,-1,1)
           By (i,-2,1) = By (i,-1,1)
           Bz (i,-2,1) = Bz (i,-1,1)
           Ex (i,-2,1) = Ex (i,-1,1)  
           Ey (i,-2,1) = Ey (i,-1,1)  
           Ez (i,-2,1) = Ez (i,-1,1)  
           Vx (i,-2,1) = Vx (i,-1,1)  
           Vy (i,-2,1) = Vy (i,-1,1)  
           Vz (i,-2,1) = Vz (i,-1,1)  
           q  (i,-2,1) = q  (i,-1,1)   
           p  (i,-2,1) = p  (i,-1,1)
           rho(i,-2,1) = rho(i,-1,1)

           psi(i,-3,1) = psi(i,-2,1)
           phi(i,-3,1) = phi(i,-2,1) 
           Bx (i,-3,1) = Bx (i,-2,1)
           By (i,-3,1) = By (i,-2,1)
           Bz (i,-3,1) = Bz (i,-2,1)
           Ex (i,-3,1) = Ex (i,-2,1)  
           Ey (i,-3,1) = Ey (i,-2,1)  
           Ez (i,-3,1) = Ez (i,-2,1)  
           Vx (i,-3,1) = Vx (i,-2,1)  
           Vy (i,-3,1) = Vy (i,-2,1)  
           Vz (i,-3,1) = Vz (i,-2,1)  
           q  (i,-3,1) = q  (i,-2,1)   
           p  (i,-3,1) = p  (i,-2,1)
           rho(i,-3,1) = rho(i,-2,1)

           psi(i,-4,1) = psi(i,-3,1)
           phi(i,-4,1) = phi(i,-3,1) 
           Bx (i,-4,1) = Bx (i,-3,1)
           By (i,-4,1) = By (i,-3,1)
           Bz (i,-4,1) = Bz (i,-3,1)
           Ex (i,-4,1) = Ex (i,-3,1)  
           Ey (i,-4,1) = Ey (i,-3,1)  
           Ez (i,-4,1) = Ez (i,-3,1)  
           Vx (i,-4,1) = Vx (i,-3,1)  
           Vy (i,-4,1) = Vy (i,-3,1)  
           Vz (i,-4,1) = Vz (i,-3,1)  
           q  (i,-4,1) = q  (i,-3,1)   
           p  (i,-4,1) = p  (i,-3,1)
           rho(i,-4,1) = rho(i,-3,1)

           psi(i,-5,1) = psi(i,-4,1)
           phi(i,-5,1) = phi(i,-4,1) 
           Bx (i,-5,1) = Bx (i,-4,1)
           By (i,-5,1) = By (i,-4,1)
           Bz (i,-5,1) = Bz (i,-4,1)
           Ex (i,-5,1) = Ex (i,-4,1)  
           Ey (i,-5,1) = Ey (i,-4,1)  
           Ez (i,-5,1) = Ez (i,-4,1)  
           Vx (i,-5,1) = Vx (i,-4,1)  
           Vy (i,-5,1) = Vy (i,-4,1)  
           Vz (i,-5,1) = Vz (i,-4,1)  
           q  (i,-5,1) = q  (i,-4,1)   
           p  (i,-5,1) = p  (i,-4,1)
           rho(i,-5,1) = rho(i,-4,1)

           psi(i,-6,1) = psi(i,-5,1)
           phi(i,-6,1) = phi(i,-5,1) 
           Bx (i,-6,1) = Bx (i,-5,1)
           By (i,-6,1) = By (i,-5,1)
           Bz (i,-6,1) = Bz (i,-5,1)
           Ex (i,-6,1) = Ex (i,-5,1)  
           Ey (i,-6,1) = Ey (i,-5,1)  
           Ez (i,-6,1) = Ez (i,-5,1)  
           Vx (i,-6,1) = Vx (i,-5,1)  
           Vy (i,-6,1) = Vy (i,-5,1)  
           Vz (i,-6,1) = Vz (i,-5,1)  
           q  (i,-6,1) = q  (i,-5,1)   
           p  (i,-6,1) = p  (i,-5,1)
           rho(i,-6,1) = rho(i,-5,1)

          
           ! Right  Boundary
           ! -----------------


           psi(i,jmax  ,1) = psi(i,jmax-1,1) 
           phi(i,jmax  ,1) = phi(i,jmax-1,1) 
           Bx (i,jmax  ,1) = Bx (i,jmax-1,1) 
           By (i,jmax  ,1) = By (i,jmax-1,1) 
           Bz (i,jmax  ,1) = Bz (i,jmax-1,1)  
           Ex (i,jmax  ,1) = Ex (i,jmax-1,1)  
           Ey (i,jmax  ,1) = Ey (i,jmax-1,1)  
           Ez (i,jmax  ,1) = Ez (i,jmax-1,1)  
           Vx (i,jmax  ,1) = Vx (i,jmax-1,1)  
           Vy (i,jmax  ,1) = Vy (i,jmax-1,1)  
           Vz (i,jmax  ,1) = Vz (i,jmax-1,1)  
           q  (i,jmax  ,1) = q  (i,jmax-1,1)   
           p  (i,jmax  ,1) = p  (i,jmax-1,1)   
           rho(i,jmax  ,1) = rho(i,jmax-1,1)   
           
           psi(i,jmax+1,1) = psi(i,jmax  ,1) 
           phi(i,jmax+1,1) = phi(i,jmax  ,1) 
           Bx (i,jmax+1,1) = Bx (i,jmax  ,1) 
           By (i,jmax+1,1) = By (i,jmax  ,1) 
           Bz (i,jmax+1,1) = Bz (i,jmax  ,1)  
           Ex (i,jmax+1,1) = Ex (i,jmax  ,1)  
           Ey (i,jmax+1,1) = Ey (i,jmax  ,1)  
           Ez (i,jmax+1,1) = Ez (i,jmax  ,1)  
           Vx (i,jmax+1,1) = Vx (i,jmax  ,1)  
           Vy (i,jmax+1,1) = Vy (i,jmax  ,1)  
           Vz (i,jmax+1,1) = Vz (i,jmax  ,1)  
           q  (i,jmax+1,1) = q  (i,jmax  ,1)   
           p  (i,jmax+1,1) = p  (i,jmax  ,1)   
           rho(i,jmax+1,1) = rho(i,jmax  ,1)   

           psi(i,jmax+2,1) = psi(i,jmax+1,1) 
           phi(i,jmax+2,1) = phi(i,jmax+1,1) 
           Bx (i,jmax+2,1) = Bx (i,jmax+1,1) 
           By (i,jmax+2,1) = By (i,jmax+1,1) 
           Bz (i,jmax+2,1) = Bz (i,jmax+1,1)  
           Ex (i,jmax+2,1) = Ex (i,jmax+1,1)  
           Ey (i,jmax+2,1) = Ey (i,jmax+1,1)  
           Ez (i,jmax+2,1) = Ez (i,jmax+1,1)  
           Vx (i,jmax+2,1) = Vx (i,jmax+1,1)  
           Vy (i,jmax+2,1) = Vy (i,jmax+1,1)  
           Vz (i,jmax+2,1) = Vz (i,jmax+1,1)  
           q  (i,jmax+2,1) = q  (i,jmax+1,1)   
           p  (i,jmax+2,1) = p  (i,jmax+1,1)   
           rho(i,jmax+2,1) = rho(i,jmax+1,1)   

           psi(i,jmax+3,1) = psi(i,jmax+2,1) 
           phi(i,jmax+3,1) = phi(i,jmax+2,1) 
           Bx (i,jmax+3,1) = Bx (i,jmax+2,1) 
           By (i,jmax+3,1) = By (i,jmax+2,1) 
           Bz (i,jmax+3,1) = Bz (i,jmax+2,1)  
           Ex (i,jmax+3,1) = Ex (i,jmax+2,1)  
           Ey (i,jmax+3,1) = Ey (i,jmax+2,1)  
           Ez (i,jmax+3,1) = Ez (i,jmax+2,1)  
           Vx (i,jmax+3,1) = Vx (i,jmax+2,1)  
           Vy (i,jmax+3,1) = Vy (i,jmax+2,1)  
           Vz (i,jmax+3,1) = Vz (i,jmax+2,1)  
           q  (i,jmax+3,1) = q  (i,jmax+2,1)   
           p  (i,jmax+3,1) = p  (i,jmax+2,1)   
           rho(i,jmax+3,1) = rho(i,jmax+2,1)   

           psi(i,jmax+4,1) = psi(i,jmax+3,1) 
           phi(i,jmax+4,1) = phi(i,jmax+3,1) 
           Bx (i,jmax+4,1) = Bx (i,jmax+3,1) 
           By (i,jmax+4,1) = By (i,jmax+3,1) 
           Bz (i,jmax+4,1) = Bz (i,jmax+3,1)  
           Ex (i,jmax+4,1) = Ex (i,jmax+3,1)  
           Ey (i,jmax+4,1) = Ey (i,jmax+3,1)  
           Ez (i,jmax+4,1) = Ez (i,jmax+3,1)  
           Vx (i,jmax+4,1) = Vx (i,jmax+3,1)  
           Vy (i,jmax+4,1) = Vy (i,jmax+3,1)  
           Vz (i,jmax+4,1) = Vz (i,jmax+3,1)  
           q  (i,jmax+4,1) = q  (i,jmax+3,1)   
           p  (i,jmax+4,1) = p  (i,jmax+3,1)   
           rho(i,jmax+4,1) = rho(i,jmax+3,1)

           psi(i,jmax+5,1) = psi(i,jmax+4,1) 
           phi(i,jmax+5,1) = phi(i,jmax+4,1) 
           Bx (i,jmax+5,1) = Bx (i,jmax+4,1) 
           By (i,jmax+5,1) = By (i,jmax+4,1) 
           Bz (i,jmax+5,1) = Bz (i,jmax+4,1)  
           Ex (i,jmax+5,1) = Ex (i,jmax+4,1)  
           Ey (i,jmax+5,1) = Ey (i,jmax+4,1)  
           Ez (i,jmax+5,1) = Ez (i,jmax+4,1)  
           Vx (i,jmax+5,1) = Vx (i,jmax+4,1)  
           Vy (i,jmax+5,1) = Vy (i,jmax+4,1)  
           Vz (i,jmax+5,1) = Vz (i,jmax+4,1)  
           q  (i,jmax+5,1) = q  (i,jmax+4,1)   
           p  (i,jmax+5,1) = p  (i,jmax+4,1)   
           rho(i,jmax+5,1) = rho(i,jmax+4,1)

           psi(i,jmax+6,1) = psi(i,jmax+5,1) 
           phi(i,jmax+6,1) = phi(i,jmax+5,1) 
           Bx (i,jmax+6,1) = Bx (i,jmax+5,1) 
           By (i,jmax+6,1) = By (i,jmax+5,1) 
           Bz (i,jmax+6,1) = Bz (i,jmax+5,1)  
           Ex (i,jmax+6,1) = Ex (i,jmax+5,1)  
           Ey (i,jmax+6,1) = Ey (i,jmax+5,1)  
           Ez (i,jmax+6,1) = Ez (i,jmax+5,1)  
           Vx (i,jmax+6,1) = Vx (i,jmax+5,1)  
           Vy (i,jmax+6,1) = Vy (i,jmax+5,1)  
           Vz (i,jmax+6,1) = Vz (i,jmax+5,1)  
           q  (i,jmax+6,1) = q  (i,jmax+5,1)   
           p  (i,jmax+6,1) = p  (i,jmax+5,1)   
           rho(i,jmax+6,1) = rho(i,jmax+5,1)   

 !        end do
     end do

!$OMP END DO


    else if (DIM == 2 .and. BOUND == 4) then

       ! Slab Jet Boundaries (Simetry outside inlet zone)

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO

    do j=-2,jmax+2
!       do k=1,kmax

           posy  = - 14.d0 + j * Dely


           if (posy .gt. -1.d0 .and. posy .le. 1.d0) then


           ! Left Boundary
           ! -----------------

           psi(0,j,1) = 0.d0
           phi(0,j,1) = 0.d0
           Bx (0,j,1) = 1.d0
           By (0,j,1) = 0.d0
           Bz (0,j,1) = 0.d0
           Vx (0,j,1) = 0.998752d0
           Vy (0,j,1) = 0.d0
           Vz (0,j,1) = 0.d0 
           Ex (0,j,1) = -(Vy(0,j,1)*Bz(0,j,1)-By(0,j,1)*Vz(0,j,1))
           Ey (0,j,1) = -(Vz(0,j,1)*Bx(0,j,1)-Bz(0,j,1)*Vx(0,j,1))
           Ez (0,j,1) = -(Vx(0,j,1)*By(0,j,1)-Bx(0,j,1)*Vy(0,j,1))
           q  (0,j,1) = 0.d0 
           rho(0,j,1) = 1.d-1
           p  (0,j,1) = 1.d-2


           psi(-1,j,1) = 0.d0
           phi(-1,j,1) = 0.d0
           Bx (-1,j,1) = 1.d0
           By (-1,j,1) = 0.d0
           Bz (-1,j,1) = 0.d0
           Vx (-1,j,1) = 0.998752d0
           Vy (-1,j,1) = 0.d0
           Vz (-1,j,1) = 0.d0 
           Ex (-1,j,1) = -(Vy(-1,j,1)*Bz(-1,j,1)-By(-1,j,1)*Vz(-1,j,1))
           Ey (-1,j,1) = -(Vz(-1,j,1)*Bx(-1,j,1)-Bz(-1,j,1)*Vx(-1,j,1))
           Ez (-1,j,1) = -(Vx(-1,j,1)*By(-1,j,1)-Bx(-1,j,1)*Vy(-1,j,1))
           q  (-1,j,1) = 0.d0 
           rho(-1,j,1) = 1.d-1
           p  (-1,j,1) = 1.d-2



           psi(-2,j,1) = 0.d0
           phi(-2,j,1) = 0.d0
           Bx (-2,j,1) = 1.d0
           By (-2,j,1) = 0.d0
           Bz (-2,j,1) = 0.d0
           Vx (-2,j,1) = 0.998752d0
           Vy (-2,j,1) = 0.d0
           Vz (-2,j,1) = 0.d0 
           Ex (-2,j,1) = -(Vy(-2,j,1)*Bz(-2,j,1)-By(-2,j,1)*Vz(-2,j,1))
           Ey (-2,j,1) = -(Vz(-2,j,1)*Bx(-2,j,1)-Bz(-2,j,1)*Vx(-2,j,1))
           Ez (-2,j,1) = -(Vx(-2,j,1)*By(-2,j,1)-Bx(-2,j,1)*Vy(-2,j,1))
           q  (-2,j,1) = 0.d0 
           rho(-2,j,1) = 1.d-1
           p  (-2,j,1) = 1.d-2

           else 

           psi(0,j,1) =  psi(1,j,1) 
           phi(0,j,1) =  phi(1,j,1) 
           Bx(0,j,1)  =  Bx(1,j,1)  
           By(0,j,1)  = -By(1,j,1)  
           Bz(0,j,1)  = -Bz(1,j,1)  
           Vx(0,j,1)  =  Vx(1,j,1)  
           Vy(0,j,1)  = -Vy(1,j,1)  
           Vz(0,j,1)  = -Vz(1,j,1)
           Ex(0,j,1)  =  Ex(1,j,1)  
           Ey(0,j,1)  = -Ey(1,j,1)  
           Ez(0,j,1)  = -Ez(1,j,1)    
           q(0,j,1)   =  q(1,j,1)   
           p(0,j,1)   =  p(1,j,1)
           rho(0,j,1) =  rho(1,j,1)

           psi(-1,j,1) = psi(0,j,1) 
           phi(-1,j,1) = phi(0,j,1) 
           Bx(-1,j,1)  = Bx(0,j,1)  
           By(-1,j,1)  = By(0,j,1)  
           Bz(-1,j,1)  = Bz(0,j,1)  
           Vx(-1,j,1)  = Vx(0,j,1)  
           Vy(-1,j,1)  = Vy(0,j,1)  
           Vz(-1,j,1)  = Vz(0,j,1)
           Ex(-1,j,1)  = Ex(0,j,1)  
           Ey(-1,j,1)  = Ey(0,j,1)  
           Ez(-1,j,1)  = Ez(0,j,1)    
           q(-1,j,1)   = q(0,j,1)   
           p(-1,j,1)   = p(0,j,1)
           rho(-1,j,1) = rho(0,j,1)

           psi(-2,j,1) = psi(-1,j,1) 
           phi(-2,j,1) = phi(-1,j,1) 
           Bx(-2,j,1)  = Bx(-1,j,1)  
           By(-2,j,1)  = By(-1,j,1)  
           Bz(-2,j,1)  = Bz(-1,j,1)  
           Vx(-2,j,1)  = Vx(-1,j,1)  
           Vy(-2,j,1)  = Vy(-1,j,1)  
           Vz(-2,j,1)  = Vz(-1,j,1)
           Ex(-2,j,1)  = Ex(-1,j,1)  
           Ey(-2,j,1)  = Ey(-1,j,1)  
           Ez(-2,j,1)  = Ez(-1,j,1)    
           q(-2,j,1)   = q(-1,j,1)   
           p(-2,j,1)   = p(-1,j,1)
           rho(-2,j,1) = rho(-1,j,1)

           end if


           ! Right  Boundary
           ! -----------------

           psi(imax,j,1) = psi(imax-1,j,1) 
           phi(imax,j,1) = phi(imax-1,j,1) 
           Bx(imax,j,1)  = Bx(imax-1,j,1) 
           By(imax,j,1)  = By(imax-1,j,1) 
           Bz(imax,j,1)  = Bz(imax-1,j,1) 
           Ex(imax,j,1)  = Ex(imax-1,j,1) 
           Ey(imax,j,1)  = Ey(imax-1,j,1) 
           Ez(imax,j,1)  = Ez(imax-1,j,1) 
           Vx(imax,j,1)  = Vx(imax-1,j,1) 
           Vy(imax,j,1)  = Vy(imax-1,j,1) 
           Vz(imax,j,1)  = Vz(imax-1,j,1) 
           q(imax,j,1)   = q(imax-1,j,1) 
           p(imax,j,1)   = p(imax-1,j,1) 
           rho(imax,j,1) = rho(imax-1,j,1) 

           psi(imax+1,j,1) = psi(imax,j,1) 
           phi(imax+1,j,1) = phi(imax,j,1) 
           Bx(imax+1,j,1)  = Bx(imax,j,1) 
           By(imax+1,j,1)  = By(imax,j,1) 
           Bz(imax+1,j,1)  = Bz(imax,j,1) 
           Ex(imax+1,j,1)  = Ex(imax,j,1) 
           Ey(imax+1,j,1)  = Ey(imax,j,1) 
           Ez(imax+1,j,1)  = Ez(imax,j,1) 
           Vx(imax+1,j,1)  = Vx(imax,j,1) 
           Vy(imax+1,j,1)  = Vy(imax,j,1) 
           Vz(imax+1,j,1)  = Vz(imax,j,1) 
           q(imax+1,j,1)   = q(imax,j,1) 
           p(imax+1,j,1)   = p(imax,j,1) 
           rho(imax+1,j,1) = rho(imax,j,1) 

           psi(imax+2,j,1) = psi(imax+1,j,1)
           phi(imax+2,j,1) = phi(imax+1,j,1)
           Bx(imax+2,j,1)  = Bx(imax+1,j,1) 
           By(imax+2,j,1)  = By(imax+1,j,1) 
           Bz(imax+2,j,1)  = Bz(imax+1,j,1) 
           Ex(imax+2,j,1)  = Ex(imax+1,j,1) 
           Ey(imax+2,j,1)  = Ey(imax+1,j,1) 
           Ez(imax+2,j,1)  = Ez(imax+1,j,1) 
           Vx(imax+2,j,1)  = Vx(imax+1,j,1) 
           Vy(imax+2,j,1)  = Vy(imax+1,j,1) 
           Vz(imax+2,j,1)  = Vz(imax+1,j,1) 
           q(imax+2,j,1)   = q(imax+1,j,1) 
           p(imax+2,j,1)   = p(imax+1,j,1) 
           rho(imax+2,j,1) = rho(imax+1,j,1) 


!        end do
     end do

!$OMP END DO 

  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO 

    do i=-2,imax+2
!       do k=1,kmax

 
           ! Left Boundary
           ! -----------------

           psi(i,0,1) = psi(i,1,1)
           phi(i,0,1) = phi(i,1,1) 
           Bx(i,0,1)  = Bx(i,1,1)
           By(i,0,1)  = By(i,1,1)
           Bz(i,0,1)  = Bz(i,1,1)
           Ex(i,0,1)  = Ex(i,1,1)  
           Ey(i,0,1)  = Ey(i,1,1)  
           Ez(i,0,1)  = Ez(i,1,1)  
           Vx(i,0,1)  = Vx(i,1,1)  
           Vy(i,0,1)  = Vy(i,1,1)  
           Vz(i,0,1)  = Vz(i,1,1)  
           q(i,0,1)   = q(i,1,1)   
           p(i,0,1)   = p(i,1,1)
           rho(i,0,1) = rho(i,1,1)

           psi(i,-1,1) = psi(i,0,1)
           phi(i,-1,1) = phi(i,0,1) 
           Bx(i,-1,1)  = Bx(i,0,1)
           By(i,-1,1)  = By(i,0,1)
           Bz(i,-1,1)  = Bz(i,0,1)
           Ex(i,-1,1)  = Ex(i,0,1)  
           Ey(i,-1,1)  = Ey(i,0,1)  
           Ez(i,-1,1)  = Ez(i,0,1)  
           Vx(i,-1,1)  = Vx(i,0,1)  
           Vy(i,-1,1)  = Vy(i,0,1)  
           Vz(i,-1,1)  = Vz(i,0,1)  
           q(i,-1,1)   = q(i,0,1)   
           p(i,-1,1)   = p(i,0,1)
           rho(i,-1,1) = rho(i,0,1)

           psi(i,-2,1) = psi(i,-1,1)
           phi(i,-2,1) = phi(i,-1,1) 
           Bx(i,-2,1)  = Bx(i,-1,1)
           By(i,-2,1)  = By(i,-1,1)
           Bz(i,-2,1)  = Bz(i,-1,1)
           Ex(i,-2,1)  = Ex(i,-1,1)  
           Ey(i,-2,1)  = Ey(i,-1,1)  
           Ez(i,-2,1)  = Ez(i,-1,1)  
           Vx(i,-2,1)  = Vx(i,-1,1)  
           Vy(i,-2,1)  = Vy(i,-1,1)  
           Vz(i,-2,1)  = Vz(i,-1,1)  
           q(i,-2,1)   = q(i,-1,1)   
           p(i,-2,1)   = p(i,-1,1)
           rho(i,-2,1) = rho(i,-1,1)

           ! Right  Boundary
           ! -----------------

           psi(i,jmax,1) = psi(i,jmax-1,1) 
           phi(i,jmax,1) = phi(i,jmax-1,1) 
           Bx(i,jmax,1)  = Bx(i,jmax-1,1) 
           By(i,jmax,1)  = By(i,jmax-1,1) 
           Bz(i,jmax,1)  = Bz(i,jmax-1,1)  
           Ex(i,jmax,1)  = Ex(i,jmax-1,1)  
           Ey(i,jmax,1)  = Ey(i,jmax-1,1)  
           Ez(i,jmax,1)  = Ez(i,jmax-1,1)  
           Vx(i,jmax,1)  = Vx(i,jmax-1,1)  
           Vy(i,jmax,1)  = Vy(i,jmax-1,1)  
           Vz(i,jmax,1)  = Vz(i,jmax-1,1)  
           q(i,jmax,1)   = q(i,jmax-1,1)   
           p(i,jmax,1)   = p(i,jmax-1,1)   
           rho(i,jmax,1) = rho(i,jmax-1,1)   

           psi(i,jmax+1,1) = psi(i,jmax,1) 
           phi(i,jmax+1,1) = phi(i,jmax,1) 
           Bx(i,jmax+1,1)  = Bx(i,jmax,1) 
           By(i,jmax+1,1)  = By(i,jmax,1) 
           Bz(i,jmax+1,1)  = Bz(i,jmax,1)  
           Ex(i,jmax+1,1)  = Ex(i,jmax,1)  
           Ey(i,jmax+1,1)  = Ey(i,jmax,1)  
           Ez(i,jmax+1,1)  = Ez(i,jmax,1)  
           Vx(i,jmax+1,1)  = Vx(i,jmax,1)  
           Vy(i,jmax+1,1)  = Vy(i,jmax,1)  
           Vz(i,jmax+1,1)  = Vz(i,jmax,1)  
           q(i,jmax+1,1)   = q(i,jmax,1)   
           p(i,jmax+1,1)   = p(i,jmax,1)   
           rho(i,jmax+1,1) = rho(i,jmax,1)   


           psi(i,jmax+2,1) = psi(i,jmax+1,1) 
           phi(i,jmax+2,1) = phi(i,jmax+1,1) 
           Bx(i,jmax+2,1)  = Bx(i,jmax+1,1) 
           By(i,jmax+2,1)  = By(i,jmax+1,1) 
           Bz(i,jmax+2,1)  = Bz(i,jmax+1,1)  
           Ex(i,jmax+2,1)  = Ex(i,jmax+1,1)  
           Ey(i,jmax+2,1)  = Ey(i,jmax+1,1)  
           Ez(i,jmax+2,1)  = Ez(i,jmax+1,1)  
           Vx(i,jmax+2,1)  = Vx(i,jmax+1,1)  
           Vy(i,jmax+2,1)  = Vy(i,jmax+1,1)  
           Vz(i,jmax+2,1)  = Vz(i,jmax+1,1)  
           q(i,jmax+2,1)   = q(i,jmax+1,1)   
           p(i,jmax+2,1)   = p(i,jmax+1,1)   
           rho(i,jmax+2,1) = rho(i,jmax+1,1)   

 
!        end do
     end do

!$OMP END DO 


    else if (DIM == 2 .and. BOUND == 5) then

       ! Cloud Shock Interaction

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO  

    do j=-2,jmax+2
!       do k=1,kmax
 

           ! Left Boundary
           ! -----------------

           psi(0,j,1) = psi(1,j,1) 
           phi(0,j,1) = phi(1,j,1) 
           Bx(0,j,1)  = Bx(1,j,1)  
           By(0,j,1)  = By(1,j,1)  
           Bz(0,j,1)  = Bz(1,j,1)  
           Ex(0,j,1)  = Ex(1,j,1)  
           Ey(0,j,1)  = Ey(1,j,1)  
           Ez(0,j,1)  = Ez(1,j,1)  
           Vx(0,j,1)  = Vx(1,j,1)  
           Vy(0,j,1)  = Vy(1,j,1)  
           Vz(0,j,1)  = Vz(1,j,1)  
           q(0,j,1)   = q(1,j,1)   
           p(0,j,1)   = p(1,j,1)   
           rho(0,j,1) = rho(1,j,1)  

           psi(-1,j,1) = psi(0,j,1) 
           phi(-1,j,1) = phi(0,j,1) 
           Bx(-1,j,1)  = Bx(0,j,1)  
           By(-1,j,1)  = By(0,j,1)  
           Bz(-1,j,1)  = Bz(0,j,1)  
           Ex(-1,j,1)  = Ex(0,j,1)  
           Ey(-1,j,1)  = Ey(0,j,1)  
           Ez(-1,j,1)  = Ez(0,j,1)  
           Vx(-1,j,1)  = Vx(0,j,1)  
           Vy(-1,j,1)  = Vy(0,j,1)  
           Vz(-1,j,1)  = Vz(0,j,1)  
           q(-1,j,1)   = q(0,j,1)   
           p(-1,j,1)   = p(0,j,1)   
           rho(-1,j,1) = rho(0,j,1)  


           psi(-2,j,1) = psi(-1,j,1) 
           phi(-2,j,1) = phi(-1,j,1) 
           Bx(-2,j,1)  = Bx(-1,j,1)  
           By(-2,j,1)  = By(-1,j,1)  
           Bz(-2,j,1)  = Bz(-1,j,1)  
           Ex(-2,j,1)  = Ex(-1,j,1)  
           Ey(-2,j,1)  = Ey(-1,j,1)  
           Ez(-2,j,1)  = Ez(-1,j,1)  
           Vx(-2,j,1)  = Vx(-1,j,1)  
           Vy(-2,j,1)  = Vy(-1,j,1)  
           Vz(-2,j,1)  = Vz(-1,j,1)  
           q(-2,j,1)   = q(-1,j,1)   
           p(-2,j,1)   = p(-1,j,1)
           rho(-2,j,1) = rho(-1,j,1)

           ! Right  Boundary
           ! -----------------

           psi(imax,j,1) = 0.d0
           phi(imax,j,1) = 0.d0
           Bx(imax,j,1)  = 0.d0
           By(imax,j,1)  = 0.564190d0
           Bz(imax,j,1)  = 0.564190d0
           Vx(imax,j,1)  = -0.995d0 
           Vy(imax,j,1)  = 0.d0
           Vz(imax,j,1)  = 0.d0
           Ex(imax,j,1)  = -(Vy(imax,j,1)*Bz(imax,j,1)-By(imax,j,1)*Vz(imax,j,1))
           Ey(imax,j,1)  = -(Vz(imax,j,1)*Bx(imax,j,1)-Bz(imax,j,1)*Vx(imax,j,1)) 
           Ez(imax,j,1)  = -(Vx(imax,j,1)*By(imax,j,1)-Bx(imax,j,1)*Vy(imax,j,1))  
           q(imax,j,1)   = 0.d0
           p(imax,j,1)   = 1.d0
           rho(imax,j,1) = 1.d0

           psi(imax+1,j,1) = 0.d0 !psi(imax,j,1) 
           phi(imax+1,j,1) = 0.d0 !phi(imax,j,1) 
           Bx(imax+1,j,1)  = 0.d0 !Bx(imax,j,1) 
           By(imax+1,j,1)  = 0.564190d0 !By(imax,j,1) 
           Bz(imax+1,j,1)  = 0.564190d0 !Bz(imax,j,1)
           Vx(imax+1,j,1)  = -0.995d0!Vx(imax,j,1) 
           Vy(imax+1,j,1)  = 0.d0 !Vy(imax,j,1) 
           Vz(imax+1,j,1)  = 0.d0 !Vz(imax,j,1) 
           Ex(imax+1,j,1)  = -(Vy(imax+1,j,1)*Bz(imax+1,j,1)-By(imax+1,j,1)*Vz(imax+1,j,1))
           Ey(imax+1,j,1)  = -(Vz(imax+1,j,1)*Bx(imax+1,j,1)-Bz(imax+1,j,1)*Vx(imax+1,j,1)) 
           Ez(imax+1,j,1)  = -(Vx(imax+1,j,1)*By(imax+1,j,1)-Bx(imax+1,j,1)*Vy(imax+1,j,1))  
           q(imax+1,j,1)   = 0.d0 !q(imax,j,1) 
           p(imax+1,j,1)   = 1.d0 !p(imax,j,1) 
           rho(imax+1,j,1) = 1.d0 !rho(imax,j,1) 

           psi(imax+2,j,1) = 0.d0 !psi(imax+1,j,1)
           phi(imax+2,j,1) = 0.d0 !phi(imax+1,j,1)
           Bx(imax+2,j,1)  = 0.d0 !Bx(imax+1,j,1) 
           By(imax+2,j,1)  = 0.564190d0 !By(imax+1,j,1) 
           Bz(imax+2,j,1)  = 0.564190d0 !Bz(imax+1,j,1) 
           Vx(imax+2,j,1)  = -0.995d0!Vx(imax+1,j,1) 
           Vy(imax+2,j,1)  = 0.d0 !Vy(imax+1,j,1) 
           Vz(imax+2,j,1)  = 0.d0 !Vz(imax+1,j,1)
           Ex(imax+2,j,1)  = -(Vy(imax+2,j,1)*Bz(imax+2,j,1)-By(imax+2,j,1)*Vz(imax+2,j,1))
           Ey(imax+2,j,1)  = -(Vz(imax+2,j,1)*Bx(imax+2,j,1)-Bz(imax+2,j,1)*Vx(imax+2,j,1)) 
           Ez(imax+2,j,1)  = -(Vx(imax+2,j,1)*By(imax+2,j,1)-Bx(imax+2,j,1)*Vy(imax+2,j,1))  
           q(imax+2,j,1)   = 0.d0 !q(imax+1,j,1) 
           p(imax+2,j,1)   = 1.d0 !p(imax+1,j,1) 
           rho(imax+2,j,1) = 1.d0 !rho(imax+1,j,1) 

 
!        end do
     end do

!$OMP END DO 

  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO  

    do i=-2,imax+2
!       do k=1,kmax

 


           ! Left Boundary
           ! -----------------

           psi(i,0,1) = psi(i,1,1)
           phi(i,0,1) = phi(i,1,1) 
           Bx(i,0,1)  = Bx(i,1,1)
           By(i,0,1)  = By(i,1,1)
           Bz(i,0,1)  = Bz(i,1,1)
           Ex(i,0,1)  = Ex(i,1,1)  
           Ey(i,0,1)  = Ey(i,1,1)  
           Ez(i,0,1)  = Ez(i,1,1)  
           Vx(i,0,1)  = Vx(i,1,1)  
           Vy(i,0,1)  = Vy(i,1,1)  
           Vz(i,0,1)  = Vz(i,1,1)  
           q(i,0,1)   = q(i,1,1)   
           p(i,0,1)   = p(i,1,1)
           rho(i,0,1) = rho(i,1,1)

           psi(i,-1,1) = psi(i,0,1)
           phi(i,-1,1) = phi(i,0,1) 
           Bx(i,-1,1)  = Bx(i,0,1)
           By(i,-1,1)  = By(i,0,1)
           Bz(i,-1,1)  = Bz(i,0,1)
           Ex(i,-1,1)  = Ex(i,0,1)  
           Ey(i,-1,1)  = Ey(i,0,1)  
           Ez(i,-1,1)  = Ez(i,0,1)  
           Vx(i,-1,1)  = Vx(i,0,1)  
           Vy(i,-1,1)  = Vy(i,0,1)  
           Vz(i,-1,1)  = Vz(i,0,1)  
           q(i,-1,1)   = q(i,0,1)   
           p(i,-1,1)   = p(i,0,1)
           rho(i,-1,1) = rho(i,0,1)

           psi(i,-2,1) = psi(i,-1,1)
           phi(i,-2,1) = phi(i,-1,1) 
           Bx(i,-2,1)  = Bx(i,-1,1)
           By(i,-2,1)  = By(i,-1,1)
           Bz(i,-2,1)  = Bz(i,-1,1)
           Ex(i,-2,1)  = Ex(i,-1,1)  
           Ey(i,-2,1)  = Ey(i,-1,1)  
           Ez(i,-2,1)  = Ez(i,-1,1)  
           Vx(i,-2,1)  = Vx(i,-1,1)  
           Vy(i,-2,1)  = Vy(i,-1,1)  
           Vz(i,-2,1)  = Vz(i,-1,1)  
           q(i,-2,1)   = q(i,-1,1)   
           p(i,-2,1)   = p(i,-1,1)
           rho(i,-2,1) = rho(i,-1,1)

           ! Right  Boundary
           ! -----------------

           psi(i,jmax,1) = psi(i,jmax-1,1) 
           phi(i,jmax,1) = phi(i,jmax-1,1) 
           Bx(i,jmax,1)  = Bx(i,jmax-1,1) 
           By(i,jmax,1)  = By(i,jmax-1,1) 
           Bz(i,jmax,1)  = Bz(i,jmax-1,1)  
           Ex(i,jmax,1)  = Ex(i,jmax-1,1)  
           Ey(i,jmax,1)  = Ey(i,jmax-1,1)  
           Ez(i,jmax,1)  = Ez(i,jmax-1,1)  
           Vx(i,jmax,1)  = Vx(i,jmax-1,1)  
           Vy(i,jmax,1)  = Vy(i,jmax-1,1)  
           Vz(i,jmax,1)  = Vz(i,jmax-1,1)  
           q(i,jmax,1)   = q(i,jmax-1,1)   
           p(i,jmax,1)   = p(i,jmax-1,1)   
           rho(i,jmax,1) = rho(i,jmax-1,1)   

           psi(i,jmax+1,1) = psi(i,jmax,1) 
           phi(i,jmax+1,1) = phi(i,jmax,1) 
           Bx(i,jmax+1,1)  = Bx(i,jmax,1) 
           By(i,jmax+1,1)  = By(i,jmax,1) 
           Bz(i,jmax+1,1)  = Bz(i,jmax,1)  
           Ex(i,jmax+1,1)  = Ex(i,jmax,1)  
           Ey(i,jmax+1,1)  = Ey(i,jmax,1)  
           Ez(i,jmax+1,1)  = Ez(i,jmax,1)  
           Vx(i,jmax+1,1)  = Vx(i,jmax,1)  
           Vy(i,jmax+1,1)  = Vy(i,jmax,1)  
           Vz(i,jmax+1,1)  = Vz(i,jmax,1)  
           q(i,jmax+1,1)   = q(i,jmax,1)   
           p(i,jmax+1,1)   = p(i,jmax,1)   
           rho(i,jmax+1,1) = rho(i,jmax,1)   


           psi(i,jmax+2,1) = psi(i,jmax+1,1) 
           phi(i,jmax+2,1) = phi(i,jmax+1,1) 
           Bx(i,jmax+2,1)  = Bx(i,jmax+1,1) 
           By(i,jmax+2,1)  = By(i,jmax+1,1) 
           Bz(i,jmax+2,1)  = Bz(i,jmax+1,1)  
           Ex(i,jmax+2,1)  = Ex(i,jmax+1,1)  
           Ey(i,jmax+2,1)  = Ey(i,jmax+1,1)  
           Ez(i,jmax+2,1)  = Ez(i,jmax+1,1)  
           Vx(i,jmax+2,1)  = Vx(i,jmax+1,1)  
           Vy(i,jmax+2,1)  = Vy(i,jmax+1,1)  
           Vz(i,jmax+2,1)  = Vz(i,jmax+1,1)  
           q(i,jmax+2,1)   = q(i,jmax+1,1)   
           p(i,jmax+2,1)   = p(i,jmax+1,1)   
           rho(i,jmax+2,1) = rho(i,jmax+1,1)   

 
!        end do
     end do

!$OMP END DO 

    else if (DIM == 2 .and. BOUND == 6) then 


!Tomek

       !stencil = 6 defined in scalar.f95

       
       ! Tearing-Mode

  ! Periodic Boundary y-direction


  !-----------------
  ! Y Direction
  !-----------------

       ! Left Boundary
       ! -----------------
       
!$OMP DO ORDERED
       
       do i=-stencil,imax+stencil
          do j = 1, stencil

!$OMP ORDERED

             psi(i,0-j,1) = psi(i,jmax-j,1)
             phi(i,0-j,1) = phi(i,jmax-j,1)
             Bx (i,0-j,1) = Bx (i,jmax-j,1)
             By (i,0-j,1) = By (i,jmax-j,1)
             Bz (i,0-j,1) = Bz (i,jmax-j,1)
             Ex (i,0-j,1) = Ex (i,jmax-j,1)
             Ey (i,0-j,1) = Ey (i,jmax-j,1)
             Ez (i,0-j,1) = Ez (i,jmax-j,1)
             Vx (i,0-j,1) = Vx (i,jmax-j,1)
             Vy (i,0-j,1) = Vy (i,jmax-j,1)
             Vz (i,0-j,1) = Vz (i,jmax-j,1)
             q  (i,0-j,1) = q  (i,jmax-j,1)
             p  (i,0-j,1) = p  (i,jmax-j,1)
             rho(i,0-j,1) = rho(i,jmax-j,1)

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
             
             psi(i,jmax+j,1) = psi(i,0+j,1)
             phi(i,jmax+j,1) = phi(i,0+j,1)
             Bx (i,jmax+j,1) = Bx (i,0+j,1)
             By (i,jmax+j,1) = By (i,0+j,1)
             Bz (i,jmax+j,1) = Bz (i,0+j,1)
             Ex (i,jmax+j,1) = Ex (i,0+j,1)
             Ey (i,jmax+j,1) = Ey (i,0+j,1)
             Ez (i,jmax+j,1) = Ez (i,0+j,1)
             Vx (i,jmax+j,1) = Vx (i,0+j,1)
             Vy (i,jmax+j,1) = Vy (i,0+j,1)
             Vz (i,jmax+j,1) = Vz (i,0+j,1)
             q  (i,jmax+j,1) = q  (i,0+j,1)
             p  (i,jmax+j,1) = p  (i,0+j,1)
             rho(i,jmax+j,1) = rho(i,0+j,1)

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

!$OMP DO ORDERED 
       
       do i = 1,stencil
          do j=-stencil,jmax+stencil

!$OMP ORDERED             

           phi(0-i,j,1) = phi(0  ,j,1) 
           psi(0-i,j,1) = psi(0  ,j,1) 

           Bx (0-i,j,1) = Bx (0  ,j,1) 
           By (0-i,j,1) = By (0  ,j,1) 
           Bz (0-i,j,1) = Bz (0  ,j,1) 

           Ex (0-i,j,1) = Ex (0  ,j,1) 
           Ey (0-i,j,1) = Ey (0  ,j,1) 
           Ez (0-i,j,1) = Ez (0  ,j,1)

           Vx (0-i,j,1) = Vx (0  ,j,1) 
           Vy (0-i,j,1) = Vy (0  ,j,1) 
           Vz (0-i,j,1) = Vz (0  ,j,1) 

           q  (0-i,j,1) = q  (0  ,j,1)  
           p  (0-i,j,1) = p  (0  ,j,1) 
           rho(0-i,j,1) = rho(0  ,j,1)

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
              
           phi(imax+i,j,1) = phi(imax  ,j,1) 
           psi(imax+i,j,1) = psi(imax  ,j,1) 

           Bx (imax+i,j,1) = Bx (imax  ,j,1)               
           By (imax+i,j,1) = By (imax  ,j,1) 
           Bz (imax+i,j,1) = Bz (imax  ,j,1) 

           Ex (imax+i,j,1) = Ex (imax  ,j,1)
           Ey (imax+i,j,1) = Ey (imax  ,j,1) 
           Ez (imax+i,j,1) = Ez (imax  ,j,1)
           
           Vx (imax+i,j,1) = Vx (imax  ,j,1) 
           Vy (imax+i,j,1) = Vy (imax  ,j,1) 
           Vz (imax+i,j,1) = Vz (imax  ,j,1)
           
           q  (imax+i,j,1) = q  (imax  ,j,1)  
           p  (imax+i,j,1) = p  (imax  ,j,1) 
           rho(imax+i,j,1) = rho(imax  ,j,1)

!$OMP END ORDERED              

        end do
     end do

!$OMP END DO 

       
!Tomek     
     

!***************************************************************************************


    else if (DIM == 2 .and. BOUND ==7) then

       ! Reconnection boundary Zenitani et al.

  !-----------------
  ! X Direction
  !-----------------

 !$OMP DO 

    do j=-2,jmax+2
!       do k=1,kmax

 


           ! Left Boundary
           ! -----------------

           psi(0,j,1) = - psi(1,j,1) 
           phi(0,j,1) = - phi(1,j,1) 
           Bx(0,j,1)  = - Bx(1,j,1)  
           By(0,j,1)  = - By(1,j,1)  
           Bz(0,j,1)  = - Bz(1,j,1)  
           Ex(0,j,1)  = - Ex(1,j,1)  
           Ey(0,j,1)  = - Ey(1,j,1)  
           Ez(0,j,1)  = - Ez(1,j,1)  
           Vx(0,j,1)  = - Vx(1,j,1)  
           Vy(0,j,1)  = - Vy(1,j,1)  
           Vz(0,j,1)  = - Vz(1,j,1)  
           q(0,j,1)   = - q(1,j,1)   
           p(0,j,1)   = - p(1,j,1)   
           rho(0,j,1) = - rho(1,j,1)  

           psi(-1,j,1) = psi(0,j,1) 
           phi(-1,j,1) = phi(0,j,1) 
           Bx(-1,j,1)  = Bx(0,j,1)  
           By(-1,j,1)  = By(0,j,1)  
           Bz(-1,j,1)  = Bz(0,j,1)  
           Ex(-1,j,1)  = Ex(0,j,1)  
           Ey(-1,j,1)  = Ey(0,j,1)  
           Ez(-1,j,1)  = Ez(0,j,1)  
           Vx(-1,j,1)  = Vx(0,j,1)  
           Vy(-1,j,1)  = Vy(0,j,1)  
           Vz(-1,j,1)  = Vz(0,j,1)  
           q(-1,j,1)   = q(0,j,1)   
           p(-1,j,1)   = p(0,j,1)   
           rho(-1,j,1) = rho(0,j,1)  


           psi(-2,j,1) = psi(-1,j,1) 
           phi(-2,j,1) = phi(-1,j,1) 
           Bx(-2,j,1)  = Bx(-1,j,1)  
           By(-2,j,1)  = By(-1,j,1)  
           Bz(-2,j,1)  = Bz(-1,j,1)  
           Ex(-2,j,1)  = Ex(-1,j,1)  
           Ey(-2,j,1)  = Ey(-1,j,1)  
           Ez(-2,j,1)  = Ez(-1,j,1)  
           Vx(-2,j,1)  = Vx(-1,j,1)  
           Vy(-2,j,1)  = Vy(-1,j,1)  
           Vz(-2,j,1)  = Vz(-1,j,1)  
           q(-2,j,1)   = q(-1,j,1)   
           p(-2,j,1)   = p(-1,j,1)
           rho(-2,j,1) = rho(-1,j,1)


           ! Right  Boundary
           ! -----------------
           psi(imax,j,1) = psi(imax-1,j,1) 
           phi(imax,j,1) = phi(imax-1,j,1) 
           Bx(imax,j,1)  = Bx(imax-1,j,1) 
           By(imax,j,1)  = By(imax-1,j,1) 
           Bz(imax,j,1)  = Bz(imax-1,j,1) 
           Ex(imax,j,1)  = Ex(imax-1,j,1) 
           Ey(imax,j,1)  = Ey(imax-1,j,1) 
           Ez(imax,j,1)  = Ez(imax-1,j,1) 
           Vx(imax,j,1)  = Vx(imax-1,j,1) 
           Vy(imax,j,1)  = Vy(imax-1,j,1) 
           Vz(imax,j,1)  = Vz(imax-1,j,1) 
           q(imax,j,1)   = q(imax-1,j,1) 
           p(imax,j,1)   = p(imax-1,j,1) 
           rho(imax,j,1) = rho(imax-1,j,1) 

           psi(imax+1,j,1) = psi(imax,j,1) 
           phi(imax+1,j,1) = phi(imax,j,1) 
           Bx(imax+1,j,1)  = Bx(imax,j,1) 
           By(imax+1,j,1)  = By(imax,j,1) 
           Bz(imax+1,j,1)  = Bz(imax,j,1) 
           Ex(imax+1,j,1)  = Ex(imax,j,1) 
           Ey(imax+1,j,1)  = Ey(imax,j,1) 
           Ez(imax+1,j,1)  = Ez(imax,j,1) 
           Vx(imax+1,j,1)  = Vx(imax,j,1) 
           Vy(imax+1,j,1)  = Vy(imax,j,1) 
           Vz(imax+1,j,1)  = Vz(imax,j,1) 
           q(imax+1,j,1)   = q(imax,j,1) 
           p(imax+1,j,1)   = p(imax,j,1) 
           rho(imax+1,j,1) = rho(imax,j,1) 

           psi(imax+2,j,1) = psi(imax,j,1)
           phi(imax+2,j,1) = phi(imax,j,1)
           Bx(imax+2,j,1)  = Bx(imax,j,1) 
           By(imax+2,j,1)  = By(imax,j,1) 
           Bz(imax+2,j,1)  = Bz(imax,j,1) 
           Ex(imax+2,j,1)  = Ex(imax,j,1) 
           Ey(imax+2,j,1)  = Ey(imax,j,1) 
           Ez(imax+2,j,1)  = Ez(imax,j,1) 
           Vx(imax+2,j,1)  = Vx(imax,j,1) 
           Vy(imax+2,j,1)  = Vy(imax,j,1) 
           Vz(imax+2,j,1)  = Vz(imax,j,1) 
           q(imax+2,j,1)   = q(imax,j,1) 
           p(imax+2,j,1)   = p(imax,j,1) 
           rho(imax+2,j,1) = rho(imax,j,1) 

 
!        end do
     end do

!$OMP END DO 

  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO  

    do i=-2,imax+2
!       do k=1,kmax


 

           posx   =  i * Delx
           posy   = - 0.5d0 * Ly + jmax * Dely

           ! Left Boundary
           ! -----------------

           psi(i,0,1) = psi(i,1,1)
           phi(i,0,1) = phi(i,1,1) 
           Bx(i,0,1)  = Bx(i,1,1)
           By(i,0,1)  = By(i,1,1) 
           Bz(i,0,1)  = Bz(i,1,1)
           Vx(i,0,1)  = Vx(i,1,1)  
           Vy(i,0,1)  = Vy(i,1,1)  
           Vz(i,0,1)  = Vz(i,1,1)
           Ex(i,0,1)  = Ex(i,1,1)
           Ey(i,0,1)  = Ey(i,1,1)
           Ez(i,0,1)  = Ez(i,1,1)
           q(i,0,1)   = q(i,1,1)   
           p(i,0,1)   = p(i,1,1)
           rho(i,0,1) = rho(i,1,1)


           psi(i,-1,1) = psi(i,0,1)
           phi(i,-1,1) = phi(i,0,1) 
           Bx(i,-1,1)  = Bx(i,0,1)
           By(i,-1,1)  = By(i,0,1) 
           Bz(i,-1,1)  = Bz(i,0,1)
           Vx(i,-1,1)  = Vx(i,0,1)  
           Vy(i,-1,1)  = Vy(i,0,1)  
           Vz(i,-1,1)  = Vz(i,0,1)
           Ex(i,-1,1)  = Ex(i,0,1)
           Ey(i,-1,1)  = Ey(i,0,1)
           Ez(i,-1,1)  = Ez(i,0,1)
           q(i,-1,1)   = q(i,0,1)   
           p(i,-1,1)   = p(i,0,1)
           rho(i,-1,1) = rho(i,0,1)

           psi(i,-2,1) = psi(i,-1,1)
           phi(i,-2,1) = phi(i,-1,1) 
           Bx(i,-2,1)  = Bx(i,-1,1)
           By(i,-2,1)  = By(i,-1,1) 
           Bz(i,-2,1)  = Bz(i,-1,1)
           Vx(i,-2,1)  = Vx(i,-1,1)  
           Vy(i,-2,1)  = Vy(i,-1,1)  
           Vz(i,-2,1)  = Vz(i,-1,1) 
           Ex(i,-2,1)  = Ex(i,-1,1)
           Ey(i,-2,1)  = Ey(i,-1,1)
           Ez(i,-2,1)  = Ez(i,-1,1)
           q(i,-2,1)   = q(i,-1,1)   
           p(i,-2,1)   = p(i,-1,1)
           rho(i,-2,1) = rho(i,-1,1)

           ! Right  Boundary
           ! -----------------

           psi(i,jmax  ,1) = psi(i,jmax-1,1) 
           phi(i,jmax  ,1) = phi(i,jmax-1,1) 
           Bx (i,jmax  ,1) = Bx (i,jmax-1,1)
           By (i,jmax  ,1) = By (i,jmax-1,1) 
           Bz (i,jmax  ,1) = Bz (i,jmax-1,1)  
           Vx (i,jmax  ,1) = Vx (i,jmax-1,1)  
           Vy (i,jmax  ,1) = Vy (i,jmax-1,1)  
           Vz (i,jmax  ,1) = Vz (i,jmax-1,1) 
           Ex (i,jmax  ,1) = Ex (i,jmax-1,1)
           Ey (i,jmax  ,1) = Ey (i,jmax-1,1)
           Ez (i,jmax  ,1) = Ez (i,jmax-1,1)
           q  (i,jmax  ,1) = q  (i,jmax-1,1)   
           p  (i,jmax  ,1) = p  (i,jmax-1,1)
           rho(i,jmax  ,1) = rho(i,jmax-1,1)

           psi(i,jmax+1,1) = psi(i,jmax  ,1) 
           phi(i,jmax+1,1) = phi(i,jmax  ,1) 
           Bx (i,jmax+1,1) = Bx (i,jmax  ,1)
           By (i,jmax+1,1) = By (i,jmax  ,1) 
           Bz (i,jmax+1,1) = Bz (i,jmax  ,1)  
           Vx (i,jmax+1,1) = Vx (i,jmax  ,1)  
           Vy (i,jmax+1,1) = Vy (i,jmax  ,1)  
           Vz (i,jmax+1,1) = Vz (i,jmax  ,1) 
           Ex (i,jmax+1,1) = Ex (i,jmax  ,1)
           Ey (i,jmax+1,1) = Ey (i,jmax  ,1)
           Ez (i,jmax+1,1) = Ez (i,jmax  ,1)
           q  (i,jmax+1,1) = q  (i,jmax  ,1)   
           p  (i,jmax+1,1) = p  (i,jmax  ,1)
           rho(i,jmax+1,1) = rho(i,jmax  ,1)


           psi(i,jmax+2,1) = psi(i,jmax+1,1) 
           phi(i,jmax+2,1) = phi(i,jmax+1,1) 
           Bx(i,jmax+2,1)  = Bx (i,jmax+1,1)
           By(i,jmax+2,1)  = By (i,jmax+1,1) 
           Bz(i,jmax+2,1)  = Bz (i,jmax+1,1)  
           Vx(i,jmax+2,1)  = Vx (i,jmax+1,1)  
           Vy(i,jmax+2,1)  = Vy (i,jmax+1,1)  
           Vz(i,jmax+2,1)  = Vz (i,jmax+1,1) 
           Ex(i,jmax+2,1)  = Ex (i,jmax+1,1)
           Ey(i,jmax+2,1)  = Ex (i,jmax+1,1)
           Ez(i,jmax+2,1)  = Ez (i,jmax+1,1)
           q(i,jmax+2,1)   = q  (i,jmax+1,1)   
           p(i,jmax+2,1)   = p  (i,jmax+1,1)
           rho(i,jmax+2,1) = rho(i,jmax+1,1)

 
!        end do
     end do

!$OMP END DO

    else

       write(*,*) "STOP: subroutine boundary d(-_-)b "
       write(*,*) "This combination of Dimension and Boundary is not correctly"
       write(*,*) "DIM =", DIM, "BOUND =", BOUND
       stop

    end if
   

    if (SLC == 1) then

       ! Copy boundary 1D

       if (DIM == 1) then
  !-----------------
  ! X Direction
  !-----------------

           ! Left Boundary
           ! -----------------

           sigma_loc( 0,1,1) = sigma_loc( 1,1,1)

           sigma_loc(-1,1,1) = sigma_loc( 0,1,1)

           sigma_loc(-2,1,1) = sigma_loc(-1,1,1)

           ! Right  Boundary
           ! -----------------

           sigma_loc(imax  ,1,1) = sigma_loc(imax-1,1,1) 

           sigma_loc(imax+1,1,1) = sigma_loc(imax  ,1,1) 

           sigma_loc(imax+2,1,1) = sigma_loc(imax+1,1,1)


        else if (DIM == 2) then

       ! Copy boundary 2D

  !-----------------
  ! X Direction
  !-----------------

!$OMP DO  

    do j=-2,jmax+2
!       do k=1,kmax

            ! Left Boundary
           ! -----------------

           sigma_loc( 0,j,1) = sigma_loc( 1,j,1) 

           sigma_loc(-1,j,1) = sigma_loc( 0,j,1) 

           sigma_loc(-2,j,1) = sigma_loc(-1,j,1) 

           ! Right  Boundary
           ! -----------------
           sigma_loc(imax  ,j,1) = sigma_loc(imax-1,j,1) 

           sigma_loc(imax+1,j,1) = sigma_loc(imax  ,j,1) 

           sigma_loc(imax+2,j,1) = sigma_loc(imax+1,j,1)


!        end do
     end do
 
!$OMP END DO


  !-----------------
  ! Y Direction
  !-----------------

!$OMP DO  

    do i=-2,imax+2
!       do k=1,kmax



           ! Left Boundary
           ! -----------------

           sigma_loc(i, 0,1) = sigma_loc(i, 1,1)

           sigma_loc(i,-1,1) = sigma_loc(i, 0,1)

           sigma_loc(i,-2,1) = sigma_loc(i,-1,1)

           ! Right  Boundary
           ! -----------------

           sigma_loc(i,jmax  ,1) = sigma_loc(i,jmax-1,1) 

           sigma_loc(i,jmax+1,1) = sigma_loc(i,jmax  ,1) 

           sigma_loc(i,jmax+2,1) = sigma_loc(i,jmax+1,1) 

 
!        end do
     end do

!$OMP END DO 
 
        end if

    end if

!$OMP END PARALLEL

   end subroutine boundary
