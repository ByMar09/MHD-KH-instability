
  module threevectors

    use parameters, only: imax, jmax, kmax
    implicit none 
    !    save

  ! ------------------
  ! " spatially localized conductivity"
  ! -----------------

  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: sigma_loc


  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: ph,hhat

  ! ------------------
  ! " Variable reconstruction"
  ! -----------------

  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: psiL,phiL,BxL,ByL,BzL,ExL,EyL,EzL,VxL,VyL,VzL,qL,pL,rhoL
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: psiR,phiR,BxR,ByR,BzR,ExR,EyR,EzR,VxR,VyR,VzR,qR,pR,rhoR


  !     Vectorial Fields
  !     ------------
  !     Electric (Ex,Ey,Ez), Magnetic (Bx,By,Bz),  
  !     Velocity (Vx,Vy,Vz), Current density (Jx,Jy,Jz), 
  !     Mass flux (FD = \rho W V), Energy flux  (Ftau),
  !      Momentum (Sx,Sy,Sz)) 

  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: psi,phi
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: Ex,Ey,Ez
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: Bx,By,Bz
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: Exn,Eyn,Ezn !used for mirk2, see subroutine 20_electricfield.f95
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: Bxn,Byn,Bzn !used for mirk1, see subroutine 20_electricfield.f95
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: Vx,Vy,Vz
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: Sx,Sy,Sz
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: Jx,Jy,Jz
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: FDx,FDy,FDz  
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: Ftaux,Ftauy,Ftauz 
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: pmag, ptotal
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: divE_E, divB_B
  
  !     Momentum  Flux Tensor
  !     ------------
  !          | FSxx FSxy FSxz |
  !    FS =  | FSyx FSyy FSyz |
  !          | FSzx FSzy FSzz |

  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: FSxx,FSxy,FSxz
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: FSyx,FSyy,FSyz
  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: FSzx,FSzy,FSzz

  !     Scalar Fields
  !     ------------
  !     Density (\rho), Pressure (p), Mass (D= \rho W), Energy
  !      density (\tau), 
  !     Specific internal energy ( \epsilon), enthalpy (h=\rho (1-
  !     \epsilon)+p), 
  !     Electric field magnitude  (E^2),  Magnetic (B^2) , Velocity
  !      (V^2), Lorentz factor (W)

  doubleprecision, dimension(-6:imax+6,-6:jmax+6,1:kmax) :: rho,D,q,p,tau,epsiln,enthpy,E2,B2,V2,W

  !     Sums of final values
  !     ------------  

  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: psifinsum, phifinsum
  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: Bxfinsum, Byfinsum, Bzfinsum
  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: Exfinsum, Eyfinsum, Ezfinsum
  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: qfinsum,  Dfinsum,taufinsum
  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: Sxfinsum, Syfinsum, Szfinsum



  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: psiintsum, phiintsum
  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: Bxintsum, Byintsum, Bzintsum
  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: Exastsum, Eyastsum, Ezastsum
  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: qintsum, Dintsum, tauintsum
  doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax)  :: Sxintsum, Syintsum, Szintsum

  end module threevectors
