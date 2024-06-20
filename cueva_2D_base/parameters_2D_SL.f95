module parameters

    implicit none


    ! If you want to use MIRK first order choose order = 1 and MIRK = 1
    ! if you want MIRK second order choose order = 2 and MIRK = 2
    ! if you want to use IMEX schemes choose any order implemented and MIRK = 0

    integer,parameter  :: order           = 2  ! 1=MIRK_1st;2=MIRK_2nd or SSP2(222);3=SSP2-LUM;4=SSP3(433);5=BPR(353)
    integer,parameter  :: MIRK            = 2  ! use 1 for MIRK_1st 2 for MIRK_2nd and 0 for IMEX
    integer,parameter  :: DIM             = 2  ! use 1 for 1D and 2 for 2D
    integer,parameter  :: fix             = 0  ! use 1 for fix loop counter in 1D tests and 0 for 2D tests
    integer,parameter  :: TEST            = 27 ! look down the number assigned for each test
    integer,parameter  :: BOUND           = 2  ! use 1 for copy boundary, 2 for periodic, 3 for reconnection special case, and 4 perfect conductor TM
    integer,parameter  :: RESTART         = 0  ! use 1 for restart from back-up file and 0 for initial conditions
    integer,parameter  :: slop            = 1  ! use 1 for MinMod 2 for MCL and 3 for Superbee (in flux reconstruction)
    integer,parameter  :: limit           = 1  ! use 1 for MinMod 2 for MCL and 3 for Superbee (in primitive variables reconstruction )

    ! Recovery and reconstruction options
    integer,parameter  :: flux_solver     = 5  ! use 1 for Simpler Lax-Friederich, 2 for HLL reconstruction flux, 3 for HLLC_1, 4 for HLLC_2 , 5 for HLL + MP5, 6 for HLLC_1 + MP5 and 7 HLLC_2 + MP5

    integer,parameter  :: REC_PRIM        = 0  ! use 1 for YES recovery primitive in intermediate steeps and 0 for NOT recovery
    integer,parameter  :: REC_ELEC        = 0  ! use 1 for YES recovery electric field in cardano subrotine and 0 for NOT recovery
    integer,parameter  :: primitive_rec   = 0  ! use 1 for YES reconstruction of primitive variables and 0 for NOT
    integer,parameter  :: recons_method   = 0  ! use 1 for MP5 limiter and 0 for MCL or MinMod limiter
    integer,parameter  :: var_source_rec  = 0  ! use 1 for YES reconstruction of intermediadet variables in source term and 0 for NOT
    integer,parameter  :: source_rec      = 0  ! use 1 for YES reconstruction of source term and 0 for NOT

                                               ! This parameter only for MP5
    integer,parameter  :: flux_rec_mp5    = 1  ! use 1 for reconstruction of fluxes using mp5 and 0 for explicit formulation of fluxes using reconstruc mp5 conserved variables
 
    integer,parameter  :: conserved_to_primitive     = 1  ! use 1 for CARDANO solver and 0 for NEWTON-RAPHSON solver
    integer,parameter  :: root_method                = 1  ! use 1 for secant method or 0 for false point method to find the star velocity in HLLC solver
    
    !Factors for slop limiters and EGLM
    doubleprecision,parameter :: FAC        = 0.5d0  ! slop factor flux reconstruction
    doubleprecision,parameter :: FAC_SR     = 0.5d0  ! slop factor source reconstruction
    doubleprecision,parameter :: FAC_EGLM   = 0.0d0  ! EGLM factor
    doubleprecision           :: FAC_SOURCE, FAC_SOURCE_1, FAC_SOURCE_2         ! factor that multiply the source terms in HLLC Riemman solver
    !MOVIE
    integer,parameter         :: movie      = 0  ! use 1 for YES movie SCS and 0 for NOT


!__________ TEST NUMBER______________

!     ONE DIMENSIONAL TEST
!----------------------------------------------------------------------------------    
!     SCS                (1D) ------> TEST 1  ; BOUND 1 ; DIM 1
!     CPAW (Palenzuela)  (1D) ------> TEST 2  ; BOUND 2 ; DIM 1
!     RST  (Palenzuela)  (1D) ------> TEST 3  ; BOUND 1 ; DIM 1
!     RST  (Bucciantini) (1D) ------> TEST 4  ; BOUND 1 ; DIM 1
!     RST  (Aloy-Miranda)(1D) ------> TEST 5  ; BOUND 1 ; DIM 1
!     Reconnection       (1D) ------> TEST 6  ; BOUND 3 ; DIM 1
!     Tearing Mode       (1D) ------> TEST 7  ; BOUND 1 ; DIM 1
!     Magnetic Diffusion (1D) ------> TEST 8  ; BOUND 2 ; DIM 1 or BOUND 4
!     Shear Layer        (1D) ------> TEST 9  ; BOUND 2 ; DIM 1 or BOUND 4
!     CPAW(Apr classical)(1D) ------> TEST 10 ; BOUND 2 ; DIM 1
!----------------------------------------------------------------------------------    
!     TWO DIMENSIONAL TEST
!----------------------------------------------------------------------------------    
!     CE                         (2D) ------> TEST 20 ; BOUND 1 ; DIM 2
!     Rotor                      (2D) ------> TEST 21 ; BOUND 1 ; DIM 2
!     CS                         (2D) ------> TEST 22 ; BOUND 1 ; DIM 2
!     SJ                         (2D) ------> TEST 23 ; BOUND 4 ; DIM 2
!     CSI                        (2D) ------> TEST 24 ; BOUND 5 ; DIM 2
!     RMI                        (2D) ------> TEST 25 ; BOUND 6 ; DIM 2
!     MD                         (2D) ----->  TEST 26 ; BOUND 2 ; DIM 2
!     SL                         (2D) ----->  TEST 27 ; BOUND 2 ; DIM 2
!     Reconnection Aloy & Mimica (2D) ------> TEST 30 ; BOUND 3 ; DIM 2 using BOUND 6
!     Reconnection Zenitani      (2D) ------> TEST 31 ; BOUND 7 ; DIM 2 using BOUND 6
!     Tearing Mode               (2D) ------> TEST 32 ; BOUND 6 ; DIM 2 
!____________________________________


    doubleprecision,parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803D0
    doubleprecision,parameter :: THIRD = 0.33333333333333333333333333333333333333333333333333333333333333333333333333333333333D0


  !     Counters

  !     time ----> h
  !     x var ---> i
  !     y var ---> j
  !     z var ---> k
  !     auxiliar var ---> l,m
  !     Variables primit ---> n
  !     NewRaph ---> nr


    integer   :: h,i,j,k,l,m,n,nr,hrec, hbck, mm, nn
    integer   :: ai, aj
    
    ! max values of

    
    integer, parameter :: hini  = 1     !use 1 for initial start
    integer, parameter :: hmax  = 32001
    integer, parameter :: imax  = 32
    integer, parameter :: jmax  = 32    ! 1 for 1D
    integer, parameter :: kmax  = 1     ! 1 for 1D
    integer, parameter :: lmax  = order
    integer, parameter :: mmax  = lmax-1
    integer, parameter :: nmax  = 40
    integer, parameter :: nrmax = 40
    integer, parameter :: imed = floor(0.5d0*imax)
    integer, parameter :: jmed = floor(0.5d0*jmax)
    integer            :: ierror

    ! ------------------
    ! " RECONNECTION "
    ! -----------------
    !   Linquist Number
    doubleprecision,parameter ::  S_l      =  1.d6 !904.53d0 !1.d2 !1.d3 !
    !   Simulation domine
    doubleprecision,parameter ::  Lr       =  2.0d-2
    doubleprecision,parameter ::  a_tm     =  S_l**(-1./3.) * Lr
    doubleprecision,parameter ::  k_tm     =  100.d0 * pi
    doubleprecision,parameter ::  k_x      =  100.d0 * pi
    doubleprecision,parameter ::  k_y      =  100.d0 * pi
    doubleprecision,parameter ::  kap_tm   =  k_tm * a_tm
    doubleprecision,parameter ::  Lx       =  2.d0 * pi / k_x !40.d0 * a_tm    !200.d0   !
    doubleprecision,parameter ::  Ly       =  2.d0 * pi / k_y !0.9975d0 !20.d0 !
    !   Upstream Temperature
    doubleprecision,parameter :: theta     =  1.0d-2
    !   Thickness of the current layer
    doubleprecision,parameter :: delta_rec =  a_tm !1.d-3 !Longx/10.d0 !1.d-5 ! sqrt(Lr/(sigma_0 * VA_rec)) !Ly/10.d0
    doubleprecision,parameter :: delta_tm  =  100.d0 * pi 
    !   Amplitud TM perturbation
    doubleprecision,parameter :: epsilon_tm=  1.d-4

    doubleprecision           :: delta_tm_t, delta_rec_t, ekin 
    doubleprecision           :: max_val_bx, EdotB, by_imed, by_prom
    doubleprecision           :: sum_psi, sum_phi, sum_ex, sum_ey, sum_ez, sum_exey, sum_exez, sum_eyez
    doubleprecision           :: sum_vx, sum_vy, sum_ekin, max_val_ekin
    doubleprecision           :: sum_bx, sum_by, sum_bz, sum_bxbybz

    
    doubleprecision,parameter :: rho_up    =  1.d0
    doubleprecision,parameter :: B0        =  1.d-2   !1.d0 !sqrt(sigma_m * rho_up)
    doubleprecision,parameter :: V0        =  1.d-5   !1.d0 !sqrt(sigma_m * rho_up)
    doubleprecision,parameter :: P_up      =  1.d-4   !theta * rho_up !multiply by constant 5 ?
    doubleprecision,parameter :: P0        =  1.d-4 !B0**2/2
    doubleprecision,parameter :: alpha_rec =  0.5d0 * pi 
    !   Magnetization
    doubleprecision,parameter :: sigma_m   =  B0**2/rho_up
    doubleprecision,parameter :: beta_m    =  2.d0 * P0 / B0**2

    !   Alfven Velocity Tearing Mode
    doubleprecision,parameter :: VA_rec_tearing =  1.d-2 !1.d0 / sqrt(1.d0/sigma_m + 2.d0 * beta_m + 1.d0)    ! B0/sqrt(rho_up) ! 
    !   Alfven Velocity
    doubleprecision,parameter :: VA_rec    =  1.d-2 !sqrt(sigma_m/(1 + sigma_m)) !Zenitani et al. (2010)

    !   Conductivity
    doubleprecision,parameter :: sigma_0   =  1.d9 !S_l/(Lr * VA_rec_tearing) !
    doubleprecision,parameter :: sigma_1   =  sigma_0/1.d1
   ! Conductivity

    doubleprecision :: sigma    = 1.d9 !S_l/(Lr * VA_rec_tearing) !


    ! Adiabatic index

    doubleprecision,parameter :: gamma = 4.d0/3.d0 !2.d0 !  5.d0/3.d0 !      monatomic gas    !  diatomic gas  !  !cp/cv  !


    ! spatially localized conductivity SLC

    integer,parameter         :: SLC       = 0     ! use 1 for YES spatially localized conductivity  and 0 for constant conductivity
    integer,parameter         :: gamma_slc = 0



!!$    doubleprecision :: sigma_0  = 1.d6  
    doubleprecision :: sigma_L, sigma_R
    ! Scalar hyperbolic divergence cleaning

    doubleprecision,parameter :: kappa     = 1.d3
    doubleprecision,parameter :: kappa_psi = 1.d3   ! clean divergence factor divE
    doubleprecision,parameter :: kappa_phi = 1.d3   ! clean divergence factor divB

    ! MIRK  Coefficients  cm_4 < 0.5 * cm_1 < 0.5  cm_1 and cm_4 € (0,1)

    doubleprecision,parameter :: cm_1 = -0.05d0 !as ISA advice this number must be negative ¡¡¡¡
    doubleprecision,parameter :: cm_4 =  0.5d0 * (1.d0 - cm_1)**2 / cm_1
    doubleprecision           :: sigma_bar, sigma_dos_bar

    ! IMEX Coefficients of the explicit matrix for the "stiff term"

    doubleprecision,parameter :: const  = 1.d0 - 1.d0/sqrt(2.d0)
    doubleprecision,parameter :: delta  = const / (1-const)
    doubleprecision,parameter :: alpha1 = 0.24169426078821d0
    doubleprecision,parameter :: beta1  = 0.06042356519705d0
    doubleprecision,parameter :: etha1  = 0.12915286960590d0


    ! ------------------
    ! " ALFVEN WAVE "
    ! -----------------

    doubleprecision           :: VA
    doubleprecision,parameter :: ethaA =  1.d0
!!$    doubleprecision,parameter :: B0    =  1.1547d0 !0053838d0

    ! ------------------
    ! " RMI Constants and  RST_AlMi Constants"
    ! -----------------

    doubleprecision           :: Vs, Ws, W_a, W_b, Vx_a, Vx_b, Vy_a, Vy_b, Vz_a, Vz_b
    doubleprecision           :: rho_a, rho_b, p_a, p_b, ph_a, ph_b, h_a, h_b, b_a
    doubleprecision           :: Ex_a, Ey_a, Ez_a, Ex_b, Ey_b, Ez_b
    doubleprecision           :: Bx_a, By_a, Bz_a, Bx_b, By_b, Bz_b
    doubleprecision           :: Ex_b1, Ey_b1, Ez_b1, Bx_b1, By_b1, Bz_b1
    doubleprecision           :: Jx_a, Jy_a, Jz_a, Jx_b,Jy_b,Jz_b 
    doubleprecision           :: hh_a, q_a, q_b,psi_a, phi_a, psi_b, phi_b
    doubleprecision           :: J_inv, eps_a,  frontera, pos_shock, Ba2, Bb2, Ea2, Eb2, Vb2
    doubleprecision           :: TOL1,  ERR_BY,  ERR_B0,  ERR_B1,  ERR_B2, EdotV_a, EdotV_b   
    doubleprecision,parameter :: x_0    = 3.d0
    doubleprecision,parameter :: lambda = 2.5d0
    doubleprecision,parameter :: a_rmi  = 0.1d0


    !  Wave sppeds S+  !*********************************************************************************************

    doubleprecision,parameter                   :: const_sp = 0.d0 ! Use 1 for low speeds and 0 for relativistic cases
    doubleprecision                             :: Sp, c_sp, Sp_x, c_sp_x, Sp_y, c_sp_y

                       !*********************************************************************************************

    ! ------------------
    ! " HLLE VELOCITIES"
    ! -----------------

    doubleprecision,parameter :: s1 = -1.d0
    doubleprecision,parameter :: s2 =  1.d0


    ! ------------------
    ! " CYLINDRICAL EXPLOTION "
    ! -----------------
    doubleprecision,parameter :: alpha = 34.5387763949d0 !Palenzuela //! 52.0715658815d0 !Komissarov //! 49.5174377627d0 !Komissarov modificado //!
    doubleprecision,parameter :: beta  = 11.512925465d0  !Palenzuela //! 23.02585093d0   !Komissarov //!
    doubleprecision           :: r, posx, posy, x, omega


        !    Intervals:
    !
    !     time   --->  Delt
    !     x var  --->  Delx
    !     y var  --->  Dely
    !     z var  --->  Delz

    doubleprecision,parameter                   :: Longx = Lx
    doubleprecision,parameter                   :: Longy = Ly
    doubleprecision,parameter                   :: Delx = Longx  / real(imax)
    doubleprecision,parameter                   :: Dely = Longy  / real(jmax)
    doubleprecision,parameter                   :: Delz = Longy  / real(kmax)
    doubleprecision,parameter                   :: CFL  = 0.5d0

    doubleprecision,parameter                   :: Delt = CFL *  min(Delx, Dely)


  end module parameters
