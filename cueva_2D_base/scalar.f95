  module scalar

    use parameters

    implicit none 


!!$!Tomek peiodic boubdary conditions (stencil refere a number of phantom cells)
    integer,parameter :: stencil = 6
!!$!Tomek

    !indices para acotar region media RDTM
    !=====================================
    
    integer           :: imed1, imed2
    !indices para fix-up density and presion subroutine 12_varprimitive.f95
    !=====================================
    integer           :: contador_r, contador_p, contador_v 
    !=====================================

    ! variables for drift velocity v_drift = (E X B) / |B^2|, ussing in subroutine 12_varprimitive.f95
    ! when appear superlumilal velocities and great values of magnetization sigma_m = |B^2|/rho > 10

    doubleprecision           :: ErotB, v_drift, B2_drift, W_drift
    
    !     ------------
    !     "NEWTON RAPHSON" 
    !     ------------

    doubleprecision                             :: pnr, pnro
    doubleprecision                             :: Exnr,Eynr,Eznr
    doubleprecision                             :: epsilonnr,rhonr,Cs2 
    doubleprecision                             :: Wnr,E2B2,V2nr
    doubleprecision                             :: CHA,CHA1,CHA2,CHA3
    doubleprecision,dimension(1:nmax)           :: Vxnr,Vynr,Vznr
    doubleprecision                             :: TOL   ! Toro E.F (4.45) originalmente 10**(-6)

    ! ------------------
    ! " CARDANO "
    ! -----------------

    integer                                     :: nrw, n_elec, n_str, l_swap, t_init
    doubleprecision                             :: gamma_1
    doubleprecision                             :: V2c, V2int, V2mp5,  V2mp5_x, V2mp5_y
    doubleprecision                             :: V2mp5_x_L, V2mp5_y_L, V2mp5_x_R, V2mp5_y_R
    doubleprecision                             :: C_1, C_2
    doubleprecision                             :: Wnro
    doubleprecision                             :: a_4, a_3, a_2, a_1, a_0
    doubleprecision                             :: b_3, b_2, b_1, b_0
    doubleprecision                             :: d_2, d_1, d_0
    doubleprecision                             :: e  , f  , g  , z_1
    doubleprecision                             :: s  , t  , y_1, y_2
    doubleprecision                             :: AA_1, BB_1, det, Q_1, Q_2, R_1, T_1
    doubleprecision                             :: s_1a, s_1b, s_1c
    doubleprecision                             :: t_1a, t_1b, t_1c
    doubleprecision                             :: z_1a, z_1b, z_1c
    doubleprecision                             :: xx1

    doubleprecision                             :: Vx_elec, Vy_elec, Vz_elec
    doubleprecision                             :: Vx_elec_0, Vy_elec_0, Vz_elec_0
    doubleprecision                             :: CHA_ELEC_1, CHA_ELEC_2, CHA_ELEC_3
    doubleprecision                             :: Edotv, vrotB_x, vrotB_y, vrotB_z, Bdotv
    doubleprecision                             :: vrotB_n_x, vrotB_n_y, vrotB_n_z         !used for mirk1 see subroutine 20.electricfield.f95
    doubleprecision                             :: Ex_elec,Ey_elec, Ez_elec 
    doubleprecision                             :: Ex_elec_0,Ey_elec_0, Ez_elec_0 
    doubleprecision                             :: Bx_elec, By_elec, Bz_elec
    doubleprecision                             :: rho_elec, D_elec, p_elec
    doubleprecision                             :: Sx_elec, Sy_elec, Sz_elec
    doubleprecision                             :: tau_elec, enthpy_elec, W_elec, W_elec_2 



    ! ------------------
    ! " MATRICES AND VECTORS BUTCHER "
    ! -----------------

    doubleprecision,dimension(1:order)       :: a, matx 
    doubleprecision,dimension(1:8)           :: wb, wbt
    doubleprecision,dimension(1:8,1:8)       :: Ab, Abt

    ! ------------------
    ! " Variable reconstruction"
    ! -----------------

    doubleprecision           :: psim1,phim1,Bxm1,Bym1,Bzm1,Exm1,Eym1,Ezm1,Vxm1,Vym1,Vzm1,qm1,pm1,rhom1
    doubleprecision           :: psim ,phim ,Bxm ,Bym ,Bzm ,Exm ,Eym ,Ezm ,Vxm ,Vym ,Vzm ,qm ,pm ,rhom 
    doubleprecision           :: psip1,phip1,Bxp1,Byp1,Bzp1,Exp1,Eyp1,Ezp1,Vxp1,Vyp1,Vzp1,qp1,pp1,rhop1
    doubleprecision           :: psip2,phip2,Bxp2,Byp2,Bzp2,Exp2,Eyp2,Ezp2,Vxp2,Vyp2,Vzp2,qp2,pp2,rhop2

    doubleprecision           :: sppsi1,spphi1,spBx1,spBy1,spBz1,spEx1,spEy1,spEz1,spVx1,spVy1,spVz1,spq1,spp1,sprho1
    doubleprecision           :: smpsi1,smphi1,smBx1,smBy1,smBz1,smEx1,smEy1,smEz1,smVx1,smVy1,smVz1,smq1,smp1,smrho1
    doubleprecision           :: sppsi2,spphi2,spBx2,spBy2,spBz2,spEx2,spEy2,spEz2,spVx2,spVy2,spVz2,spq2,spp2,sprho2
    doubleprecision           :: smpsi2,smphi2,smBx2,smBy2,smBz2,smEx2,smEy2,smEz2,smVx2,smVy2,smVz2,smq2,smp2,smrho2

    doubleprecision           :: psipp,phipp,Bxpp,Bypp,Bzpp,Expp,Eypp,Ezpp,Vxpp,Vypp,Vzpp,qpp,ppp,rhopp
    doubleprecision           :: psimm,phimm,Bxmm,Bymm,Bzmm,Exmm,Eymm,Ezmm,Vxmm,Vymm,Vzmm,qmm,pmm,rhomm


    ! ------------------
    ! " HLL Variables and Flux  reconstruction"
    ! -----------------

    doubleprecision         :: E2B2int_L, E2B2int_R, p_L, p_R, rho_L, rho_R, epsilon_L, epsilon_R
    doubleprecision         :: psiint_L, psiint_R, phiint_L, phiint_R, Bxint_L, Bxint_R, Byint_L, Byint_R, Bzint_L, Bzint_R
    doubleprecision         :: Exint_L, Exint_R, Eyint_L, Eyint_R, Ezint_L, Ezint_R, qint_L, qint_R, DDint_L, DDint_R
    doubleprecision         :: tauint_L, tauint_R, Sxint_L, Sxint_R, Syint_L, Syint_R, Szint_L, Szint_R
    doubleprecision         :: Vx_L, Vx_R, Vy_L, Vy_R, Vz_L, Vz_R, W_L, W_R, Edotv_L, Edotv_R, vrotB_x_L, vrotB_y_L, vrotB_z_L
    doubleprecision         :: Ux_L, Ux_R, Uy_L, Uy_R, Uz_L, Uz_R
    doubleprecision         :: vrotB_x_R, vrotB_y_R, vrotB_z_R, Jxint_L, Jxint_R, Jyint_L, Jyint_R, Jzint_L, Jzint_R
    doubleprecision         :: FDxint_L, FDxint_R, FDyint_L, FDyint_R, FDzint_L, FDzint_R  
    doubleprecision         :: enthpy_L, enthpy_R, Ftauxint_L, Ftauxint_R, Ftauyint_L, Ftauyint_R, Ftauzint_L, Ftauzint_R
    doubleprecision         :: FSxxint_L, FSxxint_R, FSxyint_L, FSxyint_R, FSxzint_L, FSxzint_R, FSyxint_L, FSyxint_R, FSyyint_L
    doubleprecision         :: FSyyint_R, FSyzint_L, FSyzint_R, FSzxint_L, FSzxint_R, FSzyint_L, FSzyint_R, FSzzint_L, FSzzint_R 
    doubleprecision         :: Exint_L_m, Exint_R_m, Eyint_L_m, Eyint_R_m, Ezint_L_m, Ezint_R_m
    doubleprecision         :: Bxint_L_m, Bxint_R_m, Byint_L_m, Byint_R_m, Bzint_L_m, Bzint_R_m

    doubleprecision         :: p_x_L, p_x_R, rho_x_L, rho_x_R, epsilon_x_L, epsilon_x_R
    doubleprecision         :: psiint_x_L, psiint_x_R, phiint_x_L, phiint_x_R 
    doubleprecision         :: Bxint_x_L, Bxint_x_R, Byint_x_L, Byint_x_R, Bzint_x_L, Bzint_x_R
    doubleprecision         :: Exint_x_L, Exint_x_R, Eyint_x_L, Eyint_x_R, Ezint_x_L, Ezint_x_R
    doubleprecision         :: qint_x_L, qint_x_R, DDint_x_L, DDint_x_R,tauint_x_L, tauint_x_R
    doubleprecision         :: Sxint_x_L, Sxint_x_R, Syint_x_L, Syint_x_R, Szint_x_L, Szint_x_R
    doubleprecision         :: Vx_x_L, Vx_x_R, Vy_x_L, Vy_x_R, Vz_x_L, Vz_x_R
    doubleprecision         :: Exint_x_L_m, Exint_x_R_m, Eyint_x_L_m, Eyint_x_R_m, Ezint_x_L_m, Ezint_x_R_m
    doubleprecision         :: Bxint_x_L_m, Bxint_x_R_m, Byint_x_L_m, Byint_x_R_m, Bzint_x_L_m, Bzint_x_R_m

    doubleprecision         :: p_y_L, p_y_R, rho_y_L, rho_y_R, epsilon_y_L, epsilon_y_R
    doubleprecision         :: psiint_y_L, psiint_y_R, phiint_y_L, phiint_y_R
    doubleprecision         :: Bxint_y_L, Bxint_y_R, Byint_y_L, Byint_y_R, Bzint_y_L, Bzint_y_R
    doubleprecision         :: Exint_y_L, Exint_y_R, Eyint_y_L, Eyint_y_R, Ezint_y_L, Ezint_y_R
    doubleprecision         :: qint_y_L, qint_y_R, DDint_y_L, DDint_y_R,tauint_y_L, tauint_y_R
    doubleprecision         :: Sxint_y_L, Sxint_y_R, Syint_y_L, Syint_y_R, Szint_y_L, Szint_y_R
    doubleprecision         :: Vx_y_L, Vx_y_R, Vy_y_L, Vy_y_R, Vz_y_L, Vz_y_R
    doubleprecision         :: Exint_y_L_m, Exint_y_R_m, Eyint_y_L_m, Eyint_y_R_m, Ezint_y_L_m, Ezint_y_R_m
    doubleprecision         :: Bxint_y_L_m, Bxint_y_R_m, Byint_y_L_m, Byint_y_R_m, Bzint_y_L_m, Bzint_y_R_m
    doubleprecision         :: enthpy_x_L, enthpy_x_R, enthpy_y_L, enthpy_y_R
    doubleprecision         :: Edotv_x_L, Edotv_x_R, vrotB_x_x_L, vrotB_y_x_L, vrotB_z_x_L, vrotB_x_x_R, vrotB_y_x_R, vrotB_z_x_R  
    doubleprecision         :: Edotv_y_L, Edotv_y_R, vrotB_x_y_L, vrotB_y_y_L, vrotB_z_y_L, vrotB_x_y_R, vrotB_y_y_R, vrotB_z_y_R  
    doubleprecision         :: E2B2int_x_L, E2B2int_x_R, E2B2int_y_L, E2B2int_y_R
    doubleprecision         :: W_x_L, W_x_R, W_y_L, W_y_R
    doubleprecision         :: Vx_x, Vy_x, Vz_x, Vx_y, Vy_y, Vz_y
    doubleprecision         :: epsilon_x, epsilon_y
    doubleprecision         :: Vx_rec, Vy_rec, Vz_rec 
    doubleprecision         :: p_rec, rho_rec, epsilon_rec
 
  doubleprecision           :: Ex_sour_rec_x_L, Ex_sour_rec_x_R, Ey_sour_rec_x_L, Ey_sour_rec_x_R, Ez_sour_rec_x_L, Ez_sour_rec_x_R
  doubleprecision           :: Ex_sour_rec_y_L, Ex_sour_rec_y_R, Ey_sour_rec_y_L, Ey_sour_rec_y_R, Ez_sour_rec_y_L, Ez_sour_rec_y_R
  doubleprecision           :: Ex_sour_rec_L, Ex_sour_rec_R, Ey_sour_rec_L, Ey_sour_rec_R, Ez_sour_rec_L, Ez_sour_rec_R

  doubleprecision           :: Cs_R, Cs_L, V2_R, V2_L, lambda_H_p_R, lambda_H_p_L, lambda_H_m_R, lambda_H_m_L
  doubleprecision           :: Cs_x_R, Cs_x_L, V2_x_R, V2_x_L, lambda_H_p_x_R, lambda_H_p_x_L, lambda_H_m_x_R, lambda_H_m_x_L
  doubleprecision           :: Cs_y_R, Cs_y_L, V2_y_R, V2_y_L, lambda_H_p_y_R, lambda_H_p_y_L, lambda_H_m_y_R, lambda_H_m_y_L
  doubleprecision           :: sx1, sx2, sy1, sy2, ux_p_m, ux_m_m, uy_p_m, uy_m_m
  doubleprecision           :: W_H_p_L, W_H_p_R, W_H_m_L, W_H_m_R, W2m_p, W2m_m
  doubleprecision           :: W_H_p_x_L, W_H_p_x_R, W_H_m_x_L, W_H_m_x_R
  doubleprecision           :: W_H_p_y_L, W_H_p_y_R, W_H_m_y_L, W_H_m_y_R
  
    ! ------------------
    ! " HLLC Variables and Flux  reconstruction"
    ! -----------------

  doubleprecision           :: E2B2str_L, E2B2str_R, pstr_L, pstr_R, rhostr_L, rhostr_R, epsilonstr_L, epsilonstr_R, eps, eps1, eps2
  doubleprecision           :: E2_perp, B2_perp, rt1, rt2, E2_perp_x, B2_perp_x, rt1_x, rt2_x, E2_perp_y, B2_perp_y, rt1_y, rt2_y
   
  doubleprecision           :: psistr_L, psistr_R, phistr_L, phistr_R, Bxstr_L, Bxstr_R, Bystr_L, Bystr_R, Bzstr_L, Bzstr_R
  doubleprecision           :: Exstr_L, Exstr_R, Eystr_L, Eystr_R, Ezstr_L, Ezstr_R, qstr_L, qstr_R, DDstr_L, DDstr_R
  doubleprecision           :: taustr_L, taustr_R, Sxstr_L, Sxstr_R, Systr_L, Systr_R, Szstr_L, Szstr_R

  doubleprecision           :: psistr_x_L, psistr_x_R, phistr_x_L, phistr_x_R
  doubleprecision           :: Bxstr_x_L, Bxstr_x_R, Bystr_x_L, Bystr_x_R, Bzstr_x_L, Bzstr_x_R
  doubleprecision           :: Exstr_x_L, Exstr_x_R, Eystr_x_L, Eystr_x_R, Ezstr_x_L, Ezstr_x_R
  doubleprecision           :: qstr_x_L, qstr_x_R, DDstr_x_L, DDstr_x_R
  doubleprecision           :: taustr_x_L, taustr_x_R, Sxstr_x_L, Sxstr_x_R
  doubleprecision           :: Systr_x_L, Systr_x_R, Szstr_x_L, Szstr_x_R

  doubleprecision           :: psistr_y_L, psistr_y_R, phistr_y_L, phistr_y_R
  doubleprecision           :: Bxstr_y_L, Bxstr_y_R, Bystr_y_L, Bystr_y_R, Bzstr_y_L, Bzstr_y_R
  doubleprecision           :: Exstr_y_L, Exstr_y_R, Eystr_y_L, Eystr_y_R, Ezstr_y_L, Ezstr_y_R
  doubleprecision           :: qstr_y_L, qstr_y_R, DDstr_y_L, DDstr_y_R
  doubleprecision           :: taustr_y_L, taustr_y_R, Sxstr_y_L, Sxstr_y_R
  doubleprecision           :: Systr_y_L, Systr_y_R, Szstr_y_L, Szstr_y_R
    
    
  doubleprecision           :: Vxstr_L, Vxstr_R, Vystr_L, Vystr_R, Vzstr_L, Vzstr_R, Wstr_L, Wstr_R, Edotv_str_L, Edotv_str_R 
  doubleprecision           :: Vx_str0_L, Vx_str0_R, Vy_str0_L, Vy_str0_R, Vz_str0_L, Vz_str0_R
  doubleprecision           :: vrotB_x_str_L, vrotB_y_str_L, vrotB_z_str_L, vrotB_x_str_R, vrotB_y_str_R, vrotB_z_str_R
  doubleprecision           :: ErotB_str_x_L, ErotB_str_x_R, ErotB_str_y_L, ErotB_str_y_R, ErotB_str_z_L, ErotB_str_z_R
  doubleprecision           :: ErotB_str_x_x_L, ErotB_str_x_x_R, ErotB_str_y_x_L, ErotB_str_y_x_R, ErotB_str_z_x_L, ErotB_str_z_x_R
  doubleprecision           :: ErotB_str_x_y_L, ErotB_str_x_y_R, ErotB_str_y_y_L, ErotB_str_y_y_R, ErotB_str_z_y_L, ErotB_str_z_y_R
  doubleprecision           :: Jxstr_L, Jxstr_R, Jystr_L, Jystr_R, Jzstr_L, Jzstr_R
  doubleprecision           :: FDxstr_L, FDxstr_R, FDystr_L, FDystr_R, FDzstr_L, FDzstr_R  
  doubleprecision           :: enthpystr_L, enthpystr_R, Ftauxstr_L, Ftauxstr_R, Ftauystr_L, Ftauystr_R, Ftauzstr_L, Ftauzstr_R
  doubleprecision           :: FSxxstr_L, FSxxstr_R, FSxystr_L, FSxystr_R, FSxzstr_L, FSxzstr_R, FSyxstr_L, FSyxstr_R, FSyystr_L
  doubleprecision           :: FSyystr_R, FSyzstr_L, FSyzstr_R, FSzxstr_L, FSzxstr_R, FSzystr_L, FSzystr_R, FSzzstr_L, FSzzstr_R 
  doubleprecision           :: Exstr_L_m, Exstr_R_m, Eystr_L_m, Eystr_R_m, Ezstr_L_m, Ezstr_R_m
  doubleprecision           :: Bxstr_L_m, Bxstr_R_m, Bystr_L_m, Bystr_R_m, Bzstr_L_m, Bzstr_R_m
  doubleprecision           :: CHA_Vx_L, CHA_Vy_L, CHA_Vz_L, CHA_Vx_R, CHA_Vy_R, CHA_Vz_R
  doubleprecision           :: Omega1_L, Omega2_L, Omega3_L, Omega4_L, Omega5_L
  doubleprecision           :: Omega1_R, Omega2_R, Omega3_R, Omega4_R, Omega5_R
  doubleprecision           :: Omega_1 , Omega_2 , Omega_3 , Omega_4 , Omega_5
  doubleprecision           :: DOmega_1 , DOmega_2 , DOmega_3 , DOmega_4 , DOmega_5 
    
  doubleprecision           :: a_hll, b_hll, c_hll, det_hll, q_root, q_root1, V_root_1, V_root_2, E2B2_str, E2B2_str_x, E2B2_str_y
  doubleprecision           :: a_hll_x, b_hll_x, c_hll_x, a_hll_y, b_hll_y, c_hll_y
  doubleprecision           :: psi_hll, phi_hll, Bx_hll, By_hll, Bz_hll, Ex_hll, Ey_hll, Ez_hll 
  doubleprecision           :: q_hll, DD_hll, tau_hll, Sx_hll, Sy_hll, Sz_hll, psi_hll_flux, phi_hll_flux, Bx_hll_flux
  doubleprecision           :: By_hll_flux, Bz_hll_flux, Ex_hll_flux, Ey_hll_flux, Ez_hll_flux 
  doubleprecision           :: q_hll_flux, DD_hll_flux, tau_hll_flux, Sx_hll_flux, Sy_hll_flux, Sz_hll_flux
    
  doubleprecision           :: psi_x_hll, phi_x_hll, Bx_x_hll, By_x_hll, Bz_x_hll, Ex_x_hll, Ey_x_hll, Ez_x_hll 
  doubleprecision           :: q_x_hll, DD_x_hll, tau_x_hll, Sx_x_hll, Sy_x_hll, Sz_x_hll
  doubleprecision           :: psi_x_hll_flux, phi_x_hll_flux, Bx_x_hll_flux
  doubleprecision           :: By_x_hll_flux, Bz_x_hll_flux, Ex_x_hll_flux, Ey_x_hll_flux, Ez_x_hll_flux 
  doubleprecision           :: q_x_hll_flux, DD_x_hll_flux, tau_x_hll_flux, Sx_x_hll_flux, Sy_x_hll_flux, Sz_x_hll_flux

  doubleprecision           :: psi_y_hll, phi_y_hll, Bx_y_hll, By_y_hll, Bz_y_hll, Ex_y_hll, Ey_y_hll, Ez_y_hll 
  doubleprecision           :: q_y_hll, DD_y_hll, tau_y_hll, Sx_y_hll, Sy_y_hll, Sz_y_hll
  doubleprecision           :: psi_y_hll_flux, phi_y_hll_flux, Bx_y_hll_flux
  doubleprecision           :: By_y_hll_flux, Bz_y_hll_flux, Ex_y_hll_flux, Ey_y_hll_flux, Ez_y_hll_flux 
  doubleprecision           :: q_y_hll_flux, DD_y_hll_flux, tau_y_hll_flux, Sx_y_hll_flux, Sy_y_hll_flux, Sz_y_hll_flux
    
  doubleprecision           :: pstr, pstr_x, pstr_y, Vxstr, Vystr, Vzstr, Wstr, qstr, Exstr, Eystr, Ezstr, Bxstr, Bystr, Bzstr
  doubleprecision           :: phistr, psistr, ErotB_str_x, ErotB_str_y, ErotB_str_z
  doubleprecision           :: Omega1_str, Omega2_str, Omega3_str, Omega4_str, Omega5_str
    
  doubleprecision           :: psistr_flux, phistr_flux, Exstr_flux, Eystr_flux, Ezstr_flux, Bxstr_flux, Bystr_flux, Bzstr_flux
  doubleprecision           :: psistr_flux_L, phistr_flux_L, Exstr_flux_L, Eystr_flux_L, Ezstr_flux_L
  doubleprecision           :: Bxstr_flux_L, Bystr_flux_L, Bzstr_flux_L
  doubleprecision           :: psistr_flux_R, phistr_flux_R, Exstr_flux_R, Eystr_flux_R, Ezstr_flux_R
  doubleprecision           :: Bxstr_flux_R, Bystr_flux_R, Bzstr_flux_R
    
  doubleprecision           :: psistr_x_flux, phistr_x_flux, Exstr_x_flux, Eystr_x_flux, Ezstr_x_flux
  doubleprecision           :: Bxstr_x_flux, Bystr_x_flux, Bzstr_x_flux
  doubleprecision           :: psistr_x_flux_L, phistr_x_flux_L, Exstr_x_flux_L, Eystr_x_flux_L, Ezstr_x_flux_L
  doubleprecision           :: Bxstr_x_flux_L, Bystr_x_flux_L, Bzstr_x_flux_L
  doubleprecision           :: psistr_x_flux_R, phistr_x_flux_R, Exstr_x_flux_R, Eystr_x_flux_R, Ezstr_x_flux_R
  doubleprecision           :: Bxstr_x_flux_R, Bystr_x_flux_R, Bzstr_x_flux_R

  doubleprecision           :: psistr_y_flux, phistr_y_flux, Exstr_y_flux, Eystr_y_flux, Ezstr_y_flux
  doubleprecision           :: Bxstr_y_flux, Bystr_y_flux, Bzstr_y_flux
  doubleprecision           :: psistr_y_flux_L, phistr_y_flux_L, Exstr_y_flux_L, Eystr_y_flux_L, Ezstr_y_flux_L
  doubleprecision           :: Bxstr_y_flux_L, Bystr_y_flux_L, Bzstr_y_flux_L
  doubleprecision           :: psistr_y_flux_R, phistr_y_flux_R, Exstr_y_flux_R, Eystr_y_flux_R, Ezstr_y_flux_R
  doubleprecision           :: Bxstr_y_flux_R, Bystr_y_flux_R, Bzstr_y_flux_R
    
  doubleprecision           :: Fpsi_L, Fpsi_R, Fphi_L, Fphi_R, FEx_L, FEx_R, FEy_L, FEy_R, FEz_L, FEz_R
  doubleprecision           :: FBx_L, FBx_R, FBy_L, FBy_R, FBz_L, FBz_R


    ! Sources
  doubleprecision           :: Exsource_L, Exsource_R, Eysource_L, Eysource_R, Ezsource_L, Ezsource_R
    

    ! ------------------
    ! " MP5 FLUX LEFT and RIGHT VARIABLES
    ! -----------------

 !For 1D reconstruction in x direction

    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: E2B2int_mp5_L, E2B2int_mp5_R, p_mp5_L, p_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: psiint_mp5_L, psiint_mp5_R, phiint_mp5_L, phiint_mp5_R 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Bxint_mp5_L, Bxint_mp5_R, Byint_mp5_L 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Byint_mp5_R, Bzint_mp5_L, Bzint_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Exint_mp5_L, Exint_mp5_R, Eyint_mp5_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Eyint_mp5_R, Ezint_mp5_L, Ezint_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Bxint_mp5_L_m, Bxint_mp5_R_m, Byint_mp5_L_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Byint_mp5_R_m, Bzint_mp5_L_m, Bzint_mp5_R_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Exint_mp5_L_m, Exint_mp5_R_m, Eyint_mp5_L_m 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Eyint_mp5_R_m, Ezint_mp5_L_m, Ezint_mp5_R_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: tauint_mp5_L, tauint_mp5_R, qint_mp5_L, qint_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: DDint_mp5_L, DDint_mp5_R, rho_mp5_L, rho_mp5_R 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Sxint_mp5_L, Sxint_mp5_R, Syint_mp5_L, Syint_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Szint_mp5_L, Szint_mp5_R, Vz_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Vx_mp5_L, Vx_mp5_R, Vy_mp5_L, Vy_mp5_R, Vz_mp5_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: W_mp5_L, W_mp5_R, Edotv_mp5_L, Edotv_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: vrotB_x_mp5_L, vrotB_y_mp5_L,vrotB_x_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: vrotB_y_mp5_R, vrotB_z_mp5_R, vrotB_z_mp5_L   
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Jxint_mp5_L, Jxint_mp5_R, Jyint_mp5_L, Jyint_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Jzint_mp5_L, Jzint_mp5_R  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: enthpy_mp5_L, enthpy_mp5_R, epsilon_mp5_L, epsilon_mp5_R


    

! For 2D reconstruction in x direction

    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: E2B2int_mp5_x_L, E2B2int_mp5_x_R, p_mp5_x_L, p_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: psiint_mp5_x_L, psiint_mp5_x_R, phiint_mp5_x_L, phiint_mp5_x_R 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Bxint_mp5_x_L, Bxint_mp5_x_R, Byint_mp5_x_L 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Byint_mp5_x_R, Bzint_mp5_x_L, Bzint_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Exint_mp5_x_L, Exint_mp5_x_R, Eyint_mp5_x_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Eyint_mp5_x_R, Ezint_mp5_x_L, Ezint_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Bxint_mp5_x_L_m, Bxint_mp5_x_R_m, Byint_mp5_x_L_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Byint_mp5_x_R_m, Bzint_mp5_x_L_m, Bzint_mp5_x_R_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Exint_mp5_x_L_m, Exint_mp5_x_R_m, Eyint_mp5_x_L_m 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Eyint_mp5_x_R_m, Ezint_mp5_x_L_m, Ezint_mp5_x_R_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: tauint_mp5_x_L, tauint_mp5_x_R, qint_mp5_x_L, qint_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: DDint_mp5_x_L, DDint_mp5_x_R, rho_mp5_x_L, rho_mp5_x_R 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Sxint_mp5_x_L, Sxint_mp5_x_R, Syint_mp5_x_L, Syint_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Szint_mp5_x_L, Szint_mp5_x_R, Vz_mp5_x_R, Uz_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Vx_mp5_x_L, Vx_mp5_x_R, Vy_mp5_x_L, Vy_mp5_x_R, Vz_mp5_x_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Ux_mp5_x_L, Ux_mp5_x_R, Uy_mp5_x_L, Uy_mp5_x_R, Uz_mp5_x_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: W_mp5_x_L, W_mp5_x_R, Edotv_mp5_x_L, Edotv_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: vrotB_x_mp5_x_L, vrotB_y_mp5_x_L,vrotB_x_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: vrotB_y_mp5_x_R, vrotB_z_mp5_x_R, vrotB_z_mp5_x_L   
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Jxint_mp5_x_L, Jxint_mp5_x_R, Jyint_mp5_x_L, Jyint_mp5_x_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Jzint_mp5_x_L, Jzint_mp5_x_R  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: enthpy_mp5_x_L, enthpy_mp5_x_R, epsilon_mp5_x_L, epsilon_mp5_x_R

! For 2D reconstruction in y direction

    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: E2B2int_mp5_y_L, E2B2int_mp5_y_R, p_mp5_y_L, p_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: psiint_mp5_y_L, psiint_mp5_y_R, phiint_mp5_y_L, phiint_mp5_y_R 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Bxint_mp5_y_L, Bxint_mp5_y_R, Byint_mp5_y_L 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Byint_mp5_y_R, Bzint_mp5_y_L, Bzint_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Exint_mp5_y_L, Exint_mp5_y_R, Eyint_mp5_y_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Eyint_mp5_y_R, Ezint_mp5_y_L, Ezint_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Bxint_mp5_y_L_m, Bxint_mp5_y_R_m, Byint_mp5_y_L_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Byint_mp5_y_R_m, Bzint_mp5_y_L_m, Bzint_mp5_y_R_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Exint_mp5_y_L_m, Exint_mp5_y_R_m, Eyint_mp5_y_L_m 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Eyint_mp5_y_R_m, Ezint_mp5_y_L_m, Ezint_mp5_y_R_m
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: tauint_mp5_y_L, tauint_mp5_y_R, qint_mp5_y_L, qint_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: DDint_mp5_y_L, DDint_mp5_y_R, rho_mp5_y_L, rho_mp5_y_R 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Sxint_mp5_y_L, Sxint_mp5_y_R, Syint_mp5_y_L, Syint_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Szint_mp5_y_L, Szint_mp5_y_R, Vz_mp5_y_R, Uz_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Vx_mp5_y_L, Vx_mp5_y_R, Vy_mp5_y_L, Vy_mp5_y_R, Vz_mp5_y_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Ux_mp5_y_L, Ux_mp5_y_R, Uy_mp5_y_L, Uy_mp5_y_R, Uz_mp5_y_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: W_mp5_y_L, W_mp5_y_R, Edotv_mp5_y_L, Edotv_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: vrotB_x_mp5_y_L, vrotB_y_mp5_y_L,vrotB_x_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: vrotB_y_mp5_y_R, vrotB_z_mp5_y_R, vrotB_z_mp5_y_L   
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Jxint_mp5_y_L, Jxint_mp5_y_R, Jyint_mp5_y_L, Jyint_mp5_y_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Jzint_mp5_y_L, Jzint_mp5_y_R  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: enthpy_mp5_y_L, enthpy_mp5_y_R, epsilon_mp5_y_L, epsilon_mp5_y_R


! Fluxes

    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: FDxint_mp5_L, FDxint_mp5_R, FDyint_mp5_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: FDyint_mp5_R, FDzint_mp5_L, FDzint_mp5_R  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Ftauxint_mp5_L, Ftauxint_mp5_R, Ftauyint_mp5_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Ftauyint_mp5_R, Ftauzint_mp5_L, Ftauzint_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: FSxxint_mp5_L, FSxxint_mp5_R, FSxyint_mp5_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: FSxyint_mp5_R, FSxzint_mp5_L, FSxzint_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: FSyxint_mp5_L, FSyxint_mp5_R, FSyyint_mp5_L 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: FSyyint_mp5_R, FSyzint_mp5_L, FSyzint_mp5_R
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: FSzxint_mp5_L, FSzxint_mp5_R, FSzyint_mp5_L
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: FSzyint_mp5_R, FSzzint_mp5_L, FSzzint_mp5_R


! Sources
    doubleprecision,dimension(-6:imax+6,-6:jmax+6) :: Exsour_L, Exsour_R, Eysour_L, Eysour_R, Ezsour_L, Ezsour_R
    

! max val for TM 2D

    doubleprecision,dimension(0:imax,0:jmax) ::max_bx, max_phi, max_ekin,max_vx, max_va

! Variables false-position method

    integer,parameter         ::MAXIT=1000
    integer                   ::jsec,flag
    doubleprecision           ::rtsec,rtflsp,x1,x2,xacc
    doubleprecision           ::del,dx,ftemp,fh,fl,swap,xh,xl
    doubleprecision           ::f_m, f_0, f_p
    doubleprecision           ::vx_m, vx_0, vx_p, del_v, vgg

    ! ------------------
    ! " MP5 and MP7 and MP9 "
    ! -----------------

    doubleprecision, dimension(-6:imax+6,1:68,-1:imax+1):: u, um, up
    doubleprecision, dimension(-6:imax+6,1:68,-1:imax+1):: du, d2, vorp, vorm, vmpp, vmpm
    character                                   :: rho__
    character                                   :: not__
    character                                   :: dim_x
    character                                   :: dim_y
    character*5                                 :: varname
    doubleprecision                             :: test_rho_L, test_rho_R 
    doubleprecision                             :: DMM, DM4, dm4jph, dm4jmh, vul, vav, vmd, vlc, vmin, vmax
    integer                                     :: imp
    integer                                     :: coordenate

   
!!$   @function  get_hllc_vel

    DOUBLEPRECISION :: b, c, rt, v
    DOUBLEPRECISION :: fv, dfdv, oldv, err
    INTEGER         :: iter
    
    
    
  end module scalar

