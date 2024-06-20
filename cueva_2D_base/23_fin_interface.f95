
  subroutine fin_interface

    use parameters

    implicit none


    if (DIM == 1 .and. TEST == 1) then

       !  SCS TEST

       print*, '****************** FIN ******************** '
       print*,' YOU SOLVED Similar Current Sheet test 1D  '
       print*, '           ( PALENZUELA 2009 )              ' 
       print*,'                                             '
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*,'MALLA = ', imax
       print*,' CFL = ', CFL
       print*,' \sigma =  ', sigma
       if (MIRK ==2) print*,' c_1    = ', cm_1
       print*,'                                             ' 
       print*,'CAMPO MAGNETICO EN : test_SCS_IMEX_MIRK.dat  '
       print*, '' 
       print*,'                d(-_-)b                      '
       print*, '                                            '          
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, ''


    else if (DIM == 1 .and. TEST == 2 .or. TEST == 10) then

       ! CPAW TEST

       print*, '****************** FIN ******************** '
       print*,' YOU SOLVED THE  TEST Large Amplitude CP Alfven Waves 1D  '
       if (TEST == 2 ) print*, '           ( PALENZUELA 2009 )                         '
       if (TEST == 10) print*, '           ( TOMEK 2017 )                         ' 
       print*,'                                             '
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*,'MALLA = ', imax
       print*,' CFL = ', CFL
       print*,' \sigma =  ', sigma
       print*, 'ALFVEN VELOCITY', VA  
       if (MIRK ==2) print*,' c_1    = ', cm_1
       print*,'                                             ' 
       print*,'CAMPO MAGNETICO EN : test_CPAW_IMEX_MIRK.dat  '
       print*, '' 
       print*,'                d(-_-)b                      '
       print*, '                                            ' 
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '


    else if (DIM == 1 .and. TEST == 3) then

       ! RST_Palenzuela

       print*, '****************** FIN ******************** '
       print*,' YOU ARE SOLVED THE  TEST Resistive Shock Tube 1D  '
       print*, '           ( PALENZUELA 2009 )                         ' 
       print*,'                                             '
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*,'MALLA = ', imax
       print*,' CFL = ', CFL
       print*,' \sigma =  ', sigma
       if (MIRK ==2) print*,' c_1    = ', cm_1
       print*,'                                             ' 
       print*,'CAMPO MAGNETICO EN : test_RST_Palenzuela_IMEX_MIRK.dat  '
       print*, '' 
       print*,'                d(-_-)b                      '
       print*, '                                            '     
       print*, 'Extended GLM factor = ', FAC_EGLM        
       print*, '****************** FIN ******************** '
       print*, ''   


    else if (DIM ==1 .and. TEST == 4) then 

       !RST_BZ13 

       print*, '****************** FIN ******************** '
       print*,' YOU ARE SOLVED THE  TEST Resistive Shock Tube 1D  '
       print*, '           ( BZ 2013 )                         ' 
       print*,'                                             '
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*,'MALLA = ', imax
       print*,' CFL = ', CFL
       print*,' \sigma =  ', sigma
       if (MIRK ==2) print*,' c_1    = ', cm_1
       print*,'                                             ' 
       print*,'CAMPO MAGNETICO EN : test_RST_BZ13_IMEX_MIRK.dat  '
       print*, '' 
       print*,'                d(-_-)b                      '
       print*, '                                            '                
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, ''

    else if (DIM ==1 .and. TEST == 5) then 

       !RST_BZ13 

       print*, '****************** FIN ******************** '
       print*,' YOU ARE SOLVED THE  TEST Resistive Shock Tube 1D  '
       print*, '           ( Aloy & Miranda )                         ' 
       print*,'                                             '
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       if (MIRK == 2)  print*, 'MIRK parameter c1 = ',cm_1
       print*,'MALLA = ', imax
       print*,' CFL = ', CFL
       print*,' \sigma =  ', sigma
       print*, '' 
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*, ''           
       print*,' Shock Lorentz Factor Ws = ', Ws
       print*,' Shock Velocity Vs = ', Vs
       print*, '' 
       print*, 'preshock values' 
       print*, '' 
       print*, 'Vx_a  =' , Vx_a
       print*, 'Vy_a  =' , Vy_a
       print*, 'Vz_a  =' , Vz_a
       print*, ''
       print*, 'Bx_a  =' , Bx_a
       print*, 'By_a  =' , By_a
       print*, 'Bz_a  =' , Bz_a
       print*, ''
       print*, 'Ex_a  =' , Ex_a
       print*, 'Ey_a  =' , Ey_a
       print*, 'Ez_a  =' , Ez_a
       print*, ''
       print*, 'rho_a =' , rho_a
       print*, 'p_a   =' , p_a
       print*, 'q_a   =' , q_a
       print*, ''
       print*, 'Jx_a   =', Jx_a
       print*, 'Jy_a   =', Jy_a
       print*, 'Jz_a   =', Jz_a
       print*, '' 
       print*, '' 
       print*, 'postshock values' 
       print*, '' 
       print*, 'Vx_b  =' , Vx_b
       print*, 'Vy_b  =' , Vy_b
       print*, 'Vz_b  =' , Vz_b  
       print*, ''
       print*, 'Bx_b  =' , Bx_b
       print*, 'By_b  =' , By_b
       print*, 'Bz_b  =' , Bz_b
       print*, ''
       print*, 'Ex_b  =' , Ex_b
       print*, 'Ey_b  =' , Ey_b
       print*, 'Ez_b  =' , Ez_b
       print*, ''
       print*, 'rho_b =' , rho_b
       print*, 'p_b   =' , p_b
       print*, 'q_b   =' , q_b
       print*, '' 
       print*, 'Jx_b   =', Jx_b
       print*, 'Jy_b   =', Jy_b
       print*, 'Jz_b   =', Jz_b
       print*, '' 
       print*, '' 
       print*,'CAMPO MAGNETICO EN : test_RST_AlMi_IMEX_MIRK.dat  '
       print*, '' 
       print*,'                d(-_-)b                      '
       print*, '                                            '                
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, ''


    else if (DIM ==1 .and. TEST == 6) then

       ! Reconecction TEST 1D

       print*, '****************** FIN ******************** '
       print*,'               RECONNECTION 2D  '

       print*, '' 
        print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, 'Computational Box domain (Lx X Ly)', '(', Lx, 'X', Ly , ')'
       print*, 'This test made for cut in x= ', posy
       print*, '' 
       if (SLC == 1)  print*, 'spatially localized conductivity' 
       if (SLC == 0)  print*, 'constant conductivity' 
       print*, 'Conductivity = ',sigma 
       print*, 'Length of current sheet, Lr = ',Lr
       print*, 'Thickness of current sheet, delta_rec = ',delta_rec
       print*, 'Magnetization = ', sigma_m
       print*, 'ALFVEN VELOCITY = sqrt(sigma_m/(1 + sigma_m))  =', VA_rec  
       print*, 'Linquist number = ', S_l
       print*, 'Temperature (Theta)  = ', theta  
       print*, 'rho_up  = ', rho_up
       print*, 'B0 = sqrt(sigma_m * rho_up)  = ', B0 
       print*, 'P_up = theta * rho_up  = ', P_up
       print*, 'P0 = B0**2/2  = ', P0
       print*, 'Adiabatic exponent, gamma =', gamma
       print*, '' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, '' 
       if (MIRK == 2)  print*, 'MIRK parameter c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*,'Reconnection rate        : REC_rate_Ez_1D.dat        '
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '' 
       print*, '****************** FIN ******************** '
       print*, '' 

    else if (DIM ==1 .and. TEST == 7 .or. TEST == 8 .or. TEST == 9) then

       ! Tearing Mode TEST 1D

       print*, '****************** FIN ******************** '
       if (TEST == 7) print*,'               TEARING MODE       1D   '
       if (TEST == 8) print*,'               MAGNETIC DIFFUSION 1D   '
       if (TEST == 9) print*,'               SHEAR LAYER        1D   '
       print*,'                 Tomek Setup                 '
       print*, '' 
        print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, 'Computational Box domain (Longx)', '(', Longx, ')'
       print*, '' 
       if (SLC == 1)  print*, 'spatially localized conductivity' 
       if (SLC == 0)  print*, 'constant conductivity' 
       print*, 'Conductivity = ',sigma 
       if (TEST ==7)  print*, 'Thickness of current sheet, delta_rec = ',delta_rec
       if (TEST ==8)  print*, 'Thickness of current sheet, delta_tm  = ',delta_tm
       print*, 'rho_up  = ', rho_up
       print*, 'B0 =  ', B0 
       print*, 'P_up =', P_up
       print*, 'ALFVEN VELOCITY = B0/sqrt(rho_up)  =', VA_rec_tearing  
       print*, 'Linquist number = Sl = sigma * Lx * V_A =', sigma * Longx * VA_rec_tearing
       print*, 'Adiabatic exponent, gamma =', gamma
       print*, '' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, '' 
       if (MIRK == 2)  print*, 'MIRK parameter c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*,'Reconnection rate        : REC_rate_Ez_1D.dat        '
       print*,'Temporal variation of \delta : delta_rec_t_1D.dat        '
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '' 
       print*, '****************** FIN ******************** '
       print*, '' 

    else if (DIM == 2 .and. TEST == 20) then

       ! CE TEST

       print*, '****************** FIN ******************** '
       print*,' '
       print*,' YOU ARE SOLVED CYLINDRICAL EXPLOSION 2D TEST  '
       print*,'             ( PALENZUELA 2009 ) '
       print*, '' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',  jmax
       print*, 'CFL    = ', CFL 
       print*, '\sigma = ',sigma 
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*,'CAMPO MAGNETICO EN : CE_IMEX_MIRK_Bx.dat'
       print*,'                   : CE_IMEX_MIRK_By.dat'
       print*,'PRESION            : CE_IMEX_MIRK_P.dat'
       print*,'FACTOR DE LORENTZ  : CE_IMEX_MIRK_W.dat'
       print*,'COMPARE Bx         : CE_compare_Bx.dat        '
       print*,'COMPARE P          : CE_compare_P.dat        '
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '                                            ' 
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, '' 

    else if (DIM == 2 .and. TEST == 21) then

       ! Rotor TEST 

       print*, '****************** FIN ******************** '
       print*,' YOU ARE SOLVED RESISTIVE ROTOR 2D TEST    '
       print*, '    Bucciantini & DelZanna  (2013)         ' 
       print*, '' 
       print*, '' 
       print*, ' ' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, '\sigma = ',sigma 
       print*, '\Omega = ',omega 
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, ''
       print*,'CAMPO ELECTRICO EN : Rotor_IMEX_MIRK_Ez.dat'
       print*,'DENSIDAD           : Rotor_IMEX_MIRK_rho.dat'
       print*,'PRESION            : Rotor_IMEX_MIRK_p.dat'
       print*,'FACTOR DE LORENTZ  : Rotor_IMEX_MIRK_W.dat'
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '                                            ' 
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, '' 

    else if (DIM == 2 .and. TEST == 22) then

       ! CS TEST 

       print*, '****************** FIN ******************** '
       print*,' YOU ARE SOLVED Cylindrical Star 2D TEST    '
       print*, '          ( PALENZUELA 2009 )          ' 
       print*, '' 
       print*, '' 
       print*, ' ' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, '\sigma = ',sigma 
       print*, '\Omega = ', omega 
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*,'CAMPO MAGNETICO EN : CS_IMEX_MIRK_Bz.dat'
       print*,'CAMPO ELECTRICO EN : CS_IMEX_MIRK_Bz.dat'
       print*,'DENSIDAD           : CS_IMEX_MIRK_rho.dat'
       print*,'PRESION            : CS_IMEX_MIRK_p.dat'
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '                                            ' 
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, '' 

    else if (DIM == 2 .and. TEST == 23) then

       ! SJ TEST 

       print*, '****************** FIN ******************** '
       print*,' YOU ARE SOLVED Slab Jet 2D TEST    '
       print*, '          ( Komissarov 1999 )          ' 
       print*, '' 
       print*, '' 
       print*, ' ' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, '\sigma = ',sigma 
       print*, '\Omega = ', omega 
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*,'CAMPO MAGNETICO EN : SJ_IMEX_MIRK_rho.dat'
       print*,'CAMPO ELECTRICO EN : SJ_IMEX_MIRK_Vx.dat'
       print*,'DENSIDAD           : SJ_IMEX_MIRK_P.dat'
       print*,'PRESION            : SJ_IMEX_MIRK_Pmag.dat'
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '                                            ' 
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, '' 

    else if (DIM == 2 .and. TEST == 24) then

       ! CSI TEST 

       print*, '****************** FIN ******************** '
       print*,'  YOU ARE SOLVED Cloud Shock Interaction Test    '
       print*, '   ( Carlos M. Xisto, et. al.  JCP (2014)  ) ' 
       print*, '' 
       print*, '' 
       print*, ' ' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, '\sigma = ',sigma 
       print*, '\Omega = ', omega 
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*,'DENSIDAD        : CSI_IMEX_MIRK_rho.dat'
       print*,'VELOCIDAD Vx    : CSI_IMEX_MIRK_Vx.dat'
       print*,'PRESION         : CSI_IMEX_MIRK_P.dat'
       print*,'PRESION MAG     : CSI_IMEX_MIRK_Pmag.dat'
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '                                            ' 
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, '' 

    else if (DIM == 2 .and. TEST == 25) then

       ! CSI TEST 

       print*, '****************** FIN ******************** '
       print*,'             YOU ARE SOLVED Richtmyer-Meshkov instability           '
       print*,'( O. Zanotti & M. Dumbser, Computer Physics Communications (2015)  ) ' 
       print*, '' 
       print*, '' 
       print*, ' ' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, '\sigma = ',sigma 
       print*, '\Omega = ', omega 
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*, '' 
       print*, 'Ws     = ', Ws
       print*, 'Vs     = ', Vs
       print*, 'p_b    = ', p_b
       print*, 'rho_b  = ', rho_b
       print*, '' 
       print*,'DENSIDAD  t =100%  : RMI_IMEX_MIRK_rho_t=100.dat'
       print*,'DENSIDAD  t =075%  : RMI_IMEX_MIRK_rho_t=075.dat'
       print*,'DENSIDAD  t =050%  : RMI_IMEX_MIRK_rho_t=050.dat'
       print*,'DENSIDAD  t =025%  : RMI_IMEX_MIRK_rho_t=025.dat'
       print*,'DENSIDAD  t =000%  : RMI_IMEX_MIRK_rho_t=000.dat'
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '                                            ' 
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, ''


    else if (DIM == 2 .and. TEST == 26 .or. TEST == 27) then

       ! CSI TEST 

       print*, '****************** FIN ******************** '
       if (TEST == 26) print*,'             YOU ARE SOLVED TEST 2D Magnetic Diffusion    '
       if (TEST == 27) print*,'             YOU ARE SOLVED TEST 2D Shear Layer           '
       print*,'                ( Tomek setup (2017)  )                ' 
       print*, '' 
       print*, '' 
       print*, ' ' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       if (flux_solver == 5) print*, 'HLL + MP5 '
       print*, '' 
       print*, '' 
       print*, '\sigma = ',sigma 
       if (TEST == 26)  print*, ' B0    = ', B0
       if (TEST == 27)  print*, ' V0    = ', V0
       print*, 'rho    = ', rho_up
       print*, 'p      = ', P_up
       print*, '' 
       if (TEST == 26)  print*,'t =100%  : magnetic_diffusion/MD_2D_Vx_t_100.dat'
       if (TEST == 27)  print*,'t =100%  : shear_layer/SL_2D_Vx_t_100.dat'
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '                                            ' 
       print*, 'Extended GLM factor = ', FAC_EGLM            
       print*, '****************** FIN ******************** '
       print*, ''

       

    else if (DIM ==2 .and. TEST == 30 .or. TEST ==31) then

       ! Reconecction TEST

       print*, '****************** FIN ******************** '
       print*,'               RECONNECTION 2D  '
       if (TEST ==30) print*, '                M.A. Aloy & Petar Mimica      ' 
       if (TEST ==31) print*, '                Zenitani et al. (2010)      ' 
       print*, '' 
       print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, 'Computational Box domain (Lx X Ly)', '(', Lx, 'X', Ly , ')'
       print*, '' 
       if (SLC == 1)  print*, 'spatially localized conductivity' 
       if (SLC == 1)  print*, 'sigma_0 = ',sigma_0
       if (SLC == 1)  print*, 'sigma_1 = ',sigma_1
       if (SLC == 0)  print*, 'constant conductivity' 
       if (SLC == 0)  print*, 'Conductivity = ',sigma 
       print*, 'Length of current sheet, Lr = ',Lr
       print*, 'Thickness of current sheet, delta_rec = ',delta_rec
       print*, 'Magnetization = ', sigma_m
       print*, 'ALFVEN VELOCITY = sqrt(sigma_m/(1 + sigma_m))  =', VA_rec  
       print*, 'Linquist number = ', S_l
       print*, 'Temperature (Theta)  = ', theta  
       print*, 'rho_up  = ', rho_up
       print*, 'B0 = sqrt(sigma_m * rho_up)  = ', B0 
       print*, 'P_up = theta * rho_up  = ', P_up
       print*, 'P0 = B0**2/2  = ', P0
       print*, 'Adiabatic exponent, gamma =', gamma
       print*, '' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, '' 
       if (MIRK == 2)  print*, 'MIRK parameter c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*,'CAMPO MAGNETICO EN : reconnectionBx.dat'
       print*,'                   : reconnectionBy.dat'
       print*,'PRESION            : reconnectionP.dat'
       print*,'FACTOR DE LORENTZ  : reconnectionW.dat'
       print*,'COMPARE Bx         : REC_compare_Bx.dat        '
       print*,'COMPARE P          : REC_compare_P.dat        '
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '' 
       print*, '****************** FIN ******************** '
       print*, '' 

    else if (DIM ==2 .and. TEST == 31) then

       ! Reconecction TEST

       print*, '****************** FIN ******************** '
       print*,'               RECONNECTION 2D  '

       print*, '' 
        print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, 'Computational Box domain (Lx X Ly)', '(', Lx, 'X', Ly , ')'
       print*, '' 
       if (SLC == 1)  print*, 'spatially localized conductivity' 
       if (SLC == 0)  print*, 'constant conductivity' 
       print*, 'Conductivity = ',sigma 
       print*, 'Length of current sheet, Lr = ',Lr
       print*, 'Thickness of current sheet, delta_rec = ',delta_rec
       print*, 'Magnetization = ', sigma_m
       print*, 'ALFVEN VELOCITY = sqrt(sigma_m/(1 + sigma_m))  =', VA_rec  
       print*, 'Linquist number = ', S_l
       print*, 'Temperature (Theta)  = ', theta  
       print*, 'rho_up  = ', rho_up
       print*, 'B0 = sqrt(sigma_m * rho_up)  = ', B0 
       print*, 'P_up = theta * rho_up  = ', P_up
       print*, 'P0 = B0**2/2  = ', P0
       print*, 'Adiabatic exponent, gamma =', gamma
       print*, '' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, '' 
       if (MIRK == 2)  print*, 'MIRK parameter c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       print*, '' 
       print*,'CAMPO MAGNETICO EN : reconnectionBx.dat'
       print*,'                   : reconnectionBy.dat'
       print*,'PRESION            : reconnectionP.dat'
       print*,'FACTOR DE LORENTZ  : reconnectionW.dat'
       print*,'COMPARE Bx         : REC_compare_Bx.dat        '
       print*,'COMPARE P          : REC_compare_P.dat        '
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '' 
       print*, '****************** FIN ******************** '
       print*, '' 


    else if (DIM ==2 .and. TEST == 32 .or. TEST == 33 .or. TEST == 34 .or. TEST == 36) then

       ! Reconecction TEST

       print*, '****************** FIN ******************** '
       print*,'               Tearing Mode  2D  '
       if (TEST == 32) print*, '                            ( Del Zanna et al.)            '
       if (TEST == 33) print*, '                            ( Hubert Baty  RDTM (2017))    '
       print*, '' 
       print*, 'MALLA  = ',imax,'   X',jmax
       print*, 'CFL    = ', CFL 
       print*, 'Computational Box domain (Lx X Ly)', '(', Lx, 'X', Ly , ')'
       print*, '' 
       if (SLC == 1)  print*, 'spatially localized conductivity' 
       if (SLC == 0)  print*, 'constant conductivity' 
       print*, 'Conductivity = ',sigma 
       print*, 'Linquist number Sl= ', S_l
       if (TEST == 33) print*, 'Linquist number Sa= ', S_a
       if (TEST == 33) print*, 'Length of current sheet, Lr = ',Lr
       print*, 'Thickness of current sheet, a_tm = ', a_tm
       print*, 'separation current sheest, l_tm = ', l_tm
       print*, 'Maximun wavenumber,  k_tm = ', k_tm
       print*, 'Magnetization, sigma_m = ', sigma_m
       print*, 'beta_m    = 2.d0 * P0 / B0**2 = ', beta_m
       print*, 'ALFVEN VELOCITY = 1.d0 / sqrt(1.d0/sigma_m + 2.d0 * beta_m + 1.d0) =', VA_rec_tearing
       print*, 'rho_up  = ', rho_up
       print*, 'B0 = sqrt(sigma_m * rho_up)  = ', B0 
       print*, 'P_up = theta * rho_up  = ', P_up
       print*, 'P0 = B0**2/2  = ', P0
       print*, 'Adiabatic exponent, gamma =', gamma
       print*, '' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, '' 
       if (MIRK == 2)  print*, 'MIRK parameter c1 = ',cm_1
       if (flux_solver == 1) print*, 'LFF flux solver '
       if (flux_solver == 2) print*, 'HLL flux solver '
       if (flux_solver == 5) print*, 'HLL flux solver + MP5'
       if (flux_solver == 7) print*, 'HLLC flux solver + MP5'
       print*, '' 
       print*,'CAMPO MAGNETICO EN : reconnection/TM_2D_Bx_t=100.dat'
       print*,'                   : reconnection/TM_2D_By_t=100.dat'
       print*, '' 
       print*,'                d(-_-)b          '
       print*, '' 
       print*, '****************** FIN ******************** '
       print*, '' 


    else 

       write(*,*) "STOP subroutine fin_interface"
       write(*,*) "This TEST is not implemented yet"
       stop

    end if

  end subroutine fin_interface

