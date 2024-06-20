
  !     <<<<<<<<<<<<<<<<<<<<<<<
  !     Interface
  !     <<<<<<<<<<<<<<<<<<<<<<<

  subroutine interface

    use parameters

    implicit none



    if (DIM == 1 .and. TEST == 1) then

       !  SCS TEST

       print*,  ' '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' ESTA RESOLVIENDO EL TEST SIMILAR CURRENT SHEET PROBLEM 1D '
       print*,  '              ( Palenzuela et al. (2009) )                ' 
       print*, ' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '



    else if (DIM == 1 .and. TEST == 2 .or. TEST == 10) then

       ! CPAW TEST

       print*,  ' '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' ESTA RESOLVIENDO EL TEST Large Amplitude CP Alfvèn Wave 1D  '
       if (TEST == 2 ) print*,  '              ( Palenzuela et al. (2009) )                '
       if (TEST == 10) print*,  '              ( Tomek (2017) )                ' 
       print*, ' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '


    else if (DIM == 1 .and. TEST == 3) then

       ! RST_Palenzuela

       print*,  ' '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' ESTA RESOLVIENDO EL TEST  RESISTIVE SHOCK TUBE  1D '
       print*,  '              ( Palenzuela et al. (2009) )                ' 
       print*, ' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '


    else if (DIM ==1 .and. TEST == 4) then 

       !RST_BZ13 

       print*,  ' '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' ESTA RESOLVIENDO EL TEST RESISTIVE SHOCK TUBE PROBLEM 1D  '
       print*,  '              ( Bucciantini & Del Zanna (2013) )                ' 
       print*, ' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '

    else if (DIM ==1 .and. TEST == 5) then 

       !RST_BZ13 

       print*,  ' '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' ESTA RESOLVIENDO EL TEST RESISTIVE SHOCK TUBE PROBLEM 1D  '
       print*, ' PARA COMPROBAR LAS CONDICIONES DE RANKINE-HUGONIOT DEDUCIDAS '
       print*,  '               ( Aloy & Miranda )                ' 
       print*, ' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' Shock Velocity' 
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
       print*, ' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '


    else if (DIM ==1 .and. TEST == 6) then

       ! Reconecction TEST 1D

       print*, ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '               YOU ARE SOLVING THE RECONNECTION 1D PROBLEM               '
       if (TEST == 6) print*, '                            ( M.A. Aloy & Petar Mimica)   ' 
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' PARAMETERS'
       if (SLC == 1)  print*, 'spatially localized conductivity' 
       if (SLC == 1)  print*, 'sigma_0 = ',sigma_0
       if (SLC == 1)  print*, 'sigma_1 = ',sigma_1
       if (SLC == 0)  print*, 'constant conductivity' 
       if (SLC == 0)  print*, 'Conductivity = ',sigma 
       print*, ' '
       print*, ' Alven Velocity', VA_rec
       print*, ' P_up', P_up
       print*, ' \delta ', delta_rec
       print*, ' Ly / \delta it must be great than 10.d0', Ly / delta_rec
       print*, ' '
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,'' 

    else if (DIM == 1 .and. TEST == 7 .or. TEST == 8 .or. TEST == 9) then

       ! Reconecction TEST 1D

       print*, ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       if (TEST == 7) print*, '               YOU ARE SOLVING THE TEARING MODE 1D PROBLEM        '
       if (TEST == 8) print*, '               YOU ARE SOLVING THE MAGNETIC DIFFUSION 1D PROBLEM  '
       if (TEST == 9) print*, '               YOU ARE SOLVING THE SHEAR LAYER 1D PROBLEM         '
       print*, '                            (Tomasz Rembiasz)   ' 
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' PARAMETERS'
       if (SLC == 1)  print*, 'spatially localized conductivity' 
       if (SLC == 1)  print*, 'sigma_0 = ',sigma_0
       if (SLC == 1)  print*, 'sigma_1 = ',sigma_1
       if (SLC == 0)  print*, 'constant conductivity' 
       if (SLC == 0)  print*, 'Conductivity = ',sigma 
       print*, ' '
       print*, ' Alven Velocity', VA_rec_tearing
       print*, ' P_up', P_up
       print*, ' \delta ', delta_rec
       print*, ' '
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,'' 




    else if (DIM == 2 .and. TEST == 20) then

       ! CE TEST

       print*, ' '
       print*, '       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,                      'CYLINDRICAL EXPLOSION 2D TEST '
       print*, '                      ( Palenzuela et al. (2009) )                          ' 
       print*,'' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',  jmax
       print*, '\sigma = ',sigma 
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       print*, ' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,'' 




    else if (DIM == 2 .and. TEST == 21) then

       ! Rotor TEST 

       print*,  ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '                   YOU ARE SOLVING THE TEST ROTOR 2D                       '
       print*,  '                  ( Bucciantini & Del Zanna (2013) )                ' 
       print*,'' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',  jmax
       print*, '\omega = ',omega
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       print*, ' \sigma ', sigma
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '

    else if (DIM == 2 .and. TEST == 22) then

       ! Rotor TEST 

       print*,  ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '                   YOU ARE SOLVING THE TEST Cylindrical Star 2D                       '
       print*,  '                  ( Palenzuela et al. (2009) )                 ' 
       print*,'' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',  jmax
       print*, '\omega = ', omega
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       print*, ' \sigma ', sigma
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '

    else if (DIM == 2 .and. TEST == 23) then

       ! Rotor TEST 

       print*,  ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '                   YOU ARE SOLVING THE TEST Slab Jet 2D                       '
       print*,  '                  ( Komissarov. (1999) )                 ' 
       print*,'' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',  jmax
       print*, '\omega = ', omega
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       print*, ' \sigma ', sigma
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '

    else if (DIM == 2 .and. TEST == 24) then

       ! Rotor TEST 

       print*,  ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '                   YOU ARE SOLVING THE TEST Cloud Shock Interaction        '
       print*,  '                  ( Carlos M. Xisto, et. al.  JCP (2014)  )                 ' 
       print*,'' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',  jmax
       print*, '\omega = ', omega
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       print*, ' \sigma ', sigma
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '

    else if (DIM == 2 .and. TEST == 25) then

       ! Rotor TEST 

       print*,  ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '                YOU ARE SOLVING THE TEST Richtmyer-Meshkov instability        '
       print*,  ' ( Olindo Zanotti & Michael Dumbser,  Computer Physics Communications (2015)  ) ' 
       print*,  ' ' 
       print*,'' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',  jmax
       print*, '\omega = ', omega
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       print*, ' \sigma ', sigma
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '

       
    else if (DIM == 2 .and. TEST == 26 .or. TEST == 27) then

       ! Rotor TEST 

       print*,  ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       if (TEST == 26) print*, '                YOU ARE SOLVING THE TEST 2D Magnetic Diffusion  '
       if (TEST == 27) print*, '                YOU ARE SOLVING THE TEST 2D Shear Layer         '
       print*,  '                      Tomek setup (2016)                   ' 
       print*,  ' ' 
       print*,'' 
       print*, 'ORDER ---> MIRK  = ',MIRK,'   IMEX',  order
       print*, 'MALLA  = ',imax,'   X',  jmax
       if (MIRK == 2)  print*, 'c1 = ',cm_1
       print*, ' \sigma ', sigma
       if (TEST == 26) print*, 'B0 = ', B0
       if (TEST == 27) print*, 'V0 = ', V0
       print*,  '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '


    else if (DIM ==2 .and. TEST == 30 .or. TEST == 31 .or. TEST == 32 .or. TEST == 33 .or. TEST == 34 .or. TEST == 36) then

       ! Reconecction TEST

       print*, ' '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '               YOU ARE SOLVING THE RECONNECTION 2D PROBLEM               '
       if (TEST == 30) print*, '                            ( M.A. Aloy & Petar Mimica)   ' 
       if (TEST == 31) print*, '                            ( Zenitani et al.)            ' 
       if (TEST == 32) print*, '                            ( Del Zanna et al.)            '
       if (TEST == 33) print*, '                            ( Hubert Baty  RDTM (2017))    ' 
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, ' PARAMETERS'
       if (SLC == 1)  print*, 'spatially localized conductivity' 
       if (SLC == 1)  print*, 'sigma_0 = ',sigma_0
       if (SLC == 1)  print*, 'sigma_1 = ',sigma_1
       if (SLC == 0)  print*, 'constant conductivity' 
       if (SLC == 0)  print*, 'Conductivity = ',sigma 
       print*, ' '
       if (TEST == 32 .or. TEST == 33) print*,'a_tm   = ',a_tm
       print*,'Lx     = ', Lx
       print*,'Ly     = ', Ly
       if (TEST == 32 .or. TEST == 33) print*,'k_tm   = ', k_tm
       if (TEST == 32 .or. TEST == 33) print*,'kap_tm = ', kap_tm
       print*, ' '
       if (TEST == 32 .or. TEST == 33) print*, ' Alven Velocity', VA_rec_tearing
       if (TEST == 30) print*,' Alven Velocity', VA_rec
       print*, ' P_up', P_up
       print*, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*, '     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
       print*,'' 

    else 

       write(*,*) "STOP sunroutine interface"
       write(*,*) "This TEST is not implemented yet"
       stop

    end if


  end subroutine interface
