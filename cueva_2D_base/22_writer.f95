
  !     *******************************************************************
  !     Subroutine that writes the final results
  !     *******************************************************************


  subroutine lector

    use scalar
    use threevectors

    implicit none

!------------ Reconecction TEST ------------

!      integer ii
      integer ios

      integer, parameter:: outFileUnit   = 1
      integer, parameter:: outFileUnit1  = 2
      integer, parameter:: outFileUnit2  = 3
      integer, parameter:: outFileUnit3  = 4
      integer, parameter:: outFileUnit4  = 5
      integer, parameter:: outFileUnit5  = 7
      integer, parameter:: outFileUnit6  = 8
      integer, parameter:: outFileUnit7  = 9
      integer, parameter:: outFileUnit8  = 10
      integer, parameter:: outFileUnit9  = 11
      integer, parameter:: outFileUnit10 = 12
      integer, parameter:: outFileUnit11 = 14
      integer, parameter:: outFileUnit12 = 15
      integer, parameter:: outFileUnit13 = 16
      integer, parameter:: outFileUnit14 = 17
      integer, parameter:: outFileUnit15 = 18
      integer, parameter:: outFileUnit16 = 19
            ! for back-up
      integer, parameter:: outFileUnit17 = 20
      integer, parameter:: outFileUnit18 = 21

!      Buffers to hold file names
      character(80) outputFileName
      character(80) outputFileName1
      character(80) outputFileName2
      character(80) outputFileName3
      character(80) outputFileName4
      character(80) outputFileName5
      character(80) outputFileName6
      character(80) outputFileName7
      character(80) outputFileName8      
      character(80) outputFileName9
      character(80) outputFileName10
      character(80) outputFileName11
      character(80) outputFileName12
      character(80) outputFileName13
      character(80) outputFileName14
      character(80) outputFileName15
      character(80) outputFileName16
            !for back-up
      character(80) outputFileName17
      character(80) outputFileName18

      CHARACTER(*), PARAMETER :: fileplace = "/data_local/"
 


    if (DIM == 1) then


       if ( mod(h,hloc) == 0 ) then   !( mod( floor(t_advance / Delt_loc),floor(Delt_loc) ) == 0 ) then   !

          hrec = h / hloc !ceiling(t_advance) !floor(t_advance / Delt_loc) !


       if ( hrec .lt. 10) then
          write (outputFileName , '(a, I0, a)') 'data_local/1d_data_psi_000' ,hrec,'.dat'  
          write (outputFileName1, '(a, I0, a)') 'data_local/1d_data_phi_000' ,hrec,'.dat'  
          write (outputFileName2, '(a, I0, a)') 'data_local/1d_data_vx_000'  ,hrec,'.dat'  
          write (outputFileName3, '(a, I0, a)') 'data_local/1d_data_vy_000'  ,hrec,'.dat'  
          write (outputFileName4, '(a, I0, a)') 'data_local/1d_data_vz_000'  ,hrec,'.dat'
          write (outputFileName5, '(a, I0, a)') 'data_local/1d_data_bx_000'  ,hrec,'.dat'
          write (outputFileName6, '(a, I0, a)') 'data_local/1d_data_by_000'  ,hrec,'.dat'
          write (outputFileName7, '(a, I0, a)') 'data_local/1d_data_bz_000'  ,hrec,'.dat'
          write (outputFileName8, '(a, I0, a)') 'data_local/1d_data_ex_000'  ,hrec,'.dat'
          write (outputFileName9, '(a, I0, a)') 'data_local/1d_data_ey_000'  ,hrec,'.dat'
          write (outputFileName10,'(a, I0, a)') 'data_local/1d_data_ez_000'  ,hrec,'.dat'
          write (outputFileName11,'(a, I0, a)') 'data_local/1d_data_p_000'   ,hrec,'.dat'
          write (outputFileName12,'(a, I0, a)') 'data_local/1d_data_q_000'   ,hrec,'.dat'
          write (outputFileName13,'(a, I0, a)') 'data_local/1d_data_rho_000' ,hrec,'.dat'
       else if  ( hrec .lt. 100) then
          write (outputFileName , '(a, I0, a)') 'data_local/1d_data_psi_00'  ,hrec,'.dat'  
          write (outputFileName1, '(a, I0, a)') 'data_local/1d_data_phi_00'  ,hrec,'.dat'  
          write (outputFileName2, '(a, I0, a)') 'data_local/1d_data_vx_00'   ,hrec,'.dat'  
          write (outputFileName3, '(a, I0, a)') 'data_local/1d_data_vy_00'   ,hrec,'.dat'  
          write (outputFileName4, '(a, I0, a)') 'data_local/1d_data_vz_00'   ,hrec,'.dat'
          write (outputFileName5, '(a, I0, a)') 'data_local/1d_data_bx_00'   ,hrec,'.dat'
          write (outputFileName6, '(a, I0, a)') 'data_local/1d_data_by_00'   ,hrec,'.dat'
          write (outputFileName7, '(a, I0, a)') 'data_local/1d_data_bz_00'   ,hrec,'.dat'
          write (outputFileName8, '(a, I0, a)') 'data_local/1d_data_ex_00'   ,hrec,'.dat'
          write (outputFileName9, '(a, I0, a)') 'data_local/1d_data_ey_00'   ,hrec,'.dat'
          write (outputFileName10,'(a, I0, a)') 'data_local/1d_data_ez_00'   ,hrec,'.dat'
          write (outputFileName11,'(a, I0, a)') 'data_local/1d_data_p_00'    ,hrec,'.dat'
          write (outputFileName12,'(a, I0, a)') 'data_local/1d_data_q_00'    ,hrec,'.dat'
          write (outputFileName13,'(a, I0, a)') 'data_local/1d_data_rho_00'  ,hrec,'.dat'
       else if (hrec .lt. 1000) then
          write (outputFileName , '(a, I0, a)') 'data_local/1d_data_psi_0'  ,hrec,'.dat'  
          write (outputFileName1, '(a, I0, a)') 'data_local/1d_data_phi_0'  ,hrec,'.dat'  
          write (outputFileName2, '(a, I0, a)') 'data_local/1d_data_vx_0'   ,hrec,'.dat'  
          write (outputFileName3, '(a, I0, a)') 'data_local/1d_data_vy_0'   ,hrec,'.dat'  
          write (outputFileName4, '(a, I0, a)') 'data_local/1d_data_vz_0'   ,hrec,'.dat'
          write (outputFileName5, '(a, I0, a)') 'data_local/1d_data_bx_0'   ,hrec,'.dat'
          write (outputFileName6, '(a, I0, a)') 'data_local/1d_data_by_0'   ,hrec,'.dat'
          write (outputFileName7, '(a, I0, a)') 'data_local/1d_data_bz_0'   ,hrec,'.dat'
          write (outputFileName8, '(a, I0, a)') 'data_local/1d_data_ex_0'   ,hrec,'.dat'
          write (outputFileName9, '(a, I0, a)') 'data_local/1d_data_ey_0'   ,hrec,'.dat'
          write (outputFileName10,'(a, I0, a)') 'data_local/1d_data_ez_0'   ,hrec,'.dat'
          write (outputFileName11,'(a, I0, a)') 'data_local/1d_data_p_0'    ,hrec,'.dat'
          write (outputFileName12,'(a, I0, a)') 'data_local/1d_data_q_0'    ,hrec,'.dat'
          write (outputFileName13,'(a, I0, a)') 'data_local/1d_data_rho_0'  ,hrec,'.dat'
       else
          write (outputFileName , '(a, I0, a)') 'data_local/1d_data_psi_'  ,hrec,'.dat'  
          write (outputFileName1, '(a, I0, a)') 'data_local/1d_data_phi_'  ,hrec,'.dat'  
          write (outputFileName2, '(a, I0, a)') 'data_local/1d_data_vx_'   ,hrec,'.dat'  
          write (outputFileName3, '(a, I0, a)') 'data_local/1d_data_vy_'   ,hrec,'.dat'  
          write (outputFileName4, '(a, I0, a)') 'data_local/1d_data_vz_'   ,hrec,'.dat'
          write (outputFileName5, '(a, I0, a)') 'data_local/1d_data_bx_'   ,hrec,'.dat'
          write (outputFileName6, '(a, I0, a)') 'data_local/1d_data_by_'   ,hrec,'.dat'
          write (outputFileName7, '(a, I0, a)') 'data_local/1d_data_bz_'   ,hrec,'.dat'
          write (outputFileName8, '(a, I0, a)') 'data_local/1d_data_ex_'   ,hrec,'.dat'
          write (outputFileName9, '(a, I0, a)') 'data_local/1d_data_ey_'   ,hrec,'.dat'
          write (outputFileName10,'(a, I0, a)') 'data_local/1d_data_ez_'   ,hrec,'.dat'
          write (outputFileName11,'(a, I0, a)') 'data_local/1d_data_p_'    ,hrec,'.dat'
          write (outputFileName12,'(a, I0, a)') 'data_local/1d_data_q_'    ,hrec,'.dat'
          write (outputFileName13,'(a, I0, a)') 'data_local/1d_data_rho_'  ,hrec,'.dat'
       end if



        !
        ! It will bail out rather than overwrite an existing file.
        !
        open(outFileUnit, file=outputFileName,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName)
           stop
        endif
        open(outFileUnit1, file=outputFileName1,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName1)
           stop
        endif
        open(outFileUnit2, file=outputFileName2,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName2)
           stop
        endif
        write(*, '("Opened file ", a, " for writing.")') &
                  trim(outputFileName2)
        open(outFileUnit3, file=outputFileName3,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName3)
           stop
        endif
        open(outFileUnit4, file=outputFileName4,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName4)
           stop
        endif
        open(outFileUnit5, file=outputFileName5,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName5)
           stop
        endif
        open(outFileUnit6, file=outputFileName6,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName6)
           stop
        endif
        open(outFileUnit7, file=outputFileName7,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName7)
           stop
        endif
        open(outFileUnit8, file=outputFileName8,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName8)
           stop
        endif
        open(outFileUnit9, file=outputFileName9,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName9)
           stop
        endif
        open(outFileUnit10, file=outputFileName10,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName10)
           stop
        endif
        open(outFileUnit11, file=outputFileName11,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName11)
           stop
        endif
        open(outFileUnit12, file=outputFileName12,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName12)
           stop
        endif
        open(outFileUnit13, file=outputFileName13,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName13)
           stop
        endif

        ! Here's where you read stuff from the input file
        ! and write stuff to the output file.

           
        do i=0,imax

           write (outFileUnit ,*)  - 0.5d0 * Longx + i * Delx,  psi(i,1,1)   
           write (outFileUnit1,*)  - 0.5d0 * Longx + i * Delx,  phi(i,1,1)   
           write (outFileUnit2,*)  - 0.5d0 * Longx + i * Delx,  Vx (i,1,1)   
           write (outFileUnit3,*)  - 0.5d0 * Longx + i * Delx,  Vy (i,1,1)   
           write (outFileUnit4,*)  - 0.5d0 * Longx + i * Delx,  Vz (i,1,1)
           write (outFileUnit5,*)  - 0.5d0 * Longx + i * Delx,  Bx (i,1,1)
           write (outFileUnit6,*)  - 0.5d0 * Longx + i * Delx,  By (i,1,1)
           write (outFileUnit7,*)  - 0.5d0 * Longx + i * Delx,  Bz (i,1,1)
           write (outFileUnit8,*)  - 0.5d0 * Longx + i * Delx,  Ex (i,1,1)
           write (outFileUnit9,*)  - 0.5d0 * Longx + i * Delx,  Ey (i,1,1)
           write (outFileUnit10,*) - 0.5d0 * Longx + i * Delx,  Ez (i,1,1)
           write (outFileUnit11,*) - 0.5d0 * Longx + i * Delx,  p  (i,1,1)
           write (outFileUnit12,*) - 0.5d0 * Longx + i * Delx,  q  (i,1,1)
           write (outFileUnit13,*) - 0.5d0 * Longx + i * Delx,  rho(i,1,1)

        end do ! for i

     ! Close files before next pass through the loop
     
        write(*, '("Closing files.")')
        close(outFileUnit)
        close(outFileUnit1)
        close(outFileUnit2)
        close(outFileUnit3)
        close(outFileUnit4)
        close(outFileUnit5)
        close(outFileUnit6)
        close(outFileUnit7)
        close(outFileUnit8)
        close(outFileUnit9)
        close(outFileUnit10)
        close(outFileUnit11)
        close(outFileUnit12)
        close(outFileUnit13)
        
        write(*,*) ""

     end if


!-------------------------------------------------------------
! Evolution of different quatities for 1D REC, TM, MD and SL
!-------------------------------------------------------------


     if ( h == 1 .or. mod(h,hglb) == 0 ) then

     sum_psi  = 0.d0
     sum_phi  = 0.d0
     sum_ex   = 0.d0
     sum_ey   = 0.d0
     sum_ez   = 0.d0
     sum_exey = 0.d0
     sum_exez = 0.d0
     sum_eyez = 0.d0
     sum_bx   = 0.d0
     sum_by   = 0.d0
     sum_bz   = 0.d0
     sum_bxby = 0.d0
     sum_bxbz = 0.d0
     sum_bybz = 0.d0
     sum_vx   = 0.d0
     sum_vy   = 0.d0
     sum_vz   = 0.d0
     sum_vxvy = 0.d0
     sum_vxvz = 0.d0
     sum_vyvz = 0.d0
     sum_q    = 0.d0
     sum_p    = 0.d0
     sum_rho  = 0.d0
     sum_bxbybz = 0.d0
     sum_vxvyvz = 0.d0
     ekin_c     = 0.d0
     ekin_r     = 0.d0
     eint_c     = 0.d0
     eint_r     = 0.d0
     etot       = 0.d0
     efluid_r   = 0.d0
     emag       = 0.d0
     lorentz    = 0.d0
     sum_D      = 0.d0
     sum_D_med  = 0.d0

     do i=0,imax

              sum_psi  = sum_psi  + psi(i,1,1)            * Delx 
              sum_phi  = sum_phi  + phi(i,1,1)            * Delx 
              
              sum_ex   = sum_ex   + Ex(i,1,1)**2          * Delx 
              sum_ey   = sum_ey   + Ey(i,1,1)**2          * Delx 
              sum_ez   = sum_ez   + Ez(i,1,1)**2          * Delx 
              
              sum_exey = sum_exey + Ex(i,1,1) * Ey(i,1,1) * Delx 
              sum_exez = sum_exez + Ex(i,1,1) * Ez(i,1,1) * Delx 
              sum_eyez = sum_eyez + Ey(i,1,1) * Ez(i,1,1) * Delx 

              sum_bx   = sum_bx   + Bx (i,1,1)**2         * Delx 
              sum_by   = sum_by   + By (i,1,1)**2         * Delx 
              sum_bz   = sum_bz   + Bz (i,1,1)**2         * Delx 

              max_bx(i,j)  = Bx (i,1,1) / B0
              max_phi(i,j) = phi(i,1,1)

              sum_bxby = sum_bxby + Bx(i,1,1) * By(i,1,1) * Delx 
              sum_bxbz = sum_bxbz + Bx(i,1,1) * Bz(i,1,1) * Delx 
              sum_bybz = sum_bybz + By(i,1,1) * Bz(i,1,1) * Delx 

              sum_bxbybz = sum_bxbybz + (Bx(i,1,1)**2 + By(i,1,1)**2 + Bz(i,1,1)**2) * Delx 

              sum_vxvyvz = sum_vxvyvz + (Vx(i,1,1)**2 + Vy(i,1,1)**2 + Vz(i,1,1)**2) * Delx 

              sum_vx   = sum_vx   + Vx (i,1,1)**2   * Delx 
              sum_vy   = sum_vy   + Vy (i,1,1)**2   * Delx 
              sum_vz   = sum_vz   + Vz (i,1,1)**2   * Delx 

              sum_vxvy = sum_vxvy + Vx(i,1,1) * Vy(i,1,1) * Delx 
              sum_vxvz = sum_vxvz + Vx(i,1,1) * Vz(i,1,1) * Delx 
              sum_vyvz = sum_vyvz + Vy(i,1,1) * Vz(i,1,1) * Delx 
              
              sum_q    = sum_q    + q   (i,1,1)           * Delx 
              sum_p    = sum_p    + p   (i,1,1)           * Delx 
              sum_rho  = sum_rho  + rho (i,1,1)           * Delx 
              sum_D    = sum_D    + D   (i,1,1)           * Delx 


              epsiln(i,1,1) = p(i,1,1)/((gamma-1.d0)*rho(i,1,1))
              enthpy(i,1,1)  = rho(i,1,1) * (1.d0+epsiln(i,1,1)) + p(i,1,1) 
              W(i,1,1)       = 1.d0 /sqrt(1.d0 - (Vx(i,1,1)**2 + Vy(i,1,1)**2 + Vz(i,1,1)**2))
                 
              E2(i,1,1) = (Ex(i,1,1)**2 + Ey(i,1,1)**2 + Ez(i,1,1)**2)
              B2(i,1,1) = (By(i,1,1)**2 + Bz(i,1,1)**2)

              
              ekin_c   = ekin_c   + &
              0.5d0 * rho(i,1,1) * (Vx(i,1,1)**2 + Vy(i,1,1)**2 + Vz(i,1,1)**2)   * Delx 
              ekin_r   = ekin_r   + D(i,1,1) * (W(i,1,1) - 1.d0)                  * Delx 
              eint_c   = eint_c   + rho(i,1,1) * epsiln(i,1,1)                   * Delx 
              eint_r   = eint_r   + &
              W(i,1,1)**2 * (rho(i,1,1) * epsiln(i,1,1) + p(i,1,1)) - p(i,1,1)   * Delx 
              etot     = etot     + tau(i,1,1)                                    * Delx 
              efluid_r = efluid_r + enthpy(i,1,1) * W(i,1,1)**2 - p(i,1,1)        * Delx 
              emag     = emag   + 0.5 * (E2(i,1,1) + B2(i,1,1))                   * Delx 
              lorentz  = lorentz  + W(i,1,1)                                      * Delx 
              
        end do    ! for i

        imed1 = floor(0.25d0*imax)
        imed2 = floor(0.75d0*imax)

        do i=imed1,imed2
           do j=0,jmax

              
              sum_D_med    = sum_D_med    + D   (i,1,1)           * Delx 

              
           end do ! for j
        end do    ! for i

        max_val_bx  = maxval(max_bx)
        max_val_phi = maxval(max_phi)
        
           write(80,*)  h, t_advance, sum_psi
           write(81,*)  h, t_advance, sum_phi
           write(82,*)  h, t_advance, sum_ex
           write(83,*)  h, t_advance, sum_ey
           write(84,*)  h, t_advance, sum_ez
           write(85,*)  h, t_advance, sum_exey
           write(86,*)  h, t_advance, sum_exez
           write(87,*)  h, t_advance, sum_eyez
           write(88,*)  h, t_advance, sum_bx
           write(89,*)  h, t_advance, sum_bxbybz/sum_b20

           write(90,*)  h, t_advance, sum_by
           write(91,*)  h, t_advance, sum_bz
           write(92,*)  h, t_advance, sum_bxby
           write(93,*)  h, t_advance, sum_bxbz
           write(94,*)  h, t_advance, sum_bybz
           write(95,*)  h, t_advance, sum_vxvyvz
           write(96,*)  h, t_advance, sum_bxbybz
           write(97,*)  h, t_advance, sum_vx
           write(98,*)  h, t_advance, sum_vy
           write(99,*)  h, t_advance, sum_vz
           
           write(200,*) h, t_advance, sum_vxvy
           write(201,*) h, t_advance, sum_vxvz
           write(202,*) h, t_advance, sum_vyvz
           write(203,*) h, t_advance, sum_q
           write(204,*) h, t_advance, sum_p
           write(205,*) h, t_advance, sum_rho
           write(206,*) h, t_advance, ekin_c
           write(207,*) h, t_advance, ekin_r
           write(208,*) h, t_advance, etot
           write(209,*) h, t_advance, efluid_r
           write(210,*) h, t_advance, emag
           write(211,*) h, t_advance, eint_c
           write(212,*) h, t_advance, eint_r
           write(213,*) h, t_advance, lorentz
           write(214,*) h, t_advance, sum_D
           write(215,*) h, t_advance, sum_D_med

           write(216,*) h, t_advance, max_val_bx
           write(217,*) h, t_advance, max_val_phi


           by_imed = 0
           
           do j=0,jmax

              by_imed = by_imed + By(imed+1,1,1)

           end do

        end if  ! END IF Evolution of different quatities for 1D REC, TM, MD and SL


     else if (DIM == 2) then

       !Use internal files to convert integers to string values
       ! and embed it in the file names
       ! write (inputFileName, '(a, I0, a)') 'temp.time.',i,'.sa'

       if ( mod(h,hloc) == 0 ) then

! the expression mod(n,m) gives the remainder when n is divided by m; it is meant to be applied mainly to integers. Examples are
! mod(8,3) = 2   ,   mod(27,4) = 3   ,   mod(11,2) = 1   ,   mod(20,5) = 0  .

       hrec = h / hloc

       if ( hrec .lt. 10) then
          write (outputFileName , '(a, I0, a)') 'data_local/2d_data_psi_000' ,hrec,'.dat'  
          write (outputFileName1, '(a, I0, a)') 'data_local/2d_data_phi_000' ,hrec,'.dat'  
          write (outputFileName2, '(a, I0, a)') 'data_local/2d_data_vx_000'  ,hrec,'.dat'  
          write (outputFileName3, '(a, I0, a)') 'data_local/2d_data_vy_000'  ,hrec,'.dat'  
          write (outputFileName4, '(a, I0, a)') 'data_local/2d_data_vz_000'  ,hrec,'.dat'
          write (outputFileName5, '(a, I0, a)') 'data_local/2d_data_bx_000'  ,hrec,'.dat'
          write (outputFileName6, '(a, I0, a)') 'data_local/2d_data_by_000'  ,hrec,'.dat'
          write (outputFileName7, '(a, I0, a)') 'data_local/2d_data_bz_000'  ,hrec,'.dat'
          write (outputFileName8, '(a, I0, a)') 'data_local/2d_data_ex_000'  ,hrec,'.dat'
          write (outputFileName9, '(a, I0, a)') 'data_local/2d_data_ey_000'  ,hrec,'.dat'
          write (outputFileName10,'(a, I0, a)') 'data_local/2d_data_ez_000'  ,hrec,'.dat'
          write (outputFileName11,'(a, I0, a)') 'data_local/2d_data_p_000'   ,hrec,'.dat'
          write (outputFileName12,'(a, I0, a)') 'data_local/2d_data_q_000'   ,hrec,'.dat'
          write (outputFileName13,'(a, I0, a)') 'data_local/2d_data_rho_000' ,hrec,'.dat'
          write (outputFileName14,'(a, I0, a)') 'data_local/2d_data_divE_000',hrec,'.dat'
          write (outputFileName15,'(a, I0, a)') 'data_local/2d_data_divB_000',hrec,'.dat'
       else if  ( hrec .lt. 100) then
          write (outputFileName , '(a, I0, a)') 'data_local/2d_data_psi_00'  ,hrec,'.dat'  
          write (outputFileName1, '(a, I0, a)') 'data_local/2d_data_phi_00'  ,hrec,'.dat'  
          write (outputFileName2, '(a, I0, a)') 'data_local/2d_data_vx_00'   ,hrec,'.dat'  
          write (outputFileName3, '(a, I0, a)') 'data_local/2d_data_vy_00'   ,hrec,'.dat'  
          write (outputFileName4, '(a, I0, a)') 'data_local/2d_data_vz_00'   ,hrec,'.dat'
          write (outputFileName5, '(a, I0, a)') 'data_local/2d_data_bx_00'   ,hrec,'.dat'
          write (outputFileName6, '(a, I0, a)') 'data_local/2d_data_by_00'   ,hrec,'.dat'
          write (outputFileName7, '(a, I0, a)') 'data_local/2d_data_bz_00'   ,hrec,'.dat'
          write (outputFileName8, '(a, I0, a)') 'data_local/2d_data_ex_00'   ,hrec,'.dat'
          write (outputFileName9, '(a, I0, a)') 'data_local/2d_data_ey_00'   ,hrec,'.dat'
          write (outputFileName10,'(a, I0, a)') 'data_local/2d_data_ez_00'   ,hrec,'.dat'
          write (outputFileName11,'(a, I0, a)') 'data_local/2d_data_p_00'    ,hrec,'.dat'
          write (outputFileName12,'(a, I0, a)') 'data_local/2d_data_q_00'    ,hrec,'.dat'
          write (outputFileName13,'(a, I0, a)') 'data_local/2d_data_rhoo_00'  ,hrec,'.dat'
          write (outputFileName14,'(a, I0, a)') 'data_local/2d_data_divE_00' ,hrec,'.dat'
          write (outputFileName15,'(a, I0, a)') 'data_local/2d_data_divB_00' ,hrec,'.dat'
       else if (hrec .lt. 1000) then
          write (outputFileName , '(a, I0, a)') 'data_local/2d_data_psi_0'  ,hrec,'.dat'  
          write (outputFileName1, '(a, I0, a)') 'data_local/2d_data_phi_0'  ,hrec,'.dat'  
          write (outputFileName2, '(a, I0, a)') 'data_local/2d_data_vx_0'   ,hrec,'.dat'  
          write (outputFileName3, '(a, I0, a)') 'data_local/2d_data_vy_0'   ,hrec,'.dat'  
          write (outputFileName4, '(a, I0, a)') 'data_local/2d_data_vz_0'   ,hrec,'.dat'
          write (outputFileName5, '(a, I0, a)') 'data_local/2d_data_bx_0'   ,hrec,'.dat'
          write (outputFileName6, '(a, I0, a)') 'data_local/2d_data_by_0'   ,hrec,'.dat'
          write (outputFileName7, '(a, I0, a)') 'data_local/2d_data_bz_0'   ,hrec,'.dat'
          write (outputFileName8, '(a, I0, a)') 'data_local/2d_data_ex_0'   ,hrec,'.dat'
          write (outputFileName9, '(a, I0, a)') 'data_local/2d_data_ey_0'   ,hrec,'.dat'
          write (outputFileName10,'(a, I0, a)') 'data_local/2d_data_ez_0'   ,hrec,'.dat'
          write (outputFileName11,'(a, I0, a)') 'data_local/2d_data_p_0'    ,hrec,'.dat'
          write (outputFileName12,'(a, I0, a)') 'data_local/2d_data_q_0'    ,hrec,'.dat'
          write (outputFileName13,'(a, I0, a)') 'data_local/2d_data_rho_0'  ,hrec,'.dat'
          write (outputFileName14,'(a, I0, a)') 'data_local/2d_data_divE_0' ,hrec,'.dat'
          write (outputFileName15,'(a, I0, a)') 'data_local/2d_data_divB_0' ,hrec,'.dat'
       else
          write (outputFileName , '(a, I0, a)') 'data_local/2d_data_psi_'  ,hrec,'.dat'  
          write (outputFileName1, '(a, I0, a)') 'data_local/2d_data_phi_'  ,hrec,'.dat'  
          write (outputFileName2, '(a, I0, a)') 'data_local/2d_data_vx_'   ,hrec,'.dat'  
          write (outputFileName3, '(a, I0, a)') 'data_local/2d_data_vy_'   ,hrec,'.dat'  
          write (outputFileName4, '(a, I0, a)') 'data_local/2d_data_vz_'   ,hrec,'.dat'
          write (outputFileName5, '(a, I0, a)') 'data_local/2d_data_bx_'   ,hrec,'.dat'
          write (outputFileName6, '(a, I0, a)') 'data_local/2d_data_by_'   ,hrec,'.dat'
          write (outputFileName7, '(a, I0, a)') 'data_local/2d_data_bz_'   ,hrec,'.dat'
          write (outputFileName8, '(a, I0, a)') 'data_local/2d_data_ex_'   ,hrec,'.dat'
          write (outputFileName9, '(a, I0, a)') 'data_local/2d_data_ey_'   ,hrec,'.dat'
          write (outputFileName10,'(a, I0, a)') 'data_local/2d_data_ez_'   ,hrec,'.dat'
          write (outputFileName11,'(a, I0, a)') 'data_local/2d_data_p_'    ,hrec,'.dat'
          write (outputFileName12,'(a, I0, a)') 'data_local/2d_data_q_'    ,hrec,'.dat'
          write (outputFileName13,'(a, I0, a)') 'data_local/2d_data_rho_'  ,hrec,'.dat'
          write (outputFileName14,'(a, I0, a)') 'data_local/2d_data_divE_' ,hrec,'.dat'
          write (outputFileName15,'(a, I0, a)') 'data_local/2d_data_divB_' ,hrec,'.dat'
       end if

        !
        ! It will bail out rather than overwrite an existing file.
        !
        open(outFileUnit, file=outputFileName,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName)
           stop
        endif
        open(outFileUnit1, file=outputFileName1,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName1)
           stop
        endif
        open(outFileUnit2, file=outputFileName2,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName2)
           stop
        endif
        write(*, '("Opened file ", a, " for writing.")') &
                  trim(outputFileName2)
        open(outFileUnit3, file=outputFileName3,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName3)
           stop
        endif
        open(outFileUnit4, file=outputFileName4,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName4)
           stop
        endif
        open(outFileUnit5, file=outputFileName5,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName5)
           stop
        endif
        open(outFileUnit6, file=outputFileName6,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName6)
           stop
        endif
        open(outFileUnit7, file=outputFileName7,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName7)
           stop
        endif
        open(outFileUnit8, file=outputFileName8,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName8)
           stop
        endif
        open(outFileUnit9, file=outputFileName9,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName9)
           stop
        endif
        open(outFileUnit10, file=outputFileName10,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName10)
           stop
        endif
        open(outFileUnit11, file=outputFileName11,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName11)
           stop
        endif
        open(outFileUnit12, file=outputFileName12,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName12)
           stop
        endif
        open(outFileUnit13, file=outputFileName13,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName13)
           stop
        endif
        open(outFileUnit14, file=outputFileName14,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName14)
           stop
        endif
        open(outFileUnit15, file=outputFileName15,status="new", iostat=ios)
        if (ios /= 0) then
           write(*, '("Can''t open file ", a, " for writing.")') &
                     trim(outputFileName15)
           stop
        endif
!!$        open(outFileUnit16, file=outputFileName16,status="new", iostat=ios)
!!$        if (ios /= 0) then
!!$           write(*, '("Can''t open file ", a, " for writing.")') &
!!$                     trim(outputFileName16)
!!$           stop
!!$        endif

!**********************************************************************************************
       !FOR RDTM
       do i=0,imax
           do j=0,jmax



         if (sqrt(Ex(i,j,1)**2+Ey(i,j,1)**2+Ez(i,j,1)**2) .eq. 0.d0) then
                 
            divE_E(i,j,1) = 0.d0
         else
                 
            divE_E(i,j,1) = abs( (Ex(i  ,j  ,1) - Ex(i-1,j  ,1))/Delx +   &
                                 (Ey(i  ,j  ,1) - Ey(i  ,j-1,1))/Dely +   &
                                 (Ez(i  ,j  ,1) - Ez(i  ,j  ,1))/Delz )/  &
                            (sqrt(Ex(i,j,1)**2+Ey(i,j,1)**2+Ez(i,j,1)**2))
         end if

         if (sqrt(Bx(i,j,1)**2+By(i,j,1)**2+Bz(i,j,1)**2) .eq. 0.d0) then

            divB_B(i,j,1)  = 0.d0

         else         

            divB_B(i,j,1) = abs( (Bx(i  ,j  ,1) - Bx(i-1,j  ,1))/Delx  +  &
                                 (By(i  ,j  ,1) - By(i  ,j-1,1))/Dely  +  &
                                 (Bz(i  ,j  ,1) - Bz(i  ,j  ,1))/Delz  )/ &
                            (sqrt(Bx(i,j,1)**2+By(i,j,1)**2+Bz(i,j,1)**2))

         end if

      end do
   end do



   !-----------------
   ! Y Direction
   !-----------------

   ! Left Boundary
   ! -----------------
       
       do i=-stencil,imax+stencil
          do j = 1, stencil

             divE_E(i,0-j,1) = divE_E(i,jmax-j,1)
             divB_B(i,0-j,1) = divB_B(i,jmax-j,1)

          end do
       end do

          ! Right  Boundary
          ! -----------------

       do i=-stencil,imax+stencil
          do j = 1, stencil

             
             divE_E(i,jmax+j,1) = divE_E(i,0+j,1)
             divB_B(i,jmax+j,1) = divB_B(i,0+j,1)
             
          end do
       end do

  ! Copy boundary x-direction

  !-----------------
  ! X Direction
  !-----------------

           ! Left  Boundary
           ! -----------------

       do i = 1,stencil
          do j=-stencil,jmax+stencil

             divE_E(0-i,j,1) = divE_E(0  ,j,1)
             divB_B(0-i,j,1) = divB_B(0  ,j,1) 

        end do
     end do

           ! Right  Boundary
           ! -----------------

        do i = 1,stencil
           do j=-stencil,jmax+stencil
              
              divE_E(imax+i,j,1) = divE_E(imax  ,j,1)
              divB_B(imax+i,j,1) = divB_B(imax  ,j,1) 

        end do
     end do

!**********************************************************************************************     
        
        ! Here's where you read stuff from the input file
        ! and write stuff to the output file.

        do i=-6,imax+6
           do j=-6,jmax+6


              if (j == -6) then
                 write (outFileUnit ,*)  ''
                 write (outFileUnit1,*)  ''
                 write (outFileUnit2,*)  ''
                 write (outFileUnit3,*)  ''
                 write (outFileUnit4,*)  ''
                 write (outFileUnit5,*)  ''
                 write (outFileUnit6,*)  ''
                 write (outFileUnit7,*)  ''
                 write (outFileUnit8,*)  ''
                 write (outFileUnit9,*)  ''
                 write (outFileUnit10,*) ''
                 write (outFileUnit11,*) ''
                 write (outFileUnit12,*) ''
                 write (outFileUnit13,*) ''
                 write (outFileUnit14,*) ''
                 write (outFileUnit15,*) ''
              end if

            
              write (outFileUnit ,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  psi(i,j,1)   
              write (outFileUnit1,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  phi(i,j,1)   
              write (outFileUnit2,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Vx (i,j,1)   
              write (outFileUnit3,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Vy (i,j,1)   
              write (outFileUnit4,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Vz (i,j,1)
              write (outFileUnit5,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Bx (i,j,1)
              write (outFileUnit6,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  By (i,j,1)
              write (outFileUnit7,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Bz (i,j,1)
              write (outFileUnit8,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Ex (i,j,1)
              write (outFileUnit9,*)  - 0.5d0 * Longx + i * Delx, j * Dely,  Ey (i,j,1)
              write (outFileUnit10,*) - 0.5d0 * Longx + i * Delx, j * Dely,  Ez (i,j,1)
              write (outFileUnit11,*) - 0.5d0 * Longx + i * Delx, j * Dely,  p  (i,j,1)
              write (outFileUnit12,*) - 0.5d0 * Longx + i * Delx, j * Dely,  q  (i,j,1)
              write (outFileUnit13,*) - 0.5d0 * Longx + i * Delx, j * Dely,  rho(i,j,1)
              write (outFileUnit14,*) - 0.5d0 * Longx + i * Delx, j * Dely,  divE_E(i,j,1)
              write (outFileUnit15,*) - 0.5d0 * Longx + i * Delx, j * Dely,  divB_B(i,j,1)
              


           end do
        end do

      

        ! Close files before next pass through the loop
        
        write(*, '("Closing files.")')
        close(outFileUnit)
        close(outFileUnit1)
        close(outFileUnit2)
        close(outFileUnit3)
        close(outFileUnit4)
        close(outFileUnit5)
        close(outFileUnit6)
        close(outFileUnit7)
        close(outFileUnit8)
        close(outFileUnit9)
        close(outFileUnit10)
        close(outFileUnit11)
        close(outFileUnit12)
        close(outFileUnit13)
        close(outFileUnit14)
        close(outFileUnit15)
!!$        close(outFileUnit16)        
        
        write(*,*) ""

     end if



!-------------------------------------------------------------
! Evolution of different quatities for 1D REC, TM, MD and SL
!-------------------------------------------------------------

     
     if ( h == 1 .or. mod(h,hglb) == 0 ) then

        
! Evolution of volumetric integral Bx**2 
!_________________________________________________________________________________________________________________________

     sum_psi  = 0.d0
     sum_phi  = 0.d0
     sum_ex   = 0.d0
     sum_ey   = 0.d0
     sum_ez   = 0.d0
     sum_exey = 0.d0
     sum_exez = 0.d0
     sum_eyez = 0.d0
     sum_bx   = 0.d0
     sum_by   = 0.d0
     sum_bz   = 0.d0
     sum_bxby = 0.d0
     sum_bxbz = 0.d0
     sum_bybz = 0.d0
     sum_vx   = 0.d0
     sum_vy   = 0.d0
     sum_vz   = 0.d0
     sum_vxvy = 0.d0
     sum_vxvz = 0.d0
     sum_vyvz = 0.d0
     sum_q    = 0.d0
     sum_p    = 0.d0
     sum_rho  = 0.d0
     sum_bxbybz = 0.d0
     sum_vxvyvz = 0.d0
     ekin_c     = 0.d0
     ekin_r     = 0.d0
     eint_c     = 0.d0
     eint_r     = 0.d0
     etot       = 0.d0
     efluid_r   = 0.d0
     emag       = 0.d0
     lorentz    = 0.d0
     sum_D      = 0.d0
     sum_D_med  = 0.d0

     do i=0,imax
        do j=0,jmax

              sum_psi  = sum_psi  + psi(i,j,1)            * Delx * Dely
              sum_phi  = sum_phi  + phi(i,j,1)            * Delx * Dely
              
              sum_ex   = sum_ex   + Ex(i,j,1)**2          * Delx * Dely
              sum_ey   = sum_ey   + Ey(i,j,1)**2          * Delx * Dely
              sum_ez   = sum_ez   + Ez(i,j,1)**2          * Delx * Dely
              
              sum_exey = sum_exey + Ex(i,j,1) * Ey(i,j,1) * Delx * Dely
              sum_exez = sum_exez + Ex(i,j,1) * Ez(i,j,1) * Delx * Dely
              sum_eyez = sum_eyez + Ey(i,j,1) * Ez(i,j,1) * Delx * Dely

              sum_bx   = sum_bx   + Bx (i,j,1)**2         * Delx * Dely
              sum_by   = sum_by   + By (i,j,1)**2         * Delx * Dely
              sum_bz   = sum_bz   + Bz (i,j,1)**2         * Delx * Dely

              max_bx(i,j)  = Bx (i,j,1) / B0
              max_phi(i,j) = phi(i,j,1)

              sum_bxby = sum_bxby + Bx(i,j,1) * By(i,j,1) * Delx * Dely
              sum_bxbz = sum_bxbz + Bx(i,j,1) * Bz(i,j,1) * Delx * Dely
              sum_bybz = sum_bybz + By(i,j,1) * Bz(i,j,1) * Delx * Dely

              sum_bxbybz = sum_bxbybz + (Bx(i,j,1)**2 + By(i,j,1)**2 + Bz(i,j,1)**2) *  Delx * Dely

              sum_vxvyvz = sum_vxvyvz + (Vx(i,j,1)**2 + Vy(i,j,1)**2 + Vz(i,j,1)**2) * Delx * Dely

              sum_vx   = sum_vx   + Vx (i,j,1)**2   * Delx * Dely
              sum_vy   = sum_vy   + Vy (i,j,1)**2   * Delx * Dely
              sum_vz   = sum_vz   + Vz (i,j,1)**2   * Delx * Dely

              sum_vxvy = sum_vxvy + Vx(i,j,1) * Vy(i,j,1) * Delx * Dely
              sum_vxvz = sum_vxvz + Vx(i,j,1) * Vz(i,j,1) * Delx * Dely
              sum_vyvz = sum_vyvz + Vy(i,j,1) * Vz(i,j,1) * Delx * Dely
              
              sum_q    = sum_q    + q   (i,j,1)           * Delx * Dely
              sum_p    = sum_p    + p   (i,j,1)           * Delx * Dely
              sum_rho  = sum_rho  + rho (i,j,1)           * Delx * Dely
              sum_D    = sum_D    + D   (i,j,1)           * Delx * Dely


              epsiln(i,j,1) = p(i,j,1)/((gamma-1.d0)*rho(i,j,1))
              enthpy(i,j,1)  = rho(i,j,1) * (1.d0+epsiln(i,j,1)) + p(i,j,1) 
              W(i,j,1)       = 1.d0 /sqrt(1.d0 - (Vx(i,j,1)**2 + Vy(i,j,1)**2 + Vz(i,j,1)**2))
                 
              E2(i,j,1) = (Ex(i,j,1)**2 + Ey(i,j,1)**2 + Ez(i,j,1)**2)
              B2(i,j,1) = (By(i,j,1)**2 + Bz(i,j,1)**2)

              
              ekin_c   = ekin_c   + &
              0.5d0 * rho(i,j,1) * (Vx(i,j,1)**2 + Vy(i,j,1)**2 + Vz(i,j,1)**2)   * Delx * Dely
              ekin_r   = ekin_r   + D(i,j,1) * (W(i,j,1) - 1.d0)                  * Delx * Dely
              eint_c   = eint_c   + rho(i,j,1) * epsiln(i,j,1)                   * Delx * Dely
              eint_r   = eint_r   + &
              W(i,j,1)**2 * (rho(i,j,1) * epsiln(i,j,1) + p(i,j,1)) - p(i,j,1)   * Delx * Dely
              etot     = etot     + tau(i,j,1)                                    * Delx * Dely
              efluid_r = efluid_r + enthpy(i,j,1) * W(i,j,1)**2 - p(i,j,1)        * Delx * Dely
              emag     = emag   + 0.5 * (E2(i,j,1) + B2(i,j,1))                   * Delx * Dely
              lorentz  = lorentz  + W(i,j,1)                                      * Delx * Dely
              
           end do ! for j
        end do    ! for i

        imed1 = floor(0.25d0*imax)
        imed2 = floor(0.75d0*imax)

        do i=imed1,imed2
           do j=0,jmax

              
              sum_D_med    = sum_D_med    + D   (i,j,1)           * Delx * Dely

              
           end do ! for j
        end do    ! for i

        max_val_bx  = maxval(max_bx)
        max_val_phi = maxval(max_phi)

       
           write(80,*)  h, t_advance, sum_psi
           write(81,*)  h, t_advance, sum_phi
           write(82,*)  h, t_advance, sum_ex
           write(83,*)  h, t_advance, sum_ey
           write(84,*)  h, t_advance, sum_ez
           write(85,*)  h, t_advance, sum_exey
           write(86,*)  h, t_advance, sum_exez
           write(87,*)  h, t_advance, sum_eyez
           write(88,*)  h, t_advance, sum_bx

           write(90,*)  h, t_advance, sum_by
           write(91,*)  h, t_advance, sum_bz
           write(92,*)  h, t_advance, sum_bxby
           write(93,*)  h, t_advance, sum_bxbz
           write(94,*)  h, t_advance, sum_bybz
           write(95,*)  h, t_advance, sum_vxvyvz
           write(96,*)  h, t_advance, sum_bxbybz
           write(97,*)  h, t_advance, sum_vx
           write(98,*)  h, t_advance, sum_vy
           write(99,*)  h, t_advance, sum_vz
           
           write(200,*) h, t_advance, sum_vxvy
           write(201,*) h, t_advance, sum_vxvz
           write(202,*) h, t_advance, sum_vyvz
           write(203,*) h, t_advance, sum_q
           write(204,*) h, t_advance, sum_p
           write(205,*) h, t_advance, sum_rho
           write(206,*) h, t_advance, ekin_c
           write(207,*) h, t_advance, ekin_r
           write(208,*) h, t_advance, etot
           write(209,*) h, t_advance, efluid_r
           write(210,*) h, t_advance, emag
           write(211,*) h, t_advance, eint_c
           write(212,*) h, t_advance, eint_r
           write(213,*) h, t_advance, lorentz
           write(214,*) h, t_advance, sum_D
           write(215,*) h, t_advance, sum_D_med

           write(216,*) h, t_advance, max_val_bx
           write(217,*) h, t_advance, max_val_phi

           


           by_imed = 0
           
           do j=0,jmax

              by_imed = by_imed + By(imed+1,j,1)

           end do

           by_prom = by_imed / jmax


! Evolution thickness of the current sheet in center of reconnection zone 
!_________________________________________________________________________________________________________________________
! For 2D test a_tm is equivalent to delta_rec TM V1 (TEST 7)


!!$           delta_rec_t = B0 * tanh(Delx/a_tm) * a_tm  / by_prom ! By(imed+1,0,1)
!!$
!!$           write(67,*) h, t_advance, delta_rec_t


! Evolution of the electric field in center of reconnection zone (reconection rate)
!_________________________________________________________________________________________________________________________


           posy   = - 0.5d0 * Ly + (jmed+1) * Dely
           r      =   sqrt(1.d0 - (posy**2/Lr**2))

           delta_rec_t = B0 * r * tanh(Delx/delta_rec) * delta_rec / By(imed+1,jmed+1,1)

           write(67,*) h, t_advance, delta_rec_t
           write(68,*) h, t_advance, Ez(imed,jmed,1) 


!_________________________________________________________________________________________________________________________


     end if

    else 

       write(*,*) "STOP subroutine lector"
       write(*,*) "This DIM is not implemented yet"
       stop

    end if

    if (DIM == 1 ) then
    
       if ( mod(h,hloc) == 0 ) then

          hbck = h/hloc

          if ( hbck .lt. 10) then 
             write (outputFileName17, '(a, I0, a)') 'data_restart/all_variables_000',hbck,'.dat'
             write (outputFileName18, '(a, I0, a)') 'data_print/print_backup_000',hbck,'.dat'  
          else if  ( hbck .lt. 100) then
             write (outputFileName17, '(a, I0, a)') 'data_restart/all_variables_00',hbck,'.dat'
             write (outputFileName18, '(a, I0, a)') 'data_print/print_backup_00',hbck,'.dat'  
          else if (hbck .lt. 1000) then
             write (outputFileName17, '(a, I0, a)') 'data_restart/all_variables_0',hbck,'.dat'
             write (outputFileName18, '(a, I0, a)') 'data_print/print_backup_0',hbck,'.dat'  
          else
             write (outputFileName17, '(a, I0, a)') 'data_restart/all_variables_',hbck,'.dat'
             write (outputFileName18, '(a, I0, a)') 'data_print/print_backup_',hbck,'.dat'  
          end if

        !
        ! It will bail out rather than overwrite an existing file.
        !
          open(outFileUnit17, file=outputFileName17,status="new", iostat=ios)
          if (ios /= 0) then
             write(*, '("Can''t open file ", a, " for writing.")') &
                   trim(outputFileName17)
             stop
          endif
          open(outFileUnit18, file=outputFileName18,status="new", iostat=ios)
          if (ios /= 0) then
             write(*, '("Can''t open file ", a, " for writing.")') &
                   trim(outputFileName18)
             stop
          endif

          write(*, '("Opened file ", a, " for writing.")') &
                  trim(outputFileName17)
          write(*, '("Opened file ", a, " for writing.")') &
                  trim(outputFileName18)

        ! Here's where you read stuff from the input file
        ! and write stuff to the output file.


         do i=-6,imax+6
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!::::::::: WRITE BACK UP FILE ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            write (outFileUnit17,*)   psi(i,1,1), phi(i,1,1), Bx (i,1,1), By (i,1,1), Bz (i,1,1), Vx (i,1,1), &
                                      Vy (i,1,1), Vz (i,1,1), Ex (i,1,1), Ey (i,1,1), Ez (i,1,1), q  (i,1,1), &
                                      rho(i,1,1), p  (i,1,1)
            
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!::::::::: WRITE PRINT FILES :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            write (outFileUnit18,*)   i, h, Delt, hini, hmax, hglb, hloc, CFL, sigma, kappa_psi, kappa_phi, gamma,          & !12
                                      hini * Delt, hmax * Delt, hglb * Delt, hloc * Delt, t_advance, imax, Delx, Longx,     & !20
                                      i * Delx,- 0.5d0 * Longx + i * Delx, psi(i,1,1), phi(i,1,1), Bx (i,1,1), By (i,1,1),  & !26
                                      Bz (i,1,1), Vx (i,1,1), Vy (i,1,1), Vz (i,1,1), Ex (i,1,1), Ey (i,1,1),               & !32
                                      Ez (i,1,1), q  (i,1,1), rho(i,1,1), p  (i,1,1)                                          !36

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            
     end do ! for i

        write(*,*)''
        write(*, '("Writing in backup file.")')
        write(*,*) "h=", h
        write(*,*)''

     end if

    else if (DIM ==2 ) then
    
    
       if ( mod(h,hloc)==0 ) then

          hbck = h / hloc


          if ( hbck .lt. 10) then 
             write (outputFileName17, '(a, I0, a)') 'data_restart/all_variables_000',hbck,'.dat'
             write (outputFileName18, '(a, I0, a)') 'data_print/print_backup_000',hbck,'.dat'
          else if  ( hbck .lt. 100) then
             write (outputFileName17, '(a, I0, a)') 'data_restart/all_variables_00',hbck,'.dat'
             write (outputFileName18, '(a, I0, a)') 'data_print/print_backup_00',hbck,'.dat'  
          else if (hbck .lt. 1000) then
             write (outputFileName17, '(a, I0, a)') 'data_restart/all_variables_0',hbck,'.dat'
             write (outputFileName18, '(a, I0, a)') 'data_print/print_backup_0',hbck,'.dat'  
          else
             write (outputFileName17, '(a, I0, a)') 'data_restart/all_variables_',hbck,'.dat'
             write (outputFileName18, '(a, I0, a)') 'data_print/print_backup_',hbck,'.dat'  
          end if

        !
        ! It will bail out rather than overwrite an existing file.
        !
          open(outFileUnit17, file=outputFileName17,status="new", iostat=ios)
          if (ios /= 0) then
             write(*, '("Can''t open file ", a, " for writing.")') &
                   trim(outputFileName17)
             stop
          endif
          open(outFileUnit18, file=outputFileName18,status="new", iostat=ios)
          if (ios /= 0) then
             write(*, '("Can''t open file ", a, " for writing.")') &
                   trim(outputFileName18)
             stop
          endif

          write(*, '("Opened file ", a, " for writing.")') &
                  trim(outputFileName17)
          write(*, '("Opened file ", a, " for writing.")') &
                  trim(outputFileName18)

         ! Here's where you read stuff from the input file
        ! and write stuff to the output file.


        do i=-6,imax+6
           do j=-6,jmax+6
              
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!::::::::: WRITE BACK UP FILE:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

              write (outFileUnit17,*)   psi(i,j,1), phi(i,j,1), Bx (i,j,1), By (i,j,1), Bz (i,j,1), Vx (i,j,1), &
                                       Vy (i,j,1), Vz (i,j,1), Ex (i,j,1), Ey (i,j,1), Ez (i,j,1), q  (i,j,1), &
                                       rho(i,j,1), p  (i,j,1)

! Remember when you use RESTART option, you must assignate a value to "sigma" in TMs test this value is ussually "sigma_0"
              
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!::::::::: WRITE PRINT FILES :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

             
             if (j == -6) then
                 write (outFileUnit18,*) ''
              end if

              write (outFileUnit18,*)   i,j, h, hini, hmax, hglb, hloc,CFL, sigma, kappa_psi, kappa_phi, gamma, & !12
                                        Delt, hini * Delt, hmax * Delt, hglb * Delt, hloc * Delt, t_advance,    & !18
                                        imax, jmax, Delx, Dely, Longx, Longy, i * Delx, j * Dely,               & !26
                                      - 0.5d0 * Longx + i * Delx, - 0.5d0 * Longy + j * Dely,                   & !28
                                        psi(i,j,1), phi(i,j,1), Bx (i,j,1), By (i,j,1), Bz (i,j,1), Vx (i,j,1), & !34
                                        Vy (i,j,1), Vz (i,j,1), Ex (i,j,1), Ey (i,j,1), Ez (i,j,1),             & !39
                                        q  (i,j,1), rho(i,j,1), p  (i,j,1), divB_B(i,j,1)                         !43

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
              
            end do ! for j
       end do ! for i             

        write(*,*)''
        write(*, '("Writing in backup file.")')
        write(*,*) "h=", h
        write(*,*)''
!!$        close(68)
!!$        close(69)

     end if

  end if


  
!--------FORMAT----------
120  format (E21.16,E21.16)
!--------FORMAT----------

!--------FORMAT----------
121  format (E21.16,E21.16,E21.16)
!--------FORMAT----------

!--------FORMAT----------
122 format (I15,E21.16,E21.16)
!--------FORMAT----------

!--------FORMAT----------
123  format (14E21.16)
!--------FORMAT----------

  end subroutine lector
