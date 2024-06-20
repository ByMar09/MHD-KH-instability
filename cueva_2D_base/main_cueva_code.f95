! Main CUEVA CODE

  program IMEXRK

    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

    t_init = 0
       
    call openfiles

    call open_backup
       
    if (RESTART == 0) then
       
       call init_var_primitive
       
    else if (RESTART == 1) then
       
       call init_restart

       call boundary

       close (222)

    else

       write(*,*) "STOP: main"
       write(*,*) "RESTART parameter is not valid"
       stop
       
    end if

    call interface


  !     ********************
  !     Starts Time Loop
  !     ********************

    t_advance = (hini-1) * Delt !0.d0

  do h=hini,hmax-1

              !     ------------------------------------
              !      Starts Boundary Conditions
              !     ------------------------------------

     call boundary
     
     
!::::::::: Monitorize  Initial data including phantom cells :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

     if (h == hini .or. h == hini .and. RESTART ==1) then
     
        open (unit = 333, action="write", file = "data_print/print_backup_0000.dat", IOSTAT=ierror)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!::::::::: WRITE PRINT FILES :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        if (DIM ==1) then
           
           do i=-6,imax+6
              
              write (333,*)   i, h, Delt, hini, hmax, hglb, hloc, CFL, sigma, kappa_psi, kappa_phi, gamma,          & !12
                              hini * Delt, hmax * Delt, hglb * Delt, hloc * Delt, t_advance, imax, Delx, Longx,     & !20
                              i * Delx,- 0.5d0 * Longx + i * Delx, psi(i,1,1), phi(i,1,1), Bx (i,1,1), By (i,1,1),  & !26
                              Bz (i,1,1), Vx (i,1,1), Vy (i,1,1), Vz (i,1,1), Ex (i,1,1), Ey (i,1,1),               & !32
                              Ez (i,1,1), q  (i,1,1), rho(i,1,1), p  (i,1,1)                                          !36

         end do

      else if (DIM==2) then 

        do i=-6,imax+6
           do j=-6+0+fix,jmax+6


             if (j == -6) then
                 write (333,*) ''
              end if

              write (333,*)  i,j, h, hini, hmax, hglb, hloc,CFL, sigma, kappa_psi, kappa_phi, gamma, & !12
                             Delt, hini * Delt, hmax * Delt, hglb * Delt, hloc * Delt, h * Delt,     & !18
                             imax, jmax, Delx, Dely, Longx, Longy, i * Delx, j * Dely,               & !26
                           - 0.5d0 * Longx + i * Delx, - 0.5d0 * Longy + j * Dely,                   & !28
                             psi(i,j,1), phi(i,j,1), Bx (i,j,1), By (i,j,1), Bz (i,j,1), Vx (i,j,1), & !34
                             Vy (i,j,1), Vz (i,j,1), Ex (i,j,1), Ey (i,j,1), Ez (i,j,1),             & !39
                             q  (i,j,1), rho(i,j,1), p  (i,j,1)

            end do ! for j
         end do ! for i

      end if

    end if
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        

!:::::::::::::: Changing the time step to, Delt = CFL * Deltax / (1 + max(v/v_A)) :::::::::::::::::::::::::::::::::::::::::::::

    
!!$           do j=0+fix,jmax
!!$              do i=0,imax
!!$
!!$                 enthalpy_A = rho(i,j,1) +  gamma * p(i,j,1) / (gamma - 1.d0) !Eq.2.20 Tesis
!!$
!!$                 E2B2 = Bx(i,j,1)**2+By(i,j,1)**2+Bz(i,j,1)**2 &
!!$                      - (Ex(i,j,1)**2+Ey(i,j,1)**2+Ez(i,j,1)**2)
!!$
!!$                 VA = sqrt( E2B2 / ( enthalpy_A + E2B2) )                       ! Eq.2.22 Tesis
!!$
!!$                 max_vx(i,j) = sqrt(Vx(i,j,1)**2+Vy(i,j,1)**2+Vz(i,j,1)**2) !/ VA
!!$                 
!!$                 max_va(i,j) = VA
!!$
!!$              end do
!!$           end do
!!$
!!$           j =1 
!!$
!!$           max_val_v_va = maxval(max_vx)
!!$           max_val_va   = maxval(max_va)
!!$
!!$        if (change_time_step == 1) then
!!$
!!$
!!$!           Delt  = CFL *  min(Delx, Dely) / (1.d0 + max_val_v_va)
!!$           
!!$           Wa = 1.d0 / sqrt(1.d0 - max_val_va**2)
!!$
!!$!           Delt  = sqrt(2.35d0 / sigma)
!!$!           Delt  = sqrt(Wa * sqrt(VA) / sigma)
!!$
!!$           Delt  = sqrt(Wa /(max_val_va**2 * sigma))
!!$
!!$           if (Delt .ge. Delt1) Delt = Delt1
!!$
!!$           Delt = Delt1 * Delt / sqrt(Delt1**2 + Delt**2)
!!$
!!$        end if

!!$        if (h .le. 4) then
!!$           
!!$           Delt = Delt0
!!$
!!$        else
!!$
!!$           Delt = Delt1
!!$
!!$        end if


        t_advance = t_advance + Delt


!!$        print*, "t=", t_advance, "Delt =", Delt, "Va =", VA, "max(v/va)=", max_val_v_va
    

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  
        if (primitive_rec == 1) then

           call recons_var_primitive
           call boundary

        end if

              !     ------------------------------------
              !      Starts Conserved Variables
              !     ------------------------------------

        call startconvar


        call cleanvarint

        !clean index for fix-up density and presion subroutine 12_varprimitive.f95
        !=====================================

        contador_r = 0
  
        contador_p = 0
          
        contador_v = 0

        call swap1

              !     ------------------------------------
              !      Starts Intermediate Steps
              !     ------------------------------------


     do l=1,lmax


        do m=1,l-1


           call sumint


        end do  ! for m


     call varint

     call  electricvarint

! ------------- Recovery Primitive Variables Intermediate Steps ---------

     if (REC_PRIM == 1) then

       l_swap = 0

     call varprimitivecardano

     call interconvar

     end if

!----------- Boundary conditions for intermediate conserved variables ----------     

     call boundary_conserved
     
     call boundary_electric

!-------------------------------------------------------------------------------

     call flowinit

!----------- Boundary conditions for fluxes ( not standard procedure in literature ! )  ----------     
     
     call boundary_flux_conserved

! --------------------- Flux reconstruction ---------------------


     if      (flux_solver == 1) then

        call lax_flow

     else if (flux_solver == 2) then

        call hll_flow

     else if (flux_solver == 3) then

        call hlle_flow

     else if (flux_solver == 4) then

        call hllc_flow

     else if (flux_solver == 5) then

        call hll_flow_mpx

     else if (flux_solver == 6) then

        call hlle_flow_mpx

     else if (flux_solver == 7) then

        call hllc_flow_mpx
     
     else

        write(*,*) "STOP: main"
        write(*,*) "flux_solver parameter is not valid"
        stop

     end if
     ! --------------------- Flux reconstruction ---------------------

     call fluxandsource

     call sumfin


  end do !for l


      l = lmax


     call varfin

     call electricfield

!     call swap2

     if (REC_PRIM == 1) then

        l_swap = 2

     else

        l_swap = 1
        
     end if


     if (conserved_to_primitive == 1) then

        call varprimitivecardano

     else if (conserved_to_primitive == 0) then

        call varprimitive_nr

     else

        write(*,*) "STOP: main"
        write(*,*) "conserved_to_primitive is a not valid parameter"
        stop

     end if

     call  boundary
     
     t_init = 1
     
     call lector


     if ( h == 1 .or. mod(h,hglb) == 0 ) then

        print*, ' ::::::::::: loop ::::::::::::   '
        print*, 'iteration =',h,'  ', 'time =', t_advance
!        print*, 'iteration =',h,'  ', 'time =', h * Delt
        print*, ' ::::::::::: end  ::::::::::::   '
        print*, ''
        call flush (6)

     end if


!!$     if( t_advance .ge. 201.d0 ) then
!!$
!!$        print*, "Stop at t=", t_advance
!!$        goto 073
!!$
!!$     end if



  end do ! for h


  call fin_interface

  close (26)
  close (27)
  close (28)
  close (29)
  close (30)
  close (31)
  close (32)
  close (33)
  close (34)
  close (35)
  close (36)
  close (37)
  close (68)

 end program IMEXRK
