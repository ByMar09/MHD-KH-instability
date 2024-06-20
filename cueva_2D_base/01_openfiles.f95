
    !   / --------------------\
    ! |||     OPEN FILES      ||| 
    !   \---------------------/

  subroutine openfiles 

    use parameters
    
    implicit none


    if (DIM == 1) then

       open (unit = 37, file = "data_global/test_1D_cpaw_vel_hllc.dat")

       open (unit = 80, file = "data_global/1d_data_psi_integral_vol.dat")
       open (unit = 81, file = "data_global/1d_data_phi_integral_vol.dat")
       open (unit = 82, file = "data_global/1d_data_ex_integral_vol.dat")
       open (unit = 83, file = "data_global/1d_data_ey_integral_vol.dat")
       open (unit = 84, file = "data_global/1d_data_ez_integral_vol.dat")
       open (unit = 85, file = "data_global/1d_data_exey_integral_vol.dat")
       open (unit = 86, file = "data_global/1d_data_exez_integral_vol.dat")
       open (unit = 87, file = "data_global/1d_data_eyez_integral_vol.dat")
       open (unit = 88, file = "data_global/1d_data_bx_integral_vol.dat")
       open (unit = 89, file = "data_global/1d_data_b2_b20_integral_vol.dat")

       open (unit = 90, file = "data_global/1d_data_by_integral_vol.dat")
       open (unit = 91, file = "data_global/1d_data_bz_integral_vol.dat")
       open (unit = 92, file = "data_global/1d_data_bxby_integral_vol.dat")
       open (unit = 93, file = "data_global/1d_data_bxbz_integral_vol.dat")
       open (unit = 94, file = "data_global/1d_data_bybz_integral_vol.dat")
       open (unit = 95, file = "data_global/1d_data_v2_integral_vol.dat")
       open (unit = 96, file = "data_global/1d_data_b2_integral_vol.dat")
       open (unit = 97, file = "data_global/1d_data_vx_integral_vol.dat")
       open (unit = 98, file = "data_global/1d_data_vy_integral_vol.dat")
       open (unit = 99, file = "data_global/1d_data_vz_integral_vol.dat")
       
       open (unit =200, file = "data_global/1d_data_vxvy_integral_vol.dat")
       open (unit =201, file = "data_global/1d_data_vxvz_integral_vol.dat")
       open (unit =202, file = "data_global/1d_data_vyvz_integral_vol.dat")
       open (unit =203, file = "data_global/1d_data_q_integral_vol.dat")
       open (unit =204, file = "data_global/1d_data_p_integral_vol.dat")
       open (unit =205, file = "data_global/1d_data_rho_integral_vol.dat")
       open (unit =206, file = "data_global/1d_data_ekin_c_integral_vol.dat")
       open (unit =207, file = "data_global/1d_data_ekin_r_integral_vol.dat")
       open (unit =208, file = "data_global/1d_data_etot_integral_vol.dat")
       open (unit =209, file = "data_global/1d_data_efluid_r_integral_vol.dat")
       open (unit =210, file = "data_global/1d_data_emag_integral_vol.dat")
       open (unit =211, file = "data_global/1d_data_eint_c_integral_vol.dat")
       open (unit =212, file = "data_global/1d_data_eint_r_integral_vol.dat")
       open (unit =213, file = "data_global/1d_data_lorentz_integral_vol.dat")
       open (unit =214, file = "data_global/1d_data_D_integral_vol.dat")
       open (unit =215, file = "data_global/1d_data_D_zona_media_integral_vol.dat")

       open (unit =216, file = "data_global/1d_data_bx_max_val.dat")
       open (unit =217, file = "data_global/1d_data_phi_max_val.dat")


    else if (DIM == 2) then

       open (unit = 67, file = "data_global/delta_rec_t_2D.dat")
       open (unit = 68, file = "data_global/rate_rec_ez_2D.dat")
       open (unit = 70, file = "data_global/comparative_Jz.dat")

       open (unit = 80, file = "data_global/2d_data_psi_integral_vol.dat")
       open (unit = 81, file = "data_global/2d_data_phi_integral_vol.dat")
       open (unit = 82, file = "data_global/2d_data_ex_integral_vol.dat")
       open (unit = 83, file = "data_global/2d_data_ey_integral_vol.dat")
       open (unit = 84, file = "data_global/2d_data_ez_integral_vol.dat")
       open (unit = 85, file = "data_global/2d_data_exey_integral_vol.dat")
       open (unit = 86, file = "data_global/2d_data_exez_integral_vol.dat")
       open (unit = 87, file = "data_global/2d_data_eyez_integral_vol.dat")
       open (unit = 88, file = "data_global/2d_data_bx_integral_vol.dat")

       open (unit = 90, file = "data_global/2d_data_by_integral_vol.dat")
       open (unit = 91, file = "data_global/2d_data_bz_integral_vol.dat")
       open (unit = 92, file = "data_global/2d_data_bxby_integral_vol.dat")
       open (unit = 93, file = "data_global/2d_data_bxbz_integral_vol.dat")
       open (unit = 94, file = "data_global/2d_data_bybz_integral_vol.dat")
       open (unit = 95, file = "data_global/2d_data_v2_integral_vol.dat")
       open (unit = 96, file = "data_global/2d_data_b2_integral_vol.dat")
       open (unit = 97, file = "data_global/2d_data_vx_integral_vol.dat")
       open (unit = 98, file = "data_global/2d_data_vy_integral_vol.dat")
       open (unit = 99, file = "data_global/2d_data_vz_integral_vol.dat")
       
       open (unit =200, file = "data_global/2d_data_vxvy_integral_vol.dat")
       open (unit =201, file = "data_global/2d_data_vxvz_integral_vol.dat")
       open (unit =202, file = "data_global/2d_data_vyvz_integral_vol.dat")
       open (unit =203, file = "data_global/2d_data_q_integral_vol.dat")
       open (unit =204, file = "data_global/2d_data_p_integral_vol.dat")
       open (unit =205, file = "data_global/2d_data_rho_integral_vol.dat")
       open (unit =206, file = "data_global/2d_data_ekin_c_integral_vol.dat")
       open (unit =207, file = "data_global/2d_data_ekin_r_integral_vol.dat")
       open (unit =208, file = "data_global/2d_data_etot_integral_vol.dat")
       open (unit =209, file = "data_global/2d_data_efluid_r_integral_vol.dat")
       open (unit =210, file = "data_global/2d_data_emag_integral_vol.dat")
       open (unit =211, file = "data_global/2d_data_eint_c_integral_vol.dat")
       open (unit =212, file = "data_global/2d_data_eint_r_integral_vol.dat")
       open (unit =213, file = "data_global/2d_data_lorentz_integral_vol.dat")
       open (unit =214, file = "data_global/2d_data_D_integral_vol.dat")
       open (unit =215, file = "data_global/2d_data_D_zona_media_integral_vol.dat")

       open (unit =216, file = "data_global/2d_data_bx_max_val.dat")
       open (unit =217, file = "data_global/2d_data_phi_max_val.dat")
       
    else 

       write(*,*) "STOP subroutine openfile"
       write(*,*) "This DIM is not implemented yet"
       stop

    end if

    open (unit =301, file = "fix-up_densidad.dat")
    open (unit =302, file = "fix-up_presion.dat")
    open (unit =303, file = "fix-up_velocity_drift.dat")  
    
  end subroutine openfiles
