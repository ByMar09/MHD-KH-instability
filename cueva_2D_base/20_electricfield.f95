
  subroutine electricfield

    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

!$OMP PARALLEL PRIVATE(sigma_bar,Edotv,vrotB_x,vrotB_y,vrotB_z,V2int,sigma_dos_bar)
    
    if (MIRK == 1) then

!$OMP DO 

       do j=0+fix,jmax
          do i=0,imax

              if (SLC == 1) sigma = sigma_loc(i,j,1)

    sigma_bar = Delt * sigma * W(i,j,1)

    Edotv = Ex(i,j,1) * Vx(i,j,1) +                      &
            Ey(i,j,1) * Vy(i,j,1) +                      &
            Ez(i,j,1) * Vz(i,j,1)  

! vrotB for V(t_n) and B(t_n) see subroutine 19_varfin.f95
    
    vrotB_n_x = (Bzn(i,j,1)*Vy(i,j,1)-Byn(i,j,1)*Vz(i,j,1))
    vrotB_n_y = (Bxn(i,j,1)*Vz(i,j,1)-Bzn(i,j,1)*Vx(i,j,1))
    vrotB_n_z = (Byn(i,j,1)*Vx(i,j,1)-Bxn(i,j,1)*Vy(i,j,1))

! vrotB for V(t_n) and B(t_{n+1}) see subroutine 19_varfin.f95   

    vrotB_x = (Bz(i,j,1)*Vy(i,j,1)-By(i,j,1)*Vz(i,j,1))
    vrotB_y = (Bx(i,j,1)*Vz(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
    vrotB_z = (By(i,j,1)*Vx(i,j,1)-Bx(i,j,1)*Vy(i,j,1))
    

    Ex(i,j,1)  = 1.d0 /(1.d0 + sigma_bar) * ( Ex(i,j,1) + sigma_bar * Vx(i,j,1)     * Edotv              &
                                          -   sigma_bar * ( cm_2 * vrotB_n_x + (1.d0 - cm_2) * vrotB_x ) &
                                          +   Delt      *  Exfinsum(i,j,1) )

    Ey(i,j,1)  = 1.d0 /(1.d0 + sigma_bar) * ( Ey(i,j,1) + sigma_bar * Vy(i,j,1)     * Edotv              &
                                          -   sigma_bar * ( cm_2 * vrotB_n_y + (1.d0 - cm_2) * vrotB_y ) &
                                          +   Delt      *  Eyfinsum(i,j,1) )

    Ez(i,j,1)  = 1.d0 /(1.d0 + sigma_bar) * ( Ez(i,j,1) + sigma_bar * Vz(i,j,1)     * Edotv              &
                                          -   sigma_bar * ( cm_2 * vrotB_n_z + (1.d0 - cm_2) * vrotB_z ) &
                                          +   Delt      *  Ezfinsum(i,j,1) )

      end do !i
   end do !j

!$OMP END DO

    else if (MIRK == 2) then

!$OMP DO 

       do j=0+fix,jmax
          do i=0,imax

!     Velocity module

       V2int = Vxint(i,j,1,lmax)**2 +Vyint(i,j,1,lmax)**2 +Vzint(i,j,1,lmax)**2 
       
!     Lorentz factor

       Wint(i,j,1,lmax)  = 1.d0/sqrt(1-V2int)

!    Spatially Localized Conductivity

       if (SLC == 1) sigma = sigma_loc(i,j,1)


       sigma_dos_bar = Delt * sigma * Wint(i,j,1,lmax)

              

       !normal expression for MIRK2
       Edotv = Exint(i,j,1,lmax) * Vxint(i,j,1,lmax) +                      &
               Eyint(i,j,1,lmax) * Vyint(i,j,1,lmax) +                      &
               Ezint(i,j,1,lmax) * Vzint(i,j,1,lmax)

       !ussing c3=c6=0

!!$       Edotv = Exn(i,j,1) * Vxint(i,j,1,lmax) +                      &
!!$               Eyn(i,j,1) * Vyint(i,j,1,lmax) +                      &
!!$               Ezn(i,j,1) * Vzint(i,j,1,lmax)  
    

    vrotB_x = (Bzint(i,j,1,lmax)*Vyint(i,j,1,lmax)-Byint(i,j,1,lmax)*Vzint(i,j,1,lmax))
    vrotB_y = (Bxint(i,j,1,lmax)*Vzint(i,j,1,lmax)-Bzint(i,j,1,lmax)*Vxint(i,j,1,lmax))
    vrotB_z = (Byint(i,j,1,lmax)*Vxint(i,j,1,lmax)-Bxint(i,j,1,lmax)*Vyint(i,j,1,lmax))
!!$
!!$
!!$    vrotB_x = (Bzn(i,j,1)*Vy(i,j,1)-Byn(i,j,1)*Vz(i,j,1))
!!$    vrotB_y = (Bxn(i,j,1)*Vz(i,j,1)-Bzn(i,j,1)*Vx(i,j,1))
!!$    vrotB_z = (Byn(i,j,1)*Vx(i,j,1)-Bxn(i,j,1)*Vy(i,j,1))
    

    ! Old expressions for MIRK2
!!$    Ex(i,j,1)  = 0.5d0 / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar) * (                                            &
!!$                Exn(i,j,1) * ( 1.d0 - sigma_dos_bar * (1.d0 - cm_1) ) +                                                &
!!$                Exint(i,j,1,lmax) * (1.d0 - 2.d0 * sigma_dos_bar * cm_4) + Vxint(i,j,1,lmax) * sigma_dos_bar * Edotv - &
!!$                sigma_dos_bar * vrotB_x + Delt * Exfinsum(i,j,1)  )
!!$
!!$    Ey(i,j,1)  = 0.5d0 / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar) * (                                            &
!!$                Eyn(i,j,1) * ( 1.d0 - sigma_dos_bar * (1.d0 - cm_1) ) +                                                &
!!$                Eyint(i,j,1,lmax) * (1.d0 - 2.d0 * sigma_dos_bar * cm_4) + Vyint(i,j,1,lmax) * sigma_dos_bar * Edotv - &
!!$                sigma_dos_bar * vrotB_y + Delt * Eyfinsum(i,j,1)  )
!!$
!!$    Ez(i,j,1)  = 0.5d0 / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar) * (                                            &
!!$                Ezn(i,j,1) * ( 1.d0 - sigma_dos_bar * (1.d0 - cm_1) ) +                                                &
!!$                Ezint(i,j,1,lmax) * (1.d0 - 2.d0 * sigma_dos_bar * cm_4) + Vzint(i,j,1,lmax) * sigma_dos_bar * Edotv - &
!!$                sigma_dos_bar * vrotB_z + Delt * Ezfinsum(i,j,1)  )


    !NEW expressions for MIRK2
    Ex(i,j,1)  = Exint(i,j,1,lmax)                                                                                   &
               + 0.5d0 * (-1.d0 + sigma_dos_bar * (1.d0 - cm_1) ) / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar)   &
               *(Exint(i,j,1,lmax) - Exn(i,j,1))                                                                     &
               + 0.5d0 / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar)                                              &
               *(-sigma_dos_bar * Exint(i,j,1,lmax) + sigma_dos_bar * Edotv * Vxint(i,j,1,lmax)                      &
               + Delt * Exfinsum(i,j,1) - sigma_dos_bar * vrotB_x )


    Ey(i,j,1)  = Eyint(i,j,1,lmax)                                                                                   &
               + 0.5d0 * (-1.d0 + sigma_dos_bar * (1.d0 - cm_1) ) / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar)   &
               *(Eyint(i,j,1,lmax) - Eyn(i,j,1))                                                                     &
               + 0.5d0 / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar)                                              &
               *(-sigma_dos_bar * Eyint(i,j,1,lmax) + sigma_dos_bar * Edotv * Vyint(i,j,1,lmax)                      &
               + Delt * Eyfinsum(i,j,1) - sigma_dos_bar * vrotB_y )

    Ez(i,j,1)  = Ezint(i,j,1,lmax)                                                                                   &
               + 0.5d0 * (-1.d0 + sigma_dos_bar * (1.d0 - cm_1) ) / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar)   &
               *(Ezint(i,j,1,lmax) - Ezn(i,j,1))                                                                     &
               + 0.5d0 / (1.d0 + (0.5d0 * cm_1 - cm_4) * sigma_dos_bar)                                              &
               *(-sigma_dos_bar * Ezint(i,j,1,lmax) + sigma_dos_bar * Edotv * Vzint(i,j,1,lmax)                      &
               + Delt * Ezfinsum(i,j,1) - sigma_dos_bar * vrotB_z )

    



         end do ! for i
      end do ! for j

!$OMP END DO

    else if (MIRK == 0) then

!$OMP DO 

       do j=0+fix,jmax
          do i=0,imax

    Ex(i,j,1)  = Ex(i,j,1)   + Delt *  Exfinsum(i,j,1)
    Ey(i,j,1)  = Ey(i,j,1)   + Delt *  Eyfinsum(i,j,1)
    Ez(i,j,1)  = Ez(i,j,1)   + Delt *  Ezfinsum(i,j,1)

        end do !i
     end do !j

!$OMP END DO

     else

       write(*,*) "These MIRK scheme is not implemented yet"

       stop

    end if

!$OMP END PARALLEL

  end subroutine electricfield
