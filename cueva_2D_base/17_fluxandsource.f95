  !     *******************************************************************
  !     Subroutine which construct the intermediate variables flows and sources
  !     *******************************************************************


 subroutine fluxandsource

    use scalar
    use parameters
    use threevectors
    use fourvectors
    use funciones, only: mcl

    implicit none


!$OMP   PARALLEL PRIVATE(psiint_L,psiint_R,phiint_L,phiint_R,qint_L,qint_R),   &
!$OMP & PRIVATE(Edotv,vrotB_x,vrotB_y,vrotB_z,Vx_elec,Vy_elec,Vz_elec,W_elec), &   
!$OMP & PRIVATE(Ex_sour_rec_L,Ex_sour_rec_R,Ey_sour_rec_L,Ey_sour_rec_R,Ez_sour_rec_L,Ez_sour_rec_R),             &
!$OMP & PRIVATE(Ex_sour_rec_x_L,Ex_sour_rec_x_R,Ey_sour_rec_x_L,Ey_sour_rec_x_R,Ez_sour_rec_x_L,Ez_sour_rec_x_R), &
!$OMP & PRIVATE(Ex_sour_rec_y_L,Ex_sour_rec_y_R,Ey_sour_rec_y_L,Ey_sour_rec_y_R,Ez_sour_rec_y_L,Ez_sour_rec_y_R)

    
    if (DIM == 1 ) then

       call boundary_conserved
       call boundary_electric

       if ( var_source_rec == 1 ) then

!$OMP DO 

    do i=-1,imax

       !---------------------------------
       !  Reconstruction intermediate variables ---> 1/2 (U_L + U_R)
       !---------------------------------

    !     Left - Right  intermediate values of augmented fields

    psiint_L = psiint(i  ,1,1,l) +  FAC_SR  * mcl(psiint(i+1,1,1,l) - psiint(i  ,1,1,l) ,  &
                                                  psiint(i  ,1,1,l) - psiint(i-1,1,1,l))
    psiint_R = psiint(i+1,1,1,l) -  FAC_SR  * mcl(psiint(i+2,1,1,l) - psiint(i+1,1,1,l) ,  &
                                                  psiint(i+1,1,1,l) - psiint(i  ,1,1,l))

    phiint_L = phiint(i  ,1,1,l) +  FAC_SR  * mcl(phiint(i+1,1,1,l) - phiint(i  ,1,1,l) ,  &
                                                  phiint(i  ,1,1,l) - phiint(i-1,1,1,l))
    phiint_R = phiint(i+1,1,1,l) -  FAC_SR  * mcl(phiint(i+2,1,1,l) - phiint(i+1,1,1,l) ,  &
                                                  phiint(i+1,1,1,l) - phiint(i  ,1,1,l))


    psiint(i  ,1,1,l) = 0.5d0 * (psiint_L + psiint_R)
    phiint(i  ,1,1,l) = 0.5d0 * (phiint_L + phiint_R)

    !     Left - Right intermediate values of magnetic field

!!$    Bxint_L =  Bxint(i  ,1,1,l) +  FAC_SR  * mcl(Bxint(i+1,1,1,l) - Bxint(i  ,1,1,l) ,  &
!!$                                              Bxint(i  ,1,1,l) - Bxint(i-1,1,1,l)) 
!!$    Bxint_R =  Bxint(i+1,1,1,l) -  FAC_SR  * mcl(Bxint(i+2,1,1,l) - Bxint(i+1,1,1,l) ,  &
!!$                                              Bxint(i+1,1,1,l) - Bxint(i  ,1,1,l)) 
!!$
!!$    Byint_L =  Byint(i  ,1,1,l) +  FAC_SR  * mcl(Byint(i+1,1,1,l) - Byint(i  ,1,1,l) ,  &
!!$                                              Byint(i  ,1,1,l) - Byint(i-1,1,1,l))  
!!$    Byint_R =  Byint(i+1,1,1,l) -  FAC_SR  * mcl(Byint(i+2,1,1,l) - Byint(i+1,1,1,l) ,  &
!!$                                              Byint(i+1,1,1,l) - Byint(i  ,1,1,l))  
!!$
!!$    Bzint_L =  Bzint(i  ,1,1,l) +  FAC_SR  * mcl(Bzint(i+1,1,1,l) - Bzint(i  ,1,1,l) ,  &
!!$                                              Bzint(i  ,1,1,l) - Bzint(i-1,1,1,l))  
!!$    Bzint_R =  Bzint(i+1,1,1,l) -  FAC_SR  * mcl(Bzint(i+2,1,1,l) - Bzint(i+1,1,1,l) ,  &
!!$                                              Bzint(i+1,1,1,l) - Bzint(i  ,1,1,l))  
!!$
!!$    Bxint(i  ,1,1,l) = 0.5d0 * (Bxint_L + Bxint_R)
!!$    Byint(i  ,1,1,l) = 0.5d0 * (Byint_L + Byint_R)
!!$    Bzint(i  ,1,1,l) = 0.5d0 * (Bzint_L + Bzint_R)


    !-----------------------------------------------------------------------------------

    !     Left - Right intermediate value of electric field

!!$    Exint_L =  Exint(i  ,1,1,l)  +  FAC_SR  * mcl(Exint(i+1,1,1,l) - Exint(i  ,1,1,l) ,  &
!!$                                               Exint(i  ,1,1,l) - Exint(i-1,1,1,l)) 
!!$    Exint_R =  Exint(i+1,1,1,l)  -  FAC_SR  * mcl(Exint(i+2,1,1,l) - Exint(i+1,1,1,l) ,  &
!!$                                               Exint(i+1,1,1,l) - Exint(i  ,1,1,l)) 
!!$
!!$    Eyint_L =  Eyint(i  ,1,1,l)  +  FAC_SR  * mcl(Eyint(i+1,1,1,l) - Eyint(i  ,1,1,l) ,  &
!!$                                               Eyint(i  ,1,1,l) - Eyint(i-1,1,1,l)) 
!!$    Eyint_R =  Eyint(i+1,1,1,l)  -  FAC_SR  * mcl(Eyint(i+2,1,1,l) - Eyint(i+1,1,1,l) ,  &
!!$                                               Eyint(i+1,1,1,l) - Eyint(i  ,1,1,l)) 
!!$
!!$
!!$    Ezint_L =  Ezint(i  ,1,1,l)  +  FAC_SR  * mcl(Ezint(i+1,1,1,l) - Ezint(i  ,1,1,l) ,  &
!!$                                               Ezint(i  ,1,1,l) - Ezint(i-1,1,1,l)) 
!!$    Ezint_R =  Ezint(i+1,1,1,l)  -  FAC_SR  * mcl(Ezint(i+2,1,1,l) - Ezint(i+1,1,1,l) ,  &
!!$                                               Ezint(i+1,1,1,l) - Ezint(i  ,1,1,l)) 
!!$
!!$    Exint(i  ,1,1,l) = 0.5d0 * (Exint_L + Exint_R)
!!$    Eyint(i  ,1,1,l) = 0.5d0 * (Eyint_L + Eyint_R)
!!$    Ezint(i  ,1,1,l) = 0.5d0 * (Ezint_L + Ezint_R)


    !-----------------------------------------------------------------------------------

    !     Left - Right Intermediate value of charge density 

    qint_L = qint(i  ,1,1,l)  +  FAC_SR  * mcl(qint(i+1,1,1,l) - qint(i  ,1,1,l) ,  &
                                               qint(i  ,1,1,l) - qint(i-1,1,1,l)) 
    qint_R = qint(i+1,1,1,l)  -  FAC_SR  * mcl(qint(i+2,1,1,l) - qint(i+1,1,1,l) ,  &
                                               qint(i+1,1,1,l) - qint(i  ,1,1,l)) 

    qint(i  ,1,1,l) = 0.5d0 * (qint_L + qint_R)

    !-----------------------------------------------------------------------------------

    !     Left - Right Intermediate value of Velocity

!!$          if (REC_PRIM == 0) then
!!$
!!$             Vx_L = Vx(i  ,1,1) +  FAC_SR  * mcl(Vx(i+1,1,1) - Vx(i  ,1,1) ,  &
!!$                                              Vx(i  ,1,1) - Vx(i-1,1,1))
!!$             Vx_R = Vx(i+1,1,1) -  FAC_SR  * mcl(Vx(i+2,1,1) - Vx(i+1,1,1) ,  &
!!$                                              Vx(i+1,1,1) - Vx(i  ,1,1))
!!$
!!$             Vy_L = Vy(i  ,1,1) +  FAC_SR  * mcl(Vy(i+1,1,1) - Vy(i  ,1,1) ,  &
!!$                                              Vy(i  ,1,1) - Vy(i-1,1,1))
!!$             Vy_R = Vy(i+1,1,1) -  FAC_SR  * mcl(Vy(i+2,1,1) - Vy(i+1,1,1) ,  &
!!$                                              Vy(i+1,1,1) - Vy(i  ,1,1))
!!$
!!$             Vz_L = Vz(i  ,1,1) +  FAC_SR  * mcl(Vz(i+1,1,1) - Vz(i  ,1,1) ,  &
!!$                                              Vz(i  ,1,1) - Vz(i-1,1,1))
!!$             Vz_R = Vz(i+1,1,1) -  FAC_SR  * mcl(Vz(i+2,1,1) - Vz(i+1,1,1) ,  &
!!$                                              Vz(i+1,1,1) - Vz(i  ,1,1))
!!$
!!$
!!$             Vx(i  ,1,1) = 0.5d0 * (Vx_L + Vx_R)
!!$             Vy(i  ,1,1) = 0.5d0 * (Vy_L + Vy_R)
!!$             Vz(i  ,1,1) = 0.5d0 * (Vz_L + Vz_R)
!!$
!!$          else if (REC_PRIM == 1) then
!!$
!!$             Vx_L = Vxint(i  ,1,1,l) +  FAC_SR  * mcl(Vxint(i+1,1,1,l) - Vxint(i  ,1,1,l) ,  &
!!$                                                   Vxint(i  ,1,1,l) - Vxint(i-1,1,1,l))
!!$             Vx_R = Vxint(i+1,1,1,l) -  FAC_SR  * mcl(Vxint(i+2,1,1,l) - Vxint(i+1,1,1,l) ,  &
!!$                                                   Vxint(i+1,1,1,l) - Vxint(i  ,1,1,l))
!!$
!!$             Vy_L = Vyint(i  ,1,1,l) +  FAC_SR  * mcl(Vyint(i+1,1,1,l) - Vyint(i  ,1,1,l) ,  &
!!$                                                   Vyint(i  ,1,1,l) - Vyint(i-1,1,1,l))
!!$             Vy_R = Vyint(i+1,1,1,l) -  FAC_SR  * mcl(Vyint(i+2,1,1,l) - Vyint(i+1,1,1,l) ,  &
!!$                                                   Vyint(i+1,1,1,l) - Vyint(i  ,1,1,l))
!!$
!!$             Vz_L = Vzint(i  ,1,1,l) +  FAC_SR  * mcl(Vzint(i+1,1,1,l) - Vzint(i  ,1,1,l) ,  &
!!$                                                   Vzint(i  ,1,1,l) - Vzint(i-1,1,1,l))
!!$             Vy_R = Vzint(i+1,1,1,l) -  FAC_SR  * mcl(Vzint(i+2,1,1,l) - Vzint(i+1,1,1,l) ,  &
!!$                                                   Vzint(i+1,1,1,l) - Vzint(i  ,1,1,l))
!!$
!!$
!!$             Vxint(i  ,1,1,l) = 0.5d0 * (Vx_L + Vx_R)
!!$             Vyint(i  ,1,1,l) = 0.5d0 * (Vy_L + Vy_R)
!!$             Vzint(i  ,1,1,l) = 0.5d0 * (Vz_L + Vz_R)
!!$
!!$          else
!!$
!!$             write(*,*) "STOP: subroutine fluxandsource"
!!$             write(*,*) "REC_PRIM parameter is not valid"
!!$             stop
!!$
!!$          end if


        end do

!$OMP END DO

     end if


!$OMP DO 

             do i=0,imax


          if (REC_PRIM == 0) then

             Edotv = Exint(i,1,1,l) * Vx(i,1,1) +                      &
                     Eyint(i,1,1,l) * Vy(i,1,1) +                      &
                     Ezint(i,1,1,l) * Vz(i,1,1)  

             vrotB_x = (Bzint(i,1,1,l)*Vy(i,1,1)-Byint(i,1,1,l)*Vz(i,1,1))
             vrotB_y = (Bxint(i,1,1,l)*Vz(i,1,1)-Bzint(i,1,1,l)*Vx(i,1,1))
             vrotB_z = (Byint(i,1,1,l)*Vx(i,1,1)-Bxint(i,1,1,l)*Vy(i,1,1))

!!$             vrotB_x = (Bz(i,1,1)*Vy(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$             vrotB_y = (Bx(i,1,1)*Vz(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$             vrotB_z = (By(i,1,1)*Vx(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$
!!$
!!$      vrotB_x_mp5_L(i,1) = (Bzint_mp5_L(i,1) * Vy_mp5_L(i,1) - Byint_mp5_L(i,1) * Vz_mp5_L(i,1))
!!$      vrotB_y_mp5_L(i,1) = (Bxint_mp5_L(i,1) * Vz_mp5_L(i,1) - Bzint_mp5_L(i,1) * Vx_mp5_L(i,1))
!!$      vrotB_z_mp5_L(i,1) = (Byint_mp5_L(i,1) * Vx_mp5_L(i,1) - Bxint_mp5_L(i,1) * Vy_mp5_L(i,1))
!!$
!!$      vrotB_x_mp5_R(i,1) = (Bzint_mp5_R(i,1) * Vy_mp5_R(i,1) - Byint_mp5_R(i,1) * Vz_mp5_R(i,1))
!!$      vrotB_y_mp5_R(i,1) = (Bxint_mp5_R(i,1) * Vz_mp5_R(i,1) - Bzint_mp5_R(i,1) * Vx_mp5_R(i,1))
!!$      vrotB_z_mp5_R(i,1) = (Byint_mp5_R(i,1) * Vx_mp5_R(i,1) - Bxint_mp5_R(i,1) * Vy_mp5_R(i,1))
!!$
!!$
!!$      vrotB_x = 0.5d0 * (vrotB_x_mp5_L(i,1) + vrotB_x_mp5_R(i,1))
!!$      vrotB_y = 0.5d0 * (vrotB_y_mp5_L(i,1) + vrotB_y_mp5_R(i,1))
!!$      vrotB_z = 0.5d0 * (vrotB_z_mp5_L(i,1) + vrotB_z_mp5_R(i,1))

             Vx_elec = Vx(i,1,1)
             Vy_elec = Vy(i,1,1)
             Vz_elec = Vz(i,1,1)
    
             W_elec  = 1.d0/sqrt(1.d0 - Vx_elec**2 -  Vy_elec**2 - Vz_elec**2)

          else if (REC_PRIM == 1) then

             Edotv = Exint(i,1,1,l) * Vxint(i,1,1,l) +                      &
                     Eyint(i,1,1,l) * Vyint(i,1,1,l) +                      &
                     Ezint(i,1,1,l) * Vzint(i,1,1,l)  

             vrotB_x = (Bzint(i,1,1,l)*Vyint(i,1,1,l)-Byint(i,1,1,l)*Vzint(i,1,1,l))
             vrotB_y = (Bxint(i,1,1,l)*Vzint(i,1,1,l)-Bzint(i,1,1,l)*Vxint(i,1,1,l))
             vrotB_z = (Byint(i,1,1,l)*Vxint(i,1,1,l)-Bxint(i,1,1,l)*Vyint(i,1,1,l))

             Vx_elec = Vxint(i,1,1,l)
             Vy_elec = Vyint(i,1,1,l)
             Vz_elec = Vzint(i,1,1,l)

              W_elec  = 1.d0/sqrt(1.d0 - Vx_elec**2 -  Vy_elec**2 - Vz_elec**2)

          else

             write(*,*) "STOP: subroutine fluxandsource"
             write(*,*) "REC_PRIM parameter is not valid"
             stop

          end if


              if (SLC == 1) sigma = sigma_loc(i,1,1)
                       

    psiflux(i,1,1,l) = -(Exlaxpsi(i  ,1,1,l) - Exlaxpsi(i-1,1  ,1  ,l))/Delx -   &  
                        (Eylaxpsi(i  ,1,1,l) - Eylaxpsi(i  ,1  ,1  ,l))/Dely -   &
                        (Ezlaxpsi(i  ,1,1,l) - Ezlaxpsi(i  ,1  ,1  ,l))/Delz +   &
                         qint(i,1,1,l) - kappa_psi * psiint(i,1,1,l) 

    phiflux(i,1,1,l) =  -(Bxlaxphi(i  ,1,1,l) - Bxlaxphi(i-1,1  ,1  ,l))/ Delx -  &
                         (Bylaxphi(i  ,1,1,l) - Bylaxphi(i  ,1  ,1  ,l))/Dely  -  &
                         (Bzlaxphi(i  ,1,1,l) - Bzlaxphi(i  ,1  ,1  ,l))/Delz  -  &
                          kappa_phi * phiint(i,1,1,l) 

    Bxflux(i,1,1,l)  =  -(EzlaxBx (i  ,1,1,l) - EzlaxBx (i  ,1  ,1  ,l))/ Dely +  & 
                         (EylaxBx (i  ,1,1,l) - EylaxBx (i  ,1  ,1  ,l))/ Delz -  &
                         (philaxBx(i  ,1,1,l) - philaxBx(i-1,1  ,1  ,l))/ Delx 

    Byflux(i,1,1,l)  =  -(ExlaxBy (i  ,1,1,l) - ExlaxBy (i  ,1  ,1  ,l))/ Delz +  &
                         (EzlaxBy (i  ,1,1,l) - EzlaxBy (i-1,1  ,1  ,l))/ Delx -  &
                         (philaxBy(i  ,1,1,l) - philaxBy(i  ,1  ,1  ,l))/ Dely 

    Bzflux(i,1,1,l)  =  -(EylaxBz (i  ,1,1,l) - EylaxBz (i-1,1  ,1  ,l))/ Delx +  &
                         (ExlaxBz (i  ,1,1,l) - ExlaxBz (i  ,1  ,1  ,l))/ Dely -  &
                         (philaxBz(i  ,1,1,l) - philaxBz(i  ,1  ,1  ,l))/ Delz



    Exastflux(i,1,1,l)  =  (BzlaxEx (i  ,1,1,l) - BzlaxEx (i  ,1  ,1  ,l))/ Dely -  &
                           (BylaxEx (i  ,1,1,l) - BylaxEx (i  ,1  ,1  ,l))/ Delz -  &
                           (psilaxEx(i  ,1,1,l) - psilaxEx(i-1,1  ,1  ,l))/ Delx -  &
                            qint(i,1,1,l)       * Vx_elec       

    Eyastflux(i,1,1,l)  =  (BxlaxEy (i  ,1,1,l) - BxlaxEy (i  ,1  ,1  ,l))/ Delz -  &
                           (BzlaxEy (i  ,1,1,l) - BzlaxEy (i-1,1  ,1  ,l))/ Delx -  &
                           (psilaxEy(i  ,1,1,l) - psilaxEy(i  ,1  ,1  ,l))/ Dely -  &
                            qint(i,1,1,l)       * Vy_elec 

    Ezastflux(i,1,1,l)  =  (BylaxEz (i  ,1,1,l) - BylaxEz (i-1,1  ,1  ,l))/ Delx -  &
                           (BxlaxEz (i  ,1,1,l) - BxlaxEz (i  ,1  ,1  ,l))/ Dely -  &
                           (psilaxEz(i  ,1,1,l) - psilaxEz(i  ,1  ,1  ,l))/ Delz -  &
                            qint(i,1,1,l)       * Vz_elec 


        if (source_rec_mpx == 0) then

    Exastsour(i,1,1,l)  =   sigma               * W_elec  *                       &
                           (Exint(i,1,1,l)      + vrotB_x   - Edotv * Vx_elec )

    Eyastsour(i,1,1,l)  =   sigma               * W_elec   *                      &
                           (Eyint(i,1,1,l)      + vrotB_y    - Edotv * Vy_elec )
    
    Ezastsour(i,1,1,l)  =   sigma               * W_elec   *                      &
                           (Ezint(i,1,1,l)      + vrotB_z    - Edotv * Vz_elec )

    else if (flux_solver == 5 .or. flux_solver == 6  .or. flux_solver == 7 .and. source_rec_mpx == 1) then


       Exastsour(i,1,1,l)  =  Exastsour(i,1,1,l) 
       Eyastsour(i,1,1,l)  =  Eyastsour(i,1,1,l)
       Ezastsour(i,1,1,l)  =  Ezastsour(i,1,1,l)

    else

       print*, "Revise the combination between flux_solver and source_rec_mpx parameters"
       print*, "It seems these two parameters are not valid"
       print*, "STOP subroutine 17_fluxandsource.f95"
       stop

    end if

  
    qflux(i,1,1,l)      = -(Jxlax(i,1,1,l) - Jxlax(i-1,1  ,1  ,l))/Delx -           & !ojo i-1
                           (Jylax(i,1,1,l) - Jylax(i  ,1  ,1  ,l))/Dely -           &
                           (Jzlax(i,1,1,l) - Jzlax(i  ,1  ,1  ,l))/Delz

!    print*,  i, Jxlax(i,1,1,l), Jxlax(i-1,1  ,1  ,l)

    Dflux(i,1,1,l)      = -(FDxlax(i  ,1,1,l) - FDxlax(i-1,1  ,1  ,l))/Delx -       &
                           (FDylax(i  ,1,1,l) - FDylax(i  ,1  ,1  ,l))/Dely -       &
                           (FDzlax(i  ,1,1,l) - FDzlax(i  ,1  ,1  ,l))/Delz 

    tauflux(i,1,1,l)    = -(Ftauxlax(i  ,1,1,l) - Ftauxlax(i-1,1  ,1  ,l))/Delx -   & 
                           (Ftauylax(i  ,1,1,l) - Ftauylax(i  ,1  ,1  ,l))/Dely -   &
                           (Ftauzlax(i  ,1,1,l) - Ftauzlax(i  ,1  ,1  ,l))/Delz 

! -------------------------------------- EGLM--------------------------------------

    taueglm(i,1,1,l)    = -(Exint(i,1,1,l)          *                                  & 
                           (psitauxlax(i  ,1  ,1,l) - psitauxlax(i-1,1  ,1,l))/Delx +  &
                            Eyint(i,1,1,l)          *                                  & 
                           (psitauylax(i  ,1  ,1,l) - psitauylax(i  ,1  ,1,l))/Dely +  &
                            Ezint(i,1,1,l)          *                                  & 
                           (psitauzlax(i  ,1  ,1,l) - psitauzlax(i  ,1  ,1,l))/Delz )  &
                          -(Bxint(i,1,1,l)          *                                  & 
                           (phitauxlax(i  ,1  ,1,l) - phitauxlax(i-1,1  ,1,l))/Delx +  &
                            Byint(i,1,1,l)          *                                  & 
                           (phitauylax(i  ,1  ,1,l) - phitauylax(i  ,1  ,1,l))/Dely +  &
                            Bzint(i,1,1,l)          *                                  & 
                           (phitauzlax(i  ,1  ,1,l) - phitauzlax(i  ,1  ,1,l))/Delz )  

! -------------------------------------- EGLM--------------------------------------

    Sxflux(i,1,1,l)     = -(FSxxlax(i  ,1,1,l) - FSxxlax(i-1,1  ,1  ,l))/Delx -     &
                           (FSyxlax(i  ,1,1,l) - FSyxlax(i  ,1  ,1  ,l))/Dely -     &
                           (FSzxlax(i  ,1,1,l) - FSzxlax(i  ,1  ,1  ,l))/Delz 

    Syflux(i,1,1,l)     = -(FSxylax(i  ,1,1,l) - FSxylax(i-1,1  ,1  ,l))/Delx -     &
                           (FSyylax(i  ,1,1,l) - FSyylax(i  ,1  ,1  ,l))/Dely -     &
                           (FSzylax(i  ,1,1,l) - FSzylax(i  ,1  ,1  ,l))/Delz 

    Szflux(i,1,1,l)     = -(FSxzlax(i  ,1,1,l) - FSxzlax(i-1,1  ,1  ,l))/Delx -     &
                           (FSyzlax(i  ,1,1,l) - FSyzlax(i  ,1  ,1  ,l))/Dely -     &
                           (FSzzlax(i  ,1,1,l) - FSzzlax(i  ,1  ,1  ,l))/Delz 

! -------------------------------------- EGLM--------------------------------------
    Sxeglm(i,1,1,l)     = -(Bxint(i,1,1,l)       * (                                  &
                           (BxSxlax(i  ,1  ,1,l) - BxSxlax (i-1,1  ,1  ,l))/Delx +    &
                           (BySxlax(i  ,1  ,1,l) - BySxlax (i  ,1  ,1  ,l))/Dely +    &
                           (BzSxlax(i  ,1  ,1,l) - BzSxlax (i  ,1  ,1  ,l))/Delz ))

!                           (Bxint(i  ,1  ,1,l) - Bxint (i-1,1  ,1  ,l))/Delx +    &
!                           (Byint(i  ,1  ,1,l) - Byint (i  ,1  ,1  ,l))/Dely +    &
!                           (Bzint(i  ,1  ,1,l) - Bzint (i  ,1  ,1  ,l))/Delz ))

    Syeglm(i,1,1,l)     = -(Byint(i,1,1,l)       * (                                  &
                           (BxSylax(i  ,1  ,1,l) - BxSylax (i-1,1  ,1  ,l))/Delx +    &
                           (BySylax(i  ,1  ,1,l) - BySylax (i  ,1  ,1  ,l))/Dely +    &
                           (BzSylax(i  ,1  ,1,l) - BzSylax (i  ,1  ,1  ,l))/Delz ))

!                           (Bxint(i  ,1  ,1,l) - Bxint (i-1,1  ,1  ,l))/Delx +    &
!                           (Byint(i  ,1  ,1,l) - Byint (i  ,1  ,1  ,l))/Dely +    &
!                           (Bzint(i  ,1  ,1,l) - Bzint (i  ,1  ,1  ,l))/Delz ))



    Szeglm(i,1,1,l)     = -(Bzint(i,1,1,l)       * (                                  &
                           (BxSzlax(i  ,1  ,1,l) - BxSzlax (i-1,1  ,1  ,l))/Delx +    &
                           (BySzlax(i  ,1  ,1,l) - BySzlax (i  ,1  ,1  ,l))/Dely +    &
                           (BzSzlax(i  ,1  ,1,l) - BzSzlax (i  ,1  ,1  ,l))/Delz ))

!                           (Bxint(i  ,1  ,1,l) - Bxint (i-1,1  ,1  ,l))/Delx +    &
!                           (Byint(i  ,1  ,1,l) - Byint (i  ,1-1,1  ,l))/Dely +    &
!                           (Bzint(i  ,1  ,1,l) - Bzint (i  ,1  ,1  ,l))/Delz ))



! -------------------------------------- EGLM--------------------------------------


   end do ! for i

!$OMP END DO



   if ( source_rec == 1) then

      call source_boundary

!$OMP DO 

      do i=0,imax


    Ex_sour_rec_L = Exastsour(i,1,1,l)  +  FAC_SR  * mcl(Exastsour(i+1,1,1,l) - Exastsour(i  ,1,1,l) ,  &
                                                         Exastsour(i  ,1,1,l) - Exastsour(i-1,1,1,l)) 
    Ex_sour_rec_R = Exastsour(i,1,1,l)  -  FAC_SR  * mcl(Exastsour(i+2,1,1,l) - Exastsour(i+1,1,1,l) ,  &
                                                         Exastsour(i+1,1,1,l) - Exastsour(i  ,1,1,l)) 

    Ey_sour_rec_L = Exastsour(i,1,1,l)  +  FAC_SR  * mcl(Eyastsour(i+1,1,1,l) - Eyastsour(i  ,1,1,l) ,  &
                                                         Eyastsour(i  ,1,1,l) - Eyastsour(i-1,1,1,l)) 
    Ey_sour_rec_R = Exastsour(i,1,1,l)  -  FAC_SR  * mcl(Eyastsour(i+2,1,1,l) - Eyastsour(i+1,1,1,l) ,  &
                                                         Eyastsour(i+1,1,1,l) - Eyastsour(i  ,1,1,l)) 

    Ez_sour_rec_L = Exastsour(i,1,1,l)  +  FAC_SR  * mcl(Ezastsour(i+1,1,1,l) - Ezastsour(i  ,1,1,l) ,  &
                                                         Ezastsour(i  ,1,1,l) - Ezastsour(i-1,1,1,l)) 
    Ez_sour_rec_R = Exastsour(i,1,1,l)  -  FAC_SR  * mcl(Ezastsour(i+2,1,1,l) - Ezastsour(i+1,1,1,l) ,  &
                                                         Ezastsour(i+1,1,1,l) - Ezastsour(i  ,1,1,l)) 


    Exastsour(i,1,1,l) = 0.5d0 * (Ex_sour_rec_L + Ex_sour_rec_R)
    Eyastsour(i,1,1,l) = 0.5d0 * (Ey_sour_rec_L + Ey_sour_rec_R)
    Ezastsour(i,1,1,l) = 0.5d0 * (Ez_sour_rec_L + Ez_sour_rec_R)


       end do

!$OMP END DO 

   end if


   else if (DIM == 2) then


       call boundary_conserved
       call boundary_electric

!$OMP DO 

       do j=0,jmax
          do i=0,imax

          if (REC_PRIM == 0) then

             Edotv = Exint(i,j,1,l) * Vx(i,j,1) +                      &
                     Eyint(i,j,1,l) * Vy(i,j,1) +                      &
                     Ezint(i,j,1,l) * Vz(i,j,1)  

             vrotB_x = (Bzint(i,j,1,l)*Vy(i,j,1)-Byint(i,j,1,l)*Vz(i,j,1))
             vrotB_y = (Bxint(i,j,1,l)*Vz(i,j,1)-Bzint(i,j,1,l)*Vx(i,j,1))
             vrotB_z = (Byint(i,j,1,l)*Vx(i,j,1)-Bxint(i,j,1,l)*Vy(i,j,1))

             Vx_elec = Vx(i,j,1)
             Vy_elec = Vy(i,j,1)
             Vz_elec = Vz(i,j,1)
    
             W_elec  = W(i,j,1)

          else if (REC_PRIM == 1) then

             Edotv = Exint(i,j,1,l) * Vxint(i,j,1,l) +                      &
                     Eyint(i,j,1,l) * Vyint(i,j,1,l) +                      &
                     Ezint(i,j,1,l) * Vzint(i,j,1,l)  

             vrotB_x = (Bzint(i,j,1,l)*Vyint(i,j,1,l)-Byint(i,j,1,l)*Vzint(i,j,1,l))
             vrotB_y = (Bxint(i,j,1,l)*Vzint(i,j,1,l)-Bzint(i,j,1,l)*Vxint(i,j,1,l))
             vrotB_z = (Byint(i,j,1,l)*Vxint(i,j,1,l)-Bxint(i,j,1,l)*Vyint(i,j,1,l))

             Vx_elec = Vxint(i,j,1,l)
             Vy_elec = Vyint(i,j,1,l)
             Vz_elec = Vzint(i,j,1,l)

             W_elec  = Wint(i,j,1,l)

          else

             write(*,*) "STOP: subroutine fluxandsource"
             write(*,*) "REC_PRIM parameter is not valid"
             stop

          end if


              if (SLC == 1) sigma = sigma_loc(i,j,1)


    psiflux(i,j,1,l) = -(Exlaxpsi(i  ,j  ,1,l) - Exlaxpsi(i-1,j  ,1  ,l))/Delx -   &
                        (Eylaxpsi(i  ,j  ,1,l) - Eylaxpsi(i  ,j-1,1  ,l))/Dely -   &
                        (Ezlaxpsi(i  ,j  ,1,l) - Ezlaxpsi(i  ,j  ,1  ,l))/Delz +   &
                         qint(i,j,1,l) - kappa_psi * psiint(i,j,1,l) 

    phiflux(i,j,1,l) =  -(Bxlaxphi(i  ,j  ,1,l) - Bxlaxphi(i-1,j  ,1  ,l))/Delx  -  &
                         (Bylaxphi(i  ,j  ,1,l) - Bylaxphi(i  ,j-1,1  ,l))/Dely  -  &
                         (Bzlaxphi(i  ,j  ,1,l) - Bzlaxphi(i  ,j  ,1  ,l))/Delz  -  &
                          kappa_phi * phiint(i,j,1,l) 

    Bxflux(i,j,1,l)  =  -(EzlaxBx (i  ,j  ,1,l) - EzlaxBx (i  ,j-1,1  ,l))/ Dely +  & 
                         (EylaxBx (i  ,j  ,1,l) - EylaxBx (i  ,j  ,1  ,l))/ Delz -  &
                         (philaxBx(i  ,j  ,1,l) - philaxBx(i-1,j  ,1  ,l))/ Delx 

    Byflux(i,j,1,l)  =  -(ExlaxBy (i  ,j  ,1,l) - ExlaxBy (i  ,j  ,1  ,l))/ Delz +  &
                         (EzlaxBy (i  ,j  ,1,l) - EzlaxBy (i-1,j  ,1  ,l))/ Delx -  &
                         (philaxBy(i  ,j  ,1,l) - philaxBy(i  ,j-1,1  ,l))/ Dely 

    Bzflux(i,j,1,l)  =  -(EylaxBz (i  ,j  ,1,l) - EylaxBz (i-1,j  ,1  ,l))/ Delx +  &
                         (ExlaxBz (i  ,j  ,1,l) - ExlaxBz (i  ,j-1,1  ,l))/ Dely -  &
                         (philaxBz(i  ,j  ,1,l) - philaxBz(i  ,j  ,1  ,l))/ Delz


    Exastflux(i,j,1,l)  =  (BzlaxEx (i  ,j  ,1,l) - BzlaxEx (i  ,j-1,1  ,l))/ Dely -  &
                           (BylaxEx (i  ,j  ,1,l) - BylaxEx (i  ,j  ,1  ,l))/ Delz -  &
                           (psilaxEx(i  ,j  ,1,l) - psilaxEx(i-1,j  ,1  ,l))/ Delx -  &
                            qint(i,j,1,l)       * Vx_elec       

    Eyastflux(i,j,1,l)  =  (BxlaxEy (i  ,j  ,1,l) - BxlaxEy (i  ,j  ,1  ,l))/ Delz -  &
                           (BzlaxEy (i  ,j  ,1,l) - BzlaxEy (i-1,j  ,1  ,l))/ Delx -  &
                           (psilaxEy(i  ,j  ,1,l) - psilaxEy(i  ,j-1,1  ,l))/ Dely -  &
                            qint(i,j,1,l)       * Vy_elec 

    Ezastflux(i,j,1,l)  =  (BylaxEz (i  ,j  ,1,l) - BylaxEz (i-1,j  ,1  ,l))/ Delx -  &
                           (BxlaxEz (i  ,j  ,1,l) - BxlaxEz (i  ,j-1,1  ,l))/ Dely -  &
                           (psilaxEz(i  ,j  ,1,l) - psilaxEz(i  ,j  ,1  ,l))/ Delz -  &
                           qint(i,j,1,l)       * Vz_elec


    if (source_rec_mpx == 0) then

    Exastsour(i,j,1,l)  =   sigma               * W_elec  *                       &
                           (Exint(i,j,1,l)      + vrotB_x   - Edotv * Vx_elec )

    Eyastsour(i,j,1,l)  =   sigma               * W_elec   *                      &
                           (Eyint(i,j,1,l)      + vrotB_y    - Edotv * Vy_elec )
    
    Ezastsour(i,j,1,l)  =   sigma               * W_elec   *                      &
                           (Ezint(i,j,1,l)      + vrotB_z    - Edotv * Vz_elec )

    else if (flux_solver == 5 .or. flux_solver == 6 .or. flux_solver == 7 .and. source_rec_mpx == 1) then


       Exastsour(i,j,1,l)  =  Exastsour(i,j,1,l) 
       Eyastsour(i,j,1,l)  =  Eyastsour(i,j,1,l)
       Ezastsour(i,j,1,l)  =  Ezastsour(i,j,1,l)

    else

       print*, "Revise the combination between flux_solver and source_rec_mpx parameters"
       print*, "It seems these two parameters are not valid"
       print*, "STOP subroutine 17_fluxandsource.f95"
       stop

    end if



  
    qflux(i,j,1,l)      = -(Jxlax(i,j,1,l) - Jxlax(i-1,j  ,1  ,l))/Delx -           &
                           (Jylax(i,j,1,l) - Jylax(i  ,j-1,1  ,l))/Dely -           &
                           (Jzlax(i,j,1,l) - Jzlax(i  ,j  ,1  ,l))/Delz 

    Dflux(i,j,1,l)      = -(FDxlax(i  ,j  ,1,l) - FDxlax(i-1,j  ,1  ,l))/Delx -       &
                           (FDylax(i  ,j  ,1,l) - FDylax(i  ,j-1,1  ,l))/Dely -       &
                           (FDzlax(i  ,j  ,1,l) - FDzlax(i  ,j  ,1  ,l))/Delz 



    tauflux(i,j,1,l)    = -(Ftauxlax(i  ,j  ,1,l) - Ftauxlax(i-1,j  ,1  ,l))/Delx -   & 
                           (Ftauylax(i  ,j  ,1,l) - Ftauylax(i  ,j-1,1  ,l))/Dely -   &
                           (Ftauzlax(i  ,j  ,1,l) - Ftauzlax(i  ,j  ,1  ,l))/Delz 

! -------------------------------------- EGLM--------------------------------------

    taueglm(i,j,1,l)    = -(Exint(i,j,1,l)          *                                  & 
                           (psitauxlax(i  ,j  ,1,l) - psitauxlax(i-1,j  ,1,l))/Delx +  &
                            Eyint(i,j,1,l)          *                                  & 
                           (psitauylax(i  ,j  ,1,l) - psitauylax(i  ,j-1,1,l))/Dely +  &
                            Ezint(i,j,1,l)          *                                  & 
                           (psitauzlax(i  ,j  ,1,l) - psitauzlax(i  ,j  ,1,l))/Delz )  &
                          -(Bxint(i,j,1,l)          *                                  & 
                           (phitauxlax(i  ,j  ,1,l) - phitauxlax(i-1,j  ,1,l))/Delx +  &
                            Byint(i,j,1,l)          *                                  & 
                           (phitauylax(i  ,j  ,1,l) - phitauylax(i  ,j-1,1,l))/Dely +  &
                            Bzint(i,j,1,l)          *                                  & 
                           (phitauzlax(i  ,j  ,1,l) - phitauzlax(i  ,j  ,1,l))/Delz )  

! -------------------------------------- EGLM--------------------------------------

    Sxflux(i,j,1,l)     = -(FSxxlax(i  ,j  ,1,l) - FSxxlax(i-1,j  ,1  ,l))/Delx -     &
                           (FSyxlax(i  ,j  ,1,l) - FSyxlax(i  ,j-1,1  ,l))/Dely -     &
                           (FSzxlax(i  ,j  ,1,l) - FSzxlax(i  ,j  ,1  ,l))/Delz 

    Syflux(i,j,1,l)     = -(FSxylax(i  ,j  ,1,l) - FSxylax(i-1,j  ,1  ,l))/Delx -     &
                           (FSyylax(i  ,j  ,1,l) - FSyylax(i  ,j-1,1  ,l))/Dely -     &
                           (FSzylax(i  ,j  ,1,l) - FSzylax(i  ,j  ,1  ,l))/Delz 

    Szflux(i,j,1,l)     = -(FSxzlax(i  ,j  ,1,l) - FSxzlax(i-1,j  ,1  ,l))/Delx -     &
                           (FSyzlax(i  ,j  ,1,l) - FSyzlax(i  ,j-1,1  ,l))/Dely -     &
                           (FSzzlax(i  ,j  ,1,l) - FSzzlax(i  ,j  ,1  ,l))/Delz 

! -------------------------------------- EGLM--------------------------------------
    Sxeglm(i,j,1,l)     = -(Bxint(i,j,1,l)       * (                                  &
                           (BxSxlax(i  ,j  ,1,l) - BxSxlax (i-1,j  ,1  ,l))/Delx +    &
                           (BySxlax(i  ,j  ,1,l) - BySxlax (i  ,j-1,1  ,l))/Dely +    &
                           (BzSxlax(i  ,j  ,1,l) - BzSxlax (i  ,j  ,1  ,l))/Delz ))

!                           (Bxint(i  ,j  ,1,l) - Bxint (i-1,j  ,1  ,l))/Delx +    &
!                           (Byint(i  ,j  ,1,l) - Byint (i  ,j-1,1  ,l))/Dely +    &
!                           (Bzint(i  ,j  ,1,l) - Bzint (i  ,j  ,1  ,l))/Delz ))

    Syeglm(i,j,1,l)     = -(Byint(i,j,1,l)       * (                                  &
                           (BxSylax(i  ,j  ,1,l) - BxSylax (i-1,j  ,1  ,l))/Delx +    &
                           (BySylax(i  ,j  ,1,l) - BySylax (i  ,j-1,1  ,l))/Dely +    &
                           (BzSylax(i  ,j  ,1,l) - BzSylax (i  ,j  ,1  ,l))/Delz ))

!                           (Bxint(i  ,j  ,1,l) - Bxint (i-1,j  ,1  ,l))/Delx +    &
!                           (Byint(i  ,j  ,1,l) - Byint (i  ,j-1,1  ,l))/Dely +    &
!                           (Bzint(i  ,j  ,1,l) - Bzint (i  ,j  ,1  ,l))/Delz ))



    Szeglm(i,j,1,l)     = -(Bzint(i,j,1,l)       * (                                  &
                           (BxSzlax(i  ,j  ,1,l) - BxSzlax (i-1,j  ,1  ,l))/Delx +    &
                           (BySzlax(i  ,j  ,1,l) - BySzlax (i  ,j-1,1  ,l))/Dely +    &
                           (BzSzlax(i  ,j  ,1,l) - BzSzlax (i  ,j  ,1  ,l))/Delz ))

!                           (Bxint(i  ,j  ,1,l) - Bxint (i-1,j  ,1  ,l))/Delx +    &
!                           (Byint(i  ,j  ,1,l) - Byint (i  ,j-1,1  ,l))/Dely +    &
!                           (Bzint(i  ,j  ,1,l) - Bzint (i  ,j  ,1  ,l))/Delz ))



! -------------------------------------- EGLM--------------------------------------

     end do ! for i
   end do ! for j

!$OMP END DO

   if ( source_rec == 1) then

      call source_boundary

!$OMP DO 
      
      do j=0,jmax
         do i=0,imax

    Ex_sour_rec_x_L = Exastsour(i,j,1,l)  +  FAC_SR  * mcl(Exastsour(i+1,j  ,1,l) - Exastsour(i  ,j  ,1,l) ,  &
                                                           Exastsour(i  ,j  ,1,l) - Exastsour(i-1,j  ,1,l)) 
    Ex_sour_rec_x_R = Exastsour(i,j,1,l)  -  FAC_SR  * mcl(Exastsour(i+2,j  ,1,l) - Exastsour(i+1,j  ,1,l) ,  &
                                                           Exastsour(i+1,j  ,1,l) - Exastsour(i  ,j  ,1,l)) 

    Ey_sour_rec_x_L = Exastsour(i,j,1,l)  +  FAC_SR  * mcl(Eyastsour(i+1,j  ,1,l) - Eyastsour(i  ,j  ,1,l) ,  &
                                                           Eyastsour(i  ,j  ,1,l) - Eyastsour(i-1,j  ,1,l)) 
    Ey_sour_rec_x_R = Exastsour(i,j,1,l)  -  FAC_SR  * mcl(Eyastsour(i+2,j  ,1,l) - Eyastsour(i+1,j  ,1,l) ,  &
                                                           Eyastsour(i+1,j  ,1,l) - Eyastsour(i  ,j  ,1,l)) 

    Ez_sour_rec_x_L = Exastsour(i,j,1,l)  +  FAC_SR  * mcl(Ezastsour(i+1,j  ,1,l) - Ezastsour(i  ,j  ,1,l) ,  &
                                                           Ezastsour(i  ,j  ,1,l) - Ezastsour(i-1,j  ,1,l)) 
    Ez_sour_rec_x_R = Exastsour(i,j,1,l)  -  FAC_SR  * mcl(Ezastsour(i+2,j  ,1,l) - Ezastsour(i+1,j  ,1,l) ,  &
                                                           Ezastsour(i+1,j  ,1,l) - Ezastsour(i  ,j  ,1,l)) 

    Ex_sour_rec_y_L = Exastsour(i,j,1,l)  +  FAC_SR  * mcl(Exastsour(i  ,j+1,1,l) - Exastsour(i  ,j  ,1,l) ,  &
                                                           Exastsour(i  ,j  ,1,l) - Exastsour(i  ,j-1,1,l)) 
    Ex_sour_rec_y_R = Exastsour(i,j,1,l)  -  FAC_SR  * mcl(Exastsour(i  ,j+2,1,l) - Exastsour(i  ,j+1,1,l) ,  &
                                                           Exastsour(i  ,j+1,1,l) - Exastsour(i  ,j  ,1,l)) 

    Ey_sour_rec_y_L = Exastsour(i,j,1,l)  +  FAC_SR  * mcl(Eyastsour(i  ,j+1,1,l) - Eyastsour(i  ,j  ,1,l) ,  &
                                                           Eyastsour(i  ,j  ,1,l) - Eyastsour(i  ,j-1,1,l)) 
    Ey_sour_rec_y_R = Exastsour(i,j,1,l)  -  FAC_SR  * mcl(Eyastsour(i  ,j+2,1,l) - Eyastsour(i  ,j+1,1,l) ,  &
                                                           Eyastsour(i  ,j+1,1,l) - Eyastsour(i  ,j  ,1,l)) 

    Ez_sour_rec_y_L = Exastsour(i,j,1,l)  +  FAC_SR  * mcl(Ezastsour(i  ,j+1,1,l) - Ezastsour(i  ,j  ,1,l) ,  &
                                                           Ezastsour(i  ,j  ,1,l) - Ezastsour(i  ,j-1,1,l)) 
    Ez_sour_rec_y_R = Exastsour(i,j,1,l)  -  FAC_SR  * mcl(Ezastsour(i  ,j+2,1,l) - Ezastsour(i  ,j+1,1,l) ,  &
                                                           Ezastsour(i  ,j+1,1,l) - Ezastsour(i  ,j  ,1,l)) 



    Ex_sour_rec_L   = 0.5d0 * (Ex_sour_rec_x_L + Ex_sour_rec_y_L)
    Ex_sour_rec_R   = 0.5d0 * (Ex_sour_rec_x_R + Ex_sour_rec_y_R)

    Ey_sour_rec_L   = 0.5d0 * (Ey_sour_rec_x_L + Ey_sour_rec_y_L)
    Ey_sour_rec_R   = 0.5d0 * (Ey_sour_rec_x_R + Ey_sour_rec_y_R)

    Ez_sour_rec_L   = 0.5d0 * (Ez_sour_rec_x_L + Ez_sour_rec_y_L)
    Ez_sour_rec_R   = 0.5d0 * (Ez_sour_rec_x_R + Ez_sour_rec_y_R)

    Exastsour(i,j,1,l) = 0.5d0 * (Ex_sour_rec_L + Ex_sour_rec_R)
    Eyastsour(i,j,1,l) = 0.5d0 * (Ey_sour_rec_L + Ey_sour_rec_R)
    Ezastsour(i,j,1,l) = 0.5d0 * (Ez_sour_rec_L + Ez_sour_rec_R)


          end do
       end do

!$OMP END DO
       
   end if

   else

      write(*,*) "STOP: subroutine fluxandsource"
      write(*,*) "This dimension is not implemented yet"
      stop

   end if

!$OMP END PARALLEL 

 end subroutine fluxandsource
