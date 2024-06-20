  !     *******************************************************************
  !     Subroutine which calculed  HLL flux flows 
  !     *******************************************************************

  subroutine hll_flow

    use scalar
    use parameters
    use threevectors
    use fourvectors
    use funciones


    implicit none

!$OMP   PARALLEL PRIVATE(psiint_L,psiint_R,phiint_L,phiint_R,Bxint_L,Bxint_R,Byint_L,Byint_R,Bzint_L,Bzint_R),                  &
!$OMP & PRIVATE(Bxint_L_m,Bxint_R_m,Byint_L_m,Byint_R_m,Bzint_L_m,Bzint_R_m,Exint_L,Exint_R,Eyint_L,Eyint_R,Ezint_L,Ezint_R),   &
!$OMP & PRIVATE(Exint_L_m,Exint_R_m,Eyint_L_m,Eyint_R_m,Ezint_L_m,Ezint_R_m,qint_L,qint_R,DDint_L,DDint_R),                     &
!$OMP & PRIVATE(tauint_L,tauint_R,Sxint_L,Sxint_R,Syint_L,Syint_R,Szint_L,Szint_R,E2B2int_L,E2B2int_R,Vx_L,Vx_R,Vy_L,Vy_R,Vz_L),&
!$OMP & PRIVATE(Vz_R,W_L,W_R,Edotv_L,Edotv_R,vrotB_x_L,vrotB_y_L,vrotB_z_L,vrotB_x_R,vrotB_y_R,vrotB_z_R,sigma_L,sigma_R),      &
!$OMP & PRIVATE(Jxint_L,Jxint_R,Jyint_L,Jyint_R,Jzint_L,Jzint_R,p_L,p_R,rho_L,rho_R,epsilon_L,epsilon_R,enthpy_L,enthpy_R),     &
!$OMP & PRIVATE(FDxint_L,FDxint_R,FDyint_L,FDyint_R,FDzint_L,FDzint_R),                                                         &
!$OMP & PRIVATE(Ftauxint_L,Ftauxint_R,Ftauyint_L,Ftauyint_R,Ftauzint_L,Ftauzint_R),                                             &
!$OMP & PRIVATE(FSxxint_L,FSxxint_R,FSxyint_L,FSxyint_R,FSxzint_L,FSxzint_R),                                                   &
!$OMP & PRIVATE(FSyxint_L,FSyxint_R,FSyyint_L,FSyyint_R,FSyzint_L,FSyzint_R),                                                   &
!$OMP & PRIVATE(FSzxint_L,FSzxint_R,FSzyint_L,FSzyint_R,FSzzint_L,FSzzint_R),                                                   &
!$OMP & PRIVATE(psiint_x_L,psiint_x_R,phiint_x_L,phiint_x_R,psiint_y_L,psiint_y_R,phiint_y_L,phiint_y_R),                       &
!$OMP & PRIVATE(Bxint_x_L,Bxint_x_R,Byint_x_L,Byint_x_R,Bzint_x_L,Bzint_x_R),                                                   &
!$OMP & PRIVATE(Bxint_y_L,Bxint_y_R,Byint_y_L,Byint_y_R,Bzint_y_L,Bzint_y_R),                                                   &
!$OMP & PRIVATE(Bxint_x_L_m,Bxint_x_R_m,Byint_x_L_m,Byint_x_R_m,Bzint_x_L_m,Bzint_x_R_m),                                       &
!$OMP & PRIVATE(Bxint_y_L_m,Bxint_y_R_m,Byint_y_L_m,Byint_y_R_m,Bzint_y_L_m,Bzint_y_R_m),                                       &
!$OMP & PRIVATE(Exint_x_L,Exint_x_R,Eyint_x_L,Eyint_x_R,Ezint_x_L,Ezint_x_R),                                                   &
!$OMP & PRIVATE(Exint_y_L,Exint_y_R,Eyint_y_L,Eyint_y_R,Ezint_y_L,Ezint_y_R),                                                   &
!$OMP & PRIVATE(Exint_x_L_m,Exint_x_R_m,Eyint_x_L_m,Eyint_x_R_m,Ezint_x_L_m,Ezint_x_R_m),                                       &
!$OMP & PRIVATE(Exint_y_L_m,Exint_y_R_m,Eyint_y_L_m,Eyint_y_R_m,Ezint_y_L_m,Ezint_y_R_m),                                       &
!$OMP & PRIVATE(qint_x_L,qint_x_R,qint_y_L,qint_y_R,DDint_x_L,DDint_x_R,DDint_y_L,DDint_y_R),                                   &
!$OMP & PRIVATE(tauint_x_L,tauint_x_R,tauint_y_L,tauint_y_R),                                                                   &
!$OMP & PRIVATE(Sxint_x_L,Sxint_x_R,Syint_x_L,Syint_x_R,Szint_x_L,Szint_x_R),                                                   &
!$OMP & PRIVATE(Sxint_y_L,Sxint_y_R,Syint_y_L,Syint_y_R,Szint_y_L,Szint_y_R),                                                   &
!$OMP & PRIVATE(E2B2int_x_L,E2B2int_x_R,E2B2int_y_L,E2B2int_y_R,Vx_x_L,Vx_x_R,Vy_x_L,Vy_x_R,Vz_x_L,Vz_x_R),                     &
!$OMP & PRIVATE(Edotv_x_L,Edotv_x_R,Edotv_y_L,Edotv_y_R,Vx_y_L,Vx_y_R,Vy_y_L,Vy_y_R,Vz_y_L,Vz_y_R),                             &
!$OMP & PRIVATE(Vx_x,Vx_y,Vy_x,Vy_y,Vz_x,Vz_y,Vx_rec,Vy_rec,Vz_rec,W_x_L,W_x_R,W_y_L,W_y_R),                                    &
!$OMP & PRIVATE(vrotB_x_x_L,vrotB_y_x_L,vrotB_z_x_L,vrotB_x_x_R,vrotB_y_x_R,vrotB_z_x_R),                                       &
!$OMP & PRIVATE(vrotB_x_y_L,vrotB_y_y_L,vrotB_z_y_L,vrotB_x_y_R,vrotB_y_y_R,vrotB_z_y_R),                                       &
!$OMP & PRIVATE(p_x_L,p_x_R,rho_x_L,rho_x_R,p_y_L,p_y_R,rho_y_L,rho_y_R,epsilon_x_L,epsilon_x_R,epsilon_x), & 
!$OMP & PRIVATE(epsilon_y,epsilon_y_L,epsilon_y_R,enthpy_y_L,enthpy_y_R,enthpy_x_L,enthpy_x_R,epsilon_rec,c_sp,Sp)

    
    if (DIM == 1) then

       call boundary_conserved
       call boundary_electric
       call boundary_flux_conserved

!$OMP DO 
       
    do i=-1,imax

       !---------------------------------
       ! Left - Right Conserved variables
       !---------------------------------

    !     Left - Right  intermediate values of augmented fields

    psiint_L = psiint(i  ,1,1,l) +  FAC  * mcl(psiint(i+1,1,1,l) - psiint(i  ,1,1,l) ,  &
                                               psiint(i  ,1,1,l) - psiint(i-1,1,1,l))
    psiint_R = psiint(i+1,1,1,l) -  FAC  * mcl(psiint(i+2,1,1,l) - psiint(i+1,1,1,l) ,  &
                                               psiint(i+1,1,1,l) - psiint(i  ,1,1,l))

    phiint_L = phiint(i  ,1,1,l) +  FAC  * mcl(phiint(i+1,1,1,l) - phiint(i  ,1,1,l) ,  &
                                               phiint(i  ,1,1,l) - phiint(i-1,1,1,l))
    phiint_R = phiint(i+1,1,1,l) -  FAC  * mcl(phiint(i+2,1,1,l) - phiint(i+1,1,1,l) ,  &
                                               phiint(i+1,1,1,l) - phiint(i  ,1,1,l))

    !     Left - Right intermediate values of magnetic field

    Bxint_L =  Bxint(i  ,1,1,l) +  FAC  * mcl(Bxint(i+1,1,1,l) - Bxint(i  ,1,1,l) ,  &
                                              Bxint(i  ,1,1,l) - Bxint(i-1,1,1,l)) 
    Bxint_R =  Bxint(i+1,1,1,l) -  FAC  * mcl(Bxint(i+2,1,1,l) - Bxint(i+1,1,1,l) ,  &
                                              Bxint(i+1,1,1,l) - Bxint(i  ,1,1,l)) 

    Byint_L =  Byint(i  ,1,1,l) +  FAC  * mcl(Byint(i+1,1,1,l) - Byint(i  ,1,1,l) ,  &
                                              Byint(i  ,1,1,l) - Byint(i-1,1,1,l))  
    Byint_R =  Byint(i+1,1,1,l) -  FAC  * mcl(Byint(i+2,1,1,l) - Byint(i+1,1,1,l) ,  &
                                              Byint(i+1,1,1,l) - Byint(i  ,1,1,l))  

    Bzint_L =  Bzint(i  ,1,1,l) +  FAC  * mcl(Bzint(i+1,1,1,l) - Bzint(i  ,1,1,l) ,  &
                                              Bzint(i  ,1,1,l) - Bzint(i-1,1,1,l))  
    Bzint_R =  Bzint(i+1,1,1,l) -  FAC  * mcl(Bzint(i+2,1,1,l) - Bzint(i+1,1,1,l) ,  &
                                              Bzint(i+1,1,1,l) - Bzint(i  ,1,1,l))  
    !-----------------------------------------------------------------------------------

    Bxint_L_m =  -Bxint(i  ,1,1,l) +  FAC  * mcl(Bxint(i+1,1,1,l) - Bxint(i  ,1,1,l) ,  &
                                                 Bxint(i  ,1,1,l) - Bxint(i-1,1,1,l)) 
    Bxint_R_m =  -Bxint(i+1,1,1,l) -  FAC  * mcl(Bxint(i+2,1,1,l) - Bxint(i+1,1,1,l) ,  &
                                                 Bxint(i+1,1,1,l) - Bxint(i  ,1,1,l)) 

    Byint_L_m =  -Byint(i  ,1,1,l) +  FAC  * mcl(Byint(i+1,1,1,l) - Byint(i  ,1,1,l) ,  &
                                                 Byint(i  ,1,1,l) - Byint(i-1,1,1,l))  
    Byint_R_m =  -Byint(i+1,1,1,l) -  FAC  * mcl(Byint(i+2,1,1,l) - Byint(i+1,1,1,l) ,  &
                                                 Byint(i+1,1,1,l) - Byint(i  ,1,1,l))  

    Bzint_L_m =  -Bzint(i  ,1,1,l) +  FAC  * mcl(Bzint(i+1,1,1,l) - Bzint(i  ,1,1,l) ,  &
                                                 Bzint(i  ,1,1,l) - Bzint(i-1,1,1,l))  
    Bzint_R_m =  -Bzint(i+1,1,1,l) -  FAC  * mcl(Bzint(i+2,1,1,l) - Bzint(i+1,1,1,l) ,  &
                                                 Bzint(i+1,1,1,l) - Bzint(i  ,1,1,l))  
    !-----------------------------------------------------------------------------------

    !     Left - Right intermediate value of electric field

    Exint_L =  Exint(i  ,1,1,l)  +  FAC  * mcl(Exint(i+1,1,1,l) - Exint(i  ,1,1,l) ,  &
                                               Exint(i  ,1,1,l) - Exint(i-1,1,1,l)) 
    Exint_R =  Exint(i+1,1,1,l)  -  FAC  * mcl(Exint(i+2,1,1,l) - Exint(i+1,1,1,l) ,  &
                                               Exint(i+1,1,1,l) - Exint(i  ,1,1,l)) 

    Eyint_L =  Eyint(i  ,1,1,l)  +  FAC  * mcl(Eyint(i+1,1,1,l) - Eyint(i  ,1,1,l) ,  &
                                               Eyint(i  ,1,1,l) - Eyint(i-1,1,1,l)) 
    Eyint_R =  Eyint(i+1,1,1,l)  -  FAC  * mcl(Eyint(i+2,1,1,l) - Eyint(i+1,1,1,l) ,  &
                                               Eyint(i+1,1,1,l) - Eyint(i  ,1,1,l)) 

    Ezint_L =  Ezint(i  ,1,1,l)  +  FAC  * mcl(Ezint(i+1,1,1,l) - Ezint(i  ,1,1,l) ,  &
                                               Ezint(i  ,1,1,l) - Ezint(i-1,1,1,l)) 
    Ezint_R =  Ezint(i+1,1,1,l)  -  FAC  * mcl(Ezint(i+2,1,1,l) - Ezint(i+1,1,1,l) ,  &
                                               Ezint(i+1,1,1,l) - Ezint(i  ,1,1,l)) 
    !-----------------------------------------------------------------------------------
    Exint_L_m =  - Exint(i  ,1,1,l)  +  FAC  * mcl(Exint(i+1,1,1,l) - Exint(i  ,1,1,l) ,  &
                                                   Exint(i  ,1,1,l) - Exint(i-1,1,1,l)) 
    Exint_R_m =  - Exint(i+1,1,1,l)  -  FAC  * mcl(Exint(i+2,1,1,l) - Exint(i+1,1,1,l) ,  &
                                                   Exint(i+1,1,1,l) - Exint(i  ,1,1,l)) 

    Eyint_L_m =  - Eyint(i  ,1,1,l)  +  FAC  * mcl(Eyint(i+1,1,1,l) - Eyint(i  ,1,1,l) ,  &
                                                   Eyint(i  ,1,1,l) - Eyint(i-1,1,1,l)) 
    Eyint_R_m =  - Eyint(i+1,1,1,l)  -  FAC  * mcl(Eyint(i+2,1,1,l) - Eyint(i+1,1,1,l) ,  &
                                                   Eyint(i+1,1,1,l) - Eyint(i  ,1,1,l)) 

    Ezint_L_m =  - Ezint(i  ,1,1,l)  +  FAC  * mcl(Ezint(i+1,1,1,l) - Ezint(i  ,1,1,l) ,  &
                                                   Ezint(i  ,1,1,l) - Ezint(i-1,1,1,l)) 
    Ezint_R_m =  - Ezint(i+1,1,1,l)  -  FAC  * mcl(Ezint(i+2,1,1,l) - Ezint(i+1,1,1,l) ,  &
                                                   Ezint(i+1,1,1,l) - Ezint(i  ,1,1,l)) 
    !-----------------------------------------------------------------------------------

    !     Left - Right Intermediate value of charge density 

    qint_L = qint(i  ,1,1,l)  +  FAC  * mcl(qint(i+1,1,1,l) - qint(i  ,1,1,l) ,  &
                                            qint(i  ,1,1,l) - qint(i-1,1,1,l)) 
    qint_R = qint(i+1,1,1,l)  -  FAC  * mcl(qint(i+2,1,1,l) - qint(i+1,1,1,l) ,  &
                                            qint(i+1,1,1,l) - qint(i  ,1,1,l)) 

    !     Mass Conservation Explicit intermediate value of D

    DDint_L = DDint(i  ,1,1,l) +  FAC  * mcl(DDint(i+1,1,1,l) - DDint(i  ,1,1,l) ,  &
                                             DDint(i  ,1,1,l) - DDint(i-1,1,1,l)) 
    DDint_R = DDint(i+1,1,1,l) -  FAC  * mcl(DDint(i+2,1,1,l) - DDint(i+1,1,1,l) ,  &
                                             DDint(i+1,1,1,l) - DDint(i  ,1,1,l)) 

    !     Energy Conservation Explicit intermediate value of tau

    tauint_L = tauint(i  ,1,1,l) +  FAC  * mcl(tauint(i+1,1,1,l) - tauint(i  ,1,1,l) ,  &
                                               tauint(i  ,1,1,l) - tauint(i-1,1,1,l)) 
    tauint_R = tauint(i+1,1,1,l) -  FAC  * mcl(tauint(i+2,1,1,l) - tauint(i+1,1,1,l) ,  &
                                               tauint(i+1,1,1,l) - tauint(i  ,1,1,l)) 

    !    Explicit intermediate value of S

    Sxint_L = Sxint(i  ,1,1,l) +  FAC  * mcl(Sxint(i+1,1,1,l) - Sxint(i  ,1,1,l) ,  &
                                             Sxint(i  ,1,1,l) - Sxint(i-1,1,1,l)) 
    Sxint_R = Sxint(i+1,1,1,l) -  FAC  * mcl(Sxint(i+2,1,1,l) - Sxint(i+1,1,1,l) ,  &
                                             Sxint(i+1,1,1,l) - Sxint(i  ,1,1,l)) 

    Syint_L = Syint(i  ,1,1,l) +  FAC  * mcl(Syint(i+1,1,1,l) - Syint(i  ,1,1,l) ,  &
                                             Syint(i  ,1,1,l) - Syint(i-1,1,1,l))
    Syint_R = Syint(i+1,1,1,l) -  FAC  * mcl(Syint(i+2,1,1,l) - Syint(i+1,1,1,l) ,  &
                                             Syint(i+1,1,1,l) - Syint(i  ,1,1,l))


    Szint_L = Szint(i  ,1,1,l) +  FAC  * mcl(Szint(i+1,1,1,l) - Szint(i  ,1,1,l) ,  &
                                             Szint(i  ,1,1,l) - Szint(i-1,1,1,l))
    Szint_R = Szint(i+1,1,1,l) -  FAC  * mcl(Szint(i+2,1,1,l) - Szint(i+1,1,1,l) ,  &
                                             Szint(i+1,1,1,l) - Szint(i  ,1,1,l))


       !---------------------------------
       ! Left - Right Fluxes
       !---------------------------------

    !     Electric and Magnetic field modules 


    E2B2int_L = Exint_L**2 + Eyint_L**2 + Ezint_L**2 + &
                Bxint_L**2 + Byint_L**2 + Bzint_L**2

    E2B2int_R = Exint_R**2 + Eyint_R**2 + Ezint_R**2 + &
                Bxint_R**2 + Byint_R**2 + Bzint_R**2




          if (REC_PRIM == 0) then

             Vx_L = Vx(i  ,1,1) +  FAC  * mcl(Vx(i+1,1,1) - Vx(i  ,1,1) ,  &
                                              Vx(i  ,1,1) - Vx(i-1,1,1))
             Vx_R = Vx(i+1,1,1) -  FAC  * mcl(Vx(i+2,1,1) - Vx(i+1,1,1) ,  &
                                              Vx(i+1,1,1) - Vx(i  ,1,1))

             Vy_L = Vy(i  ,1,1) +  FAC  * mcl(Vy(i+1,1,1) - Vy(i  ,1,1) ,  &
                                              Vy(i  ,1,1) - Vy(i-1,1,1))
             Vy_R = Vy(i+1,1,1) -  FAC  * mcl(Vy(i+2,1,1) - Vy(i+1,1,1) ,  &
                                              Vy(i+1,1,1) - Vy(i  ,1,1))

             Vz_L = Vz(i  ,1,1) +  FAC  * mcl(Vz(i+1,1,1) - Vz(i  ,1,1) ,  &
                                              Vz(i  ,1,1) - Vz(i-1,1,1))
             Vz_R = Vz(i+1,1,1) -  FAC  * mcl(Vz(i+2,1,1) - Vz(i+1,1,1) ,  &
                                              Vz(i+1,1,1) - Vz(i  ,1,1))

             W_L = 1.d0/sqrt(1.d0 - Vx_L**2 -  Vy_L**2 - Vz_L**2)
             W_R = 1.d0/sqrt(1.d0 - Vx_R**2 -  Vy_R**2 - Vz_R**2)


          else if (REC_PRIM == 1) then

             Vx_L = Vxint(i  ,1,1,l) +  FAC  * mcl(Vxint(i+1,1,1,l) - Vxint(i  ,1,1,l) ,  &
                                                   Vxint(i  ,1,1,l) - Vxint(i-1,1,1,l))
             Vx_R = Vxint(i+1,1,1,l) -  FAC  * mcl(Vxint(i+2,1,1,l) - Vxint(i+1,1,1,l) ,  &
                                                   Vxint(i+1,1,1,l) - Vxint(i  ,1,1,l))

             Vy_L = Vyint(i  ,1,1,l) +  FAC  * mcl(Vyint(i+1,1,1,l) - Vyint(i  ,1,1,l) ,  &
                                                   Vyint(i  ,1,1,l) - Vyint(i-1,1,1,l))
             Vy_R = Vyint(i+1,1,1,l) -  FAC  * mcl(Vyint(i+2,1,1,l) - Vyint(i+1,1,1,l) ,  &
                                                   Vyint(i+1,1,1,l) - Vyint(i  ,1,1,l))

             Vz_L = Vzint(i  ,1,1,l) +  FAC  * mcl(Vzint(i+1,1,1,l) - Vzint(i  ,1,1,l) ,  &
                                                   Vzint(i  ,1,1,l) - Vzint(i-1,1,1,l))
             Vy_R = Vzint(i+1,1,1,l) -  FAC  * mcl(Vzint(i+2,1,1,l) - Vzint(i+1,1,1,l) ,  &
                                                   Vzint(i+1,1,1,l) - Vzint(i  ,1,1,l))

             W_L = 1.d0/sqrt(1.d0 - Vx_L**2 -  Vy_L**2 - Vz_L**2)
             W_R = 1.d0/sqrt(1.d0 - Vx_R**2 -  Vy_R**2 - Vz_R**2)


          else

             write(*,*) "STOP: subroutine hllflux"
             write(*,*) "REC_PRIM parameter is not valid"
             stop

          end if



             Edotv_L = Exint_L * Vx_L + Eyint_L * Vy_L + Ezint_L * Vz_L
             Edotv_R = Exint_R * Vx_R + Eyint_R * Vy_R + Ezint_R * Vz_R

             vrotB_x_L = (Bzint_L * Vy_L - Byint_L * Vz_L)
             vrotB_y_L = (Bxint_L * Vz_L - Bzint_L * Vx_L)
             vrotB_z_L = (Byint_L * Vx_L - Bxint_L * Vy_L)

             vrotB_x_R = (Bzint_R * Vy_R - Byint_R * Vz_R)
             vrotB_y_R = (Bxint_R * Vz_R - Bzint_R * Vx_R)
             vrotB_z_R = (Byint_R * Vx_R - Bxint_R * Vy_R)


              if (SLC == 1) then
                 
                 sigma_L = sigma_0 * DDint_L**gamma_slc 
                 sigma_R = sigma_0 * DDint_R**gamma_slc 

              else if (SLC == 0) then

                 sigma_L = sigma
                 sigma_R = sigma

              end if


    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxint_L  = sigma_L * W_L * ( Exint_L +  vrotB_x_L  - Edotv_L * Vx_L )   &
                           +   qint_L  *  Vx_L
    Jxint_R  = sigma_R * W_R * ( Exint_R +  vrotB_x_R  - Edotv_R * Vx_R )   &
                           +   qint_R  *  Vx_R


    Jyint_L = sigma_L * W_L  * ( Eyint_L +  vrotB_y_L  - Edotv_L * Vy_L )   &
                           +   qint_L  *  Vy_L
    Jyint_R = sigma_R * W_R  * ( Eyint_R +  vrotB_y_R  - Edotv_R * Vy_R )   &
                           +   qint_R  *  Vy_R

    Jzint_L = sigma_L * W_L  * ( Ezint_L +  vrotB_z_L  - Edotv_L * Vz_L )   &
                           +   qint_L  *  Vz_L
    Jzint_R = sigma_R * W_R  * ( Ezint_R +  vrotB_z_R  - Edotv_R * Vz_R )   &
                           +   qint_R  *  Vz_R

    !     Density flux (FD = \rho W V)

    FDxint_L = DDint_L * Vx_L
    FDxint_R = DDint_R * Vx_R

    FDyint_L = DDint_L * Vy_L
    FDyint_R = DDint_R * Vy_R

    FDzint_L = DDint_L * Vz_L
    FDzint_R = DDint_R * Vz_R

              !     Entalphy 

    p_L   = p(i,1,1)     +  FAC  * mcl(p(i+1,1,1) - p(i  ,1,1) ,  &
                                       p(i  ,1,1) - p(i-1,1,1))
    p_R   = p(i+1,1,1)   -  FAC  * mcl(p(i+2,1,1) - p(i+1,1,1) ,  &
                                       p(i+1,1,1) - p(i  ,1,1))

    rho_L = rho(i,1,1)   + FAC  * mcl(rho(i+1,1,1) - rho(i  ,1,1) ,  &
                                      rho(i  ,1,1) - rho(i-1,1,1))
    rho_R = rho(i+1,1,1) - FAC  * mcl(rho(i+2,1,1) - rho(i+1,1,1) ,  &
                                      rho(i+1,1,1) - rho(i  ,1,1))

!!$    if (p_L .lt. 0 .or. p_R .lt. 0) then 
!!$
!!$       print*, "negative pressure in HLL ---> p_L =", p_L, "p_R =", p_R, "i =", i
!!$
!!$    end if

!!$
!!$    p_L   = pint(i  ,1,1,l)     +  FAC  * mcl(pint(i+1,1,1,l) - pint(i  ,1,1,l) ,  &
!!$                                              pint(i  ,1,1,l) - pint(i-1,1,1,l))
!!$    p_R   = pint(i+1,1,1,l)     -  FAC  * mcl(pint(i+2,1,1,l) - pint(i+1,1,1,l) ,  &
!!$                                              pint(i+1,1,1,l) - pint(i  ,1,1,l))
!!$
!!$    rho_L   = rhoint(i  ,1,1,l) +  FAC  * mcl(rhoint(i+1,1,1,l) - rhoint(i  ,1,1,l) ,  &
!!$                                              rhoint(i  ,1,1,l) - rhoint(i-1,1,1,l))
!!$    rho_R   = rhoint(i+1,1,1,l) -  FAC  * mcl(rhoint(i+2,1,1,l) - rhoint(i+1,1,1,l) ,  &
!!$                                              rhoint(i+1,1,1,l) - rhoint(i  ,1,1,l))

!!$    if (i == 100) then
!!$
!!$       print*, p_L
!!$
!!$    end if

    epsilon_L = p_L / ((gamma-1.d0) * rho_L)
    epsilon_R = p_R / ((gamma-1.d0) * rho_R)

    enthpy_L   =  rho_L * ( 1.d0 + epsilon_L ) + p_L 
    enthpy_R   =  rho_R * ( 1.d0 + epsilon_R ) + p_R 
!!$    
!!$       enthpy_L   =  enthpyint(i  ,1,1,l) +  FAC  * mcl(enthpyint(i+1,1,1,l) - enthpyint(i  ,1,1,l) ,  &
!!$                                                        enthpyint(i  ,1,1,l) - enthpyint(i-1,1,1,l))
!!$       enthpy_R   =  enthpyint(i+1,1,1,l) -  FAC  * mcl(enthpyint(i+2,1,1,l) - enthpyint(i+1,1,1,l) ,  &
!!$                                                        enthpyint(i+1,1,1,l) - enthpyint(i  ,1,1,l))



    !     Flujos conservados de energía

    Ftauxint_L = (Bzint_L  * Eyint_L - Byint_L * Ezint_L) +   &
                  enthpy_L * W_L**2  * Vx_L   
    Ftauxint_R = (Bzint_R  * Eyint_R - Byint_R * Ezint_R) +   &
                  enthpy_R * W_R**2  * Vx_R   

    Ftauyint_L = (Bxint_L  * Ezint_L - Bzint_L * Exint_L) +   &
                  enthpy_L * W_L**2  * Vy_L
    Ftauyint_R = (Bxint_R  * Ezint_R - Bzint_R * Exint_R) +   &
                  enthpy_R * W_R**2  * Vy_R

    Ftauzint_L = (Byint_L  * Exint_L - Bxint_L * Eyint_L) +   &
                  enthpy_L * W_L**2  * Vz_L
    Ftauzint_R = (Byint_R  * Exint_R - Bxint_R * Eyint_R) +   &
                  enthpy_R * W_R**2  * Vz_R

    !     Components of the flux momentum tensor 

    FSxxint_L =  - Exint_L**2 - Bxint_L**2 + enthpy_L * W_L**2 * Vx_L**2 + &
                   0.5d0 * E2B2int_L + p_L  
    FSxxint_R =  - Exint_R**2 - Bxint_R**2 + enthpy_R * W_R**2 * Vx_R**2 + &
                   0.5d0 * E2B2int_R + p_R  

    FSxyint_L = -  Exint_L  * Eyint_L - Bxint_L * Byint_L +    &
                   enthpy_L * W_L**2  * Vx_L    * Vy_L
    FSxyint_R = -  Exint_R  * Eyint_R - Bxint_R * Byint_R +    &
                   enthpy_R * W_R**2  * Vx_R    * Vy_R

    FSxzint_L = -  Exint_L  * Ezint_L - Bxint_L * Bzint_L +    &
                   enthpy_L * W_L**2  * Vx_L    * Vz_L
    FSxzint_R = -  Exint_R  * Ezint_R - Bxint_R * Bzint_R +    &
                   enthpy_R * W_R**2  * Vx_R    * Vz_R

    FSyxint_L = -  Exint_L  * Eyint_L - Bxint_L * Byint_L +    &
                   enthpy_L * W_L**2  * Vx_L    * Vy_L
    FSyxint_R = -  Exint_R  * Eyint_R - Bxint_R * Byint_R +    &
                   enthpy_R * W_R**2  * Vx_R    * Vy_R

    FSyyint_L = -  Eyint_L**2 - Byint_L**2 + enthpy_L * W_L**2 * Vy_L**2 + &
                   0.5d0 * E2B2int_L  + p_L
    FSyyint_R = -  Eyint_R**2 - Byint_R**2 + enthpy_R * W_R**2 * Vy_R**2 + &
                   0.5d0 * E2B2int_R  + p_R

    FSyzint_L = -  Eyint_L  * Ezint_L - Byint_L * Bzint_L +    &
                   enthpy_L * W_L**2  * Vy_L    * Vz_L
    FSyzint_R = -  Eyint_R  * Ezint_R - Byint_R * Bzint_R +    &
                   enthpy_R * W_R**2  * Vy_R    * Vz_R

    FSzxint_L = -  Exint_L  * Ezint_L - Bxint_L * Bzint_L +    &
                   enthpy_L * W_L**2  * Vx_L    * Vz_L       
    FSzxint_R = -  Exint_R  * Ezint_R - Bxint_R * Bzint_R +    &
                   enthpy_R * W_R**2  * Vx_R    * Vz_R             

    FSzyint_L = -  Eyint_L  * Ezint_L - Byint_L * Bzint_L +    &
                   enthpy_L * W_L**2  * Vy_L    * Vz_L
    FSzyint_R = -  Eyint_R  * Ezint_R - Byint_R * Bzint_R +    &
                   enthpy_R * W_R**2  * Vy_R    * Vz_R

    FSzzint_L = -  Ezint_L**2 - Bzint_L**2 + enthpy_L * W_L**2 * Vz_L**2 + &
                   0.5d0 * E2B2int_L + p_L
    FSzzint_R = -  Ezint_R**2 - Bzint_R**2 + enthpy_R * W_R**2 * Vz_R**2 + &
                   0.5d0 * E2B2int_R + p_R

! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,1,1,l)  =   0.5d0 * (Exint_L + Exint_R - ( psiint_R - psiint_L ) )  
       Eylaxpsi(i,1,1,l)  =   0.5d0 * (Eyint_L + Eyint_R - ( psiint_R - psiint_L ) ) 
       Ezlaxpsi(i,1,1,l)  =   0.5d0 * (Ezint_L + Ezint_R - ( psiint_R - psiint_L ) ) 

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,1,1,l)  =   0.5d0 * (Bxint_L + Bxint_R - ( phiint_R - phiint_L ) )
       Bylaxphi(i,1,1,l)  =   0.5d0 * (Byint_L + Byint_R - ( phiint_R - phiint_L ) )
       Bzlaxphi(i,1,1,l)  =   0.5d0 * (Bzint_L + Bzint_R - ( phiint_R - phiint_L ) )

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,1,1,l)  =   0.5d0 * (phiint_L   + phiint_R   - ( Bxint_R - Bxint_L ) )
       EzlaxBx (i,1,1,l)  =   0.5d0 * (Ezint_L    + Ezint_R    - ( Bxint_R - Bxint_L ) )
       EylaxBx (i,1,1,l)  = - 0.5d0 * (Eyint_L_m  + Eyint_R_m  - ( Bxint_R - Bxint_L ) )

! Flujos By

       EzlaxBy (i,1,1,l)  = - 0.5d0 * (Ezint_L_m  + Ezint_R_m  - ( Byint_R - Byint_L ) )
       philaxBy(i,1,1,l)  =   0.5d0 * (phiint_L   + phiint_R   - ( Byint_R - Byint_L ) )
       ExlaxBy (i,1,1,l)  =   0.5d0 * (Exint_L    + Exint_R    - ( Byint_R - Byint_L ) )

! Flujos Bz

       EylaxBz (i,1,1,l)  =   0.5d0 * (Eyint_L    + Eyint_R    - ( Bzint_R - Bzint_L ) )
       ExlaxBz (i,1,1,l)  = - 0.5d0 * (Exint_L_m  + Exint_R_m  - ( Bzint_R - Bzint_L ) )
       philaxBz(i,1,1,l)  =   0.5d0 * (phiint_L   + phiint_R   - ( Bzint_R - Bzint_L ) )

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,1,1,l)  =   0.5d0 * (psiint_L   + psiint_R   - ( Exint_R - Exint_L ) )
       BzlaxEx(i,1,1,l)   = - 0.5d0 * (Bzint_L_m  + Bzint_R_m  - ( Exint_R - Exint_L ) )
       BylaxEx(i,1,1,l)   =   0.5d0 * (Byint_L    + Byint_R    - ( Exint_R - Exint_L ) )

! Flujos Ey


       BzlaxEy (i,1,1,l)  =   0.5d0 * (Bzint_L    + Bzint_R    - ( Eyint_R - Eyint_L ) )
       psilaxEy(i,1,1,l)  =   0.5d0 * (psiint_L   + psiint_R   - ( Eyint_R - Eyint_L ) )
       BxlaxEy (i,1,1,l)  = - 0.5d0 * (Bxint_L_m  + Bxint_R_m  - ( Eyint_R - Eyint_L ) )

! Flujos Ez


       BylaxEz (i,1,1,l)  = - 0.5d0 * (Byint_L_m  + Byint_R_m  - ( Ezint_R - Ezint_L ) )
       BxlaxEz (i,1,1,l)  =   0.5d0 * (Bxint_L    + Bxint_R    - ( Ezint_R - Ezint_L ) )
       psilaxEz(i,1,1,l)  =   0.5d0 * (psiint_L   + psiint_R   - ( Ezint_R - Ezint_L ) )

! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,1,1,l)     =   0.5d0 * (Jxint_L  + Jxint_R - ( qint_R - qint_L ) )
       Jylax(i,1,1,l)     =   0.5d0 * (Jyint_L  + Jyint_R - ( qint_R - qint_L ) )
       Jzlax(i,1,1,l)     =   0.5d0 * (Jzint_L  + Jzint_R - ( qint_R - qint_L ) )

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,1,1,l)    =   0.5d0 * (FDxint_L  + FDxint_R - ( DDint_R - DDint_L ) )
       FDylax(i,1,1,l)    =   0.5d0 * (FDyint_L  + FDyint_R - ( DDint_R - DDint_L ) )
       FDzlax(i,1,1,l)    =   0.5d0 * (FDzint_L  + FDzint_R - ( DDint_R - DDint_L ) )

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,1,1,l)  =   0.5d0 * (Ftauxint_L + Ftauxint_R - ( tauint_R - tauint_L  ) )
       Ftauylax(i,1,1,l)  =   0.5d0 * (Ftauyint_L + Ftauyint_R - ( tauint_R - tauint_L  ) )
       Ftauzlax(i,1,1,l)  =   0.5d0 * (Ftauzint_L + Ftauzint_R - ( tauint_R - tauint_L  ) )

! -------------------------------------- EGLM-------------------------------------- ?????? REVISAR ??????

       psitauxlax(i,1,1,l)  =   0.5d0 * (psiint_L + psiint_R - ( tauint_R - tauint_L  ) )
       psitauylax(i,1,1,l)  =   0.5d0 * (psiint_L + psiint_R - ( tauint_R - tauint_L  ) )
       psitauzlax(i,1,1,l)  =   0.5d0 * (psiint_L + psiint_R - ( tauint_R - tauint_L  ) )


       phitauxlax(i,1,1,l)  =   0.5d0 * (phiint_L + phiint_R - ( tauint_R - tauint_L  ) )
       phitauylax(i,1,1,l)  =   0.5d0 * (phiint_L + phiint_R - ( tauint_R - tauint_L  ) )
       phitauzlax(i,1,1,l)  =   0.5d0 * (phiint_L + phiint_R - ( tauint_R - tauint_L  ) )

! -------------------------------------- EGLM-------------------------------------- ?????? REVISAR ??????

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,1,1,l)   =   0.5d0 * (FSxxint_L + FSxxint_R - (Sxint_R - Sxint_L ) ) 
       FSxylax(i,1,1,l)   =   0.5d0 * (FSxyint_L + FSxyint_R - (Syint_R - Syint_L ) ) 
       FSxzlax(i,1,1,l)   =   0.5d0 * (FSxzint_L + FSxzint_R - (Szint_R - Szint_L ) ) 

 
       FSyxlax(i,1,1,l)   =   0.5d0 * (FSyxint_L + FSyxint_R - (Sxint_R - Sxint_L ) ) 
       FSyylax(i,1,1,l)   =   0.5d0 * (FSyyint_L + FSyyint_R - (Syint_R - Syint_L ) )
       FSyzlax(i,1,1,l)   =   0.5d0 * (FSyzint_L + FSyzint_R - (Szint_R - Szint_L ) )


       FSzxlax(i,1,1,l)  =   0.5d0 * (FSzxint_L + FSzxint_R - (Sxint_R - Sxint_L ) )
       FSzylax(i,1,1,l)  =   0.5d0 * (FSzyint_L + FSzyint_R - (Syint_R - Syint_L ) )
       FSzzlax(i,1,1,l)  =   0.5d0 * (FSzzint_L + FSzzint_R - (Szint_R - Szint_L ) )

! -------------------------------------- EGLM-------------------------------------- ¡¡¡¡ seems ok just to revise the index slop for uni and multi dimensional case


       BxSxlax(i,1,1,l)   = 0.5d0 * (Bxint_L  + Bxint_R - (Sxint_R - Sxint_L ) )
       BySxlax(i,1,1,l)   = 0.5d0 * (Byint_L  + Byint_R - (Sxint_R - Sxint_L ) )
       BzSxlax(i,1,1,l)   = 0.5d0 * (Bzint_L  + Bzint_R - (Sxint_R - Sxint_L ) )



       BxSylax(i,1,1,l)   = 0.5d0 * (Bxint_L  + Bxint_R - (Syint_R - Syint_L ) )
       BySylax(i,1,1,l)   = 0.5d0 * (Byint_L  + Byint_R - (Syint_R - Syint_L ) )
       BzSylax(i,1,1,l)   = 0.5d0 * (Bzint_L  + Bzint_R - (Syint_R - Syint_L ) )


       BxSzlax(i,1,1,l)   = 0.5d0 * (Bxint_L  + Bxint_R - (Szint_R - Szint_L ) )
       BySzlax(i,1,1,l)   = 0.5d0 * (Byint_L  + Byint_R - (Szint_R - Szint_L ) )
       BzSzlax(i,1,1,l)   = 0.5d0 * (Bzint_L  + Bzint_R - (Szint_R - Szint_L ) )

! -------------------------------------- EGLM-------------------------------------- 

    end do

!$OMP END DO    

    else if (DIM == 2) then

       call boundary_conserved
       call boundary_electric
       call boundary_flux_conserved


!$OMP DO        

       do j=-1,jmax
          do i=-1,imax

       !---------------------------------
       ! Left - Right Conserved variables
       !---------------------------------

    !     Left - Right  intermediate values of augmented fields

    psiint_x_L = psiint(i  ,j  ,1,l) +  FAC  * mcl(psiint(i+1,j  ,1,l) - psiint(i  ,j  ,1,l) ,  &
                                                   psiint(i  ,j  ,1,l) - psiint(i-1,j  ,1,l))
    psiint_x_R = psiint(i+1,j  ,1,l) -  FAC  * mcl(psiint(i+2,j  ,1,l) - psiint(i+1,j  ,1,l) ,  &
                                                   psiint(i+1,j  ,1,l) - psiint(i  ,j  ,1,l))

    phiint_x_L = phiint(i  ,j  ,1,l) +  FAC  * mcl(phiint(i+1,j  ,1,l) - phiint(i  ,j  ,1,l) ,  &
                                                   phiint(i  ,j  ,1,l) - phiint(i-1,j  ,1,l))
    phiint_x_R = phiint(i+1,j  ,1,l) -  FAC  * mcl(phiint(i+2,j  ,1,l) - phiint(i+1,j  ,1,l) ,  &
                                                   phiint(i+1,j  ,1,l) - phiint(i  ,j  ,1,l))

    psiint_y_L = psiint(i  ,j  ,1,l) +  FAC  * mcl(psiint(i  ,j+1,1,l) - psiint(i  ,j  ,1,l) ,  &
                                                   psiint(i  ,j  ,1,l) - psiint(i  ,j-1,1,l))
    psiint_y_R = psiint(i  ,j+1,1,l) -  FAC  * mcl(psiint(i  ,j+2,1,l) - psiint(i  ,j+1,1,l) ,  &
                                                   psiint(i  ,j+1,1,l) - psiint(i  ,j  ,1,l))

    phiint_y_L = phiint(i  ,j  ,1,l) +  FAC  * mcl(phiint(i  ,j+1,1,l) - phiint(i  ,j  ,1,l) ,  &
                                                   phiint(i  ,j  ,1,l) - phiint(i  ,j-1,1,l))
    phiint_y_R = phiint(i  ,j+1,1,l) -  FAC  * mcl(phiint(i  ,j+2,1,l) - phiint(i  ,j+1,1,l) ,  &
                                                   phiint(i  ,j+1,1,l) - phiint(i  ,j  ,1,l))

    !     Left - Right intermediate values of magnetic field

    Bxint_x_L =  Bxint(i  ,j  ,1,l) +  FAC_SR  * mcl(Bxint(i+1,j  ,1,l) - Bxint(i  ,j  ,1,l) ,  &
                                                     Bxint(i  ,j  ,1,l) - Bxint(i-1,j  ,1,l)) 
    Bxint_x_R =  Bxint(i+1,j  ,1,l) -  FAC_SR  * mcl(Bxint(i+2,j  ,1,l) - Bxint(i+1,j  ,1,l) ,  &
                                                     Bxint(i+1,j  ,1,l) - Bxint(i  ,j  ,1,l)) 

    Byint_x_L =  Byint(i  ,j  ,1,l) +  FAC_SR  * mcl(Byint(i+1,j  ,1,l) - Byint(i  ,j  ,1,l) ,  &
                                                     Byint(i  ,j  ,1,l) - Byint(i-1,j  ,1,l))  
    Byint_x_R =  Byint(i+1,j  ,1,l) -  FAC_SR  * mcl(Byint(i+2,j  ,1,l) - Byint(i+1,j  ,1,l) ,  &
                                                     Byint(i+1,j  ,1,l) - Byint(i  ,j  ,1,l))  

    Bzint_x_L =  Bzint(i  ,j  ,1,l) +  FAC_SR  * mcl(Bzint(i+1,j  ,1,l) - Bzint(i  ,j  ,1,l) ,  &
                                                     Bzint(i  ,j  ,1,l) - Bzint(i-1,j  ,1,l))  
    Bzint_x_R =  Bzint(i+1,j  ,1,l) -  FAC_SR  * mcl(Bzint(i+2,j  ,1,l) - Bzint(i+1,j  ,1,l) ,  &
                                                     Bzint(i+1,j  ,1,l) - Bzint(i  ,j  ,1,l))  

    Bxint_y_L =  Bxint(i  ,j  ,1,l) +  FAC_SR  * mcl(Bxint(i  ,j+1,1,l) - Bxint(i  ,j  ,1,l) ,  &
                                                     Bxint(i  ,j  ,1,l) - Bxint(i  ,j-1,1,l)) 
    Bxint_y_R =  Bxint(i  ,j+1,1,l) -  FAC_SR  * mcl(Bxint(i  ,j+2,1,l) - Bxint(i  ,j+1,1,l) ,  &
                                                     Bxint(i  ,j+1,1,l) - Bxint(i  ,j  ,1,l)) 

    Byint_y_L =  Byint(i  ,j  ,1,l) +  FAC_SR  * mcl(Byint(i  ,j+1,1,l) - Byint(i  ,j  ,1,l) ,  &
                                                     Byint(i  ,j  ,1,l) - Byint(i  ,j-1,1,l))  
    Byint_y_R =  Byint(i  ,j+1,1,l) -  FAC_SR  * mcl(Byint(i  ,j+2,1,l) - Byint(i  ,j+1,1,l) ,  &
                                                     Byint(i  ,j+1,1,l) - Byint(i  ,j  ,1,l))  

    Bzint_y_L =  Bzint(i  ,j  ,1,l) +  FAC_SR  * mcl(Bzint(i  ,j+1,1,l) - Bzint(i  ,j  ,1,l) ,  &
                                                     Bzint(i  ,j  ,1,l) - Bzint(i  ,j-1,1,l))  
    Bzint_y_R =  Bzint(i  ,j+1,1,l) -  FAC_SR  * mcl(Bzint(i  ,j+2,1,l) - Bzint(i  ,j+1,1,l) ,  &
                                                     Bzint(i  ,j+1,1,l) - Bzint(i  ,j  ,1,l))  

    Bxint_L     = 0.5d0 * (Bxint_x_L + Bxint_y_L)
    Byint_L     = 0.5d0 * (Byint_x_L + Byint_y_L)
    Bzint_L     = 0.5d0 * (Bzint_x_L + Bzint_y_L)

    Bxint_R     = 0.5d0 * (Bxint_x_R + Bxint_y_R)
    Byint_R     = 0.5d0 * (Byint_x_R + Byint_y_R)
    Bzint_R     = 0.5d0 * (Bzint_x_R + Bzint_y_R)

    !-----------------------------------------------------------------------------------

    Bxint_x_L_m =  - Bxint(i  ,j  ,1,l) +  FAC  * mcl(Bxint(i+1,j  ,1,l) - Bxint(i  ,j  ,1,l) ,  &
                                                      Bxint(i  ,j  ,1,l) - Bxint(i-1,j  ,1,l)) 
    Bxint_x_R_m =  - Bxint(i+1,j  ,1,l) -  FAC  * mcl(Bxint(i+2,j  ,1,l) - Bxint(i+1,j  ,1,l) ,  &
                                                      Bxint(i+1,j  ,1,l) - Bxint(i  ,j  ,1,l)) 

    Byint_x_L_m =  - Byint(i  ,j  ,1,l) +  FAC  * mcl(Byint(i+1,j  ,1,l) - Byint(i  ,j  ,1,l) ,  &
                                                      Byint(i  ,j  ,1,l) - Byint(i-1,j  ,1,l))  
    Byint_x_R_m =  - Byint(i+1,j  ,1,l) -  FAC  * mcl(Byint(i+2,j  ,1,l) - Byint(i+1,j  ,1,l) ,  &
                                                      Byint(i+1,j  ,1,l) - Byint(i  ,j  ,1,l))  

    Bzint_x_L_m =  - Bzint(i  ,j  ,1,l) +  FAC  * mcl(Bzint(i+1,j  ,1,l) - Bzint(i  ,j  ,1,l) ,  &
                                                      Bzint(i  ,j  ,1,l) - Bzint(i-1,j  ,1,l))  
    Bzint_x_R_m =  - Bzint(i+1,j  ,1,l) -  FAC  * mcl(Bzint(i+2,j  ,1,l) - Bzint(i+1,j  ,1,l) ,  &
                                                      Bzint(i+1,j  ,1,l) - Bzint(i  ,j  ,1,l))  

    Bxint_y_L_m =  - Bxint(i  ,j  ,1,l) +  FAC  * mcl(Bxint(i  ,j+1,1,l) - Bxint(i  ,j  ,1,l) ,  &
                                                      Bxint(i  ,j  ,1,l) - Bxint(i  ,j-1,1,l)) 
    Bxint_y_R_m =  - Bxint(i  ,j+1,1,l) -  FAC  * mcl(Bxint(i  ,j+2,1,l) - Bxint(i  ,j+1,1,l) ,  &
                                                      Bxint(i  ,j+1,1,l) - Bxint(i  ,j  ,1,l)) 

    Byint_y_L_m =  - Byint(i  ,j  ,1,l) +  FAC  * mcl(Byint(i  ,j+1,1,l) - Byint(i  ,j  ,1,l) ,  &
                                                      Byint(i  ,j  ,1,l) - Byint(i  ,j-1,1,l))  
    Byint_y_R_m =  - Byint(i  ,j+1,1,l) -  FAC  * mcl(Byint(i  ,j+2,1,l) - Byint(i  ,j+1,1,l) ,  &
                                                      Byint(i  ,j+1,1,l) - Byint(i  ,j  ,1,l))  

    Bzint_y_L_m =  - Bzint(i  ,j  ,1,l) +  FAC  * mcl(Bzint(i  ,j+1,1,l) - Bzint(i  ,j  ,1,l) ,  &
                                                      Bzint(i  ,j  ,1,l) - Bzint(i  ,j-1,1,l))  
    Bzint_y_R_m =  - Bzint(i  ,j+1,1,l) -  FAC  * mcl(Bzint(i  ,j+2,1,l) - Bzint(i  ,j+1,1,l) ,  &
                                                      Bzint(i  ,j+1,1,l) - Bzint(i  ,j  ,1,l))  

    !-----------------------------------------------------------------------------------

    !     Left - Right intermediate value of electric field

    Exint_x_L =  Exint(i  ,j  ,1,l)  +  FAC_SR  * mcl(Exint(i+1,j  ,1,l) - Exint(i  ,j  ,1,l) ,  &
                                                      Exint(i  ,j  ,1,l) - Exint(i-1,j  ,1,l)) 
    Exint_x_R =  Exint(i+1,j  ,1,l)  -  FAC_SR  * mcl(Exint(i+2,j  ,1,l) - Exint(i+1,j  ,1,l) ,  &
                                                      Exint(i+1,j  ,1,l) - Exint(i  ,j  ,1,l)) 

    Eyint_x_L =  Eyint(i  ,j  ,1,l)  +  FAC_SR  * mcl(Eyint(i+1,j  ,1,l) - Eyint(i  ,j  ,1,l) ,  &
                                                      Eyint(i  ,j  ,1,l) - Eyint(i-1,j  ,1,l)) 
    Eyint_x_R =  Eyint(i+1,j  ,1,l)  -  FAC_SR  * mcl(Eyint(i+2,j  ,1,l) - Eyint(i+1,j  ,1,l) ,  &
                                                      Eyint(i+1,j  ,1,l) - Eyint(i  ,j  ,1,l)) 

    Ezint_x_L =  Ezint(i  ,j  ,1,l)  +  FAC_SR  * mcl(Ezint(i+1,j  ,1,l) - Ezint(i  ,j  ,1,l) ,  &
                                                      Ezint(i  ,j  ,1,l) - Ezint(i-1,j  ,1,l)) 
    Ezint_x_R =  Ezint(i+1,j  ,1,l)  -  FAC_SR  * mcl(Ezint(i+2,j  ,1,l) - Ezint(i+1,j  ,1,l) ,  &
                                                      Ezint(i+1,j  ,1,l) - Ezint(i  ,j  ,1,l)) 

    Exint_y_L =  Exint(i  ,j  ,1,l)  +  FAC_SR  * mcl(Exint(i  ,j+1,1,l) - Exint(i  ,j  ,1,l) ,  &
                                                      Exint(i  ,j  ,1,l) - Exint(i  ,j-1,1,l)) 
    Exint_y_R =  Exint(i  ,j+1,1,l)  -  FAC_SR  * mcl(Exint(i  ,j+2,1,l) - Exint(i  ,j+1,1,l) ,  &
                                                      Exint(i  ,j+1,1,l) - Exint(i  ,j  ,1,l)) 

    Eyint_y_L =  Eyint(i  ,j  ,1,l)  +  FAC_SR  * mcl(Eyint(i  ,j+1,1,l) - Eyint(i  ,j  ,1,l) ,  &
                                                      Eyint(i  ,j  ,1,l) - Eyint(i  ,j-1,1,l)) 
    Eyint_y_R =  Eyint(i  ,j+1,1,l)  -  FAC_SR  * mcl(Eyint(i  ,j+2,1,l) - Eyint(i  ,j+1,1,l) ,  &
                                                      Eyint(i  ,j+1,1,l) - Eyint(i  ,j  ,1,l)) 

    Ezint_y_L =  Ezint(i  ,j  ,1,l)  +  FAC_SR  * mcl(Ezint(i  ,j+1,1,l) - Ezint(i  ,j  ,1,l) ,  &
                                                      Ezint(i  ,j  ,1,l) - Ezint(i  ,j-1,1,l)) 
    Ezint_y_R =  Ezint(i  ,j+1,1,l)  -  FAC_SR  * mcl(Ezint(i  ,j+2,1,l) - Ezint(i  ,j+1,1,l) ,  &
                                                      Ezint(i  ,j+1,1,l) - Ezint(i  ,j  ,1,l)) 

    Exint_L     = 0.5d0 * (Exint_x_L + Exint_y_L)
    Eyint_L     = 0.5d0 * (Eyint_x_L + Eyint_y_L)
    Ezint_L     = 0.5d0 * (Ezint_x_L + Ezint_y_L)

    Exint_R     = 0.5d0 * (Exint_x_R + Exint_y_R)
    Eyint_R     = 0.5d0 * (Eyint_x_R + Eyint_y_R)
    Ezint_R     = 0.5d0 * (Ezint_x_R + Ezint_y_R)

    !-----------------------------------------------------------------------------------

    Exint_x_L_m =  - Exint(i  ,j  ,1,l)  +  FAC  * mcl(Exint(i+1,j  ,1,l) - Exint(i  ,j  ,1,l) ,  &
                                                       Exint(i  ,j  ,1,l) - Exint(i-1,j  ,1,l)) 
    Exint_x_R_m =  - Exint(i+1,j  ,1,l)  -  FAC  * mcl(Exint(i+2,j  ,1,l) - Exint(i+1,j  ,1,l) ,  &
                                                       Exint(i+1,j  ,1,l) - Exint(i  ,j  ,1,l)) 

    Eyint_x_L_m =  - Eyint(i  ,j  ,1,l)  +  FAC  * mcl(Eyint(i+1,j  ,1,l) - Eyint(i  ,j  ,1,l) ,  &
                                                       Eyint(i  ,j  ,1,l) - Eyint(i-1,j  ,1,l)) 
    Eyint_x_R_m =  - Eyint(i+1,j  ,1,l)  -  FAC  * mcl(Eyint(i+2,j  ,1,l) - Eyint(i+1,j  ,1,l) ,  &
                                                       Eyint(i+1,j  ,1,l) - Eyint(i  ,j  ,1,l)) 

    Ezint_x_L_m =  - Ezint(i  ,j  ,1,l)  +  FAC  * mcl(Ezint(i+1,j  ,1,l) - Ezint(i  ,j  ,1,l) ,  &
                                                       Ezint(i  ,j  ,1,l) - Ezint(i-1,j  ,1,l)) 
    Ezint_x_R_m =  - Ezint(i+1,j  ,1,l)  -  FAC  * mcl(Ezint(i+2,j  ,1,l) - Ezint(i+1,j  ,1,l) ,  &
                                                       Ezint(i+1,j  ,1,l) - Ezint(i  ,j  ,1,l)) 

    Exint_y_L_m =  - Exint(i  ,j  ,1,l)  +  FAC  * mcl(Exint(i  ,j+1,1,l) - Exint(i  ,j  ,1,l) ,  &
                                                       Exint(i  ,j  ,1,l) - Exint(i  ,j-1,1,l)) 
    Exint_y_R_m =  - Exint(i  ,j+1,1,l)  -  FAC  * mcl(Exint(i  ,j+2,1,l) - Exint(i  ,j+1,1,l) ,  &
                                                       Exint(i  ,j+1,1,l) - Exint(i  ,j  ,1,l)) 

    Eyint_y_L_m =  - Eyint(i  ,j  ,1,l)  +  FAC  * mcl(Eyint(i  ,j+1,1,l) - Eyint(i  ,j  ,1,l) ,  &
                                                       Eyint(i  ,j  ,1,l) - Eyint(i  ,j-1,1,l)) 
    Eyint_y_R_m =  - Eyint(i  ,j+1,1,l)  -  FAC  * mcl(Eyint(i  ,j+2,1,l) - Eyint(i  ,j+1,1,l) ,  &
                                                       Eyint(i  ,j+1,1,l) - Eyint(i  ,j  ,1,l)) 

    Ezint_y_L_m =  - Ezint(i  ,j  ,1,l)  +  FAC  * mcl(Ezint(i  ,j+1,1,l) - Ezint(i  ,j  ,1,l) ,  &
                                                       Ezint(i  ,j  ,1,l) - Ezint(i  ,j-1,1,l)) 
    Ezint_y_R_m =  - Ezint(i  ,j+1,1,l)  -  FAC  * mcl(Ezint(i  ,j+2,1,l) - Ezint(i  ,j+1,1,l) ,  &
                                                       Ezint(i  ,j+1,1,l) - Ezint(i  ,j  ,1,l)) 


    !-----------------------------------------------------------------------------------

    !     Left - Right Intermediate value of charge density 

    qint_x_L = qint(i  ,j  ,1,l)  +  FAC_SR  * mcl(qint(i+1,j  ,1,l) - qint(i  ,j  ,1,l) ,  &
                                                   qint(i  ,j  ,1,l) - qint(i-1,j  ,1,l)) 
    qint_x_R = qint(i+1,j  ,1,l)  -  FAC_SR  * mcl(qint(i+2,j  ,1,l) - qint(i+1,j  ,1,l) ,  &
                                                   qint(i+1,j  ,1,l) - qint(i  ,j  ,1,l)) 

    qint_y_L = qint(i  ,j  ,1,l)  +  FAC_SR  * mcl(qint(i  ,j+1,1,l) - qint(i  ,j  ,1,l) ,  &
                                                   qint(i  ,j  ,1,l) - qint(i  ,j-1,1,l)) 
    qint_y_R = qint(i  ,j+1,1,l)  -  FAC_SR  * mcl(qint(i  ,j+2,1,l) - qint(i  ,j+1,1,l) ,  &
                                                   qint(i  ,j+1,1,l) - qint(i  ,j  ,1,l))

    qint_L = 0.5d0 * (qint_x_L + qint_y_L)
    qint_R = 0.5d0 * (qint_x_R + qint_y_R)
 

    !     Mass Conservation Explicit intermediate value of D

    DDint_x_L = DDint(i  ,j  ,1,l) +  FAC  * mcl(DDint(i+1,j  ,1,l) - DDint(i  ,j  ,1,l) ,  &
                                                 DDint(i  ,j  ,1,l) - DDint(i-1,j  ,1,l)) 
    DDint_x_R = DDint(i+1,j  ,1,l) -  FAC  * mcl(DDint(i+2,j  ,1,l) - DDint(i+1,j  ,1,l) ,  &
                                                 DDint(i+1,j  ,1,l) - DDint(i  ,j  ,1,l)) 

    DDint_y_L = DDint(i  ,j  ,1,l) +  FAC  * mcl(DDint(i  ,j+1,1,l) - DDint(i  ,j  ,1,l) ,  &
                                                 DDint(i  ,j  ,1,l) - DDint(i  ,j-1,1,l)) 
    DDint_y_R = DDint(i  ,j+1,1,l) -  FAC  * mcl(DDint(i  ,j+2,1,l) - DDint(i  ,j+1,1,l) ,  &
                                                 DDint(i  ,j+1,1,l) - DDint(i  ,j  ,1,l)) 

    !     Energy Conservation Explicit intermediate value of tau

    tauint_x_L = tauint(i  ,j  ,1,l) +  FAC  * mcl(tauint(i+1,j  ,1,l) - tauint(i  ,j  ,1,l) ,  &
                                                   tauint(i  ,j  ,1,l) - tauint(i-1,j  ,1,l)) 
    tauint_x_R = tauint(i+1,j  ,1,l) -  FAC  * mcl(tauint(i+2,j  ,1,l) - tauint(i+1,j  ,1,l) ,  &
                                                   tauint(i+1,j  ,1,l) - tauint(i  ,j  ,1,l)) 

    tauint_y_L = tauint(i  ,j  ,1,l) +  FAC  * mcl(tauint(i  ,j+1,1,l) - tauint(i  ,j  ,1,l) ,  &
                                                   tauint(i  ,j  ,1,l) - tauint(i  ,j-1,1,l)) 
    tauint_y_R = tauint(i  ,j+1,1,l) -  FAC  * mcl(tauint(i  ,j+2,1,l) - tauint(i  ,j+1,1,l) ,  &
                                                   tauint(i  ,j+1,1,l) - tauint(i  ,j  ,1,l)) 

    !    Explicit intermediate value of S

    Sxint_x_L = Sxint(i  ,j  ,1,l) +  FAC  * mcl(Sxint(i+1,j  ,1,l) - Sxint(i  ,j  ,1,l) ,  &
                                                 Sxint(i  ,j  ,1,l) - Sxint(i-1,j  ,1,l)) 
    Sxint_x_R = Sxint(i+1,j  ,1,l) -  FAC  * mcl(Sxint(i+2,j  ,1,l) - Sxint(i+1,j  ,1,l) ,  &
                                                 Sxint(i+1,j  ,1,l) - Sxint(i  ,j  ,1,l)) 

    Syint_x_L = Syint(i  ,j  ,1,l) +  FAC  * mcl(Syint(i+1,j  ,1,l) - Syint(i  ,j  ,1,l) ,  &
                                                 Syint(i  ,j  ,1,l) - Syint(i-1,j  ,1,l))
    Syint_x_R = Syint(i+1,j  ,1,l) -  FAC  * mcl(Syint(i+2,j  ,1,l) - Syint(i+1,j  ,1,l) ,  &
                                                 Syint(i+1,j  ,1,l) - Syint(i  ,j  ,1,l))

    Szint_x_L = Szint(i  ,j  ,1,l) +  FAC  * mcl(Szint(i+1,j  ,1,l) - Szint(i  ,j  ,1,l) ,  &
                                                 Szint(i  ,j  ,1,l) - Szint(i-1,j  ,1,l))
    Szint_x_R = Szint(i+1,j  ,1,l) -  FAC  * mcl(Szint(i+2,j  ,1,l) - Szint(i+1,j  ,1,l) ,  &
                                                 Szint(i+1,j  ,1,l) - Szint(i  ,j  ,1,l))

    Sxint_y_L = Sxint(i  ,j  ,1,l) +  FAC  * mcl(Sxint(i  ,j+1,1,l) - Sxint(i  ,j  ,1,l) ,  &
                                                 Sxint(i  ,j  ,1,l) - Sxint(i  ,j-1,1,l)) 
    Sxint_y_R = Sxint(i  ,j+1,1,l) -  FAC  * mcl(Sxint(i  ,j+2,1,l) - Sxint(i  ,j+1,1,l) ,  &
                                                 Sxint(i  ,j+1,1,l) - Sxint(i  ,j  ,1,l)) 

    Syint_y_L = Syint(i  ,j  ,1,l) +  FAC  * mcl(Syint(i  ,j+1,1,l) - Syint(i  ,j  ,1,l) ,  &
                                                 Syint(i  ,j  ,1,l) - Syint(i  ,j-1,1,l))
    Syint_y_R = Syint(i  ,j+1,1,l) -  FAC  * mcl(Syint(i  ,j+2,1,l) - Syint(i  ,j+1,1,l) ,  &
                                                 Syint(i  ,j+1,1,l) - Syint(i  ,j  ,1,l))

    Szint_y_L = Szint(i  ,j  ,1,l) +  FAC  * mcl(Szint(i  ,j+1,1,l) - Szint(i  ,j  ,1,l) ,  &
                                                 Szint(i  ,j  ,1,l) - Szint(i  ,j-1,1,l))
    Szint_y_R = Szint(i  ,j+1,1,l) -  FAC  * mcl(Szint(i  ,j+2,1,l) - Szint(i  ,j+1,1,l) ,  &
                                                 Szint(i  ,j+1,1,l) - Szint(i  ,j  ,1,l))


       !---------------------------------
       ! Left - Right Fluxes
       !---------------------------------

    !     Electric and Magnetic field modules 


    E2B2int_x_L = Exint_x_L**2 + Eyint_x_L**2 + Ezint_x_L**2 + &
                  Bxint_x_L**2 + Byint_x_L**2 + Bzint_x_L**2

    E2B2int_x_R = Exint_x_R**2 + Eyint_x_R**2 + Ezint_x_R**2 + &
                  Bxint_x_R**2 + Byint_x_R**2 + Bzint_x_R**2

    E2B2int_y_L = Exint_y_L**2 + Eyint_y_L**2 + Ezint_y_L**2 + &
                  Bxint_y_L**2 + Byint_y_L**2 + Bzint_y_L**2

    E2B2int_y_R = Exint_y_R**2 + Eyint_y_R**2 + Ezint_y_R**2 + &
                  Bxint_y_R**2 + Byint_y_R**2 + Bzint_y_R**2


    E2B2int_L = Exint_L**2 + Eyint_L**2 + Ezint_L**2 + &
                Bxint_L**2 + Byint_L**2 + Bzint_L**2

    E2B2int_R = Exint_R**2 + Eyint_R**2 + Ezint_R**2 + &
                Bxint_R**2 + Byint_R**2 + Bzint_R**2






          if (REC_PRIM == 0) then

             Vx_x_L = Vx(i  ,j  ,1) +  FAC_SR  * mcl(Vx(i+1,j  ,1) - Vx(i  ,j  ,1) ,  &
                                                     Vx(i  ,j  ,1) - Vx(i-1,j  ,1))
             Vx_x_R = Vx(i+1,j  ,1) -  FAC_SR  * mcl(Vx(i+2,j  ,1) - Vx(i+1,j  ,1) ,  &
                                                     Vx(i+1,j  ,1) - Vx(i  ,j  ,1))

             Vy_x_L = Vy(i  ,j  ,1) +  FAC_SR  * mcl(Vy(i+1,j  ,1) - Vy(i  ,j  ,1) ,  &
                                                     Vy(i  ,j  ,1) - Vy(i-1,j  ,1))
             Vy_x_R = Vy(i+1,j  ,1) -  FAC_SR  * mcl(Vy(i+2,j  ,1) - Vy(i+1,j  ,1) ,  &
                                                     Vy(i+1,j  ,1) - Vy(i  ,j  ,1))

             Vz_x_L = Vz(i  ,j  ,1) +  FAC_SR  * mcl(Vz(i+1,j  ,1) - Vz(i  ,j  ,1) ,  &
                                                     Vz(i  ,j  ,1) - Vz(i-1,j  ,1))
             Vz_x_R = Vz(i+1,j  ,1) -  FAC_SR  * mcl(Vz(i+2,j  ,1) - Vz(i+1,j  ,1) ,  &
                                                     Vz(i+1,j  ,1) - Vz(i  ,j  ,1))

             Vx_y_L = Vx(i  ,j  ,1) +  FAC_SR  * mcl(Vx(i  ,j+1,1) - Vx(i  ,j  ,1) ,  &
                                                     Vx(i  ,j  ,1) - Vx(i  ,j-1,1))
             Vx_y_R = Vx(i  ,j+1,1) -  FAC_SR  * mcl(Vx(i  ,j+2,1) - Vx(i  ,j+1,1) ,  &
                                                     Vx(i  ,j+1,1) - Vx(i  ,j  ,1))

             Vy_y_L = Vy(i  ,j  ,1) +  FAC_SR  * mcl(Vy(i  ,j+1,1) - Vy(i  ,j  ,1) ,  &
                                                     Vy(i  ,j  ,1) - Vy(i  ,j-1,1))
             Vy_y_R = Vy(i  ,j+1,1) -  FAC_SR  * mcl(Vy(i  ,j+2,1) - Vy(i  ,j+1,1) ,  &
                                                     Vy(i  ,j+1,1) - Vy(i  ,j  ,1))

             Vz_y_L = Vz(i  ,j  ,1) +  FAC_SR  * mcl(Vz(i  ,j+1,1) - Vz(i  ,j  ,1) ,  &
                                                     Vz(i  ,j  ,1) - Vz(i  ,j-1,1))
             Vz_y_R = Vz(i  ,j+1,1) -  FAC_SR  * mcl(Vz(i  ,j+2,1) - Vz(i  ,j+1,1) ,  &
                                                     Vz(i  ,j+1,1) - Vz(i  ,j  ,1))


             Vx_x = 0.5d0 * (Vx_x_L + Vx_x_R)
             Vx_y = 0.5d0 * (Vx_y_L + Vx_y_R)

             Vy_x = 0.5d0 * (Vy_x_L + Vy_x_R)
             Vy_y = 0.5d0 * (Vy_y_L + Vy_y_R)

             Vz_x = 0.5d0 * (Vz_x_L + Vz_x_R)
             Vz_y = 0.5d0 * (Vz_y_L + Vz_y_R)

             Vx_L = 0.5d0 * (Vx_x_L + Vx_y_L)
             Vy_L = 0.5d0 * (Vy_x_L + Vy_y_L)
             Vz_L = 0.5d0 * (Vz_x_L + Vz_y_L)

             Vx_R = 0.5d0 * (Vx_x_R + Vx_y_R)
             Vy_R = 0.5d0 * (Vy_x_R + Vy_y_R)
             Vz_R = 0.5d0 * (Vz_x_R + Vz_y_R)

             Vx_rec = 0.5d0 * (Vx_L + Vx_R)
             Vy_rec = 0.5d0 * (Vy_L + Vy_R)
             Vz_rec = 0.5d0 * (Vz_L + Vz_R)


             W_x_L = 1.d0/sqrt(1.d0 - Vx_x_L**2 -  Vy_x_L**2 - Vz_x_L**2)
             W_x_R = 1.d0/sqrt(1.d0 - Vx_x_R**2 -  Vy_x_R**2 - Vz_x_R**2)


             W_y_L = 1.d0/sqrt(1.d0 - Vx_y_L**2 -  Vy_y_L**2 - Vz_y_L**2)
             W_y_R = 1.d0/sqrt(1.d0 - Vx_y_R**2 -  Vy_y_R**2 - Vz_y_R**2)

             W_L = 1.d0/sqrt(1.d0 - Vx_L**2 -  Vy_L**2 - Vz_L**2)
             W_R = 1.d0/sqrt(1.d0 - Vx_R**2 -  Vy_R**2 - Vz_R**2)


          else if (REC_PRIM == 1) then

             Vx_x_L = Vxint(i  ,j  ,1,l) +  FAC_SR  * mcl(Vxint(i+1,j  ,1,l) - Vxint(i  ,j  ,1,l) ,  &
                                                          Vxint(i  ,j  ,1,l) - Vxint(i-1,j  ,1,l))
             Vx_x_R = Vxint(i+1,j  ,1,l) -  FAC_SR  * mcl(Vxint(i+2,j  ,1,l) - Vxint(i+1,j  ,1,l) ,  &
                                                          Vxint(i+1,j  ,1,l) - Vxint(i  ,j  ,1,l))

             Vy_x_L = Vyint(i  ,j  ,1,l) +  FAC_SR  * mcl(Vyint(i+1,j  ,1,l) - Vyint(i  ,j  ,1,l) ,  &
                                                          Vyint(i  ,j  ,1,l) - Vyint(i-1,j  ,1,l))
             Vy_x_R = Vyint(i+1,j  ,1,l) -  FAC_SR  * mcl(Vyint(i+2,j  ,1,l) - Vyint(i+1,j  ,1,l) ,  &
                                                          Vyint(i+1,j  ,1,l) - Vyint(i  ,j  ,1,l))

             Vz_x_L = Vzint(i  ,j  ,1,l) +  FAC_SR  * mcl(Vzint(i+1,j  ,1,l) - Vzint(i  ,j  ,1,l) ,  &
                                                          Vzint(i  ,j  ,1,l) - Vzint(i-1,j  ,1,l))
             Vy_x_R = Vzint(i+1,j  ,1,l) -  FAC_SR  * mcl(Vzint(i+2,j  ,1,l) - Vzint(i+1,j  ,1,l) ,  &
                                                          Vzint(i+1,j  ,1,l) - Vzint(i  ,j  ,1,l))


             Vx_y_L = Vxint(i  ,j  ,1,l) +  FAC_SR  * mcl(Vxint(i  ,j+1,1,l) - Vxint(i  ,j  ,1,l) ,  &
                                                          Vxint(i  ,j  ,1,l) - Vxint(i  ,j-1,1,l))
             Vx_y_R = Vxint(i  ,j+1,1,l) -  FAC_SR  * mcl(Vxint(i  ,j+2,1,l) - Vxint(i  ,j+1,1,l) ,  &
                                                          Vxint(i  ,j+1,1,l) - Vxint(i  ,j  ,1,l))

             Vy_y_L = Vyint(i  ,j  ,1,l) +  FAC_SR  * mcl(Vyint(i  ,j+1,1,l) - Vyint(i  ,j  ,1,l) ,  &
                                                          Vyint(i  ,j  ,1,l) - Vyint(i  ,j-1,1,l))
             Vy_y_R = Vyint(i  ,j+1,1,l) -  FAC_SR  * mcl(Vyint(i  ,j+2,1,l) - Vyint(i  ,j+1,1,l) ,  &
                                                          Vyint(i  ,j+1,1,l) - Vyint(i  ,j  ,1,l))

             Vz_y_L = Vzint(i  ,j  ,1,l) +  FAC_SR  * mcl(Vzint(i  ,j+1,1,l) - Vzint(i  ,j  ,1,l) ,  &
                                                          Vzint(i  ,j  ,1,l) - Vzint(i  ,j-1,1,l))
             Vy_y_R = Vzint(i  ,j+1,1,l) -  FAC_SR  * mcl(Vzint(i  ,j+2,1,l) - Vzint(i  ,j+1,1,l) ,  &
                                                          Vzint(i  ,j+1,1,l) - Vzint(i  ,j  ,1,l))



             Vx_x = 0.5d0 * (Vx_x_L + Vx_x_R)
             Vx_y = 0.5d0 * (Vx_y_L + Vx_y_R)

             Vy_x = 0.5d0 * (Vy_x_L + Vy_x_R)
             Vy_y = 0.5d0 * (Vy_y_L + Vy_y_R)

             Vz_x = 0.5d0 * (Vz_x_L + Vz_x_R)
             Vz_y = 0.5d0 * (Vz_y_L + Vz_y_R)

             Vx_L = 0.5d0 * (Vx_x_L + Vx_y_L)
             Vy_L = 0.5d0 * (Vy_x_L + Vy_y_L)
             Vz_L = 0.5d0 * (Vz_x_L + Vz_y_L)

             Vx_R = 0.5d0 * (Vx_x_R + Vx_y_R)
             Vy_R = 0.5d0 * (Vy_x_R + Vy_y_R)
             Vz_R = 0.5d0 * (Vz_x_R + Vz_y_R)

             Vx_L = 0.5d0 * (Vx_x_L + Vx_y_L)
             Vy_L = 0.5d0 * (Vy_x_L + Vy_y_L)
             Vz_L = 0.5d0 * (Vz_x_L + Vz_y_L)

             Vx_R = 0.5d0 * (Vx_x_R + Vx_y_R)
             Vy_R = 0.5d0 * (Vy_x_R + Vy_y_R)
             Vz_R = 0.5d0 * (Vz_x_R + Vz_y_R)


             W_x_L = 1.d0/sqrt(1.d0 - Vx_x_L**2 -  Vy_x_L**2 - Vz_x_L**2)
             W_x_R = 1.d0/sqrt(1.d0 - Vx_x_R**2 -  Vy_x_R**2 - Vz_x_R**2)


             W_y_L = 1.d0/sqrt(1.d0 - Vx_y_L**2 -  Vy_y_L**2 - Vz_y_L**2)
             W_y_R = 1.d0/sqrt(1.d0 - Vx_y_R**2 -  Vy_y_R**2 - Vz_y_R**2)


             W_L = 1.d0/sqrt(1.d0 - Vx_L**2 -  Vy_L**2 - Vz_L**2)
             W_R = 1.d0/sqrt(1.d0 - Vx_R**2 -  Vy_R**2 - Vz_R**2)


          else

             write(*,*) "STOP: subroutine hllflux"
             write(*,*) "REC_PRIM parameter is not valid"
             stop

          end if

             Edotv_L = Exint_L * Vx_L + Eyint_L * Vy_L + Ezint_L * Vz_L
             Edotv_R = Exint_R * Vx_R + Eyint_R * Vy_R + Ezint_R * Vz_R

             Edotv_x_L = Exint_x_L * Vx_x_L + Eyint_x_L * Vy_x_L + Ezint_x_L * Vz_x_L
             Edotv_x_R = Exint_x_R * Vx_x_R + Eyint_x_R * Vy_x_R + Ezint_x_R * Vz_x_R

             Edotv_y_L = Exint_y_L * Vx_y_L + Eyint_y_L * Vy_y_L + Ezint_y_L * Vz_y_L
             Edotv_y_R = Exint_y_R * Vx_y_R + Eyint_y_R * Vy_y_R + Ezint_y_R * Vz_y_R


             vrotB_x_L = (Bzint_L * Vy_L - Byint_L * Vz_L)
             vrotB_y_L = (Bxint_L * Vz_L - Bzint_L * Vx_L)
             vrotB_z_L = (Byint_L * Vx_L - Bxint_L * Vy_L)

             vrotB_x_R = (Bzint_R * Vy_R - Byint_R * Vz_R)
             vrotB_y_R = (Bxint_R * Vz_R - Bzint_R * Vx_R)
             vrotB_z_R = (Byint_R * Vx_R - Bxint_R * Vy_R)

             vrotB_x_x_L = (Bzint_x_L * Vy_x_L - Byint_x_L * Vz_x_L)
             vrotB_y_x_L = (Bxint_x_L * Vz_x_L - Bzint_x_L * Vx_x_L)
             vrotB_z_x_L = (Byint_x_L * Vx_x_L - Bxint_x_L * Vy_x_L)

             vrotB_x_x_R = (Bzint_x_R * Vy_x_R - Byint_x_R * Vz_x_R)
             vrotB_y_x_R = (Bxint_x_R * Vz_x_R - Bzint_x_R * Vx_x_R)
             vrotB_z_x_R = (Byint_x_R * Vx_x_R - Bxint_x_R * Vy_x_R)

             vrotB_x_y_L = (Bzint_y_L * Vy_y_L - Byint_y_L * Vz_y_L)
             vrotB_y_y_L = (Bxint_y_L * Vz_y_L - Bzint_y_L * Vx_y_L)
             vrotB_z_y_L = (Byint_y_L * Vx_y_L - Bxint_y_L * Vy_y_L)

             vrotB_x_y_R = (Bzint_y_R * Vy_y_R - Byint_y_R * Vz_y_R)
             vrotB_y_y_R = (Bxint_y_R * Vz_y_R - Bzint_y_R * Vx_y_R)
             vrotB_z_y_R = (Byint_y_R * Vx_y_R - Bxint_y_R * Vy_y_R)



              if (SLC == 1) sigma = sigma_loc(i,j,1)


    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxint_L  = sigma * W_x_L * ( Exint_x_L +  vrotB_x_x_L  - Edotv_x_L * Vx_x_L )   &
                             +   qint_x_L  *  Vx_x_L
    Jxint_R  = sigma * W_x_R * ( Exint_x_R +  vrotB_x_x_R  - Edotv_x_R * Vx_x_R )   &
                             +   qint_x_R  *  Vx_x_R

    Jyint_L = sigma * W_y_L  * ( Eyint_y_L +  vrotB_y_y_L  - Edotv_y_L * Vy_y_L )   &
                             +   qint_y_L  *  Vy_y_L
    Jyint_R = sigma * W_y_R  * ( Eyint_y_R +  vrotB_y_y_R  - Edotv_y_R * Vy_y_R )   &
                             +   qint_y_R  *  Vy_y_R

    Jzint_L = sigma * W_x_L  * ( Ezint_x_L +  vrotB_z_x_L  - Edotv_x_L * Vz_x_L )   &
                             +   qint_x_L  *  Vz_x_L
    Jzint_R = sigma * W_x_R  * ( Ezint_x_R +  vrotB_z_x_R  - Edotv_x_R * Vz_x_R )   &
                             +   qint_x_R  *  Vz_x_R

    !     Density flux (FD = \rho W V)

    FDxint_L = DDint_x_L * Vx_x_L
    FDxint_R = DDint_x_R * Vx_x_R

    FDyint_L = DDint_y_L * Vy_y_L
    FDyint_R = DDint_y_R * Vy_y_R

    FDzint_L = DDint_x_L * Vz_x_L
    FDzint_R = DDint_x_R * Vz_x_R !---> in 2D we are using x reconstruction component
                                  !---> this must be zero when realize de flux diference in z component

              !     Entalphy 


    if (REC_PRIM == 0 ) then

    p_x_L   = p  (i  ,j  ,1) + FAC  * mcl(p  (i+1,j  ,1) - p  (i  ,j  ,1) ,  &
                                          p  (i  ,j  ,1) - p  (i-1,j  ,1))
    p_x_R   = p  (i+1,j  ,1) - FAC  * mcl(p  (i+2,j  ,1) - p  (i+1,j  ,1) ,  &
                                          p  (i+1,j  ,1) - p  (i  ,j  ,1))

    rho_x_L = rho(i  ,j  ,1) + FAC  * mcl(rho(i+1,j  ,1) - rho(i  ,j  ,1) ,  &
                                          rho(i  ,j  ,1) - rho(i-1,j  ,1))
    rho_x_R = rho(i+1,j  ,1) - FAC  * mcl(rho(i+2,j  ,1) - rho(i+1,j  ,1) ,  &
                                          rho(i+1,j  ,1) - rho(i  ,j  ,1))

    p_y_L   = p  (i  ,j  ,1) + FAC  * mcl(p  (i  ,j+1,1) - p  (i  ,j  ,1) ,  &
                                          p  (i  ,j  ,1) - p  (i  ,j-1,1))
    p_y_R   = p  (i  ,j+1,1) - FAC  * mcl(p  (i  ,j+2,1) - p  (i  ,j+1,1) ,  &
                                          p  (i  ,j+1,1) - p  (i  ,j  ,1))

    rho_y_L = rho(i  ,j  ,1) + FAC  * mcl(rho(i  ,j+1,1) - rho(i  ,j  ,1) ,  &
                                          rho(i  ,j  ,1) - rho(i  ,j-1,1))
    rho_y_R = rho(i  ,j+1,1) - FAC  * mcl(rho(i  ,j+2,1) - rho(i  ,j+1,1) ,  &
                                          rho(i  ,j+1,1) - rho(i  ,j  ,1))

    else if (REC_PRIM == 1 ) then

    p_x_L   = pint(i  ,j  ,1,l)     +  FAC  * mcl(pint(i+1,j  ,1,l) - pint(i  ,j  ,1,l) ,  &
                                                  pint(i  ,j  ,1,l) - pint(i-1,j  ,1,l))
    p_x_R   = pint(i+1,j  ,1,l)     -  FAC  * mcl(pint(i+2,j  ,1,l) - pint(i+1,j  ,1,l) ,  &
                                                  pint(i+1,j  ,1,l) - pint(i  ,j  ,1,l))
    p_y_L   = pint(i  ,j  ,1,l)     +  FAC  * mcl(pint(i  ,j+1,1,l) - pint(i  ,j  ,1,l) ,  &
                                                  pint(i  ,j  ,1,l) - pint(i  ,j-1,1,l))
    p_y_R   = pint(i  ,j+1,1,l)     -  FAC  * mcl(pint(i  ,j+2,1,l) - pint(i  ,j+1,1,l) ,  &
                                                  pint(i  ,j+1,1,l) - pint(i  ,j  ,1,l))

    rho_x_L   = rhoint(i  ,j  ,1,l) +  FAC  * mcl(rhoint(i+1,j  ,1,l) - rhoint(i  ,j  ,1,l) ,  &
                                                  rhoint(i  ,j  ,1,l) - rhoint(i-1,j  ,1,l))
    rho_x_R   = rhoint(i+1,j  ,1,l) -  FAC  * mcl(rhoint(i+2,j  ,1,l) - rhoint(i+1,j  ,1,l) ,  &
                                                  rhoint(i+1,j  ,1,l) - rhoint(i  ,j  ,1,l))
    rho_y_L   = rhoint(i  ,j  ,1,l) +  FAC  * mcl(rhoint(i  ,j+1,1,l) - rhoint(i  ,j  ,1,l) ,  &
                                                  rhoint(i  ,j  ,1,l) - rhoint(i  ,j-1,1,l))
    rho_y_R   = rhoint(i  ,j+1,1,l) -  FAC  * mcl(rhoint(i  ,j+2,1,l) - rhoint(i  ,j+1,1,l) ,  &
                                                  rhoint(i  ,j+1,1,l) - rhoint(i  ,j  ,1,l))

     else

        write(*,*) "STOP: subroutine hllflux"
        write(*,*) "REC_PRIM parameter is not valid"
        stop

     end if

    p_L   = 0.5d0 * (p_x_L   + p_y_L  )
    p_R   = 0.5d0 * (p_x_R   + p_y_R  )

!!$    p_rec = 0.5d0 * (p_L + p_R)

    rho_L = 0.5d0 * (rho_x_L + rho_y_L)
    rho_R = 0.5d0 * (rho_x_R + rho_y_R)

!!$    rho_rec = 0.5d0 * (rho_L + rho_R)


    epsilon_x_L  = p_x_L   / ((gamma-1.d0) * rho_x_L)
    epsilon_x_R  = p_x_R   / ((gamma-1.d0) * rho_x_R)

    epsilon_x    = 0.5d0 * (epsilon_x_L + epsilon_x_R)

    enthpy_x_L   = rho_x_L * ( 1.d0 + epsilon_x_L ) + p_x_L 
    enthpy_x_R   = rho_x_R * ( 1.d0 + epsilon_x_R ) + p_x_R 

    epsilon_y_L  = p_y_L   / ((gamma-1.d0) * rho_y_L)
    epsilon_y_R  = p_y_R   / ((gamma-1.d0) * rho_y_R)

    epsilon_y    = 0.5d0 * (epsilon_y_L + epsilon_y_R)

    enthpy_y_L   = rho_y_L * ( 1.d0 + epsilon_y_L ) + p_y_L 
    enthpy_y_R   = rho_y_R * ( 1.d0 + epsilon_y_R ) + p_y_R 


    epsilon_L   = p_L      / ((gamma-1.d0) * rho_L)
    epsilon_R   = p_R      / ((gamma-1.d0) * rho_R)

    epsilon_rec = 0.5d0 * (epsilon_L + epsilon_R)

    enthpy_L    =  rho_L   * ( 1.d0 + epsilon_L ) + p_L 
    enthpy_R    =  rho_R   * ( 1.d0 + epsilon_R ) + p_R 


!!$    if (p_x_L .lt. 0 .or. p_x_R .lt. 0) then 
!!$
!!$       print*, "negative pressure in HLL ---> p_L =", p_x_L, "p_R =", p_x_R, "i,j =", i, j
!!$
!!$    end if

!!$    if (rho_x_L .lt. 0 .or. rho_x_R .lt. 0) then 
!!$
!!$       print*, "negative density in HLL ---> rho_L =", rho_L, "rho_R =", rho_R, "i =", i, j
!!$
!!$    end if
    


!!$       enthpy_x_L   =  enthpyint(i  ,j  ,1,l) +  FAC_SR  * mcl(enthpyint(i+1,j  ,1,l) - enthpyint(i  ,j  ,1,l) ,  &
!!$                                                               enthpyint(i  ,j  ,1,l) - enthpyint(i-1,j  ,1,l))
!!$       enthpy_x_R   =  enthpyint(i+1,j  ,1,l) -  FAC_SR  * mcl(enthpyint(i+2,j  ,1,l) - enthpyint(i+1,j  ,1,l) ,  &
!!$                                                               enthpyint(i+1,j  ,1,l) - enthpyint(i  ,j  ,1,l))
!!$
!!$       enthpy_y_L   =  enthpyint(i  ,j  ,1,l) +  FAC_SR  * mcl(enthpyint(i  ,j+1,1,l) - enthpyint(i  ,j  ,1,l) ,  &
!!$                                                               enthpyint(i  ,j  ,1,l) - enthpyint(i  ,j-1,1,l))
!!$       enthpy_y_R   =  enthpyint(i  ,j+1,1,l) -  FAC_SR  * mcl(enthpyint(i  ,j+2,1,l) - enthpyint(i  ,j+1,1,l) ,  &
!!$                                                               enthpyint(i  ,j+1,1,l) - enthpyint(i  ,j  ,1,l))


       enthpy_L   =  0.5d0 * ( enthpy_x_L + enthpy_y_L )
       enthpy_R   =  0.5d0 * ( enthpy_x_R + enthpy_y_R )


    !     Flujos conservados de energía

    Ftauxint_L = (Bzint_x_L  * Eyint_x_L - Byint_x_L * Ezint_x_L) +   &
                  enthpy_x_L * W_x_L**2  * Vx_x_L   
    Ftauxint_R = (Bzint_x_R  * Eyint_x_R - Byint_x_R * Ezint_x_R) +   &
                  enthpy_x_R * W_x_R**2  * Vx_x_R   

    Ftauyint_L = (Bxint_y_L  * Ezint_y_L - Bzint_y_L * Exint_y_L) +   &
                  enthpy_y_L * W_y_L**2  * Vy_y_L
    Ftauyint_R = (Bxint_y_R  * Ezint_y_R - Bzint_y_R * Exint_y_R) +   &
                  enthpy_y_R * W_y_R**2  * Vy_y_R

    Ftauzint_L = (Byint_x_L  * Exint_x_L - Bxint_x_L * Eyint_x_L) +   &
                  enthpy_x_L * W_x_L**2  * Vz_x_L
    Ftauzint_R = (Byint_x_R  * Exint_x_R - Bxint_x_R * Eyint_x_R) +   &
                  enthpy_x_R * W_x_R**2  * Vz_x_R

    !     Components of the flux momentum tensor 

    FSxxint_L =  - Exint_x_L**2 - Bxint_x_L**2 + enthpy_x_L * W_x_L**2 * Vx_x_L**2 + &
                   0.5d0 * E2B2int_x_L + p_x_L  
    FSxxint_R =  - Exint_x_R**2 - Bxint_x_R**2 + enthpy_x_R * W_x_R**2 * Vx_x_R**2 + &
                   0.5d0 * E2B2int_x_R + p_x_R  

    FSxyint_L = -  Exint_x_L  * Eyint_x_L - Bxint_x_L * Byint_x_L +    &
                   enthpy_x_L * W_x_L**2  * Vx_x_L    * Vy_x_L
    FSxyint_R = -  Exint_x_R  * Eyint_x_R - Bxint_x_R * Byint_x_R +    &
                   enthpy_x_R * W_x_R**2  * Vx_x_R    * Vy_x_R

    FSxzint_L = -  Exint_x_L  * Ezint_x_L - Bxint_x_L * Bzint_x_L +    &
                   enthpy_x_L * W_x_L**2  * Vx_x_L    * Vz_x_L
    FSxzint_R = -  Exint_x_R  * Ezint_x_R - Bxint_x_R * Bzint_x_R +    &
                   enthpy_x_R * W_x_R**2  * Vx_x_R    * Vz_x_R

    FSyxint_L = -  Exint_y_L  * Eyint_y_L - Bxint_y_L * Byint_y_L +    &
                   enthpy_y_L * W_y_L**2  * Vx_y_L    * Vy_y_L
    FSyxint_R = -  Exint_y_R  * Eyint_y_R - Bxint_y_R * Byint_y_R +    &
                   enthpy_y_R * W_y_R**2  * Vx_y_R    * Vy_y_R

    FSyyint_L = -  Eyint_y_L**2 - Byint_y_L**2 + enthpy_y_L * W_y_L**2 * Vy_y_L**2 + &
                   0.5d0 * E2B2int_y_L  + p_y_L
    FSyyint_R = -  Eyint_y_R**2 - Byint_y_R**2 + enthpy_y_R * W_y_R**2 * Vy_y_R**2 + &
                   0.5d0 * E2B2int_y_R  + p_y_R

    FSyzint_L = -  Eyint_y_L  * Ezint_y_L - Byint_y_L * Bzint_y_L +    &
                   enthpy_y_L * W_y_L**2  * Vy_y_L    * Vz_y_L
    FSyzint_R = -  Eyint_y_R  * Ezint_y_R - Byint_y_R * Bzint_y_R +    &
                   enthpy_y_R * W_y_R**2  * Vy_y_R    * Vz_y_R

    FSzxint_L = -  Exint_x_L  * Ezint_x_L - Bxint_x_L * Bzint_x_L +    &
                   enthpy_x_L * W_x_L**2  * Vx_x_L    * Vz_x_L       
    FSzxint_R = -  Exint_x_R  * Ezint_x_R - Bxint_x_R * Bzint_x_R +    &
                   enthpy_x_R * W_x_R**2  * Vx_x_R    * Vz_x_R             

    FSzyint_L = -  Eyint_x_L  * Ezint_x_L - Byint_x_L * Bzint_x_L +    &
                   enthpy_x_L * W_x_L**2  * Vy_x_L    * Vz_x_L
    FSzyint_R = -  Eyint_x_R  * Ezint_x_R - Byint_x_R * Bzint_x_R +    &
                   enthpy_x_R * W_x_R**2  * Vy_x_R    * Vz_x_R

    FSzzint_L = -  Ezint_x_L**2 - Bzint_x_L**2 + enthpy_x_L * W_x_L**2 * Vz_x_L**2 + &
                   0.5d0 * E2B2int_x_L + p_x_L
    FSzzint_R = -  Ezint_x_R**2 - Bzint_x_R**2 + enthpy_x_R * W_x_R**2 * Vz_x_R**2 + &
                   0.5d0 * E2B2int_x_R + p_x_R


          if (const_sp .eq. 1.d0) then

             c_sp = (gamma * (gamma -1.d0) * epsiln(i,j,1))/(1.d0 + gamma * epsiln(i,j,1))
             Sp   = const_sp * abs((sqrt(Vx(i,j,1)**2 + Vy(i,j,1)**2 + Vz(i,j,1)**2) + c_sp))

!!$             c_sp = (gamma * (gamma -1.d0) * epsiln_rec)/(1.d0 + gamma * epsiln_rec)
!!$             Sp   = const_sp * abs((sqrt(Vx_rec**2 + Vy_rec**2 + Vz_rec**2) + c_sp))

!!$             c_sp_x = (gamma * (gamma -1.d0) * epsiln_x)/(1.d0 + gamma * epsiln_x)
!!$             Sp_x   = const_sp * abs((sqrt(Vx_x**2 + Vy_x**2 + Vz_x**2) + c_sp_x))
!!$
!!$             c_sp_y = (gamma * (gamma -1.d0) * epsiln_y)/(1.d0 + gamma * epsiln_y)
!!$             Sp_y   = const_sp * abs((sqrt(Vx_y**2 + Vy_y**2 + Vz_y**2) + c_sp_y))
          else

             Sp   = 1.d0
!!$             Sp_x = 1.d0
!!$             Sp_y = 1.d0

          end if


! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,j,1,l)  =   0.5d0 * (Exint_x_L + Exint_x_R - Sp * ( psiint_x_R - psiint_x_L ) )  
       Eylaxpsi(i,j,1,l)  =   0.5d0 * (Eyint_y_L + Eyint_y_R - Sp * ( psiint_y_R - psiint_y_L ) ) 
       Ezlaxpsi(i,j,1,l)  =   0.5d0 * (Ezint_x_L + Ezint_x_R - Sp * ( psiint_x_R - psiint_x_L ) ) 

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,j,1,l)  =   0.5d0 * (Bxint_x_L + Bxint_x_R - Sp * ( phiint_x_R - phiint_x_L ) )
       Bylaxphi(i,j,1,l)  =   0.5d0 * (Byint_y_L + Byint_y_R - Sp * ( phiint_y_R - phiint_y_L ) )
       Bzlaxphi(i,j,1,l)  =   0.5d0 * (Bzint_x_L + Bzint_x_R - Sp * ( phiint_x_R - phiint_x_L ) )

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,j,1,l)  =   0.5d0 * (phiint_x_L   + phiint_x_R   - Sp * ( Bxint_x_R - Bxint_x_L ) )
       EzlaxBx (i,j,1,l)  =   0.5d0 * (Ezint_y_L    + Ezint_y_R    - Sp * ( Bxint_y_R - Bxint_y_L ) )
       EylaxBx (i,j,1,l)  = - 0.5d0 * (Eyint_x_L_m  + Eyint_x_R_m  - Sp * ( Bxint_x_R - Bxint_x_L ) )

! Flujos By

       EzlaxBy (i,j,1,l)  = - 0.5d0 * (Ezint_x_L_m  + Ezint_x_R_m  - Sp * ( Byint_x_R - Byint_x_L ) )
       philaxBy(i,j,1,l)  =   0.5d0 * (phiint_y_L   + phiint_y_R   - Sp * ( Byint_y_R - Byint_y_L ) )
       ExlaxBy (i,j,1,l)  =   0.5d0 * (Exint_x_L    + Exint_x_R    - Sp * ( Byint_x_R - Byint_x_L ) )

! Flujos Bz

       EylaxBz (i,j,1,l)  =   0.5d0 * (Eyint_x_L    + Eyint_x_R    - Sp * ( Bzint_x_R - Bzint_x_L ) )
       ExlaxBz (i,j,1,l)  = - 0.5d0 * (Exint_y_L_m  + Exint_y_R_m  - Sp * ( Bzint_y_R - Bzint_y_L ) )
       philaxBz(i,j,1,l)  =   0.5d0 * (phiint_x_L   + phiint_x_R   - Sp * ( Bzint_x_R - Bzint_x_L ) )

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,j,1,l)  =   0.5d0 * (psiint_x_L   + psiint_x_R   - Sp * ( Exint_x_R - Exint_x_L ) )
       BzlaxEx(i,j,1,l)   = - 0.5d0 * (Bzint_y_L_m  + Bzint_y_R_m  - Sp * ( Exint_y_R - Exint_y_L ) )
       BylaxEx(i,j,1,l)   =   0.5d0 * (Byint_x_L    + Byint_x_R    - Sp * ( Exint_x_R - Exint_x_L ) )

! Flujos Ey


       BzlaxEy (i,j,1,l)  =   0.5d0 * (Bzint_x_L    + Bzint_x_R    - Sp * ( Eyint_x_R - Eyint_x_L ) )
       psilaxEy(i,j,1,l)  =   0.5d0 * (psiint_y_L   + psiint_y_R   - Sp * ( Eyint_y_R - Eyint_y_L ) )
       BxlaxEy (i,j,1,l)  = - 0.5d0 * (Bxint_x_L_m  + Bxint_x_R_m  - Sp * ( Eyint_x_R - Eyint_x_L ) )

! Flujos Ez


       BylaxEz (i,j,1,l)  = - 0.5d0 * (Byint_x_L_m  + Byint_x_R_m  - Sp * ( Ezint_x_R - Ezint_x_L ) )
       BxlaxEz (i,j,1,l)  =   0.5d0 * (Bxint_y_L    + Bxint_y_R    - Sp * ( Ezint_y_R - Ezint_y_L ) )
       psilaxEz(i,j,1,l)  =   0.5d0 * (psiint_x_L   + psiint_x_R   - Sp * ( Ezint_x_R - Ezint_x_L ) )

! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,j,1,l)     =   0.5d0 * (Jxint_L  + Jxint_R - Sp * ( qint_x_R - qint_x_L ) )
       Jylax(i,j,1,l)     =   0.5d0 * (Jyint_L  + Jyint_R - Sp * ( qint_y_R - qint_y_L ) )
       Jzlax(i,j,1,l)     =   0.5d0 * (Jzint_L  + Jzint_R - Sp * ( qint_x_R - qint_x_L ) )  ! REVISAR Jx, Jy and Jz

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,j,1,l)    =   0.5d0 * (FDxint_L  + FDxint_R - Sp * ( DDint_x_R - DDint_x_L ) )
       FDylax(i,j,1,l)    =   0.5d0 * (FDyint_L  + FDyint_R - Sp * ( DDint_y_R - DDint_y_L ) )
       FDzlax(i,j,1,l)    =   0.5d0 * (FDzint_L  + FDzint_R - Sp * ( DDint_x_R - DDint_x_L ) )

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,j,1,l)  =   0.5d0 * (Ftauxint_L + Ftauxint_R - Sp * ( tauint_x_R - tauint_x_L  ) )
       Ftauylax(i,j,1,l)  =   0.5d0 * (Ftauyint_L + Ftauyint_R - Sp * ( tauint_y_R - tauint_y_L  ) )
       Ftauzlax(i,j,1,l)  =   0.5d0 * (Ftauzint_L + Ftauzint_R - Sp * ( tauint_x_R - tauint_x_L  ) )

! -------------------------------------- EGLM-------------------------------------- ?????? REVISAR ??????

       psitauxlax(i,j,1,l)  =   0.5d0 * (psiint_x_L + psiint_x_R - Sp * ( tauint_x_R - tauint_x_L  ) )
       psitauylax(i,j,1,l)  =   0.5d0 * (psiint_y_L + psiint_y_R - Sp * ( tauint_y_R - tauint_y_L  ) )
       psitauzlax(i,j,1,l)  =   0.5d0 * (psiint_x_L + psiint_x_R - Sp * ( tauint_x_R - tauint_x_L  ) )


       phitauxlax(i,j,1,l)  =   0.5d0 * (phiint_x_L + phiint_x_R - Sp * ( tauint_x_R - tauint_x_L  ) )
       phitauylax(i,j,1,l)  =   0.5d0 * (phiint_y_L + phiint_y_R - Sp * ( tauint_y_R - tauint_y_L  ) )
       phitauzlax(i,j,1,l)  =   0.5d0 * (phiint_x_L + phiint_x_R - Sp * ( tauint_x_R - tauint_x_L  ) )

! -------------------------------------- EGLM-------------------------------------- ?????? REVISAR ??????

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,j,1,l)   =   0.5d0 * (FSxxint_L + FSxxint_R - Sp * (Sxint_x_R - Sxint_x_L ) ) 
       FSxylax(i,j,1,l)   =   0.5d0 * (FSxyint_L + FSxyint_R - Sp * (Syint_x_R - Syint_x_L ) ) 
       FSxzlax(i,j,1,l)   =   0.5d0 * (FSxzint_L + FSxzint_R - Sp * (Szint_x_R - Szint_x_L ) ) 

 
       FSyxlax(i,j,1,l)   =   0.5d0 * (FSyxint_L + FSyxint_R - Sp * (Sxint_y_R - Sxint_y_L ) ) 
       FSyylax(i,j,1,l)   =   0.5d0 * (FSyyint_L + FSyyint_R - Sp * (Syint_y_R - Syint_y_L ) )
       FSyzlax(i,j,1,l)   =   0.5d0 * (FSyzint_L + FSyzint_R - Sp * (Szint_y_R - Szint_y_L ) )


       FSzxlax(i,j,1,l)   =   0.5d0 * (FSzxint_L + FSzxint_R - Sp * (Sxint_x_R - Sxint_x_L ) )
       FSzylax(i,j,1,l)   =   0.5d0 * (FSzyint_L + FSzyint_R - Sp * (Syint_x_R - Syint_x_L ) )
       FSzzlax(i,j,1,l)   =   0.5d0 * (FSzzint_L + FSzzint_R - Sp * (Szint_x_R - Szint_x_L ) )

! -------------------------------------- EGLM-------------------------------------- ¡¡¡¡ seems ok just to revise ¡¡¡¡¡¡¡


       BxSxlax(i,j,1,l)   = 0.5d0 * (Bxint_x_L  + Bxint_x_R - Sp * (Sxint_x_R - Sxint_x_L ) )
       BySxlax(i,j,1,l)   = 0.5d0 * (Byint_y_L  + Byint_y_R - Sp * (Sxint_y_R - Sxint_y_L ) )
       BzSxlax(i,j,1,l)   = 0.5d0 * (Bzint_x_L  + Bzint_x_R - Sp * (Sxint_x_R - Sxint_x_L ) )



       BxSylax(i,j,1,l)   = 0.5d0 * (Bxint_x_L  + Bxint_x_R - Sp * (Syint_x_R - Syint_x_L ) )
       BySylax(i,j,1,l)   = 0.5d0 * (Byint_y_L  + Byint_y_R - Sp * (Syint_y_R - Syint_y_L ) )
       BzSylax(i,j,1,l)   = 0.5d0 * (Bzint_x_L  + Bzint_x_R - Sp * (Syint_x_R - Syint_x_L ) )


       BxSzlax(i,j,1,l)   = 0.5d0 * (Bxint_x_L  + Bxint_x_R - Sp * (Szint_x_R - Szint_x_L ) )
       BySzlax(i,j,1,l)   = 0.5d0 * (Byint_y_L  + Byint_y_R - Sp * (Szint_y_R - Szint_y_L ) )
       BzSzlax(i,j,1,l)   = 0.5d0 * (Bzint_x_L  + Bzint_x_R - Sp * (Szint_x_R - Szint_x_L ) )

! -------------------------------------- EGLM-------------------------------------- 

        end do !i
     end do  !j

!$OMP END DO     

!************************************************************************************
!************************************************************************************
   else

      write(*,*) "STOP: subroutine hllflux"
      write(*,*) "This dimension is not implemented yet"
      stop

   end if

!$OMP END PARALLEL    

 end subroutine hll_flow
