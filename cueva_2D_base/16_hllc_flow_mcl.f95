  !     *******************************************************************
  !     Subroutine which calculed  HLL flux flows 
  !     *******************************************************************

   subroutine hllc_flow

    use scalar
    use parameters
    use threevectors
    use fourvectors
    use funciones, only: mcl

    implicit none


   
    if (DIM == 1) then


       call boundary_conserved
       call boundary_electric
       call boundary_flux_conserved
     

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

             Vx_L = Vx (i  ,1,1) +  FAC  * mcl(Vx(i+1,1,1) - Vx(i  ,1,1) ,  &
                                               Vx(i  ,1,1) - Vx(i-1,1,1))
             Vx_R = Vx (i+1,1,1) -  FAC  * mcl(Vx(i+2,1,1) - Vx(i+1,1,1) ,  &
                                               Vx(i+1,1,1) - Vx(i  ,1,1))

             Vy_L = Vy (i  ,1,1) +  FAC  * mcl(Vy(i+1,1,1) - Vy(i  ,1,1) ,  &
                                               Vy(i  ,1,1) - Vy(i-1,1,1))
             Vy_R = Vy (i+1,1,1) -  FAC  * mcl(Vy(i+2,1,1) - Vy(i+1,1,1) ,  &
                                               Vy(i+1,1,1) - Vy(i  ,1,1))

             Vz_L = Vz (i  ,1,1) +  FAC  * mcl(Vz(i+1,1,1) - Vz(i  ,1,1) ,  &
                                               Vz(i  ,1,1) - Vz(i-1,1,1))
             Vz_R = Vz (i+1,1,1) -  FAC  * mcl(Vz(i+2,1,1) - Vz(i+1,1,1) ,  &
                                               Vz(i+1,1,1) - Vz(i  ,1,1))

             p_L  = p  (i  ,1,1)  +  FAC  * mcl(p(i+1,1,1) - p(i  ,1,1) ,   &
                                                p(i  ,1,1) - p(i-1,1,1))
             p_R  = p  (i+1,1,1)  -  FAC  * mcl(p(i+2,1,1) - p(i+1,1,1) ,   &
                                                p(i+1,1,1) - p(i  ,1,1))

             rho_L = rho(i,1,1)   +  FAC  * mcl(rho(i+1,1,1) - rho(i  ,1,1) ,  &
                                                rho(i  ,1,1) - rho(i-1,1,1))
             rho_R = rho(i+1,1,1) -  FAC  * mcl(rho(i+2,1,1) - rho(i+1,1,1) ,  &
                                                rho(i+1,1,1) - rho(i  ,1,1))

             W_L = 1.d0/sqrt(1.d0 - Vx_L**2 -  Vy_L**2 - Vz_L**2)
             W_R = 1.d0/sqrt(1.d0 - Vx_R**2 -  Vy_R**2 - Vz_R**2)


          else if (REC_PRIM == 1) then

             Vx_L = Vxint(i  ,1,1,l)   +  FAC  * mcl(Vxint(i+1,1,1,l) - Vxint(i  ,1,1,l) ,  &
                                                     Vxint(i  ,1,1,l) - Vxint(i-1,1,1,l))
             Vx_R = Vxint(i+1,1,1,l)   -  FAC  * mcl(Vxint(i+2,1,1,l) - Vxint(i+1,1,1,l) ,  &
                                                     Vxint(i+1,1,1,l) - Vxint(i  ,1,1,l))

             Vy_L = Vyint(i  ,1,1,l)   +  FAC  * mcl(Vyint(i+1,1,1,l) - Vyint(i  ,1,1,l) ,  &
                                                     Vyint(i  ,1,1,l) - Vyint(i-1,1,1,l))
             Vy_R = Vyint(i+1,1,1,l)   -  FAC  * mcl(Vyint(i+2,1,1,l) - Vyint(i+1,1,1,l) ,  &
                                                     Vyint(i+1,1,1,l) - Vyint(i  ,1,1,l))

             Vz_L = Vzint(i  ,1,1,l)   +  FAC  * mcl(Vzint(i+1,1,1,l) - Vzint(i  ,1,1,l) ,  &
                                                     Vzint(i  ,1,1,l) - Vzint(i-1,1,1,l))
             Vz_R = Vzint(i+1,1,1,l)   -  FAC  * mcl(Vzint(i+2,1,1,l) - Vzint(i+1,1,1,l) ,  &
                                                     Vzint(i+1,1,1,l) - Vzint(i  ,1,1,l))


!!$             p_L  = p  (i  ,1,1)  +  FAC  * mcl(p(i+1,1,1) - p(i  ,1,1) ,   &
!!$                                                p(i  ,1,1) - p(i-1,1,1))
!!$             p_R  = p  (i+1,1,1)  -  FAC  * mcl(p(i+2,1,1) - p(i+1,1,1) ,   &
!!$                                                p(i+1,1,1) - p(i  ,1,1))
!!$
!!$             rho_L = rho(i,1,1)   +  FAC  * mcl(rho(i+1,1,1) - rho(i  ,1,1) ,  &
!!$                                                rho(i  ,1,1) - rho(i-1,1,1))
!!$             rho_R = rho(i+1,1,1) -  FAC  * mcl(rho(i+2,1,1) - rho(i+1,1,1) ,  &
!!$                                                rho(i+1,1,1) - rho(i  ,1,1))
!!$
!!$             
             p_L  = pint (i  ,1,1,l)   +  FAC  * mcl(pint(i+1,1,1,l) - pint(i  ,1,1,l) ,  &
                                                     pint(i  ,1,1,l) - pint(i-1,1,1,l))
             p_R  = pint (i+1,1,1,l)   -  FAC  * mcl(pint(i+2,1,1,l) - pint(i+1,1,1,l) ,  &
                                                     pint(i+1,1,1,l) - pint(i  ,1,1,l))

             rho_L = rhoint(i  ,1,1,l) +  FAC  * mcl(rhoint(i+1,1,1,l) - rhoint(i  ,1,1,l) ,  &
                                                     rhoint(i  ,1,1,l) - rhoint(i-1,1,1,l))
             rho_R = rhoint(i+1,1,1,l) -  FAC  * mcl(rhoint(i+2,1,1,l) - rhoint(i+1,1,1,l) ,  &
                                                     rhoint(i+1,1,1,l) - rhoint(i  ,1,1,l))

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


!     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

!!$    Jxint_L = sigma * W_L * ( Exint_L +  vrotB_x_L  - Edotv_L * Vx_L )   &
!!$                          +   qint_L  *  Vx_L
!!$    Jxint_R = sigma * W_R * ( Exint_R +  vrotB_x_R  - Edotv_R * Vx_R )   &
!!$                          +   qint_R  *  Vx_R
!!$
!!$    Jyint_L = sigma * W_L * ( Eyint_L +  vrotB_y_L  - Edotv_L * Vy_L )   &
!!$                          +   qint_L  *  Vy_L
!!$    Jyint_R = sigma * W_R * ( Eyint_R +  vrotB_y_R  - Edotv_R * Vy_R )   &
!!$                          +   qint_R  *  Vy_R
!!$
!!$    Jzint_L = sigma * W_L * ( Ezint_L +  vrotB_z_L  - Edotv_L * Vz_L )   &
!!$                          +   qint_L  *  Vz_L
!!$    Jzint_R = sigma * W_R * ( Ezint_R +  vrotB_z_R  - Edotv_R * Vz_R )   &
!!$                          +   qint_R  *  Vz_R


!!$    Jxint_L = Jxlaxp(i  ,1,1,l) +  FAC  * mcl(Jxlaxp(i+1,1,1,l) - Jxlaxp(i  ,1,1,l) ,  &
!!$                                              Jxlaxp(i  ,1,1,l) - Jxlaxp(i-1,1,1,l)) 
!!$    Jxint_R = Jxlaxm(i+1,1,1,l) -  FAC  * mcl(Jxlaxm(i+2,1,1,l) - Jxlaxm(i+1,1,1,l) ,  &
!!$                                              Jxlaxm(i+1,1,1,l) - Jxlaxm(i  ,1,1,l)) 
!!$
!!$    Jyint_L = Jylaxp(i  ,1,1,l) +  FAC  * mcl(Jylaxp(i+1,1,1,l) - Jylaxp(i  ,1,1,l) ,  &
!!$                                              Jylaxp(i  ,1,1,l) - Jylaxp(i-1,1,1,l)) 
!!$    Jyint_R = Jylaxm(i+1,1,1,l) -  FAC  * mcl(Jylaxm(i+2,1,1,l) - Jylaxm(i+1,1,1,l) ,  &
!!$                                              Jylaxm(i+1,1,1,l) - Jylaxm(i  ,1,1,l)) 
!!$
!!$    Jzint_L = Jzlaxp(i  ,1,1,l) +  FAC  * mcl(Jzlaxp(i+1,1,1,l) - Jzlaxp(i  ,1,1,l) ,  &
!!$                                              Jzlaxp(i  ,1,1,l) - Jzlaxp(i-1,1,1,l)) 
!!$    Jzint_R = Jzlaxm(i+1,1,1,l) -  FAC  * mcl(Jzlaxm(i+2,1,1,l) - Jzlaxm(i+1,1,1,l) ,  &
!!$                                              Jzlaxm(i+1,1,1,l) - Jzlaxm(i  ,1,1,l)) 

    Jxint_L = Jxint(i  ,1,1,l)  +  FAC  * mcl(Jxint(i+1,1,1,l) - Jxint(i  ,1,1,l) ,  &
                                              Jxint(i  ,1,1,l) - Jxint(i-1,1,1,l)) 
    Jxint_R = Jxint(i+1,1,1,l)  -  FAC  * mcl(Jxint(i+2,1,1,l) - Jxint(i+1,1,1,l) ,  &
                                              Jxint(i+1,1,1,l) - Jxint(i  ,1,1,l)) 
    Jyint_L = Jyint(i  ,1,1,l)  +  FAC  * mcl(Jyint(i+1,1,1,l) - Jyint(i  ,1,1,l) ,  &
                                              Jyint(i  ,1,1,l) - Jyint(i-1,1,1,l)) 
    Jyint_R = Jyint(i+1,1,1,l)  -  FAC  * mcl(Jyint(i+2,1,1,l) - Jyint(i+1,1,1,l) ,  &
                                              Jyint(i+1,1,1,l) - Jyint(i  ,1,1,l)) 
    Jzint_L = Jzint(i  ,1,1,l)  +  FAC  * mcl(Jzint(i+1,1,1,l) - Jzint(i  ,1,1,l) ,  &
                                              Jzint(i  ,1,1,l) - Jzint(i-1,1,1,l))
    Jzint_R = Jzint(i+1,1,1,l)  -  FAC  * mcl(Jzint(i+2,1,1,l) - Jzint(i+1,1,1,l) ,  &
                                              Jzint(i+1,1,1,l) - Jzint(i  ,1,1,l)) 

!!$            Jxint_R   = Jxint_L
!!$              
!!$            Jyint_R   = Jyint_L
!!$             
!!$            Jzint_R   = Jzint_L


    !     Density flux (FD = \rho W V)

    FDxint_L = DDint_L * Vx_L
    FDxint_R = DDint_R * Vx_R

    FDyint_L = DDint_L * Vy_L
    FDyint_R = DDint_R * Vy_R

    FDzint_L = DDint_L * Vz_L
    FDzint_R = DDint_R * Vz_R

              !     Entalphy 

    epsilon_L = p_L / ((gamma-1.d0) * rho_L)
    epsilon_R = p_R / ((gamma-1.d0) * rho_R)

    enthpy_L   =  rho_L * ( 1.d0 + epsilon_L ) + p_L 
    enthpy_R   =  rho_R * ( 1.d0 + epsilon_R ) + p_R

!    print*, enthpy_R,enthpy_L,rho_L



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



!/////////////////////// SOURCE TERMS ////////////////////////////////////


    Omega1_L =  qint_L - kappa *  psiint_L
    Omega1_R =  qint_R - kappa *  psiint_R
    DOmega_1 =  (Omega1_R - Omega1_L)

    Omega2_L = - kappa *  phiint_L
    Omega2_R = - kappa *  phiint_R
    DOmega_2 =  (Omega2_R - Omega2_L)

    Omega3_L = - Jxint_L
    Omega3_R = - Jxint_R
    DOmega_3 =  (Omega3_R - Omega3_L)

    Omega4_L = - Jyint_L
    Omega4_R = - Jyint_R
    DOmega_4 =  (Omega4_R - Omega4_L)

    Omega5_L = - Jzint_L
    Omega5_R = - Jzint_R
    DOmega_5 =  (Omega5_R - Omega5_L)


!/////////////////////// HLL Variables  ////////////////////////////////////

! we difine de characteristic velocitis in diferent form from notes as -----> s2 = lambda_r  and s1 = lambda_l

    !     hll  values of augmented fields

    psi_hll = ( s2 * psiint_R - s1 * psiint_L + Exint_L  - Exint_R   ) / (s2 -s1)
    phi_hll = ( s2 * phiint_R - s1 * phiint_L + Bxint_L  - Bxint_R   ) / (s2 -s1)

    !     hll  values of electric field

    Ex_hll  = ( s2 * Exint_R  - s1 * Exint_L  + psiint_L  - psiint_R  ) / (s2 -s1)
    Ey_hll  = ( s2 * Eyint_R  - s1 * Eyint_L  + Bzint_L   - Bzint_R   ) / (s2 -s1)
    Ez_hll  = ( s2 * Ezint_R  - s1 * Ezint_L  + Byint_L_m - Byint_R_m   ) / (s2 -s1)  !ojo minus sign ¡¡¡


    !     hll  values of magnetic field

    Bx_hll  = ( s2 * Bxint_R  - s1 * Bxint_L  + phiint_L  - phiint_R  ) / (s2 -s1)
    By_hll  = ( s2 * Byint_R  - s1 * Byint_L  + Ezint_L_m - Ezint_R_m   ) / (s2 -s1)  !ojo minus sign ¡¡¡
    Bz_hll  = ( s2 * Bzint_R  - s1 * Bzint_L  + Eyint_L   - Eyint_R   ) / (s2 -s1)


    !     hll  value of charge density 

    q_hll   = ( s2 * qint_R   - s1 * qint_L   + Jxint_L   - Jxint_R   ) / (s2 -s1)

    !     hll  value of mass conservation 

    DD_hll  = ( s2 * DDint_R  - s1 * DDint_L  + FDxint_L  - FDxint_R  ) / (s2 -s1)

    !     hll  value of energy conservation 

    tau_hll = ( s2 * tauint_R - s1 * tauint_L + Ftauxint_L - Ftauxint_R ) / (s2 -s1)

    !    hll  value of momentum S

    Sx_hll  = ( s2 * Sxint_R  - s1 * Sxint_L  + FSxxint_L  - FSxxint_R  ) / (s2 -s1)
    Sy_hll  = ( s2 * Syint_R  - s1 * Syint_L  + FSxyint_L  - FSxyint_R  ) / (s2 -s1)
    Sz_hll  = ( s2 * Szint_R  - s1 * Szint_L  + FSxzint_L  - FSxzint_R  ) / (s2 -s1)


!/////////////////////// HLL Fluxes  ////////////////////////////////////

! Electric Gauss HLL flux
!_________________________________________________________________________

      psi_hll_flux  =   0.5d0 * (Exint_L + Exint_R - ( psiint_R - psiint_L ) )  

! Magnetic Gauss HLL flux
!_________________________________________________________________________


      phi_hll_flux  =   0.5d0 * (Bxint_L + Bxint_R - ( phiint_R - phiint_L ) )

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       Bx_hll_flux  =   0.5d0 * (phiint_L   + phiint_R   - ( Bxint_R - Bxint_L ) )

! Flujos By

       By_hll_flux  =   0.5d0 * (Ezint_L_m   + Ezint_R_m    - ( Byint_R - Byint_L ) ) ! minus sign

! Flujos Bz

       Bz_hll_flux  =   0.5d0 * (Eyint_L    + Eyint_R    - ( Bzint_R - Bzint_L ) )

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       Ex_hll_flux  =   0.5d0 * (psiint_L   + psiint_R   - ( Exint_R - Exint_L ) )

! Flujos Ey


       Ey_hll_flux  =   0.5d0 * (Bzint_L    + Bzint_R    - ( Eyint_R - Eyint_L ) )

! Flujos Ez


       Ez_hll_flux  =   0.5d0 * (Byint_L_m  + Byint_R_m    - ( Ezint_R - Ezint_L ) ) ! minus sign

! Conserved current HLL flux
!_________________________________________________________________________


       q_hll_flux   =   0.5d0 * (Jxint_L  + Jxint_R - ( qint_R - qint_L ) )

! Conserved Mass HLL flux
!_________________________________________________________________________


       DD_hll_flux  =   0.5d0 * (FDxint_L  + FDxint_R - ( DDint_R - DDint_L ) )

! Conserved Energy HLL flux
!_________________________________________________________________________


       tau_hll_flux =   0.5d0 * (Ftauxint_L + Ftauxint_R - ( tauint_R - tauint_L  ) )


! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       Sx_hll_flux  =   0.5d0 * (FSxxint_L + FSxxint_R - (Sxint_R - Sxint_L ) ) 
       Sy_hll_flux  =   0.5d0 * (FSxyint_L + FSxyint_R - (Syint_R - Syint_L ) ) 
       Sz_hll_flux  =   0.5d0 * (FSxzint_L + FSxzint_R - (Szint_R - Szint_L ) )



! Dedudction of HLLC solver M.A Aloy 08-2017

! 1) First calculate the continuos conserved variables in start region (q, E, B) as its avarage HLL values

   Exstr_R  = Ex_hll  
   Exstr_L  = Ex_hll
   
   Eystr_R  = Ey_hll
   Eystr_L  = Ey_hll 

   Ezstr_R  = Ez_hll 
   Ezstr_L  = Ez_hll 

   Bxstr_R  = Bx_hll
   Bxstr_L  = Bx_hll 
   
   Bystr_R  = By_hll 
   Bystr_L  = By_hll 

   Bzstr_R  = Bz_hll 
   Bzstr_L  = Bz_hll 

   qstr_R   = q_hll  
   qstr_L   = q_hll



   !  2) then find the contac wave velocitu vx^* by means of quadratic equation
   !                    a (vx^*)^2 + b vx^* + c = 0
   !  with:
   !           a = tau_hll_flux - ErotB_str_x
   !           b = E2_perp + B2_perp - tau_hll - Sx_hll_flux
   !           c = Sx_hll - ErotB_str_x
       
   E2_perp     = (Eystr_R**2 + Ezstr_R**2)
   B2_perp     = (Bystr_R**2 + Bzstr_R**2)

   E2B2_str    = 0.5d0 * ( (Exstr_R**2 + Eystr_R**2 + Ezstr_R**2) + (Bxstr_R**2 + Bystr_R**2 + Bzstr_R**2) )

   ErotB_str_x = (Bzstr_R * Eystr_R - Bystr_R * Ezstr_R)
   ErotB_str_y = (Bxstr_R * Ezstr_R - Bzstr_R * Exstr_R)
   ErotB_str_z = (Bystr_R * Exstr_R - Bxstr_R * Eystr_R)

   

   a_hll       = tau_hll_flux - ErotB_str_x
   b_hll       = E2_perp + B2_perp - tau_hll - Sx_hll_flux
   c_hll       = Sx_hll - ErotB_str_x

   ! Analitic solution quadratic equations

   
   eps  = 1.d-10
   eps1 = 1.d-6
   eps2 = 1.d-3
   

   det_hll     = b_hll**2 - 4.d0 * a_hll * c_hll

   q_root1 = - 0.5d0 * ( b_hll + sign(1.d0,b_hll) * abs(sqrt( det_hll )) )

   rt1    = q_root1 / a_hll
   rt2    = c_hll  / q_root1 ! this is the negative root with has physical meaning

   if (det_hll .lt. 0.d0) then

      print*, "discriminante imaginario en la ecuacion cuadratica para obtener Vxstr en HLLC Riemman solver"
      print*, "det_hll =", det_hll
      rt2    = c_hll  / b_hll
      
   end if


!!$   if ( abs(rt2) .le. 1.d0 ) then
!!$
!!$      if ( abs(rt2) .le. 1.d-6 ) rt2 = 0.d0
!!$
!!$      Vxstr = rt2
!!$
!!$   else if ( abs(rt2) .gt. 1.d0 ) then
!!$
!!$      print*, "superluminal velocity ( -1 > Vxstr > 1) in  quadratic equation HLLC Solver"
!!$      print*, "rt1 =", rt1
!!$      print*, "rt2 =", rt2
   
!!$      Vxstr= max(-1.d0 + eps1 , min( 1.d0 - eps1, Vxstr ))
   
!!$   end if


   

!!$Para que \lambda^* sea igual a cero en la ecuación cuadrática (42), debe suceder que c=0.
!!$
!!$Sin embargo, incluso en el caso c=0, la ecuación cuadrática (42) tiene dos soluciones:
!!$
!!$\lambda^*_a =0
!!$\lambda^*_b = -b/a (suponiendo que a \ne 0)
!!$
!!$Nótese que las soluciones generales de la cuadrática son:
!!$
!!$\lambda^*_1 = [-b - \sqrt( b^2 - 4ac ) ] / (2a)
!!$\lambda^*_2 = [-b + \sqrt( b^2 - 4ac ) ] / (2a)
!!$
!!$De esas dos soluciones, se toma sólo \lambda^*_1, pues como razonan Mignone & Bodo (2005, 2006), esa es la única que garantiza que
!!$
!!$\lambda_l \le \lambda^* \le \lambda_r.  [A]
!!$
!!$Fijémonos ahora que, cuando c ~> 0, la única solución que se anula de las dos en el caso general es la segunda. Es decir, que en caso c~>0, tenemos:
!!$
!!$
!!$\lambda^*_1 -> \lambda^*_b = -b/a
!!$\lambda^*_2 -> \lambda^*_a = 0
!!$
!!$Como hemos dicho que la solución \lambda^*_2, no garantiza que se cumpla la condición [A], concluimos que: [eps es un valor positivo pequeño, p.e., eps = 1e-13]

!!$   eps = 1d-13   
!!$
!!$
!!$   if ( abs(c_hll) .le. eps ) then
!!$
!!$      if ( abs(a_hll) .ge. eps ) then
!!$
!!$        Vxstr  = - b_hll / a_hll 
!!$
!!$     else
!!$!*
!!$!*       Este caso sólo se da cuando a~0 y c~0
!!$!*
!!$          Vxstr  = 0.d0
!!$
!!$       end if
!!$
!!$    else

!!$       [solución numérica de la ecuación cuadrática]
!!$
!!$   This function obtains the root of the quadratic equation
!!$
!!$       a*v^2 - b*v + c = 0
!!$
!!$   in the physical range [-1,1] usign a Newton-Rapson iterative method.


   Vxstr       = get_hllc_vel(a_hll,-b_hll,c_hll,rt2)


   if ( abs(Vxstr) .le. 1.d-9 ) then

      Vxstr = 0.d0

   else if ( abs(Vxstr) .gt. 1.d0 ) then

      print*, "superluminal velocity ( -1 > Vxstr > 1) in  quadratic equation HLLC Solver"
      print*, "Vxstr =", Vxstr

      Vxstr= max(-1.d0 + eps1 , min( 1.d0 - eps1, Vxstr ))
      
   end if

!!$end if



      ! 3)  Find the total pressure (P = p_{gas} + p_{em}) in start region


   pstr = Sx_hll_flux + (Exstr_R**2 + Bxstr_R**2 ) - ( tau_hll_flux - ErotB_str_x ) * Vxstr

      ! 4) Calculate the remainder conserved variables in its right and left star values

      
   psistr_R =  psi_hll
   psistr_L =  psi_hll 

   phistr_R =  phi_hll
   phistr_L =  phi_hll

!!$!=======================================================================================================================   
!!$! characteristic sound velocity, with cs= sqrt(gamma * p/(rho * h), observe that you use h = rho * h
!!$!=======================================================================================================================      
!!$
!!$   Cs_R      = sqrt(gamma * p_R / (enthpy_R))
!!$   Cs_L      = sqrt(gamma * p_L / (enthpy_L))
!!$
!!$   V2_R      =  Vx_R**2 +  Vy_R**2 + Vz_R**2
!!$   V2_L      =  Vx_L**2 +  Vy_L**2 + Vz_L**2
!!$
!!$   lambda_H_p_R = ( Vx_R * ( 1.d0 - Cs_R**2) + Cs_R * sqrt( (1.d0 - V2_R) * (1.d0 - V2_R * Cs_R**2 - Vx_R * (1.d0 - Cs_R**2) )) ) &
!!$                / ( 1.d0 - Cs_R**2 * V2_R)
!!$
!!$   lambda_H_p_L = ( Vx_L * ( 1.d0 - Cs_L**2) + Cs_L * sqrt( (1.d0 - V2_L) * (1.d0 - V2_L * Cs_L**2 - Vx_L * (1.d0 - Cs_L**2) )) ) &
!!$                / ( 1.d0 - Cs_L**2 * V2_L)
!!$
!!$   
!!$   lambda_H_m_R = ( Vx_R * ( 1.d0 - Cs_R**2) - Cs_R * sqrt( (1.d0 - V2_R) * (1.d0 - V2_R * Cs_R**2 - Vx_R * (1.d0 - Cs_R**2) )) ) &
!!$                / ( 1.d0 - Cs_R**2 * V2_R)
!!$
!!$   lambda_H_m_L = ( Vx_L * ( 1.d0 - Cs_L**2) - Cs_L * sqrt( (1.d0 - V2_L) * (1.d0 - V2_L * Cs_L**2 - Vx_L * (1.d0 - Cs_L**2) )) ) &
!!$                / ( 1.d0 - Cs_L**2 * V2_L)
!!$
!!$
!!$   s1 = min(lambda_H_m_L, lambda_H_m_R, 0.d0)
!!$   s2 = max(lambda_H_p_L, lambda_H_p_R, 0.d0)
!!$
!!$!   print*, "s1", s1
!!$!   print*, "s2", s2

!=======================================================================================================================         
!=======================================================================================================================      
   

!!$   Sxstr_R  = (  s2 * Sx_hll - Sx_hll_flux + pstr - Bxstr_R**2 - Exstr_R**2 - ErotB_str_x * Vxstr ) &
!!$              / (s2 - Vxstr) 
!!$   Sxstr_L  = ( -s1 * Sx_hll + Sx_hll_flux - pstr + Bxstr_R**2 + Exstr_R**2 + ErotB_str_x * Vxstr ) &
!!$              / (Vxstr - s1) 
!!$
   Sxstr_R  = (  s2 * Sx_hll - tau_hll_flux * Vxstr ) / (s2 - Vxstr) 
   Sxstr_L  = ( -s1 * Sx_hll + tau_hll_flux * Vxstr ) / (Vxstr - s1) 
   
   Systr_R  = (  s2 * Sy_hll - Sy_hll_flux - ErotB_str_y * Vxstr - Exstr_R * Eystr_R - Bxstr_R * Bystr_R )               &
              / (s2 - Vxstr)
   Systr_L  = ( -s1 * Sy_hll + Sy_hll_flux + ErotB_str_y * Vxstr + Exstr_L * Eystr_L + Bxstr_L * Bystr_L )               &
              / (Vxstr - s1)

   Szstr_R  = (  s2 * Sz_hll - Sz_hll_flux - ErotB_str_z * Vxstr - Exstr_R * Ezstr_R - Bxstr_R * Bzstr_R )               &
              / (s2 - Vxstr)
   Szstr_L  = ( -s1 * Sz_hll + Sz_hll_flux + ErotB_str_z * Vxstr + Exstr_L * Ezstr_L + Bxstr_L * Bzstr_L )               &
              / (Vxstr - s1)


   taustr_R = (Sxstr_R - ErotB_str_x) / Vxstr - pstr + E2B2_str 
   taustr_L = (Sxstr_L - ErotB_str_x) / Vxstr - pstr + E2B2_str
!!$
!!$   taustr_R = (  s2 * tau_hll - tau_hll_flux + (pstr - E2B2_str) * Vxstr + ErotB_str_x ) &
!!$              / (s2 - Vxstr)
!!$   taustr_L = ( -s1 * tau_hll + tau_hll_flux - (pstr - E2B2_str) * Vxstr - ErotB_str_x ) &
!!$              / (Vxstr - s1)

   DDstr_R  = (  s2 * DD_hll - DD_hll_flux ) / (s2 - Vxstr)
   DDstr_L  = ( -s1 * DD_hll + DD_hll_flux ) / (Vxstr - s1)


   ! 5) Finally find the numerical fluxes fot the cases Vxstr=0 and Vxstr =\ 0

   !**************************** STAR FLUXES ******************************************


   ! Using Rankine-Hugoniot conditions and consistency flux conditions, we can calculate
   ! Right and Left star flux values as:



       if ( Vxstr == 0.d0 ) then

          
    ! EM FLUXES

    psistr_flux_R =  psi_hll_flux 
    psistr_flux_L =  psi_hll_flux 
    
    phistr_flux_R =  phi_hll_flux 
    phistr_flux_L =  phi_hll_flux 

    Exstr_flux_R  =  Ex_hll_flux  
    Exstr_flux_L  =  Ex_hll_flux  

    Eystr_flux_R  =  Ey_hll_flux  
    Eystr_flux_L  =  Ey_hll_flux  

    Ezstr_flux_R  =  Ez_hll_flux  
    Ezstr_flux_L  =  Ez_hll_flux  

    Bxstr_flux_R  =  Bx_hll_flux  
    Bxstr_flux_L  =  Bx_hll_flux  

    Bystr_flux_R  =  By_hll_flux  
    Bystr_flux_L  =  By_hll_flux  

    Bzstr_flux_R  =  Bz_hll_flux  
    Bzstr_flux_L  =  Bz_hll_flux  
    

    !     Density flux (FD = \rho W V)

    FDxstr_R = DD_hll_flux  
    FDxstr_L = DD_hll_flux  

    ! these values must be canceled when we made diference of fluxes --->

    FDystr_R = 1.d0
    FDystr_L = 1.d0
    FDzstr_R = 1.d0
    FDzstr_L = 1.d0


    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxstr_L  = q_hll_flux
    Jxstr_R  = q_hll_flux  

    ! these values must be canceled when we made diference of fluxes --->

    Jystr_R = 1.d0
    Jystr_L = 1.d0
    Jzstr_R = 1.d0
    Jzstr_L = 1.d0

    !     Flujos conservados de energía

    Ftauxstr_R = tau_hll_flux  
    Ftauxstr_L = tau_hll_flux  

! these values must be canceled when we made diference of fluxes --->

    Ftauystr_R = 1.d0
    Ftauystr_L = 1.d0
    Ftauzstr_R = 1.d0
    Ftauzstr_L = 1.d0


    !     Components of the flux momentum tensor 

    FSxxstr_R =  Sx_hll_flux  
    FSxxstr_L =  Sx_hll_flux  

    FSxystr_R =  Sy_hll_flux  
    FSxystr_L =  Sy_hll_flux  

    FSxzstr_R =  Sz_hll_flux  
    FSxzstr_L =  Sz_hll_flux  

! these values must be canceled when we made diference of fluxes --->

    FSyxstr_R = 1.d0
    FSyxstr_L = 1.d0
    FSyystr_R = 1.d0
    FSyystr_L = 1.d0
    FSyzstr_R = 1.d0
    FSyzstr_L = 1.d0

    FSzxstr_R = 1.d0
    FSzxstr_L = 1.d0
    FSzystr_R = 1.d0
    FSzystr_L = 1.d0
    FSzzstr_R = 1.d0
    FSzzstr_L = 1.d0


       else

          ! We use the RH conditions
          !  F^*_l = F_l + s1 (U^*_l - U_l)
          !  F^*_r = F_r + s1 (U^*_r - U_r)

!   EM FLUXES       

    psistr_flux_R =  psi_hll_flux  
    psistr_flux_L =  psi_hll_flux  
    
    phistr_flux_R =  phi_hll_flux  
    phistr_flux_L =  phi_hll_flux  

    Exstr_flux_R  =  Ex_hll_flux  
    Exstr_flux_L  =  Ex_hll_flux  

    Eystr_flux_R  =  Ey_hll_flux  
    Eystr_flux_L  =  Ey_hll_flux  

    Ezstr_flux_R  =  Ez_hll_flux  
    Ezstr_flux_L  =  Ez_hll_flux  

    Bxstr_flux_R  =  Bx_hll_flux  
    Bxstr_flux_L  =  Bx_hll_flux  

    Bystr_flux_R  =  By_hll_flux  
    Bystr_flux_L  =  By_hll_flux  

    Bzstr_flux_R  =  Bz_hll_flux  
    Bzstr_flux_L  =  Bz_hll_flux  

!!$    Bxstr_flux_R  =  phiint_R + s2 * ( Bxstr_R  -  Bxint_R) 
!!$    Bxstr_flux_L  =  phiint_L + s1 * ( Bxstr_L  -  Bxint_L) 
!!$
!!$    Bystr_flux_R  =  Ezint_R  + s2 * ( Bystr_R  -  Byint_R)
!!$    Bystr_flux_L  =  Ezint_L  + s1 * ( Bystr_L  -  Byint_L)
!!$
!!$    Bzstr_flux_R  =  Eyint_R  + s2 * ( Bzstr_R  -  Bzint_R)
!!$    Bzstr_flux_L  =  Eyint_L  + s1 * ( Bzstr_L  -  Bzint_L)
!!$    
!!$    Exstr_flux_R  =  psiint_R + s2 * ( Exstr_R  -  Exint_R)
!!$    Exstr_flux_L  =  psiint_L + s1 * ( Exstr_L  -  Exint_L)
!!$
!!$    Eystr_flux_R  =  Bzint_R  + s2 * ( Eystr_R  -  Eyint_R)
!!$    Eystr_flux_L  =  Bzint_L  + s1 * ( Eystr_L  -  Eyint_L)
!!$
!!$    Ezstr_flux_R  =  Byint_R  + s2 * ( Ezstr_R  -  Ezint_R)
!!$    Ezstr_flux_L  =  Byint_L  + s1 * ( Ezstr_L  -  Ezint_L)

!   Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxstr_L       = q_hll_flux ! Jxint_R   + s2 * ( qstr_R   -   qint_R)
    Jxstr_R       = q_hll_flux ! Jxint_L   + s1 * ( qstr_L   -   qint_L)

!   these values must be canceled when we made diference of fluxes --->

    Jystr_R       = 1.d0
    Jystr_L       = 1.d0
    Jzstr_R       = 1.d0
    Jzstr_L       = 1.d0
    
!   Density flux (FD = \rho W V)


!!$    FDxstr_R      = FDxint_R  + s2 * ( DDstr_R  -  DDint_R) ! DDstr_R * Vxstr !
!!$    FDxstr_L      = FDxint_L  + s1 * ( DDstr_L  -  DDint_L) ! DDstr_L * Vxstr !


    ! PDF NOTES (pag 80) VERSION
    
    FDxstr_L = DD_hll_flux  - s1 * (s2 - Vxstr) * (DDstr_R - DDstr_L) / (s2 -s1)
    FDxstr_R = DD_hll_flux  + s2 * (Vxstr - s1) * (DDstr_R - DDstr_L) / (s2 -s1)


!   these values must be canceled when we made diference of fluxes --->

    FDystr_R      = 1.d0
    FDystr_L      = 1.d0
    FDzstr_R      = 1.d0
    FDzstr_L      = 1.d0

!   Flujos conservados de energía

!!$    Ftauxstr_R    = Sxstr_R !Ftauxint_R + s2 * (taustr_R - tauint_R) ! 
!!$    Ftauxstr_L    = Sxstr_L !Ftauxint_L + s1 * (taustr_L - tauint_L) ! 

    ! PDF NOTES (pag 80) VERSION

    Ftauxstr_L = tau_hll_flux  - s1 * (s2 - Vxstr) * (taustr_R - taustr_L) / (s2 -s1)
    Ftauxstr_R = tau_hll_flux  + s2 * (Vxstr - s1) * (taustr_R - taustr_L) / (s2 -s1)


!   these values must be canceled when we made diference of fluxes --->

    Ftauystr_R    = 1.d0
    Ftauystr_L    = 1.d0
    Ftauzstr_R    = 1.d0
    Ftauzstr_L    = 1.d0



!   Components of the flux momentum tensor 

!!$    FSxxstr_R     = - Exstr_R**2 - Bxstr_R**2 + ( Sxstr_R - ErotB_str_x ) * Vxstr + pstr !FSxxint_R + s2 * ( Sxstr_R - Sxint_R ) ! 
!!$    FSxxstr_L     = - Exstr_L**2 - Bxstr_L**2 + ( Sxstr_L - ErotB_str_x ) * Vxstr + pstr !FSxxint_L + s1 * ( Sxstr_L - Sxint_L ) ! 
!!$
!!$    FSxystr_R     = - Exstr_R * Eystr_R - Bxstr_R * Bystr_R + (Systr_R - ErotB_str_y) * Vxstr !FSxyint_R + s2 * ( Systr_R - Syint_R ) ! 
!!$    FSxystr_L     = - Exstr_L * Eystr_L - Bxstr_L * Bystr_L + (Systr_L - ErotB_str_y) * Vxstr !FSxyint_L + s1 * ( Systr_L - Syint_L ) ! 
!!$
!!$    FSxzstr_R     = - Exstr_R * Ezstr_R - Bxstr_R * Bzstr_R + (Szstr_R - ErotB_str_z) * Vxstr !FSxzint_R + s2 * ( Szstr_R - Szint_R ) ! 
!!$    FSxzstr_L     = - Exstr_L * Ezstr_L - Bxstr_L * Bzstr_L + (Szstr_L - ErotB_str_z) * Vxstr !FSxzint_L + s1 * ( Szstr_L - Szint_L ) ! 


!     PDF NOTES (pag 80) VERSION

    FSxxstr_L =  Sx_hll_flux  - s1 * (s2 - Vxstr) * (Sxstr_R - Sxstr_L) / (s2 -s1)
    FSxxstr_R =  Sx_hll_flux  + s2 * (Vxstr - s1) * (Sxstr_R - Sxstr_L) / (s2 -s1)

    FSxystr_L =  Sy_hll_flux  - s1 * (s2 - Vxstr) * (Systr_R - Systr_L) / (s2 -s1)
    FSxystr_R =  Sy_hll_flux  + s2 * (Vxstr - s1) * (Systr_R - Systr_L) / (s2 -s1)

    FSxzstr_L =  Sz_hll_flux  - s1 * (s2 - Vxstr) * (Szstr_R - Szstr_L) / (s2 -s1)
    FSxzstr_R =  Sz_hll_flux  + s2 * (Vxstr - s1) * (Szstr_R - Szstr_L) / (s2 -s1)

!   these values must be canceled when we made diference of fluxes --->

    FSyxstr_R = 1.d0
    FSyxstr_L = 1.d0
    FSyystr_R = 1.d0
    FSyystr_L = 1.d0
    FSyzstr_R = 1.d0
    FSyzstr_L = 1.d0

    FSzxstr_R = 1.d0
    FSzxstr_L = 1.d0
    FSzystr_R = 1.d0
    FSzystr_L = 1.d0
    FSzzstr_R = 1.d0
    FSzzstr_L = 1.d0    


    end if


!******************************************************************************************************************

    ! Finally we use the HLLC criterion for flux


    if ( s1 .ge. 0.d0 ) then

! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,1,1,l)  =   Exint_L
       Eylaxpsi(i,1,1,l)  =   1.d0
       Ezlaxpsi(i,1,1,l)  =   1.d0

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,1,1,l)  =   Bxint_L
       Bylaxphi(i,1,1,l)  =   1.d0
       Bzlaxphi(i,1,1,l)  =   1.d0

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,1,1,l)  =   phiint_L
       EzlaxBx (i,1,1,l)  =   1.d0
       EylaxBx (i,1,1,l)  =   1.d0

! Flujos By

       EzlaxBy (i,1,1,l)  =  -Ezint_L ! Ezint_L_m  minus sign
       philaxBy(i,1,1,l)  =   1.d0
       ExlaxBy (i,1,1,l)  =   1.d0

! Flujos Bz

       EylaxBz (i,1,1,l)  =   Eyint_L
       ExlaxBz (i,1,1,l)  =   1.d0
       philaxBz(i,1,1,l)  =   1.d0

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,1,1,l)  =   psiint_L
       BzlaxEx(i,1,1,l)   =   1.d0
       BylaxEx(i,1,1,l)   =   1.d0

! Flujos Ey


       BzlaxEy (i,1,1,l)  =   Bzint_L
       psilaxEy(i,1,1,l)  =   1.d0
       BxlaxEy (i,1,1,l)  =   1.d0

! Flujos Ez


       BylaxEz (i,1,1,l)  =  - Byint_L ! Byint_L_m  minus sign
       BxlaxEz (i,1,1,l)  =   1.d0
       psilaxEz(i,1,1,l)  =   1.d0


! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,1,1,l)     =   Jxint_L  
       Jylax(i,1,1,l)     =   1.d0
       Jzlax(i,1,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,1,1,l)    =   FDxint_L  
       FDylax(i,1,1,l)    =   1.d0
       FDzlax(i,1,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,1,1,l)  =   Ftauxint_L 
       Ftauylax(i,1,1,l)  =   1.d0
       Ftauzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,1,1,l)  =   psiint_L 
       psitauylax(i,1,1,l)  =   1.d0
       psitauzlax(i,1,1,l)  =   1.d0


       phitauxlax(i,1,1,l)  =   phiint_L 
       phitauylax(i,1,1,l)  =   1.d0
       phitauzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,1,1,l)   =   FSxxint_L 
       FSxylax(i,1,1,l)   =   FSxyint_L 
       FSxzlax(i,1,1,l)   =   FSxzint_L 

 
       FSyxlax(i,1,1,l)   =   1.d0
       FSyylax(i,1,1,l)   =   1.d0
       FSyzlax(i,1,1,l)   =   1.d0


       FSzxlax(i,1,1,l)  =   1.d0
       FSzylax(i,1,1,l)  =   1.d0
       FSzzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,1,1,l)   = Bxint_L  
       BySxlax(i,1,1,l)   = 1.d0
       BzSxlax(i,1,1,l)   = 1.d0



       BxSylax(i,1,1,l)   = Bxint_L  
       BySylax(i,1,1,l)   = 1.d0
       BzSylax(i,1,1,l)   = 1.d0


       BxSzlax(i,1,1,l)   = Bxint_L  
       BySzlax(i,1,1,l)   = 1.d0
       BzSzlax(i,1,1,l)   = 1.d0

! -------------------------------------- EGLM-------------------------------------- 


    else if ( s1 .lt. 0.d0 .and. 0.d0 .le. Vxstr ) then

! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,1,1,l)  =   psistr_flux_L
       Eylaxpsi(i,1,1,l)  =   1.d0
       Ezlaxpsi(i,1,1,l)  =   1.d0

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,1,1,l)  =   phistr_flux_L
       Bylaxphi(i,1,1,l)  =   1.d0
       Bzlaxphi(i,1,1,l)  =   1.d0

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,1,1,l)  =   Bxstr_flux_L
       EzlaxBx (i,1,1,l)  =   1.d0
       EylaxBx (i,1,1,l)  =   1.d0

! Flujos By

       EzlaxBy (i,1,1,l)  =  - Bystr_flux_L ! minus sign
       philaxBy(i,1,1,l)  =   1.d0
       ExlaxBy (i,1,1,l)  =   1.d0

! Flujos Bz

       EylaxBz (i,1,1,l)  =   Bzstr_flux_L
       ExlaxBz (i,1,1,l)  =   1.d0
       philaxBz(i,1,1,l)  =   1.d0

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,1,1,l)  =   Exstr_flux_L
       BzlaxEx(i,1,1,l)   =   1.d0
       BylaxEx(i,1,1,l)   =   1.d0

! Flujos Ey


       BzlaxEy (i,1,1,l)  =   Eystr_flux_L
       psilaxEy(i,1,1,l)  =   1.d0
       BxlaxEy (i,1,1,l)  =   1.d0

! Flujos Ez


       BylaxEz (i,1,1,l)  =  - Ezstr_flux_L ! minus sign
       BxlaxEz (i,1,1,l)  =   1.d0
       psilaxEz(i,1,1,l)  =   1.d0


! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,1,1,l)     =   Jxstr_L  
       Jylax(i,1,1,l)     =   1.d0
       Jzlax(i,1,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,1,1,l)    =   FDxstr_L  
       FDylax(i,1,1,l)    =   1.d0
       FDzlax(i,1,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,1,1,l)  =   Ftauxstr_L 
       Ftauylax(i,1,1,l)  =   1.d0
       Ftauzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,1,1,l)  =   psistr_L 
       psitauylax(i,1,1,l)  =   1.d0
       psitauzlax(i,1,1,l)  =   1.d0


       phitauxlax(i,1,1,l)  =   phistr_L 
       phitauylax(i,1,1,l)  =   1.d0
       phitauzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,1,1,l)   =   FSxxstr_L 
       FSxylax(i,1,1,l)   =   FSxystr_L 
       FSxzlax(i,1,1,l)   =   FSxzstr_L 

 
       FSyxlax(i,1,1,l)   =   1.d0
       FSyylax(i,1,1,l)   =   1.d0
       FSyzlax(i,1,1,l)   =   1.d0


       FSzxlax(i,1,1,l)  =   1.d0
       FSzylax(i,1,1,l)  =   1.d0
       FSzzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,1,1,l)   = Bxstr_L  
       BySxlax(i,1,1,l)   = 1.d0
       BzSxlax(i,1,1,l)   = 1.d0



       BxSylax(i,1,1,l)   = Bxstr_L  
       BySylax(i,1,1,l)   = 1.d0
       BzSylax(i,1,1,l)   = 1.d0


       BxSzlax(i,1,1,l)   = Bxstr_L  
       BySzlax(i,1,1,l)   = 1.d0
       BzSzlax(i,1,1,l)   = 1.d0

! -------------------------------------- EGLM-------------------------------------- 

    else if (Vxstr .lt. 0.d0 .and. 0.d0 .le. s2) then


! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,1,1,l)  =   psistr_flux_R
       Eylaxpsi(i,1,1,l)  =   1.d0
       Ezlaxpsi(i,1,1,l)  =   1.d0

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,1,1,l)  =   phistr_flux_R 
       Bylaxphi(i,1,1,l)  =   1.d0
       Bzlaxphi(i,1,1,l)  =   1.d0

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,1,1,l)  =   Bxstr_flux_R
       EzlaxBx (i,1,1,l)  =   1.d0
       EylaxBx (i,1,1,l)  =   1.d0

! Flujos By

       EzlaxBy (i,1,1,l)  =  - Bystr_flux_R ! minus sign
       philaxBy(i,1,1,l)  =   1.d0
       ExlaxBy (i,1,1,l)  =   1.d0

! Flujos Bz

       EylaxBz (i,1,1,l)  =   Bzstr_flux_R
       ExlaxBz (i,1,1,l)  =   1.d0
       philaxBz(i,1,1,l)  =   1.d0

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,1,1,l)  =   Exstr_flux_R
       BzlaxEx(i,1,1,l)   =   1.d0
       BylaxEx(i,1,1,l)   =   1.d0

! Flujos Ey


       BzlaxEy (i,1,1,l)  =   Eystr_flux_R
       psilaxEy(i,1,1,l)  =   1.d0
       BxlaxEy (i,1,1,l)  =   1.d0

! Flujos Ez


       BylaxEz (i,1,1,l)  =  - Ezstr_flux_R ! minus sign
       BxlaxEz (i,1,1,l)  =   1.d0
       psilaxEz(i,1,1,l)  =   1.d0


! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,1,1,l)     =   Jxstr_R  
       Jylax(i,1,1,l)     =   1.d0
       Jzlax(i,1,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,1,1,l)    =   FDxstr_R  
       FDylax(i,1,1,l)    =   1.d0
       FDzlax(i,1,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,1,1,l)  =   Ftauxstr_R 
       Ftauylax(i,1,1,l)  =   1.d0
       Ftauzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,1,1,l)  =   psistr_R 
       psitauylax(i,1,1,l)  =   1.d0
       psitauzlax(i,1,1,l)  =   1.d0


       phitauxlax(i,1,1,l)  =   phistr_R 
       phitauylax(i,1,1,l)  =   1.d0
       phitauzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,1,1,l)   =   FSxxstr_R 
       FSxylax(i,1,1,l)   =   FSxystr_R 
       FSxzlax(i,1,1,l)   =   FSxzstr_R 

 
       FSyxlax(i,1,1,l)   =   1.d0
       FSyylax(i,1,1,l)   =   1.d0
       FSyzlax(i,1,1,l)   =   1.d0


       FSzxlax(i,1,1,l)  =   1.d0
       FSzylax(i,1,1,l)  =   1.d0
       FSzzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,1,1,l)   = Bxstr_R  
       BySxlax(i,1,1,l)   = 1.d0
       BzSxlax(i,1,1,l)   = 1.d0



       BxSylax(i,1,1,l)   = Bxstr_R  
       BySylax(i,1,1,l)   = 1.d0
       BzSylax(i,1,1,l)   = 1.d0


       BxSzlax(i,1,1,l)   = Bxstr_R  
       BySzlax(i,1,1,l)   = 1.d0
       BzSzlax(i,1,1,l)   = 1.d0

! -------------------------------------- EGLM--------------------------------------

    else if ( s2 .lt. 0.d0 ) then

! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,1,1,l)  =   Exint_R
       Eylaxpsi(i,1,1,l)  =   1.d0
       Ezlaxpsi(i,1,1,l)  =   1.d0

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,1,1,l)  =   Bxint_R
       Bylaxphi(i,1,1,l)  =   1.d0
       Bzlaxphi(i,1,1,l)  =   1.d0

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,1,1,l)  =   phiint_R
       EzlaxBx (i,1,1,l)  =   1.d0
       EylaxBx (i,1,1,l)  =   1.d0

! Flujos By

       EzlaxBy (i,1,1,l)  =  -Ezint_R ! Ezint_R_m  minus sign
       philaxBy(i,1,1,l)  =   1.d0
       ExlaxBy (i,1,1,l)  =   1.d0

! Flujos Bz

       EylaxBz (i,1,1,l)  =   Eyint_R
       ExlaxBz (i,1,1,l)  =   1.d0
       philaxBz(i,1,1,l)  =   1.d0

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,1,1,l)  =   psiint_R
       BzlaxEx(i,1,1,l)   =   1.d0
       BylaxEx(i,1,1,l)   =   1.d0

! Flujos Ey


       BzlaxEy (i,1,1,l)  =   Bzint_R
       psilaxEy(i,1,1,l)  =   1.d0
       BxlaxEy (i,1,1,l)  =   1.d0

! Flujos Ez


       BylaxEz (i,1,1,l)  =  - Byint_R ! Byint_R_m  minus sign
       BxlaxEz (i,1,1,l)  =   1.d0
       psilaxEz(i,1,1,l)  =   1.d0


! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,1,1,l)     =   Jxint_R  
       Jylax(i,1,1,l)     =   1.d0
       Jzlax(i,1,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,1,1,l)    =   FDxint_R  
       FDylax(i,1,1,l)    =   1.d0
       FDzlax(i,1,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,1,1,l)  =   Ftauxint_R 
       Ftauylax(i,1,1,l)  =   1.d0
       Ftauzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,1,1,l)  =   psiint_R 
       psitauylax(i,1,1,l)  =   1.d0
       psitauzlax(i,1,1,l)  =   1.d0


       phitauxlax(i,1,1,l)  =   phiint_R 
       phitauylax(i,1,1,l)  =   1.d0
       phitauzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,1,1,l)   =   FSxxint_R 
       FSxylax(i,1,1,l)   =   FSxyint_R 
       FSxzlax(i,1,1,l)   =   FSxzint_R 

 
       FSyxlax(i,1,1,l)   =   1.d0
       FSyylax(i,1,1,l)   =   1.d0
       FSyzlax(i,1,1,l)   =   1.d0


       FSzxlax(i,1,1,l)  =   1.d0
       FSzylax(i,1,1,l)  =   1.d0
       FSzzlax(i,1,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,1,1,l)   = Bxint_R  
       BySxlax(i,1,1,l)   = 1.d0
       BzSxlax(i,1,1,l)   = 1.d0



       BxSylax(i,1,1,l)   = Bxint_R  
       BySylax(i,1,1,l)   = 1.d0
       BzSylax(i,1,1,l)   = 1.d0


       BxSzlax(i,1,1,l)   = Bxint_R  
       BySzlax(i,1,1,l)   = 1.d0
       BzSzlax(i,1,1,l)   = 1.d0

! -------------------------------------- EGLM-------------------------------------- 
     

!!$   else if (Vxstr .eq. 0.d0 ) then
!!$
!!$
!!$! Electric Gauss HLL flux
!!$!_________________________________________________________________________
!!$
!!$       Exlaxpsi(i,1,1,l)  =   0.5d0 * (psistr_flux_R + psistr_flux_L)
!!$       Eylaxpsi(i,1,1,l)  =   1.d0
!!$       Ezlaxpsi(i,1,1,l)  =   1.d0
!!$
!!$! Magnetic Gauss HLL flux
!!$!_________________________________________________________________________
!!$
!!$
!!$       Bxlaxphi(i,1,1,l)  =   0.5d0 * (phistr_flux_R + phistr_flux_L)
!!$       Bylaxphi(i,1,1,l)  =   1.d0
!!$       Bzlaxphi(i,1,1,l)  =   1.d0
!!$
!!$ ! Faraday HLL flux
!!$!_________________________________________________________________________
!!$
!!$! Flujos Bx 
!!$
!!$       philaxBx(i,1,1,l)  =   0.5d0 * (Bxstr_flux_R + Bxstr_flux_L)
!!$       EzlaxBx (i,1,1,l)  =   1.d0
!!$       EylaxBx (i,1,1,l)  =   1.d0
!!$
!!$! Flujos By
!!$
!!$       EzlaxBy (i,1,1,l)  =  - 0.5d0 * (Bystr_flux_R + Bystr_flux_L) ! minus sign
!!$       philaxBy(i,1,1,l)  =   1.d0
!!$       ExlaxBy (i,1,1,l)  =   1.d0
!!$
!!$! Flujos Bz
!!$
!!$       EylaxBz (i,1,1,l)  =   0.5d0 * (Bzstr_flux_R + Bzstr_flux_L)
!!$       ExlaxBz (i,1,1,l)  =   1.d0
!!$       philaxBz(i,1,1,l)  =   1.d0
!!$
!!$! Ampere Law flux
!!$!_________________________________________________________________________
!!$
!!$! Flujos Ex
!!$
!!$       psilaxEx(i,1,1,l)  =   0.5d0 * (Exstr_flux_R + Exstr_flux_L)
!!$       BzlaxEx(i,1,1,l)   =   1.d0
!!$       BylaxEx(i,1,1,l)   =   1.d0
!!$
!!$! Flujos Ey
!!$
!!$
!!$       BzlaxEy (i,1,1,l)  =   0.5d0 * (Eystr_flux_R + Eystr_flux_L)
!!$       psilaxEy(i,1,1,l)  =   1.d0
!!$       BxlaxEy (i,1,1,l)  =   1.d0
!!$
!!$! Flujos Ez
!!$
!!$
!!$       BylaxEz (i,1,1,l)  =  - 0.5d0 * (Ezstr_flux_R + Ezstr_flux_L) ! minus sign
!!$       BxlaxEz (i,1,1,l)  =   1.d0
!!$       psilaxEz(i,1,1,l)  =   1.d0
!!$
!!$! Conserved current HLL flux
!!$!_________________________________________________________________________
!!$
!!$
!!$       Jxlax(i,1,1,l)     =   0.5d0 * (Jxstr_R + Jxstr_L)
!!$       Jylax(i,1,1,l)     =   1.d0
!!$       Jzlax(i,1,1,l)     =   1.d0
!!$
!!$! Conserved Mass HLL flux
!!$!_________________________________________________________________________
!!$
!!$
!!$       FDxlax(i,1,1,l)    =   0.5d0 * (FDxstr_R + FDxstr_L)
!!$       FDylax(i,1,1,l)    =   1.d0
!!$       FDzlax(i,1,1,l)    =   1.d0
!!$
!!$! Conserved Energy HLL flux
!!$!_________________________________________________________________________
!!$
!!$
!!$       Ftauxlax(i,1,1,l)  =   0.5d0 * (Ftauxstr_R + Ftauxstr_L)
!!$       Ftauylax(i,1,1,l)  =   1.d0
!!$       Ftauzlax(i,1,1,l)  =   1.d0
!!$
!!$! -------------------------------------- EGLM-------------------------------------- 
!!$
!!$       psitauxlax(i,1,1,l)  =   0.5d0 * (psistr_R + psistr_L)
!!$       psitauylax(i,1,1,l)  =   1.d0
!!$       psitauzlax(i,1,1,l)  =   1.d0
!!$
!!$
!!$       phitauxlax(i,1,1,l)  =   0.5d0 * (phistr_R + phistr_L)
!!$       phitauylax(i,1,1,l)  =   1.d0
!!$       phitauzlax(i,1,1,l)  =   1.d0
!!$
!!$! -------------------------------------- EGLM-------------------------------------- 
!!$
!!$! Conserved Momentum HLL flux tensor
!!$!_________________________________________________________________________
!!$
!!$
!!$       FSxxlax(i,1,1,l)   =   0.5d0 * (FSxxstr_R + FSxxstr_L)
!!$       FSxylax(i,1,1,l)   =   0.5d0 * (FSxystr_R + FSxystr_L)
!!$       FSxzlax(i,1,1,l)   =   0.5d0 * (FSxzstr_R + FSxzstr_L)
!!$
!!$ 
!!$       FSyxlax(i,1,1,l)   =   1.d0
!!$       FSyylax(i,1,1,l)   =   1.d0
!!$       FSyzlax(i,1,1,l)   =   1.d0
!!$
!!$
!!$       FSzxlax(i,1,1,l)  =   1.d0
!!$       FSzylax(i,1,1,l)  =   1.d0
!!$       FSzzlax(i,1,1,l)  =   1.d0
!!$
!!$! -------------------------------------- EGLM-------------------------------------- 
!!$
!!$
!!$       BxSxlax(i,1,1,l)   = 0.5d0 * (Bxstr_R  + Bxstr_L)
!!$       BySxlax(i,1,1,l)   = 1.d0
!!$       BzSxlax(i,1,1,l)   = 1.d0
!!$
!!$
!!$
!!$       BxSylax(i,1,1,l)   = 0.5d0 * (Bxstr_R  + Bxstr_L)
!!$       BySylax(i,1,1,l)   = 1.d0
!!$       BzSylax(i,1,1,l)   = 1.d0
!!$
!!$
!!$       BxSzlax(i,1,1,l)   = 0.5d0 * (Bxstr_R  + Bxstr_L)
!!$       BySzlax(i,1,1,l)   = 1.d0
!!$       BzSzlax(i,1,1,l)   = 1.d0
!!$
!!$! -------------------------------------- EGLM-------------------------------------- 
!!$

    else 

       print*, "Revise Vxstar velocity it could be imaginary subroutine hllc_flow "
       print*, "Vxstar =", Vxstr
       print*, "s1=", s1
       print*, "s2=", s2
       stop

    end if

  end do
 
      !///////////////////////////////////////////////////////////////////////////////////////////////////////
      ! 2D
      !///////////////////////////////////////////////////////////////////////////////////////////////////////


    else if (DIM == 2) then


       call boundary_conserved
       call boundary_electric
       call boundary_flux_conserved

             
 
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

             Vx_x_L = Vx(i  ,j  ,1) +  FAC  * mcl(Vx(i+1,j  ,1) - Vx(i  ,j  ,1) ,  &
                                                  Vx(i  ,j  ,1) - Vx(i-1,j  ,1))
             Vx_x_R = Vx(i+1,j  ,1) -  FAC  * mcl(Vx(i+2,j  ,1) - Vx(i+1,j  ,1) ,  &
                                                  Vx(i+1,j  ,1) - Vx(i  ,j  ,1))

             Vy_x_L = Vy(i  ,j  ,1) +  FAC  * mcl(Vy(i+1,j  ,1) - Vy(i  ,j  ,1) ,  &
                                                  Vy(i  ,j  ,1) - Vy(i-1,j  ,1))
             Vy_x_R = Vy(i+1,j  ,1) -  FAC  * mcl(Vy(i+2,j  ,1) - Vy(i+1,j  ,1) ,  &
                                                  Vy(i+1,j  ,1) - Vy(i  ,j  ,1))

             Vz_x_L = Vz(i  ,j  ,1) +  FAC  * mcl(Vz(i+1,j  ,1) - Vz(i  ,j  ,1) ,  &
                                                  Vz(i  ,j  ,1) - Vz(i-1,j  ,1))
             Vz_x_R = Vz(i+1,j  ,1) -  FAC  * mcl(Vz(i+2,j  ,1) - Vz(i+1,j  ,1) ,  &
                                                  Vz(i+1,j  ,1) - Vz(i  ,j  ,1))

             Vx_y_L = Vx(i  ,j  ,1) +  FAC  * mcl(Vx(i  ,j+1,1) - Vx(i  ,j  ,1) ,  &
                                                  Vx(i  ,j  ,1) - Vx(i  ,j-1,1))
             Vx_y_R = Vx(i  ,j+1,1) -  FAC  * mcl(Vx(i  ,j+2,1) - Vx(i  ,j+1,1) ,  &
                                                  Vx(i  ,j+1,1) - Vx(i  ,j  ,1))

             Vy_y_L = Vy(i  ,j  ,1) +  FAC  * mcl(Vy(i  ,j+1,1) - Vy(i  ,j  ,1) ,  &
                                                  Vy(i  ,j  ,1) - Vy(i  ,j-1,1))
             Vy_y_R = Vy(i  ,j+1,1) -  FAC  * mcl(Vy(i  ,j+2,1) - Vy(i  ,j+1,1) ,  &
                                                  Vy(i  ,j+1,1) - Vy(i  ,j  ,1))

             Vz_y_L = Vz(i  ,j  ,1) +  FAC  * mcl(Vz(i  ,j+1,1) - Vz(i  ,j  ,1) ,  &
                                                  Vz(i  ,j  ,1) - Vz(i  ,j-1,1))
             Vz_y_R = Vz(i  ,j+1,1) -  FAC  * mcl(Vz(i  ,j+2,1) - Vz(i  ,j+1,1) ,  &
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

             Vx_x_L = Vxint(i  ,j  ,1,l) +  FAC  * mcl(Vxint(i+1,j  ,1,l) - Vxint(i  ,j  ,1,l) ,  &
                                                       Vxint(i  ,j  ,1,l) - Vxint(i-1,j  ,1,l))
             Vx_x_R = Vxint(i+1,j  ,1,l) -  FAC  * mcl(Vxint(i+2,j  ,1,l) - Vxint(i+1,j  ,1,l) ,  &
                                                       Vxint(i+1,j  ,1,l) - Vxint(i  ,j  ,1,l))

             Vy_x_L = Vyint(i  ,j  ,1,l) +  FAC  * mcl(Vyint(i+1,j  ,1,l) - Vyint(i  ,j  ,1,l) ,  &
                                                       Vyint(i  ,j  ,1,l) - Vyint(i-1,j  ,1,l))
             Vy_x_R = Vyint(i+1,j  ,1,l) -  FAC  * mcl(Vyint(i+2,j  ,1,l) - Vyint(i+1,j  ,1,l) ,  &
                                                       Vyint(i+1,j  ,1,l) - Vyint(i  ,j  ,1,l))

             Vz_x_L = Vzint(i  ,j  ,1,l) +  FAC  * mcl(Vzint(i+1,j  ,1,l) - Vzint(i  ,j  ,1,l) ,  &
                                                       Vzint(i  ,j  ,1,l) - Vzint(i-1,j  ,1,l))
             Vy_x_R = Vzint(i+1,j  ,1,l) -  FAC  * mcl(Vzint(i+2,j  ,1,l) - Vzint(i+1,j  ,1,l) ,  &
                                                       Vzint(i+1,j  ,1,l) - Vzint(i  ,j  ,1,l))


             Vx_y_L = Vxint(i  ,j  ,1,l) +  FAC  * mcl(Vxint(i  ,j+1,1,l) - Vxint(i  ,j  ,1,l) ,  &
                                                       Vxint(i  ,j  ,1,l) - Vxint(i  ,j-1,1,l))
             Vx_y_R = Vxint(i  ,j+1,1,l) -  FAC  * mcl(Vxint(i  ,j+2,1,l) - Vxint(i  ,j+1,1,l) ,  &
                                                       Vxint(i  ,j+1,1,l) - Vxint(i  ,j  ,1,l))

             Vy_y_L = Vyint(i  ,j  ,1,l) +  FAC  * mcl(Vyint(i  ,j+1,1,l) - Vyint(i  ,j  ,1,l) ,  &
                                                       Vyint(i  ,j  ,1,l) - Vyint(i  ,j-1,1,l))
             Vy_y_R = Vyint(i  ,j+1,1,l) -  FAC  * mcl(Vyint(i  ,j+2,1,l) - Vyint(i  ,j+1,1,l) ,  &
                                                       Vyint(i  ,j+1,1,l) - Vyint(i  ,j  ,1,l))

             Vz_y_L = Vzint(i  ,j  ,1,l) +  FAC  * mcl(Vzint(i  ,j+1,1,l) - Vzint(i  ,j  ,1,l) ,  &
                                                       Vzint(i  ,j  ,1,l) - Vzint(i  ,j-1,1,l))
             Vy_y_R = Vzint(i  ,j+1,1,l) -  FAC  * mcl(Vzint(i  ,j+2,1,l) - Vzint(i  ,j+1,1,l) ,  &
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
!!$
!!$
!!$
!!$    Jxint_L = Jxint(i  ,j  ,1,l)  +  FAC * mcl(Jxint(i+1,j  ,1,l) - Jxint(i  ,j  ,1,l) ,  &
!!$                                               Jxint(i  ,j  ,1,l) - Jxint(i-1,j  ,1,l)) 
!!$    Jxint_R = Jxint(i+1,j  ,1,l)  -  FAC * mcl(Jxint(i+2,j  ,1,l) - Jxint(i+1,j  ,1,l) ,  &
!!$                                               Jxint(i+1,j  ,1,l) - Jxint(i  ,j  ,1,l)) 
!!$    Jyint_L = Jyint(i  ,j  ,1,l)  +  FAC * mcl(Jyint(i  ,j+1,1,l) - Jyint(i  ,j  ,1,l) ,  &
!!$                                               Jyint(i  ,j  ,1,l) - Jyint(i  ,j-1,1,l))
!!$    Jyint_R = Jyint(i  ,j+1,1,l)  -  FAC * mcl(Jyint(i  ,j+2,1,l) - Jyint(i  ,j+1,1,l) ,  &
!!$                                               Jyint(i  ,j+1,1,l) - Jyint(i  ,j  ,1,l))
!!$    Jzint_L = Jzint(i  ,j  ,1,l)  +  FAC * mcl(Jzint(i+1,j  ,1,l) - Jzint(i  ,j  ,1,l) ,  &
!!$                                               Jzint(i  ,j  ,1,l) - Jzint(i-1,j  ,1,l))
!!$    Jzint_R = Jzint(i+1,j  ,1,l)  -  FAC * mcl(Jzint(i+2,j  ,1,l) - Jzint(i+1,j  ,1,l) ,  &
!!$                                               Jzint(i+1,j  ,1,l) - Jzint(i  ,j  ,1,l)) 

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

!!$    p_x_L   = p  (i  ,j  ,1) + FAC  * mcl(p  (i+1,j  ,1) - p  (i  ,j  ,1) ,  &
!!$                                          p  (i  ,j  ,1) - p  (i-1,j  ,1))
!!$    p_x_R   = p  (i+1,j  ,1) - FAC  * mcl(p  (i+2,j  ,1) - p  (i+1,j  ,1) ,  &
!!$                                          p  (i+1,j  ,1) - p  (i  ,j  ,1))
!!$
!!$    rho_x_L = rho(i  ,j  ,1) + FAC  * mcl(rho(i+1,j  ,1) - rho(i  ,j  ,1) ,  &
!!$                                          rho(i  ,j  ,1) - rho(i-1,j  ,1))
!!$    rho_x_R = rho(i+1,j  ,1) - FAC  * mcl(rho(i+2,j  ,1) - rho(i+1,j  ,1) ,  &
!!$                                          rho(i+1,j  ,1) - rho(i  ,j  ,1))
!!$
!!$    p_y_L   = p  (i  ,j  ,1) + FAC  * mcl(p  (i  ,j+1,1) - p  (i  ,j  ,1) ,  &
!!$                                          p  (i  ,j  ,1) - p  (i  ,j-1,1))
!!$    p_y_R   = p  (i  ,j+1,1) - FAC  * mcl(p  (i  ,j+2,1) - p  (i  ,j+1,1) ,  &
!!$                                          p  (i  ,j+1,1) - p  (i  ,j  ,1))
!!$
!!$    rho_y_L = rho(i  ,j  ,1) + FAC  * mcl(rho(i  ,j+1,1) - rho(i  ,j  ,1) ,  &
!!$                                          rho(i  ,j  ,1) - rho(i  ,j-1,1))
!!$    rho_y_R = rho(i  ,j+1,1) - FAC  * mcl(rho(i  ,j+2,1) - rho(i  ,j+1,1) ,  &
!!$                                          rho(i  ,j+1,1) - rho(i  ,j  ,1))
!!$    
!!$
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


    if (p_x_L .lt. 0 .or. p_x_R .lt. 0) then 

       print*, "negative pressure in HLL ---> p_L =", p_x_L, "p_R =", p_x_R, "i,j =", i, j

    end if

    if (rho_x_L .lt. 0 .or. rho_x_R .lt. 0) then 

       print*, "negative density in HLL ---> rho_L =", rho_L, "rho_R =", rho_R, "i =", i, j

    end if
    


!!$       enthpy_x_L   =  enthpyint(i  ,j  ,1,l) +  FAC  * mcl(enthpyint(i+1,j  ,1,l) - enthpyint(i  ,j  ,1,l) ,  &
!!$                                                               enthpyint(i  ,j  ,1,l) - enthpyint(i-1,j  ,1,l))
!!$       enthpy_x_R   =  enthpyint(i+1,j  ,1,l) -  FAC  * mcl(enthpyint(i+2,j  ,1,l) - enthpyint(i+1,j  ,1,l) ,  &
!!$                                                               enthpyint(i+1,j  ,1,l) - enthpyint(i  ,j  ,1,l))
!!$
!!$       enthpy_y_L   =  enthpyint(i  ,j  ,1,l) +  FAC  * mcl(enthpyint(i  ,j+1,1,l) - enthpyint(i  ,j  ,1,l) ,  &
!!$                                                               enthpyint(i  ,j  ,1,l) - enthpyint(i  ,j-1,1,l))
!!$       enthpy_y_R   =  enthpyint(i  ,j+1,1,l) -  FAC  * mcl(enthpyint(i  ,j+2,1,l) - enthpyint(i  ,j+1,1,l) ,  &
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


!/////////////////////// HLL Variables  X DIRECTION ////////////////////////////////////

! we difine de characteristic velocitis in diferent form from notes as -----> s2 = lambda_r  and s1 = lambda_l

    !     hll  values of augmented fields
    psi_x_hll = ( s2 * psiint_x_R - s1 * psiint_x_L + Exint_x_L  - Exint_x_R   ) / (s2 -s1)
    phi_x_hll = ( s2 * phiint_x_R - s1 * phiint_x_L + Bxint_x_L  - Bxint_x_R   ) / (s2 -s1)

    !     hll  values of electric field

    Ex_x_hll  = ( s2 * Exint_x_R  - s1 * Exint_x_L  + psiint_x_L  - psiint_x_R  ) / (s2 -s1)
    Ey_x_hll  = ( s2 * Eyint_x_R  - s1 * Eyint_x_L  + Bzint_x_L   - Bzint_x_R   ) / (s2 -s1)
    Ez_x_hll  = ( s2 * Ezint_x_R  - s1 * Ezint_x_L  + Byint_x_L_m - Byint_x_R_m ) / (s2 -s1)  !ojo minus sign ¡¡¡


    !     hll  values of magnetic field

    Bx_x_hll  = ( s2 * Bxint_x_R  - s1 * Bxint_x_L  + phiint_x_L  - phiint_x_R  ) / (s2 -s1)
    By_x_hll  = ( s2 * Byint_x_R  - s1 * Byint_x_L  + Ezint_x_L_m - Ezint_x_R_m ) / (s2 -s1)  !ojo minus sign ¡¡¡
    Bz_x_hll  = ( s2 * Bzint_x_R  - s1 * Bzint_x_L  + Eyint_x_L   - Eyint_x_R   ) / (s2 -s1)


    !     hll  value of charge density 

    q_x_hll   = ( s2 * qint_x_R   - s1 * qint_x_L   + Jxint_L   - Jxint_R   ) / (s2 -s1)

    !     hll  value of mass conservation 

    DD_x_hll  = ( s2 * DDint_x_R  - s1 * DDint_x_L  + FDxint_L  - FDxint_R  ) / (s2 -s1)

    !     hll  value of energy conservation 

    tau_x_hll = ( s2 * tauint_x_R - s1 * tauint_x_L + Ftauxint_L - Ftauxint_R ) / (s2 -s1)

    !    hll  value of momentum S

    Sx_x_hll  = ( s2 * Sxint_x_R  - s1 * Sxint_x_L  + FSxxint_L  - FSxxint_R  ) / (s2 -s1)
    Sy_x_hll  = ( s2 * Syint_x_R  - s1 * Syint_x_L  + FSxyint_L  - FSxyint_R  ) / (s2 -s1)
    Sz_x_hll  = ( s2 * Szint_x_R  - s1 * Szint_x_L  + FSxzint_L  - FSxzint_R  ) / (s2 -s1)


!/////////////////////// HLL Fluxes  X DIRECTION ////////////////////////////////////

! Electric Gauss HLL flux
!_________________________________________________________________________

      psi_x_hll_flux  =   0.5d0 * (Exint_x_L + Exint_x_R - ( psiint_x_R - psiint_x_L ) )  

! Magnetic Gauss HLL flux
!_________________________________________________________________________


      phi_x_hll_flux  =   0.5d0 * (Bxint_x_L + Bxint_x_R - ( phiint_x_R - phiint_x_L ) )

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       Bx_x_hll_flux  =   0.5d0 * (phiint_x_L   + phiint_x_R   - ( Bxint_x_R - Bxint_x_L ) )

! Flujos By

       By_x_hll_flux  =   0.5d0 * (Ezint_x_L_m  + Ezint_x_R_m  - ( Byint_x_R - Byint_x_L ) ) ! minus sign

! Flujos Bz

       Bz_x_hll_flux  =   0.5d0 * (Eyint_x_L    + Eyint_x_R    - ( Bzint_x_R - Bzint_x_L ) )

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       Ex_x_hll_flux  =   0.5d0 * (psiint_x_L   + psiint_x_R   - ( Exint_x_R - Exint_x_L ) )

! Flujos Ey


       Ey_x_hll_flux  =   0.5d0 * (Bzint_x_L    + Bzint_x_R    - ( Eyint_x_R - Eyint_x_L ) )

! Flujos Ez


       Ez_x_hll_flux  =   0.5d0 * (Byint_x_L_m  + Byint_x_R_m    - ( Ezint_x_R - Ezint_x_L ) ) ! minus sign

! Conserved current HLL flux
!_________________________________________________________________________


       q_x_hll_flux   =   0.5d0 * (Jxint_L  + Jxint_R - ( qint_x_R - qint_x_L ) )

! Conserved Mass HLL flux
!_________________________________________________________________________


       DD_x_hll_flux  =   0.5d0 * (FDxint_L  + FDxint_R - ( DDint_x_R - DDint_x_L ) )

! Conserved Energy HLL flux
!_________________________________________________________________________


       tau_x_hll_flux =   0.5d0 * (Ftauxint_L + Ftauxint_R - ( tauint_x_R - tauint_x_L  ) )


! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       Sx_x_hll_flux  =   0.5d0 * (FSxxint_L + FSxxint_R - (Sxint_x_R - Sxint_x_L ) ) 
       Sy_x_hll_flux  =   0.5d0 * (FSxyint_L + FSxyint_R - (Syint_x_R - Syint_x_L ) )  
       Sz_x_hll_flux  =   0.5d0 * (FSxzint_L + FSxzint_R - (Szint_x_R - Szint_x_L ) )


       
!/////////////////////// HLL Variables  Y DIRECTION ////////////////////////////////////

! we difine de characteristic velocitis in diferent form from notes as -----> s2 = lambda_r  and s1 = lambda_l

    !     hll  values of augmented fields
    psi_y_hll = ( s2 * psiint_y_R - s1 * psiint_y_L + Eyint_y_L  - Eyint_y_R   ) / (s2 -s1)
    phi_y_hll = ( s2 * phiint_y_R - s1 * phiint_y_L + Byint_y_L  - Byint_y_R   ) / (s2 -s1)

    !     hll  values of electric field

    Ey_y_hll  = ( s2 * Eyint_y_R  - s1 * Eyint_y_L  + psiint_y_L  - psiint_y_R  ) / (s2 -s1)
    Ez_y_hll  = ( s2 * Ezint_y_R  - s1 * Ezint_y_L  + Bxint_y_L   - Bxint_y_R   ) / (s2 -s1)
    Ex_y_hll  = ( s2 * Exint_y_R  - s1 * Exint_y_L  + Bzint_y_L_m - Bzint_y_R_m ) / (s2 -s1)  !ojo minus sign ¡¡¡


    !     hll  values of magnetic field

    By_y_hll  = ( s2 * Byint_y_R  - s1 * Byint_y_L  + phiint_y_L  - phiint_y_R  ) / (s2 -s1)
    Bz_y_hll  = ( s2 * Bzint_y_R  - s1 * Bzint_y_L  + Exint_y_L_m - Exint_y_R_m ) / (s2 -s1)  !ojo minus sign ¡¡¡
    Bx_y_hll  = ( s2 * Bxint_y_R  - s1 * Bxint_y_L  + Ezint_y_L   - Ezint_y_R   ) / (s2 -s1)


    !     hll  value of charge density 

    q_y_hll   = ( s2 * qint_y_R   - s1 * qint_y_L   + Jyint_L   - Jyint_R   ) / (s2 -s1)

    !     hll  value of mass conservation 

    DD_y_hll  = ( s2 * DDint_y_R  - s1 * DDint_y_L  + FDyint_L  - FDyint_R  ) / (s2 -s1)

    !     hll  value of energy conservation 

    tau_y_hll = ( s2 * tauint_y_R - s1 * tauint_y_L + Ftauyint_L - Ftauyint_R ) / (s2 -s1)

    !    hll  value of momentum S

    Sy_y_hll  = ( s2 * Syint_y_R  - s1 * Syint_y_L  + FSyyint_L  - FSyyint_R  ) / (s2 -s1)
    Sz_y_hll  = ( s2 * Szint_y_R  - s1 * Szint_y_L  + FSyzint_L  - FSyzint_R  ) / (s2 -s1)
    Sx_y_hll  = ( s2 * Sxint_y_R  - s1 * Sxint_y_L  + FSyxint_L  - FSyxint_R  ) / (s2 -s1)


!/////////////////////// HLL Fluxes  Y DIRECTION ////////////////////////////////////

! Electric Gauss HLL flux
!_________________________________________________________________________

      psi_y_hll_flux  =   0.5d0 * (Eyint_y_L + Eyint_y_R - ( psiint_y_R - psiint_y_L ) )  

! Magnetic Gauss HLL flux
!_________________________________________________________________________


      phi_y_hll_flux  =   0.5d0 * (Byint_y_L + Byint_y_R - ( phiint_y_R - phiint_y_L ) )

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos By 

       By_y_hll_flux  =   0.5d0 * (phiint_y_L   + phiint_y_R   - ( Byint_y_R - Byint_y_L ) )

! Flujos Bz

       Bz_y_hll_flux  =   0.5d0 * (Exint_y_L_m  + Exint_y_R_m  - ( Bzint_y_R - Bzint_y_L ) ) ! minus sign

! Flujos Bx

       Bx_y_hll_flux  =   0.5d0 * (Ezint_y_L    + Ezint_y_R    - ( Bxint_y_R - Bxint_y_L ) )

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ey

       Ey_y_hll_flux  =   0.5d0 * (psiint_y_L   + psiint_y_R   - ( Eyint_y_R - Eyint_y_L ) )

! Flujos Ez


       Ez_y_hll_flux  =   0.5d0 * (Bxint_y_L    + Bxint_y_R    - ( Ezint_y_R - Ezint_y_L ) )

! Flujos Ex


       Ex_y_hll_flux  =   0.5d0 * (Bzint_y_L_m  + Bzint_y_R_m    - ( Exint_y_R - Exint_y_L ) ) ! minus sign

! Conserved current HLL flux
!_________________________________________________________________________


       q_y_hll_flux   =   0.5d0 * (Jyint_L  + Jyint_R - ( qint_y_R - qint_y_L ) )

! Conserved Mass HLL flux
!_________________________________________________________________________


       DD_y_hll_flux  =   0.5d0 * (FDyint_L  + FDyint_R - ( DDint_y_R - DDint_y_L ) )

! Conserved Energy HLL flux
!_________________________________________________________________________


       tau_y_hll_flux =   0.5d0 * (Ftauyint_L + Ftauyint_R - ( tauint_y_R - tauint_y_L  ) )


! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       Sy_y_hll_flux  =   0.5d0 * (FSyyint_L + FSyyint_R - (Syint_y_R - Syint_y_L ) ) 
       Sz_y_hll_flux  =   0.5d0 * (FSyzint_L + FSyzint_R - (Szint_y_R - Szint_y_L ) ) 
       Sx_y_hll_flux  =   0.5d0 * (FSyxint_L + FSyxint_R - (Sxint_y_R - Sxint_y_L ) )

! Dedudction of HLLC solver M.A Aloy 08-2017

! 1) First calculate the continuos conserved variables in start region (q, E, B) as its avarage HLL values


       ! X direction 
       
   Exstr_x_R  = Ex_x_hll  
   Exstr_x_L  = Ex_x_hll
   
   Eystr_x_R  = Ey_x_hll
   Eystr_x_L  = Ey_x_hll 

   Ezstr_x_R  = Ez_x_hll 
   Ezstr_x_L  = Ez_x_hll 

   Bxstr_x_R  = Bx_x_hll
   Bxstr_x_L  = Bx_x_hll 
   
   Bystr_x_R  = By_x_hll 
   Bystr_x_L  = By_x_hll 

   Bzstr_x_R  = Bz_x_hll 
   Bzstr_x_L  = Bz_x_hll 

   qstr_x_R   = q_x_hll  
   qstr_x_L   = q_x_hll

          ! Y direction 
       
   Exstr_y_R  = Ex_y_hll  
   Exstr_y_L  = Ex_y_hll
   
   Eystr_y_R  = Ey_y_hll
   Eystr_y_L  = Ey_y_hll 

   Ezstr_y_R  = Ez_y_hll 
   Ezstr_y_L  = Ez_y_hll 

   Bxstr_y_R  = Bx_y_hll
   Bxstr_y_L  = Bx_y_hll 
   
   Bystr_y_R  = By_y_hll 
   Bystr_y_L  = By_y_hll 

   Bzstr_y_R  = Bz_y_hll 
   Bzstr_y_L  = Bz_y_hll 

   qstr_y_R   = q_y_hll  
   qstr_y_L   = q_y_hll

   !  2) then find the contac wave velocitu vx^* by means of quadratic equation
   !                    a (vx^*)^2 + b vx^* + c = 0
   !  with:
   !           a = tau_hll_flux - ErotB_str_x
   !           b = E2_perp + B2_perp - tau_hll - Sx_hll_flux
   !           c = Sx_hll - ErotB_str_x
       
   ErotB_str_x_x_L = (Bzstr_x_L * Eystr_x_L - Bystr_x_L * Ezstr_x_L)
   ErotB_str_y_x_L = (Bxstr_x_L * Ezstr_x_L - Bzstr_x_L * Exstr_x_L)
   ErotB_str_z_x_L = (Bystr_x_L * Exstr_x_L - Bxstr_x_L * Eystr_x_L)

   ErotB_str_x_x_R = (Bzstr_x_R * Eystr_x_R - Bystr_x_R * Ezstr_x_R)
   ErotB_str_y_x_R = (Bxstr_x_R * Ezstr_x_R - Bzstr_x_R * Exstr_x_R)
   ErotB_str_z_x_R = (Bystr_x_R * Exstr_x_R - Bxstr_x_R * Eystr_x_R)

   ErotB_str_x_y_L = (Bzstr_y_L * Eystr_y_L - Bystr_y_L * Ezstr_y_L)
   ErotB_str_y_y_L = (Bxstr_y_L * Ezstr_y_L - Bzstr_y_L * Exstr_y_L)
   ErotB_str_z_y_L = (Bystr_y_L * Exstr_y_L - Bxstr_y_L * Eystr_y_L)

   ErotB_str_x_y_R = (Bzstr_y_R * Eystr_y_R - Bystr_y_R * Ezstr_y_R)
   ErotB_str_y_y_R = (Bxstr_y_R * Ezstr_y_R - Bzstr_y_R * Exstr_y_R)
   ErotB_str_z_y_R = (Bystr_y_R * Exstr_y_R - Bxstr_y_R * Eystr_y_R)



   ! FIND ROOT Vxstr !

   E2_perp_x   = (Eystr_x_R**2 + Ezstr_x_R**2)
   B2_perp_x   = (Bystr_x_R**2 + Bzstr_x_R**2)

   E2B2_str_x  = 0.5d0 * ( (Exstr_x_R**2 + Eystr_x_R**2 + Ezstr_x_R**2) &
                         + (Bxstr_x_R**2 + Bystr_x_R**2 + Bzstr_x_R**2) )

   a_hll_x     = tau_x_hll_flux - ErotB_str_x_x_R
   b_hll_x     = E2_perp_x      + B2_perp_x - tau_x_hll - Sx_x_hll_flux
   c_hll_x     = Sx_x_hll       - ErotB_str_x_x_R


   ! Analitic solution quadratic equations


   eps  = 1.d-10
   eps1 = 1.d-6
   eps2 = 1.d-3

   
   det_hll     = b_hll_x**2 - 4.d0 * a_hll_x * c_hll_x

   q_root1  = - 0.5d0 * ( b_hll_x + sign(1.d0,b_hll_x) * abs(sqrt( det_hll )) )

   rt1_x    = q_root1 / a_hll_x
   rt2_x    = c_hll_x  / q_root1 ! this is the negative root with has physical meaning

   if (det_hll .lt. 0.d0) then

      print*, "discriminante imaginario en la ecuacion cuadratica para obtener Vxstr en HLLC Riemman solver (x direction)"
      print*, "det_hll =", det_hll

      rt2_x    = c_hll_x  / b_hll_x

   end if


!!$   if ( abs(rt2_x) .le. 1.d0 ) then
!!$
!!$      if ( abs(rt2_x) .le. 1.d-6 ) rt2_x = 0.d0
!!$
!!$      Vxstr = rt2_x
!!$
!!$   else if ( abs(rt2_x) .gt. 1.d0 ) then
!!$
!!$      print*, "superluminal velocity ( -1 > Vxstr > 1) in  quadratic equation HLLC Solver"
!!$      print*, "rt1_x =", rt1_x
!!$      print*, "rt2_x =", rt2_x
!!$
!!$   end if

   
   ! Numerical Solution quadratic equations

!!$Para que \lambda^* sea igual a cero en la ecuación cuadrática (42), debe suceder que c=0.
!!$
!!$Sin embargo, incluso en el caso c=0, la ecuación cuadrática (42) tiene dos soluciones:
!!$
!!$\lambda^*_a =0
!!$\lambda^*_b = -b/a (suponiendo que a \ne 0)
!!$
!!$Nótese que las soluciones generales de la cuadrática son:
!!$
!!$\lambda^*_1 = [-b - \sqrt( b^2 - 4ac ) ] / (2a)
!!$\lambda^*_2 = [-b + \sqrt( b^2 - 4ac ) ] / (2a)
!!$
!!$De esas dos soluciones, se toma sólo \lambda^*_1, pues como razonan Mignone & Bodo (2005, 2006), esa es la única que garantiza que
!!$
!!$\lambda_l \le \lambda^* \le \lambda_r.  [A]
!!$
!!$Fijémonos ahora que, cuando c -> 0, la única solución que se anula de las dos en el caso general es la segunda. Es decir, que en caso c->0, tenemos:
!!$
!!$
!!$\lambda^*_1 -> \lambda^*_b = -b/a
!!$\lambda^*_2 -> \lambda^*_a = 0
!!$
!!$Como hemos dicho que la solución \lambda^*_2, no garantiza que se cumpla la condición [A], concluimos que:

!!$   eps  = 1.d-13   
!!$
!!$   if ( abs(c_hll_x) .le. eps ) then
!!$
!!$      if ( abs(a_hll_x) .ge. eps ) then
!!$
!!$        Vxstr  = - b_hll_x / a_hll_x 
!!$
!!$   if ( abs(Vxstr) .le. 1.d-9 ) then
!!$
!!$      Vxstr = 0.d0
!!$
!!$   else if ( abs(Vxstr) .gt. 1.d0 ) then
!!$
!!$      print*, "superluminal velocity ( -1 > Vxstr > 1) in  quadratic equation HLLC Solver"
!!$      print*, "Vxstr =", Vxstr
!!$
!!$      Vxstr= max(-1.d0 + eps1 , min( 1.d0 - eps1, Vxstr ))
!!$
!!$   end if   
!!$
!!$  else
!!$!*
!!$!*       Este caso sólo se da cuando a~0 y c~0
!!$!*
!!$          Vxstr  = 0.d0
!!$
!!$       end if
!!$
!!$    else

!!$       [solución numérica de la ecuación cuadrática]
!!$
!!$   This function obtains the root of the quadratic equation
!!$
!!$       a*v^2 - b*v + c = 0
!!$
!!$   in the physical range [-1,1] usign a Newton-Rapson iterative method.


   Vxstr       = get_hllc_vel(a_hll_x,-b_hll_x,c_hll_x,rt2_x)


   if ( abs(Vxstr) .le. 1.d-6 ) then

      Vxstr = 0.d0

   else if ( abs(Vxstr) .gt. 1.d0 ) then

      print*, "superluminal velocity ( -1 > Vxstr > 1) in  quadratic equation HLLC Solver"
      print*, "Vxstr =", Vxstr

      Vxstr= max(-1.d0 + eps1 , min( 1.d0 - eps1, Vxstr ))


   end if

!!$   end if   


   ! FIND ROOT Vystr !


   E2_perp_y   = (Exstr_y_R**2 + Ezstr_y_R**2)
   B2_perp_y   = (Bxstr_y_R**2 + Bzstr_y_R**2)

   E2B2_str_y  = 0.5d0 * ( (Exstr_y_R**2 + Eystr_y_R**2 + Ezstr_y_R**2) &
                         + (Bxstr_y_R**2 + Bystr_y_R**2 + Bzstr_y_R**2) )

   a_hll_y     = tau_y_hll_flux - ErotB_str_y_y_R
   b_hll_y     = E2_perp_y      + B2_perp_y - tau_y_hll - Sy_y_hll_flux
   c_hll_y     = Sy_y_hll       - ErotB_str_y_y_R

      ! Analitic solution quadratic equations
   
   det_hll     = b_hll_y**2 - 4.d0 * a_hll_y * c_hll_y

   q_root1  = - 0.5d0 * ( b_hll_y + sign(1.d0,b_hll_y) * abs(sqrt( det_hll )) )

   rt1_y    = q_root1 / a_hll_y
   rt2_y    = c_hll_y / q_root1 ! this is the negative root with has physical meaning


   if (det_hll .lt. 0.d0) then

      print*, "discriminante imaginario en la ecuacion cuadratica para obtener Vxstr en HLLC Riemman solver (y direction)"
      print*, "det_hll =", det_hll
      
      rt2_y    = c_hll_y / b_hll_y

    end if



!!$   if ( abs(rt2_y) .le. 1.d0 ) then
!!$
!!$      if ( abs(rt2_y) .le. 1.d-6 ) rt2_y = 0.d0
!!$
!!$      Vystr = rt2_y
!!$
!!$   else if ( abs(rt2_y) .gt. 1.d0 ) then
!!$
!!$      print*, "superluminal velocity ( -1 > Vystr > 1) in  quadratic equation HLLC Solver"
!!$      print*, "rt1_y =", rt1_y
!!$      print*, "rt2_y =", rt2_y
!!$
!!$   end if
 
   ! Numerical Solution quadratic equations

!!$   eps = 1d-13   
!!$
!!$
!!$   if ( abs(c_hll_y) .le. eps ) then
!!$
!!$      if ( abs(a_hll_y) .ge. eps ) then
!!$
!!$        Vystr  = - b_hll_y / a_hll_y 
!!$
!!$   if ( abs(Vystr) .le. 1.d-9 ) then
!!$
!!$      Vystr = 0.d0
!!$
!!$   else if ( abs(Vystr) .gt. 1.d0 ) then
!!$
!!$      print*, "superluminal velocity ( -1 > Vystr > 1) in  quadratic equation HLLC Solver"
!!$      print*, "Vystr =", Vystr
!!$
!!$      Vystr= max(-1.d0 + eps1 , min( 1.d0 - eps1, Vystr ))
!!$
!!$   end if   
!!$
!!$    else
!!$!*
!!$!*       Este caso sólo se da cuando a~0 y c~0
!!$!*
!!$          Vystr  = 0.d0
!!$
!!$       end if
!!$
!!$    else

!!$       [solución numérica de la ecuación cuadrática]
!!$
!!$   This function obtains the root of the quadratic equation
!!$
!!$       a*v^2 - b*v + c = 0
!!$
!!$   in the physical range [-1,1] usign a Newton-Rapson iterative method.


   Vystr       = get_hllc_vel(a_hll_y,-b_hll_y,c_hll_y,rt2_y)


   if ( abs(Vystr) .le. 1.d-6 ) then

      Vystr = 0.d0

   else if ( abs(Vystr) .gt. 1.d0 ) then

      print*, "superluminal velocity ( -1 > Vystr > 1) in  quadratic equation HLLC Solver"
      print*, "Vystr =", Vystr

      Vystr= max(-1.d0 + eps1 , min( 1.d0 - eps1, Vystr ))

   end if   

!!$ end if

   ! 3)  Find the total pressure (P = p_{gas} + p_{em}) in start region

   ! X direction

   pstr_x = Sx_x_hll_flux + (Exstr_x_R**2 + Bxstr_x_R**2 ) - ( tau_x_hll_flux - ErotB_str_x_x_R ) * Vxstr

   ! Y direction

   pstr_y = Sy_y_hll_flux + (Eystr_y_R**2 + Bystr_y_R**2 ) - ( tau_y_hll_flux - ErotB_str_y_y_R ) * Vystr

   ! 4) Calculate the remainder conserved variables in its right and left star values

   ! X direction
      
   psistr_x_R = ( s2 * psi_x_hll - Vxstr * Ex_x_hll_flux) / (s2 - Vxstr) 
   psistr_x_L = (-s1 * psi_x_hll + Vxstr * Ex_x_hll_flux) / (Vxstr - s1) 
   
   phistr_x_R = ( s2 * phi_x_hll - Vxstr * Bx_x_hll_flux) / (s2 - Vxstr) 
   phistr_x_L = (-s1 * phi_x_hll + Vxstr * Bx_x_hll_flux) / (Vxstr - s1)

   Sxstr_x_R  = ( s2 * Sx_x_hll - tau_x_hll_flux * Vxstr ) / (s2 - Vxstr) 
   Sxstr_x_L  = (-s1 * Sx_x_hll + tau_x_hll_flux * Vxstr ) / (Vxstr - s1)

!!$   Sxstr_x_R  = (  s2 * Sx_x_hll - Sx_x_hll_flux + pstr_x - Bxstr_x_R**2 - Exstr_x_R**2 - ErotB_str_x_x_R * Vxstr ) &
!!$              / (s2 - Vxstr) 
!!$   Sxstr_x_L  = ( -s1 * Sx_x_hll + Sx_x_hll_flux - pstr_x + Bxstr_x_R**2 + Exstr_x_R**2 + ErotB_str_x_x_R * Vxstr ) &
!!$              / (Vxstr - s1) 
    
   Systr_x_R  = ( s2 * Sy_x_hll - Sy_x_hll_flux - ErotB_str_y_x_R * Vxstr - Exstr_x_R * Eystr_x_R - Bxstr_x_R * Bystr_x_R )  &
              / ( s2 - Vxstr)
   Systr_x_L  = (-s1 * Sy_x_hll + Sy_x_hll_flux + ErotB_str_y_x_R * Vxstr + Exstr_x_L * Eystr_x_L + Bxstr_x_L * Bystr_x_L ) &
              / ( Vxstr - s1)

   Szstr_x_R  = ( s2 * Sz_x_hll - Sz_x_hll_flux - ErotB_str_z_x_R * Vxstr - Exstr_x_R * Ezstr_x_R - Bxstr_x_R * Bzstr_x_R )  &
              / ( s2 - Vxstr)
   Szstr_x_L  = (-s1 * Sz_x_hll + Sz_x_hll_flux + ErotB_str_z_x_R * Vxstr + Exstr_x_L * Ezstr_x_L + Bxstr_x_L * Bzstr_x_L ) &
              / ( Vxstr - s1)

   taustr_x_R = (Sxstr_x_R - ErotB_str_x_x_R) / Vxstr - pstr_x + E2B2_str_x 
   taustr_x_L = (Sxstr_x_L - ErotB_str_x_x_R) / Vxstr - pstr_x + E2B2_str_x

!!$
!!$   taustr_x_R = (  s2 * tau_x_hll - tau_x_hll_flux + (pstr_x - E2B2_str_x) * Vxstr + ErotB_str_x_x_R ) &
!!$              / (s2 - Vxstr)
!!$   taustr_x_L = ( -s1 * tau_x_hll + tau_x_hll_flux - (pstr_x - E2B2_str_x) * Vxstr - ErotB_str_x_x_R ) &
!!$              / (Vxstr - s1)

   DDstr_x_R  = ( s2 * DD_x_hll - DD_x_hll_flux ) / (s2 - Vxstr)
   DDstr_x_L  = (-s1 * DD_x_hll + DD_x_hll_flux ) / (Vxstr - s1)

   ! Y direction
      
   psistr_y_R = ( s2 * psi_y_hll - Vystr * Ey_y_hll_flux) / (s2 - Vystr) 
   psistr_y_L = (-s1 * psi_y_hll + Vystr * Ey_y_hll_flux) / (Vystr - s1) 
   
   phistr_y_R = ( s2 * phi_y_hll - Vystr * By_y_hll_flux) / (s2 - Vystr) 
   phistr_y_L = (-s1 * phi_y_hll + Vystr * By_y_hll_flux) / (Vystr - s1)

   Systr_y_R  = ( s2 * Sy_y_hll - tau_y_hll_flux * Vystr ) / (s2 - Vystr) 
   Systr_y_L  = (-s1 * Sy_y_hll + tau_y_hll_flux * Vystr ) / (Vystr - s1)
!!$
!!$   Systr_y_R  = (  s2 * Sy_y_hll - Sy_y_hll_flux + pstr_y - Bystr_y_R**2 - Eystr_y_R**2 - ErotB_str_y_y_R * Vystr ) &
!!$              / (s2 - Vystr) 
!!$   Systr_y_L  = ( -s1 * Sy_y_hll + Sy_y_hll_flux - pstr_y + Bystr_y_R**2 + Eystr_y_R**2 + ErotB_str_y_y_R * Vystr ) &
!!$              / (Vystr - s1) 
   
   Szstr_y_R  = ( s2 * Sz_y_hll - Sz_y_hll_flux - ErotB_str_z_y_R * Vystr - Eystr_y_R * Ezstr_y_R - Bystr_y_R * Bzstr_y_R )  &
              / ( s2 - Vystr)
   Szstr_y_L  = (-s1 * Sz_y_hll + Sz_y_hll_flux + ErotB_str_z_y_R * Vystr + Eystr_y_L * Ezstr_y_L + Bystr_y_L * Bzstr_y_L ) &
              / (Vystr - s1)

   Sxstr_y_R  = ( s2 * Sx_y_hll - Sx_y_hll_flux - ErotB_str_x_y_R * Vystr - Eystr_y_R * Exstr_y_R - Bystr_y_R * Bxstr_y_R )  &
              / ( s2 - Vystr)
   Sxstr_y_L  = (-s1 * Sx_y_hll + Sx_y_hll_flux + ErotB_str_x_y_R * Vystr + Eystr_y_L * Exstr_y_L + Bystr_y_L * Bxstr_y_L )  &
              / (Vystr - s1)

   taustr_y_R = (Systr_y_R - ErotB_str_y_y_R) / Vystr - pstr_y + E2B2_str_y 
   taustr_y_L = (Systr_y_L - ErotB_str_y_y_R) / Vystr - pstr_y + E2B2_str_y
!!$
!!$   taustr_y_R = (  s2 * tau_y_hll - tau_y_hll_flux + (pstr_y - E2B2_str_y) * Vystr + ErotB_str_y_y_R ) &
!!$              / (s2 - Vystr)
!!$   taustr_y_L = ( -s1 * tau_y_hll + tau_y_hll_flux - (pstr_y - E2B2_str_y) * Vystr - ErotB_str_y_y_R ) &
!!$              / (Vystr - s1)

   DDstr_y_R  = ( s2 * DD_y_hll - DD_y_hll_flux ) / (s2 - Vystr)
   DDstr_y_L  = (-s1 * DD_y_hll + DD_y_hll_flux ) / (Vystr - s1)

!!$   Vxstr = 0.d0
!!$   Vystr = 0.d0

!   print*, Vxstr, Vystr
 

   ! 5) Finally find the numerical fluxes fot the cases Vxstr=0 and Vxstr =\ 0

   !**************************** STAR FLUXES ******************************************


   ! Using Rankine-Hugoniot conditions and consistency flux conditions, we can calculate
   ! Right and Left star flux values as:

   ! X Direction

       if ( Vxstr == 0.d0 ) then

 ! EM FLUXES

    psistr_x_flux_R =  psi_x_hll_flux 
    psistr_x_flux_L =  psi_x_hll_flux 
    
    phistr_x_flux_R =  phi_x_hll_flux 
    phistr_x_flux_L =  phi_x_hll_flux 

    Exstr_x_flux_R  =  Ex_x_hll_flux  
    Exstr_x_flux_L  =  Ex_x_hll_flux  

    Eystr_x_flux_R  =  Ey_x_hll_flux  
    Eystr_x_flux_L  =  Ey_x_hll_flux  

    Ezstr_x_flux_R  =  Ez_x_hll_flux  
    Ezstr_x_flux_L  =  Ez_x_hll_flux  

    Bxstr_x_flux_R  =  Bx_x_hll_flux  
    Bxstr_x_flux_L  =  Bx_x_hll_flux  

    Bystr_x_flux_R  =  By_x_hll_flux  
    Bystr_x_flux_L  =  By_x_hll_flux  

    Bzstr_x_flux_R  =  Bz_x_hll_flux  
    Bzstr_x_flux_L  =  Bz_x_hll_flux  
    

    !     Density flux (FD = \rho W V)

    FDxstr_R = DD_x_hll_flux  
    FDxstr_L = DD_x_hll_flux  

! these values must be canceled when we made diference of fluxes --->

    FDzstr_R = 1.d0
    FDzstr_L = 1.d0


!     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxstr_L  = q_x_hll_flux  
    Jxstr_R  = q_x_hll_flux  

! these values must be canceled when we made diference of fluxes --->

    Jzstr_R = 1.d0
    Jzstr_L = 1.d0

    !     Flujos conservados de energía

    Ftauxstr_R = tau_x_hll_flux  
    Ftauxstr_L = tau_x_hll_flux  

! these values must be canceled when we made diference of fluxes --->

    Ftauzstr_R = 1.d0
    Ftauzstr_L = 1.d0


    !     Components of the flux momentum tensor 

    FSxxstr_R =  Sx_x_hll_flux  
    FSxxstr_L =  Sx_x_hll_flux  

    FSxystr_R =  Sy_x_hll_flux  
    FSxystr_L =  Sy_x_hll_flux  

    FSxzstr_R =  Sz_x_hll_flux  
    FSxzstr_L =  Sz_x_hll_flux  

! these values must be canceled when we made diference of fluxes --->

    FSzxstr_R = 1.d0
    FSzxstr_L = 1.d0
    FSzystr_R = 1.d0
    FSzystr_L = 1.d0
    FSzzstr_R = 1.d0
    FSzzstr_L = 1.d0


       else

 ! EM FLUXES

    psistr_x_flux_R =  psi_x_hll_flux + (s2 * (Vxstr - s1) * (psistr_x_R - psistr_x_L) )/(s2 -s1) 
    psistr_x_flux_L =  psi_x_hll_flux - (s1 * (s2 - Vxstr) * (psistr_x_R - psistr_x_L) )/(s2 -s1) 
    
    phistr_x_flux_R =  phi_x_hll_flux + (s2 * (Vxstr - s1) * (phistr_x_R - phistr_x_L) )/(s2 -s1) 
    phistr_x_flux_L =  phi_x_hll_flux - (s1 * (s2 - Vxstr) * (phistr_x_R - phistr_x_L) )/(s2 -s1) 

    Exstr_x_flux_R  =  Ex_x_hll_flux  
    Exstr_x_flux_L  =  Ex_x_hll_flux  

    Eystr_x_flux_R  =  Ey_x_hll_flux  
    Eystr_x_flux_L  =  Ey_x_hll_flux  

    Ezstr_x_flux_R  =  Ez_x_hll_flux  
    Ezstr_x_flux_L  =  Ez_x_hll_flux  

    Bxstr_x_flux_R  =  Bx_x_hll_flux  
    Bxstr_x_flux_L  =  Bx_x_hll_flux  

    Bystr_x_flux_R  =  By_x_hll_flux  
    Bystr_x_flux_L  =  By_x_hll_flux  

    Bzstr_x_flux_R  =  Bz_x_hll_flux  
    Bzstr_x_flux_L  =  Bz_x_hll_flux  
    

    !     Density flux (FD = \rho W V)

    FDxstr_R = DD_x_hll_flux  + s2 * (Vxstr - s1) * (DDstr_x_R - DDstr_x_L) / (s2 -s1)
    FDxstr_L = DD_x_hll_flux  - s1 * (s2 - Vxstr) * (DDstr_x_R - DDstr_x_L) / (s2 -s1)

! these values must be canceled when we made diference of fluxes --->

    FDzstr_R = 1.d0
    FDzstr_L = 1.d0


!     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxstr_L  = q_x_hll_flux  
    Jxstr_R  = q_x_hll_flux  

! these values must be canceled when we made diference of fluxes --->

    Jzstr_R = 1.d0
    Jzstr_L = 1.d0

    !     Flujos conservados de energía

    Ftauxstr_R = tau_x_hll_flux  + s2 * (Vxstr - s1) * (taustr_x_R - taustr_x_L) / (s2 -s1)
    Ftauxstr_L = tau_x_hll_flux  - s1 * (s2 - Vxstr) * (taustr_x_R - taustr_x_L) / (s2 -s1)

! these values must be canceled when we made diference of fluxes --->

    Ftauzstr_R = 1.d0
    Ftauzstr_L = 1.d0


    !     Components of the flux momentum tensor 

    FSxxstr_R =  Sx_x_hll_flux  + s2 * (Vxstr - s1) * (Sxstr_x_R - Sxstr_x_L) / (s2 -s1)
    FSxxstr_L =  Sx_x_hll_flux  - s1 * (s2 - Vxstr) * (Sxstr_x_R - Sxstr_x_L) / (s2 -s1)

    FSxystr_R =  Sy_x_hll_flux  + s2 * (Vxstr - s1) * (Systr_x_R - Systr_x_L) / (s2 -s1)
    FSxystr_L =  Sy_x_hll_flux  - s1 * (s2 - Vxstr) * (Systr_x_R - Systr_x_L) / (s2 -s1)

    FSxzstr_R =  Sz_x_hll_flux  + s2 * (Vxstr - s1) * (Szstr_x_R - Szstr_x_L) / (s2 -s1)
    FSxzstr_L =  Sz_x_hll_flux  - s1 * (s2 - Vxstr) * (Szstr_x_R - Szstr_x_L) / (s2 -s1)

! these values must be canceled when we made diference of fluxes --->

    FSzxstr_R = 1.d0
    FSzxstr_L = 1.d0
    FSzystr_R = 1.d0
    FSzystr_L = 1.d0
    FSzzstr_R = 1.d0
    FSzzstr_L = 1.d0

    end if

!////////////////////////////// Y DIRECTION ///////////////////////////////////////////////////

       if ( Vystr == 0.d0 ) then

          
 ! EM FLUXES

    psistr_y_flux_R =  psi_y_hll_flux 
    psistr_y_flux_L =  psi_y_hll_flux 
    
    phistr_y_flux_R =  phi_y_hll_flux 
    phistr_y_flux_L =  phi_y_hll_flux 

    Exstr_y_flux_R  =  Ex_y_hll_flux  
    Exstr_y_flux_L  =  Ex_y_hll_flux  

    Eystr_y_flux_R  =  Ey_y_hll_flux  
    Eystr_y_flux_L  =  Ey_y_hll_flux  

    Ezstr_y_flux_R  =  Ez_y_hll_flux  
    Ezstr_y_flux_L  =  Ez_y_hll_flux  

    Bxstr_y_flux_R  =  Bx_y_hll_flux  
    Bxstr_y_flux_L  =  Bx_y_hll_flux  

    Bystr_y_flux_R  =  By_y_hll_flux  
    Bystr_y_flux_L  =  By_y_hll_flux  

    Bzstr_y_flux_R  =  Bz_y_hll_flux  
    Bzstr_y_flux_L  =  Bz_y_hll_flux  
    

    !     Density flux (FD = \rho W V)

    FDystr_R = DD_y_hll_flux  
    FDystr_L = DD_y_hll_flux  

!     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jystr_L  = q_y_hll_flux  
    Jystr_R  = q_y_hll_flux  

    !     Flujos conservados de energía

    Ftauystr_R = tau_y_hll_flux  
    Ftauystr_L = tau_y_hll_flux  

    !     Components of the flux momentum tensor 

    FSyxstr_R =  Sx_y_hll_flux  
    FSyxstr_L =  Sx_y_hll_flux  

    FSyystr_R =  Sy_y_hll_flux  
    FSyystr_L =  Sy_y_hll_flux  

    FSyzstr_R =  Sz_y_hll_flux  
    FSyzstr_L =  Sz_y_hll_flux  

       else

 ! EM FLUXES

    psistr_y_flux_R =  psi_y_hll_flux + (s2 * (Vystr - s1) * (psistr_y_R - psistr_y_L) )/(s2 -s1) 
    psistr_y_flux_L =  psi_y_hll_flux - (s1 * (s2 - Vystr) * (psistr_y_R - psistr_y_L) )/(s2 -s1) 
    
    phistr_y_flux_R =  phi_y_hll_flux + (s2 * (Vystr - s1) * (phistr_y_R - phistr_y_L) )/(s2 -s1) 
    phistr_y_flux_L =  phi_y_hll_flux - (s1 * (s2 - Vystr) * (phistr_y_R - phistr_y_L) )/(s2 -s1) 

    Exstr_y_flux_R  =  Ex_y_hll_flux  
    Exstr_y_flux_L  =  Ex_y_hll_flux  

    Eystr_y_flux_R  =  Ey_y_hll_flux  
    Eystr_y_flux_L  =  Ey_y_hll_flux  

    Ezstr_y_flux_R  =  Ez_y_hll_flux  
    Ezstr_y_flux_L  =  Ez_y_hll_flux  

    Bxstr_y_flux_R  =  Bx_y_hll_flux  
    Bxstr_y_flux_L  =  Bx_y_hll_flux  

    Bystr_y_flux_R  =  By_y_hll_flux  
    Bystr_y_flux_L  =  By_y_hll_flux  

    Bzstr_y_flux_R  =  Bz_y_hll_flux  
    Bzstr_y_flux_L  =  Bz_y_hll_flux  
    

    !     Density flux (FD = \rho W V)

    FDystr_R = DD_y_hll_flux  + s2 * (Vystr - s1) * (DDstr_y_R - DDstr_y_L) / (s2 -s1)
    FDystr_L = DD_y_hll_flux  - s1 * (s2 - Vystr) * (DDstr_y_R - DDstr_y_L) / (s2 -s1)

    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jystr_L  = q_y_hll_flux 
    Jystr_R  = q_y_hll_flux 

    !     Flujos conservados de energía

    Ftauystr_R = tau_y_hll_flux  + s2 * (Vystr - s1) * (taustr_y_R - taustr_y_L) / (s2 -s1)
    Ftauystr_L = tau_y_hll_flux  - s1 * (s2 - Vystr) * (taustr_y_R - taustr_y_L) / (s2 -s1)



    !     Components of the flux momentum tensor 

    FSyxstr_R =  Sx_y_hll_flux  + s2 * (Vystr - s1) * (Sxstr_y_R - Sxstr_y_L) / (s2 -s1)
    FSyxstr_L =  Sx_y_hll_flux  - s1 * (s2 - Vystr) * (Sxstr_y_R - Sxstr_y_L) / (s2 -s1)

    FSyystr_R =  Sy_y_hll_flux  + s2 * (Vystr - s1) * (Systr_y_R - Systr_y_L) / (s2 -s1)
    FSyystr_L =  Sy_y_hll_flux  - s1 * (s2 - Vystr) * (Systr_y_R - Systr_y_L) / (s2 -s1)

    FSyzstr_R =  Sz_y_hll_flux  + s2 * (Vystr - s1) * (Szstr_y_R - Szstr_y_L) / (s2 -s1)
    FSyzstr_L =  Sz_y_hll_flux  - s1 * (s2 - Vystr) * (Szstr_y_R - Szstr_y_L) / (s2 -s1)


    end if

!******************************************************************************************************************

!/////////////////////// Finally we use the HLLC criterion for flux //////////////////////////////////////////////


    if ( s1 .lt. 0.d0 .and. 0.d0 .lt. Vxstr ) then

! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,j,1,l)  =   psistr_x_flux_L
       Ezlaxpsi(i,j,1,l)  =   1.d0

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,j,1,l)  =   phistr_x_flux_L
       Bzlaxphi(i,j,1,l)  =   1.d0

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,j,1,l)  =   Bxstr_x_flux_L
       EylaxBx (i,j,1,l)  = - 1.d0            ! minus sign

! Flujos By

       EzlaxBy (i,j,1,l)  = - Bystr_x_flux_L ! minus sign
       ExlaxBy (i,j,1,l)  =   1.d0

! Flujos Bz

       EylaxBz (i,j,1,l)  =   Bzstr_x_flux_L
       philaxBz(i,j,1,l)  =   1.d0

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,j,1,l)  =   Exstr_x_flux_L
       BylaxEx(i,j,1,l)   =   1.d0

! Flujos Ey


       BzlaxEy (i,j,1,l)  =   Eystr_x_flux_L
       BxlaxEy (i,j,1,l)  = - 1.d0           ! minus sign

! Flujos Ez


       BylaxEz (i,j,1,l)  = -  Ezstr_x_flux_L ! minus sign
       psilaxEz(i,j,1,l)  =    1.d0


! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,j,1,l)     =   Jxstr_L  
       Jzlax(i,j,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,j,1,l)    =   FDxstr_L  
       FDzlax(i,j,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,j,1,l)  =   Ftauxstr_L 
       Ftauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,j,1,l)  =   psistr_x_L 
       psitauzlax(i,j,1,l)  =   1.d0


       phitauxlax(i,j,1,l)  =   phistr_x_L 
       phitauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,j,1,l)   =   FSxxstr_L 
       FSxylax(i,j,1,l)   =   FSxystr_L 
       FSxzlax(i,j,1,l)   =   FSxzstr_L 

 
       FSzxlax(i,j,1,l)  =   1.d0
       FSzylax(i,j,1,l)  =   1.d0
       FSzzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,j,1,l)   = Bxstr_x_L  
       BzSxlax(i,j,1,l)   = 1.d0



       BxSylax(i,j,1,l)   = Bxstr_x_L  
       BzSylax(i,j,1,l)   = 1.d0


       BxSzlax(i,j,1,l)   = Bxstr_x_L  
       BzSzlax(i,j,1,l)   = 1.d0

! -------------------------------------- EGLM-------------------------------------- 

    else if (Vxstr .lt. 0.d0 .and. 0.d0 .lt. s2) then


! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,j,1,l)  =   psistr_x_flux_R
       Ezlaxpsi(i,j,1,l)  =   1.d0

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,j,1,l)  =   phistr_x_flux_R 
       Bzlaxphi(i,j,1,l)  =   1.d0

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,j,1,l)  =   Bxstr_x_flux_R
       EylaxBx (i,j,1,l)  = - 1.d0            ! minus sign

! Flujos By

       EzlaxBy (i,j,1,l)  = -  Bystr_x_flux_R ! minus sign
       ExlaxBy (i,j,1,l)  =    1.d0

! Flujos Bz

       EylaxBz (i,j,1,l)  =   Bzstr_x_flux_R
       philaxBz(i,j,1,l)  =   1.d0

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,j,1,l)  =   Exstr_x_flux_R
       BylaxEx(i,j,1,l)   =   1.d0

! Flujos Ey


       BzlaxEy (i,j,1,l)  =   Eystr_x_flux_R
       BxlaxEy (i,j,1,l)  = - 1.d0            ! minus sign

! Flujos Ez


       BylaxEz (i,j,1,l)  = -  Ezstr_x_flux_R ! minus sign
       psilaxEz(i,j,1,l)  =    1.d0


! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,j,1,l)     =   Jxstr_R  
       Jzlax(i,j,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,j,1,l)    =   FDxstr_R  
       FDzlax(i,j,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,j,1,l)  =   Ftauxstr_R 
       Ftauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,j,1,l)  =   psistr_x_R 
       psitauzlax(i,j,1,l)  =   1.d0


       phitauxlax(i,j,1,l)  =   phistr_x_R 
       phitauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,j,1,l)   =   FSxxstr_R 
       FSxylax(i,j,1,l)   =   FSxystr_R 
       FSxzlax(i,j,1,l)   =   FSxzstr_R 

       FSzxlax(i,j,1,l)  =   1.d0
       FSzylax(i,j,1,l)  =   1.d0
       FSzzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,j,1,l)   = Bxstr_x_R  
       BzSxlax(i,j,1,l)   = 1.d0



       BxSylax(i,j,1,l)   = Bxstr_x_R  
       BzSylax(i,j,1,l)   = 1.d0


       BxSzlax(i,j,1,l)   = Bxstr_x_R  
       BzSzlax(i,j,1,l)   = 1.d0

! -------------------------------------- EGLM-------------------------------------- 

   else if (Vxstr .eq. 0.d0 ) then


! Electric Gauss HLL flux
!_________________________________________________________________________

       Exlaxpsi(i,j,1,l)  =   0.5d0 * (psistr_x_flux_R + psistr_x_flux_L)
       Ezlaxpsi(i,j,1,l)  =   1.d0

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bxlaxphi(i,j,1,l)  =   0.5d0 * (phistr_x_flux_R + phistr_x_flux_L)
       Bzlaxphi(i,j,1,l)  =   1.d0

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBx(i,j,1,l)  =   0.5d0 * (Bxstr_x_flux_R + Bxstr_x_flux_L)
       EylaxBx (i,j,1,l)  = - 1.d0                                      ! minus sign

! Flujos By

       EzlaxBy (i,j,1,l)  = - 0.5d0 * (Bystr_x_flux_R + Bystr_x_flux_L) ! minus sign
       ExlaxBy (i,j,1,l)  =   1.d0

! Flujos Bz

       EylaxBz (i,j,1,l)  =   0.5d0 * (Bzstr_x_flux_R + Bzstr_x_flux_L)
       philaxBz(i,j,1,l)  =   1.d0

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxEx(i,j,1,l)  =   0.5d0 * (Exstr_x_flux_R + Exstr_x_flux_L)
       BylaxEx(i,j,1,l)   =   1.d0

! Flujos Ey


       BzlaxEy (i,j,1,l)  =   0.5d0 * (Eystr_x_flux_R + Eystr_x_flux_L)
       BxlaxEy (i,j,1,l)  = - 1.d0                                      ! minus sign

! Flujos Ez


       BylaxEz (i,j,1,l)  = - 0.5d0 * (Ezstr_x_flux_R + Ezstr_x_flux_L) ! minus sign
       psilaxEz(i,j,1,l)  =   1.d0

! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,j,1,l)     =   0.5d0 * (Jxstr_R + Jxstr_L)
       Jzlax(i,j,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,j,1,l)    =   0.5d0 * (FDxstr_R + FDxstr_L)
       FDzlax(i,j,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,j,1,l)  =   0.5d0 * (Ftauxstr_R + Ftauxstr_L)
       Ftauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,j,1,l)  =   0.5d0 * (psistr_x_R + psistr_x_L)
       psitauzlax(i,j,1,l)  =   1.d0


       phitauxlax(i,j,1,l)  =   0.5d0 * (phistr_x_R + phistr_x_L)
       phitauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,j,1,l)   =   0.5d0 * (FSxxstr_R + FSxxstr_L)
       FSxylax(i,j,1,l)   =   0.5d0 * (FSxystr_R + FSxystr_L)
       FSxzlax(i,j,1,l)   =   0.5d0 * (FSxzstr_R + FSxzstr_L)

 
       FSzxlax(i,j,1,l)  =   1.d0
       FSzylax(i,j,1,l)  =   1.d0
       FSzzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,j,1,l)   = 0.5d0 * (Bxstr_x_R  + Bxstr_x_L)
       BzSxlax(i,j,1,l)   = 1.d0



       BxSylax(i,j,1,l)   = 0.5d0 * (Bxstr_x_R  + Bxstr_x_L)
       BzSylax(i,j,1,l)   = 1.d0


       BxSzlax(i,j,1,l)   = 0.5d0 * (Bxstr_x_R  + Bxstr_x_L)
       BzSzlax(i,j,1,l)   = 1.d0

! -------------------------------------- EGLM-------------------------------------- 

    else 

       print*, "Revise Vxstar velocity it could be imaginary subroutine hllc_flow "
       stop

    end if


   !******************************************************************************************************************

!/////////////////////// Finally we use the HLLC criterion for flux //////////////////////////////////////////////


    if ( s1 .lt. 0.d0 .and. 0.d0 .lt. Vystr ) then

! Electric Gauss HLL flux
!_________________________________________________________________________

       Eylaxpsi(i,j,1,l)  =   psistr_y_flux_L

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bylaxphi(i,j,1,l)  =   phistr_y_flux_L

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       EzlaxBx (i,j,1,l)  =   Bxstr_y_flux_L

! Flujos By

       philaxBy(i,j,1,l)  =   Bystr_y_flux_L

! Flujos Bz

       ExlaxBz (i,j,1,l)  = - Bzstr_y_flux_L ! minus sign

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       BzlaxEx(i,j,1,l)   = - Exstr_y_flux_L ! minus sign

! Flujos Ey

       psilaxEy(i,j,1,l)  =   Eystr_y_flux_L

! Flujos Ez

       BxlaxEz (i,j,1,l)  =   Ezstr_y_flux_L

! Conserved current HLL flux
!_________________________________________________________________________


       Jylax(i,j,1,l)     =   Jystr_L 

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDylax(i,j,1,l)    =   FDystr_L

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauylax(i,j,1,l)  =   Ftauystr_L 

! -------------------------------------- EGLM-------------------------------------- 

       psitauylax(i,j,1,l)  =   psistr_y_L 

       phitauylax(i,j,1,l)  =   phistr_y_L 

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSyxlax(i,j,1,l)   =   FSyxstr_L
       FSyylax(i,j,1,l)   =   FSyystr_L
       FSyzlax(i,j,1,l)   =   FSyzstr_L


! -------------------------------------- EGLM-------------------------------------- 


       BySxlax(i,j,1,l)   = Bxstr_y_L  

       BySylax(i,j,1,l)   = Bxstr_y_L 

       BySzlax(i,j,1,l)   = Bxstr_y_L  

! -------------------------------------- EGLM-------------------------------------- 

    else if (Vystr .lt. 0.d0 .and. 0.d0 .lt. s2) then


! Electric Gauss HLL flux
!_________________________________________________________________________

       Eylaxpsi(i,j,1,l)  =   psistr_y_flux_R

! Magnetic Gauss HLL flux
!_________________________________________________________________________


       Bylaxphi(i,j,1,l)  =   phistr_y_flux_R 

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       EzlaxBx (i,j,1,l)  =   Bxstr_y_flux_R

! Flujos By

       philaxBy(i,j,1,l)  =    Bystr_y_flux_R 

! Flujos Bz

       ExlaxBz (i,j,1,l)  = - Bzstr_y_flux_R ! minus sign

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       BzlaxEx(i,j,1,l)   = - Exstr_y_flux_R ! minus sign

! Flujos Ey

       psilaxEy(i,j,1,l)  =   Eystr_y_flux_R

! Flujos Ez

       BxlaxEz (i,j,1,l)  =    Ezstr_y_flux_R

! Conserved current HLL flux
!_________________________________________________________________________

       Jylax(i,j,1,l)     =   Jystr_R 

! Conserved Mass HLL flux
!_________________________________________________________________________

       FDylax(i,j,1,l)    =   FDystr_R  

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauylax(i,j,1,l)  =   Ftauystr_R 

! -------------------------------------- EGLM-------------------------------------- 

       psitauylax(i,j,1,l)  =   psistr_y_R 

       phitauylax(i,j,1,l)  =   phistr_y_R 

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSyxlax(i,j,1,l)   =   FSyxstr_R 
       FSyylax(i,j,1,l)   =   FSyystr_R 
       FSyzlax(i,j,1,l)   =   FSyzstr_R 

! -------------------------------------- EGLM-------------------------------------- 


       BySxlax(i,j,1,l)   = Bystr_y_R

       BySylax(i,j,1,l)   = Bystr_y_R  

       BySzlax(i,j,1,l)   = Bystr_y_R  

! -------------------------------------- EGLM-------------------------------------- 

   else if (Vystr .eq. 0.d0 ) then


! Electric Gauss HLL flux
!_________________________________________________________________________

       Eylaxpsi(i,j,1,l)  =   0.5d0 * (psistr_y_flux_R + psistr_y_flux_L)

! Magnetic Gauss HLL flux
!_________________________________________________________________________

       Bylaxphi(i,j,1,l)  =   0.5d0 * (phistr_y_flux_R + phistr_y_flux_L)

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       EzlaxBx (i,j,1,l)  =   0.5d0 * (Bxstr_y_flux_R + Bxstr_y_flux_L)

! Flujos By

       philaxBy(i,j,1,l)  =    0.5d0 * (Bystr_y_flux_R + Bystr_y_flux_L)

! Flujos Bz

       ExlaxBz (i,j,1,l)  = -  0.5d0 * (Bzstr_y_flux_R + Bzstr_y_flux_L) ! minus sign

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       BzlaxEx(i,j,1,l)   = - 0.5d0 * (Exstr_y_flux_R + Exstr_y_flux_L) ! minus sign

! Flujos Ey


       psilaxEy(i,j,1,l)  =   0.5d0 * (Eystr_y_flux_R + Eystr_y_flux_L)

! Flujos Ez


       BxlaxEz (i,j,1,l)  =   0.5d0 * (Ezstr_y_flux_R + Ezstr_y_flux_L)

! Conserved current HLL flux
!_________________________________________________________________________


       Jylax(i,j,1,l)     =   0.5d0 * (Jystr_R + Jystr_L)

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDylax(i,j,1,l)    =   0.5d0 * (FDystr_R + FDystr_L)

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauylax(i,j,1,l)  =   0.5d0 * (Ftauystr_R + Ftauystr_L)

! -------------------------------------- EGLM-------------------------------------- 

       psitauylax(i,j,1,l)  =   0.5d0 * (psistr_y_R + psistr_y_L)

       phitauylax(i,j,1,l)  =   0.5d0 * (phistr_y_R + phistr_y_L)

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSyxlax(i,j,1,l)   =   0.5d0 * (FSyxstr_R + FSyxstr_L)
       FSyylax(i,j,1,l)   =   0.5d0 * (FSyystr_R + FSyystr_L)
       FSyzlax(i,j,1,l)   =   0.5d0 * (FSyzstr_R + FSyzstr_L)

! -------------------------------------- EGLM-------------------------------------- 


       BySxlax(i,j,1,l)   = 0.5d0 * (Bystr_y_R  + Bystr_y_L)

       BySylax(i,j,1,l)   = 0.5d0 * (Bystr_y_R  + Bystr_y_L)

       BySzlax(i,j,1,l)   = 0.5d0 * (Bystr_y_R  + Bystr_y_L)

! -------------------------------------- EGLM-------------------------------------- 


    else 

       print*, "Revise Vystar velocity it could be imaginary subroutine hllc_flow "
       stop

    end if



        end do
     end do


   else

      write(*,*) "STOP: subroutine hllflux"
      write(*,*) "This dimension is not implemented yet"
      stop

   end if


   CONTAINS



         doubleprecision function get_hllc_vel(a,b,c,rt)

        
!!$/*@@
!!$   @function  get_hllc_vel
!!$   @date      Tuesday 18 July 2006
!!$   @author    Miguel Angel Aloy Toras
!!$   @desc
!!$   This function obtains the root of the quadratic equation
!!$
!!$       a*v^2 - b*v + c = 0
!!$
!!$   in the physical range [-1,1] usign a Newton-Rapson iterative method.
!!$  @enddesc
!!$@@*/

        
      implicit none

!     Coeffcients of the quadratic equation a*v^2 - b*v + c = 0
      DOUBLEPRECISION a, b, c, rt

      DOUBLEPRECISION v

      INTEGER, PARAMETER :: maxiter = 100
      INTEGER               iter

      DOUBLEPRECISION, PARAMETER:: minerr = 1.d-12
      DOUBLEPRECISION, PARAMETER:: diferr = 1.d-3
      DOUBLEPRECISION fv, dfdv, oldv, err
     
!*
!*     Initial seed (as if a=0) but limited to the range [-1,1]
!*
      oldv = max(-1.d0, min( 1.d0, rt ))
!c      print*,'Initial seed', oldv, a,b,c
      iter = 0
!*
!*     Start Newton-Rapson loop
!*

333   fv = a*v*v - b*v + c
      dfdv = 2.d0*a*v - b

      v = max(-1.d0, min( 1.d0, oldv - fv/dfdv ))
      err = abs( oldv - v )
!c      print'(i3,1x,5e12.5)',iter,oldv,v,fv,dfdv,err
      oldv = v
      iter = iter+1
!*
!*     Maximum number of iterations reached?
!*
      if (iter.ge.maxiter) then
         print'(a,i3,1x,a)',                                            &
              '[GET_HLLC_VEL]: Non convergence in velocity after',iter, &
              'iterations'
         stop
      endif
!*
!*     Error larger than tolerance?
!*
      if (err.gt.minerr) goto 333

      get_hllc_vel = max(-1.d0 + diferr, min( 1.d0 - diferr, v )) ! v
!*
!*
!*
      return
!*
!*     end of get_hllc_vel
!*
    end function get_hllc_vel   



  end subroutine hllc_flow
