  !     *******************************************************************
  !     Subroutine which calculed  HLL flux flows 
  !     *******************************************************************

   subroutine hlle_flow_mpx

    use scalar
    use parameters
    use threevectors
    use fourvectors
    use funciones

    implicit none


    if (DIM == 1) then

       call boundary_conserved
       call boundary_electric
       call boundary_flux_conserved
       


         do i=-6,imax+6


            
            u(i, 1,1) =  psiint(i,1,1,l)
            u(i, 2,1) =  phiint(i,1,1,l)

            u(i, 3,1) =  Exint (i,1,1,l)
            u(i, 4,1) =  Eyint (i,1,1,l)
            u(i, 5,1) =  Ezint (i,1,1,l)

            u(i, 6,1) =- Exint (i,1,1,l)
            u(i, 7,1) =- Eyint (i,1,1,l)
            u(i, 8,1) =- Ezint (i,1,1,l)

            u(i, 9,1) =  Bxint (i,1,1,l)
            u(i,10,1) =  Byint (i,1,1,l)
            u(i,11,1) =  Bzint (i,1,1,l)

            u(i,12,1) =- Bxint (i,1,1,l)
            u(i,13,1) =- Byint (i,1,1,l)
            u(i,14,1) =- Bzint (i,1,1,l)

            u(i,15,1) =  qint  (i,1,1,l)
            u(i,16,1) =  DDint (i,1,1,l)

            u(i,17,1) =  tauint(i,1,1,l)

            u(i,18,1) =  Sxint (i,1,1,l)
            u(i,19,1) =  Syint (i,1,1,l)
            u(i,20,1) =  Szint (i,1,1,l)

            if (REC_PRIM == 0) then

               u(i,21,1) = Vx(i,1,1)
               u(i,22,1) = Vy(i,1,1)
               u(i,23,1) = Vz(i,1,1)

               u(i,24,1) = p  (i,1,1)
               u(i,25,1) = rho(i,1,1)

            else if (REC_PRIM == 1) then

               u(i,21,1) = Vxint(i,1,1,l)
               u(i,22,1) = Vyint(i,1,1,l)
               u(i,23,1) = Vzint(i,1,1,l)


!!$               u(i,24,1) = p  (i,1,1)
!!$               u(i,25,1) = rho(i,1,1)
!!$               
               u(i,24,1) = pint  (i,1,1,l)
               u(i,25,1) = rhoint(i,1,1,l)
               

            else

               write(*,*) "STOP: subroutine hllflux"
               write(*,*) "REC_PRIM parameter is not valid"
               stop

            end if

!////////////////////////////  FLUXES    ///////////////////////////////
!*********** Condicional para reconstruir flujos usando mp5 o valores a left y right 

!!$            u(i,26,1) = Jx(i,1,1) 
!!$            u(i,27,1) = Jy(i,1,1) 
!!$            u(i,28,1) = Jz(i,1,1) 
!!$
            u(i,26,1) = Jxint(i,1,1,l)
            u(i,27,1) = Jyint(i,1,1,l)
            u(i,28,1) = Jzint(i,1,1,l)


            if ( flux_rec_mp5 == 1 ) then

               u(i,29,1) = FDxint(i,1,1,l)
               u(i,30,1) = FDyint(i,1,1,l)
               u(i,31,1) = FDzint(i,1,1,l)

               u(i,32,1) = Ftauxint(i,1,1,l)
               u(i,33,1) = Ftauyint(i,1,1,l)
               u(i,34,1) = Ftauzint(i,1,1,l)

               u(i,35,1) = FSxxint(i,1,1,l)
               u(i,36,1) = FSxyint(i,1,1,l)
               u(i,37,1) = FSxzint(i,1,1,l)

               u(i,38,1) = FSyxint(i,1,1,l)
               u(i,39,1) = FSyyint(i,1,1,l)
               u(i,40,1) = FSyzint(i,1,1,l)

               u(i,41,1) = FSzxint(i,1,1,l)
               u(i,42,1) = FSzyint(i,1,1,l)
               u(i,43,1) = FSzzint(i,1,1,l)

            end if
            
 
               
         end do


         
         coordenate = 1
         imp = 1

         if (rec_mp == 0) then
            
            call MP5(imp)

         else if (rec_mp == 1) then

            call MP7(imp)

         else if (rec_mp == 2) then

            call MP9(imp)
            
         else

            print*, " rec_mp parameter not valid! "
            stop

         end if
         
         varname = "not__"


         
         do i=-1,imax+1

              

            psiint_L  = up(i  , 1,1) 
            psiint_R  = um(i+1, 1,1)

            phiint_L  = up(i  , 2,1) 
            phiint_R  = um(i+1, 2,1)

            Exint_L   = up(i  , 3,1) 
            Exint_R   = um(i+1, 3,1)
            Eyint_L   = up(i  , 4,1) 
            Eyint_R   = um(i+1, 4,1)
            Ezint_L   = up(i  , 5,1) 
            Ezint_R   = um(i+1, 5,1)

            Exint_L_m = up(i  , 6,1) 
            Exint_R_m = um(i+1, 6,1)
            Eyint_L_m = up(i  , 7,1) 
            Eyint_R_m = um(i+1, 7,1)
            Ezint_L_m = up(i  , 8,1) 
            Ezint_R_m = um(i+1, 8,1)

            Bxint_L   = up(i  , 9,1) 
            Bxint_R   = um(i+1, 9,1)
            Byint_L   = up(i  ,10,1) 
            Byint_R   = um(i+1,10,1)
            Bzint_L   = up(i  ,11,1) 
            Bzint_R   = um(i+1,11,1)

            Bxint_L_m = up(i  ,12,1) 
            Bxint_R_m = um(i+1,12,1)
            Byint_L_m = up(i  ,13,1) 
            Byint_R_m = um(i+1,13,1)
            Bzint_L_m = up(i  ,14,1) 
            Bzint_R_m = um(i+1,14,1)

            qint_L    = up(i  ,15,1) 
            qint_R    = um(i+1,15,1)
            
            DDint_L   = up(i  ,16,1) 
            DDint_R   = um(i+1,16,1)

            tauint_L  = up(i  ,17,1) 
            tauint_R  = um(i+1,17,1)

            Sxint_L   = up(i  ,18,1) 
            Sxint_R   = um(i+1,18,1)
            Syint_L   = up(i  ,19,1) 
            Syint_R   = um(i+1,19,1)
            Szint_L   = up(i  ,20,1) 
            Szint_R   = um(i+1,20,1)

            Vx_L      = up(i  ,21,1) 
            Vx_R      = um(i+1,21,1)
            Vy_L      = up(i  ,22,1) 
            Vy_R      = um(i+1,22,1)
            Vz_L      = up(i  ,23,1) 
            Vz_R      = um(i+1,23,1)

            p_L       = up(i  ,24,1) 
            p_R       = um(i+1,24,1)

            rho_L     = up(i  ,25,1) 
            rho_R     = um(i+1,25,1)

            Jxint_L   = up(i  ,26,1) 
            Jxint_R   = um(i+1,26,1)
            Jyint_L   = up(i  ,27,1) 
            Jyint_R   = um(i+1,27,1)
            Jzint_L   = up(i  ,28,1) 
            Jzint_R   = um(i+1,28,1)


!!$            Jxint_R   = Jxint_L
!!$              
!!$            Jyint_R   = Jyint_L
!!$             
!!$            Jzint_R   = Jzint_L


!     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV


    E2B2int_L = Exint_L**2 + Eyint_L**2 + Ezint_L**2 + &
                Bxint_L**2 + Byint_L**2 + Bzint_L**2

    E2B2int_R = Exint_R**2 + Eyint_R**2 + Ezint_R**2 + &
                Bxint_R**2 + Byint_R**2 + Bzint_R**2
    
    Edotv_L = Exint_L * Vx_L + Eyint_L * Vy_L + Ezint_L * Vz_L
    Edotv_R = Exint_R * Vx_R + Eyint_R * Vy_R + Ezint_R * Vz_R

    vrotB_x_L = (Bzint_L * Vy_L - Byint_L * Vz_L)
    vrotB_y_L = (Bxint_L * Vz_L - Bzint_L * Vx_L)
    vrotB_z_L = (Byint_L * Vx_L - Bxint_L * Vy_L)

    vrotB_x_R = (Bzint_R * Vy_R - Byint_R * Vz_R)
    vrotB_y_R = (Bxint_R * Vz_R - Bzint_R * Vx_R)
    vrotB_z_R = (Byint_R * Vx_R - Bxint_R * Vy_R)

    W_L       = 1.d0/sqrt(1.d0 - (Vx_L**2+Vy_L**2+Vz_L**2))
    W_R       = 1.d0/sqrt(1.d0 - (Vx_R**2+Vy_R**2+Vz_R**2))



    if (source_rec_mpx == 1) then

    Exsource_L  = sigma * W_L * ( Exint_L +  vrotB_x_L  - Edotv_L * Vx_L )
    Exsource_R  = sigma * W_R * ( Exint_R +  vrotB_x_R  - Edotv_R * Vx_R )

    Eysource_L = sigma * W_L  * ( Eyint_L +  vrotB_y_L  - Edotv_L * Vy_L )
    Eysource_R = sigma * W_R  * ( Eyint_R +  vrotB_y_R  - Edotv_R * Vy_R )

    Ezsource_L = sigma * W_L  * ( Ezint_L +  vrotB_z_L  - Edotv_L * Vz_L )
    Ezsource_R = sigma * W_R  * ( Ezint_R +  vrotB_z_R  - Edotv_R * Vz_R )


!!$    Exsource_L = Jxint_L - qint_L  *  Vx_L
!!$    Exsource_R = Jxint_R - qint_R  *  Vx_R
!!$
!!$    Eysource_L = Jyint_L - qint_L  *  Vy_L
!!$    Eysource_R = Jyint_R - qint_R  *  Vy_R
!!$
!!$    Ezsource_L = Jzint_L - qint_L  *  Vz_L
!!$    Ezsource_R = Jzint_R - qint_R  *  Vz_R
    

    end if

            

!!$    Jxint_L  = sigma * W_L * ( Exint_L +  vrotB_x_L  - Edotv_L * Vx_L )   &
!!$                           +   qint_L  *  Vx_L
!!$    Jxint_R  = sigma * W_R * ( Exint_R +  vrotB_x_R  - Edotv_R * Vx_R )   &
!!$                           +   qint_R  *  Vx_R
!!$
!!$    Jyint_L = sigma * W_L  * ( Eyint_L +  vrotB_y_L  - Edotv_L * Vy_L )   &
!!$                           +   qint_L  *  Vy_L
!!$    Jyint_R = sigma * W_R  * ( Eyint_R +  vrotB_y_R  - Edotv_R * Vy_R )   &
!!$                           +   qint_R  *  Vy_R
!!$
!!$    Jzint_L = sigma * W_L  * ( Ezint_L +  vrotB_z_L  - Edotv_L * Vz_L )   &
!!$                           +   qint_L  *  Vz_L
!!$    Jzint_R = sigma * W_R  * ( Ezint_R +  vrotB_z_R  - Edotv_R * Vz_R )   &
!!$                           +   qint_R  *  Vz_R

!!$    Jxint_R   = Jxint_L
!!$
!!$    Jyint_R   = Jyint_L
!!$
!!$    Jzint_R   = Jzint_L

    

!*********** Condicional para reconstruir flujos usando mp5 o valores a left y right
            
            if ( flux_rec_mp5 == 1 ) then
               
               FDxint_L   = up(i  ,29,1) 
               FDxint_R   = um(i+1,29,1)
               FDyint_L   = up(i  ,30,1) 
               FDyint_R   = um(i+1,30,1)
               FDzint_L   = up(i  ,31,1) 
               FDzint_R   = um(i+1,31,1)

               Ftauxint_L = up(i  ,32,1) 
               Ftauxint_R = um(i+1,32,1)
               Ftauyint_L = up(i  ,33,1) 
               Ftauyint_R = um(i+1,33,1)
               Ftauzint_L = up(i  ,34,1) 
               Ftauzint_R = um(i+1,34,1)

               FSxxint_L  = up(i  ,35,1) 
               FSxxint_R  = um(i+1,35,1)
               FSxyint_L  = up(i  ,36,1) 
               FSxyint_R  = um(i+1,36,1)
               FSxzint_L  = up(i  ,37,1) 
               FSxzint_R  = um(i+1,37,1)

               FSyxint_L  = up(i  ,38,1) 
               FSyxint_R  = um(i+1,38,1)
               FSyyint_L  = up(i  ,39,1) 
               FSyyint_R  = um(i+1,39,1)
               FSyzint_L  = up(i  ,40,1) 
               FSyzint_R  = um(i+1,40,1)

               FSzxint_L  = up(i  ,41,1) 
               FSzxint_R  = um(i+1,41,1)
               FSzyint_L  = up(i  ,42,1) 
               FSzyint_R  = um(i+1,42,1)
               FSzzint_L  = up(i  ,43,1) 
               FSzzint_R  = um(i+1,43,1)

               
            else if (flux_rec_mp5 == 0) then
               

       !---------------------------------
       ! Left - Right Fluxes
       !---------------------------------

    !     Electric and Magnetic field modules 


    E2B2int_L = Exint_L**2 + Eyint_L**2 + Ezint_L**2 + &
                Bxint_L**2 + Byint_L**2 + Bzint_L**2

    E2B2int_R = Exint_R**2 + Eyint_R**2 + Ezint_R**2 + &
                Bxint_R**2 + Byint_R**2 + Bzint_R**2
    
    Edotv_L = Exint_L * Vx_L + Eyint_L * Vy_L + Ezint_L * Vz_L
    Edotv_R = Exint_R * Vx_R + Eyint_R * Vy_R + Ezint_R * Vz_R

    vrotB_x_L = (Bzint_L * Vy_L - Byint_L * Vz_L)
    vrotB_y_L = (Bxint_L * Vz_L - Bzint_L * Vx_L)
    vrotB_z_L = (Byint_L * Vx_L - Bxint_L * Vy_L)

    vrotB_x_R = (Bzint_R * Vy_R - Byint_R * Vz_R)
    vrotB_y_R = (Bxint_R * Vz_R - Bzint_R * Vx_R)
    vrotB_z_R = (Byint_R * Vx_R - Bxint_R * Vy_R)

    W_L       = 1.d0/sqrt(1.d0 - (Vx_L**2+Vy_L**2+Vz_L**2))
    W_R       = 1.d0/sqrt(1.d0 - (Vx_R**2+Vy_R**2+Vz_R**2))
    
    
    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV
    

!!$    Jxint_L = sigma * W_L * ( Exint_L +  vrotB_x_L  - Edotv_L * Vx_L )   &
!!$                             +   qint_L  *  Vx_L
!!$    Jxint_R = sigma * W_R * ( Exint_R +  vrotB_x_R  - Edotv_R * Vx_R )   &
!!$                             +   qint_R  *  Vx_R
!!$
!!$    Jyint_L = sigma * W_L  * ( Eyint_L +  vrotB_y_L  - Edotv_L * Vy_L )   &
!!$                             +   qint_L  *  Vy_L
!!$    Jyint_R = sigma * W_R  * ( Eyint_R +  vrotB_y_R  - Edotv_R * Vy_R )   &
!!$                             +   qint_R  *  Vy_R
!!$
!!$    Jzint_L = sigma * W_L  * ( Ezint_L +  vrotB_z_L  - Edotv_L * Vz_L )   &
!!$                             +   qint_L  *  Vz_L
!!$    Jzint_R = sigma * W_R  * ( Ezint_R +  vrotB_z_R  - Edotv_R * Vz_R )   &
!!$                             +   qint_R  *  Vz_R


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


    else 

           print*, " flux_rec_mp5 parameter not valid in subroutine hllc_flow "
           stop

     end if


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

    ! since for HLLE Riemann solver the values of s1 and s2 change, we must force to its values of -1 and 1 for EM fields
    !-----------------------------------------------------------------------------------------------------
    s1 = -1.d0
    s2 =  1.d0
    !-----------------------------------------------------------------------------------------------------

    psi_hll = ( s2 * psiint_R - s1 * psiint_L + Exint_L  - Exint_R   ) / (s2 -s1)
    phi_hll = ( s2 * phiint_R - s1 * phiint_L + Bxint_L  - Bxint_R   ) / (s2 -s1)

    !     hll  values of electric field

    Ex_hll  = ( s2 * Exint_R  - s1 * Exint_L  + psiint_L  - psiint_R  ) / (s2 -s1)
    Ey_hll  = ( s2 * Eyint_R  - s1 * Eyint_L  + Bzint_L   - Bzint_R   ) / (s2 -s1)
    Ez_hll  = ( s2 * Ezint_R  - s1 * Ezint_L  + Byint_L_m - Byint_R_m ) / (s2 -s1)  !ojo minus sign ¡¡¡


    !     hll  values of magnetic field

    Bx_hll  = ( s2 * Bxint_R  - s1 * Bxint_L  + phiint_L  - phiint_R  ) / (s2 -s1)
    By_hll  = ( s2 * Byint_R  - s1 * Byint_L  + Ezint_L_m - Ezint_R_m ) / (s2 -s1)  !ojo minus sign ¡¡¡
    Bz_hll  = ( s2 * Bzint_R  - s1 * Bzint_L  + Eyint_L   - Eyint_R   ) / (s2 -s1)


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

       By_hll_flux  =   0.5d0 * (Ezint_L_m   + Ezint_R_m - ( Byint_R - Byint_L ) ) ! minus sign

! Flujos Bz

       Bz_hll_flux  =   0.5d0 * (Eyint_L    + Eyint_R    - ( Bzint_R - Bzint_L ) )

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       Ex_hll_flux  =   0.5d0 * (psiint_L   + psiint_R   - ( Exint_R - Exint_L ) )

! Flujos Ey


       Ey_hll_flux  =   0.5d0 * (Bzint_L    + Bzint_R    - ( Eyint_R - Eyint_L ) )

! Flujos Ez


       Ez_hll_flux  =   0.5d0 * (Byint_L_m  + Byint_R_m  - ( Ezint_R - Ezint_L ) ) ! minus sign

!=======================================================================================================================   
! characteristic sound velocity, with cs= sqrt(gamma * p/(rho * h), observe that you use h = rho * h
!=======================================================================================================================      

   Cs_R      = sqrt(gamma * p_R / (enthpy_R))
   Cs_L      = sqrt(gamma * p_L / (enthpy_L))

   V2_R      =  Vx_R**2 +  Vy_R**2 + Vz_R**2
   V2_L      =  Vx_L**2 +  Vy_L**2 + Vz_L**2

   lambda_H_p_R = ( Vx_R * ( 1.d0 - Cs_R**2) + Cs_R * sqrt( (1.d0 - V2_R) * (1.d0 - V2_R * Cs_R**2 - Vx_R * (1.d0 - Cs_R**2) )) ) &
                / ( 1.d0 - Cs_R**2 * V2_R)

   lambda_H_p_L = ( Vx_L * ( 1.d0 - Cs_L**2) + Cs_L * sqrt( (1.d0 - V2_L) * (1.d0 - V2_L * Cs_L**2 - Vx_L * (1.d0 - Cs_L**2) )) ) &
                / ( 1.d0 - Cs_L**2 * V2_L)

   
   lambda_H_m_R = ( Vx_R * ( 1.d0 - Cs_R**2) - Cs_R * sqrt( (1.d0 - V2_R) * (1.d0 - V2_R * Cs_R**2 - Vx_R * (1.d0 - Cs_R**2) )) ) &
                / ( 1.d0 - Cs_R**2 * V2_R)

   lambda_H_m_L = ( Vx_L * ( 1.d0 - Cs_L**2) - Cs_L * sqrt( (1.d0 - V2_L) * (1.d0 - V2_L * Cs_L**2 - Vx_L * (1.d0 - Cs_L**2) )) ) &
                / ( 1.d0 - Cs_L**2 * V2_L)


   !   These estimates are not recommended for practical computations, Eleuterio Toro Sec.10.5 Eq.10.48 3ed (2009)

   !   If appear NaN numbers. By definition, NAN is not equal to anything, even itself. Simply compare the variable to itself: if(my_var /= my_var)


!!$   eps2 = 8.d-2
!!$
!!$
!!$   if ( isnan(lambda_H_m_L) .or. isnan(lambda_H_m_R) .or. &
!!$        isnan(lambda_H_p_L) .or. isnan(lambda_H_p_R)    ) then
!!$
!!$      if ( isnan(lambda_H_m_L) ) lambda_H_m_L = -1.d0
!!$      if ( isnan(lambda_H_m_R) ) lambda_H_m_R = -1.d0
!!$
!!$      if ( isnan(lambda_H_p_L) ) lambda_H_p_L =  1.d0
!!$      if ( isnan(lambda_H_p_R) ) lambda_H_p_R =  1.d0
!!$
!!$      s1 = min(lambda_H_m_L, lambda_H_m_R)
!!$      s2 = max(lambda_H_p_L, lambda_H_p_R)
!!$
!!$      !  proposed average eigenvalues for the left and right non–linear waves, E. Toro (2009) 3ed
!!$      ! we made avarage over cuadri-vectors ux = 0.5 (Wl vl + ur vl)
!!$      ! find the lorentz factor -W^2 + ux^2 = -1 ---> W^2 = ux^2 + 1
!!$      ! and the three velocity as vx = sqrt( 1 - 1/W^2)
!!$
!!$   else if ( abs(lambda_H_m_L) .le. 1.d0 - eps2 .and. abs(lambda_H_m_R) .le. 1.d0 - eps2 .and. &
!!$             abs(lambda_H_p_L) .le. 1.d0 - eps2 .and. abs(lambda_H_p_R) .le. 1.d0 - eps2     ) then
!!$
!!$      W_H_p_L = 1.d0 / sqrt(1.d0 - lambda_H_p_L**2)
!!$      W_H_p_R = 1.d0 / sqrt(1.d0 - lambda_H_p_R**2)
!!$
!!$      W_H_m_L = 1.d0 / sqrt(1.d0 - lambda_H_m_L**2)
!!$      W_H_m_R = 1.d0 / sqrt(1.d0 - lambda_H_m_R**2)
!!$
!!$      ux_p_m    = 0.5d0 * (W_H_p_L * lambda_H_p_L + W_H_p_R * lambda_H_p_R)
!!$      ux_m_m    = 0.5d0 * (W_H_m_L * lambda_H_m_L + W_H_m_R * lambda_H_m_R)
!!$
!!$      W2m_p     = ux_p_m**2  + 1.d0
!!$      W2m_m     = ux_m_m**2  + 1.d0
!!$
!!$
!!$      s2       = sqrt( 1.d0 - 1.d0/W2m_p)
!!$      s1       = sqrt( 1.d0 - 1.d0/W2m_m)
!!$
!!$   else

   s1 = min(lambda_H_m_L, lambda_H_m_R)
   s2 = max(lambda_H_p_L, lambda_H_p_R)

!!$   end if


   if ( isnan(s1) .or. isnan(s2)) then

      print*, " "
      print*, "s1 ==", s1, l, i, h, lambda_H_m_L, lambda_H_m_R, W_H_m_L, W_H_m_R
      print*, "s2 ==", s1, l, i, h, lambda_H_p_L, lambda_H_p_R, W_H_p_L, W_H_p_R
      print*, "NaN value for characteristic velocities s1 and s1. HLLC soubroutine"
      print*, " "
      stop

   end if

   eps2 = 1.d-9

   
   if ( abs(s1) .gt. 1.d0 ) then

      s1 = max(-1.d0 + eps2 , min( 1.d0 - eps2, s1 ))
      
   end if

   if ( abs(s2) .gt. 1.d0 ) then

      s2 = max(-1.d0 + eps2 , min( 1.d0 - eps2, s2 ))
      
   end if
    
!=======================================================================================================================

   !/////////////////////// HLL Variables  ////////////////////////////////////

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

! Conserved current HLL flux
!_________________________________________________________________________


       q_hll_flux   =   (s2 * Jxint_L  - s1 * Jxint_R + s1 * s2 * ( qint_R - qint_L ) ) / (s2 -s1)

! Conserved Mass HLL flux
!_________________________________________________________________________


       DD_hll_flux  =   (s2 * FDxint_L  - s1 * FDxint_R + s1 * s2 * (DDint_R - DDint_L) ) / (s2 -s1)

! Conserved Energy HLL flux
!_________________________________________________________________________


       tau_hll_flux =   (s2 * Ftauxint_L - s1 *  Ftauxint_R + s1 * s2 * (tauint_R - tauint_L) ) / (s2 -s1)


! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       Sx_hll_flux  =   (s2 * FSxxint_L - s1 * FSxxint_R + s1 * s2 * (Sxint_R - Sxint_L) ) /  (s2 -s1)
       Sy_hll_flux  =   (s2 * FSxyint_L - s1 * FSxyint_R + s1 * s2 * (Syint_R - Syint_L) ) /  (s2 -s1)
       Sz_hll_flux  =   (s2 * FSxzint_L - s1 * FSxzint_R + s1 * s2 * (Szint_R - Szint_L) ) /  (s2 -s1)



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

   psistr_R =  psi_hll
   psistr_L =  psi_hll 

   phistr_R =  phi_hll
   phistr_L =  phi_hll



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

      

!=======================================================================================================================   
! characteristic sound velocity, with cs= sqrt(gamma * p/(rho * h), observe that you use h = rho * h
!=======================================================================================================================      

   Cs_R      = sqrt(gamma * p_R / (enthpy_R))
   Cs_L      = sqrt(gamma * p_L / (enthpy_L))

   V2_R      =  Vx_R**2 +  Vy_R**2 + Vz_R**2
   V2_L      =  Vx_L**2 +  Vy_L**2 + Vz_L**2

   lambda_H_p_R = ( Vx_R * ( 1.d0 - Cs_R**2) + Cs_R * sqrt( (1.d0 - V2_R) * (1.d0 - V2_R * Cs_R**2 - Vx_R * (1.d0 - Cs_R**2) )) ) &
                / ( 1.d0 - Cs_R**2 * V2_R)

   lambda_H_p_L = ( Vx_L * ( 1.d0 - Cs_L**2) + Cs_L * sqrt( (1.d0 - V2_L) * (1.d0 - V2_L * Cs_L**2 - Vx_L * (1.d0 - Cs_L**2) )) ) &
                / ( 1.d0 - Cs_L**2 * V2_L)

   
   lambda_H_m_R = ( Vx_R * ( 1.d0 - Cs_R**2) - Cs_R * sqrt( (1.d0 - V2_R) * (1.d0 - V2_R * Cs_R**2 - Vx_R * (1.d0 - Cs_R**2) )) ) &
                / ( 1.d0 - Cs_R**2 * V2_R)

   lambda_H_m_L = ( Vx_L * ( 1.d0 - Cs_L**2) - Cs_L * sqrt( (1.d0 - V2_L) * (1.d0 - V2_L * Cs_L**2 - Vx_L * (1.d0 - Cs_L**2) )) ) &
                / ( 1.d0 - Cs_L**2 * V2_L)


   s1 = min(lambda_H_m_L, lambda_H_m_R, 0.d0)
   s2 = max(lambda_H_p_L, lambda_H_p_R, 0.d0)

!!$   print*, "s1", s1, i
!!$   print*, "s2", s2, i

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

    psistr_flux_R =  Exstr_R  !psi_hll_flux  
    psistr_flux_L =  Exstr_L  !psi_hll_flux  
    
    phistr_flux_R =  Bxstr_R  !phi_hll_flux  
    phistr_flux_L =  Bxstr_L  !phi_hll_flux  

    Exstr_flux_R  =  psistr_R !Ex_hll_flux  
    Exstr_flux_L  =  psistr_L !Ex_hll_flux  

    Eystr_flux_R  =  Bzstr_R  !Ey_hll_flux  
    Eystr_flux_L  =  Bzstr_L  !Ey_hll_flux  

    Ezstr_flux_R  = -Bystr_R  !Ez_hll_flux  
    Ezstr_flux_L  = -Bystr_L  !Ez_hll_flux  

    Bxstr_flux_R  =  phistr_R !Bx_hll_flux  
    Bxstr_flux_L  =  phistr_L !Bx_hll_flux  

    Bystr_flux_R  = -Ezstr_R  !By_hll_flux  
    Bystr_flux_L  = -Ezstr_L  ! By_hll_flux  

    Bzstr_flux_R  =  Eystr_R  !Bz_hll_flux  
    Bzstr_flux_L  =  Eystr_L  !Bz_hll_flux

    

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


    FDxstr_R      = FDxint_R  + s2 * ( DDstr_R  -  DDint_R) ! DDstr_R * Vxstr !
    FDxstr_L      = FDxint_L  + s1 * ( DDstr_L  -  DDint_L) ! DDstr_L * Vxstr !


    ! PDF NOTES (pag 80) VERSION
    
!!$    FDxstr_L = DD_hll_flux  - s1 * (s2 - Vxstr) * (DDstr_R - DDstr_L) / (s2 -s1)
!!$    FDxstr_R = DD_hll_flux  + s2 * (Vxstr - s1) * (DDstr_R - DDstr_L) / (s2 -s1)


!   these values must be canceled when we made diference of fluxes --->

    FDystr_R      = 1.d0
    FDystr_L      = 1.d0
    FDzstr_R      = 1.d0
    FDzstr_L      = 1.d0

!   Flujos conservados de energía

    Ftauxstr_R    = Sxstr_R !Ftauxint_R + s2 * (taustr_R - tauint_R) ! 
    Ftauxstr_L    = Sxstr_L !Ftauxint_L + s1 * (taustr_L - tauint_L) ! 

    ! PDF NOTES (pag 80) VERSION

!!$    Ftauxstr_L = tau_hll_flux  - s1 * (s2 - Vxstr) * (taustr_R - taustr_L) / (s2 -s1)
!!$    Ftauxstr_R = tau_hll_flux  + s2 * (Vxstr - s1) * (taustr_R - taustr_L) / (s2 -s1)


!   these values must be canceled when we made diference of fluxes --->

    Ftauystr_R    = 1.d0
    Ftauystr_L    = 1.d0
    Ftauzstr_R    = 1.d0
    Ftauzstr_L    = 1.d0



!   Components of the flux momentum tensor 

    FSxxstr_R     = FSxxint_R + s2 * ( Sxstr_R - Sxint_R ) ! - Exstr_R**2 - Bxstr_R**2 + ( Sxstr_R - ErotB_str_x ) * Vxstr + pstr !
    FSxxstr_L     = FSxxint_L + s1 * ( Sxstr_L - Sxint_L ) ! - Exstr_L**2 - Bxstr_L**2 + ( Sxstr_L - ErotB_str_x ) * Vxstr + pstr !

    FSxystr_R     = FSxyint_R + s2 * ( Systr_R - Syint_R ) ! - Exstr_R * Eystr_R - Bxstr_R * Bystr_R + (Systr_R - ErotB_str_y) * Vxstr !
    FSxystr_L     = FSxyint_L + s1 * ( Systr_L - Syint_L ) ! - Exstr_L * Eystr_L - Bxstr_L * Bystr_L + (Systr_L - ErotB_str_y) * Vxstr !

    FSxzstr_R     = FSxzint_R + s2 * ( Szstr_R - Szint_R ) ! - Exstr_R * Ezstr_R - Bxstr_R * Bzstr_R + (Szstr_R - ErotB_str_z) * Vxstr !
    FSxzstr_L     = FSxzint_L + s1 * ( Szstr_L - Szint_L ) ! - Exstr_L * Ezstr_L - Bxstr_L * Bzstr_L + (Szstr_L - ErotB_str_z) * Vxstr !


!!$    print*, "s1=", s1, i
!!$    print*, "s2=", s2, i

!     PDF NOTES (pag 80) VERSION

!!$    FSxxstr_L =  Sx_hll_flux  - s1 * (s2 - Vxstr) * (Sxstr_R - Sxstr_L) / (s2 -s1)
!!$    FSxxstr_R =  Sx_hll_flux  + s2 * (Vxstr - s1) * (Sxstr_R - Sxstr_L) / (s2 -s1)
!!$
!!$    FSxystr_L =  Sy_hll_flux  - s1 * (s2 - Vxstr) * (Systr_R - Systr_L) / (s2 -s1)
!!$    FSxystr_R =  Sy_hll_flux  + s2 * (Vxstr - s1) * (Systr_R - Systr_L) / (s2 -s1)
!!$
!!$    FSxzstr_L =  Sz_hll_flux  - s1 * (s2 - Vxstr) * (Szstr_R - Szstr_L) / (s2 -s1)
!!$    FSxzstr_R =  Sz_hll_flux  + s2 * (Vxstr - s1) * (Szstr_R - Szstr_L) / (s2 -s1)

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



!******************************************     EM FIELDS       ************************************************************

    ! Finally we use the HLLC criterion for flux


    if ( -1.d0 .lt. 0.d0 .and. 0.d0 .le. Vxstr ) then

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

! Electric field sources
!_________________________________________________________________________

       if (source_rec_mpx == 1) then
   
          Exastsour(i,j,1,l) = Exsource_L

          Eyastsour(i,j,1,l) = Eysource_L

          Ezastsour(i,j,1,l) = Ezsource_L

       end if


       

    else if (Vxstr .lt. 0.d0 .and. 0.d0 .le. 1.d0) then


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
       

! Electric field sources
!_________________________________________________________________________

       if (source_rec_mpx == 1) then
   
          Exastsour(i,j,1,l) = Exsource_R

          Eyastsour(i,j,1,l) = Eysource_R

          Ezastsour(i,j,1,l) = Ezsource_R

       end if

    else 

       print*, "Revise Vxstar velocity it could be imaginary subroutine hllc_flow "
       print*, "Vxstar =", Vxstr
       print*, "s1=", s1
       print*, "s2=", s2
       stop

    end if




    
!******************************************  HD FIELDS + CHARGE ************************************************************

    ! Finally we use the HLLC criterion for flux


    if ( s1 .gt. 0.d0 ) then

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


    else if ( s1 .le. 0.d0 .and. 0.d0 .le. Vxstr ) then

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


       
!$OMP   PARALLEL DEFAULT(SHARED) &
!$OMP & PRIVATE(dim_x,dim_y,rho__,not__,V2mp5,V2mp5_x,V2mp5_y,sigma_L,sigma_R,test_rho_L,test_rho_R),            &
!$OMP & PRIVATE(dm4jph, dm4jmh, vul, vav, vmd, vlc, vmin, vmax, s1, s2, sx1, sx2, sy1,sy2, eps, eps1, eps2),     &
!$OMP & PRIVATE(Cs_x_R, Cs_x_L, Vx_x_R, Vx_x_L, lambda_H_p_x_R, lambda_H_p_x_L, lambda_H_m_x_R, lambda_H_m_x_L), &
!$OMP & PRIVATE(Cs_y_R, Cs_y_L, Vy_y_R, Vy_y_L, lambda_H_p_y_R, lambda_H_p_y_L, lambda_H_m_y_R, lambda_H_m_y_L)

      
!//////////////////////////// X - Direction ///////////////////////////////

       varname = "not__"

         
!$OMP DO   

       do imp=-1,jmax+1
          

            u(-6:imax+6, 1,imp) =  psiint(-6:imax+6,imp,1,l)
            u(-6:imax+6, 2,imp) =  phiint(-6:imax+6,imp,1,l)
            
            u(-6:imax+6, 3,imp) =  Exint (-6:imax+6,imp,1,l)
            u(-6:imax+6, 4,imp) =  Eyint (-6:imax+6,imp,1,l)
            u(-6:imax+6, 5,imp) =  Ezint (-6:imax+6,imp,1,l)

            u(-6:imax+6, 6,imp) =- Exint (-6:imax+6,imp,1,l)
            u(-6:imax+6, 7,imp) =- Eyint (-6:imax+6,imp,1,l)
            u(-6:imax+6, 8,imp) =- Ezint (-6:imax+6,imp,1,l)

            u(-6:imax+6, 9,imp) =  Bxint (-6:imax+6,imp,1,l)
            u(-6:imax+6,10,imp) =  Byint (-6:imax+6,imp,1,l)
            u(-6:imax+6,11,imp) =  Bzint (-6:imax+6,imp,1,l)

            u(-6:imax+6,12,imp) =- Bxint (-6:imax+6,imp,1,l)
            u(-6:imax+6,13,imp) =- Byint (-6:imax+6,imp,1,l)
            u(-6:imax+6,14,imp) =- Bzint (-6:imax+6,imp,1,l)

            u(-6:imax+6,15,imp) =  qint  (-6:imax+6,imp,1,l)

            u(-6:imax+6,16,imp) =  DDint (-6:imax+6,imp,1,l)

            u(-6:imax+6,17,imp) =  tauint(-6:imax+6,imp,1,l)

            u(-6:imax+6,18,imp) =  Sxint (-6:imax+6,imp,1,l)
            u(-6:imax+6,19,imp) =  Syint (-6:imax+6,imp,1,l)
            u(-6:imax+6,20,imp) =  Szint (-6:imax+6,imp,1,l)


            if (REC_PRIM == 0) then

               u(-6:imax+6,21,imp) = Vx (-6:imax+6,imp,1)
               u(-6:imax+6,22,imp) = Vy (-6:imax+6,imp,1)
               u(-6:imax+6,23,imp) = Vz (-6:imax+6,imp,1)

               u(-6:imax+6,24,imp) = p  (-6:imax+6,imp,1)
               u(-6:imax+6,25,imp) = rho(-6:imax+6,imp,1)

            else if (REC_PRIM == 1) then

               u(-6:imax+6,21,imp) = Vxint (-6:imax+6,imp,1,l)
               u(-6:imax+6,22,imp) = Vyint (-6:imax+6,imp,1,l)
               u(-6:imax+6,23,imp) = Vzint (-6:imax+6,imp,1,l)

               u(-6:imax+6,24,imp) = pint  (-6:imax+6,imp,1,l)
               u(-6:imax+6,25,imp) = rhoint(-6:imax+6,imp,1,l)

            else

               write(*,*) "STOP: subroutine hllflux"
               write(*,*) "REC_PRIM parameter is not valid"
               stop

            end if

!*** reconstruction only in x direction

               u(-6:imax+6,26,imp) = Jxint(-6:imax+6,imp,1,l)
               u(-6:imax+6,27,imp) = Jzint(-6:imax+6,imp,1,l)
!!$
!!$               u(-6:imax+6,26,imp) = Jx(-6:imax+6,imp,1)
!!$               u(-6:imax+6,27,imp) = Jz(-6:imax+6,imp,1)

!*** Condicional para reconstruir flujos usando mp5 o valores a left y right 

            if ( flux_rec_mp5 == 1 ) then


               u(-6:imax+6,28,imp) = FDxint  (-6:imax+6,imp,1,l)
               u(-6:imax+6,29,imp) = FDzint  (-6:imax+6,imp,1,l)

               u(-6:imax+6,30,imp) = Ftauxint(-6:imax+6,imp,1,l)
               u(-6:imax+6,31,imp) = Ftauzint(-6:imax+6,imp,1,l)
               
               u(-6:imax+6,32,imp) = FSxxint (-6:imax+6,imp,1,l)
               u(-6:imax+6,33,imp) = FSxyint (-6:imax+6,imp,1,l)
               u(-6:imax+6,34,imp) = FSxzint (-6:imax+6,imp,1,l)

               u(-6:imax+6,35,imp) = FSzxint (-6:imax+6,imp,1,l)
               u(-6:imax+6,36,imp) = FSzyint (-6:imax+6,imp,1,l)
               u(-6:imax+6,37,imp) = FSzzint (-6:imax+6,imp,1,l)

            end if
            
            coordenate = 1

          if (rec_mp == 0) then
            
            call MP5(imp)

         else if (rec_mp == 1) then

            call MP7(imp)

         else if (rec_mp == 2) then

            call MP9(imp)
            
         else

            print*, " rec_mp parameter not valid! "
            stop

         end if
            

            psiint_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1, 1,imp) 
            psiint_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2, 1,imp)

            phiint_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1, 2,imp) 
            phiint_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2, 2,imp)

            Exint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1, 3,imp) 
            Exint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2, 3,imp)
            Eyint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1, 4,imp) 
            Eyint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2, 4,imp)
            Ezint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1, 5,imp) 
            Ezint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2, 5,imp)

            Exint_mp5_x_L_m(-1:imax+1,imp) = up(-1:imax+1, 6,imp) 
            Exint_mp5_x_R_m(-1:imax+1,imp) = um( 0:imax+2, 6,imp)
            Eyint_mp5_x_L_m(-1:imax+1,imp) = up(-1:imax+1, 7,imp) 
            Eyint_mp5_x_R_m(-1:imax+1,imp) = um( 0:imax+2, 7,imp)
            Ezint_mp5_x_L_m(-1:imax+1,imp) = up(-1:imax+1, 8,imp) 
            Ezint_mp5_x_R_m(-1:imax+1,imp) = um( 0:imax+2, 8,imp)

            Bxint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1, 9,imp) 
            Bxint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2, 9,imp)
            Byint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,10,imp) 
            Byint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,10,imp)
            Bzint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,11,imp) 
            Bzint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,11,imp)
            
            Bxint_mp5_x_L_m(-1:imax+1,imp) = up(-1:imax+1,12,imp) 
            Bxint_mp5_x_R_m(-1:imax+1,imp) = um( 0:imax+2,12,imp)
            Byint_mp5_x_L_m(-1:imax+1,imp) = up(-1:imax+1,13,imp) 
            Byint_mp5_x_R_m(-1:imax+1,imp) = um( 0:imax+2,13,imp)
            Bzint_mp5_x_L_m(-1:imax+1,imp) = up(-1:imax+1,14,imp) 
            Bzint_mp5_x_R_m(-1:imax+1,imp) = um( 0:imax+2,14,imp)

            qint_mp5_x_L   (-1:imax+1,imp) = up(-1:imax+1,15,imp) 
            qint_mp5_x_R   (-1:imax+1,imp) = um( 0:imax+2,15,imp)
            
            DDint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,16,imp) 
            DDint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,16,imp)

            tauint_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,17,imp) 
            tauint_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,17,imp)

            Sxint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,18,imp) 
            Sxint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,18,imp)
            Syint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,19,imp) 
            Syint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,19,imp)
            Szint_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,20,imp) 
            Szint_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,20,imp) 

            if (REC_PRIM == 0) then

               Vx_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,21,imp) 
               Vx_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,21,imp)
               Vy_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,22,imp) 
               Vy_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,22,imp)
               Vz_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,23,imp) 
               Vz_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,23,imp)

               p_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,24,imp) 
               p_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,24,imp)
               rho_mp5_x_L(-1:imax+1,imp) = up(-1:imax+1,25,imp) 
               rho_mp5_x_R(-1:imax+1,imp) = um( 0:imax+2,25,imp) 

            else if (REC_PRIM == 1) then

               Vx_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,21,imp) 
               Vx_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,21,imp)
               Vy_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,22,imp) 
               Vy_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,22,imp)
               Vz_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,23,imp) 
               Vz_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,23,imp)

               p_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,24,imp) 
               p_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,24,imp)
               rho_mp5_x_L(-1:imax+1,imp) = up(-1:imax+1,25,imp) 
               rho_mp5_x_R(-1:imax+1,imp) = um( 0:imax+2,25,imp) 

            else

               write(*,*) "STOP: subroutine hllflux"
               write(*,*) "REC_PRIM parameter is not valid"
               stop

            end if


!*** reconstruction only in x direction

            Jxint_mp5_L(-1:imax+1,imp) = up(-1:imax+1,26,imp) 
            Jxint_mp5_R(-1:imax+1,imp) = um( 0:imax+2,26,imp)
            Jzint_mp5_L(-1:imax+1,imp) = up(-1:imax+1,27,imp) 
            Jzint_mp5_R(-1:imax+1,imp) = um( 0:imax+2,27,imp) 

!*** Condicional para reconstruir flujos usando mp5 o valores a left y right 

            if ( flux_rec_mp5 == 1 ) then

               
               FDxint_mp5_L  (-1:imax+1,imp) = up(-1:imax+1,28,imp) 
               FDxint_mp5_R  (-1:imax+1,imp) = um( 0:imax+2,28,imp)
               FDzint_mp5_L  (-1:imax+1,imp) = up(-1:imax+1,29,imp) 
               FDzint_mp5_R  (-1:imax+1,imp) = um( 0:imax+2,29,imp)

               Ftauxint_mp5_L(-1:imax+1,imp) = up(-1:imax+1,30,imp) 
               Ftauxint_mp5_R(-1:imax+1,imp) = um( 0:imax+2,30,imp)
               Ftauzint_mp5_L(-1:imax+1,imp) = up(-1:imax+1,31,imp) 
               Ftauzint_mp5_R(-1:imax+1,imp) = um( 0:imax+2,31,imp)

               FSxxint_mp5_L (-1:imax+1,imp) = up(-1:imax+1,32,imp) 
               FSxxint_mp5_R (-1:imax+1,imp) = um( 0:imax+2,32,imp)
               FSxyint_mp5_L (-1:imax+1,imp) = up(-1:imax+1,33,imp) 
               FSxyint_mp5_R (-1:imax+1,imp) = um( 0:imax+2,33,imp)
               FSxzint_mp5_L (-1:imax+1,imp) = up(-1:imax+1,34,imp) 
               FSxzint_mp5_R (-1:imax+1,imp) = um( 0:imax+2,34,imp)

               FSzxint_mp5_L (-1:imax+1,imp) = up(-1:imax+1,35,imp) 
               FSzxint_mp5_R (-1:imax+1,imp) = um( 0:imax+2,35,imp)
               FSzyint_mp5_L (-1:imax+1,imp) = up(-1:imax+1,36,imp) 
               FSzyint_mp5_R (-1:imax+1,imp) = um( 0:imax+2,36,imp)
               FSzzint_mp5_L (-1:imax+1,imp) = up(-1:imax+1,37,imp) 
               FSzzint_mp5_R (-1:imax+1,imp) = um( 0:imax+2,37,imp) 
               
               

            end if
            
         end do

!$OMP END DO


 
!//////////////////////////// Y - Direction ///////////////////////////////

        
                
!$OMP DO  

         do imp=-1,imax+1

            
            u(-6:jmax+6,38,imp) =  psiint(imp,-6:jmax+6,1,l)
            u(-6:jmax+6,39,imp) =  phiint(imp,-6:jmax+6,1,l)

            u(-6:jmax+6,40,imp) =  Exint (imp,-6:jmax+6,1,l)
            u(-6:jmax+6,41,imp) =  Eyint (imp,-6:jmax+6,1,l)
            u(-6:jmax+6,42,imp) =  Ezint (imp,-6:jmax+6,1,l)

            u(-6:jmax+6,43,imp) =- Exint (imp,-6:jmax+6,1,l)
            u(-6:jmax+6,44,imp) =- Eyint (imp,-6:jmax+6,1,l)
            u(-6:jmax+6,45,imp) =- Ezint (imp,-6:jmax+6,1,l)

            u(-6:jmax+6,46,imp) =  Bxint (imp,-6:jmax+6,1,l)
            u(-6:jmax+6,47,imp) =  Byint (imp,-6:jmax+6,1,l)
            u(-6:jmax+6,48,imp) =  Bzint (imp,-6:jmax+6,1,l)

            u(-6:jmax+6,49,imp) =- Bxint (imp,-6:jmax+6,1,l)
            u(-6:jmax+6,50,imp) =- Byint (imp,-6:jmax+6,1,l)
            u(-6:jmax+6,51,imp) =- Bzint (imp,-6:jmax+6,1,l)

            u(-6:jmax+6,52,imp) =  qint  (imp,-6:jmax+6,1,l)

            u(-6:jmax+6,53,imp) =  DDint (imp,-6:jmax+6,1,l)

            u(-6:jmax+6,54,imp) =  tauint(imp,-6:jmax+6,1,l)

            u(-6:jmax+6,55,imp) =  Sxint(imp,-6:jmax+6,1,l)
            u(-6:jmax+6,56,imp) =  Syint(imp,-6:jmax+6,1,l)
            u(-6:jmax+6,57,imp) =  Szint(imp,-6:jmax+6,1,l)

            if (REC_PRIM == 0) then

               u(-6:jmax+6,58,imp) = Vx (imp,-6:jmax+6,1)
               u(-6:jmax+6,59,imp) = Vy (imp,-6:jmax+6,1)
               u(-6:jmax+6,60,imp) = Vz (imp,-6:jmax+6,1)

               u(-6:jmax+6,61,imp) = p  (imp,-6:jmax+6,1)
               u(-6:jmax+6,62,imp) = rho(imp,-6:jmax+6,1)

            else if (REC_PRIM == 1) then

               u(-6:jmax+6,58,imp) = Vxint (imp,-6:jmax+6,1,l)
               u(-6:jmax+6,59,imp) = Vyint (imp,-6:jmax+6,1,l)
               u(-6:jmax+6,60,imp) = Vzint (imp,-6:jmax+6,1,l)

               u(-6:jmax+6,61,imp) = pint  (imp,-6:jmax+6,1,l)
               u(-6:jmax+6,62,imp) = rhoint(imp,-6:jmax+6,1,l)
               
            else

               write(*,*) "STOP: subroutine hllflux"
               write(*,*) "REC_PRIM parameter is not valid"
               stop

            end if


!*** reconstruction only in y direction

 
            u(-6:jmax+6,63,imp) = Jyint(imp,-6:jmax+6,1,l)
!!$
!!$            u(-6:jmax+6,63,imp) = Jy(imp,-6:jmax+6,1)      

!*** Condicional para reconstruir flujos usando mp5 o valores a left y right 

            if ( flux_rec_mp5 == 1 ) then

               u(-6:jmax+6,64,imp) = FDyint  (imp,-6:jmax+6,1,l)

               u(-6:jmax+6,65,imp) = Ftauyint(imp,-6:jmax+6,1,l)

               u(-6:jmax+6,66,imp) = FSyxint (imp,-6:jmax+6,1,l)
               u(-6:jmax+6,67,imp) = FSyyint (imp,-6:jmax+6,1,l)
               u(-6:jmax+6,68,imp) = FSyzint (imp,-6:jmax+6,1,l)

            end if

!___________________________________


            coordenate = 2

          if (rec_mp == 0) then
            
            call MP5(imp)

         else if (rec_mp == 1) then

            call MP7(imp)

         else if (rec_mp == 2) then

            call MP9(imp)
            
         else

            print*, " rec_mp parameter not valid! "
            stop

         end if


            psiint_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,38,imp) 
            psiint_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,38,imp)

            phiint_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,39,imp) 
            phiint_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,39,imp)

            Exint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,40,imp) 
            Exint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,40,imp)
            Eyint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,41,imp) 
            Eyint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,41,imp)
            Ezint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,42,imp) 
            Ezint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,42,imp)
            
            Exint_mp5_y_L_m(imp,-1:jmax+1) = up(-1:jmax+1,43,imp) 
            Exint_mp5_y_R_m(imp,-1:jmax+1) = um( 0:jmax+2,43,imp)
            Eyint_mp5_y_L_m(imp,-1:jmax+1) = up(-1:jmax+1,44,imp) 
            Eyint_mp5_y_R_m(imp,-1:jmax+1) = um( 0:jmax+2,44,imp)
            Ezint_mp5_y_L_m(imp,-1:jmax+1) = up(-1:jmax+1,45,imp) 
            Ezint_mp5_y_R_m(imp,-1:jmax+1) = um( 0:jmax+2,45,imp)

            Bxint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,46,imp) 
            Bxint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,46,imp)
            Byint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,47,imp) 
            Byint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,47,imp)
            Bzint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,48,imp) 
            Bzint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,48,imp)

            Bxint_mp5_y_L_m(imp,-1:jmax+1) = up(-1:jmax+1,49,imp) 
            Bxint_mp5_y_R_m(imp,-1:jmax+1) = um( 0:jmax+2,49,imp)
            Byint_mp5_y_L_m(imp,-1:jmax+1) = up(-1:jmax+1,50,imp) 
            Byint_mp5_y_R_m(imp,-1:jmax+1) = um( 0:jmax+2,50,imp)
            Bzint_mp5_y_L_m(imp,-1:jmax+1) = up(-1:jmax+1,51,imp) 
            Bzint_mp5_y_R_m(imp,-1:jmax+1) = um( 0:jmax+2,51,imp)

            qint_mp5_y_L   (imp,-1:jmax+1) = up(-1:jmax+1,52,imp) 
            qint_mp5_y_R   (imp,-1:jmax+1) = um( 0:jmax+2,52,imp)

            DDint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,53,imp) 
            DDint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,53,imp)

            tauint_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,54,imp) 
            tauint_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,54,imp)

            Sxint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,55,imp) 
            Sxint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,55,imp)
            Syint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,56,imp) 
            Syint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,56,imp)
            Szint_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,57,imp) 
            Szint_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,57,imp) 
            
            if (REC_PRIM == 0) then

               Vx_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,58,imp) 
               Vx_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,58,imp)
               Vy_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,59,imp) 
               Vy_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,59,imp)
               Vz_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,60,imp) 
               Vz_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,60,imp)

               p_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,61,imp) 
               p_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,61,imp)
               rho_mp5_y_L(imp,-1:jmax+1) = up(-1:jmax+1,62,imp) 
               rho_mp5_y_R(imp,-1:jmax+1) = um( 0:jmax+2,62,imp) 

            else if (REC_PRIM == 1) then

               Vx_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,58,imp) 
               Vx_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,58,imp)
               Vy_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,59,imp) 
               Vy_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,59,imp)
               Vz_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,60,imp) 
               Vz_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,60,imp)

               p_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,61,imp) 
               p_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,61,imp)
               rho_mp5_y_L(imp,-1:jmax+1) = up(-1:jmax+1,62,imp) 
               rho_mp5_y_R(imp,-1:jmax+1) = um( 0:jmax+2,62,imp) 
               

            else

               write(*,*) "STOP: subroutine hllflux"
               write(*,*) "REC_PRIM parameter is not valid"
               stop

            end if

!*** reconstruction only in y direction
            
               Jyint_mp5_L(imp,-1:jmax+1) = up(-1:jmax+1,63,imp) 
               Jyint_mp5_R(imp,-1:jmax+1) = um( 0:jmax+2,63,imp) 
            
!*** Condicional para reconstruir flujos usando mp5 o valores a left y right 

            if ( flux_rec_mp5 == 1 ) then

               FDyint_mp5_L  (imp,-1:jmax+1) = up(-1:jmax+1,64,imp) 
               FDyint_mp5_R  (imp,-1:jmax+1) = um( 0:jmax+2,64,imp)

               Ftauyint_mp5_L(imp,-1:jmax+1) = up(-1:jmax+1,65,imp) 
               Ftauyint_mp5_R(imp,-1:jmax+1) = um( 0:jmax+2,65,imp)

               FSyxint_mp5_L (imp,-1:jmax+1) = up(-1:jmax+1,66,imp) 
               FSyxint_mp5_R (imp,-1:jmax+1) = um( 0:jmax+2,66,imp)
               FSyyint_mp5_L (imp,-1:jmax+1) = up(-1:jmax+1,67,imp) 
               FSyyint_mp5_R (imp,-1:jmax+1) = um( 0:jmax+2,67,imp)
               FSyzint_mp5_L (imp,-1:jmax+1) = up(-1:jmax+1,68,imp) 
               FSyzint_mp5_R (imp,-1:jmax+1) = um( 0:jmax+2,68,imp) 
               

            end if

            
         end do

!$OMP END DO         

     
             
!$OMP DO   ORDERED       
 
          do j=-1,jmax
             do i=-1,imax

!$OMP ORDERED
                

!//////////////////////////// LORENTZ FACTOR ///////////////////////////////
                

             V2mp5_x = Vx_mp5_x_L(i,j)**2 +  Vy_mp5_x_L(i,j)**2 + Vz_mp5_x_L(i,j)**2
             V2mp5_y = Vx_mp5_y_L(i,j)**2 +  Vy_mp5_y_L(i,j)**2 + Vz_mp5_y_L(i,j)**2

            if (V2mp5_x .ge. 1.d0) then

               print*, "super luminal left velocity in mp5 soubroutine, V2_x =", V2mp5_x, "i=", i

            else if (V2mp5_y .ge. 1.d0) then

               print*, "super luminal left velocity in mp5 soubroutine, V2_Y =", V2mp5_y, "i=", i

            end if

            V2mp5_x = Vx_mp5_x_R(i,j)**2 +  Vy_mp5_x_R(i,j)**2 + Vz_mp5_x_R(i,j)**2
            V2mp5_y = Vx_mp5_y_R(i,j)**2 +  Vy_mp5_y_R(i,j)**2 + Vz_mp5_y_R(i,j)**2

            if (V2mp5_x .ge. 1.d0) then

             print*, "super luminal right velocity in mp5 soubroutine, V2_x =", V2mp5, "i=", i

            else if (V2mp5_y .ge. 1.d0) then

               print*, "super luminal left velocity in mp5 soubroutine, V2_Y =", V2mp5_y, "i=", i

            end if

             W_mp5_x_L(i,j) = 1.d0/sqrt(1.d0 - Vx_mp5_x_L(i,j)**2 -  Vy_mp5_x_L(i,j)**2 - Vz_mp5_x_L(i,j)**2)
             W_mp5_x_R(i,j) = 1.d0/sqrt(1.d0 - Vx_mp5_x_R(i,j)**2 -  Vy_mp5_x_R(i,j)**2 - Vz_mp5_x_R(i,j)**2)

             W_mp5_y_L(i,j) = 1.d0/sqrt(1.d0 - Vx_mp5_y_L(i,j)**2 -  Vy_mp5_y_L(i,j)**2 - Vz_mp5_y_L(i,j)**2)
             W_mp5_y_R(i,j) = 1.d0/sqrt(1.d0 - Vx_mp5_y_R(i,j)**2 -  Vy_mp5_y_R(i,j)**2 - Vz_mp5_y_R(i,j)**2)



!//////////////////////////// SQUARE EM FIELDS ///////////////////////////////


      E2B2int_mp5_x_L(i,j) = Exint_mp5_x_L(i,j)**2 + Eyint_mp5_x_L(i,j)**2 + Ezint_mp5_x_L(i,j)**2 + &
                             Bxint_mp5_x_L(i,j)**2 + Byint_mp5_x_L(i,j)**2 + Bzint_mp5_x_L(i,j)**2

      E2B2int_mp5_x_R(i,j) = Exint_mp5_x_R(i,j)**2 + Eyint_mp5_x_R(i,j)**2 + Ezint_mp5_x_R(i,j)**2 + &
                             Bxint_mp5_x_R(i,j)**2 + Byint_mp5_x_R(i,j)**2 + Bzint_mp5_x_R(i,j)**2

      E2B2int_mp5_y_L(i,j) = Exint_mp5_y_L(i,j)**2 + Eyint_mp5_y_L(i,j)**2 + Ezint_mp5_y_L(i,j)**2 + &
                             Bxint_mp5_y_L(i,j)**2 + Byint_mp5_y_L(i,j)**2 + Bzint_mp5_y_L(i,j)**2

      E2B2int_mp5_y_R(i,j) = Exint_mp5_y_R(i,j)**2 + Eyint_mp5_y_R(i,j)**2 + Ezint_mp5_y_R(i,j)**2 + &
                             Bxint_mp5_y_R(i,j)**2 + Byint_mp5_y_R(i,j)**2 + Bzint_mp5_y_R(i,j)**2


      Edotv_mp5_x_L(i,j) = Exint_mp5_x_L(i,j) * Vx_mp5_x_L(i,j) + Eyint_mp5_x_L(i,j) * Vy_mp5_x_L(i,j) &
                         + Ezint_mp5_x_L(i,j) * Vz_mp5_x_L(i,j)
      Edotv_mp5_x_R(i,j) = Exint_mp5_x_R(i,j) * Vx_mp5_x_R(i,j) + Eyint_mp5_x_R(i,j) * Vy_mp5_x_R(i,j) &
                         + Ezint_mp5_x_R(i,j) * Vz_mp5_x_R(i,j)

      Edotv_mp5_y_L(i,j) = Exint_mp5_y_L(i,j) * Vx_mp5_y_L(i,j) + Eyint_mp5_y_L(i,j) * Vy_mp5_y_L(i,j) &
                         + Ezint_mp5_y_L(i,j) * Vz_mp5_y_L(i,j)
      Edotv_mp5_y_R(i,j) = Exint_mp5_y_R(i,j) * Vx_mp5_y_R(i,j) + Eyint_mp5_y_R(i,j) * Vy_mp5_y_R(i,j) &
                         + Ezint_mp5_y_R(i,j) * Vz_mp5_y_R(i,j)



      vrotB_x_mp5_x_L(i,j) = (Bzint_mp5_x_L(i,j) * Vy_mp5_x_L(i,j) - Byint_mp5_x_L(i,j) * Vz_mp5_x_L(i,j))
      vrotB_y_mp5_x_L(i,j) = (Bxint_mp5_x_L(i,j) * Vz_mp5_x_L(i,j) - Bzint_mp5_x_L(i,j) * Vx_mp5_x_L(i,j))
      vrotB_z_mp5_x_L(i,j) = (Byint_mp5_x_L(i,j) * Vx_mp5_x_L(i,j) - Bxint_mp5_x_L(i,j) * Vy_mp5_x_L(i,j))

      vrotB_x_mp5_x_R(i,j) = (Bzint_mp5_x_R(i,j) * Vy_mp5_x_R(i,j) - Byint_mp5_x_R(i,j) * Vz_mp5_x_R(i,j))
      vrotB_y_mp5_x_R(i,j) = (Bxint_mp5_x_R(i,j) * Vz_mp5_x_R(i,j) - Bzint_mp5_x_R(i,j) * Vx_mp5_x_R(i,j))
      vrotB_z_mp5_x_R(i,j) = (Byint_mp5_x_R(i,j) * Vx_mp5_x_R(i,j) - Bxint_mp5_x_R(i,j) * Vy_mp5_x_R(i,j))

      vrotB_x_mp5_y_L(i,j) = (Bzint_mp5_y_L(i,j) * Vy_mp5_y_L(i,j) - Byint_mp5_y_L(i,j) * Vz_mp5_y_L(i,j))
      vrotB_y_mp5_y_L(i,j) = (Bxint_mp5_y_L(i,j) * Vz_mp5_y_L(i,j) - Bzint_mp5_y_L(i,j) * Vx_mp5_y_L(i,j))
      vrotB_z_mp5_y_L(i,j) = (Byint_mp5_y_L(i,j) * Vx_mp5_y_L(i,j) - Bxint_mp5_y_L(i,j) * Vy_mp5_y_L(i,j))

      vrotB_x_mp5_y_R(i,j) = (Bzint_mp5_y_R(i,j) * Vy_mp5_y_R(i,j) - Byint_mp5_y_R(i,j) * Vz_mp5_y_R(i,j))
      vrotB_y_mp5_y_R(i,j) = (Bxint_mp5_y_R(i,j) * Vz_mp5_y_R(i,j) - Bzint_mp5_y_R(i,j) * Vx_mp5_y_R(i,j))
      vrotB_z_mp5_y_R(i,j) = (Byint_mp5_y_R(i,j) * Vx_mp5_y_R(i,j) - Bxint_mp5_y_R(i,j) * Vy_mp5_y_R(i,j))



              if (SLC == 1) then
                 
                 sigma_L = sigma_0 * DDint_mp5_x_L(i,j)**gamma_slc !You must revise this we reconstruct only in x direction
                 sigma_R = sigma_0 * DDint_mp5_x_R(i,j)**gamma_slc 

              else if (SLC == 0) then

                 sigma_L = sigma
                 sigma_R = sigma

              end if

!Source reconstruction
!------------------------------------------------------------------------------------------------------              
      Exsour_L(i,j)    = sigma * W_mp5_x_L(i,j) * ( Exint_mp5_x_L(i,j) + vrotB_x_mp5_x_L(i,j) - &
                         Edotv_mp5_x_L(i,j) * Vx_mp5_x_L(i,j)  ) 
      Exsour_R(i,j)    = sigma * W_mp5_x_R(i,j) * ( Exint_mp5_x_R(i,j) + vrotB_x_mp5_x_R(i,j) - &
                         Edotv_mp5_x_R(i,j) * Vx_mp5_x_R(i,j)  )

      Eysour_L(i,j)    = sigma * W_mp5_y_L(i,j) * ( Eyint_mp5_y_L(i,j) + vrotB_y_mp5_y_L(i,j) - &
                         Edotv_mp5_y_L(i,j) * Vy_mp5_y_L(i,j)  )
      Eysour_R(i,j)    = sigma * W_mp5_y_R(i,j) * ( Eyint_mp5_y_R(i,j) + vrotB_y_mp5_y_R(i,j) - &
                         Edotv_mp5_y_R(i,j) * Vy_mp5_y_R(i,j)  )

      Ezsour_L(i,j)    = sigma * W_mp5_x_L(i,j) * ( Ezint_mp5_x_L(i,j) + vrotB_z_mp5_x_L(i,j) - &
                         Edotv_mp5_x_L(i,j) * Vz_mp5_x_L(i,j)  ) 
      Ezsour_R(i,j)    = sigma * W_mp5_x_R(i,j) * ( Ezint_mp5_x_R(i,j) + vrotB_z_mp5_x_R(i,j) - &
                         Edotv_mp5_x_R(i,j) * Vz_mp5_x_R(i,j)  ) 

              
 

!!$      Jxint_mp5_L(i,j) = sigma * W_mp5_x_L(i,j) * ( Exint_mp5_x_L(i,j) + vrotB_x_mp5_x_L(i,j) - &
!!$                         Edotv_mp5_x_L(i,j) * Vx_mp5_x_L(i,j)  )+ qint_mp5_x_L(i,j) * Vx_mp5_x_L(i,j) 
!!$      Jxint_mp5_R(i,j) = sigma * W_mp5_x_R(i,j) * ( Exint_mp5_x_R(i,j) + vrotB_x_mp5_x_R(i,j) - &
!!$                         Edotv_mp5_x_R(i,j) * Vx_mp5_x_R(i,j)  )+ qint_mp5_x_R(i,j) * Vx_mp5_x_R(i,j)
!!$
!!$      Jyint_mp5_L(i,j) = sigma * W_mp5_y_L(i,j) * ( Eyint_mp5_y_L(i,j) + vrotB_y_mp5_y_L(i,j) - &
!!$                         Edotv_mp5_y_L(i,j) * Vy_mp5_y_L(i,j)  )+ qint_mp5_y_L(i,j) * Vy_mp5_y_L(i,j)
!!$      Jyint_mp5_R(i,j) = sigma * W_mp5_y_R(i,j) * ( Eyint_mp5_y_R(i,j) + vrotB_y_mp5_y_R(i,j) - &
!!$                         Edotv_mp5_y_R(i,j) * Vy_mp5_y_R(i,j)  )+ qint_mp5_y_R(i,j) * Vy_mp5_y_R(i,j)
!!$
!!$      Jzint_mp5_L(i,j) = sigma * W_mp5_x_L(i,j) * ( Ezint_mp5_x_L(i,j) + vrotB_z_mp5_x_L(i,j) - &
!!$                         Edotv_mp5_x_L(i,j) * Vz_mp5_x_L(i,j)  )+  qint_mp5_x_L(i,j) * Vz_mp5_x_L(i,j) 
!!$      Jzint_mp5_R(i,j) = sigma * W_mp5_x_R(i,j) * ( Ezint_mp5_x_R(i,j) + vrotB_z_mp5_x_R(i,j) - &
!!$                         Edotv_mp5_x_R(i,j) * Vz_mp5_x_R(i,j)  )+  qint_mp5_x_R(i,j) * Vz_mp5_x_R(i,j) 


             test_rho_L = rho_mp5_y_L(i,j)
             test_rho_R = rho_mp5_y_R(i,j)

            if (test_rho_L .lt. 0.d0) then

               print*, "Negative density mp5 test_rho_L=", rho_mp5_y_L(i,j), "i=", i, "j=", j

            else if (test_rho_R .lt. 0.d0) then

              print*, "Negative density mp5 test_rho_R=", rho_mp5_y_R(i,j), "i=", i, "j=", j

            end if

 
!//////////////////////////// ENTHALPY ///////////////////////////////

 
    epsilon_mp5_x_L(i,j) = p_mp5_x_L(i,j) / ((gamma-1.d0) * rho_mp5_x_L(i,j))
    epsilon_mp5_x_R(i,j) = p_mp5_x_R(i,j) / ((gamma-1.d0) * rho_mp5_x_R(i,j))

    enthpy_mp5_x_L(i,j)  = rho_mp5_x_L(i,j) * ( 1.d0 + epsilon_mp5_x_L(i,j) ) + p_mp5_x_L(i,j) 
    enthpy_mp5_x_R(i,j)  = rho_mp5_x_R(i,j) * ( 1.d0 + epsilon_mp5_x_R(i,j) ) + p_mp5_x_R(i,j) 

    epsilon_mp5_y_L(i,j) = p_mp5_y_L(i,j) / ((gamma-1.d0) * rho_mp5_y_L(i,j))
    epsilon_mp5_y_R(i,j) = p_mp5_y_R(i,j) / ((gamma-1.d0) * rho_mp5_y_R(i,j))

    enthpy_mp5_y_L(i,j)  = rho_mp5_y_L(i,j) * ( 1.d0 + epsilon_mp5_y_L(i,j) ) + p_mp5_y_L(i,j) 
    enthpy_mp5_y_R(i,j)  = rho_mp5_y_R(i,j) * ( 1.d0 + epsilon_mp5_y_R(i,j) ) + p_mp5_y_R(i,j) 

 

!*********** ELSE del Condicional para reconstruir flujos usando mp5 o valores a left y right 

         if ( flux_rec_mp5 == 0 ) then
 

    !     Density flux (FD = \rho W V)

    FDxint_mp5_L(i,j) = DDint_mp5_x_L(i,j) * Vx_mp5_x_L(i,j)
    FDxint_mp5_R(i,j) = DDint_mp5_x_R(i,j) * Vx_mp5_x_R(i,j)

    FDyint_mp5_L(i,j) = DDint_mp5_y_L(i,j) * Vy_mp5_y_L(i,j)
    FDyint_mp5_R(i,j) = DDint_mp5_y_R(i,j) * Vy_mp5_y_R(i,j)

    FDzint_mp5_L(i,j) = DDint_mp5_x_L(i,j) * Vz_mp5_x_L(i,j)
    FDzint_mp5_R(i,j) = DDint_mp5_x_R(i,j) * Vz_mp5_x_R(i,j)!---> in 2D we are using x reconstruction component
                                                            !---> this must be zero when realize de flux diference in z component


    !     Flujos conservados de energía

    Ftauxint_mp5_L(i,j) = (Bzint_mp5_x_L(i,j)  * Eyint_mp5_x_L(i,j) - Byint_mp5_x_L(i,j) * Ezint_mp5_x_L(i,j)) +   &
                           enthpy_mp5_x_L(i,j) * W_mp5_x_L(i,j)**2  * Vx_mp5_x_L(i,j)   
    Ftauxint_mp5_R(i,j) = (Bzint_mp5_x_R(i,j)  * Eyint_mp5_x_R(i,j) - Byint_mp5_x_R(i,j) * Ezint_mp5_x_R(i,j)) +   &
                           enthpy_mp5_x_R(i,j) * W_mp5_x_R(i,j)**2  * Vx_mp5_x_R(i,j)   

    Ftauyint_mp5_L(i,j) = (Bxint_mp5_y_L(i,j)  * Ezint_mp5_y_L(i,j) - Bzint_mp5_y_L(i,j) * Exint_mp5_y_L(i,j)) +   &
                           enthpy_mp5_y_L(i,j) * W_mp5_y_L(i,j)**2  * Vy_mp5_y_L(i,j)
    Ftauyint_mp5_R(i,j) = (Bxint_mp5_y_R(i,j)  * Ezint_mp5_y_R(i,j) - Bzint_mp5_y_R(i,j) * Exint_mp5_y_R(i,j)) +   &
                           enthpy_mp5_y_R(i,j) * W_mp5_y_R(i,j)**2  * Vy_mp5_y_R(i,j)

    Ftauzint_mp5_L(i,j) = (Byint_mp5_x_L(i,j)  * Exint_mp5_x_L(i,j) - Bxint_mp5_x_L(i,j) * Eyint_mp5_x_L(i,j)) +   &
                           enthpy_mp5_x_L(i,j) * W_mp5_x_L(i,j)**2  * Vz_mp5_x_L(i,j)
    Ftauzint_mp5_R(i,j) = (Byint_mp5_x_R(i,j)  * Exint_mp5_x_R(i,j) - Bxint_mp5_x_R(i,j) * Eyint_mp5_x_R(i,j)) +   &
                           enthpy_mp5_x_R(i,j) * W_mp5_x_R(i,j)**2  * Vz_mp5_x_R(i,j)

    !     Components of the flux momentum tensor 

    FSxxint_mp5_L(i,j) =  - Exint_mp5_x_L(i,j)**2 - Bxint_mp5_x_L(i,j)**2 +                &
                            enthpy_mp5_x_L(i,j) * W_mp5_x_L(i,j)**2 * Vx_mp5_x_L(i,j)**2 + &
                            0.5d0 * E2B2int_mp5_x_L(i,j) + p_mp5_x_L(i,j)  
    FSxxint_mp5_R(i,j) =  - Exint_mp5_x_R(i,j)**2 - Bxint_mp5_x_R(i,j)**2 +                &
                            enthpy_mp5_x_R(i,j) * W_mp5_x_R(i,j)**2 * Vx_mp5_x_R(i,j)**2 + &
                            0.5d0 * E2B2int_mp5_x_R(i,j) + p_mp5_x_R(i,j)  

    FSxyint_mp5_L(i,j) = -  Exint_mp5_x_L(i,j)  * Eyint_mp5_x_L(i,j) - Bxint_mp5_x_L(i,j) * Byint_mp5_x_L(i,j) +    &
                            enthpy_mp5_x_L(i,j) * W_mp5_x_L(i,j)**2  * Vx_mp5_x_L(i,j)    * Vy_mp5_x_L(i,j)
    FSxyint_mp5_R(i,j) = -  Exint_mp5_x_R(i,j)  * Eyint_mp5_x_R(i,j) - Bxint_mp5_x_R(i,j) * Byint_mp5_x_R(i,j) +    &
                            enthpy_mp5_x_R(i,j) * W_mp5_x_R(i,j)**2  * Vx_mp5_x_R(i,j)    * Vy_mp5_x_R(i,j)

    FSxzint_mp5_L(i,j) = -  Exint_mp5_x_L(i,j)  * Ezint_mp5_x_L(i,j) - Bxint_mp5_x_L(i,j) * Bzint_mp5_x_L(i,j) +    &
                            enthpy_mp5_x_L(i,j) * W_mp5_x_L(i,j)**2  * Vx_mp5_x_L(i,j)    * Vz_mp5_x_L(i,j)
    FSxzint_mp5_R(i,j) = -  Exint_mp5_x_R(i,j)  * Ezint_mp5_x_R(i,j) - Bxint_mp5_x_R(i,j) * Bzint_mp5_x_R(i,j) +    &
                            enthpy_mp5_x_R(i,j) * W_mp5_x_R(i,j)**2  * Vx_mp5_x_R(i,j)    * Vz_mp5_x_R(i,j)

    FSyxint_mp5_L(i,j) = -  Exint_mp5_y_L(i,j)  * Eyint_mp5_y_L(i,j) - Bxint_mp5_y_L(i,j) * Byint_mp5_y_L(i,j) +    &
                            enthpy_mp5_y_L(i,j) * W_mp5_y_L(i,j)**2  * Vx_mp5_y_L(i,j)    * Vy_mp5_y_L(i,j)
    FSyxint_mp5_R(i,j) = -  Exint_mp5_y_R(i,j)  * Eyint_mp5_y_R(i,j) - Bxint_mp5_y_R(i,j) * Byint_mp5_y_R(i,j) +    &
                            enthpy_mp5_y_R(i,j) * W_mp5_y_R(i,j)**2  * Vx_mp5_y_R(i,j)    * Vy_mp5_y_R(i,j)

    FSyyint_mp5_L(i,j) = -  Eyint_mp5_y_L(i,j)**2 - Byint_mp5_y_L(i,j)**2 +                &
                            enthpy_mp5_y_L(i,j) * W_mp5_y_L(i,j)**2 * Vy_mp5_y_L(i,j)**2 + &
                            0.5d0 * E2B2int_mp5_y_L(i,j)  + p_mp5_y_L(i,j)
    FSyyint_mp5_R(i,j) = -  Eyint_mp5_y_R(i,j)**2 - Byint_mp5_y_R(i,j)**2 +                &
                            enthpy_mp5_y_R(i,j) * W_mp5_y_R(i,j)**2 * Vy_mp5_y_R(i,j)**2 + &
                            0.5d0 * E2B2int_mp5_y_R(i,j)  + p_mp5_y_R(i,j)

    FSyzint_mp5_L(i,j) = -  Eyint_mp5_y_L(i,j)  * Ezint_mp5_y_L(i,j) - Byint_mp5_y_L(i,j) * Bzint_mp5_y_L(i,j) +    &
                            enthpy_mp5_y_L(i,j) * W_mp5_y_L(i,j)**2  * Vy_mp5_y_L(i,j)    * Vz_mp5_y_L(i,j)
    FSyzint_mp5_R(i,j) = -  Eyint_mp5_y_R(i,j)  * Ezint_mp5_y_R(i,j) - Byint_mp5_y_R(i,j) * Bzint_mp5_y_R(i,j) +    &
                            enthpy_mp5_y_R(i,j) * W_mp5_y_R(i,j)**2  * Vy_mp5_y_R(i,j)    * Vz_mp5_y_R(i,j)

    FSzxint_mp5_L(i,j) = -  Exint_mp5_x_L(i,j)  * Ezint_mp5_x_L(i,j) - Bxint_mp5_x_L(i,j) * Bzint_mp5_x_L(i,j) +    &
                            enthpy_mp5_x_L(i,j) * W_mp5_x_L(i,j)**2  * Vx_mp5_x_L(i,j)    * Vz_mp5_x_L(i,j)       
    FSzxint_mp5_R(i,j) = -  Exint_mp5_x_R(i,j)  * Ezint_mp5_x_R(i,j) - Bxint_mp5_x_R(i,j) * Bzint_mp5_x_R(i,j) +    &
                            enthpy_mp5_x_R(i,j) * W_mp5_x_R(i,j)**2  * Vx_mp5_x_R(i,j)    * Vz_mp5_x_R(i,j)             

    FSzyint_mp5_L(i,j) = -  Eyint_mp5_x_L(i,j)  * Ezint_mp5_x_L(i,j) - Byint_mp5_x_L(i,j) * Bzint_mp5_x_L(i,j) +    &
                            enthpy_mp5_x_L(i,j) * W_mp5_x_L(i,j)**2  * Vy_mp5_x_L(i,j)    * Vz_mp5_x_L(i,j)
    FSzyint_mp5_R(i,j) = -  Eyint_mp5_x_R(i,j)  * Ezint_mp5_x_R(i,j) - Byint_mp5_x_R(i,j) * Bzint_mp5_x_R(i,j) +    &
                            enthpy_mp5_x_R(i,j) * W_mp5_x_R(i,j)**2  * Vy_mp5_x_R(i,j)    * Vz_mp5_x_R(i,j)

    FSzzint_mp5_L(i,j) = -  Ezint_mp5_x_L(i,j)**2 - Bzint_mp5_x_L(i,j)**2 +                &
                            enthpy_mp5_x_L(i,j) * W_mp5_x_L(i,j)**2 * Vz_mp5_x_L(i,j)**2 + &
                            0.5d0 * E2B2int_mp5_x_L(i,j) + p_mp5_x_L(i,j)
    FSzzint_mp5_R(i,j) = -  Ezint_mp5_x_R(i,j)**2 - Bzint_mp5_x_R(i,j)**2 +                &
                            enthpy_mp5_x_R(i,j) * W_mp5_x_R(i,j)**2 * Vz_mp5_x_R(i,j)**2 + &
                            0.5d0 * E2B2int_mp5_x_R(i,j) + p_mp5_x_R(i,j)


        else   if (flux_rec_mp5 .gt. 1) then 

       print*, " flux_rec_mp5 parameter not valid in subroutine hllc_flow 2D"
       stop

    end if
        


!/////////////////////// HLL Variables  X DIRECTION ////////////////////////////////////

    ! we difine de characteristic velocitis in diferent form from notes as -----> s2 = lambda_r  and s1 = lambda_l

    ! since for HLLE Riemann solver the values of s1 and s2 change, we must force to its values of -1 and 1 for EM fields
    !-----------------------------------------------------------------------------------------------------
    s1 = -1.d0
    s2 =  1.d0
    !-----------------------------------------------------------------------------------------------------

    !     hll  values of augmented fields
    psi_x_hll = ( s2 * psiint_mp5_x_R(i,j) - s1 * psiint_mp5_x_L(i,j) &
         + Exint_mp5_x_L(i,j)  - Exint_mp5_x_R(i,j)   ) / (s2 -s1)
    phi_x_hll = ( s2 * phiint_mp5_x_R(i,j) - s1 * phiint_mp5_x_L(i,j) &
         + Bxint_mp5_x_L(i,j)  - Bxint_mp5_x_R(i,j)   ) / (s2 -s1)

    !     hll  values of electric field

    Ex_x_hll  = ( s2 * Exint_mp5_x_R(i,j)  - s1 * Exint_mp5_x_L(i,j)  &
         + psiint_mp5_x_L(i,j)  - psiint_mp5_x_R(i,j)  ) / (s2 -s1)
    Ey_x_hll  = ( s2 * Eyint_mp5_x_R(i,j)  - s1 * Eyint_mp5_x_L(i,j)  &
         + Bzint_mp5_x_L(i,j)   - Bzint_mp5_x_R(i,j)   ) / (s2 -s1)
    Ez_x_hll  = ( s2 * Ezint_mp5_x_R(i,j)  - s1 * Ezint_mp5_x_L(i,j)  &
         + Byint_mp5_x_L_m(i,j) - Byint_mp5_x_R_m(i,j)   ) / (s2 -s1)  !ojo minus sign ¡¡¡


    !     hll  values of magnetic field

    Bx_x_hll  = ( s2 * Bxint_mp5_x_R(i,j)  - s1 * Bxint_mp5_x_L(i,j)  &
         + phiint_mp5_x_L(i,j)  - phiint_mp5_x_R(i,j)  ) / (s2 -s1)
    By_x_hll  = ( s2 * Byint_mp5_x_R(i,j)  - s1 * Byint_mp5_x_L(i,j)  &
         + Ezint_mp5_x_L_m(i,j) - Ezint_mp5_x_R_m(i,j)   ) / (s2 -s1)  !ojo minus sign ¡¡¡
    Bz_x_hll  = ( s2 * Bzint_mp5_x_R(i,j)  - s1 * Bzint_mp5_x_L(i,j)  &
         + Eyint_mp5_x_L(i,j)   - Eyint_mp5_x_R(i,j)   ) / (s2 -s1)

!/////////////////////// HLL Variables  Y DIRECTION ////////////////////////////////////

    ! we difine de characteristic velocitis in diferent form from notes as -----> s2 = lambda_r  and s1 = lambda_l

    !     hll  values of augmented fields
    psi_y_hll = ( s2 * psiint_mp5_y_R(i,j) - s1 * psiint_mp5_y_L(i,j) &
         + Eyint_mp5_y_L(i,j)  - Eyint_mp5_y_R(i,j)   ) / (s2 -s1)
    phi_y_hll = ( s2 * phiint_mp5_y_R(i,j) - s1 * phiint_mp5_y_L(i,j) &
         + Byint_mp5_y_L(i,j)  - Byint_mp5_y_R(i,j)   ) / (s2 -s1)

    !     hll  values of electric field

    Ey_y_hll  = ( s2 * Eyint_mp5_y_R(i,j)  - s1 * Eyint_mp5_y_L(i,j)  &
         + psiint_mp5_y_L(i,j)  - psiint_mp5_y_R(i,j)  ) / (s2 -s1)
    Ez_y_hll  = ( s2 * Ezint_mp5_y_R(i,j)  - s1 * Ezint_mp5_y_L(i,j)  &
         + Bxint_mp5_y_L(i,j)   - Bxint_mp5_y_R(i,j)   ) / (s2 -s1)
    Ex_y_hll  = ( s2 * Exint_mp5_y_R(i,j)  - s1 * Exint_mp5_y_L(i,j)  &
         + Bzint_mp5_y_L_m(i,j) - Bzint_mp5_y_R_m(i,j)   ) / (s2 -s1)  !ojo minus sign ¡¡¡


    !     hll  values of magnetic field

    By_y_hll  = ( s2 * Byint_mp5_y_R(i,j)  - s1 * Byint_mp5_y_L(i,j)  &
         + phiint_mp5_y_L(i,j)  - phiint_mp5_y_R(i,j)  ) / (s2 -s1)
    Bz_y_hll  = ( s2 * Bzint_mp5_y_R(i,j)  - s1 * Bzint_mp5_y_L(i,j)  &
         + Exint_mp5_y_L_m(i,j) - Exint_mp5_y_R_m(i,j)   ) / (s2 -s1)  !ojo minus sign ¡¡¡
    Bx_y_hll  = ( s2 * Bxint_mp5_y_R(i,j)  - s1 * Bxint_mp5_y_L(i,j)  &
         + Ezint_mp5_y_L(i,j)   - Ezint_mp5_y_R(i,j)   ) / (s2 -s1)

!/////////////////////// HLL Fluxes  X DIRECTION ////////////////////////////////////

! Electric Gauss HLL flux
!_________________________________________________________________________

      psi_x_hll_flux  =   0.5d0 * (Exint_mp5_x_L(i,j) + Exint_mp5_x_R(i,j) &
           - ( psiint_mp5_x_R(i,j) - psiint_mp5_x_L(i,j) ) )  

! Magnetic Gauss HLL flux
!_________________________________________________________________________


      phi_x_hll_flux  =   0.5d0 * (Bxint_mp5_x_L(i,j) + Bxint_mp5_x_R(i,j) &
           - ( phiint_mp5_x_R(i,j) - phiint_mp5_x_L(i,j) ) )

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

       Bx_x_hll_flux  =   0.5d0 * (phiint_mp5_x_L(i,j)   + phiint_mp5_x_R(i,j)  &
            - ( Bxint_mp5_x_R(i,j) - Bxint_mp5_x_L(i,j) ) )

! Flujos By

       By_x_hll_flux  =   0.5d0 * (Ezint_mp5_x_L_m(i,j)  + Ezint_mp5_x_R_m(i,j) &
            - ( Byint_mp5_x_R(i,j) - Byint_mp5_x_L(i,j) ) ) ! minus sign

! Flujos Bz

       Bz_x_hll_flux  =   0.5d0 * (Eyint_mp5_x_L(i,j)    + Eyint_mp5_x_R(i,j)  &
            - ( Bzint_mp5_x_R(i,j) - Bzint_mp5_x_L(i,j) ) )

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       Ex_x_hll_flux  =   0.5d0 * (psiint_mp5_x_L(i,j)   + psiint_mp5_x_R(i,j) &
            - ( Exint_mp5_x_R(i,j) - Exint_mp5_x_L(i,j) ) )

! Flujos Ey


       Ey_x_hll_flux  =   0.5d0 * (Bzint_mp5_x_L(i,j)    + Bzint_mp5_x_R(i,j)  &
            - ( Eyint_mp5_x_R(i,j) - Eyint_mp5_x_L(i,j) ) )

! Flujos Ez


       Ez_x_hll_flux  =   0.5d0 * (Byint_mp5_x_L_m(i,j)  + Byint_mp5_x_R_m(i,j) &
            - ( Ezint_mp5_x_R(i,j) - Ezint_mp5_x_L(i,j) ) ) ! minus sign


!Source
!_________________________________________________________________________
       
      Exsource_L      = Exsour_L(i,j)    
      Exsource_R      = Exsour_R(i,j)   
                      
      Eysource_L      = Eysour_L(i,j)   
      Eysource_R      = Eysour_R(i,j)   

      Ezsource_L      = Ezsour_L(i,j)   
      Ezsource_R      = Ezsour_R(i,j)   
    

!/////////////////////// HLL Fluxes  Y DIRECTION ////////////////////////////////////

! Electric Gauss HLL flux
!_________________________________________________________________________

      psi_y_hll_flux  =   0.5d0 * (Eyint_mp5_y_L(i,j) + Eyint_mp5_y_R(i,j) &
           - ( psiint_mp5_y_R(i,j) - psiint_mp5_y_L(i,j) ) )  

! Magnetic Gauss HLL flux
!_________________________________________________________________________


      phi_y_hll_flux  =   0.5d0 * (Byint_mp5_y_L(i,j) + Byint_mp5_y_R(i,j) &
           - ( phiint_mp5_y_R(i,j) - phiint_mp5_y_L(i,j) ) )

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos By 

       By_y_hll_flux  =   0.5d0 * (phiint_mp5_y_L(i,j)   + phiint_mp5_y_R(i,j) &
            - ( Byint_mp5_y_R(i,j) - Byint_mp5_y_L(i,j) ) )

! Flujos Bz

       Bz_y_hll_flux  =   0.5d0 * (Exint_mp5_y_L_m(i,j)  + Exint_mp5_y_R_m(i,j) &
            - ( Bzint_mp5_y_R(i,j) - Bzint_mp5_y_L(i,j) ) ) ! minus sign

! Flujos Bx

       Bx_y_hll_flux  =   0.5d0 * (Ezint_mp5_y_L(i,j)    + Ezint_mp5_y_R(i,j) &
            - ( Bxint_mp5_y_R(i,j) - Bxint_mp5_y_L(i,j) ) )

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ey

       Ey_y_hll_flux  =   0.5d0 * (psiint_mp5_y_L(i,j)   + psiint_mp5_y_R(i,j) &
            - ( Eyint_mp5_y_R(i,j) - Eyint_mp5_y_L(i,j) ) )

! Flujos Ez

       Ez_y_hll_flux  =   0.5d0 * (Bxint_mp5_y_L(i,j)    + Bxint_mp5_y_R(i,j) &
            - ( Ezint_mp5_y_R(i,j) - Ezint_mp5_y_L(i,j) ) )

! Flujos Ex

       Ex_y_hll_flux  =   0.5d0 * (Bzint_mp5_y_L_m(i,j)  + Bzint_mp5_y_R_m(i,j) &
            - ( Exint_mp5_y_R(i,j) - Exint_mp5_y_L(i,j) ) ) ! minus sign


!********************************** HIDRODYNAMICS  X DIRECTION *********************************************************

!=======================================================================================================================   
! characteristic sound velocity, with cs= sqrt(gamma * p/(rho * h), observe that you use h = rho * h
!=======================================================================================================================      
       
   Cs_x_R      = sqrt(gamma * p_mp5_x_R(i,j) / (enthpy_mp5_x_R(i,j) ))
   Cs_x_L      = sqrt(gamma * p_mp5_x_L(i,j) / (enthpy_mp5_x_L(i,j) ))

   Vx_x_R      = Vx_mp5_x_R(i,j)
   Vx_x_L      = Vx_mp5_x_L(i,j)
   

   V2_x_R      =  Vx_mp5_x_R(i,j)**2 +  Vy_mp5_x_R(i,j)**2 + Vz_mp5_x_R(i,j)**2
   V2_x_L      =  Vx_mp5_x_L(i,j)**2 +  Vy_mp5_x_L(i,j)**2 + Vz_mp5_x_L(i,j)**2

   lambda_H_p_x_R = ( Vx_x_R * ( 1.d0 - Cs_x_R**2) &
                  + Cs_x_R * sqrt( (1.d0 - V2_x_R) * (1.d0 - V2_x_R * Cs_x_R**2 - Vx_x_R * (1.d0 - Cs_x_R**2) )) ) &
                  / ( 1.d0 - Cs_x_R**2 * V2_x_R)

   lambda_H_p_x_L = ( Vx_x_L * ( 1.d0 - Cs_x_L**2) &
                  + Cs_x_L * sqrt( (1.d0 - V2_x_L) * (1.d0 - V2_x_L * Cs_x_L**2 - Vx_x_L * (1.d0 - Cs_x_L**2) )) ) &
                  / ( 1.d0 - Cs_x_L**2 * V2_x_L)

   
   lambda_H_m_x_R = ( Vx_x_R * ( 1.d0 - Cs_x_R**2) &
                  - Cs_x_R * sqrt( (1.d0 - V2_x_R) * (1.d0 - V2_x_R * Cs_x_R**2 - Vx_x_R * (1.d0 - Cs_x_R**2) )) ) &
                  / ( 1.d0 - Cs_x_R**2 * V2_x_R)

   lambda_H_m_x_L = ( Vx_x_L * ( 1.d0 - Cs_x_L**2) &
                  - Cs_x_L * sqrt( (1.d0 - V2_x_L) * (1.d0 - V2_x_L * Cs_x_L**2 - Vx_x_L * (1.d0 - Cs_x_L**2) )) ) &
                  / ( 1.d0 - Cs_x_L**2 * V2_x_L)


   !   These estimates are not recommended for practical computations, Eleuterio Toro Sec.10.5 Eq.10.48 3ed (2009)

   !   If appear NaN numbers. By definition, NAN is not equal to anything, even itself. Simply compare the variable to itself: if(my_var /= my_var)


!!$   eps = 1.d-1
!!$
!!$
!!$   if ( isnan(lambda_H_m_x_L) .or. isnan(lambda_H_m_x_R) .or. &
!!$        isnan(lambda_H_p_x_L) .or. isnan(lambda_H_p_x_R)    ) then
!!$
!!$      if ( isnan(lambda_H_m_x_L) ) lambda_H_m_x_L = -1.d0
!!$      if ( isnan(lambda_H_m_x_R) ) lambda_H_m_x_R = -1.d0
!!$
!!$      if ( isnan(lambda_H_p_x_L) ) lambda_H_p_x_L =  1.d0
!!$      if ( isnan(lambda_H_p_x_R) ) lambda_H_p_x_R =  1.d0
!!$
!!$      sx1 = min(lambda_H_m_x_L, lambda_H_m_x_R)
!!$      sx2 = max(lambda_H_p_x_L, lambda_H_p_x_R)
!!$
!!$      !  proposed average eigenvalues for the left and right non–linear waves, E. Toro (2009) 3ed
!!$      ! we made avarage over cuadri-vectors ux = 0.5 (Wl vl + ur vl)
!!$      ! find the lorentz factor -W^2 + ux^2 = -1 ---> W^2 = ux^2 + 1
!!$      ! and the three velocity as vx = sqrt( 1 - 1/W^2)
!!$
!!$   else if ( abs(lambda_H_m_x_L) .le. 1.d0 - eps .and. abs(lambda_H_m_x_R) .le. 1.d0 - eps .and. &
!!$             abs(lambda_H_p_x_L) .le. 1.d0 - eps .and. abs(lambda_H_p_x_R) .le. 1.d0 - eps     ) then
!!$
!!$      W_H_p_x_L = 1.d0 / sqrt(1.d0 - lambda_H_p_x_L**2)
!!$      W_H_p_x_R = 1.d0 / sqrt(1.d0 - lambda_H_p_x_R**2)
!!$
!!$      W_H_m_x_L = 1.d0 / sqrt(1.d0 - lambda_H_m_x_L**2)
!!$      W_H_m_x_R = 1.d0 / sqrt(1.d0 - lambda_H_m_x_R**2)
!!$
!!$      ux_p_m    = 0.5d0 * (W_H_p_x_L * lambda_H_p_x_L + W_H_p_x_R * lambda_H_p_x_R)
!!$      ux_m_m    = 0.5d0 * (W_H_m_x_L * lambda_H_m_x_L + W_H_m_x_R * lambda_H_m_x_R)
!!$
!!$      W2m_p     = ux_p_m**2  + 1.d0
!!$      W2m_m     = ux_m_m**2  + 1.d0
!!$
!!$
!!$      sx2       = sqrt( 1.d0 - 1.d0/W2m_p)
!!$      sx1       = sqrt( 1.d0 - 1.d0/W2m_m)
!!$
!!$   else

      sx1 = min(lambda_H_m_x_L, lambda_H_m_x_R)
      sx2 = max(lambda_H_p_x_L, lambda_H_p_x_R)

!!$   end if


   if ( isnan(sx1) .or. isnan(sx2)) then

      print*, " "
      print*, "sx1 ==", sx1, l, i, j, h, lambda_H_m_x_L, lambda_H_m_x_R, W_H_m_L, W_H_m_x_R
      print*, "sx2 ==", sx2, l, i, j, h, lambda_H_p_x_L, lambda_H_p_x_R, W_H_p_L, W_H_p_x_R
      print*, "NaN value for characteristic velocities sx1 and sx2. HLLC soubroutine"
      print*, " "

      stop

   end if

   eps2 = 1.d-9

   
   if ( abs(sx1) .gt. 1.d0 ) then

      sx1 = max(-1.d0 + eps2 , min( 1.d0 - eps2, sx1 ))
      
   end if

   if ( abs(sx2) .gt. 1.d0 ) then

      sx2 = max(-1.d0 + eps2 , min( 1.d0 - eps2, sx2 ))
      
   end if
    
!!$!===================
!!$   sx1 = -1.d0
!!$   sx2 =  1.d0
!!$!===================
   
!/////////////////////// HLL Variables  X DIRECTION ////////////////////////////////////
    
    !     hll  value of charge density 

    q_x_hll   = ( sx2 * qint_mp5_x_R(i,j)   - sx1 * qint_mp5_x_L(i,j)   + Jxint_mp5_L(i,j)   - Jxint_mp5_R(i,j)   ) / (sx2 -sx1)

    !     hll  value of mass conservation 

    DD_x_hll  = ( sx2 * DDint_mp5_x_R(i,j)  - sx1 * DDint_mp5_x_L(i,j)  + FDxint_mp5_L(i,j)  - FDxint_mp5_R(i,j)  ) / (sx2 -sx1)

    !     hll  value of energy conservation 

    tau_x_hll = ( sx2 * tauint_mp5_x_R(i,j) - sx1 * tauint_mp5_x_L(i,j) + Ftauxint_mp5_L(i,j) - Ftauxint_mp5_R(i,j) ) / (sx2 -sx1)

    !    hll  value of momentum S

    Sx_x_hll  = ( sx2 * Sxint_mp5_x_R(i,j)  - sx1 * Sxint_mp5_x_L(i,j)  + FSxxint_mp5_L(i,j)  - FSxxint_mp5_R(i,j)  ) / (sx2 -sx1)
    Sy_x_hll  = ( sx2 * Syint_mp5_x_R(i,j)  - sx1 * Syint_mp5_x_L(i,j)  + FSxyint_mp5_L(i,j)  - FSxyint_mp5_R(i,j)  ) / (sx2 -sx1)
    Sz_x_hll  = ( sx2 * Szint_mp5_x_R(i,j)  - sx1 * Szint_mp5_x_L(i,j)  + FSxzint_mp5_L(i,j)  - FSxzint_mp5_R(i,j)  ) / (sx2 -sx1)


!/////////////////////// HLL Fluxes  X DIRECTION ////////////////////////////////////

! Conserved current HLL flux
!_________________________________________________________________________


    q_x_hll_flux = (sx2 * Jxint_mp5_L(i,j) - sx1 * Jxint_mp5_R(i,j) &
                 +  sx1 * sx2 * (qint_mp5_x_R(i,j) - qint_mp5_x_L(i,j))) / (sx2 -sx1)

! Conserved Mass HLL flux
!_________________________________________________________________________


    DD_x_hll_flux = (sx2 * FDxint_mp5_L(i,j) - sx1 * FDxint_mp5_R(i,j) &
                  +  sx1 * sx2 * (DDint_mp5_x_R(i,j) - DDint_mp5_x_L(i,j))) / (sx2 -sx1)

! Conserved Energy HLL flux
!_________________________________________________________________________


    tau_x_hll_flux = (sx2 * Ftauxint_mp5_L(i,j) - sx1 * Ftauxint_mp5_R(i,j) &
                   +  sx1 * sx2 * (tauint_mp5_x_R(i,j) - tauint_mp5_x_L(i,j))) / (sx2 -sx1)


! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


    Sx_x_hll_flux  = (sx2 * FSxxint_mp5_L(i,j) - sx1 * FSxxint_mp5_R(i,j) &
                   +  sx1 * sx2 * (Sxint_mp5_x_R(i,j) - Sxint_mp5_x_L(i,j))) / (sx2 -sx1)
    
    Sy_x_hll_flux  = (sx2 * FSxyint_mp5_L(i,j) - sx1 * FSxyint_mp5_R(i,j) &
                   +  sx1 * sx2 * (Syint_mp5_x_R(i,j) - Syint_mp5_x_L(i,j))) / (sx2 -sx1)
    
    Sz_x_hll_flux  = (sx2 * FSxzint_mp5_L(i,j) - sx1 * FSxzint_mp5_R(i,j) &
                   +  sx1 * sx2 * (Szint_mp5_x_R(i,j) - Szint_mp5_x_L(i,j))) / (sx2 -sx1) 


!********************************** HIDRODYNAMICS  Y DIRECTION *********************************************************       

!=======================================================================================================================   
! characteristic sound velocity, with cs= sqrt(gamma * p/(rho * h), observe that you use h = rho * h
!=======================================================================================================================


   Cs_y_R      = sqrt(gamma * p_mp5_y_R(i,j) / (enthpy_mp5_y_R(i,j) ))
   Cs_y_L      = sqrt(gamma * p_mp5_y_L(i,j) / (enthpy_mp5_y_L(i,j) ))

   Vy_y_R      = Vy_mp5_y_R(i,j)
   Vy_y_L      = Vy_mp5_y_L(i,j)
   
   V2_y_R      =  Vx_mp5_y_R(i,j)**2 +  Vy_mp5_y_R(i,j)**2 + Vz_mp5_y_R(i,j)**2
   V2_y_L      =  Vx_mp5_y_L(i,j)**2 +  Vy_mp5_y_L(i,j)**2 + Vz_mp5_y_L(i,j)**2
 

   lambda_H_p_y_R = ( Vy_y_R * ( 1.d0 - Cs_y_R**2) &
                  + Cs_y_R * sqrt( (1.d0 - V2_y_R) * (1.d0 - V2_y_R * Cs_y_R**2 - Vy_y_R * (1.d0 - Cs_y_R**2) )) ) &
                  / ( 1.d0 - Cs_y_R**2 * V2_y_R)

   lambda_H_p_y_L = ( Vy_y_L * ( 1.d0 - Cs_y_L**2) &
                  + Cs_y_L * sqrt( (1.d0 - V2_y_L) * (1.d0 - V2_y_L * Cs_y_L**2 - Vy_y_L * (1.d0 - Cs_y_L**2) )) ) &
                  / ( 1.d0 - Cs_y_L**2 * V2_y_L)

   
   lambda_H_m_y_R = ( Vy_y_R * ( 1.d0 - Cs_y_R**2) &
                  - Cs_y_R * sqrt( (1.d0 - V2_y_R) * (1.d0 - V2_y_R * Cs_y_R**2 - Vy_y_R * (1.d0 - Cs_y_R**2) )) ) &
                  / ( 1.d0 - Cs_y_R**2 * V2_y_R)

   lambda_H_m_y_L = ( Vy_y_L * ( 1.d0 - Cs_y_L**2) &
                  - Cs_y_L * sqrt( (1.d0 - V2_y_L) * (1.d0 - V2_y_L * Cs_y_L**2 - Vy_y_L * (1.d0 - Cs_y_L**2) )) ) &
                  / ( 1.d0 - Cs_y_L**2 * V2_y_L)


   !   These estimates are not recommended for practical computations, Eleuterio Toro Sec.10.5 Eq.10.48 3ed (2009)

   !   If appear NaN numbers. By definition, NAN is not equal to anything, even itself. Simply compare the variable to itself: if(my_var /= my_var)


!!$   eps = 1.d-1
!!$
!!$
!!$   if ( isnan(lambda_H_m_y_L) .or. isnan(lambda_H_m_y_R) .or. &
!!$        isnan(lambda_H_p_y_L) .or. isnan(lambda_H_p_y_R)    ) then
!!$
!!$      if ( isnan(lambda_H_m_y_L) ) lambda_H_m_y_L = -1.d0
!!$      if ( isnan(lambda_H_m_y_R) ) lambda_H_m_y_R = -1.d0
!!$
!!$      if ( isnan(lambda_H_p_y_L) ) lambda_H_p_y_L =  1.d0
!!$      if ( isnan(lambda_H_p_y_R) ) lambda_H_p_y_R =  1.d0
!!$
!!$      sy1 = min(lambda_H_m_y_L, lambda_H_m_y_R)
!!$      sy2 = max(lambda_H_p_y_L, lambda_H_p_y_R)
!!$
!!$      ! proposed average eigenvalues for the left and right non–linear waves, E. Toro (2009) 3ed
!!$      ! we made avarage over cuadri-vectors ux = 0.5 (Wl vl + ur vl)
!!$      ! find the lorentz factor -W^2 + ux^2 = -1 ---> W^2 = ux^2 + 1
!!$      ! and the three velocity as vx = sqrt( 1 - 1/W^2)
!!$
!!$   else if ( abs(lambda_H_m_y_L) .le. 1.d0 - eps .and. abs(lambda_H_m_y_R) .le. 1.d0 - eps .and. &
!!$             abs(lambda_H_p_y_L) .le. 1.d0 - eps .and. abs(lambda_H_p_y_R) .le. 1.d0 - eps     ) then
!!$
!!$      W_H_p_y_L = 1.d0 / sqrt(1.d0 - lambda_H_p_y_L**2)
!!$      W_H_p_y_R = 1.d0 / sqrt(1.d0 - lambda_H_p_y_R**2)
!!$
!!$      W_H_m_y_L = 1.d0 / sqrt(1.d0 - lambda_H_m_y_L**2)
!!$      W_H_m_y_R = 1.d0 / sqrt(1.d0 - lambda_H_m_y_R**2)
!!$
!!$      uy_p_m    = 0.5d0 * (W_H_p_y_L * lambda_H_p_y_L + W_H_p_y_R * lambda_H_p_y_R)
!!$      uy_m_m    = 0.5d0 * (W_H_m_y_L * lambda_H_m_y_L + W_H_m_y_R * lambda_H_m_y_R)
!!$
!!$      W2m_p     = uy_p_m**2  + 1.d0
!!$      W2m_m     = uy_m_m**2  + 1.d0
!!$
!!$
!!$      sy2       = sqrt( 1.d0 - 1.d0/W2m_p)
!!$      sy1       = sqrt( 1.d0 - 1.d0/W2m_m)
!!$ 
!!$   else

      sy1 = min(lambda_H_m_y_L, lambda_H_m_y_R)
      sy2 = max(lambda_H_p_y_L, lambda_H_p_y_R)

!!$   end if


   if ( isnan(sy1) .or. isnan(sy2)) then

     print*, " "
      print*, "sy1 ==", sy1, l, i, j, h, lambda_H_m_x_L, lambda_H_m_x_R, W_H_m_L, W_H_m_x_R
      print*, "sy2 ==", sy2, l, i, j, h, lambda_H_p_x_L, lambda_H_p_x_R, W_H_p_L, W_H_p_x_R
      print*, "NaN value for characteristic velocities sy1 and sy2. HLLC soubroutine"
      print*, " "

      stop

   end if

   eps2 = 1.d-9

   
   if ( abs(sy1) .gt. 1.d0 ) then

      sy1 = max(-1.d0 + eps2 , min( 1.d0 - eps2, sy1 ))
      
   end if

   if ( abs(sy2) .gt. 1.d0 ) then

      sy2 = max(-1.d0 + eps2 , min( 1.d0 - eps2, sy2 ))
      
   end if

!!$!===================
!!$   sy1 = -1.d0
!!$   sy2 =  1.d0
!!$!===================   

       
!/////////////////////// HLL Variables  Y DIRECTION ////////////////////////////////////


    !     hll  value of charge density 

    q_y_hll   = ( sy2 * qint_mp5_y_R(i,j)   - sy1 * qint_mp5_y_L(i,j)   + Jyint_mp5_L(i,j)   - Jyint_mp5_R(i,j)   ) / (sy2 -sy1)

    !     hll  value of mass conservation 

    DD_y_hll  = ( sy2 * DDint_mp5_y_R(i,j)  - sy1 * DDint_mp5_y_L(i,j)  + FDyint_mp5_L(i,j)  - FDyint_mp5_R(i,j)  ) / (sy2 -sy1)

    !     hll  value of energy conservation 

    tau_y_hll = ( sy2 * tauint_mp5_y_R(i,j) - sy1 * tauint_mp5_y_L(i,j) + Ftauyint_mp5_L(i,j) - Ftauyint_mp5_R(i,j) ) / (sy2 -sy1)

    !    hll  value of momentum S

    Sy_y_hll  = ( sy2 * Syint_mp5_y_R(i,j)  - sy1 * Syint_mp5_y_L(i,j)  + FSyyint_mp5_L(i,j)  - FSyyint_mp5_R(i,j)  ) / (sy2 -sy1)
    Sz_y_hll  = ( sy2 * Szint_mp5_y_R(i,j)  - sy1 * Szint_mp5_y_L(i,j)  + FSyzint_mp5_L(i,j)  - FSyzint_mp5_R(i,j)  ) / (sy2 -sy1)
    Sx_y_hll  = ( sy2 * Sxint_mp5_y_R(i,j)  - sy1 * Sxint_mp5_y_L(i,j)  + FSyxint_mp5_L(i,j)  - FSyxint_mp5_R(i,j)  ) / (sy2 -sy1)


!/////////////////////// HLL Fluxes  Y DIRECTION ////////////////////////////////////

! Conserved current HLL flux
!_________________________________________________________________________


    q_y_hll_flux = (sy2 * Jyint_mp5_L(i,j) - sy1 * Jyint_mp5_R(i,j) &
                 +  sy1 * sy2 * (qint_mp5_y_R(i,j) - qint_mp5_y_L(i,j))) / (sy2 -sy1)

! Conserved Mass HLL flux
!_________________________________________________________________________


    DD_y_hll_flux = (sy2 * FDyint_mp5_L(i,j) - sy1 * FDyint_mp5_R(i,j) &
                  +  sy1 * sy2 * (DDint_mp5_y_R(i,j) - DDint_mp5_y_L(i,j))) / (sy2 -sy1)

! Conserved Energy HLL flux
!_________________________________________________________________________


    tau_y_hll_flux = (sy2 * Ftauyint_mp5_L(i,j) - sy1 * Ftauyint_mp5_R(i,j) &
                   +  sy1 * sy2 * (tauint_mp5_y_R(i,j) - tauint_mp5_y_L(i,j))) / (sy2 -sy1)


! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


    Sy_y_hll_flux = (sy2 * FSyyint_mp5_L(i,j) - sy1 * FSyyint_mp5_R(i,j) &
                  +  sy1 * sy2 * (Syint_mp5_y_R(i,j) - Syint_mp5_y_L(i,j))) / (sy2 -sy1)
    
    Sz_y_hll_flux = (sy2 * FSyzint_mp5_L(i,j) - sy1 * FSyzint_mp5_R(i,j) &
                  +  sy1 * sy2 * (Szint_mp5_y_R(i,j) - Szint_mp5_y_L(i,j))) / (sy2 -sy1)
    
    Sx_y_hll_flux = (sy2 * FSyxint_mp5_L(i,j) - sy1 * FSyxint_mp5_R(i,j) &
                  +  sy1 * sy2 * (Sxint_mp5_y_R(i,j) - Sxint_mp5_y_L(i,j))) / (sy2 -sy1)  


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

   psistr_x_R = psi_x_hll 
   psistr_x_L = psi_x_hll 
   
   phistr_x_R = phi_x_hll 
   phistr_x_L = phi_x_hll

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

   psistr_y_R = psi_y_hll 
   psistr_y_L = psi_y_hll 
   
   phistr_y_R = phi_y_hll 
   phistr_y_L = phi_y_hll
   

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


!   eps  = 1.d-10
   eps1 = 1.d-6
!   eps2 = 1.d-3

   
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


   if ( abs(Vxstr) .le. 1.d-9 ) then

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

      print*, "discriminante imaginario en la ecuacion cuadratica para obtener Vystr en HLLC Riemman solver (y direction)"
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


   if ( abs(Vystr) .le. 1.d-9 ) then

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
   
   Sxstr_x_R  = ( sx2 * Sx_x_hll - tau_x_hll_flux * Vxstr ) / (sx2 - Vxstr) 
   Sxstr_x_L  = (-sx1 * Sx_x_hll + tau_x_hll_flux * Vxstr ) / (Vxstr - sx1)

!!$   Sxstr_x_R  = (  sx2 * Sx_x_hll - Sx_x_hll_flux + pstr_x - Bxstr_x_R**2 - Exstr_x_R**2 - ErotB_str_x_x_R * Vxstr ) &
!!$              / (sx2 - Vxstr) 
!!$   Sxstr_x_L  = ( -sx1 * Sx_x_hll + Sx_x_hll_flux - pstr_x + Bxstr_x_R**2 + Exstr_x_R**2 + ErotB_str_x_x_R * Vxstr ) &
!!$              / (Vxstr - sx1) 
    
   Systr_x_R  = ( sx2 * Sy_x_hll - Sy_x_hll_flux - ErotB_str_y_x_R * Vxstr - Exstr_x_R * Eystr_x_R - Bxstr_x_R * Bystr_x_R )  &
              / ( sx2 - Vxstr)
   Systr_x_L  = (-sx1 * Sy_x_hll + Sy_x_hll_flux + ErotB_str_y_x_R * Vxstr + Exstr_x_L * Eystr_x_L + Bxstr_x_L * Bystr_x_L ) &
              / ( Vxstr - sx1)

   Szstr_x_R  = ( sx2 * Sz_x_hll - Sz_x_hll_flux - ErotB_str_z_x_R * Vxstr - Exstr_x_R * Ezstr_x_R - Bxstr_x_R * Bzstr_x_R )  &
              / ( sx2 - Vxstr)
   Szstr_x_L  = (-sx1 * Sz_x_hll + Sz_x_hll_flux + ErotB_str_z_x_R * Vxstr + Exstr_x_L * Ezstr_x_L + Bxstr_x_L * Bzstr_x_L ) &
              / ( Vxstr - sx1)

   taustr_x_R = (Sxstr_x_R - ErotB_str_x_x_R) / Vxstr - pstr_x + E2B2_str_x 
   taustr_x_L = (Sxstr_x_L - ErotB_str_x_x_R) / Vxstr - pstr_x + E2B2_str_x

!!$
!!$   taustr_x_R = (  sx2 * tau_x_hll - tau_x_hll_flux + (pstr_x - E2B2_str_x) * Vxstr + ErotB_str_x_x_R ) &
!!$              / (sx2 - Vxstr)
!!$   taustr_x_L = ( -sx1 * tau_x_hll + tau_x_hll_flux - (pstr_x - E2B2_str_x) * Vxstr - ErotB_str_x_x_R ) &
!!$              / (Vxstr - sx1)

   DDstr_x_R  = ( sx2 * DD_x_hll - DD_x_hll_flux ) / (sx2 - Vxstr)
   DDstr_x_L  = (-sx1 * DD_x_hll + DD_x_hll_flux ) / (Vxstr - sx1)

   ! Y direction
      
   Systr_y_R  = ( sy2 * Sy_y_hll - tau_y_hll_flux * Vystr ) / (sy2 - Vystr) 
   Systr_y_L  = (-sy1 * Sy_y_hll + tau_y_hll_flux * Vystr ) / (Vystr - sy1)
!!$
!!$   Systr_y_R  = (  sy2 * Sy_y_hll - Sy_y_hll_flux + pstr_y - Bystr_y_R**2 - Eystr_y_R**2 - ErotB_str_y_y_R * Vystr ) &
!!$              / (sy2 - Vystr) 
!!$   Systr_y_L  = ( -sy1 * Sy_y_hll + Sy_y_hll_flux - pstr_y + Bystr_y_R**2 + Eystr_y_R**2 + ErotB_str_y_y_R * Vystr ) &
!!$              / (Vystr - sy1) 
   
   Szstr_y_R  = ( sy2 * Sz_y_hll - Sz_y_hll_flux - ErotB_str_z_y_R * Vystr - Eystr_y_R * Ezstr_y_R - Bystr_y_R * Bzstr_y_R )  &
              / ( sy2 - Vystr)
   Szstr_y_L  = (-sy1 * Sz_y_hll + Sz_y_hll_flux + ErotB_str_z_y_R * Vystr + Eystr_y_L * Ezstr_y_L + Bystr_y_L * Bzstr_y_L ) &
              / (Vystr - sy1)

   Sxstr_y_R  = ( sy2 * Sx_y_hll - Sx_y_hll_flux - ErotB_str_x_y_R * Vystr - Eystr_y_R * Exstr_y_R - Bystr_y_R * Bxstr_y_R )  &
              / ( sy2 - Vystr)
   Sxstr_y_L  = (-sy1 * Sx_y_hll + Sx_y_hll_flux + ErotB_str_x_y_R * Vystr + Eystr_y_L * Exstr_y_L + Bystr_y_L * Bxstr_y_L )  &
              / (Vystr - sy1)

   taustr_y_R = (Systr_y_R - ErotB_str_y_y_R) / Vystr - pstr_y + E2B2_str_y 
   taustr_y_L = (Systr_y_L - ErotB_str_y_y_R) / Vystr - pstr_y + E2B2_str_y
!!$
!!$   taustr_y_R = (  sy2 * tau_y_hll - tau_y_hll_flux + (pstr_y - E2B2_str_y) * Vystr + ErotB_str_y_y_R ) &
!!$              / (sy2 - Vystr)
!!$   taustr_y_L = ( -sy1 * tau_y_hll + tau_y_hll_flux - (pstr_y - E2B2_str_y) * Vystr - ErotB_str_y_y_R ) &
!!$              / (Vystr - sy1)

   DDstr_y_R  = ( sy2 * DD_y_hll - DD_y_hll_flux ) / (sy2 - Vystr)
   DDstr_y_L  = (-sy1 * DD_y_hll + DD_y_hll_flux ) / (Vystr - sy1)


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

!!$    psistr_x_flux_R =  Exstr_x_R  !psi_x_hll_flux 
!!$    psistr_x_flux_L =  Exstr_x_L  !psi_x_hll_flux 
!!$    
!!$    phistr_x_flux_R =  Bxstr_x_R  !phi_x_hll_flux 
!!$    phistr_x_flux_L =  Bxstr_x_L  !phi_x_hll_flux 
!!$
!!$    Exstr_x_flux_R  =  psistr_x_R !Ex_x_hll_flux  
!!$    Exstr_x_flux_L  =  psistr_x_L !Ex_x_hll_flux  
!!$
!!$    Eystr_x_flux_R  =  Bzstr_x_R  !Ey_x_hll_flux  
!!$    Eystr_x_flux_L  =  Bzstr_x_L  !Ey_x_hll_flux  
!!$
!!$    Ezstr_x_flux_R  = -Bystr_x_R  ! Ez_x_hll_flux  
!!$    Ezstr_x_flux_L  = -Bystr_x_L  ! Ez_x_hll_flux  
!!$
!!$    Bxstr_x_flux_R  =  phistr_x_R !Bx_x_hll_flux  
!!$    Bxstr_x_flux_L  =  phistr_x_L !Bx_x_hll_flux  
!!$
!!$    Bystr_x_flux_R  = -Ezstr_x_R  !By_x_hll_flux  
!!$    Bystr_x_flux_L  = -Ezstr_x_R  !By_x_hll_flux  
!!$
!!$    Bzstr_x_flux_R  =  Eystr_x_R  !Bz_x_hll_flux  
!!$    Bzstr_x_flux_L  =  Eystr_x_L  !Bz_x_hll_flux


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

!!$    FDxstr_R      = FDxint_R  + sx2 * ( DDstr_x_R  -  DDint_x_R) ! DDstr_x_R * Vxstr !
!!$    FDxstr_L      = FDxint_L  + sx1 * ( DDstr_x_L  -  DDint_x_L) ! DDstr_x_L * Vxstr !


    ! PDF NOTES (pag 80) VERSION
    
    FDxstr_R = DD_x_hll_flux  + sx2 * (Vxstr - sx1) * (DDstr_x_R - DDstr_x_L) / (sx2 -sx1)
    FDxstr_L = DD_x_hll_flux  - sx1 * (sx2 - Vxstr) * (DDstr_x_R - DDstr_x_L) / (sx2 -sx1)

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

!!$    Ftauxstr_R    = Sxstr_x_R !Ftauxint_x_R + sx2 * (taustr_x_R - tauint_x_R) ! 
!!$    Ftauxstr_L    = Sxstr_x_L !Ftauxint_x_L + sx1 * (taustr_x_L - tauint_x_L) !

!!$    Ftauxstr_R    = Ftauxint_R + sx2 * (taustr_x_R - tauint_x_R) ! 
!!$    Ftauxstr_L    = Ftauxint_L + sx1 * (taustr_x_L - tauint_x_L) ! 

    ! PDF NOTES (pag 80) VERSION
    
    Ftauxstr_R = tau_x_hll_flux  + sx2 * (Vxstr - sx1) * (taustr_x_R - taustr_x_L) / (sx2 -sx1)
    Ftauxstr_L = tau_x_hll_flux  - sx1 * (sx2 - Vxstr) * (taustr_x_R - taustr_x_L) / (sx2 -sx1)

! these values must be canceled when we made diference of fluxes --->

    Ftauzstr_R = 1.d0
    Ftauzstr_L = 1.d0


    !     Components of the flux momentum tensor

!!$    FSxxstr_R     = FSxxint_R + sx2 * ( Sxstr_x_R - Sxint_x_R ) ! - Exstr_x_R**2 - Bxstr_x_R**2 + ( Sxstr_x_R - ErotB_str_x_x_R ) * Vxstr + pstr_x !
!!$    FSxxstr_L     = FSxxint_L + sx1 * ( Sxstr_x_L - Sxint_x_L ) ! - Exstr_x_L**2 - Bxstr_x_L**2 + ( Sxstr_x_L - ErotB_str_x_x_L ) * Vxstr + pstr_x !
!!$
!!$    FSxystr_R     = FSxyint_R + sx2 * ( Systr_x_R - Syint_x_R ) ! - Exstr_x_R * Eystr_x_R - Bxstr_x_R * Bystr_x_R + (Systr_x_R - ErotB_str_y_x_R) * Vxstr !
!!$    FSxystr_L     = FSxyint_L + sx1 * ( Systr_x_L - Syint_x_L ) ! - Exstr_x_L * Eystr_x_L - Bxstr_x_L * Bystr_x_L + (Systr_x_L - ErotB_str_y_x_L) * Vxstr !
!!$
!!$    FSxzstr_R     = FSxzint_R + sx2 * ( Szstr_x_R - Szint_x_R ) ! - Exstr_x_R * Ezstr_x_R - Bxstr_x_R * Bzstr_x_R + (Szstr_x_R - ErotB_str_z_x_R) * Vxstr !
!!$    FSxzstr_L     = FSxzint_L + sx1 * ( Szstr_x_L - Szint_x_L ) ! - Exstr_x_L * Ezstr_x_L - Bxstr_x_L * Bzstr_x_L + (Szstr_x_L - ErotB_str_z_x_L) * Vxstr !

!     PDF NOTES (pag 80) VERSION
    
    FSxxstr_R =  Sx_x_hll_flux  + sx2 * (Vxstr - sx1) * (Sxstr_x_R - Sxstr_x_L) / (sx2 -sx1)
    FSxxstr_L =  Sx_x_hll_flux  - sx1 * (sx2 - Vxstr) * (Sxstr_x_R - Sxstr_x_L) / (sx2 -sx1)

    FSxystr_R =  Sy_x_hll_flux  + sx2 * (Vxstr - sx1) * (Systr_x_R - Systr_x_L) / (sx2 -sx1)
    FSxystr_L =  Sy_x_hll_flux  - sx1 * (sx2 - Vxstr) * (Systr_x_R - Systr_x_L) / (sx2 -sx1)

    FSxzstr_R =  Sz_x_hll_flux  + sx2 * (Vxstr - sx1) * (Szstr_x_R - Szstr_x_L) / (sx2 -sx1)
    FSxzstr_L =  Sz_x_hll_flux  - sx1 * (sx2 - Vxstr) * (Szstr_x_R - Szstr_x_L) / (sx2 -sx1)

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
      

!!$    psistr_y_flux_R =  Eystr_y_R  !psi_y_hll_flux 
!!$    psistr_y_flux_L =  Eystr_y_L  !psi_y_hll_flux 
!!$    
!!$    phistr_y_flux_R =  Bystr_y_R  !phi_y_hll_flux 
!!$    phistr_y_flux_L =  Bystr_y_L  !phi_y_hll_flux 
!!$
!!$    Exstr_y_flux_R  = -Bzstr_y_R  !Ex_y_hll_flux  
!!$    Exstr_y_flux_L  = -Bzstr_y_L  ! Ex_y_hll_flux  
!!$
!!$    Eystr_y_flux_R  =  psistr_y_R !Ey_y_hll_flux  
!!$    Eystr_y_flux_L  =  psistr_y_L !Ey_y_hll_flux  
!!$
!!$    Ezstr_y_flux_R  =  Bxstr_y_R  !Ez_y_hll_flux  
!!$    Ezstr_y_flux_L  =  Bxstr_y_L  !Ez_y_hll_flux  
!!$
!!$    Bxstr_y_flux_R  =  Ezstr_y_R  !Bx_y_hll_flux  
!!$    Bxstr_y_flux_L  =  Ezstr_y_L  !Bx_y_hll_flux  
!!$
!!$    Bystr_y_flux_R  =  phistr_y_R !By_y_hll_flux  
!!$    Bystr_y_flux_L  =  phistr_y_L !By_y_hll_flux  
!!$
!!$    Bzstr_y_flux_R  = -Exstr_y_R  !Bz_y_hll_flux  
!!$    Bzstr_y_flux_L  = -Exstr_y_L  !Bz_y_hll_flux


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

!!$    FDystr_R      = FDyint_R  + sy2 * ( DDstr_y_R  -  DDint_y_R) ! DDstr_R * Vxstr !
!!$    FDystr_L      = FDyint_L  + sy1 * ( DDstr_y_L  -  DDint_y_L) ! DDstr_L * Vxstr !


    ! PDF NOTES (pag 80) VERSION

    FDystr_R = DD_y_hll_flux  + sy2 * (Vystr - sy1) * (DDstr_y_R - DDstr_y_L) / (sy2 -sy1)
    FDystr_L = DD_y_hll_flux  - sy1 * (sy2 - Vystr) * (DDstr_y_R - DDstr_y_L) / (sy2 -sy1)

    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jystr_L  = q_y_hll_flux 
    Jystr_R  = q_y_hll_flux 

    !     Flujos conservados de energía

!!$    Ftauystr_R    = Systr_R !Ftauxint_R + sy2 * (taustr_y_R - tauint_y_R) ! 
!!$    Ftauystr_L    = Systr_L !Ftauxint_L + sy1 * (taustr_y_L - tauint_y_L) !

!!$    Ftauystr_R    = Ftauyint_R + sy2 * (taustr_y_R - tauint_y_R) ! 
!!$    Ftauystr_L    = Ftauyint_L + sy1 * (taustr_y_L - tauint_y_L) ! 

    ! PDF NOTES (pag 80) VERSION
    
    Ftauystr_R = tau_y_hll_flux  + sy2 * (Vystr - sy1) * (taustr_y_R - taustr_y_L) / (sy2 -sy1)
    Ftauystr_L = tau_y_hll_flux  - sy1 * (sy2 - Vystr) * (taustr_y_R - taustr_y_L) / (sy2 -sy1)

    !     Components of the flux momentum tensor 

!!$    FSyxstr_R =  FSyxint_R + sy2 * ( Sxstr_y_R - Sxint_y_R )
!!$    FSyxstr_L =  FSyxint_L + sy1 * ( Sxstr_y_L - Sxint_y_L )
!!$
!!$    FSyystr_R =  FSyyint_R + sy2 * ( Systr_y_R - Syint_y_R )
!!$    FSyystr_L =  FSyyint_L + sy1 * ( Systr_y_L - Syint_y_L )
!!$
!!$    FSyzstr_R =  FSyzint_R + sy2 * ( Szstr_y_R - Szint_y_R )
!!$    FSyzstr_L =  FSyzint_L + sy1 * ( Szstr_y_L - Szint_y_L )

!     PDF NOTES (pag 80) VERSION

    FSyxstr_R =  Sx_y_hll_flux  + sy2 * (Vystr - sy1) * (Sxstr_y_R - Sxstr_y_L) / (sy2 -sy1)
    FSyxstr_L =  Sx_y_hll_flux  - sy1 * (sy2 - Vystr) * (Sxstr_y_R - Sxstr_y_L) / (sy2 -sy1)

    FSyystr_R =  Sy_y_hll_flux  + sy2 * (Vystr - sy1) * (Systr_y_R - Systr_y_L) / (sy2 -sy1)
    FSyystr_L =  Sy_y_hll_flux  - sy1 * (sy2 - Vystr) * (Systr_y_R - Systr_y_L) / (sy2 -sy1)

    FSyzstr_R =  Sz_y_hll_flux  + sy2 * (Vystr - sy1) * (Szstr_y_R - Szstr_y_L) / (sy2 -sy1)
    FSyzstr_L =  Sz_y_hll_flux  - sy1 * (sy2 - Vystr) * (Szstr_y_R - Szstr_y_L) / (sy2 -sy1)


    end if

!******************************************************************************************************************

    !/////////////////////// Finally we use the HLLC criterion for flux X DIECTION  //////////////////////////////////////////////



    !******************************************     EM FIELDS       ************************************************************


        if ( -1.d0 .lt. 0.d0 .and. 0.d0 .le. Vxstr ) then
    
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


! Electric field sources
!_________________________________________________________________________

       if (source_rec_mpx == 1) then
   
          Exastsour(i,j,1,l) = 0.5d0 * (Exsource_L - Exsource_R)

          Ezastsour(i,j,1,l) = 0.5d0 * (Ezsource_L - Ezsource_R)

       end if



   else if (Vxstr .lt. 0.d0 .and. 0.d0 .lt. 1.d0) then



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


! Electric field sources
!_________________________________________________________________________

       if (source_rec_mpx == 1) then
   
          Exastsour(i,j,1,l) = 0.5d0 * (Exsource_L - Exsource_R)

          Ezastsour(i,j,1,l) = 0.5d0 * (Ezsource_L - Ezsource_R)

       end if
       

    end if


!******************************************  HD FIELDS + CHARGE ************************************************************

    if ( sx1 .ge. 0.d0 ) then


! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,j,1,l)     =   Jxint_L  
       Jzlax(i,j,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,j,1,l)    =   FDxint_L  
       FDzlax(i,j,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,j,1,l)  =   Ftauxint_L 
       Ftauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,j,1,l)  =   psiint_x_L 
       psitauzlax(i,j,1,l)  =   1.d0


       phitauxlax(i,j,1,l)  =   phiint_x_L 
       phitauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,j,1,l)   =   FSxxint_L 
       FSxylax(i,j,1,l)   =   FSxyint_L 
       FSxzlax(i,j,1,l)   =   FSxzint_L 

 
       FSzxlax(i,j,1,l)  =   1.d0
       FSzylax(i,j,1,l)  =   1.d0
       FSzzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,j,1,l)   = Bxint_x_L  
       BzSxlax(i,j,1,l)   = 1.d0



       BxSylax(i,j,1,l)   = Bxint_x_L  
       BzSylax(i,j,1,l)   = 1.d0


       BxSzlax(i,j,1,l)   = Bxint_x_L  
       BzSzlax(i,j,1,l)   = 1.d0

! -------------------------------------- EGLM-------------------------------------- 

    else if ( sx1 .le. 0.d0 .and. 0.d0 .le. Vxstr ) then



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

    else if (Vxstr .le. 0.d0 .and. 0.d0 .le. sx2) then



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

       else if ( sx2 .lt. 0.d0 ) then

! Conserved current HLL flux
!_________________________________________________________________________


       Jxlax(i,j,1,l)     =   Jxint_R  
       Jzlax(i,j,1,l)     =   1.d0

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDxlax(i,j,1,l)    =   FDxint_R  
       FDzlax(i,j,1,l)    =   1.d0

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauxlax(i,j,1,l)  =   Ftauxint_R 
       Ftauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

       psitauxlax(i,j,1,l)  =   psiint_x_R 
       psitauzlax(i,j,1,l)  =   1.d0


       phitauxlax(i,j,1,l)  =   phiint_x_R 
       phitauzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSxxlax(i,j,1,l)   =   FSxxint_R 
       FSxylax(i,j,1,l)   =   FSxyint_R 
       FSxzlax(i,j,1,l)   =   FSxzint_R 

 
       FSzxlax(i,j,1,l)  =   1.d0
       FSzylax(i,j,1,l)  =   1.d0
       FSzzlax(i,j,1,l)  =   1.d0

! -------------------------------------- EGLM-------------------------------------- 


       BxSxlax(i,j,1,l)   = Bxint_x_R  
       BzSxlax(i,j,1,l)   = 1.d0



       BxSylax(i,j,1,l)   = Bxint_x_R  
       BzSylax(i,j,1,l)   = 1.d0


       BxSzlax(i,j,1,l)   = Bxint_x_R  
       BzSzlax(i,j,1,l)   = 1.d0

! -------------------------------------- EGLM-------------------------------------- 

          
    else 

       print*, "Revise Vxstar velocity it could be imaginary subroutine hlle_flow "
       stop

    end if


    !******************************************************************************************************************

    !/////////////////////// Finally we use the HLLC criterion for flux Y DIECTION  //////////////////////////////////////////////

    !******************************************     EM FIELDS       ************************************************************


        if ( -1.d0 .lt. 0.d0 .and. 0.d0 .le. Vystr ) then


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


! Electric field sources
!_________________________________________________________________________

       if (source_rec_mpx == 1) then
   
          Eyastsour(i,j,1,l) = 0.5d0 * (Eysource_L - Eysource_R)

       end if       
       

       
    else if (Vystr .lt. 0.d0 .and. 0.d0 .lt. 1.d0) then



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

! Electric field sources
!_________________________________________________________________________

       if (source_rec_mpx == 1) then
   
          Eyastsour(i,j,1,l) = 0.5d0 * (Eysource_L - Eysource_R)

       end if     
       

    end if


!******************************************  HD FIELDS + CHARGE ************************************************************

    if ( sy1 .ge. 0.d0 ) then



! Conserved current HLL flux
!_________________________________________________________________________


       Jylax(i,j,1,l)     =   Jyint_L 

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDylax(i,j,1,l)    =   FDyint_L

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauylax(i,j,1,l)  =   Ftauyint_L 

! -------------------------------------- EGLM-------------------------------------- 

       psitauylax(i,j,1,l)  =   psiint_y_L 

       phitauylax(i,j,1,l)  =   phiint_y_L 

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSyxlax(i,j,1,l)   =   FSyxint_L
       FSyylax(i,j,1,l)   =   FSyyint_L
       FSyzlax(i,j,1,l)   =   FSyzint_L


! -------------------------------------- EGLM-------------------------------------- 


       BySxlax(i,j,1,l)   = Bxint_y_L  

       BySylax(i,j,1,l)   = Bxint_y_L 

       BySzlax(i,j,1,l)   = Bxint_y_L  

! -------------------------------------- EGLM-------------------------------------- 
       

    else if ( sy1 .le. 0.d0 .and. 0.d0 .le. Vystr ) then



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

    else if (Vystr .le. 0.d0 .and. 0.d0 .le. sy2) then


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

    else if ( sy2 .lt. 0.d0 ) then


! Conserved current HLL flux
!_________________________________________________________________________


       Jylax(i,j,1,l)     =   Jyint_R 

! Conserved Mass HLL flux
!_________________________________________________________________________


       FDylax(i,j,1,l)    =   FDyint_R

! Conserved Energy HLL flux
!_________________________________________________________________________


       Ftauylax(i,j,1,l)  =   Ftauyint_R 

! -------------------------------------- EGLM-------------------------------------- 

       psitauylax(i,j,1,l)  =   psiint_y_R 

       phitauylax(i,j,1,l)  =   phiint_y_R 

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


       FSyxlax(i,j,1,l)   =   FSyxint_R
       FSyylax(i,j,1,l)   =   FSyyint_R
       FSyzlax(i,j,1,l)   =   FSyzint_R


! -------------------------------------- EGLM-------------------------------------- 


       BySxlax(i,j,1,l)   = Bxint_y_R  

       BySylax(i,j,1,l)   = Bxint_y_R 

       BySzlax(i,j,1,l)   = Bxint_y_R  

! -------------------------------------- EGLM-------------------------------------- 


    else 

       print*, "Revise Vystar velocity it could be imaginary subroutine hllc_flow "
       stop

    end if


!$OMP END ORDERED
    

        end do
     end do

!$OMP END DO

     
!$OMP END PARALLEL


   else

      write(*,*) "STOP: subroutine hllflux"
      write(*,*) "This dimension is not implemented yet"
      stop

   end if


  end subroutine hlle_flow_mpx
