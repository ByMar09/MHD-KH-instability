  !     *******************************************************************
  !     Subroutine which calculed  HLL flux flows 
  !     *******************************************************************

  subroutine hll_flow_mpx

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

       !---------------------------------
       ! Left - Right Conserved variables
       !---------------------------------

!//////////////////////////// PSI ///////////////////////////////


!$OMP   PARALLEL        
!$OMP DO 

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

            else if (REC_PRIM == 1) then

               u(i,21,1) = Vxint(i,1,1,l)
               u(i,22,1) = Vyint(i,1,1,l)
               u(i,23,1) = Vzint(i,1,1,l)

            else

               write(*,*) "STOP: subroutine hllflux"
               write(*,*) "REC_PRIM parameter is not valid"
               stop

            end if


            u(i,24,1) = p  (i,1,1)
            u(i,25,1) = rho(i,1,1)


!////////////////////////////  FLUXES    ///////////////////////////////
!*********** Condicional para reconstruir flujos usando mp5 o valores a left y right 

!!$            u(i,26,1) = Jx(i,1,1) !Jxint(i,1,1,l)
!!$            u(i,27,1) = Jy(i,1,1) !Jyint(i,1,1,l)1
!!$            u(i,28,1) = Jz(i,1,1) !Jzint(i,1,1,l)

            u(i,26,1) = Jxint(i,1,1,l)
            u(i,27,1) = Jyint(i,1,1,l)
            u(i,28,1) = Jzint(i,1,1,l)

            
!        print*, i, l, Jxint(i,1,1,l)

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
!$OMP END DO
!$OMP END PARALLEL
         
         
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

!$OMP PARALLEL        
!$OMP DO 
         do i=-1,imax+1

            psiint_mp5_L (i,1) = up(i  , 1,1) 
            psiint_mp5_R (i,1) = um(i+1, 1,1)

            phiint_mp5_L (i,1) = up(i  , 2,1) 
            phiint_mp5_R (i,1) = um(i+1, 2,1)

            Exint_mp5_L  (i,1) = up(i  , 3,1) 
            Exint_mp5_R  (i,1) = um(i+1, 3,1)
            Eyint_mp5_L  (i,1) = up(i  , 4,1) 
            Eyint_mp5_R  (i,1) = um(i+1, 4,1)
            Ezint_mp5_L  (i,1) = up(i  , 5,1) 
            Ezint_mp5_R  (i,1) = um(i+1, 5,1)

            Exint_mp5_L_m(i,1) = up(i  , 6,1) 
            Exint_mp5_R_m(i,1) = um(i+1, 6,1)
            Eyint_mp5_L_m(i,1) = up(i  , 7,1) 
            Eyint_mp5_R_m(i,1) = um(i+1, 7,1)
            Ezint_mp5_L_m(i,1) = up(i  , 8,1) 
            Ezint_mp5_R_m(i,1) = um(i+1, 8,1)

            Bxint_mp5_L  (i,1) = up(i  , 9,1) 
            Bxint_mp5_R  (i,1) = um(i+1, 9,1)
            Byint_mp5_L  (i,1) = up(i  ,10,1) 
            Byint_mp5_R  (i,1) = um(i+1,10,1)
            Bzint_mp5_L  (i,1) = up(i  ,11,1) 
            Bzint_mp5_R  (i,1) = um(i+1,11,1)

            Bxint_mp5_L_m(i,1) = up(i  ,12,1) 
            Bxint_mp5_R_m(i,1) = um(i+1,12,1)
            Byint_mp5_L_m(i,1) = up(i  ,13,1) 
            Byint_mp5_R_m(i,1) = um(i+1,13,1)
            Bzint_mp5_L_m(i,1) = up(i  ,14,1) 
            Bzint_mp5_R_m(i,1) = um(i+1,14,1)

            qint_mp5_L   (i,1) = up(i  ,15,1) 
            qint_mp5_R   (i,1) = um(i+1,15,1)
            
            DDint_mp5_L  (i,1) = up(i  ,16,1) 
            DDint_mp5_R  (i,1) = um(i+1,16,1)

            tauint_mp5_L (i,1) = up(i  ,17,1) 
            tauint_mp5_R (i,1) = um(i+1,17,1)

            Sxint_mp5_L  (i,1) = up(i  ,18,1) 
            Sxint_mp5_R  (i,1) = um(i+1,18,1)
            Syint_mp5_L  (i,1) = up(i  ,19,1) 
            Syint_mp5_R  (i,1) = um(i+1,19,1)
            Szint_mp5_L  (i,1) = up(i  ,20,1) 
            Szint_mp5_R  (i,1) = um(i+1,20,1)

            Vx_mp5_L     (i,1) = up(i  ,21,1) 
            Vx_mp5_R     (i,1) = um(i+1,21,1)
            Vy_mp5_L     (i,1) = up(i  ,22,1) 
            Vy_mp5_R     (i,1) = um(i+1,22,1)
            Vz_mp5_L     (i,1) = up(i  ,23,1) 
            Vz_mp5_R     (i,1) = um(i+1,23,1)

            p_mp5_L      (i,1) = up(i  ,24,1) 
            p_mp5_R      (i,1) = um(i+1,24,1)

            rho_mp5_L    (i,1) = up(i  ,25,1) 
            rho_mp5_R    (i,1) = um(i+1,25,1)

            Jxint_mp5_L  (i,1) = up(i  ,26,1) 
            Jxint_mp5_R  (i,1) = um(i+1,26,1) 
            Jyint_mp5_L  (i,1) = up(i  ,27,1) 
            Jyint_mp5_R  (i,1) = um(i+1,27,1)
            Jzint_mp5_L  (i,1) = up(i  ,28,1) 
            Jzint_mp5_R  (i,1) = um(i+1,28,1)


!*********** Condicional para reconstruir flujos usando mp5 o valores a left y right
            
            if ( flux_rec_mp5 == 1 ) then
               
               FDxint_mp5_L  (i,1) = up(i  ,29,1) 
               FDxint_mp5_R  (i,1) = um(i+1,29,1)
               FDyint_mp5_L  (i,1) = up(i  ,30,1) 
               FDyint_mp5_R  (i,1) = um(i+1,30,1)
               FDzint_mp5_L  (i,1) = up(i  ,31,1) 
               FDzint_mp5_R  (i,1) = um(i+1,31,1)

               Ftauxint_mp5_L(i,1) = up(i  ,32,1) 
               Ftauxint_mp5_R(i,1) = um(i+1,32,1)
               Ftauyint_mp5_L(i,1) = up(i  ,33,1) 
               Ftauyint_mp5_R(i,1) = um(i+1,33,1)
               Ftauzint_mp5_L(i,1) = up(i  ,34,1) 
               Ftauzint_mp5_R(i,1) = um(i+1,34,1)

               FSxxint_mp5_L (i,1) = up(i  ,35,1) 
               FSxxint_mp5_R (i,1) = um(i+1,35,1)
               FSxyint_mp5_L (i,1) = up(i  ,36,1) 
               FSxyint_mp5_R (i,1) = um(i+1,36,1)
               FSxzint_mp5_L (i,1) = up(i  ,37,1) 
               FSxzint_mp5_R (i,1) = um(i+1,37,1)

               FSyxint_mp5_L (i,1) = up(i  ,38,1) 
               FSyxint_mp5_R (i,1) = um(i+1,38,1)
               FSyyint_mp5_L (i,1) = up(i  ,39,1) 
               FSyyint_mp5_R (i,1) = um(i+1,39,1)
               FSyzint_mp5_L (i,1) = up(i  ,40,1) 
               FSyzint_mp5_R (i,1) = um(i+1,40,1)

               FSzxint_mp5_L (i,1) = up(i  ,41,1) 
               FSzxint_mp5_R (i,1) = um(i+1,41,1)
               FSzyint_mp5_L (i,1) = up(i  ,42,1) 
               FSzyint_mp5_R (i,1) = um(i+1,42,1)
               FSzzint_mp5_L (i,1) = up(i  ,43,1) 
               FSzzint_mp5_R (i,1) = um(i+1,43,1)

               
            end if
               
         end do

!$OMP END DO
!$OMP END PARALLEL         

!//////////////////////////// LORENTZ FACTOR ///////////////////////////////

!$OMP PARALLEL         
!$OMP DO 

         do i=-1,imax+1
 

             V2mp5 = Vx_mp5_L(i,1)**2 +  Vy_mp5_L(i,1)**2 + Vz_mp5_L(i,1)**2

            if (V2mp5 .ge. 1.d0) then

             print*, "super luminal left velocity in mp5 soubroutine, V2 =", V2mp5, "i=", i
!!$             stop

            end if

            V2mp5 = Vx_mp5_R(i,1)**2 +  Vy_mp5_R(i,1)**2 + Vz_mp5_R(i,1)**2

            if (V2mp5 .ge. 1.d0) then

             print*, "super luminal right velocity in mp5 soubroutine, V2 =", V2mp5, "i=", i
!!$             stop

            end if

             W_mp5_L(i,1) = 1.d0/sqrt(1.d0 - Vx_mp5_L(i,1)**2 -  Vy_mp5_L(i,1)**2 - Vz_mp5_L(i,1)**2)
             W_mp5_R(i,1) = 1.d0/sqrt(1.d0 - Vx_mp5_R(i,1)**2 -  Vy_mp5_R(i,1)**2 - Vz_mp5_R(i,1)**2)

 
          end do
          
!$OMP END DO
!$OMP END PARALLEL
          
!//////////////////////////// SQUARE EM FIELDS ///////////////////////////////


!$OMP PARALLEL          
!$OMP DO 

         do i=-1,imax+1
 

      E2B2int_mp5_L(i,1) = Exint_mp5_L(i,1)**2 + Eyint_mp5_L(i,1)**2 + Ezint_mp5_L(i,1)**2 + &
                           Bxint_mp5_L(i,1)**2 + Byint_mp5_L(i,1)**2 + Bzint_mp5_L(i,1)**2

      E2B2int_mp5_R(i,1) = Exint_mp5_R(i,1)**2 + Eyint_mp5_R(i,1)**2 + Ezint_mp5_R(i,1)**2 + &
                           Bxint_mp5_R(i,1)**2 + Byint_mp5_R(i,1)**2 + Bzint_mp5_R(i,1)**2


      Edotv_mp5_L(i,1) = Exint_mp5_L(i,1) * Vx_mp5_L(i,1) + Eyint_mp5_L(i,1) * Vy_mp5_L(i,1) + Ezint_mp5_L(i,1) * Vz_mp5_L(i,1)
      Edotv_mp5_R(i,1) = Exint_mp5_R(i,1) * Vx_mp5_R(i,1) + Eyint_mp5_R(i,1) * Vy_mp5_R(i,1) + Ezint_mp5_R(i,1) * Vz_mp5_R(i,1)

      vrotB_x_mp5_L(i,1) = (Bzint_mp5_L(i,1) * Vy_mp5_L(i,1) - Byint_mp5_L(i,1) * Vz_mp5_L(i,1))
      vrotB_y_mp5_L(i,1) = (Bxint_mp5_L(i,1) * Vz_mp5_L(i,1) - Bzint_mp5_L(i,1) * Vx_mp5_L(i,1))
      vrotB_z_mp5_L(i,1) = (Byint_mp5_L(i,1) * Vx_mp5_L(i,1) - Bxint_mp5_L(i,1) * Vy_mp5_L(i,1))

      vrotB_x_mp5_R(i,1) = (Bzint_mp5_R(i,1) * Vy_mp5_R(i,1) - Byint_mp5_R(i,1) * Vz_mp5_R(i,1))
      vrotB_y_mp5_R(i,1) = (Bxint_mp5_R(i,1) * Vz_mp5_R(i,1) - Bzint_mp5_R(i,1) * Vx_mp5_R(i,1))
      vrotB_z_mp5_R(i,1) = (Byint_mp5_R(i,1) * Vx_mp5_R(i,1) - Bxint_mp5_R(i,1) * Vy_mp5_R(i,1))


!      print*, i, vrotB_x_mp5_L(i,1), vrotB_x_mp5_R(i,1), Bzint_mp5_L(i,1), Vz_mp5_L(i,1), Byint_mp5_L(i,1), Vy_mp5_L(i,1)


 
         end do

!$OMP END DO
!$OMP END PARALLEL         

              if (SLC == 1) then
                 
                 sigma_L = sigma_0 * DDint_mp5_L(i,1)**gamma_slc 
                 sigma_R = sigma_0 * DDint_mp5_R(i,1)**gamma_slc 

              else if (SLC == 0) then

                 sigma_L = sigma
                 sigma_R = sigma

              end if


!SOURCE RECONSTRUCTION

      if (source_rec_mpx == 1) then

!$OMP PARALLEL
!$OMP DO 

         do i=-1,imax+1
 

    !     stiff source expression S = sigma * W * (E + V X B - (E.V)V) 

!!$    Exsour_L(i,1) = Jxint_mp5_L(i,1) - qint_mp5_L(i,1) *  Vx_mp5_L(i,1)
!!$    Exsour_R(i,1) = Jxint_mp5_R(i,1) - qint_mp5_R(i,1) *  Vx_mp5_R(i,1)
!!$
!!$    Eysour_L(i,1) = Jyint_mp5_L(i,1) - qint_mp5_L(i,1) *  Vy_mp5_L(i,1)
!!$    Eysour_R(i,1) = Jyint_mp5_R(i,1) - qint_mp5_R(i,1)  * Vy_mp5_R(i,1)
!!$
!!$    Ezsour_L(i,1) = Jzint_mp5_L(i,1) - qint_mp5_L(i,1) *  Vz_mp5_L(i,1)
!!$    Ezsour_R(i,1) = Jzint_mp5_R(i,1) - qint_mp5_R(i,1) *  Vz_mp5_R(i,1)
!!$
!!$
!!$
    Exsour_L(i,1) = sigma_L * W_mp5_L(i,1)  * ( Exint_mp5_L(i,1) +  vrotB_x_mp5_L(i,1)  - Edotv_mp5_L(i,1) * Vx_mp5_L(i,1) ) 
    Exsour_R(i,1) = sigma_R * W_mp5_R(i,1)  * ( Exint_mp5_R(i,1) +  vrotB_x_mp5_R(i,1)  - Edotv_mp5_R(i,1) * Vx_mp5_R(i,1) ) 

    Eysour_L(i,1) = sigma_L * W_mp5_L(i,1)  * ( Eyint_mp5_L(i,1) +  vrotB_y_mp5_L(i,1)  - Edotv_mp5_L(i,1) * Vy_mp5_L(i,1) ) 
    Eysour_R(i,1) = sigma_R * W_mp5_R(i,1)  * ( Eyint_mp5_R(i,1) +  vrotB_y_mp5_R(i,1)  - Edotv_mp5_R(i,1) * Vy_mp5_R(i,1) )

    Ezsour_L(i,1) = sigma_L * W_mp5_L(i,1)  * ( Ezint_mp5_L(i,1) +  vrotB_z_mp5_L(i,1)  - Edotv_mp5_L(i,1) * Vz_mp5_L(i,1) )
    Ezsour_R(i,1) = sigma_R * W_mp5_R(i,1)  * ( Ezint_mp5_R(i,1) +  vrotB_z_mp5_R(i,1)  - Edotv_mp5_R(i,1) * Vz_mp5_R(i,1) ) 
    


!        print*, i, vrotB_x_mp5_L(i,1), vrotB_x_mp5_R(i,1), Jxint_mp5_L(i,1), Jxint_mp5_R(i,1)
 
         end do
 
!$OMP END DO
!$OMP END PARALLEL

      end if



              

!!$!$OMP PARALLEL
!!$!$OMP DO 
!!$
!!$         do i=-1,imax+1
!!$ 
!!$
!!$    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV
!!$
!!$    Jxint_mp5_L(i,1) = sigma_L * W_mp5_L(i,1)  * ( Exint_mp5_L(i,1) +  vrotB_x_mp5_L(i,1)  - Edotv_mp5_L(i,1) * Vx_mp5_L(i,1) )   &
!!$                           + qint_mp5_L(i,1) *  Vx_mp5_L(i,1)
!!$    Jxint_mp5_R(i,1) = sigma_R * W_mp5_R(i,1)  * ( Exint_mp5_R(i,1) +  vrotB_x_mp5_R(i,1)  - Edotv_mp5_R(i,1) * Vx_mp5_R(i,1) )   &
!!$                           + qint_mp5_R(i,1) *  Vx_mp5_R(i,1)
!!$
!!$
!!$    Jyint_mp5_L(i,1) = sigma_L * W_mp5_L(i,1)  * ( Eyint_mp5_L(i,1) +  vrotB_y_mp5_L(i,1)  - Edotv_mp5_L(i,1) * Vy_mp5_L(i,1) )   &
!!$                           + qint_mp5_L(i,1) *  Vy_mp5_L(i,1)
!!$    Jyint_mp5_R(i,1) = sigma_R * W_mp5_R(i,1)  * ( Eyint_mp5_R(i,1) +  vrotB_y_mp5_R(i,1)  - Edotv_mp5_R(i,1) * Vy_mp5_R(i,1) )   &
!!$                           +   qint_mp5_R(i,1)  *  Vy_mp5_R(i,1)
!!$
!!$    Jzint_mp5_L(i,1) = sigma_L * W_mp5_L(i,1)  * ( Ezint_mp5_L(i,1) +  vrotB_z_mp5_L(i,1)  - Edotv_mp5_L(i,1) * Vz_mp5_L(i,1) )   &
!!$                           + qint_mp5_L(i,1) *  Vz_mp5_L(i,1)
!!$    Jzint_mp5_R(i,1) = sigma_R * W_mp5_R(i,1)  * ( Ezint_mp5_R(i,1) +  vrotB_z_mp5_R(i,1)  - Edotv_mp5_R(i,1) * Vz_mp5_R(i,1) )   &
!!$                           + qint_mp5_R(i,1) *  Vz_mp5_R(i,1)
!!$
!!$
!!$!        print*, i, vrotB_x_mp5_L(i,1), vrotB_x_mp5_R(i,1), Jxint_mp5_L(i,1), Jxint_mp5_R(i,1)
!!$ 
!!$         end do
!!$ 
!!$!$OMP END DO
!!$!$OMP END PARALLEL
         

!//////////////////////////// ENTHALPY ///////////////////////////////

!$OMP PARALLEL
!$OMP DO   
         do i=-1,imax+1
 

    epsilon_mp5_L(i,1) = p_mp5_L(i,1) / ((gamma-1.d0) * rho_mp5_L(i,1))
    epsilon_mp5_R(i,1) = p_mp5_R(i,1) / ((gamma-1.d0) * rho_mp5_R(i,1))

    enthpy_mp5_L(i,1)  = rho_mp5_L(i,1) * ( 1.d0 + epsilon_mp5_L(i,1) ) + p_mp5_L(i,1) 
    enthpy_mp5_R(i,1)  = rho_mp5_R(i,1) * ( 1.d0 + epsilon_mp5_R(i,1) ) + p_mp5_R(i,1) 

 
         end do
!$OMP END DO
!$OMP END PARALLEL

!*********** ELSE del Condicional para reconstruir flujos usando mp5 o valores a left y right 

      if ( flux_rec_mp5 == 0 ) then
 
!$OMP PARALLEL      
!$OMP DO   
         do i=-1,imax+1
 

    !     Density flux (FD = \rho W V)

    FDxint_mp5_L(i,1) = DDint_mp5_L(i,1) * Vx_mp5_L(i,1)
    FDxint_mp5_R(i,1) = DDint_mp5_R(i,1) * Vx_mp5_R(i,1)

    FDyint_mp5_L(i,1) = DDint_mp5_L(i,1) * Vy_mp5_L(i,1)
    FDyint_mp5_R(i,1) = DDint_mp5_R(i,1) * Vy_mp5_R(i,1)

    FDzint_mp5_L(i,1) = DDint_mp5_L(i,1) * Vz_mp5_L(i,1)
    FDzint_mp5_R(i,1) = DDint_mp5_R(i,1) * Vz_mp5_R(i,1)


    !     Flujos conservados de energía

    Ftauxint_mp5_L(i,1) = (Bzint_mp5_L(i,1)  * Eyint_mp5_L(i,1) - Byint_mp5_L(i,1) * Ezint_mp5_L(i,1)) +   &
                  enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vx_mp5_L(i,1)   
    Ftauxint_mp5_R(i,1) = (Bzint_mp5_R(i,1)  * Eyint_mp5_R(i,1) - Byint_mp5_R(i,1) * Ezint_mp5_R(i,1)) +   &
                  enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vx_mp5_R(i,1)   

    Ftauyint_mp5_L(i,1) = (Bxint_mp5_L(i,1)  * Ezint_mp5_L(i,1) - Bzint_mp5_L(i,1) * Exint_mp5_L(i,1)) +   &
                  enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vy_mp5_L(i,1)
    Ftauyint_mp5_R(i,1) = (Bxint_mp5_R(i,1)  * Ezint_mp5_R(i,1) - Bzint_mp5_R(i,1) * Exint_mp5_R(i,1)) +   &
                  enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vy_mp5_R(i,1)

    Ftauzint_mp5_L(i,1) = (Byint_mp5_L(i,1)  * Exint_mp5_L(i,1) - Bxint_mp5_L(i,1) * Eyint_mp5_L(i,1)) +   &
                  enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vz_mp5_L(i,1)
    Ftauzint_mp5_R(i,1) = (Byint_mp5_R(i,1)  * Exint_mp5_R(i,1) - Bxint_mp5_R(i,1) * Eyint_mp5_R(i,1)) +   &
                  enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vz_mp5_R(i,1)

    !     Components of the flux momentum tensor 

    FSxxint_mp5_L(i,1) =  - Exint_mp5_L(i,1)**2 - Bxint_mp5_L(i,1)**2 + enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2 * Vx_mp5_L(i,1)**2 + &
                   0.5d0 * E2B2int_mp5_L(i,1) + p_mp5_L(i,1)  
    FSxxint_mp5_R(i,1) =  - Exint_mp5_R(i,1)**2 - Bxint_mp5_R(i,1)**2 + enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2 * Vx_mp5_R(i,1)**2 + &
                   0.5d0 * E2B2int_mp5_R(i,1) + p_mp5_R(i,1)  

    FSxyint_mp5_L(i,1) = -  Exint_mp5_L(i,1)  * Eyint_mp5_L(i,1) - Bxint_mp5_L(i,1) * Byint_mp5_L(i,1) +    &
                   enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vx_mp5_L(i,1)    * Vy_mp5_L(i,1)
    FSxyint_mp5_R(i,1) = -  Exint_mp5_R(i,1)  * Eyint_mp5_R(i,1) - Bxint_mp5_R(i,1) * Byint_mp5_R(i,1) +    &
                   enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vx_mp5_R(i,1)    * Vy_mp5_R(i,1)

    FSxzint_mp5_L(i,1) = -  Exint_mp5_L(i,1)  * Ezint_mp5_L(i,1) - Bxint_mp5_L(i,1) * Bzint_mp5_L(i,1) +    &
                   enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vx_mp5_L(i,1)    * Vz_mp5_L(i,1)
    FSxzint_mp5_R(i,1) = -  Exint_mp5_R(i,1)  * Ezint_mp5_R(i,1) - Bxint_mp5_R(i,1) * Bzint_mp5_R(i,1) +    &
                   enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vx_mp5_R(i,1)    * Vz_mp5_R(i,1)

    FSyxint_mp5_L(i,1) = -  Exint_mp5_L(i,1)  * Eyint_mp5_L(i,1) - Bxint_mp5_L(i,1) * Byint_mp5_L(i,1) +    &
                   enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vx_mp5_L(i,1)    * Vy_mp5_L(i,1)
    FSyxint_mp5_R(i,1) = -  Exint_mp5_R(i,1)  * Eyint_mp5_R(i,1) - Bxint_mp5_R(i,1) * Byint_mp5_R(i,1) +    &
                   enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vx_mp5_R(i,1)    * Vy_mp5_R(i,1)

    FSyyint_mp5_L(i,1) = -  Eyint_mp5_L(i,1)**2 - Byint_mp5_L(i,1)**2 + enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2 * Vy_mp5_L(i,1)**2 + &
                   0.5d0 * E2B2int_mp5_L(i,1)  + p_mp5_L(i,1)
    FSyyint_mp5_R(i,1) = -  Eyint_mp5_R(i,1)**2 - Byint_mp5_R(i,1)**2 + enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2 * Vy_mp5_R(i,1)**2 + &
                   0.5d0 * E2B2int_mp5_R(i,1)  + p_mp5_R(i,1)

    FSyzint_mp5_L(i,1) = -  Eyint_mp5_L(i,1)  * Ezint_mp5_L(i,1) - Byint_mp5_L(i,1) * Bzint_mp5_L(i,1) +    &
                   enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vy_mp5_L(i,1)    * Vz_mp5_L(i,1)
    FSyzint_mp5_R(i,1) = -  Eyint_mp5_R(i,1)  * Ezint_mp5_R(i,1) - Byint_mp5_R(i,1) * Bzint_mp5_R(i,1) +    &
                   enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vy_mp5_R(i,1)    * Vz_mp5_R(i,1)

    FSzxint_mp5_L(i,1) = -  Exint_mp5_L(i,1)  * Ezint_mp5_L(i,1) - Bxint_mp5_L(i,1) * Bzint_mp5_L(i,1) +    &
                   enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vx_mp5_L(i,1)    * Vz_mp5_L(i,1)       
    FSzxint_mp5_R(i,1) = -  Exint_mp5_R(i,1)  * Ezint_mp5_R(i,1) - Bxint_mp5_R(i,1) * Bzint_mp5_R(i,1) +    &
                   enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vx_mp5_R(i,1)    * Vz_mp5_R(i,1)             

    FSzyint_mp5_L(i,1) = -  Eyint_mp5_L(i,1)  * Ezint_mp5_L(i,1) - Byint_mp5_L(i,1) * Bzint_mp5_L(i,1) +    &
                   enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2  * Vy_mp5_L(i,1)    * Vz_mp5_L(i,1)
    FSzyint_mp5_R(i,1) = -  Eyint_mp5_R(i,1)  * Ezint_mp5_R(i,1) - Byint_mp5_R(i,1) * Bzint_mp5_R(i,1) +    &
                   enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2  * Vy_mp5_R(i,1)    * Vz_mp5_R(i,1)

    FSzzint_mp5_L(i,1) = -  Ezint_mp5_L(i,1)**2 - Bzint_mp5_L(i,1)**2 + enthpy_mp5_L(i,1) * W_mp5_L(i,1)**2 * Vz_mp5_L(i,1)**2 + &
                   0.5d0 * E2B2int_mp5_L(i,1) + p_mp5_L(i,1)
    FSzzint_mp5_R(i,1) = -  Ezint_mp5_R(i,1)**2 - Bzint_mp5_R(i,1)**2 + enthpy_mp5_R(i,1) * W_mp5_R(i,1)**2 * Vz_mp5_R(i,1)**2 + &
                   0.5d0 * E2B2int_mp5_R(i,1) + p_mp5_R(i,1)


 
         end do

!$OMP END DO
!$OMP END PARALLEL

      end if
      


!$OMP PARALLEL
!$OMP DO 

    do i=-1,imax+1


!////////////////////////////  HLL FLUXES    ///////////////////////////////



! Electric Gauss HLL flux
!_________________________________________________________________________

   Exlaxpsi(i,1,1,l) = (s2* Exint_mp5_L(i,1) - s1* Exint_mp5_R(i,1) + s1*s2* ( psiint_mp5_R(i,1) - psiint_mp5_L(i,1) ))/(s2 -s1) 
   Eylaxpsi(i,1,1,l) = (s2* Eyint_mp5_L(i,1) - s1* Eyint_mp5_R(i,1) + s1*s2* ( psiint_mp5_R(i,1) - psiint_mp5_L(i,1) ))/(s2 -s1) 
   Ezlaxpsi(i,1,1,l) = (s2* Ezint_mp5_L(i,1) - s1* Ezint_mp5_R(i,1) + s1*s2* ( psiint_mp5_R(i,1) - psiint_mp5_L(i,1) ))/(s2 -s1)

! Magnetic Gauss HLL flux
!_________________________________________________________________________


   Bxlaxphi(i,1,1,l) = (s2* Bxint_mp5_L(i,1) - s1* Bxint_mp5_R(i,1) + s1*s2* ( phiint_mp5_R(i,1) - phiint_mp5_L(i,1) ))/(s2 -s1)
   Bylaxphi(i,1,1,l) = (s2* Byint_mp5_L(i,1) - s1* Byint_mp5_R(i,1) + s1*s2* ( phiint_mp5_R(i,1) - phiint_mp5_L(i,1) ))/(s2 -s1)
   Bzlaxphi(i,1,1,l) = (s2* Bzint_mp5_L(i,1) - s1* Bzint_mp5_R(i,1) + s1*s2* ( phiint_mp5_R(i,1) - phiint_mp5_L(i,1) ))/(s2 -s1)

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

   philaxBx(i,1,1,l) = (s2* phiint_mp5_L (i,1) - s1* phiint_mp5_R (i,1) + s1*s2* ( Bxint_mp5_R(i,1) - Bxint_mp5_L(i,1) ))/(s2 -s1)
   EzlaxBx (i,1,1,l) = (s2* Ezint_mp5_L  (i,1) - s1* Ezint_mp5_R  (i,1) + s1*s2* ( Bxint_mp5_R(i,1) - Bxint_mp5_L(i,1) ))/(s2 -s1)
   EylaxBx (i,1,1,l) =-(s2* Eyint_mp5_L_m(i,1) - s1* Eyint_mp5_R_m(i,1) + s1*s2* ( Bxint_mp5_R(i,1) - Bxint_mp5_L(i,1) ))/(s2 -s1) !sign

! Flujos By

   EzlaxBy (i,1,1,l) =-(s2* Ezint_mp5_L_m(i,1) - s1* Ezint_mp5_R_m(i,1) + s1*s2* ( Byint_mp5_R(i,1) - Byint_mp5_L(i,1) ))/(s2 -s1) !sign
   philaxBy(i,1,1,l) = (s2* phiint_mp5_L (i,1) - s1* phiint_mp5_R (i,1) + s1*s2* ( Byint_mp5_R(i,1) - Byint_mp5_L(i,1) ))/(s2 -s1) 
   ExlaxBy (i,1,1,l) = (s2* Exint_mp5_L  (i,1) - s1* Exint_mp5_R  (i,1) + s1*s2* ( Byint_mp5_R(i,1) - Byint_mp5_L(i,1) ))/(s2 -s1) 

! Flujos Bz

   EylaxBz (i,1,1,l) = (s2* Eyint_mp5_L  (i,1) - s1* Eyint_mp5_R  (i,1) + s1*s2* ( Bzint_mp5_R(i,1) - Bzint_mp5_L(i,1) ))/(s2 -s1) 
   ExlaxBz (i,1,1,l) =-(s2* Exint_mp5_L_m(i,1) - s1* Exint_mp5_R_m(i,1) + s1*s2* ( Bzint_mp5_R(i,1) - Bzint_mp5_L(i,1) ))/(s2 -s1) !sign
   philaxBz(i,1,1,l) = (s2* phiint_mp5_L (i,1) - s1* phiint_mp5_R (i,1) + s1*s2* ( Bzint_mp5_R(i,1) - Bzint_mp5_L(i,1) ))/(s2 -s1)

! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

   psilaxEx(i,1,1,l) = (s2* psiint_mp5_L (i,1) - s1* psiint_mp5_R (i,1) + s1*s2* ( Exint_mp5_R(i,1) - Exint_mp5_L(i,1) ))/(s2 -s1)
   BzlaxEx(i,1,1,l)  =-(s2* Bzint_mp5_L_m(i,1) - s1* Bzint_mp5_R_m(i,1) + s1*s2* ( Exint_mp5_R(i,1) - Exint_mp5_L(i,1) ))/(s2 -s1) !sign
   BylaxEx(i,1,1,l)  = (s2* Byint_mp5_L  (i,1) - s1* Byint_mp5_R  (i,1) + s1*s2* ( Exint_mp5_R(i,1) - Exint_mp5_L(i,1) ))/(s2 -s1)

!!$   if ( l == lmax .and. h == hmax-1) then
!!$
!!$      open (unit = 500, file = "numerical_flux_ex.dat")
!!$
!!$      write (500,*)  - 0.5d0 * Longx + i * Delx, psilaxEx(i,1,1,l)
!!$
!!$   end if
   

!   print*, i,l, psilaxEx(i,1,1,l), psiint_mp5_L (i,1), psiint_mp5_R (i,1), Exint_mp5_R(i,1), Exint_mp5_L(i,1), s1, s2

! Flujos Ey


   BzlaxEy (i,1,1,l) = (s2* Bzint_mp5_L  (i,1) - s1* Bzint_mp5_R  (i,1) + s1*s2* ( Eyint_mp5_R(i,1) - Eyint_mp5_L(i,1) ))/(s2 -s1)
   psilaxEy(i,1,1,l) = (s2* psiint_mp5_L (i,1) - s1* psiint_mp5_R (i,1) + s1*s2* ( Eyint_mp5_R(i,1) - Eyint_mp5_L(i,1) ))/(s2 -s1)
   BxlaxEy (i,1,1,l) =-(s2* Bxint_mp5_L_m(i,1) - s1* Bxint_mp5_R_m(i,1) + s1*s2* ( Eyint_mp5_R(i,1) - Eyint_mp5_L(i,1) ))/(s2 -s1) !sign

! Flujos Ez


   BylaxEz (i,1,1,l) =-(s2* Byint_mp5_L_m(i,1) - s1* Byint_mp5_R_m(i,1) + s1*s2* ( Ezint_mp5_R(i,1) - Ezint_mp5_L(i,1) ))/(s2 -s1) !sign
   BxlaxEz (i,1,1,l) = (s2* Bxint_mp5_L  (i,1) - s1* Bxint_mp5_R  (i,1) + s1*s2* ( Ezint_mp5_R(i,1) - Ezint_mp5_L(i,1) ))/(s2 -s1)
   psilaxEz(i,1,1,l) = (s2* psiint_mp5_L (i,1) - s1* psiint_mp5_R (i,1) + s1*s2* ( Ezint_mp5_R(i,1) - Ezint_mp5_L(i,1) ))/(s2 -s1)


! Electric field sources
!_________________________________________________________________________
   

      if (source_rec_mpx == 1) then
   
   Exastsour(i,1,1,l) = 0.5d0 * (Exsour_L(i,1) + Exsour_R(i,1))

   Eyastsour(i,1,1,l) = 0.5d0 * (Eysour_L(i,1) + Eysour_R(i,1))

   Ezastsour(i,1,1,l) = 0.5d0 * (Ezsour_L(i,1) + Ezsour_R(i,1))


      end if

! Conserved current HLL flux
!_________________________________________________________________________


   Jxlax(i,1,1,l)    = (s2* Jxint_mp5_L(i,1) - s1* Jxint_mp5_R(i,1) + s1*s2* ( qint_mp5_R(i,1) - qint_mp5_L(i,1) ))/ (s2 -s1)
   Jylax(i,1,1,l)    = (s2* Jyint_mp5_L(i,1) - s1* Jyint_mp5_R(i,1) + s1*s2* ( qint_mp5_R(i,1) - qint_mp5_L(i,1) ))/ (s2 -s1)
   Jzlax(i,1,1,l)    = (s2* Jzint_mp5_L(i,1) - s1* Jzint_mp5_R(i,1) + s1*s2* ( qint_mp5_R(i,1) - qint_mp5_L(i,1) ))/ (s2 -s1)


! Conserved Mass HLL flux
!_________________________________________________________________________


   FDxlax(i,1,1,l)   = (s2* FDxint_mp5_L(i,1) - s1* FDxint_mp5_R(i,1) + s1*s2* ( DDint_mp5_R(i,1) - DDint_mp5_L(i,1) ))/ (s2 -s1)
   FDylax(i,1,1,l)   = (s2* FDyint_mp5_L(i,1) - s1* FDyint_mp5_R(i,1) + s1*s2* ( DDint_mp5_R(i,1) - DDint_mp5_L(i,1) ))/ (s2 -s1)
   FDzlax(i,1,1,l)   = (s2* FDzint_mp5_L(i,1) - s1* FDzint_mp5_R(i,1) + s1*s2* ( DDint_mp5_R(i,1) - DDint_mp5_L(i,1) ))/ (s2 -s1)

! Conserved Energy HLL flux
!_________________________________________________________________________


   Ftauxlax(i,1,1,l) = (s2* Ftauxint_mp5_L(i,1) - s1* Ftauxint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) -tauint_mp5_L(i,1) ))/(s2-s1)
   Ftauylax(i,1,1,l) = (s2* Ftauyint_mp5_L(i,1) - s1* Ftauyint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) -tauint_mp5_L(i,1) ))/(s2-s1)
   Ftauzlax(i,1,1,l) = (s2* Ftauzint_mp5_L(i,1) - s1* Ftauzint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) -tauint_mp5_L(i,1) ))/(s2-s1) 

! -------------------------------------- EGLM-------------------------------------- 

   psitauxlax(i,1,1,l)= (s2* psiint_mp5_L(i,1) - s1* psiint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) - tauint_mp5_L(i,1) ))/(s2-s1) 
   psitauylax(i,1,1,l)= (s2* psiint_mp5_L(i,1) - s1* psiint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) - tauint_mp5_L(i,1) ))/(s2-s1) 
   psitauzlax(i,1,1,l)= (s2* psiint_mp5_L(i,1) - s1* psiint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) - tauint_mp5_L(i,1) ))/(s2-s1)


   phitauxlax(i,1,1,l)= (s2* phiint_mp5_L(i,1) - s1* phiint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) - tauint_mp5_L(i,1) ))/(s2-s1)
   phitauylax(i,1,1,l)= (s2* phiint_mp5_L(i,1) - s1* phiint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) - tauint_mp5_L(i,1) ))/(s2-s1)
   phitauzlax(i,1,1,l)= (s2* phiint_mp5_L(i,1) - s1* phiint_mp5_R(i,1) + s1*s2* (tauint_mp5_R(i,1) - tauint_mp5_L(i,1) ))/(s2-s1)

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


   FSxxlax(i,1,1,l) = (s2* FSxxint_mp5_L(i,1) - s1* FSxxint_mp5_R(i,1) + s1*s2* (Sxint_mp5_R(i,1) - Sxint_mp5_L(i,1) )) / (s2 -s1)
   FSxylax(i,1,1,l) = (s2* FSxyint_mp5_L(i,1) - s1* FSxyint_mp5_R(i,1) + s1*s2* (Syint_mp5_R(i,1) - Syint_mp5_L(i,1) )) / (s2 -s1)
   FSxzlax(i,1,1,l) = (s2* FSxzint_mp5_L(i,1) - s1* FSxzint_mp5_R(i,1) + s1*s2* (Szint_mp5_R(i,1) - Szint_mp5_L(i,1) )) / (s2 -s1) 

 
   FSyxlax(i,1,1,l) = (s2* FSyxint_mp5_L(i,1) - s1* FSyxint_mp5_R(i,1) + s1*s2* (Sxint_mp5_R(i,1) - Sxint_mp5_L(i,1) ))/ (s2 -s1)  
   FSyylax(i,1,1,l) = (s2* FSyyint_mp5_L(i,1) - s1* FSyyint_mp5_R(i,1) + s1*s2* (Syint_mp5_R(i,1) - Syint_mp5_L(i,1) ))/ (s2 -s1) 
   FSyzlax(i,1,1,l) = (s2* FSyzint_mp5_L(i,1) - s1* FSyzint_mp5_R(i,1) + s1*s2* (Szint_mp5_R(i,1) - Szint_mp5_L(i,1) ))/ (s2 -s1)


   FSzxlax(i,1,1,l) = (s2* FSzxint_mp5_L(i,1) - s1* FSzxint_mp5_R(i,1) + s1*s2* (Sxint_mp5_R(i,1) - Sxint_mp5_L(i,1) ))/ (s2 -s1)
   FSzylax(i,1,1,l) = (s2* FSzyint_mp5_L(i,1) - s1* FSzyint_mp5_R(i,1) + s1*s2* (Syint_mp5_R(i,1) - Syint_mp5_L(i,1) ))/ (s2 -s1)
   FSzzlax(i,1,1,l) = (s2* FSzzint_mp5_L(i,1) - s1* FSzzint_mp5_R(i,1) + s1*s2* (Szint_mp5_R(i,1) - Szint_mp5_L(i,1) ))/ (s2 -s1) 

! -------------------------------------- EGLM-------------------------------------- 


   BxSxlax(i,1,1,l) = (s2* Bxint_mp5_L(i,1) - s1* Bxint_mp5_R(i,1) + s1*s2* (Sxint_mp5_R(i,1) - Sxint_mp5_L(i,1) ) ) / (s2 -s1) 
   BySxlax(i,1,1,l) = (s2* Byint_mp5_L(i,1) - s1* Byint_mp5_R(i,1) + s1*s2* (Sxint_mp5_R(i,1) - Sxint_mp5_L(i,1) ) ) / (s2 -s1)
   BzSxlax(i,1,1,l) = (s2* Bzint_mp5_L(i,1) - s1* Bzint_mp5_R(i,1) + s1*s2* (Sxint_mp5_R(i,1) - Sxint_mp5_L(i,1) ) ) / (s2 -s1)



   BxSylax(i,1,1,l) = (s2* Bxint_mp5_L(i,1) - s1* Bxint_mp5_R(i,1) + s1*s2* (Syint_mp5_R(i,1) - Syint_mp5_L(i,1) ) ) / (s2 -s1)
   BySylax(i,1,1,l) = (s2* Byint_mp5_L(i,1) - s1* Byint_mp5_R(i,1) + s1*s2* (Syint_mp5_R(i,1) - Syint_mp5_L(i,1) ) ) / (s2 -s1)
   BzSylax(i,1,1,l) = (s2* Bzint_mp5_L(i,1) - s1* Bzint_mp5_R(i,1) + s1*s2* (Syint_mp5_R(i,1) - Syint_mp5_L(i,1) ) ) / (s2 -s1)


   BxSzlax(i,1,1,l) = (s2* Bxint_mp5_L(i,1) - s1* Bxint_mp5_R(i,1) + s1*s2* (Szint_mp5_R(i,1) - Szint_mp5_L(i,1) ) ) / (s2 -s1)
   BySzlax(i,1,1,l) = (s2* Byint_mp5_L(i,1) - s1* Byint_mp5_R(i,1) + s1*s2* (Szint_mp5_R(i,1) - Szint_mp5_L(i,1) ) ) / (s2 -s1)
   BzSzlax(i,1,1,l) = (s2* Bzint_mp5_L(i,1) - s1* Bzint_mp5_R(i,1) + s1*s2* (Szint_mp5_R(i,1) - Szint_mp5_L(i,1) ) ) / (s2 -s1)

! -------------------------------------- EGLM-------------------------------------- 

 

      end do

!$OMP END DO
!$OMP END PARALLEL

    else if (DIM == 2) then

!************************************************************************************

       call boundary_conserved
       call boundary_electric
       call boundary_flux_conserved


!$OMP   PARALLEL DEFAULT(SHARED) &
!$OMP & PRIVATE(dim_x,dim_y,rho__,not__,V2mp5,V2mp5_x,V2mp5_y,sigma_L,sigma_R,test_rho_L,test_rho_R),                             &
!$OMP & PRIVATE(dm4jph, dm4jmh, vul, vav, vmd, vlc, vmin, vmax,V2mp5_x_L,V2mp5_x_R,V2mp5_y_L,V2mp5_y_R),                          &
!$OMP & PRIVATE(psi_x_hll, phi_x_hll, Ex_x_hll, Ey_x_hll, Ez_x_hll, Bx_x_hll, By_x_hll, Bz_x_hll, q_x_hll, DD_x_hll, tau_x_hll),  &
!$OMP & PRIVATE(Ex_x_hll_flux, Sx_x_hll, Sy_x_hll, Sz_x_hll, psi_x_hll_flux, phi_x_hll_flux, Bx_x_hll_flux, By_x_hll_flux),       &      
!$OMP & PRIVATE(Bz_x_hll_flux, Ey_x_hll_flux, Ez_x_hll_flux, Exsource_L, Exsource_R, Eysource_L, Eysource_R, Ezsource_L),         &
!$OMP & PRIVATE(Ezsource_R, q_x_hll_flux, DD_x_hll_flux, tau_x_hll_flux, Sx_x_hll_flux, Sy_x_hll_flux, Sz_x_hll_flux, psi_y_hll), &
!$OMP & PRIVATE(phi_y_hll, Ey_y_hll, Ez_y_hll, Ex_y_hll, By_y_hll, Bz_y_hll, Bx_y_hll, q_y_hll, DD_y_hll, tau_y_hll, Sy_y_hll),   &
!$OMP & PRIVATE(Sz_y_hll, Sx_y_hll, psi_y_hll_flux, phi_y_hll_flux, By_y_hll_flux, Bz_y_hll_flux, Bx_y_hll_flux, Ey_y_hll_flux),  &
!$OMP & PRIVATE(Ez_y_hll_flux, Ex_y_hll_flux, q_y_hll_flux, DD_y_hll_flux, tau_y_hll_flux, Sy_y_hll_flux, Sz_y_hll_flux),         &
!$OMP & PRIVATE(Sx_y_hll_flux, Exstr_x_R, Exstr_x_L, Eystr_x_R, Eystr_x_L, Ezstr_x_R, Ezstr_x_L, Bxstr_x_R, Bxstr_x_L),           &
!$OMP & PRIVATE(Bystr_x_R, Bystr_x_L, Bzstr_x_R, Bzstr_x_L, qstr_x_R, qstr_x_L, Exstr_y_R, Exstr_y_L, Eystr_y_R, Eystr_y_L),      &
!$OMP & PRIVATE(Ezstr_y_R, Ezstr_y_L, Bxstr_y_R, Bxstr_y_L, Bystr_y_R, Bystr_y_L, Bzstr_y_R, Bzstr_y_L, qstr_y_R, qstr_y_L),      &
!$OMP & PRIVATE(ErotB_str_x_x_L, ErotB_str_y_x_L, ErotB_str_z_x_L, ErotB_str_x_x_R, ErotB_str_y_x_R, ErotB_str_z_x_R),            &
!$OMP & PRIVATE(ErotB_str_x_y_L, ErotB_str_y_y_L, ErotB_str_z_y_L, ErotB_str_x_y_R, ErotB_str_y_y_R, ErotB_str_z_y_R, E2_perp_x), &
!$OMP & PRIVATE(B2_perp_x, E2B2_str_x, p_elec, Vxstr, a_hll_x, b_hll_x, c_hll_x, det_hll, q_root1, rt1_x, rt2_x, E2_perp_y),      &
!$OMP & PRIVATE(B2_perp_y, E2B2_str_y, Vystr, a_hll_y, b_hll_y, c_hll_y, rt1_y, rt2_y, pstr_x, pstr_y, psistr_x_R, psistr_x_L),   &
!$OMP & PRIVATE(phistr_x_R, phistr_x_L, Sxstr_x_R, Sxstr_x_L, Systr_x_R, Systr_x_L, Szstr_x_R, Szstr_x_L, taustr_x_R,taustr_x_L), &
!$OMP & PRIVATE(DDstr_x_R, DDstr_x_L, psistr_y_R, psistr_y_L, phistr_y_R, phistr_y_L, Systr_y_R, Systr_y_L, Szstr_y_R),           &
!$OMP & PRIVATE(Szstr_y_L, Sxstr_y_R, Sxstr_y_L, taustr_y_R, taustr_y_L, DDstr_y_R, DDstr_y_L, psistr_x_flux_R, psistr_x_flux_L), &
!$OMP & PRIVATE(phistr_x_flux_R, phistr_x_flux_L, Exstr_x_flux_R, Exstr_x_flux_L, Eystr_x_flux_R, Eystr_x_flux_L,Ezstr_x_flux_R), &
!$OMP & PRIVATE(Ezstr_x_flux_L, Bxstr_x_flux_R, Bxstr_x_flux_L, Bystr_x_flux_R, Bystr_x_flux_L, Bzstr_x_flux_R, Bzstr_x_flux_L),  &
!$OMP & PRIVATE(FDxstr_R, FDxstr_L, FDzstr_R, FDzstr_L, Jxstr_L, Jxstr_R, Jzstr_R, Jzstr_L, Ftauxstr_R, Ftauxstr_L, Ftauzstr_R),  &
!$OMP & PRIVATE(Ftauzstr_L, FSxxstr_R, FSxxstr_L, FSxystr_R, FSxystr_L, FSxzstr_R, FSxzstr_L, FSzxstr_R, FSzxstr_L, FSzystr_R),   &
!$OMP & PRIVATE(FSzystr_L, FSzzstr_R, FSzzstr_L, psistr_y_flux_R, psistr_y_flux_L, phistr_y_flux_R, phistr_y_flux_L),             &
!$OMP & PRIVATE(Exstr_y_flux_R, Exstr_y_flux_L, Eystr_y_flux_R, Eystr_y_flux_L, Ezstr_y_flux_R, Ezstr_y_flux_L, Bxstr_y_flux_R),  &
!$OMP & PRIVATE(Bxstr_y_flux_L, Bystr_y_flux_R, Bystr_y_flux_L, Bzstr_y_flux_R, Bzstr_y_flux_L, FDystr_R, FDystr_L, Jystr_L),     &
!$OMP & PRIVATE(Jystr_R, Ftauystr_R, Ftauystr_L, FSyxstr_R, FSyxstr_L, FSyystr_R, FSyystr_L, FSyzstr_R, FSyzstr_L)
  
     
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

               if (REC_4VECTOR == 0) then

                  u(-6:imax+6,21,imp) = Vx (-6:imax+6,imp,1)
                  u(-6:imax+6,22,imp) = Vy (-6:imax+6,imp,1)
                  u(-6:imax+6,23,imp) = Vz (-6:imax+6,imp,1)

                  u(-6:imax+6,24,imp) = p  (-6:imax+6,imp,1)
                  u(-6:imax+6,25,imp) = rho(-6:imax+6,imp,1)

               else if (REC_4VECTOR == 1) then

                  W(-6:imax+6,imp,1) = 1.d0 / sqrt( 1.d0 -      &
                                       Vx(-6:imax+6,imp,1)**2 + &
                                       Vy(-6:imax+6,imp,1)**2 + &
                                       Vz(-6:imax+6,imp,1)**2 )

                  u(-6:imax+6,21,imp) = W(-6:imax+6,imp,1) * Vx (-6:imax+6,imp,1)
                  u(-6:imax+6,22,imp) = W(-6:imax+6,imp,1) * Vy (-6:imax+6,imp,1)
                  u(-6:imax+6,23,imp) = W(-6:imax+6,imp,1) * Vz (-6:imax+6,imp,1)

                  u(-6:imax+6,24,imp) = p  (-6:imax+6,imp,1)
                  u(-6:imax+6,25,imp) = rho(-6:imax+6,imp,1)
                  

               else

                  print*, "This REC_4VECTOR parameter is not valid. STOP subroutine 16_hllc_flow_mpx.f95"
                  stop

               end if

            else if (REC_PRIM == 1) then


               if (REC_4VECTOR == 0) then

                  u(-6:imax+6,21,imp) = Vxint (-6:imax+6,imp,1,l)
                  u(-6:imax+6,22,imp) = Vyint (-6:imax+6,imp,1,l)
                  u(-6:imax+6,23,imp) = Vzint (-6:imax+6,imp,1,l)

                  u(-6:imax+6,24,imp) = pint  (-6:imax+6,imp,1,l)
                  u(-6:imax+6,25,imp) = rhoint(-6:imax+6,imp,1,l)

              
               else if (REC_4VECTOR == 1) then

                  Wint(-6:imax+6,imp,1,l) = 1.d0 / sqrt( 1.d0 -      &
                                       Vxint(-6:imax+6,imp,1,l)**2 + &
                                       Vyint(-6:imax+6,imp,1,l)**2 + &
                                       Vzint(-6:imax+6,imp,1,l)**2 )

                  u(-6:imax+6,21,imp) = Wint(-6:imax+6,imp,1,l) * Vxint (-6:imax+6,imp,1,l)
                  u(-6:imax+6,22,imp) = Wint(-6:imax+6,imp,1,l) * Vyint (-6:imax+6,imp,1,l)
                  u(-6:imax+6,23,imp) = Wint(-6:imax+6,imp,1,l) * Vzint (-6:imax+6,imp,1,l)

                  u(-6:imax+6,24,imp) = p  (-6:imax+6,imp,1)
                  u(-6:imax+6,25,imp) = rho(-6:imax+6,imp,1)
                  

               else

                  print*, "This REC_4VECTOR parameter is not valid. STOP subroutine 16_hllc_flow_mpx.f95"
                  stop

               end if

 
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


            if (REC_4VECTOR == 0) then

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


            else if (REC_4VECTOR == 1) then


               Ux_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,21,imp) 
               Ux_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,21,imp)
               Uy_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,22,imp) 
               Uy_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,22,imp)
               Uz_mp5_x_L (-1:imax+1,imp) = up(-1:imax+1,23,imp) 
               Uz_mp5_x_R (-1:imax+1,imp) = um( 0:imax+2,23,imp)

               p_mp5_x_L  (-1:imax+1,imp) = up(-1:imax+1,24,imp) 
               p_mp5_x_R  (-1:imax+1,imp) = um( 0:imax+2,24,imp)
               rho_mp5_x_L(-1:imax+1,imp) = up(-1:imax+1,25,imp) 
               rho_mp5_x_R(-1:imax+1,imp) = um( 0:imax+2,25,imp) 

            else

                  print*, "This REC_4VECTOR parameter is not valid. STOP subroutine 16_hllc_flow_mpx.f95"
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

               if (REC_4VECTOR == 0) then

                  u(-6:jmax+6,58,imp) = Vx (imp,-6:jmax+6,1)
                  u(-6:jmax+6,59,imp) = Vy (imp,-6:jmax+6,1)
                  u(-6:jmax+6,60,imp) = Vz (imp,-6:jmax+6,1)

                  u(-6:jmax+6,61,imp) = p  (imp,-6:jmax+6,1)
                  u(-6:jmax+6,62,imp) = rho(imp,-6:jmax+6,1)

               else if (REC_4VECTOR == 1) then

                  W(imp,-6:jmax+6,1) = 1.d0 / sqrt( 1.d0 -      &
                                       Vx(imp,-6:jmax+6,1)**2 + &
                                       Vy(imp,-6:jmax+6,1)**2 + &
                                       Vz(imp,-6:jmax+6,1)**2 )
                  
                  u(-6:jmax+6,58,imp) = W(imp,-6:jmax+6,1) * Vx (imp,-6:jmax+6,1)
                  u(-6:jmax+6,59,imp) = W(imp,-6:jmax+6,1) * Vy (imp,-6:jmax+6,1)
                  u(-6:jmax+6,60,imp) = W(imp,-6:jmax+6,1) * Vz (imp,-6:jmax+6,1)

                  u(-6:jmax+6,61,imp) = p  (imp,-6:jmax+6,1)
                  u(-6:jmax+6,62,imp) = rho(imp,-6:jmax+6,1)

                
               else

                  print*, "This REC_4VECTOR parameter is not valid. STOP subroutine 16_hllc_flow_mpx.f95"
                  stop

               end if

            else if (REC_PRIM == 1) then


               if (REC_4VECTOR == 0) then


                  u(-6:jmax+6,58,imp) = Vxint (imp,-6:jmax+6,1,l)
                  u(-6:jmax+6,59,imp) = Vyint (imp,-6:jmax+6,1,l)
                  u(-6:jmax+6,60,imp) = Vzint (imp,-6:jmax+6,1,l)

                  u(-6:jmax+6,61,imp) = pint  (imp,-6:jmax+6,1,l)
                  u(-6:jmax+6,62,imp) = rhoint(imp,-6:jmax+6,1,l)
               
             
               else if (REC_4VECTOR == 1) then

                  Wint(imp,-6:jmax+6,1,l) = 1.d0 / sqrt( 1.d0 -      &
                                       Vxint(imp,-6:jmax+6,1,l)**2 + &
                                       Vyint(imp,-6:jmax+6,1,l)**2 + &
                                       Vzint(imp,-6:jmax+6,1,l)**2 )
                  
                  u(-6:jmax+6,58,imp) = Wint(imp,-6:jmax+6,1,l) * Vxint (imp,-6:jmax+6,1,l)
                  u(-6:jmax+6,59,imp) = Wint(imp,-6:jmax+6,1,l) * Vyint (imp,-6:jmax+6,1,l)
                  u(-6:jmax+6,60,imp) = Wint(imp,-6:jmax+6,1,l) * Vzint (imp,-6:jmax+6,1,l)

                  u(-6:jmax+6,61,imp) = pint  (imp,-6:jmax+6,1,l)
                  u(-6:jmax+6,62,imp) = rhoint(imp,-6:jmax+6,1,l)

               else

                  print*, "This REC_4VECTOR parameter is not valid. STOP subroutine 16_hllc_flow_mpx.f95"
                  stop

               end if

 
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


            if (REC_4VECTOR == 0) then

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
                  
            else if (REC_4VECTOR == 1) then

               Ux_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,58,imp) 
               Ux_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,58,imp)
               Uy_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,59,imp) 
               Uy_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,59,imp)
               Uz_mp5_y_L (imp,-1:jmax+1) = up(-1:jmax+1,60,imp) 
               Uz_mp5_y_R (imp,-1:jmax+1) = um( 0:jmax+2,60,imp)

               p_mp5_y_L  (imp,-1:jmax+1) = up(-1:jmax+1,61,imp) 
               p_mp5_y_R  (imp,-1:jmax+1) = um( 0:jmax+2,61,imp)
               rho_mp5_y_L(imp,-1:jmax+1) = up(-1:jmax+1,62,imp) 
               rho_mp5_y_R(imp,-1:jmax+1) = um( 0:jmax+2,62,imp) 

            else

               print*, "This REC_4VECTOR parameter is not valid. STOP subroutine 16_hllc_flow_mpx.f95"
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

     
             
!$OMP DO          
 
          do j=-1,jmax
             do i=-1,imax




                !//////////////////////////// LORENTZ FACTOR ///////////////////////////////


     ! Fix-up over velocities


     if (REC_4VECTOR == 0) then


            eps1    = 1.d-6

!!$            if (Vx_mp5_x_L(i,j) .ge. 1.d0) Vx_mp5_x_L(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vx_mp5_x_L(i,j) ))
!!$            if (Vx_mp5_x_R(i,j) .ge. 1.d0) Vx_mp5_x_R(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vx_mp5_x_R(i,j) ))
!!$            if (Vy_mp5_x_L(i,j) .ge. 1.d0) Vy_mp5_x_L(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vy_mp5_x_L(i,j) ))
!!$            if (Vy_mp5_x_R(i,j) .ge. 1.d0) Vy_mp5_x_R(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vy_mp5_x_R(i,j) ))
!!$            if (Vz_mp5_x_L(i,j) .ge. 1.d0) Vz_mp5_x_L(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vz_mp5_x_L(i,j) ))
!!$            if (Vz_mp5_x_R(i,j) .ge. 1.d0) Vz_mp5_x_R(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vz_mp5_x_R(i,j) ))
!!$
!!$            if (Vx_mp5_y_L(i,j) .ge. 1.d0) Vx_mp5_y_L(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vx_mp5_y_L(i,j) ))
!!$            if (Vx_mp5_y_R(i,j) .ge. 1.d0) Vx_mp5_y_R(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vx_mp5_y_R(i,j) ))
!!$            if (Vy_mp5_y_L(i,j) .ge. 1.d0) Vy_mp5_y_L(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vy_mp5_y_L(i,j) ))
!!$            if (Vy_mp5_y_R(i,j) .ge. 1.d0) Vy_mp5_y_R(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vy_mp5_y_R(i,j) ))
!!$            if (Vz_mp5_y_L(i,j) .ge. 1.d0) Vz_mp5_y_L(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vz_mp5_y_L(i,j) ))
!!$            if (Vz_mp5_y_R(i,j) .ge. 1.d0) Vz_mp5_y_R(i,j) = max(-1.d0 + eps1 , min( 1.d0 - eps1, Vz_mp5_y_R(i,j) ))
            

   
            V2mp5_x_L = Vx_mp5_x_L(i,j)**2 +  Vy_mp5_x_L(i,j)**2 + Vz_mp5_x_L(i,j)**2
            V2mp5_x_R = Vx_mp5_x_R(i,j)**2 +  Vy_mp5_x_R(i,j)**2 + Vz_mp5_x_R(i,j)**2
            
            V2mp5_y_L = Vx_mp5_y_L(i,j)**2 +  Vy_mp5_y_L(i,j)**2 + Vz_mp5_y_L(i,j)**2
            V2mp5_y_R = Vx_mp5_y_R(i,j)**2 +  Vy_mp5_y_R(i,j)**2 + Vz_mp5_y_R(i,j)**2


            if ( V2mp5_x_L .ge. 1.d0 .or. V2mp5_x_R .ge. 1.d0) then

!               print*, "super luminal left velocity in mp5 soubroutine, V2_x =", V2mp5_x_L, "i=", i, "j=", j

!               V2mp5_x_L = max(-1.d0 + eps1 , min( 1.d0 - eps1, V2mp5_x_L  ))

               Vx_mp5_x_L(i,j) = Vx(i  ,j  ,1)
               Vx_mp5_x_R(i,j) = Vx(i-1,j  ,1)
               Vy_mp5_x_L(i,j) = Vy(i  ,j  ,1)
               Vy_mp5_x_R(i,j) = Vy(i-1,j  ,1)
               Vz_mp5_x_L(i,j) = Vz(i  ,j  ,1)
               Vz_mp5_x_R(i,j) = Vz(i-1,j  ,1)

               Vx_mp5_y_L(i,j) = Vx(i  ,j  ,1)
               Vx_mp5_y_R(i,j) = Vx(i  ,j-1,1)
               Vy_mp5_y_L(i,j) = Vy(i  ,j  ,1)
               Vy_mp5_y_R(i,j) = Vy(i  ,j-1,1)
               Vz_mp5_y_L(i,j) = Vz(i  ,j  ,1)
               Vz_mp5_y_R(i,j) = Vz(i  ,j-1,1)

               V2mp5_y_L = Vx_mp5_y_L(i,j)**2 +  Vy_mp5_y_L(i,j)**2 + Vz_mp5_y_L(i,j)**2
               V2mp5_y_R = Vx_mp5_y_R(i,j)**2 +  Vy_mp5_y_R(i,j)**2 + Vz_mp5_y_R(i,j)**2

               V2mp5_x_L = Vx_mp5_x_L(i,j)**2 +  Vy_mp5_x_L(i,j)**2 + Vz_mp5_x_L(i,j)**2
               V2mp5_x_R = Vx_mp5_x_R(i,j)**2 +  Vy_mp5_x_R(i,j)**2 + Vz_mp5_x_R(i,j)**2

               print*, "FIXED velocity mp5 soubroutine, V2_X =", V2mp5_x_L,V2mp5_x_R,V2mp5_y_L,V2mp5_y_R,"i=", i, "j=", j

            end if


            if (V2mp5_y_L .ge. 1.d0 .or. V2mp5_y_R .ge. 1.d0) then


!               print*, "super luminal left velocity in mp5 soubroutine, V2_Y =", V2mp5_y_L, "i=", i, "j=", j

!               V2mp5_y_L = max(-1.d0 + eps1 , min( 1.d0 - eps1, V2mp5_y_L  ))

               Vx_mp5_x_L(i,j) = Vx(i  ,j  ,1)
               Vx_mp5_x_R(i,j) = Vx(i-1,j  ,1)
               Vy_mp5_x_L(i,j) = Vy(i  ,j  ,1)
               Vy_mp5_x_R(i,j) = Vy(i-1,j  ,1)
               Vz_mp5_x_L(i,j) = Vz(i  ,j  ,1)
               Vz_mp5_x_R(i,j) = Vz(i-1,j  ,1)

               Vx_mp5_y_L(i,j) = Vx(i  ,j  ,1)
               Vx_mp5_y_R(i,j) = Vx(i  ,j-1,1)
               Vy_mp5_y_L(i,j) = Vy(i  ,j  ,1)
               Vy_mp5_y_R(i,j) = Vy(i  ,j-1,1)
               Vz_mp5_y_L(i,j) = Vz(i  ,j  ,1)
               Vz_mp5_y_R(i,j) = Vz(i  ,j-1,1)

               V2mp5_y_L = Vx_mp5_y_L(i,j)**2 +  Vy_mp5_y_L(i,j)**2 + Vz_mp5_y_L(i,j)**2
               V2mp5_y_R = Vx_mp5_y_R(i,j)**2 +  Vy_mp5_y_R(i,j)**2 + Vz_mp5_y_R(i,j)**2


               V2mp5_x_L = Vx_mp5_x_L(i,j)**2 +  Vy_mp5_x_L(i,j)**2 + Vz_mp5_x_L(i,j)**2
               V2mp5_x_R = Vx_mp5_x_R(i,j)**2 +  Vy_mp5_x_R(i,j)**2 + Vz_mp5_x_R(i,j)**2


               print*, "FIXED velocity mp5 soubroutine, V2_Y =", V2mp5_x_L,V2mp5_x_R,V2mp5_y_L,V2mp5_y_R,"i=", i, "j=", j

            end if

             W_mp5_x_L(i,j) = 1.d0/sqrt(1.d0 - V2mp5_x_L)
             W_mp5_x_R(i,j) = 1.d0/sqrt(1.d0 - V2mp5_x_R)

             W_mp5_y_L(i,j) = 1.d0/sqrt(1.d0 - V2mp5_y_L)
             W_mp5_y_R(i,j) = 1.d0/sqrt(1.d0 - V2mp5_y_R)


          else if (REC_4VECTOR == 1) then

             W_mp5_x_L(i,j) = sqrt( 1.d0 + Ux_mp5_x_L(i,j)**2 &
                                         + Uy_mp5_x_L(i,j)**2 &
                                         + Uz_mp5_x_L(i,j)**2 )


             Vx_mp5_x_L(i,j) = Ux_mp5_x_L(i,j) / W_mp5_x_L(i,j)
             Vy_mp5_x_L(i,j) = Uy_mp5_x_L(i,j) / W_mp5_x_L(i,j)
             Vy_mp5_x_L(i,j) = Uy_mp5_x_L(i,j) / W_mp5_x_L(i,j)


             W_mp5_x_R(i,j) = sqrt( 1.d0 + Ux_mp5_x_R(i,j)**2 &
                                         + Uy_mp5_x_R(i,j)**2 &
                                         + Uz_mp5_x_R(i,j)**2 )
             

             Vx_mp5_x_R(i,j) = Ux_mp5_x_R(i,j) / W_mp5_x_R(i,j)
             Vy_mp5_x_R(i,j) = Uy_mp5_x_R(i,j) / W_mp5_x_R(i,j)
             Vy_mp5_x_R(i,j) = Uy_mp5_x_R(i,j) / W_mp5_x_R(i,j)


             W_mp5_y_L(i,j) = sqrt( 1.d0 + Ux_mp5_y_L(i,j)**2 &
                                         + Uy_mp5_y_L(i,j)**2 &
                                         + Uz_mp5_y_L(i,j)**2 )


             Vx_mp5_y_L(i,j) = Ux_mp5_y_L(i,j) / W_mp5_y_L(i,j)
             Vy_mp5_y_L(i,j) = Uy_mp5_y_L(i,j) / W_mp5_y_L(i,j)
             Vy_mp5_y_L(i,j) = Uy_mp5_y_L(i,j) / W_mp5_y_L(i,j)


             W_mp5_y_R(i,j) = sqrt( 1.d0 + Ux_mp5_y_R(i,j)**2 &
                                         + Uy_mp5_y_R(i,j)**2 &
                                         + Uz_mp5_y_R(i,j)**2 )
             

             Vx_mp5_y_R(i,j) = Ux_mp5_y_R(i,j) / W_mp5_y_R(i,j)
             Vy_mp5_y_R(i,j) = Uy_mp5_y_R(i,j) / W_mp5_y_R(i,j)
             Vy_mp5_y_R(i,j) = Uy_mp5_y_R(i,j) / W_mp5_y_R(i,j)

             
                
          else

             print*, "This REC_4VECTOR parameter is not valid. STOP subroutine 16_hllc_flow_mpx.f95"
             stop

          end if



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
      
      print*, "parametro flux_rec_mp5 no valido "
      stop

   end if


!!$         do i=-1,imax+1
!!$            do j=-1,jmax+1
!!$
!!$           print*, Bxint_mp5_x_L(i,j)
!!$
!!$        end do
!!$     end do


!////////////////////////////  HLL FLUXES    ///////////////////////////////


   ! Electric Gauss HLL flux

!_________________________________________________________________________

   Exlaxpsi(i,j,1,l) = (s2* Exint_mp5_x_L(i,j)- s1* Exint_mp5_x_R(i,j) &
        + s1*s2* (psiint_mp5_x_R(i,j)- psiint_mp5_x_L(i,j) ))/(s2-s1) 
   Eylaxpsi(i,j,1,l) = (s2* Eyint_mp5_y_L(i,j)- s1* Eyint_mp5_y_R(i,j) &
        + s1*s2* (psiint_mp5_y_R(i,j)- psiint_mp5_y_L(i,j) ))/(s2-s1)  
   Ezlaxpsi(i,j,1,l) = (s2* Ezint_mp5_x_L(i,j)- s1* Ezint_mp5_x_R(i,j) &
        + s1*s2* (psiint_mp5_x_R(i,j)- psiint_mp5_x_L(i,j) ))/(s2-s1) 

! Magnetic Gauss HLL flux
!_________________________________________________________________________


   Bxlaxphi(i,j,1,l) = (s2* Bxint_mp5_x_L(i,j)- s1* Bxint_mp5_x_R(i,j) &
        + s1*s2* (phiint_mp5_x_R(i,j)- phiint_mp5_x_L(i,j) ))/(s2-s1)
   Bylaxphi(i,j,1,l) = (s2* Byint_mp5_y_L(i,j)- s1* Byint_mp5_y_R(i,j) &
        + s1*s2* (phiint_mp5_y_R(i,j)- phiint_mp5_y_L(i,j) ))/(s2-s1)
   Bzlaxphi(i,j,1,l) = (s2* Bzint_mp5_x_L(i,j)- s1* Bzint_mp5_x_R(i,j) &
        + s1*s2* (phiint_mp5_x_R(i,j)- phiint_mp5_x_L(i,j) ))/(s2-s1)

 ! Faraday HLL flux
!_________________________________________________________________________

! Flujos Bx 

   philaxBx(i,j,1,l) = (s2* phiint_mp5_x_L (i,j)-s1* phiint_mp5_x_R (i,j) &
        +s1*s2* (Bxint_mp5_x_R(i,j)- Bxint_mp5_x_L(i,j) ))/(s2-s1)
   EzlaxBx (i,j,1,l) = (s2* Ezint_mp5_y_L  (i,j)-s1* Ezint_mp5_y_R  (i,j) &
        +s1*s2* (Bxint_mp5_y_R(i,j)- Bxint_mp5_y_L(i,j) ))/(s2-s1)
   EylaxBx (i,j,1,l) =-(s2* Eyint_mp5_x_L_m(i,j)-s1* Eyint_mp5_x_R_m(i,j) &
        +s1*s2* (Bxint_mp5_x_R(i,j)- Bxint_mp5_x_L(i,j) ))/(s2-s1) !sign

! Flujos By

   EzlaxBy (i,j,1,l) =-(s2* Ezint_mp5_x_L_m(i,j)-s1* Ezint_mp5_x_R_m(i,j) &
        +s1*s2* (Byint_mp5_x_R(i,j)- Byint_mp5_x_L(i,j) ))/(s2-s1) !sign
   philaxBy(i,j,1,l) = (s2* phiint_mp5_y_L (i,j)-s1* phiint_mp5_y_R (i,j) &
        +s1*s2* (Byint_mp5_y_R(i,j)- Byint_mp5_y_L(i,j) ))/(s2-s1)
   ExlaxBy (i,j,1,l) = (s2* Exint_mp5_x_L  (i,j)-s1* Exint_mp5_x_R  (i,j) &
        +s1*s2* (Byint_mp5_x_R(i,j)- Byint_mp5_x_L(i,j) ))/(s2-s1)

! Flujos Bz

   EylaxBz (i,j,1,l) = (s2* Eyint_mp5_x_L  (i,j)-s1* Eyint_mp5_x_R  (i,j) &
        +s1*s2* (Bzint_mp5_x_R(i,j)- Bzint_mp5_x_L(i,j) ))/(s2-s1)
   ExlaxBz (i,j,1,l) =-(s2* Exint_mp5_y_L_m(i,j)-s1* Exint_mp5_y_R_m(i,j) &
        +s1*s2* (Bzint_mp5_y_R(i,j)- Bzint_mp5_y_L(i,j) ))/(s2-s1) !sign
   philaxBz(i,j,1,l) = (s2* phiint_mp5_x_L (i,j)-s1* phiint_mp5_x_R (i,j) &
        +s1*s2* (Bzint_mp5_x_R(i,j)- Bzint_mp5_x_L(i,j) ))/(s2-s1)
   
! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

   psilaxEx(i,j,1,l) = (s2* psiint_mp5_x_L (i,j)-s1* psiint_mp5_x_R (i,j) &
        +s1*s2* (Exint_mp5_x_R(i,j)- Exint_mp5_x_L(i,j) ))/(s2-s1)
   BzlaxEx(i,j,1,l)  =-(s2* Bzint_mp5_y_L_m(i,j)-s1* Bzint_mp5_y_R_m(i,j) &
        +s1*s2* (Exint_mp5_y_R(i,j)- Exint_mp5_y_L(i,j) ))/(s2-s1) !sign
   BylaxEx(i,j,1,l)  = (s2* Byint_mp5_x_L  (i,j)-s1* Byint_mp5_x_R  (i,j) &
        +s1*s2* (Exint_mp5_x_R(i,j)- Exint_mp5_x_L(i,j) ))/(s2-s1)

! Flujos Ey


   BzlaxEy (i,j,1,l) = (s2* Bzint_mp5_x_L  (i,j)-s1* Bzint_mp5_x_R  (i,j) &
        +s1*s2* (Eyint_mp5_x_R(i,j)- Eyint_mp5_x_L(i,j) ))/(s2-s1)
   psilaxEy(i,j,1,l) = (s2* psiint_mp5_y_L (i,j)-s1* psiint_mp5_y_R (i,j) &
        +s1*s2* (Eyint_mp5_y_R(i,j)- Eyint_mp5_y_L(i,j) ))/(s2-s1)
   BxlaxEy (i,j,1,l) =-(s2* Bxint_mp5_x_L_m(i,j)-s1* Bxint_mp5_x_R_m(i,j) &
        +s1*s2* (Eyint_mp5_x_R(i,j)- Eyint_mp5_x_L(i,j) ))/(s2-s1) !sign

! Flujos Ez


   BylaxEz (i,j,1,l) =-(s2* Byint_mp5_x_L_m(i,j)-s1* Byint_mp5_x_R_m(i,j) &
        +s1*s2* (Ezint_mp5_x_R(i,j)- Ezint_mp5_x_L(i,j) ))/(s2-s1) !sign
   BxlaxEz (i,j,1,l) = (s2* Bxint_mp5_y_L  (i,j)-s1* Bxint_mp5_y_R  (i,j) &
        +s1*s2* (Ezint_mp5_y_R(i,j)- Ezint_mp5_y_L(i,j) ))/(s2-s1)
   psilaxEz(i,j,1,l) = (s2* psiint_mp5_x_L (i,j)-s1* psiint_mp5_x_R (i,j) &
        +s1*s2* (Ezint_mp5_x_R(i,j)- Ezint_mp5_x_L(i,j) ))/(s2-s1)


! Electric field sources
!_________________________________________________________________________
   

      if (source_rec_mpx == 1) then
   
   Exastsour(i,j,1,l) = 0.5d0 * (Exsour_L(i,j) + Exsour_R(i,j))

   Eyastsour(i,j,1,l) = 0.5d0 * (Eysour_L(i,j) + Eysour_R(i,j))

   Ezastsour(i,j,1,l) = 0.5d0 * (Ezsour_L(i,j) + Ezsour_R(i,j))


      end if   

! Conserved current HLL flux
!_________________________________________________________________________


   Jxlax(i,j,1,l)   =   (s2* Jxint_mp5_L(i,j) -s1* Jxint_mp5_R(i,j) &
        +s1*s2* ( qint_mp5_x_R(i,j) - qint_mp5_x_L(i,j) ) )/(s2-s1)
   Jylax(i,j,1,l)   =   (s2* Jyint_mp5_L(i,j) -s1* Jyint_mp5_R(i,j) &
        +s1*s2* ( qint_mp5_y_R(i,j) - qint_mp5_y_L(i,j) ) )/(s2-s1)
   Jzlax(i,j,1,l)   =   (s2* Jzint_mp5_L(i,j) -s1* Jzint_mp5_R(i,j) &
        +s1*s2* ( qint_mp5_x_R(i,j) - qint_mp5_x_L(i,j) ) )/(s2-s1)

! Conserved Mass HLL flux
!_________________________________________________________________________


   FDxlax(i,j,1,l)  =  (s2* FDxint_mp5_L(i,j) -s1* FDxint_mp5_R(i,j) &
        +s1*s2* ( DDint_mp5_x_R(i,j) - DDint_mp5_x_L(i,j) ) )/(s2-s1)
   FDylax(i,j,1,l)  =  (s2* FDyint_mp5_L(i,j) -s1* FDyint_mp5_R(i,j) &
        +s1*s2* ( DDint_mp5_y_R(i,j) - DDint_mp5_y_L(i,j) ) )/(s2-s1)
   FDzlax(i,j,1,l)  =  (s2* FDzint_mp5_L(i,j) -s1* FDzint_mp5_R(i,j) &
        +s1*s2* ( DDint_mp5_x_R(i,j) - DDint_mp5_x_L(i,j) ) )/(s2-s1)

! Conserved Energy HLL flux
!_________________________________________________________________________


   Ftauxlax(i,j,1,l) = (s2* Ftauxint_mp5_L(i,j)-s1* Ftauxint_mp5_R(i,j) &
        +s1*s2* (tauint_mp5_x_R(i,j)- tauint_mp5_x_L(i,j)) )/(s2-s1)
   Ftauylax(i,j,1,l) = (s2* Ftauyint_mp5_L(i,j)-s1* Ftauyint_mp5_R(i,j) &
        +s1*s2* (tauint_mp5_y_R(i,j)- tauint_mp5_y_L(i,j)) )/(s2-s1)
   Ftauzlax(i,j,1,l) = (s2* Ftauzint_mp5_L(i,j)-s1* Ftauzint_mp5_R(i,j) &
        +s1*s2* (tauint_mp5_x_R(i,j)- tauint_mp5_x_L(i,j)) )/(s2-s1)

! -------------------------------------- EGLM-------------------------------------- 

   psitauxlax(i,j,1,l)=(s2* psiint_mp5_x_L(i,j)-s1* psiint_mp5_x_R(i,j) &
        +s1*s2* (tauint_mp5_x_R(i,j)- tauint_mp5_x_L(i,j) ))/(s2-s1)
   psitauylax(i,j,1,l)=(s2* psiint_mp5_y_L(i,j)-s1* psiint_mp5_y_R(i,j) &
        +s1*s2* (tauint_mp5_y_R(i,j)- tauint_mp5_y_L(i,j) ))/(s2-s1)
   psitauzlax(i,j,1,l)=(s2* psiint_mp5_x_L(i,j)-s1* psiint_mp5_x_R(i,j) &
        +s1*s2* (tauint_mp5_x_R(i,j)- tauint_mp5_x_L(i,j) ))/(s2-s1)

   phitauxlax(i,j,1,l)=(s2* phiint_mp5_x_L(i,j)-s1* phiint_mp5_x_R(i,j) &
        +s1*s2* (tauint_mp5_x_R(i,j)- tauint_mp5_x_L(i,j) ))/(s2-s1)
   phitauylax(i,j,1,l)=(s2* phiint_mp5_y_L(i,j)-s1* phiint_mp5_y_R(i,j) &
        +s1*s2* (tauint_mp5_y_R(i,j)- tauint_mp5_y_L(i,j) ))/(s2-s1)
   phitauzlax(i,j,1,l)=(s2* phiint_mp5_x_L(i,j)-s1* phiint_mp5_x_R(i,j) &
        +s1*s2* (tauint_mp5_x_R(i,j)- tauint_mp5_x_L(i,j) ))/(s2-s1)

! -------------------------------------- EGLM-------------------------------------- 

! Conserved Momentum HLL flux tensor
!_________________________________________________________________________


   FSxxlax(i,j,1,l) = (s2* FSxxint_mp5_L(i,j)-s1* FSxxint_mp5_R(i,j) &
        +s1*s2* (Sxint_mp5_x_R(i,j) - Sxint_mp5_x_L(i,j) ) ) /(s2-s1)
   FSxylax(i,j,1,l) = (s2* FSxyint_mp5_L(i,j)-s1* FSxyint_mp5_R(i,j) &
        +s1*s2* (Syint_mp5_x_R(i,j) - Syint_mp5_x_L(i,j) ) ) /(s2-s1)
   FSxzlax(i,j,1,l) = (s2* FSxzint_mp5_L(i,j)-s1* FSxzint_mp5_R(i,j) &
        +s1*s2* (Szint_mp5_x_R(i,j) - Szint_mp5_x_L(i,j) ) ) /(s2-s1)

 
   FSyxlax(i,j,1,l) = (s2* FSyxint_mp5_L(i,j)-s1* FSyxint_mp5_R(i,j) &
        +s1*s2* (Sxint_mp5_y_R(i,j) - Sxint_mp5_y_L(i,j) ) ) /(s2-s1)
   FSyylax(i,j,1,l) = (s2* FSyyint_mp5_L(i,j)-s1* FSyyint_mp5_R(i,j) &
        +s1*s2* (Syint_mp5_y_R(i,j) - Syint_mp5_y_L(i,j) ) ) /(s2-s1)
   FSyzlax(i,j,1,l) = (s2* FSyzint_mp5_L(i,j)-s1* FSyzint_mp5_R(i,j) &
        +s1*s2* (Szint_mp5_y_R(i,j) - Szint_mp5_y_L(i,j) ) ) /(s2-s1)


   FSzxlax(i,j,1,l) = (s2* FSzxint_mp5_L(i,j)-s1* FSzxint_mp5_R(i,j) &
        +s1*s2* (Sxint_mp5_x_R(i,j) - Sxint_mp5_x_L(i,j) ) ) /(s2-s1)
   FSzylax(i,j,1,l) = (s2* FSzyint_mp5_L(i,j)-s1* FSzyint_mp5_R(i,j) &
        +s1*s2* (Syint_mp5_x_R(i,j) - Syint_mp5_x_L(i,j) ) ) /(s2-s1)
   FSzzlax(i,j,1,l) = (s2* FSzzint_mp5_L(i,j)-s1* FSzzint_mp5_R(i,j) &
        +s1*s2* (Szint_mp5_x_R(i,j) - Szint_mp5_x_L(i,j) ) ) /(s2-s1)

! -------------------------------------- EGLM-------------------------------------- 


   BxSxlax(i,j,1,l) = (s2* Bxint_mp5_x_L(i,j)-s1* Bxint_mp5_x_R(i,j) &
        +s1*s2* (Sxint_mp5_x_R(i,j) - Sxint_mp5_x_L(i,j) ) ) /(s2-s1) 
   BySxlax(i,j,1,l) = (s2* Byint_mp5_y_L(i,j)-s1* Byint_mp5_y_R(i,j) &
        +s1*s2* (Sxint_mp5_y_R(i,j) - Sxint_mp5_y_L(i,j) ) ) /(s2-s1)
   BzSxlax(i,j,1,l) = (s2* Bzint_mp5_x_L(i,j)-s1* Bzint_mp5_x_R(i,j) &
        +s1*s2* (Sxint_mp5_x_R(i,j) - Sxint_mp5_x_L(i,j) ) ) /(s2-s1)


   BxSylax(i,j,1,l) = (s2* Bxint_mp5_x_L(i,j)-s1* Bxint_mp5_x_R(i,j) &
        +s1*s2* (Syint_mp5_x_R(i,j) - Syint_mp5_x_L(i,j) ) ) /(s2-s1)
   BySylax(i,j,1,l) = (s2* Byint_mp5_y_L(i,j)-s1* Byint_mp5_y_R(i,j) &
        +s1*s2* (Syint_mp5_y_R(i,j) - Syint_mp5_y_L(i,j) ) ) /(s2-s1)
   BzSylax(i,j,1,l) = (s2* Bzint_mp5_x_L(i,j)-s1* Bzint_mp5_x_R(i,j) &
        +s1*s2* (Syint_mp5_x_R(i,j) - Syint_mp5_x_L(i,j) ) ) /(s2-s1)


   BxSzlax(i,j,1,l) = (s2* Bxint_mp5_x_L(i,j)-s1* Bxint_mp5_x_R(i,j) &
        +s1*s2* (Szint_mp5_x_R(i,j) - Szint_mp5_x_L(i,j) ) ) /(s2-s1)
   BySzlax(i,j,1,l) = (s2* Byint_mp5_y_L(i,j)-s1* Byint_mp5_y_R(i,j) &
        +s1*s2* (Szint_mp5_y_R(i,j) - Szint_mp5_y_L(i,j) ) ) /(s2-s1)
   BzSzlax(i,j,1,l) = (s2* Bzint_mp5_x_L(i,j)-s1* Bzint_mp5_x_R(i,j) &
        +s1*s2* (Szint_mp5_x_R(i,j) - Szint_mp5_x_L(i,j) ) ) /(s2-s1)

! -------------------------------------- EGLM-------------------------------------- 

 
        end do
     end do

!$OMP END DO
!$OMP END PARALLEL


     

   else

      write(*,*) "STOP: subroutine hllflux"
      write(*,*) "This dimension is not implemented yet"
      stop

   end if


 end subroutine hll_flow_mpx
