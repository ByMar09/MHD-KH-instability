  !     *******************************************************************
  !     Subroutine which initialized  flows of the intermediate variables
  !     *******************************************************************

   subroutine flowinit

    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

!$OMP PARALLEL PRIVATE(Edotv,vrotB_x,vrotB_y,vrotB_z)
    
    if (DIM == 1) then

       if(REC_PRIM == 0) then

!$OMP DO 

        do i=0,imax

           if (SLC == 1) sigma = sigma_loc(i,1,1)

    !     Electric and Magnetic field modules 


    E2int(i,1,1,l) = Exint(i,1,1,l)**2 + Eyint(i,1,1,l)**2   &
                   + Ezint(i,1,1,l)**2

    B2int(i,1,1,l) = Bxint(i,1,1,l)**2 + Byint(i,1,1,l)**2   &
                   + Bzint(i,1,1,l)**2


    Edotv = Exint(i,1,1,l) * Vx(i,1,1) +                      &
            Eyint(i,1,1,l) * Vy(i,1,1) +                      &
            Ezint(i,1,1,l) * Vz(i,1,1)  

    vrotB_x = (Bzint(i,1,1,l)*Vy(i,1,1)-Byint(i,1,1,l)*Vz(i,1,1))
    vrotB_y = (Bxint(i,1,1,l)*Vz(i,1,1)-Bzint(i,1,1,l)*Vx(i,1,1))
    vrotB_z = (Byint(i,1,1,l)*Vx(i,1,1)-Bxint(i,1,1,l)*Vy(i,1,1))
!!$
!!$    vrotB_x = (Bz(i,1,1)*Vy(i,1,1)-By(i,1,1)*Vz(i,1,1))
!!$    vrotB_y = (Bx(i,1,1)*Vz(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
!!$    vrotB_z = (By(i,1,1)*Vx(i,1,1)-Bx(i,1,1)*Vy(i,1,1))

    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxint(i,1,1,l) = sigma * W(i,1,1) *                       &
         ( Exint(i,1,1,l) +  vrotB_x  - Edotv * Vx(i,1,1) )   &
         + qint(i,1,1,l) *Vx(i,1,1)

!    print*, i, l, Jxint(i,1,1,l), vrotB_x, Bzint(i,1,1,l), Byint(i,1,1,l), Vy(i,1,1), "see subrotine 14_flowinit.f95"


    Jyint(i,1,1,l) = sigma * W(i,1,1) *                       &
         ( Eyint(i,1,1,l) +  vrotB_y  - Edotv * Vy(i,1,1) )   &
         + qint(i,1,1,l) * Vy(i,1,1)

    Jzint(i,1,1,l) = sigma * W(i,1,1) *                       &       
         ( Ezint(i,1,1,l) +  vrotB_z  - Edotv * Vz(i,1,1) )   &
         + qint(i,1,1,l) * Vz(i,1,1)

    !     Density flux (FD = \rho W V)

    FDxint(i,1,1,l) = DDint(i,1,1,l) * Vx(i,1,1)

    FDyint(i,1,1,l) = DDint(i,1,1,l) * Vy(i,1,1)

    FDzint(i,1,1,l) = DDint(i,1,1,l) * Vz(i,1,1)

    !     Flujos conservados de energía

    Ftauxint(i,1,1,l) = (Bzint(i,1,1,l)*Eyint(i,1,1,l)-       &
         Byint(i,1,1,l) * Ezint(i,1,1,l) ) +                  &
         enthpy(i,1,1)*W(i,1,1)**2 * Vx(i,1,1)   

    Ftauyint(i,1,1,l) = (Bxint(i,1,1,l)*Ezint(i,1,1,l)-       &
         Bzint(i,1,1,l) * Exint(i,1,1,l) ) +                  &
         enthpy(i,1,1)*W(i,1,1)**2 * Vy(i,1,1)

    Ftauzint(i,1,1,l) = (Byint(i,1,1,l)*Exint(i,1,1,l)-       &
         Bxint(i,1,1,l) * Eyint(i,1,1,l) ) +                  &
         enthpy(i,1,1)*W(i,1,1)**2 * Vz(i,1,1)

    !     Components of the flux momentum tensor 

    FSxxint(i,1,1,l) =  - Exint(i,1,1,l)**2 -                &
         Bxint(i,1,1,l)**2 + enthpy(i,1,1) *                   &  
         W(i,1,1)**2 * Vx(i,1,1)**2 +                         &
         0.5*(E2int(i,1,1,l)+B2int(i,1,1,l))+p(i,1,1)  

    FSxyint(i,1,1,l)  = -Exint(i,1,1,l)*Eyint(i,1,1,l) -      &
         Bxint(i,1,1,l)*Byint(i,1,1,l) +                      &
         enthpy(i,1,1) * W(i,1,1)**2*Vx(i,1,1) *              &
         Vy(i,1,1)

    FSxzint(i,1,1,l) = -Exint(i,1,1,l)*Ezint(i,1,1,l) -       &
         Bxint(i,1,1,l)*Bzint(i,1,1,l)+enthpy(i,1,1)*         &
         W(i,1,1)**2 * Vx(i,1,1)*Vz(i,1,1)

    FSyxint(i,1,1,l) = -Exint(i,1,1,l)*Eyint(i,1,1,l) -       &
         Bxint(i,1,1,l)*Byint(i,1,1,l) +                      &
         enthpy(i,1,1) * W(i,1,1)**2 * Vx(i,1,1) *              &
         Vy(i,1,1)

    FSyyint(i,1,1,l) = -Eyint(i,1,1,l)**2 -                   &
         Byint(i,1,1,l)**2 +enthpy(i,1,1) *                   &
         W(i,1,1)**2 * Vy(i,1,1)**2 +                         &
         0.5*(E2int(i,1,1,l)+B2int(i,1,1,l))+p(i,1,1)

    FSyzint(i,1,1,l) = -Eyint(i,1,1,l)*Ezint(i,1,1,l) -       &
         Byint(i,1,1,l)*Bzint(i,1,1,l) +enthpy(i,1,1)*        &
         W(i,1,1)**2*Vy(i,1,1)*Vz(i,1,1)

    FSzxint(i,1,1,l) = -Exint(i,1,1,l)*Ezint(i,1,1,l) -       &
         Bxint(i,1,1,l)*Bzint(i,1,1,l)+enthpy(i,1,1)*         &
         W(i,1,1)**2 * Vx(i,1,1)*Vz(i,1,1)                   

    FSzyint(i,1,1,l) = -Eyint(i,1,1,l)*Ezint(i,1,1,l) -       &
         Byint(i,1,1,l)*Bzint(i,1,1,l) +enthpy(i,1,1)*        &
         W(i,1,1)**2*Vy(i,1,1)*Vz(i,1,1)

    FSzzint(i,1,1,l) = - Ezint(i,1,1,l)**2 -                  &
         Bzint(i,1,1,l)**2 +enthpy(i,1,1) *                   &
         W(i,1,1)**2 * Vz(i,1,1)**2 +                         &
         0.5*(E2int(i,1,1,l)+B2int(i,1,1,l))+p(i,1,1)



        end do ! for i

!$OMP END DO

       else if (REC_PRIM ==1) then

!$OMP DO 

        do i=0,imax

           if (SLC == 1) sigma = sigma_loc(i,1,1)

    !     Electric and Magnetic field modules 


    E2int(i,1,1,l) = Exint(i,1,1,l)**2 + Eyint(i,1,1,l)**2   &
                   + Ezint(i,1,1,l)**2

    B2int(i,1,1,l) = Bxint(i,1,1,l)**2 + Byint(i,1,1,l)**2   &
                   + Bzint(i,1,1,l)**2


    Edotv = Exint(i,1,1,l) * Vxint(i,1,1,l) +                      &
            Eyint(i,1,1,l) * Vyint(i,1,1,l) +                      &
            Ezint(i,1,1,l) * Vzint(i,1,1,l)  

    vrotB_x = (Bzint(i,1,1,l)*Vyint(i,1,1,l)-Byint(i,1,1,l)*Vzint(i,1,1,l))
    vrotB_y = (Bxint(i,1,1,l)*Vzint(i,1,1,l)-Bzint(i,1,1,l)*Vxint(i,1,1,l))
    vrotB_z = (Byint(i,1,1,l)*Vxint(i,1,1,l)-Bxint(i,1,1,l)*Vyint(i,1,1,l))

    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxint(i,1,1,l) = sigma * Wint(i,1,1,l) *                       &
         ( Exint(i,1,1,l) +  vrotB_x  - Edotv * Vxint(i,1,1,l) )   &
         + qint(i,1,1,l) *Vxint(i,1,1,l)


    Jyint(i,1,1,l) = sigma * Wint(i,1,1,l) *                       &
         ( Eyint(i,1,1,l) +  vrotB_y  - Edotv * Vyint(i,1,1,l) )   &
         + qint(i,1,1,l) * Vyint(i,1,1,l)

    Jzint(i,1,1,l) = sigma * Wint(i,1,1,l) *                       &       
         ( Ezint(i,1,1,l) +  vrotB_z  - Edotv * Vzint(i,1,1,l) )   &
         + qint(i,1,1,l) * Vzint(i,1,1,l)

    !     Density flux (FD = \rho W V)

    FDxint(i,1,1,l) = DDint(i,1,1,l) * Vxint(i,1,1,l)

    FDyint(i,1,1,l) = DDint(i,1,1,l) * Vyint(i,1,1,l)

    FDzint(i,1,1,l) = DDint(i,1,1,l) * Vzint(i,1,1,l)

    !     Flujos conservados de energía

    Ftauxint(i,1,1,l) = (Bzint(i,1,1,l)*Eyint(i,1,1,l)-       &
         Byint(i,1,1,l) * Ezint(i,1,1,l) ) +                  &
         enthpyint(i,1,1,l)*Wint(i,1,1,l)**2 * Vxint(i,1,1,l)   

    Ftauyint(i,1,1,l) = (Bxint(i,1,1,l)*Ezint(i,1,1,l)-       &
         Bzint(i,1,1,l) * Exint(i,1,1,l) ) +                  &
         enthpyint(i,1,1,l)*Wint(i,1,1,l)**2 * Vyint(i,1,1,l)

    Ftauzint(i,1,1,l) = (Byint(i,1,1,l)*Exint(i,1,1,l)-       &
         Bxint(i,1,1,l) * Eyint(i,1,1,l) ) +                  &
         enthpyint(i,1,1,l)*Wint(i,1,1,l)**2 * Vzint(i,1,1,l)

    !     Components of the flux momentum tensor 

    FSxxint(i,1,1,l) =  - Exint(i,1,1,l)**2 -                &
         Bxint(i,1,1,l)**2 + enthpyint(i,1,1,l) *            &  
         Wint(i,1,1,l)**2 * Vxint(i,1,1,l)**2 +              &
         0.5*(E2int(i,1,1,l)+B2int(i,1,1,l))+pint(i,1,1,l)  

    FSxyint(i,1,1,l)  = -Exint(i,1,1,l)*Eyint(i,1,1,l) -        &
         Bxint(i,1,1,l)*Byint(i,1,1,l) +                        &
         enthpyint(i,1,1,l) * Wint(i,1,1,l)**2*Vxint(i,1,1,l) * &
         Vyint(i,1,1,l)

    FSxzint(i,1,1,l) = -Exint(i,1,1,l)*Ezint(i,1,1,l) -       &
         Bxint(i,1,1,l)*Bzint(i,1,1,l)+enthpyint(i,1,1,l)*    &
         Wint(i,1,1,l)**2 * Vxint(i,1,1,l)*Vzint(i,1,1,l)

    FSyxint(i,1,1,l) = -Exint(i,1,1,l)*Eyint(i,1,1,l) -           &
         Bxint(i,1,1,l)*Byint(i,1,1,l) +                          &
         enthpyint(i,1,1,l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l) * &
         Vyint(i,1,1,l)

    FSyyint(i,1,1,l) = -Eyint(i,1,1,l)**2 -                   &
         Byint(i,1,1,l)**2 +enthpyint(i,1,1,l) *              &
         Wint(i,1,1,l)**2 * Vyint(i,1,1,l)**2 +               &
         0.5*(E2int(i,1,1,l)+B2int(i,1,1,l))+pint(i,1,1,l)

    FSyzint(i,1,1,l) = -Eyint(i,1,1,l)*Ezint(i,1,1,l) -       &
         Byint(i,1,1,l)*Bzint(i,1,1,l) +enthpyint(i,1,1,l)*   &
         Wint(i,1,1,l)**2*Vyint(i,1,1,l)*Vzint(i,1,1,l)

    FSzxint(i,1,1,l) = -Exint(i,1,1,l)*Ezint(i,1,1,l) -       &
         Bxint(i,1,1,l)*Bzint(i,1,1,l)+enthpyint(i,1,1,l)*    &
         Wint(i,1,1,l)**2 * Vxint(i,1,1,l)*Vzint(i,1,1,l)                   

    FSzyint(i,1,1,l) = -Eyint(i,1,1,l)*Ezint(i,1,1,l) -       &
         Byint(i,1,1,l)*Bzint(i,1,1,l) +enthpyint(i,1,1,l)*   &
         Wint(i,1,1,l)**2*Vyint(i,1,1,l)*Vzint(i,1,1,l)

    FSzzint(i,1,1,l) = - Ezint(i,1,1,l)**2 -                  &
         Bzint(i,1,1,l)**2 +enthpyint(i,1,1,l) *              &
         Wint(i,1,1,l)**2 * Vzint(i,1,1,l)**2 +               &
         0.5*(E2int(i,1,1,l)+B2int(i,1,1,l))+pint(i,1,1,l)



        end do ! for i

!$OMP END DO

       else

          write(*,*) "STOP: subroutine flowinit"
          write(*,*) "REC_PRIM parameter is not valid"
          stop

       end if



    else if (DIM ==2) then

       if (REC_PRIM == 0) then

!$OMP DO 
 
          do j=0,jmax
             do i=0,imax


              if (SLC == 1) sigma = sigma_loc(i,j,1)

    !     Electric and Magnetic field modules 


    E2int(i,j,1,l) = Exint(i,j,1,l)**2 + Eyint(i,j,1,l)**2   &
                   + Ezint(i,j,1,l)**2

    B2int(i,j,1,l) = Bxint(i,j,1,l)**2 + Byint(i,j,1,l)**2   &
                   + Bzint(i,j,1,l)**2


    Edotv = Exint(i,j,1,l) * Vx(i,j,1) +                      &
            Eyint(i,j,1,l) * Vy(i,j,1) +                      &
            Ezint(i,j,1,l) * Vz(i,j,1)  

    vrotB_x = (Bzint(i,j,1,l)*Vy(i,j,1)-Byint(i,j,1,l)*Vz(i,j,1))
    vrotB_y = (Bxint(i,j,1,l)*Vz(i,j,1)-Bzint(i,j,1,l)*Vx(i,j,1))
    vrotB_z = (Byint(i,j,1,l)*Vx(i,j,1)-Bxint(i,j,1,l)*Vy(i,j,1))

    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV


    Jxint(i,j,1,l) = sigma * W(i,j,1) *                       &
         ( Exint(i,j,1,l) +  vrotB_x  - Edotv * Vx(i,j,1) )   &
         + qint(i,j,1,l) *Vx(i,j,1)


    Jyint(i,j,1,l) = sigma * W(i,j,1) *                       &
         ( Eyint(i,j,1,l) +  vrotB_y  - Edotv * Vy(i,j,1) )   &
         + qint(i,j,1,l) * Vy(i,j,1)

    Jzint(i,j,1,l) = sigma * W(i,j,1) *                       &       
         ( Ezint(i,j,1,l) +  vrotB_z  - Edotv * Vz(i,j,1) )   &
         + qint(i,j,1,l) * Vz(i,j,1)

    !     Density flux (FD = \rho W V)


    FDxint(i,j,1,l) = DDint(i,j,1,l) * Vx(i,j,1)

    FDyint(i,j,1,l) = DDint(i,j,1,l) * Vy(i,j,1)

    FDzint(i,j,1,l) = DDint(i,j,1,l) * Vz(i,j,1)

    !     Flujos conservados de energía

    Ftauxint(i,j,1,l) = (Bzint(i,j,1,l)*Eyint(i,j,1,l)-       &
         Byint(i,j,1,l) * Ezint(i,j,1,l) ) +                  &
         enthpy(i,j,1)*W(i,j,1)**2 * Vx(i,j,1)   

    Ftauyint(i,j,1,l) = (Bxint(i,j,1,l)*Ezint(i,j,1,l)-       &
         Bzint(i,j,1,l) * Exint(i,j,1,l) ) +                  &
         enthpy(i,j,1)*W(i,j,1)**2 * Vy(i,j,1)

    Ftauzint(i,j,1,l) = (Byint(i,j,1,l)*Exint(i,j,1,l)-       &
         Bxint(i,j,1,l) * Eyint(i,j,1,l) ) +                  &
         enthpy(i,j,1)*W(i,j,1)**2 * Vz(i,j,1)

    !     Components of the flux momentum tensor 

    FSxxint(i,j,1,l) = - Exint(i,j,1,l)**2 -                  &
         Bxint(i,j,1,l)**2 + enthpy(i,j,1) *                  &
         W(i,j,1)**2 * Vx(i,j,1)**2 +                         &
         0.5*(E2int(i,j,1,l)+B2int(i,j,1,l))+p(i,j,1)  

    FSxyint(i,j,1,l)  = -Exint(i,j,1,l)*Eyint(i,j,1,l) -      &
         Bxint(i,j,1,l)*Byint(i,j,1,l) +                      &
         enthpy(i,j,1) * W(i,j,1)**2*Vx(i,j,1) *              &
         Vy(i,j,1)

    FSxzint(i,j,1,l) = -Exint(i,j,1,l)*Ezint(i,j,1,l) -       &
         Bxint(i,j,1,l)*Bzint(i,j,1,l)+enthpy(i,j,1)*         &
         W(i,j,1)**2 * Vx(i,j,1)*Vz(i,j,1)

    FSyxint(i,j,1,l) = -Exint(i,j,1,l)*Eyint(i,j,1,l) -       &
         Bxint(i,j,1,l)*Byint(i,j,1,l) +                      &
         enthpy(i,j,1) * W(i,j,1)**2 * Vx(i,j,1) *            &
         Vy(i,j,1)

    FSyyint(i,j,1,l) = -Eyint(i,j,1,l)**2 -                   &
         Byint(i,j,1,l)**2 +enthpy(i,j,1) *                   &
         W(i,j,1)**2 * Vy(i,j,1)**2 +                         &
         0.5*(E2int(i,j,1,l)+B2int(i,j,1,l))+p(i,j,1)

    FSyzint(i,j,1,l) = -Eyint(i,j,1,l)*Ezint(i,j,1,l) -       &
         Byint(i,j,1,l)*Bzint(i,j,1,l) +enthpy(i,j,1)*        &
         W(i,j,1)**2*Vy(i,j,1)*Vz(i,j,1)

    FSzxint(i,j,1,l) = -Exint(i,j,1,l)*Ezint(i,j,1,l) -       &
         Bxint(i,j,1,l)*Bzint(i,j,1,l)+enthpy(i,j,1)*         &
         W(i,j,1)**2 * Vx(i,j,1)*Vz(i,j,1)                   

    FSzyint(i,j,1,l) = -Eyint(i,j,1,l)*Ezint(i,j,1,l) -       &
         Byint(i,j,1,l)*Bzint(i,j,1,l) +enthpy(i,j,1)*        &
         W(i,j,1)**2*Vy(i,j,1)*Vz(i,j,1)

    FSzzint(i,j,1,l) = - Ezint(i,j,1,l)**2 -                  &
         Bzint(i,j,1,l)**2 +enthpy(i,j,1) *                   &
         W(i,j,1)**2 * Vz(i,j,1)**2 +                         &
         0.5*(E2int(i,j,1,l)+B2int(i,j,1,l))+p(i,j,1)


           end do ! for i
        end do ! for j

!$OMP END DO

     else if (REC_PRIM == 1) then

!$OMP DO 

        do j=0,jmax
           do i=0,imax

              if (SLC == 1) sigma = sigma_loc(i,j,1)


    E2int(i,j,1,l) = Exint(i,j,1,l)**2 + Eyint(i,j,1,l)**2   &
                   + Ezint(i,j,1,l)**2

    B2int(i,j,1,l) = Bxint(i,j,1,l)**2 + Byint(i,j,1,l)**2   &
                   + Bzint(i,j,1,l)**2


    Edotv = Exint(i,j,1,l) * Vxint(i,j,1,l) +                      &
            Eyint(i,j,1,l) * Vyint(i,j,1,l) +                      &
            Ezint(i,j,1,l) * Vzint(i,j,1,l)  

    vrotB_x = (Bzint(i,j,1,l)*Vyint(i,j,1,l)-Byint(i,j,1,l)*Vzint(i,j,1,l))
    vrotB_y = (Bxint(i,j,1,l)*Vzint(i,j,1,l)-Bzint(i,j,1,l)*Vxint(i,j,1,l))
    vrotB_z = (Byint(i,j,1,l)*Vxint(i,j,1,l)-Bxint(i,j,1,l)*Vyint(i,j,1,l))

    !     Current density expression J = sigma * W * (E + V X B - (E.V)V) + qV

    Jxint(i,j,1,l) = sigma * Wint(i,j,1,l) *                       &
         ( Exint(i,j,1,l) +  vrotB_x  - Edotv * Vxint(i,j,1,l) )   &
         + qint(i,j,1,l) *Vxint(i,j,1,l)


    Jyint(i,j,1,l) = sigma * Wint(i,j,1,l) *                       &
         ( Eyint(i,j,1,l) +  vrotB_y  - Edotv * Vyint(i,j,1,l) )   &
         + qint(i,j,1,l) * Vyint(i,j,1,l)

    Jzint(i,j,1,l) = sigma * Wint(i,j,1,l) *                       &       
         ( Ezint(i,j,1,l) +  vrotB_z  - Edotv * Vzint(i,j,1,l) )   &
         + qint(i,j,1,l) * Vzint(i,j,1,l)

    !     Density flux (FD = \rho W V)

    FDxint(i,j,1,l) = DDint(i,j,1,l) * Vxint(i,j,1,l)

    FDyint(i,j,1,l) = DDint(i,j,1,l) * Vyint(i,j,1,l)

    FDzint(i,j,1,l) = DDint(i,j,1,l) * Vzint(i,j,1,l)

    !     Flujos conservados de energía

    Ftauxint(i,j,1,l) = (Bzint(i,j,1,l)*Eyint(i,j,1,l)-       &
         Byint(i,j,1,l) * Ezint(i,j,1,l) ) +                  &
         enthpyint(i,j,1,l)*Wint(i,j,1,l)**2 * Vxint(i,j,1,l)   

    Ftauyint(i,j,1,l) = (Bxint(i,j,1,l)*Ezint(i,j,1,l)-       &
         Bzint(i,j,1,l) * Exint(i,j,1,l) ) +                  &
         enthpyint(i,j,1,l)*Wint(i,j,1,l)**2 * Vyint(i,j,1,l)

    Ftauzint(i,j,1,l) = (Byint(i,j,1,l)*Exint(i,j,1,l)-       &
         Bxint(i,j,1,l) * Eyint(i,j,1,l) ) +                  &
         enthpyint(i,j,1,l)*Wint(i,j,1,l)**2 * Vzint(i,j,1,l)

    !     Components of the flux momentum tensor 

    FSxxint(i,j,1,l) =  - Exint(i,j,1,l)**2 -                &
         Bxint(i,j,1,l)**2 + enthpyint(i,j,1,l) *            &  
         Wint(i,j,1,l)**2 * Vxint(i,j,1,l)**2 +              &
         0.5d0*(E2int(i,j,1,l)+B2int(i,j,1,l))+pint(i,j,1,l)  

    FSxyint(i,j,1,l)  = -Exint(i,j,1,l)*Eyint(i,j,1,l) -        &
         Bxint(i,j,1,l)*Byint(i,j,1,l) +                        &
         enthpyint(i,j,1,l) * Wint(i,j,1,l)**2*Vxint(i,j,1,l) * &
         Vyint(i,j,1,l)

    FSxzint(i,j,1,l) = -Exint(i,j,1,l)*Ezint(i,j,1,l) -          &
         Bxint(i,j,1,l)*Bzint(i,j,1,l)+enthpyint(i,j,1,l)*       &
         Wint(i,j,1,l)**2 * Vxint(i,j,1,l)*Vzint(i,j,1,l)

    FSyxint(i,j,1,l) = -Exint(i,j,1,l)*Eyint(i,j,1,l) -           &
         Bxint(i,j,1,l)*Byint(i,j,1,l) +                          &
         enthpyint(i,j,1,l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l) * &
         Vyint(i,j,1,l)

    FSyyint(i,j,1,l) = -Eyint(i,j,1,l)**2 -                   &
         Byint(i,j,1,l)**2 +enthpyint(i,j,1,l) *              &
         Wint(i,j,1,l)**2 * Vyint(i,j,1,l)**2 +               &
         0.5d0*(E2int(i,j,1,l)+B2int(i,j,1,l))+pint(i,j,1,l)

    FSyzint(i,j,1,l) = -Eyint(i,j,1,l)*Ezint(i,j,1,l) -       &
         Byint(i,j,1,l)*Bzint(i,j,1,l) +enthpyint(i,j,1,l)*   &
         Wint(i,j,1,l)**2*Vyint(i,j,1,l)*Vzint(i,j,1,l)

    FSzxint(i,j,1,l) = -Exint(i,j,1,l)*Ezint(i,j,1,l) -       &
         Bxint(i,j,1,l)*Bzint(i,j,1,l)+enthpyint(i,j,1,l)*    &
         Wint(i,j,1,l)**2 * Vxint(i,j,1,l)*Vzint(i,j,1,l)                   

    FSzyint(i,j,1,l) = -Eyint(i,j,1,l)*Ezint(i,j,1,l) -       &
         Byint(i,j,1,l)*Bzint(i,j,1,l) +enthpyint(i,j,1,l)*   &
         Wint(i,j,1,l)**2*Vyint(i,j,1,l)*Vzint(i,j,1,l)

    FSzzint(i,j,1,l) = - Ezint(i,j,1,l)**2 -                  &
         Bzint(i,j,1,l)**2 +enthpyint(i,j,1,l) *              &
         Wint(i,j,1,l)**2 * Vzint(i,j,1,l)**2 +               &
         0.5d0*(E2int(i,j,1,l)+B2int(i,j,1,l))+pint(i,j,1,l)


           end do ! for i
        end do ! for j

!$OMP END DO

     else

          write(*,*) "STOP: subroutine flowinit"
          write(*,*) "REC_PRIM parameter is not valid"
          stop

       end if

    else

      write(*,*) "STOP: subroutine flowinit"
      write(*,*) "This dimension is not implemented yet"
      stop


   end if

!$OMP END PARALLEL

  end subroutine flowinit
