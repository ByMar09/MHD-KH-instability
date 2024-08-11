!     Auxiliar intermediate value for E (solution of the implicit part) IMEX - MIRK

  subroutine  electricvarint

    use scalar
    use parameters
    use threevectors
    use fourvectors
    implicit none

    call imex_schemes

!$OMP PARALLEL PRIVATE(sigma_bar,Edotv,vrotB_x,vrotB_y,vrotB_z)


    if (DIM == 1) then

       if (MIRK == 1) then


!$OMP DO 

       do i=0,imax

       Exint(i,1,1,l)  =  Ex(i,1,1) !Exast(i,1,1,l) !

       Eyint(i,1,1,l)  =  Ey(i,1,1) !Eyast(i,1,1,l) !

       Ezint(i,1,1,l)  =  Ez(i,1,1) !Ezast(i,1,1,l) !


       end do !i

!$OMP END DO

    else if (MIRK == 2) then

!$OMP DO 

       do i=0,imax


     if ( l== 1) then 

    Exint(i,1,1,l) = Ex(i,1,1) !Exast(i,1,1,l)
    Eyint(i,1,1,l) = Ey(i,1,1) !Eyast(i,1,1,l)
    Ezint(i,1,1,l) = Ez(i,1,1) !Ezast(i,1,1,l)


    else

    
    if(SLC == 1) sigma = sigma_loc(i,1,1)

    sigma_bar = Delt * sigma * W(i,1,1)


    Edotv = Ex(i,1,1) * Vx(i,1,1) +                      &
            Ey(i,1,1) * Vy(i,1,1) +                      &
            Ez(i,1,1) * Vz(i,1,1)  

    vrotB_x = (Bz(i,1,1)*Vy(i,1,1)-By(i,1,1)*Vz(i,1,1))
    vrotB_y = (Bx(i,1,1)*Vz(i,1,1)-Bz(i,1,1)*Vx(i,1,1))
    vrotB_z = (By(i,1,1)*Vx(i,1,1)-Bx(i,1,1)*Vy(i,1,1))
!!$
!!$    vrotB_x = (Bzint(i,1,1,l)*Vy(i,1,1)-Byint(i,1,1,l)*Vz(i,1,1))
!!$    vrotB_y = (Bxint(i,1,1,l)*Vz(i,1,1)-Bzint(i,1,1,l)*Vx(i,1,1))
!!$    vrotB_z = (Byint(i,1,1,l)*Vx(i,1,1)-Bxint(i,1,1,l)*Vy(i,1,1))



    Exint(i,1,1,l) = 1.d0 / (1.d0 + sigma_bar * ( 1.d0 - cm_1))  *  (                          &
                     Ex(i,1,1) * ( 1.d0 - sigma_bar * cm_1) + sigma_bar *  Vx(i,1,1) * Edotv   & 
                     - sigma_bar * vrotB_x + Delt * Exastsum(i,1,1) )

    Eyint(i,1,1,l) = 1.d0 / (1.d0 + sigma_bar * ( 1.d0 - cm_1))  *  (                          &
                     Ey(i,1,1) * ( 1.d0 - sigma_bar * cm_1) + sigma_bar *  Vy(i,1,1) * Edotv   & 
                     - sigma_bar * vrotB_y + Delt * Eyastsum(i,1,1) )

    Ezint(i,1,1,l) = 1.d0 / (1.d0 + sigma_bar * ( 1.d0 - cm_1))  *  (                          &
                     Ez(i,1,1) * ( 1.d0 - sigma_bar * cm_1) + sigma_bar *  Vz(i,1,1) * Edotv   & 
                     - sigma_bar * vrotB_z + Delt * Ezastsum(i,1,1) )

    end if

     end do !i

!$OMP END DO

     !     Auxiliar intermediate value for E (solution of the implicit part) IMEX

    else if (MIRK == 0 .and. REC_PRIM == 0) then

!$OMP DO 

    do i=0,imax


          if(SLC == 1) sigma = sigma_loc(i,1,1)

    a(l)    = Abt(l,l) * sigma * Delt

    matx(l) = W(i,1,1)**2*a(l) + W(i,1,1)*a(l)**2      &
            + W(i,1,1)  + a(l)

    vrotB_x = (Bzint(i,1,1,l)*Vy(i,1,1)-Byint(i,1,1,l)*Vz(i,1,1))
    vrotB_y = (Bxint(i,1,1,l)*Vz(i,1,1)-Bzint(i,1,1,l)*Vx(i,1,1))
    vrotB_z = (Byint(i,1,1,l)*Vx(i,1,1)-Bxint(i,1,1,l)*Vy(i,1,1))

    Exint(i,1,1,l) = (1./matx(l)) * (                                               &
         (a(l) + W(i,1,1)    + a(l) * W(i,1,1)**2 * Vx(i,1,1)**2 )* Exast(i,1,1,l)+ &
          a(l) * W(i,1,1)**2 * Vx(i,1,1)          * Vy(i,1,1)     * Eyast(i,1,1,l)+ &
          a(l) * W(i,1,1)**2 * Vx(i,1,1)          * Vz(i,1,1)     * Ezast(i,1,1,l)- &
          Abt(l,l) * sigma   * W(i,1,1)           * Delt          * (               &
         (a(l) + W(i,1,1)    + a(l) * W(i,1,1)**2 * Vx(i,1,1)**2 )* vrotB_x       + &
          a(l) * W(i,1,1)**2 * Vx(i,1,1)          * Vy(i,1,1)     * vrotB_y       + &
          a(l) * W(i,1,1)**2 * Vx(i,1,1)          * Vz(i,1,1)     * vrotB_z      ))

    Eyint(i,1,1,l) = (1./matx(l)) * (                                               &
          a(l) * W(i,1,1)**2 * Vx(i,1,1)          * Vy(i,1,1)     * Exast(i,1,1,l)+ &
         (a(l) + W(i,1,1)    + a(l) * W(i,1,1)**2 * Vy(i,1,1)**2 )* Eyast(i,1,1,l)+ &
          a(l) * W(i,1,1)**2 * Vy(i,1,1)          * Vz(i,1,1)     * Ezast(i,1,1,l)- &
          Abt(l,l) * sigma   * W(i,1,1)           * Delt * (                        &
          a(l) * W(i,1,1)**2 * Vx(i,1,1)          * Vy(i,1,1)     * vrotB_x       + &
         (a(l) + W(i,1,1)    + a(l) * W(i,1,1)**2 * Vy(i,1,1)**2 )* vrotB_y       + &
          a(l) * W(i,1,1)**2 * Vy(i,1,1)          * Vz(i,1,1)     * vrotB_z      ))

    Ezint(i,1,1,l) = (1./matx(l)) * (                                               &
          a(l) * W(i,1,1)**2 * Vx(i,1,1)          * Vz(i,1,1)     * Exast(i,1,1,l)+ &
          a(l) * W(i,1,1)**2 * Vy(i,1,1)          * Vz(i,1,1)     * Eyast(i,1,1,l)+ &
         (a(l) + W(i,1,1)    + a(l) * W(i,1,1)**2 * Vz(i,1,1)**2 )* Ezast(i,1,1,l)- &
          Abt(l,l) * sigma   * W(i,1,1)           * Delt * (                        &
          a(l) * W(i,1,1)**2 * Vx(i,1,1)          * Vz(i,1,1)     * vrotB_x       + &
          a(l) * W(i,1,1)**2 * Vy(i,1,1)          * Vz(i,1,1)     * vrotB_y       + &
         (a(l) + W(i,1,1)    + a(l) * W(i,1,1)**2 * Vz(i,1,1)**2 )* vrotB_z      ))


     end do !i

!$OMP END DO

    else if (MIRK == 0 .and. REC_PRIM == 1) then

!$OMP DO 

    do i=0,imax


          if(SLC == 1) sigma = sigma_loc(i,1,1)


    a(l)    = Abt(l,l) * sigma * Delt

    matx(l) = Wint(i,1,1,l)**2*a(l) + Wint(i,1,1,l)*a(l)**2      &
            + Wint(i,1,1,l)  + a(l)

    vrotB_x = (Bzint(i,1,1,l)*Vyint(i,1,1,l)-Byint(i,1,1,l)*Vzint(i,1,1,l))
    vrotB_y = (Bxint(i,1,1,l)*Vzint(i,1,1,l)-Bzint(i,1,1,l)*Vxint(i,1,1,l))
    vrotB_z = (Byint(i,1,1,l)*Vxint(i,1,1,l)-Bxint(i,1,1,l)*Vyint(i,1,1,l))

    Exint(i,1,1,l) = (1./matx(l)) * (                                               &
         (a(l) + Wint(i,1,1,l)    + a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)**2 )* Exast(i,1,1,l)+ &
          a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)          * Vyint(i,1,1,l)     * Eyast(i,1,1,l)+ &
          a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)          * Vzint(i,1,1,l)     * Ezast(i,1,1,l)- &
          Abt(l,l) * sigma   * Wint(i,1,1,l)           * Delt          * (               &
         (a(l) + Wint(i,1,1,l)    + a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)**2 )* vrotB_x       + &
          a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)          * Vyint(i,1,1,l)     * vrotB_y       + &
          a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)          * Vzint(i,1,1,l)     * vrotB_z      ))

    Eyint(i,1,1,l) = (1./matx(l)) * (                                               &
          a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)          * Vyint(i,1,1,l)     * Exast(i,1,1,l)+ &
         (a(l) + Wint(i,1,1,l)    + a(l) * Wint(i,1,1,l)**2 * Vyint(i,1,1,l)**2 )* Eyast(i,1,1,l)+ &
          a(l) * Wint(i,1,1,l)**2 * Vyint(i,1,1,l)          * Vzint(i,1,1,l)     * Ezast(i,1,1,l)- &
          Abt(l,l) * sigma   * Wint(i,1,1,l)           * Delt * (                        &
          a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)          * Vyint(i,1,1,l)     * vrotB_x       + &
         (a(l) + Wint(i,1,1,l)    + a(l) * Wint(i,1,1,l)**2 * Vyint(i,1,1,l)**2 )* vrotB_y       + &
          a(l) * Wint(i,1,1,l)**2 * Vyint(i,1,1,l)          * Vzint(i,1,1,l)     * vrotB_z      ))

    Ezint(i,1,1,l) = (1./matx(l)) * (                                               &
          a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)          * Vzint(i,1,1,l)     * Exast(i,1,1,l)+ &
          a(l) * Wint(i,1,1,l)**2 * Vyint(i,1,1,l)          * Vzint(i,1,1,l)     * Eyast(i,1,1,l)+ &
         (a(l) + Wint(i,1,1,l)    + a(l) * Wint(i,1,1,l)**2 * Vzint(i,1,1,l)**2 )* Ezast(i,1,1,l)- &
          Abt(l,l) * sigma   * Wint(i,1,1,l)           * Delt * (                        &
          a(l) * Wint(i,1,1,l)**2 * Vxint(i,1,1,l)          * Vzint(i,1,1,l)     * vrotB_x       + &
          a(l) * Wint(i,1,1,l)**2 * Vyint(i,1,1,l)          * Vzint(i,1,1,l)     * vrotB_y       + &
         (a(l) + Wint(i,1,1,l)    + a(l) * Wint(i,1,1,l)**2 * Vzint(i,1,1,l)**2 )* vrotB_z      ))


     end do !i

!$OMP END DO

    else


       write(*,*) "STOP: subroutine electricfield"
       write(*,*) "These MIRK scheme is not implemented yet"
       write(*,*) "Or recovery primitive intemediate (REC_PRIM) parameter is not valid"
       write(*,*) "MIRK =", MIRK, "REC_PRIM =", REC_PRIM
       stop

     stop

     end if


    else if (DIM == 2) then


    if (MIRK == 1) then

!$OMP DO 
       
       do j=0,jmax
          do i=0,imax
 


       Exint(i,j,1,l)  =  Ex(i,j,1) !Exast(i,j,1,l) !

       Eyint(i,j,1,l)  =  Ey(i,j,1) !Eyast(i,j,1,l) !

       Ezint(i,j,1,l)  =  Ez(i,j,1) !Ezast(i,j,1,l) !


           end do !i
        end do ! j

!$OMP END DO

    else if (MIRK == 2) then

!$OMP DO 

       do j=0,jmax
          do i=0,imax
 
          if(SLC == 1) sigma = sigma_loc(i,j,1)

    sigma_bar = Delt * sigma * W(i,j,1)


    Edotv = Ex(i,j,1) * Vx(i,j,1) +                      &
            Ey(i,j,1) * Vy(i,j,1) +                      &
            Ez(i,j,1) * Vz(i,j,1)  

    vrotB_x = (Bz(i,j,1)*Vy(i,j,1)-By(i,j,1)*Vz(i,j,1))
    vrotB_y = (Bx(i,j,1)*Vz(i,j,1)-Bz(i,j,1)*Vx(i,j,1))
    vrotB_z = (By(i,j,1)*Vx(i,j,1)-Bx(i,j,1)*Vy(i,j,1))

    if ( l== 1) then 

    Exint(i,j,1,l) = Exast(i,j,1,l)
    Eyint(i,j,1,l) = Eyast(i,j,1,l)
    Ezint(i,j,1,l) = Ezast(i,j,1,l)

    else


    Exint(i,j,1,l) = 1.d0 / (1.d0 + sigma_bar * ( 1.d0 - cm_1))  *  (                          &
                     Ex(i,j,1) * ( 1.d0 - sigma_bar * cm_1) + sigma_bar *  Vx(i,j,1) * Edotv   & 
                     - sigma_bar * vrotB_x + Delt * Exastsum(i,j,1) )

    Eyint(i,j,1,l) = 1.d0 / (1.d0 + sigma_bar * ( 1.d0 - cm_1))  *  (                          &
                     Ey(i,j,1) * ( 1.d0 - sigma_bar * cm_1) + sigma_bar *  Vy(i,j,1) * Edotv   & 
                     - sigma_bar * vrotB_y + Delt * Eyastsum(i,j,1) )

    Ezint(i,j,1,l) = 1.d0 / (1.d0 + sigma_bar * ( 1.d0 - cm_1))  *  (                          &
                     Ez(i,j,1) * ( 1.d0 - sigma_bar * cm_1) + sigma_bar *  Vz(i,j,1) * Edotv   & 
                     - sigma_bar * vrotB_z + Delt * Ezastsum(i,j,1) )

    end if

       end do !i
    end do !j

!$OMP END DO

    else if (MIRK == 0 .and. REC_PRIM == 0) then

!$OMP DO 

       do j=0,jmax
          do i=0,imax

          if(SLC == 1) sigma = sigma_loc(i,j,1)


    a(l)    = Abt(l,l) * sigma * Delt


    matx(l) = W(i,j,1)**2*a(l) + W(i,j,1)*a(l)**2 + W(i,j,1)  + a(l)

    vrotB_x = (Bzint(i,j,1,l)*Vy(i,j,1)-Byint(i,j,1,l)*Vz(i,j,1))
    vrotB_y = (Bxint(i,j,1,l)*Vz(i,j,1)-Bzint(i,j,1,l)*Vx(i,j,1))
    vrotB_z = (Byint(i,j,1,l)*Vx(i,j,1)-Bxint(i,j,1,l)*Vy(i,j,1))

    Exint(i,j,1,l) = (1./matx(l)) * (                                               &
         (a(l) + W(i,j,1)    + a(l) * W(i,j,1)**2 * Vx(i,j,1)**2 )* Exast(i,j,1,l)+ &
          a(l) * W(i,j,1)**2 * Vx(i,j,1)          * Vy(i,j,1)     * Eyast(i,j,1,l)+ &
          a(l) * W(i,j,1)**2 * Vx(i,j,1)          * Vz(i,j,1)     * Ezast(i,j,1,l)- &
          Abt(l,l) * sigma   * W(i,j,1)           * Delt          * (               &
         (a(l) + W(i,j,1)    + a(l) * W(i,j,1)**2 * Vx(i,j,1)**2 )* vrotB_x       + &
          a(l) * W(i,j,1)**2 * Vx(i,j,1)          * Vy(i,j,1)     * vrotB_y       + &
          a(l) * W(i,j,1)**2 * Vx(i,j,1)          * Vz(i,j,1)     * vrotB_z      ))

    Eyint(i,j,1,l) = (1./matx(l)) * (                                               &
          a(l) * W(i,j,1)**2 * Vx(i,j,1)          * Vy(i,j,1)     * Exast(i,j,1,l)+ &
         (a(l) + W(i,j,1)    + a(l) * W(i,j,1)**2 * Vy(i,j,1)**2 )* Eyast(i,j,1,l)+ &
          a(l) * W(i,j,1)**2 * Vy(i,j,1)          * Vz(i,j,1)     * Ezast(i,j,1,l)- &
          Abt(l,l) * sigma   * W(i,j,1)           * Delt * (                        &
          a(l) * W(i,j,1)**2 * Vx(i,j,1)          * Vy(i,j,1)     * vrotB_x       + &
         (a(l) + W(i,j,1)    + a(l) * W(i,j,1)**2 * Vy(i,j,1)**2 )* vrotB_y       + &
          a(l) * W(i,j,1)**2 * Vy(i,j,1)          * Vz(i,j,1)     * vrotB_z      ))

    Ezint(i,j,1,l) = (1./matx(l)) * (                                               &
          a(l) * W(i,j,1)**2 * Vx(i,j,1)          * Vz(i,j,1)     * Exast(i,j,1,l)+ &
          a(l) * W(i,j,1)**2 * Vy(i,j,1)          * Vz(i,j,1)     * Eyast(i,j,1,l)+ &
         (a(l) + W(i,j,1)    + a(l) * W(i,j,1)**2 * Vz(i,j,1)**2 )* Ezast(i,j,1,l)- &
          Abt(l,l) * sigma   * W(i,j,1)           * Delt * (                        &
          a(l) * W(i,j,1)**2 * Vx(i,j,1)          * Vz(i,j,1)     * vrotB_x       + &
          a(l) * W(i,j,1)**2 * Vy(i,j,1)          * Vz(i,j,1)     * vrotB_y       + &
         (a(l) + W(i,j,1)    + a(l) * W(i,j,1)**2 * Vz(i,j,1)**2 )* vrotB_z      ))

        end do !i
     end do !j

!$OMP END DO

    else if (MIRK == 0 .and. REC_PRIM == 1) then

!$OMP DO 

       do j=0,jmax
          do i=0,imax

          if(SLC == 1) sigma = sigma_loc(i,j,1)


    a(l)    = Abt(l,l) * sigma * Delt

    matx(l) = Wint(i,j,1,l)**2*a(l) + Wint(i,j,1,l)*a(l)**2      &
            + Wint(i,j,1,l)  + a(l)

    vrotB_x = (Bzint(i,j,1,l)*Vyint(i,j,1,l)-Byint(i,j,1,l)*Vzint(i,j,1,l))
    vrotB_y = (Bxint(i,j,1,l)*Vzint(i,j,1,l)-Bzint(i,j,1,l)*Vxint(i,j,1,l))
    vrotB_z = (Byint(i,j,1,l)*Vxint(i,j,1,l)-Bxint(i,j,1,l)*Vyint(i,j,1,l))

    Exint(i,j,1,l) = (1./matx(l)) * (                                               &
         (a(l) + Wint(i,j,1,l)    + a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)**2 )* Exast(i,j,1,l)+ &
          a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)          * Vyint(i,j,1,l)     * Eyast(i,j,1,l)+ &
          a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)          * Vzint(i,j,1,l)     * Ezast(i,j,1,l)- &
          Abt(l,l) * sigma   * Wint(i,j,1,l)           * Delt          * (               &
         (a(l) + Wint(i,j,1,l)    + a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)**2 )* vrotB_x       + &
          a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)          * Vyint(i,j,1,l)     * vrotB_y       + &
          a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)          * Vzint(i,j,1,l)     * vrotB_z      ))

    Eyint(i,j,1,l) = (1./matx(l)) * (                                               &
          a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)          * Vyint(i,j,1,l)     * Exast(i,j,1,l)+ &
         (a(l) + Wint(i,j,1,l)    + a(l) * Wint(i,j,1,l)**2 * Vyint(i,j,1,l)**2 )* Eyast(i,j,1,l)+ &
          a(l) * Wint(i,j,1,l)**2 * Vyint(i,j,1,l)          * Vzint(i,j,1,l)     * Ezast(i,j,1,l)- &
          Abt(l,l) * sigma   * Wint(i,j,1,l)           * Delt * (                        &
          a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)          * Vyint(i,j,1,l)     * vrotB_x       + &
         (a(l) + Wint(i,j,1,l)    + a(l) * Wint(i,j,1,l)**2 * Vyint(i,j,1,l)**2 )* vrotB_y       + &
          a(l) * Wint(i,j,1,l)**2 * Vyint(i,j,1,l)          * Vzint(i,j,1,l)     * vrotB_z      ))

    Ezint(i,j,1,l) = (1./matx(l)) * (                                               &
          a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)          * Vzint(i,j,1,l)     * Exast(i,j,1,l)+ &
          a(l) * Wint(i,j,1,l)**2 * Vyint(i,j,1,l)          * Vzint(i,j,1,l)     * Eyast(i,j,1,l)+ &
         (a(l) + Wint(i,j,1,l)    + a(l) * Wint(i,j,1,l)**2 * Vzint(i,j,1,l)**2 )* Ezast(i,j,1,l)- &
          Abt(l,l) * sigma   * Wint(i,j,1,l)           * Delt * (                        &
          a(l) * Wint(i,j,1,l)**2 * Vxint(i,j,1,l)          * Vzint(i,j,1,l)     * vrotB_x       + &
          a(l) * Wint(i,j,1,l)**2 * Vyint(i,j,1,l)          * Vzint(i,j,1,l)     * vrotB_y       + &
         (a(l) + Wint(i,j,1,l)    + a(l) * Wint(i,j,1,l)**2 * Vzint(i,j,1,l)**2 )* vrotB_z      ))


        end do !i
     end do !j

!$OMP END DO

    else

      write(*,*) "STOP: subroutine electricfield"
       write(*,*) "These MIRK scheme is not implemented yet"
       write(*,*) "Or recovery primitive intemediate (REC_PRIM) parameter is not valid"
       write(*,*) "MIRK =", MIRK, "REC_PRIM =", REC_PRIM
       stop

    end if

   else

      write(*,*) "STOP: subroutine electricvarint"
      write(*,*) "This dimension is not implemented yet"
      stop

   end if

!$OMP END PARALLEL   
   
   call boundary_electric


  end subroutine electricvarint
