  !     *******************************************************************
  !     Subroutine which make the sum of flow for the final variables
  !     *******************************************************************

  subroutine sumfin


    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

    call imex_schemes

!$OMP PARALLEL
!$OMP DO 

    do j=0+fix,jmax
       do i=0,imax

    psifinsum(i,j,1) = psifinsum(i,j,1) +  wb(l) * psiflux(i,j,1,l)
 
    phifinsum(i,j,1) = phifinsum(i,j,1) +  wb(l) * phiflux(i,j,1,l)

    Bxfinsum(i,j,1)  = Bxfinsum(i,j,1)  +  wb(l) * Bxflux(i,j,1,l)

    Byfinsum(i,j,1)  = Byfinsum(i,j,1)  +  wb(l) * Byflux(i,j,1,l)


    Bzfinsum(i,j,1)  =  Bzfinsum(i,j,1) +  wb(l) * Bzflux(i,j,1,l)


    Exfinsum(i,j,1)  =  Exfinsum(i,j,1) + (wb (l) * Exastflux(i,j,1,l)   &
                                        -  wbt(l) * Exastsour(i,j,1,l) )  

    Eyfinsum(i,j,1)  =  Eyfinsum(i,j,1) + (wb (l) * Eyastflux(i,j,1,l)   &
                                        -  wbt(l) * Eyastsour(i,j,1,l) ) 

    Ezfinsum(i,j,1)  =  Ezfinsum(i,j,1) + (wb (l) * Ezastflux(i,j,1,l)   &
                                        -  wbt(l) * Ezastsour(i,j,1,l) ) 

    qfinsum(i,j,1)   =  qfinsum(i,j,1)  +  wb(l)  * qflux(i,j,1,l)

    Dfinsum(i,j,1)   =  Dfinsum(i,j,1)  +  wb(l)  * Dflux(i,j,1,l)

! -------------------------------------- EGLM--------------------------------------

    taufinsum(i,j,1) = taufinsum(i,j,1) +  wb(l)    * tauflux(i,j,1,l)    &
                                        +  FAC_EGLM * wb(l)  * taueglm(i,j,1,l)

    Sxfinsum(i,j,1)  = Sxfinsum(i,j,1)  +  wb(l)    * Sxflux(i,j,1,l)     &
                                        +  FAC_EGLM * wb(l)  * Sxeglm(i,j,1,l)

    Syfinsum(i,j,1)  = Syfinsum(i,j,1)  +  wb(l)    * Syflux(i,j,1,l)     &
                                        +  FAC_EGLM * wb(l)  * Syeglm(i,j,1,l)

    Szfinsum(i,j,1)  = Szfinsum(i,j,1)  +  wb(l)    * Szflux(i,j,1,l)      &
                                        +  FAC_EGLM * wb(l)  * Szeglm(i,j,1,l)

    ! -------------------------------------- EGLM--------------------------------------


            end do ! for i
         end do ! for j

!$OMP END DO
!$OMP END PARALLEL

  end subroutine sumfin
