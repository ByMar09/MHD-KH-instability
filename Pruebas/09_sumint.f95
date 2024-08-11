  !     *******************************************************************
  !     Subroutine which make the sum of flows to the intermediate variables
  !     *******************************************************************

 subroutine sumint 

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


    psiintsum(i,j,1) =  psiintsum(i,j,1) +  Ab(l,m)  *  psiflux(i,j,1,m) 

    phiintsum(i,j,1) =  phiintsum(i,j,1) +  Ab(l,m)  *  phiflux(i,j,1,m)

    Bxintsum(i,j,1)  =  Bxintsum(i,j,1)  +  Ab(l,m)  *  Bxflux (i,j,1,m)


    Byintsum(i,j,1)  =  Byintsum(i,j,1)  +  Ab(l,m)  *  Byflux (i,j,1,m)  

    Bzintsum(i,j,1)  =  Bzintsum(i,j,1)  +  Ab(l,m)  *  Bzflux (i,j,1,m) 


    Exastsum(i,j,1)  =  Exastsum(i,j,1)  + ( Ab (l,m) *  Exastflux(i,j,1,m)  &
                                           - Abt(l,m) *  Exastsour(i,j,1,m) )


    Eyastsum(i,j,1)  =  Eyastsum(i,j,1)  + ( Ab (l,m) *  Eyastflux(i,j,1,m)  &
                                           - Abt(l,m) *  Eyastsour(i,j,1,m) )

    Ezastsum(i,j,1)  =  Ezastsum(i,j,1)  + ( Ab (l,m) *  Ezastflux(i,j,1,m)  &
                                           - Abt(l,m) *  Ezastsour(i,j,1,m) )

  
    qintsum(i,j,1)   =  qintsum(i,j,1)   +  Ab(l,m)  *  qflux(i,j,1,m)

    Dintsum(i,j,1)   =  Dintsum(i,j,1)   +  Ab(l,m)  *  Dflux(i,j,1,m)

! -------------------------------------- EGLM--------------------------------------

    tauintsum(i,j,1) =  tauintsum(i,j,1) +  Ab(l,m)  *  tauflux(i,j,1,m)     &
                                         +  FAC_EGLM * Ab(l,m)  *  taueglm(i,j,1,m)

    Sxintsum(i,j,1)  =  Sxintsum(i,j,1)  +  Ab(l,m)  *  Sxflux(i,j,1,m)      &
                                         +  FAC_EGLM * Ab(l,m)  *  Sxeglm(i,j,1,m)

    Syintsum(i,j,1)  =  Syintsum(i,j,1)  +  Ab(l,m)  *  Syflux(i,j,1,m)      &
                                         +  FAC_EGLM * Ab(l,m)  *  Syeglm(i,j,1,m)

    Szintsum(i,j,1)  =  Szintsum(i,j,1)  +  Ab(l,m)  *  Szflux(i,j,1,m)      &
                                         +  FAC_EGLM * Ab(l,m)  *  Szeglm(i,j,1,m)

    ! -------------------------------------- EGLM--------------------------------------


            end do ! for i
         end do ! for j

!$OMP END DO
!$OMP END PARALLEL

  end subroutine sumint
