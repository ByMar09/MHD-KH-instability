  !     *******************************************************************
  !     Subroutine which calculed  simpler Lax-Friedrich flows 
  !     *******************************************************************

  subroutine laxflow

    use scalar
    use parameters
    use threevectors
    use fourvectors
    use funciones, only: mcl

    implicit none


!$OMP PARALLEL PRIVATE(c_sp,Sp)

    
    if (DIM == 1) then

!$OMP DO 

    do i=-2,imax+2


! Electric Gauss Law flux
!_________________________________________________________________________
       Exlaxpsip(i,1,1,l) =  Exint (i  ,1  ,1  ,l) + psiint(i  ,1  ,1  ,l)
       Exlaxpsim(i,1,1,l) =  Exint (i  ,1  ,1  ,l) - psiint(i  ,1  ,1  ,l)
       Eylaxpsip(i,1,1,l) =  Eyint (i  ,1  ,1  ,l) + psiint(i  ,1  ,1  ,l)
       Eylaxpsim(i,1,1,l) =  Eyint (i  ,1  ,1  ,l) - psiint(i  ,1  ,1  ,l)
       Ezlaxpsip(i,1,1,l) =  Ezint (i  ,1  ,1  ,l) + psiint(i  ,1  ,1  ,l)
       Ezlaxpsim(i,1,1,l) =  Ezint (i  ,1  ,1  ,l) - psiint(i  ,1  ,1  ,l)
! Magnetic Gauss Law flux
!_________________________________________________________________________

       Bxlaxphip(i,1,1,l) =  Bxint (i  ,1  ,1  ,l) + phiint(i  ,1  ,1  ,l) 
       Bxlaxphim(i,1,1,l) =  Bxint (i  ,1  ,1  ,l) - phiint(i  ,1  ,1  ,l) 
       Bylaxphip(i,1,1,l) =  Byint (i  ,1  ,1  ,l) + phiint(i  ,1  ,1  ,l) 
       Bylaxphim(i,1,1,l) =  Byint (i  ,1  ,1  ,l) - phiint(i  ,1  ,1  ,l) 
       Bzlaxphip(i,1,1,l) =  Bzint (i  ,1  ,1  ,l) + phiint(i  ,1  ,1  ,l) 
       Bzlaxphim(i,1,1,l) =  Bzint (i  ,1  ,1  ,l) - phiint(i  ,1  ,1  ,l) 

! Faraday Law flux
!_________________________________________________________________________

! Flujos Bx

       philaxBxp(i,1,1,l) =  phiint(i  ,1  ,1  ,l) + Bxint (i  ,1  ,1  ,l) 
       philaxBxm(i,1,1,l) =  phiint(i  ,1  ,1  ,l) - Bxint (i  ,1  ,1  ,l) 
       EzlaxBxp (i,1,1,l) =  Ezint (i  ,1  ,1  ,l) + Bxint (i  ,1  ,1  ,l) 
       EzlaxBxm (i,1,1,l) =  Ezint (i  ,1  ,1  ,l) - Bxint (i  ,1  ,1  ,l) 
       EylaxBxp (i,1,1,l) = -Eyint (i  ,1  ,1  ,l) + Bxint (i  ,1  ,1  ,l)  
       EylaxBxm (i,1,1,l) = -Eyint (i  ,1  ,1  ,l) - Bxint (i  ,1  ,1  ,l)
! Flujos By

       EzlaxByp (i,1,1,l) = -Ezint (i  ,1  ,1  ,l) + Byint (i  ,1  ,1  ,l) 
       EzlaxBym (i,1,1,l) = -Ezint (i  ,1  ,1  ,l) - Byint (i  ,1  ,1  ,l) 
       philaxByp(i,1,1,l) =  phiint(i  ,1  ,1  ,l) + Byint (i  ,1  ,1  ,l) 
       philaxBym(i,1,1,l) =  phiint(i  ,1  ,1  ,l) - Byint (i  ,1  ,1  ,l) 
       ExlaxByp (i,1,1,l) =  Exint (i  ,1  ,1  ,l) + Byint (i  ,1  ,1  ,l) 
       ExlaxBym (i,1,1,l) =  Exint (i  ,1  ,1  ,l) - Byint (i  ,1  ,1  ,l) 


! Flujos Bz

       EylaxBzp (i,1,1,l) =  Eyint (i  ,1  ,1  ,l) + Bzint (i  ,1  ,1  ,l)   
       EylaxBzm (i,1,1,l) =  Eyint (i  ,1  ,1  ,l) - Bzint (i  ,1  ,1  ,l)   
       ExlaxBzp (i,1,1,l) = -Exint (i  ,1  ,1  ,l) + Bzint (i  ,1  ,1  ,l) 
       ExlaxBzm (i,1,1,l) = -Exint (i  ,1  ,1  ,l) - Bzint (i  ,1  ,1  ,l) 
       philaxBzp(i,1,1,l) =  phiint(i  ,1  ,1  ,l) + Bzint (i  ,1  ,1  ,l)  
       philaxBzm(i,1,1,l) =  phiint(i  ,1  ,1  ,l) - Bzint (i  ,1  ,1  ,l)  


! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxExp(i,1,1,l) =  psiint(i  ,1  ,1  ,l) + Exint (i  ,1  ,1  ,l) 
       psilaxExm(i,1,1,l) =  psiint(i  ,1  ,1  ,l) - Exint (i  ,1  ,1  ,l) 
       BzlaxExp (i,1,1,l) = -Bzint (i  ,1  ,1  ,l) + Exint (i  ,1  ,1  ,l) 
       BzlaxExm (i,1,1,l) = -Bzint (i  ,1  ,1  ,l) - Exint (i  ,1  ,1  ,l) 
       BylaxExp (i,1,1,l) =  Byint (i  ,1  ,1  ,l) + Exint (i  ,1  ,1  ,l) 
       BylaxExm (i,1,1,l) =  Byint (i  ,1  ,1  ,l) - Exint (i  ,1  ,1  ,l) 

! Flujos Ey

       BzlaxEyp (i,1,1,l) =  Bzint (i  ,1  ,1  ,l) + Eyint (i  ,1  ,1  ,l) 
       BzlaxEym (i,1,1,l) =  Bzint (i  ,1  ,1  ,l) - Eyint (i  ,1  ,1  ,l) 
       psilaxEyp(i,1,1,l) =  psiint(i  ,1  ,1  ,l) + Eyint (i  ,1  ,1  ,l) 
       psilaxEym(i,1,1,l) =  psiint(i  ,1  ,1  ,l) - Eyint (i  ,1  ,1  ,l) 
       BxlaxEyp (i,1,1,l) = -Bxint (i  ,1  ,1  ,l) + Eyint (i  ,1  ,1  ,l) 
       BxlaxEym (i,1,1,l) = -Bxint (i  ,1  ,1  ,l) - Eyint (i  ,1  ,1  ,l) 

! Flujos Ez

       BylaxEzp (i,1,1,l) = -Byint (i  ,1  ,1  ,l) + Ezint (i  ,1  ,1  ,l) 
       BylaxEzm (i,1,1,l) = -Byint (i  ,1  ,1  ,l) - Ezint (i  ,1  ,1  ,l)  
       BxlaxEzp (i,1,1,l) =  Bxint (i  ,1  ,1  ,l) + Ezint (i  ,1  ,1  ,l) 
       BxlaxEzm (i,1,1,l) =  Bxint (i  ,1  ,1  ,l) - Ezint (i  ,1  ,1  ,l) 
       psilaxEzp(i,1,1,l) =  psiint(i  ,1  ,1  ,l) + Ezint (i  ,1  ,1  ,l) 
       psilaxEzm(i,1,1,l) =  psiint(i  ,1  ,1  ,l) - Ezint (i  ,1  ,1  ,l) 

! Conserved current flux
!_________________________________________________________________________

       Jxlaxp(i,1,1,l)    =  Jxint(i  ,1  ,1  ,l)  + qint (i  ,1  ,1  ,l) 
       Jxlaxm(i,1,1,l)    =  Jxint(i  ,1  ,1  ,l)  - qint (i  ,1  ,1  ,l) 
       Jylaxp(i,1,1,l)    =  Jyint(i  ,1  ,1  ,l)  + qint (i  ,1  ,1  ,l) 
       Jylaxm(i,1,1,l)    =  Jyint(i  ,1  ,1  ,l)  - qint (i  ,1  ,1  ,l) 
       Jzlaxp(i,1,1,l)    =  Jzint(i  ,1  ,1  ,l)  + qint (i  ,1  ,1  ,l) 
       Jzlaxm(i,1,1,l)    =  Jzint(i  ,1  ,1  ,l)  - qint (i  ,1  ,1  ,l) 

! Conserved Mass flux
!_________________________________________________________________________

       FDxlaxp(i,1,1,l)   =  FDxint(i  ,1  ,1  ,l) + DDint (i  ,1  ,1  ,l)  
       FDxlaxm(i,1,1,l)   =  FDxint(i  ,1  ,1  ,l) - DDint (i  ,1  ,1  ,l)  
       FDylaxp(i,1,1,l)   =  FDyint(i  ,1  ,1  ,l) + DDint (i  ,1  ,1  ,l) 
       FDylaxm(i,1,1,l)   =  FDyint(i  ,1  ,1  ,l) - DDint (i  ,1  ,1  ,l) 
       FDzlaxp(i,1,1,l)   =  FDzint(i  ,1  ,1  ,l) + DDint (i  ,1  ,1  ,l) 
       FDzlaxm(i,1,1,l)   =  FDzint(i  ,1  ,1  ,l) - DDint (i  ,1  ,1  ,l) 

! Conserved Energy flux
!_________________________________________________________________________

       Ftauxlaxp(i,1,1,l) =  Ftauxint(i  ,1  ,1  ,l) + tauint  (i  ,1  ,1  ,l) 
       Ftauxlaxm(i,1,1,l) =  Ftauxint(i  ,1  ,1  ,l) - tauint  (i  ,1  ,1  ,l) 
       Ftauylaxp(i,1,1,l) =  Ftauyint(i  ,1  ,1  ,l) + tauint  (i  ,1  ,1  ,l) 
       Ftauylaxm(i,1,1,l) =  Ftauyint(i  ,1  ,1  ,l) - tauint  (i  ,1  ,1  ,l) 
       Ftauzlaxp(i,1,1,l) =  Ftauzint(i  ,1  ,1  ,l) + tauint  (i  ,1  ,1  ,l) 
       Ftauzlaxm(i,1,1,l) =  Ftauzint(i  ,1  ,1  ,l) - tauint  (i  ,1  ,1  ,l) 

! -------------------------------------- EGLM--------------------------------------
       psilaxtaup(i,1,1,l) = psiint  (i  ,1  ,1  ,l) + tauint  (i  ,1  ,1  ,l) 
       psilaxtaum(i,1,1,l) = psiint  (i  ,1  ,1  ,l) - tauint  (i  ,1  ,1  ,l) 
       philaxtaup(i,1,1,l) = phiint  (i  ,1  ,1  ,l) + tauint  (i  ,1  ,1  ,l)
       philaxtaum(i,1,1,l) = phiint  (i  ,1  ,1  ,l) - tauint  (i  ,1  ,1  ,l)  
! -------------------------------------- EGLM--------------------------------------

! Conserved Momentum flux tensor
!_________________________________________________________________________

       FSxxlaxp(i,1,1,l)  =  FSxxint(i  ,1  ,1  ,l) + Sxint (i  ,1  ,1  ,l) 
       FSxxlaxm(i,1,1,l)  =  FSxxint(i  ,1  ,1  ,l) - Sxint (i  ,1  ,1  ,l) 
       FSxylaxp(i,1,1,l)  =  FSxyint(i  ,1  ,1  ,l) + Syint (i  ,1  ,1  ,l)
       FSxylaxm(i,1,1,l)  =  FSxyint(i  ,1  ,1  ,l) - Syint (i  ,1  ,1  ,l)
       FSxzlaxp(i,1,1,l)  =  FSxzint(i  ,1  ,1  ,l) + Szint (i  ,1  ,1  ,l) 
       FSxzlaxm(i,1,1,l)  =  FSxzint(i  ,1  ,1  ,l) - Szint (i  ,1  ,1  ,l) 

       FSyxlaxp(i,1,1,l)  =  FSyxint(i  ,1  ,1  ,l) + Sxint (i  ,1  ,1  ,l) 
       FSyxlaxm(i,1,1,l)  =  FSyxint(i  ,1  ,1  ,l) - Sxint (i  ,1  ,1  ,l) 
       FSyylaxp(i,1,1,l)  =  FSyyint(i  ,1  ,1  ,l) + Syint (i  ,1  ,1  ,l) 
       FSyylaxm(i,1,1,l)  =  FSyyint(i  ,1  ,1  ,l) - Syint (i  ,1  ,1  ,l) 
       FSyzlaxp(i,1,1,l)  =  FSyzint(i  ,1  ,1  ,l) + Szint (i  ,1  ,1  ,l) 
       FSyzlaxm(i,1,1,l)  =  FSyzint(i  ,1  ,1  ,l) - Szint (i  ,1  ,1  ,l) 

       FSzxlaxp(i,1,1,l)  =  FSzxint(i  ,1  ,1  ,l) + Sxint (i  ,1  ,1  ,l) 
       FSzxlaxm(i,1,1,l)  =  FSzxint(i  ,1  ,1  ,l) - Sxint (i  ,1  ,1  ,l) 
       FSzylaxp(i,1,1,l)  =  FSzyint(i  ,1  ,1  ,l) + Syint (i  ,1  ,1  ,l) 
       FSzylaxm(i,1,1,l)  =  FSzyint(i  ,1  ,1  ,l) - Syint (i  ,1  ,1  ,l) 
       FSzzlaxp(i,1,1,l)  =  FSzzint(i  ,1  ,1  ,l) + Szint (i  ,1  ,1  ,l)   
       FSzzlaxm(i,1,1,l)  =  FSzzint(i  ,1  ,1  ,l) - Szint (i  ,1  ,1  ,l) 

! -------------------------------------- EGLM--------------------------------------

       BxlaxSxp(i,1,1,l)   =  Bxint  (i  ,1  ,1  ,l) + Sxint (i  ,1  ,1  ,l)
       BxlaxSxm(i,1,1,l)   =  Bxint  (i  ,1  ,1  ,l) - Sxint (i  ,1  ,1  ,l)
       BylaxSxp(i,1,1,l)   =  Byint  (i  ,1  ,1  ,l) + Sxint (i  ,1  ,1  ,l)
       BylaxSxm(i,1,1,l)   =  Byint  (i  ,1  ,1  ,l) - Sxint (i  ,1  ,1  ,l)
       BzlaxSxp(i,1,1,l)   =  Bzint  (i  ,1  ,1  ,l) + Sxint (i  ,1  ,1  ,l)
       BzlaxSxm(i,1,1,l)   =  Bzint  (i  ,1  ,1  ,l) - Sxint (i  ,1  ,1  ,l)

       BxlaxSyp(i,1,1,l)   =  Bxint  (i  ,1  ,1  ,l) + Syint (i  ,1  ,1  ,l)
       BxlaxSym(i,1,1,l)   =  Bxint  (i  ,1  ,1  ,l) - Syint (i  ,1  ,1  ,l)
       BylaxSyp(i,1,1,l)   =  Byint  (i  ,1  ,1  ,l) + Syint (i  ,1  ,1  ,l)
       BylaxSym(i,1,1,l)   =  Byint  (i  ,1  ,1  ,l) - Syint (i  ,1  ,1  ,l)
       BzlaxSyp(i,1,1,l)   =  Bzint  (i  ,1  ,1  ,l) + Syint (i  ,1  ,1  ,l)
       BzlaxSym(i,1,1,l)   =  Bzint  (i  ,1  ,1  ,l) - Syint (i  ,1  ,1  ,l)

       BxlaxSzp(i,1,1,l)   =  Bxint  (i  ,1  ,1  ,l) + Szint (i  ,1  ,1  ,l)
       BxlaxSzm(i,1,1,l)   =  Bxint  (i  ,1  ,1  ,l) - Szint (i  ,1  ,1  ,l)
       BylaxSzp(i,1,1,l)   =  Byint  (i  ,1  ,1  ,l) + Szint (i  ,1  ,1  ,l)
       BylaxSzm(i,1,1,l)   =  Byint  (i  ,1  ,1  ,l) - Szint (i  ,1  ,1  ,l)
       BzlaxSzp(i,1,1,l)   =  Bzint  (i  ,1  ,1  ,l) + Szint (i  ,1  ,1  ,l)
       BzlaxSzm(i,1,1,l)   =  Bzint  (i  ,1  ,1  ,l) - Szint (i  ,1  ,1  ,l)

! -------------------------------------- EGLM--------------------------------------


   end do ! for i

!$OMP END DO

!$OMP   PARALLEL        
!$OMP DO 

         do i=-4,imax+4

            u(i, 1,1) =  Exlaxpsip(i,1,1,l)
            u(i, 2,1) =  Exlaxpsim(i,1,1,l)
            u(i, 3,1) =  Eylaxpsip(i,1,1,l)
            u(i, 4,1) =  Eylaxpsim(i,1,1,l)
            u(i, 5,1) =  Ezlaxpsip(i,1,1,l)
            u(i, 6,1) =  Ezlaxpsim(i,1,1,l)
            
            u(i, 7,1) =  Bxlaxphip(i,1,1,l)
            u(i, 8,1) =  Bxlaxphim(i,1,1,l)Ezint (i,1,1,l)
            u(i, 9,1) =  Bylaxphip(i,1,1,l)Bxint (i,1,1,l)
            u(i,10,1) =  Bylaxphim(i,1,1,l)Byint (i,1,1,l)
            u(i,11,1) =  Bzlaxphip(i,1,1,l)Bzint (i,1,1,l)
            u(i,12,1) =  Bzlaxphim(i,1,1,l)Bxint (i,1,1,l)

       Bxlaxphip(i,1,1,l) =  Bxint (i  ,1  ,1  ,l) + phiint(i  ,1  ,1  ,l) 
       Bxlaxphim(i,1,1,l) =  Bxint (i  ,1  ,1  ,l) - phiint(i  ,1  ,1  ,l) 
       Bylaxphip(i,1,1,l) =  Byint (i  ,1  ,1  ,l) + phiint(i  ,1  ,1  ,l) 
       Bylaxphim(i,1,1,l) =  Byint (i  ,1  ,1  ,l) - phiint(i  ,1  ,1  ,l) 
       Bzlaxphip(i,1,1,l) =  Bzint (i  ,1  ,1  ,l) + phiint(i  ,1  ,1  ,l) 
       Bzlaxphim(i,1,1,l) =  Bzint (i  ,1  ,1  ,l) - phiint(i  ,1  ,1  ,l) 



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
!$OMP END DO
!$OMP END PARALLEL
   

!$OMP DO 
         
    do i=-1,imax


! Electric Gauss Law flux
!_________________________________________________________________________

       ExlaxpsiL(i,1,1,l) =     Exlaxpsip(i  ,1,1,l) + FAC  *                  &
                            mcl(Exlaxpsip(i+1,1,1,l) - Exlaxpsip(i  ,1,1,l) ,  &
                                Exlaxpsip(i  ,1,1,l) - Exlaxpsip(i-1,1,1,l))
       ExlaxpsiR(i,1,1,l) =     Exlaxpsim(i+1,1,1,l) - FAC  *                  &
                            mcl(Exlaxpsim(i+2,1,1,l) - Exlaxpsim(i+1,1,1,l) ,  &
                                Exlaxpsim(i+1,1,1,l) - Exlaxpsim(i  ,1,1,l))

       EylaxpsiL(i,1,1,l) =     Eylaxpsip(i  ,1,1,l) + FAC  *                  &
                            mcl(Eylaxpsip(i  ,1,1,l) - Eylaxpsip(i  ,1,1,l) ,  &
                                Eylaxpsip(i  ,1,1,l) - Eylaxpsip(i  ,1,1,l))
       EylaxpsiR(i,1,1,l) =     Eylaxpsim(i  ,1,1,l) - FAC  *                  &
                            mcl(Eylaxpsim(i  ,1,1,l) - Eylaxpsim(i  ,1,1,l) ,  &
                                Eylaxpsim(i  ,1,1,l) - Eylaxpsim(i  ,1,1,l))

       EzlaxpsiL(i,1,1,l) =     Ezlaxpsip(i  ,1,1,l) + FAC  *                  &
                            mcl(Ezlaxpsip(i  ,1,1,l) - Ezlaxpsip(i  ,1,1,l) ,  &
                                Ezlaxpsip(i  ,1,1,l) - Ezlaxpsip(i  ,1,1,l))
       EzlaxpsiR(i,1,1,l) =     Ezlaxpsim(i  ,1,1,l) - FAC  *                  &
                            mcl(Ezlaxpsim(i  ,1,1,l) - Ezlaxpsim(i  ,1,1,l) ,  &
                                Ezlaxpsim(i  ,1,1,l) - Ezlaxpsim(i  ,1,1,l))

       Exlaxpsi(i,1,1,l)  =   0.5d0 * (ExlaxpsiL(i,1,1,l) + ExlaxpsiR(i,1,1,l))  
       Eylaxpsi(i,1,1,l)  =   0.5d0 * (EylaxpsiL(i,1,1,l) + EylaxpsiR(i,1,1,l))
       Ezlaxpsi(i,1,1,l)  =   0.5d0 * (EzlaxpsiL(i,1,1,l) + EzlaxpsiR(i,1,1,l))

! Magnetic Gauss Law flux
!_________________________________________________________________________

       BxlaxphiL(i,1,1,l) =     Bxlaxphip(i  ,1,1,l) + FAC  *                  &
                            mcl(Bxlaxphip(i+1,1,1,l) - Bxlaxphip(i  ,1,1,l) ,  &
                                Bxlaxphip(i  ,1,1,l) - Bxlaxphip(i-1,1,1,l))
       BxlaxphiR(i,1,1,l) =     Bxlaxphim(i+1,1,1,l) - FAC  *                  &
                            mcl(Bxlaxphim(i+2,1,1,l) - Bxlaxphim(i+1,1,1,l) ,  &
                                Bxlaxphim(i+1,1,1,l) - Bxlaxphim(i  ,1,1,l))

       BylaxphiL(i,1,1,l) =     Bylaxphip(i  ,1,1,l) + FAC  *                  &
                            mcl(Bylaxphip(i  ,1,1,l) - Bylaxphip(i  ,1,1,l) ,  &
                                Bylaxphip(i  ,1,1,l) - Bylaxphip(i  ,1,1,l))
       BylaxphiR(i,1,1,l) =     Bylaxphim(i  ,1,1,l) - FAC  *                  &
                            mcl(Bylaxphim(i  ,1,1,l) - Bylaxphim(i  ,1,1,l) ,  &
                                Bylaxphim(i  ,1,1,l) - Bylaxphim(i  ,1,1,l))

       BzlaxphiL(i,1,1,l) =     Bzlaxphip(i  ,1,1,l) + FAC  *                  &
                            mcl(Bzlaxphip(i  ,1,1,l) - Bzlaxphip(i  ,1,1,l) ,  &
                                Bzlaxphip(i  ,1,1,l) - Bzlaxphip(i  ,1,1,l))
       BzlaxphiR(i,1,1,l) =     Bzlaxphim(i  ,1,1,l) - FAC  *                  &
                            mcl(Bzlaxphim(i  ,1,1,l) - Bzlaxphim(i  ,1,1,l) ,  &
                                Bzlaxphim(i  ,1,1,l) - Bzlaxphim(i  ,1,1,l))

       Bxlaxphi(i,1,1,l)  =   0.5d0 * (BxlaxphiL(i,1,1,l) + BxlaxphiR(i,1,1,l))
       Bylaxphi(i,1,1,l)  =   0.5d0 * (BylaxphiL(i,1,1,l) + BylaxphiR(i,1,1,l))
       Bzlaxphi(i,1,1,l)  =   0.5d0 * (BzlaxphiL(i,1,1,l) + BzlaxphiR(i,1,1,l))

! Faraday Law flux
!_________________________________________________________________________

! Flujos Bx 

       philaxBxL(i,1,1,l) =     philaxBxp(i  ,1,1,l) + FAC  *                 &
                            mcl(philaxBxp(i+1,1,1,l) - philaxBxp(i  ,1,1,l),  &
                                philaxBxp(i  ,1,1,l) - philaxBxp(i-1,1,1,l))
       philaxBxR(i,1,1,l) =     philaxBxm(i+1,1,1,l) - FAC  *                 &
                            mcl(philaxBxm(i+2,1,1,l) - philaxBxm(i+1,1,1,l),  &
                                philaxBxm(i+1,1,1,l) - philaxBxm(i  ,1,1,l)) 

       EzlaxBxL(i,1,1,l)  =     EzlaxBxp(i  ,1,1,l)  + FAC  *                 &
                            mcl(EzlaxBxp(i  ,1,1,l) - EzlaxBxp(i  ,1,1,l),    &
                                EzlaxBxp(i  ,1,1,l) - EzlaxBxp(i  ,1,1,l))
       EzlaxBxR(i,1,1,l)  =     EzlaxBxm(i  ,1,1,l) - FAC  *                  &
                            mcl(EzlaxBxm(i  ,1,1,l) - EzlaxBxm(i  ,1,1,l),    &
                                EzlaxBxm(i  ,1,1,l) - EzlaxBxm(i  ,1,1,l))  

       EylaxBxL(i,1,1,l)  =     EylaxBxp(i  ,1,1,l) + FAC  *                  &
                            mcl(EylaxBxp(i  ,1,1,l) - EylaxBxp(i  ,1,1,l),    &
                                EylaxBxp(i  ,1,1,l) - EylaxBxp(i  ,1,1,l))
       EylaxBxR(i,1,1,l)  =     EylaxBxm(i  ,1,1,l) - FAC  *                  &
                            mcl(EylaxBxm(i  ,1,1,l) - EylaxBxm(i  ,1,1,l),    &
                                EylaxBxm(i  ,1,1,l) - EylaxBxm(i  ,1,1,l))   


       philaxBx(i,1,1,l)  =   0.5d0 * (philaxBxL(i,1,1,l) + philaxBxR(i,1,1,l))
       EzlaxBx (i,1,1,l)  =   0.5d0 * (EzlaxBxL (i,1,1,l) + EzlaxBxR (i,1,1,l))
       EylaxBx (i,1,1,l)  = - 0.5d0 * (EylaxBxL (i,1,1,l) + EylaxBxR (i,1,1,l))   

! Flujos By

       EzlaxByL(i,1,1,l)  =     EzlaxByp(i  ,1,1,l) + FAC  *                  &
                            mcl(EzlaxByp(i+1,1,1,l) - EzlaxByp(i  ,1,1,l),    &
                                EzlaxByp(i  ,1,1,l) - EzlaxByp(i-1,1,1,l))
       EzlaxByR(i,1,1,l)  =     EzlaxBym(i+1,1,1,l) - FAC  *                  &
                            mcl(EzlaxBym(i+2,1,1,l) - EzlaxBym(i+1,1,1,l),    &
                                EzlaxBym(i+1,1,1,l) - EzlaxBym(i  ,1,1,l))

       philaxByL(i,1,1,l) =     philaxByp(i  ,1,1,l) + FAC  *                 &
                            mcl(philaxByp(i  ,1,1,l) - philaxByp(i  ,1,1,l),  &
                                philaxByp(i  ,1,1,l) - philaxByp(i  ,1,1,l))
       philaxByR(i,1,1,l) =     philaxBym(i  ,1,1,l) - FAC  *                 &
                            mcl(philaxBym(i  ,1,1,l) - philaxBym(i  ,1,1,l),  &
                                philaxBym(i  ,1,1,l) - philaxBym(i  ,1,1,l))

       ExlaxByL(i,1,1,l)  =     ExlaxByp(i  ,1,1,l) + FAC  *                  &
                            mcl(ExlaxByp(i  ,1,1,l) - ExlaxByp(i  ,1,1,l),    &
                                ExlaxByp(i  ,1,1,l) - ExlaxByp(i  ,1,1,l))
       ExlaxByR(i,1,1,l)  =     ExlaxBym(i  ,1,1,l) - FAC  *                  &
                            mcl(ExlaxBym(i  ,1,1,l) - ExlaxBym(i  ,1,1,l),    &
                                ExlaxBym(i  ,1,1,l) - ExlaxBym(i  ,1,1,l))

       EzlaxBy (i,1,1,l)  = - 0.5d0 * (EzlaxByL (i,1,1,l) + EzlaxByR (i,1,1,l))
       philaxBy(i,1,1,l)  =   0.5d0 * (philaxByL(i,1,1,l) + philaxByR(i,1,1,l)) 
       ExlaxBy (i,1,1,l)  =   0.5d0 * (ExlaxByL (i,1,1,l) + ExlaxByR (i,1,1,l))  

! Flujos Bz

       EylaxBzL(i,1,1,l)  =     EylaxBzp(i  ,1,1,l) + FAC  *                  &
                            mcl(EylaxBzp(i+1,1,1,l) - EylaxBzp(i  ,1,1,l),    &
                                EylaxBzp(i  ,1,1,l) - EylaxBzp(i-1,1,1,l))  
       EylaxBzR(i,1,1,l)  =     EylaxBzm(i+1,1,1,l) - FAC  *                  &
                            mcl(EylaxBzm(i+2,1,1,l) - EylaxBzm(i+1,1,1,l),    &
                                EylaxBzm(i+1,1,1,l) - EylaxBzm(i  ,1,1,l))

       ExlaxBzL(i,1,1,l)  =     ExlaxBzp(i  ,1,1,l) + FAC  *                  &
                            mcl(ExlaxBzp(i  ,1,1,l) - ExlaxBzp(i  ,1,1,l),    &
                                ExlaxBzp(i  ,1,1,l) - ExlaxBzp(i  ,1,1,l))
       ExlaxBzR(i,1,1,l)  =     ExlaxBzm(i  ,1,1,l) - FAC  *                  &
                            mcl(ExlaxBzm(i  ,1,1,l) - ExlaxBzm(i  ,1,1,l),    &
                                ExlaxBzm(i  ,1,1,l) - ExlaxBzm(i  ,1,1,l))

       philaxBzL(i,1,1,l) =     philaxBzp(i  ,1,1,l) + FAC  *                 &
                            mcl(philaxBzp(i  ,1,1,l) - philaxBzp(i  ,1,1,l),  &
                                philaxBzp(i  ,1,1,l) - philaxBzp(i  ,1,1,l))
       philaxBzR(i,1,1,l) =     philaxBzm(i  ,1,1,l) - FAC  *                 &
                            mcl(philaxBzm(i  ,1,1,l) - philaxBzm(i  ,1,1,l),  &
                                philaxBzm(i  ,1,1,l) - philaxBzm(i  ,1,1,l))

       EylaxBz (i,1,1,l)  =   0.5d0 * (EylaxBzL (i,1,1,l)  + EylaxBzR (i,1,1,l))
       ExlaxBz (i,1,1,l)  = - 0.5d0 * (ExlaxBzL (i,1,1,l)  + ExlaxBzR (i,1,1,l))
       philaxBz(i,1,1,l)  =   0.5d0 * (philaxBzL(i,1,1,l)  + philaxBzR(i,1,1,l))


! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxExL(i,1,1,l) =     psilaxExp(i  ,1,1,l) + FAC  *                 &
                            mcl(psilaxExp(i+1,1,1,l) - psilaxExp(i  ,1,1,l),  &
                                psilaxExp(i  ,1,1,l) - psilaxExp(i-1,1,1,l))  
       psilaxExR(i,1,1,l) =     psilaxExm(i+1,1,1,l) - FAC  *                 &
                            mcl(psilaxExm(i+2,1,1,l) - psilaxExm(i+1,1,1,l),  &
                                psilaxExm(i+1,1,1,l) - psilaxExm(i  ,1,1,l))

       BzlaxExL(i,1,1,l)  =     BzlaxExp(i  ,1,1,l) + FAC  *                  &
                            mcl(BzlaxExp(i  ,1,1,l) - BzlaxExp(i  ,1,1,l),    &
                                BzlaxExp(i  ,1,1,l) - BzlaxExp(i  ,1,1,l))
       BzlaxExR(i,1,1,l)  =     BzlaxExm(i  ,1,1,l) - FAC  *                  &
                            mcl(BzlaxExm(i  ,1,1,l) - BzlaxExm(i  ,1,1,l),    &
                                BzlaxExm(i  ,1,1,l) - BzlaxExm(i  ,1,1,l))

       BylaxExL(i,1,1,l)  =     BylaxExp(i  ,1,1,l) + FAC  *                  &
                            mcl(BylaxExp(i  ,1,1,l) - BylaxExp(i  ,1,1,l),    &
                                BylaxExp(i  ,1,1,l) - BylaxExp(i  ,1,1,l))
       BylaxExR(i,1,1,l)  =     BylaxExm(i  ,1,1,l) - FAC  *                  &
                            mcl(BylaxExm(i  ,1,1,l) - BylaxExm(i  ,1,1,l),    &
                                BylaxExm(i  ,1,1,l) - BylaxExm(i  ,1,1,l))

       psilaxEx(i,1,1,l)  =   0.5d0 * (psilaxExL(i,1,1,l) + psilaxExR(i,1,1,l))
       BzlaxEx(i,1,1,l)   = - 0.5d0 * (BzlaxExL (i,1,1,l) + BzlaxExR (i,1,1,l))  
       BylaxEx(i,1,1,l)   =   0.5d0 * (BylaxExL (i,1,1,l) + BylaxExR (i,1,1,l))  

! Flujos Ey

       BzlaxEyL(i,1,1,l)  =     BzlaxEyp(i  ,1,1,l) + FAC  *                  &
                            mcl(BzlaxEyp(i+1,1,1,l) - BzlaxEyp(i  ,1,1,l),    &
                                BzlaxEyp(i  ,1,1,l) - BzlaxEyp(i-1,1,1,l))  
       BzlaxEyR(i,1,1,l)  =     BzlaxEym(i+1,1,1,l) - FAC  *                  &
                            mcl(BzlaxEym(i+2,1,1,l) - BzlaxEym(i+1,1,1,l),    &
                                BzlaxEym(i+1,1,1,l) - BzlaxEym(i  ,1,1,l))

       psilaxEyL(i,1,1,l) =     psilaxEyp(i  ,1,1,l) + FAC  *                 &
                            mcl(psilaxEyp(i  ,1,1,l) - psilaxEyp(i  ,1,1,l),  &
                                psilaxEyp(i  ,1,1,l) - psilaxEyp(i  ,1,1,l))
       psilaxEyR(i,1,1,l) =     psilaxEym(i  ,1,1,l) - FAC  *                 &
                            mcl(psilaxEym(i  ,1,1,l) - psilaxEym(i  ,1,1,l),  &
                                psilaxEym(i  ,1,1,l) - psilaxEym(i  ,1,1,l))

       BxlaxEyL(i,1,1,l)  =     BxlaxEyp(i  ,1,1,l) + FAC  *                 &
                            mcl(BxlaxEyp(i  ,1,1,l) - BxlaxEyp(i  ,1,1,l),   &
                                BxlaxEyp(i  ,1,1,l) - BxlaxEyp(i  ,1,1,l))
       BxlaxEyR(i,1,1,l)  =     BxlaxEym(i  ,1,1,l) - FAC  *                 &
                            mcl(BxlaxEym(i  ,1,1,l) - BxlaxEym(i  ,1,1,l),   &
                                BxlaxEym(i  ,1,1,l) - BxlaxEym(i  ,1,1,l))

       BzlaxEy (i,1,1,l)  =   0.5d0 * (BzlaxEyL (i,1,1,l) + BzlaxEyR (i,1,1,l))
       psilaxEy(i,1,1,l)  =   0.5d0 * (psilaxEyL(i,1,1,l) + psilaxEyR(i,1,1,l))
       BxlaxEy (i,1,1,l)  = - 0.5d0 * (BxlaxEyL (i,1,1,l) + BxlaxEyR (i,1,1,l))

! Flujos Ez

       BylaxEzL(i,1,1,l)  =     BylaxEzp(i  ,1,1,l) + FAC  *                  &
                            mcl(BylaxEzp(i+1,1,1,l) - BylaxEzp(i  ,1,1,l),    &
                                BylaxEzp(i  ,1,1,l) - BylaxEzp(i-1,1,1,l)) 
       BylaxEzR(i,1,1,l)  =     BylaxEzm(i+1,1,1,l) -  FAC  *                 &
                            mcl(BylaxEzm(i+2,1,1,l) - BylaxEzm(i+1,1,1,l),    &
                                BylaxEzm(i+1,1,1,l) - BylaxEzm(i  ,1,1,l))

       BxlaxEzL(i,1,1,l)  =     BxlaxEzp(i  ,1,1,l) + FAC  *                  &
                            mcl(BxlaxEzp(i  ,1,1,l) - BxlaxEzp(i  ,1,1,l),    &
                                BxlaxEzp(i  ,1,1,l) - BxlaxEzp(i  ,1,1,l))
       BxlaxEzR(i,1,1,l)  =     BxlaxEzm(i  ,1,1,l) - FAC  *                  &
                            mcl(BxlaxEzm(i  ,1,1,l) - BxlaxEzm(i  ,1,1,l),    &
                                BxlaxEzm(i  ,1,1,l) - BxlaxEzm(i  ,1,1,l))

       psilaxEzL(i,1,1,l) =     psilaxEzp(i  ,1,1,l) + FAC  *                 &
                            mcl(psilaxEzp(i  ,1,1,l) - psilaxEzp(i  ,1,1,l),  &
                                psilaxEzp(i  ,1,1,l) - psilaxEzp(i  ,1,1,l))
       psilaxEzR(i,1,1,l) =     psilaxEzm(i  ,1,1,l) - FAC  *                 &
                            mcl(psilaxEzm(i  ,1,1,l) - psilaxEzm(i  ,1,1,l),  &
                                psilaxEzm(i  ,1,1,l) - psilaxEzm(i  ,1,1,l))

       BylaxEz (i,1,1,l)  = - 0.5d0 * (BylaxEzL (i,1,1,l) + BylaxEzR (i,1,1,l))
       BxlaxEz (i,1,1,l)  =   0.5d0 * (BxlaxEzL (i,1,1,l) + BxlaxEzR (i,1,1,l))
       psilaxEz(i,1,1,l)  =   0.5d0 * (psilaxEzL(i,1,1,l) + psilaxEzR(i,1,1,l))

! Conserved current flux
!_________________________________________________________________________

       JxlaxL(i,1,1,l)    =     Jxlaxp(i  ,1,1,l) + FAC  *                    &
                            mcl(Jxlaxp(i+1,1,1,l) - Jxlaxp(i  ,1,1,l),        &
                                Jxlaxp(i  ,1,1,l) - Jxlaxp(i-1,1,1,l))  
       JxlaxR(i,1,1,l)    =     Jxlaxm(i+1,1,1,l) -  FAC  *                   &
                            mcl(Jxlaxm(i+2,1,1,l) - Jxlaxm(i+1,1,1,l),        &
                                Jxlaxm(i+1,1,1,l) - Jxlaxm(i  ,1,1,l))

       JylaxL(i,1,1,l)    =     Jylaxp(i  ,1,1,l) +  FAC  *                   &
                            mcl(Jylaxp(i  ,1,1,l) - Jylaxp(i  ,1,1,l),        &
                                Jylaxp(i  ,1,1,l) - Jylaxp(i  ,1,1,l))
       JylaxR(i,1,1,l)    =     Jylaxm(i  ,1,1,l) -  FAC  *                   &
                            mcl(Jylaxm(i  ,1,1,l) - Jylaxm(i  ,1,1,l),        &
                                Jylaxm(i  ,1,1,l) - Jylaxm(i  ,1,1,l))

       JzlaxL(i,1,1,l)    =     Jzlaxp(i  ,1,1,l) +  FAC  *                   &
                            mcl(Jzlaxp(i  ,1,1,l) - Jzlaxp(i  ,1,1,l),        &
                                Jzlaxp(i  ,1,1,l) - Jzlaxp(i  ,1,1,l))
       JzlaxR(i,1,1,l)    =     Jzlaxm(i  ,1,1,l) -  FAC  *                   &
                            mcl(Jzlaxm(i  ,1,1,l) - Jzlaxm(i  ,1,1,l),        &
                                Jzlaxm(i  ,1,1,l) - Jzlaxm(i  ,1,1,l))

       Jxlax(i,1,1,l)     =   0.5d0 * (JxlaxL(i,1,1,l) + JxlaxR(i,1,1,l))
       Jylax(i,1,1,l)     =   0.5d0 * (JylaxL(i,1,1,l) + JylaxR(i,1,1,l))
       Jzlax(i,1,1,l)     =   0.5d0 * (JzlaxL(i,1,1,l) + JzlaxR(i,1,1,l))

! Conserved Mass flux
!_________________________________________________________________________

       FDxlaxL(i,1,1,l)   =     FDxlaxp(i  ,1,1,l) + FAC  *                   &
                            mcl(FDxlaxp(i+1,1,1,l) - FDxlaxp(i  ,1,1,l),      &
                                FDxlaxp(i  ,1,1,l) - FDxlaxp(i-1,1,1,l))   
       FDxlaxR(i,1,1,l)   =     FDxlaxm(i+1,1,1,l) -  FAC  *                  &
                            mcl(FDxlaxm(i+2,1,1,l) - FDxlaxm(i+1,1,1,l),      &
                                FDxlaxm(i+1,1,1,l) - FDxlaxm(i  ,1,1,l))    

       FDylaxL(i,1,1,l)   =     FDylaxp(i  ,1,1,l) +  FAC  *                  &
                            mcl(FDylaxp(i  ,1,1,l) - FDylaxp(i  ,1,1,l),      &
                                FDylaxp(i  ,1,1,l) - FDylaxp(i  ,1,1,l))
       FDylaxR(i,1,1,l)   =     FDylaxm(i  ,1,1,l) -  FAC  *                  &
                            mcl(FDylaxm(i  ,1,1,l) - FDylaxm(i  ,1,1,l),      &
                                FDylaxm(i  ,1,1,l) - FDylaxm(i  ,1,1,l))

       FDzlaxL(i,1,1,l)   =     FDzlaxp(i  ,1,1,l) +  FAC  *                  &
                            mcl(FDzlaxp(i  ,1,1,l) - FDzlaxp(i  ,1,1,l),      &
                                FDzlaxp(i  ,1,1,l) - FDzlaxp(i  ,1,1,l))
       FDzlaxR(i,1,1,l)   =     FDzlaxm(i  ,1,1,l) -  FAC  *                  &
                            mcl(FDzlaxm(i  ,1,1,l) - FDzlaxm(i  ,1,1,l),      &
                                FDzlaxm(i  ,1,1,l) - FDzlaxm(i  ,1,1,l))

       FDxlax(i,1,1,l)    =   0.5d0 * (FDxlaxL(i,1,1,l) + FDxlaxR(i,1,1,l))
       FDylax(i,1,1,l)    =   0.5d0 * (FDylaxL(i,1,1,l) + FDylaxR(i,1,1,l))
       FDzlax(i,1,1,l)    =   0.5d0 * (FDzlaxL(i,1,1,l) + FDzlaxR(i,1,1,l))

! Conserved Energy flux
!_________________________________________________________________________

       FtauxlaxL(i,1,1,l) =     Ftauxlaxp(i  ,1,1,l) + FAC *                  &
                            mcl(Ftauxlaxp(i+1,1,1,l) - Ftauxlaxp(i  ,1,1,l),  &
                                Ftauxlaxp(i  ,1,1,l) - Ftauxlaxp(i-1,1,1,l)) 
       FtauxlaxR(i,1,1,l) =     Ftauxlaxm(i+1,1,1,l) - FAC *                  &
                            mcl(Ftauxlaxm(i+2,1,1,l) - Ftauxlaxm(i+1,1,1,l),  &
                                Ftauxlaxm(i+1,1,1,l) - Ftauxlaxm(i  ,1,1,l)) 

       FtauylaxL(i,1,1,l) =     Ftauylaxp(i  ,1,1,l) +  FAC  *                &
                            mcl(Ftauylaxp(i  ,1,1,l) - Ftauylaxp(i  ,1,1,l),  &
                                Ftauylaxp(i  ,1,1,l) - Ftauylaxp(i  ,1,1,l))
       FtauylaxR(i,1,1,l) =     Ftauylaxm(i  ,1,1,l) -  FAC  *                &
                            mcl(Ftauylaxm(i  ,1,1,l) - Ftauylaxm(i  ,1,1,l),  &
                                Ftauylaxm(i  ,1,1,l) - Ftauylaxm(i  ,1,1,l))

       FtauzlaxL(i,1,1,l) =     Ftauzlaxp(i  ,1,1,l) +  FAC  *                &
                            mcl(Ftauzlaxp(i  ,1,1,l) - Ftauzlaxp(i  ,1,1,l),  &
                                Ftauzlaxp(i  ,1,1,l) - Ftauzlaxp(i  ,1,1,l))
       FtauzlaxR(i,1,1,l) =     Ftauzlaxm(i  ,1,1,l) -  FAC  *                &
                            mcl(Ftauzlaxm(i  ,1,1,l) - Ftauzlaxm(i  ,1,1,l),  &
                                Ftauzlaxm(i  ,1,1,l) - Ftauzlaxm(i  ,1,1,l))


       Ftauxlax(i,1,1,l)  =   0.5d0 * (FtauxlaxL(i,1,1,l) + FtauxlaxR(i,1,1,l))
       Ftauylax(i,1,1,l)  =   0.5d0 * (FtauylaxL(i,1,1,l) + FtauylaxR(i,1,1,l))
       Ftauzlax(i,1,1,l)  =   0.5d0 * (FtauzlaxL(i,1,1,l) + FtauzlaxR(i,1,1,l))

! -------------------------------------- EGLM--------------------------------------

       psilaxtauxL(i,1,1,l) =     psilaxtaup(i  ,1  ,1,l) +  FAC  *                   &
                              mcl(psilaxtaup(i+1,1  ,1,l) - psilaxtaup(i  ,1  ,1,l),  &
                                  psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i-1,1  ,1,l))
       psilaxtauxR(i,1,1,l) =     psilaxtaum(i+1,1  ,1,l) -  FAC  *                   &
                              mcl(psilaxtaup(i+2,1  ,1,l) - psilaxtaup(i+1,1  ,1,l),  &
                                  psilaxtaup(i+1,1  ,1,l) - psilaxtaup(i  ,1  ,1,l))

       psilaxtauyL(i,1,1,l) =     psilaxtaup(i  ,1  ,1,l) +  FAC  *                   &
                              mcl(psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i  ,1  ,1,l),  &
                                  psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i  ,1  ,1,l))
       psilaxtauyR(i,1,1,l) =     psilaxtaum(i  ,1  ,1,l) -  FAC  *                   &
                              mcl(psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i  ,1  ,1,l),  &
                                  psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i  ,1  ,1,l))

       psilaxtauzL(i,1,1,l) =     psilaxtaup(i  ,1  ,1,l) +  FAC  *                   &
                              mcl(psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i  ,1  ,1,l),  &
                                  psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i  ,1  ,1,l))
       psilaxtauzR(i,1,1,l) =     psilaxtaum(i  ,1  ,1,l) -  FAC  *                   &
                              mcl(psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i  ,1  ,1,l),  &
                                  psilaxtaup(i  ,1  ,1,l) - psilaxtaup(i  ,1  ,1,l))

       psitauxlax(i,1,1,l)  =   0.5d0 * (psilaxtauxL(i,1,1,l) + psilaxtauxR(i,1,1,l))
       psitauylax(i,1,1,l)  =   0.5d0 * (psilaxtauyL(i,1,1,l) + psilaxtauyR(i,1,1,l))
       psitauzlax(i,1,1,l)  =   0.5d0 * (psilaxtauzL(i,1,1,l) + psilaxtauzR(i,1,1,l))


       philaxtauxL(i,1,1,l) =     philaxtaup(i  ,1  ,1,l) +  FAC  *                   &
                              mcl(philaxtaup(i+1,1  ,1,l) - philaxtaup(i  ,1  ,1,l),  &
                                  philaxtaup(i  ,1  ,1,l) - philaxtaup(i-1,1  ,1,l))
       philaxtauxR(i,1,1,l) =     philaxtaum(i+1,1  ,1,l) -  FAC  *                   &
                              mcl(philaxtaup(i+2,1  ,1,l) - philaxtaup(i+1,1  ,1,l),  &
                                  philaxtaup(i+1,1  ,1,l) - philaxtaup(i  ,1  ,1,l))

       philaxtauyL(i,1,1,l) =     philaxtaup(i  ,1  ,1,l) +  FAC  *                   &
                              mcl(philaxtaup(i  ,1  ,1,l) - philaxtaup(i  ,1  ,1,l),  &
                                  philaxtaup(i  ,1  ,1,l) - philaxtaup(i  ,1  ,1,l))
       philaxtauyR(i,1,1,l) =     philaxtaum(i  ,1  ,1,l) -  FAC  *                   &
                              mcl(philaxtaup(i  ,1  ,1,l) - philaxtaup(i  ,1  ,1,l),  &
                                  philaxtaup(i  ,1  ,1,l) - philaxtaup(i  ,1  ,1,l))

       philaxtauzL(i,1,1,l) =     philaxtaup(i  ,1  ,1,l) +  FAC  *                   &
                              mcl(philaxtaup(i  ,1  ,1,l) - philaxtaup(i  ,1  ,1,l),  &
                                  philaxtaup(i  ,1  ,1,l) - philaxtaup(i  ,1  ,1,l))
       philaxtauzR(i,1,1,l) =     philaxtaum(i  ,1  ,1,l) -  FAC  *                   &
                              mcl(philaxtaup(i  ,1  ,1,l) - philaxtaup(i  ,1  ,1,l),  &
                                  philaxtaup(i  ,1  ,1,l) - philaxtaup(i  ,1  ,1,l))


       phitauxlax(i,1,1,l)  =   0.5d0 * (philaxtauxL(i,1,1,l) + philaxtauxR(i,1,1,l))
       phitauylax(i,1,1,l)  =   0.5d0 * (philaxtauyL(i,1,1,l) + philaxtauyR(i,1,1,l))
       phitauzlax(i,1,1,l)  =   0.5d0 * (philaxtauzL(i,1,1,l) + philaxtauzR(i,1,1,l))

! -------------------------------------- EGLM--------------------------------------



! Conserved Momentum flux tensor
!_________________________________________________________________________

       FSxxlaxL(i,1,1,l)  =     FSxxlaxp(i  ,1,1,l) + FAC  *                  &
                            mcl(FSxxlaxp(i+1,1,1,l) - FSxxlaxp(i  ,1,1,l),    &
                                FSxxlaxp(i  ,1,1,l) - FSxxlaxp(i-1,1,1,l))  
       FSxxlaxR(i,1,1,l)  =     FSxxlaxm(i+1,1,1,l) - FAC  *                  &
                            mcl(FSxxlaxm(i+2,1,1,l) - FSxxlaxm(i+1,1,1,l),    &
                                FSxxlaxm(i+1,1,1,l) - FSxxlaxm(i  ,1,1,l))
       FSxylaxL(i,1,1,l)  =     FSxylaxp(i  ,1,1,l) + FAC  *                  &
                            mcl(FSxylaxp(i+1,1,1,l) - FSxylaxp(i  ,1,1,l),    &
                                FSxylaxp(i  ,1,1,l) - FSxylaxp(i-1,1,1,l))  
       FSxylaxR(i,1,1,l)  =     FSxylaxm(i+1,1,1,l) - FAC  *                  &
                            mcl(FSxylaxm(i+2,1,1,l) - FSxylaxm(i+1,1,1,l),    &
                                FSxylaxm(i+1,1,1,l) - FSxylaxm(i  ,1,1,l))
       FSxzlaxL(i,1,1,l)  =     FSxzlaxp(i  ,1,1,l) + FAC  *                  &
                            mcl(FSxzlaxp(i+1,1,1,l) - FSxzlaxp(i  ,1,1,l),    &
                                FSxzlaxp(i  ,1,1,l) - FSxzlaxp(i-1,1,1,l))  
       FSxzlaxR(i,1,1,l)  =     FSxzlaxm(i+1,1,1,l) - FAC  *                  &
                            mcl(FSxzlaxm(i+2,1,1,l) - FSxzlaxm(i+1,1,1,l),    &
                                FSxzlaxm(i+1,1,1,l) - FSxzlaxm(i  ,1,1,l))


       FSxxlax(i,1,1,l)   =   0.5d0 * (FSxxlaxL(i,1,1,l)  + FSxxlaxR(i,1,1,l)  ) 
       FSxylax(i,1,1,l)   =   0.5d0 * (FSxylaxL(i,1,1,l)  + FSxylaxR(i,1,1,l)  )
       FSxzlax(i,1,1,l)   =   0.5d0 * (FSxzlaxL(i,1,1,l)  + FSxzlaxR(i,1,1,l)  )

       FSyxlaxL(i,1,1,l)  =     FSyxlaxp(i  ,1,1,l) +  FAC  *                 &
                            mcl(FSyxlaxp(i  ,1,1,l) - FSyxlaxp(i  ,1,1,l),    &
                                FSyxlaxp(i  ,1,1,l) - FSyxlaxp(i  ,1,1,l))
       FSyxlaxR(i,1,1,l)  =     FSyxlaxm(i  ,1,1,l) -  FAC  *                 &
                            mcl(FSyxlaxm(i  ,1,1,l) - FSyxlaxm(i  ,1,1,l),    &
                                FSyxlaxm(i  ,1,1,l) - FSyxlaxm(i  ,1,1,l))
       FSyylaxL(i,1,1,l)  =     FSyylaxp(i  ,1,1,l) +  FAC  *                 &
                            mcl(FSyylaxp(i  ,1,1,l) - FSyylaxp(i  ,1,1,l),    &
                                FSyylaxp(i  ,1,1,l) - FSyylaxp(i  ,1,1,l))
       FSyylaxR(i,1,1,l)  =     FSyylaxm(i  ,1,1,l) -  FAC  *                 &
                            mcl(FSyylaxm(i  ,1,1,l) - FSyylaxm(i  ,1,1,l),    &
                                FSyylaxm(i  ,1,1,l) - FSyylaxm(i  ,1,1,l))
       FSyzlaxL(i,1,1,l)  =     FSyzlaxp(i  ,1,1,l) +  FAC  *                 &
                            mcl(FSyzlaxp(i  ,1,1,l) - FSyzlaxp(i  ,1,1,l),    &
                                FSyzlaxp(i  ,1,1,l) - FSyzlaxp(i  ,1,1,l))
       FSyzlaxR(i,1,1,l)  =     FSyzlaxm(i  ,1,1,l) -  FAC  *                 &
                            mcl(FSyzlaxm(i  ,1,1,l) - FSyzlaxm(i  ,1,1,l),    &
                                FSyzlaxm(i  ,1,1,l) - FSyzlaxm(i  ,1,1,l))


 
       FSyxlax(i,1,1,l)   =   0.5d0 * (FSyxlaxL(i,1,1,l)  + FSyxlaxR(i,1,1,l)  )
       FSyylax(i,1,1,l)   =   0.5d0 * (FSyylaxL(i,1,1,l)  + FSyylaxR(i,1,1,l)  )
       FSyzlax(i,1,1,l)   =   0.5d0 * (FSyzlaxL(i,1,1,l)  + FSyzlaxR(i,1,1,l)  )


       FSzxlaxL(i,1,1,l)  =     FSzxlaxp(i  ,1,1,l) +  FAC  *                 &
                            mcl(FSzxlaxp(i  ,1,1,l) - FSzxlaxp(i  ,1,1,l),    &
                                FSzxlaxp(i  ,1,1,l) - FSzxlaxp(i  ,1,1,l))
       FSzxlaxR(i,1,1,l)  =     FSzxlaxm(i  ,1,1,l) -  FAC  *                 &
                            mcl(FSzxlaxm(i  ,1,1,l) - FSzxlaxm(i  ,1,1,l),    &
                                FSzxlaxm(i  ,1,1,l) - FSzxlaxm(i  ,1,1,l))
       FSzylaxL(i,1,1,l)  =     FSzylaxp(i  ,1,1,l) +  FAC  *                 &
                            mcl(FSzylaxp(i  ,1,1,l) - FSzylaxp(i  ,1,1,l),    &
                                FSzylaxp(i  ,1,1,l) - FSzylaxp(i  ,1,1,l))
       FSzylaxR(i,1,1,l)  =     FSzylaxm(i  ,1,1,l) -  FAC  *                 &
                            mcl(FSzylaxm(i  ,1,1,l) - FSzylaxm(i  ,1,1,l),    &
                                FSzylaxm(i  ,1,1,l) - FSzylaxm(i  ,1,1,l))
       FSzzlaxL(i,1,1,l)  =     FSzzlaxp(i  ,1,1,l) +  FAC  *                 &
                            mcl(FSzzlaxp(i  ,1,1,l) - FSzzlaxp(i  ,1,1,l),    &
                                FSzzlaxp(i  ,1,1,l) - FSzzlaxp(i  ,1,1,l))
       FSzzlaxR(i,1,1,l)  =     FSzzlaxm(i  ,1,1,l) -  FAC  *                 &
                            mcl(FSzzlaxm(i  ,1,1,l) - FSzzlaxm(i  ,1,1,l),    &
                                FSzzlaxm(i  ,1,1,l) - FSzzlaxm(i  ,1,1,l))

       FSzxlax(i,1,1,l)  =   0.5d0 * (FSzxlaxL(i,1,1,l)  + FSzxlaxR(i,1,1,l)  )
       FSzylax(i,1,1,l)  =   0.5d0 * (FSzylaxL(i,1,1,l)  + FSzylaxR(i,1,1,l)  )
       FSzzlax(i,1,1,l)  =   0.5d0 * (FSzzlaxL(i,1,1,l)  + FSzzlaxR(i,1,1,l)  )

! -------------------------------------- EGLM--------------------------------------

       BxSxlaxL(i,1,1,l)  =     BxlaxSxp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BxlaxSxp(i+1,1  ,1,l) - BxlaxSxp(i  ,1  ,1,l),    &
                                BxlaxSxp(i  ,1  ,1,l) - BxlaxSxp(i-1,1  ,1,l))  
       BxSxlaxR(i,1,1,l)  =     BxlaxSxm(i+1,1  ,1,l) - FAC  *                    &
                            mcl(BxlaxSxp(i+2,1  ,1,l) - BxlaxSxp(i+1,1  ,1,l),    &
                                BxlaxSxp(i+1,1  ,1,l) - BxlaxSxp(i  ,1  ,1,l))
       BySxlaxL(i,1,1,l)  =     BylaxSxp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BylaxSxp(i  ,1  ,1,l) - BylaxSxp(i  ,1  ,1,l),    &
                                BylaxSxp(i  ,1  ,1,l) - BylaxSxp(i  ,1  ,1,l))  
       BySxlaxR(i,1,1,l)  =     BylaxSxm(i  ,1  ,1,l) - FAC  *                    &
                            mcl(BylaxSxp(i  ,1  ,1,l) - BylaxSxp(i  ,1  ,1,l),    &
                                BylaxSxp(i  ,1  ,1,l) - BylaxSxp(i  ,1  ,1,l))
       BzSxlaxL(i,1,1,l)  =     BzlaxSxp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BzlaxSxp(i  ,1  ,1,l) - BzlaxSxp(i  ,1  ,1,l),    &
                                BzlaxSxp(i  ,1  ,1,l) - BzlaxSxp(i  ,1  ,1,l))  
       BzSxlaxR(i,1,1,l)  =     BzlaxSxm(i  ,1  ,1,l) - FAC  *                    &
                            mcl(BzlaxSxp(i  ,1  ,1,l) - BzlaxSxp(i  ,1  ,1,l),    &
                                BzlaxSxp(i  ,1  ,1,l) - BzlaxSxp(i  ,1  ,1,l))

       BxSxlax(i,1,1,l)   = 0.5d0 * (BxSxlaxL(i,1,1,l)  + BxSxlaxR(i,1,1,l)    )
       BySxlax(i,1,1,l)   = 0.5d0 * (BySxlaxL(i,1,1,l)  + BySxlaxR(i,1,1,l)    )
       BzSxlax(i,1,1,l)   = 0.5d0 * (BzSxlaxL(i,1,1,l)  + BzSxlaxR(i,1,1,l)    )


       BxSylaxL(i,1,1,l)  =     BxlaxSyp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BxlaxSyp(i+1,1  ,1,l) - BxlaxSyp(i  ,1  ,1,l),    &
                                BxlaxSyp(i  ,1  ,1,l) - BxlaxSyp(i-1,1  ,1,l))  
       BxSylaxR(i,1,1,l)  =     BxlaxSym(i+1,1  ,1,l) - FAC  *                    &
                            mcl(BxlaxSyp(i+2,1  ,1,l) - BxlaxSyp(i+1,1  ,1,l),    &
                                BxlaxSyp(i+1,1  ,1,l) - BxlaxSyp(i  ,1  ,1,l))
       BySylaxL(i,1,1,l)  =     BylaxSyp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BylaxSyp(i  ,1  ,1,l) - BylaxSyp(i  ,1  ,1,l),    &
                                BylaxSyp(i  ,1  ,1,l) - BylaxSyp(i  ,1  ,1,l))  
       BySylaxR(i,1,1,l)  =     BylaxSym(i  ,1  ,1,l) - FAC  *                    &
                            mcl(BylaxSyp(i  ,1  ,1,l) - BylaxSyp(i  ,1  ,1,l),    &
                                BylaxSyp(i  ,1  ,1,l) - BylaxSyp(i  ,1  ,1,l))
       BzSylaxL(i,1,1,l)  =     BzlaxSyp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BzlaxSyp(i  ,1  ,1,l) - BzlaxSyp(i  ,1  ,1,l),    &
                                BzlaxSyp(i  ,1  ,1,l) - BzlaxSyp(i  ,1  ,1,l))  
       BzSylaxR(i,1,1,l)  =     BzlaxSym(i  ,1  ,1,l) - FAC  *                    &
                            mcl(BzlaxSyp(i  ,1  ,1,l) - BzlaxSyp(i  ,1  ,1,l),    &
                                BzlaxSyp(i  ,1  ,1,l) - BzlaxSyp(i  ,1  ,1,l))

       BxSylax(i,1,1,l)   = 0.5d0 * (BxSylaxL(i,1,1,l)  + BxSylaxR(i,1,1,l)    )
       BySylax(i,1,1,l)   = 0.5d0 * (BySylaxL(i,1,1,l)  + BySylaxR(i,1,1,l)    )
       BzSylax(i,1,1,l)   = 0.5d0 * (BzSylaxL(i,1,1,l)  + BzSylaxR(i,1,1,l)    )


       BxSzlaxL(i,1,1,l)  =     BxlaxSzp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BxlaxSzp(i+1,1  ,1,l) - BxlaxSzp(i  ,1  ,1,l),    &
                                BxlaxSzp(i  ,1  ,1,l) - BxlaxSzp(i-1,1  ,1,l))  
       BxSzlaxR(i,1,1,l)  =     BxlaxSzm(i+1,1  ,1,l) - FAC  *                    &
                            mcl(BxlaxSzp(i+2,1  ,1,l) - BxlaxSzp(i+1,1  ,1,l),    &
                                BxlaxSzp(i+1,1  ,1,l) - BxlaxSzp(i  ,1  ,1,l))
       BySzlaxL(i,1,1,l)  =     BylaxSzp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BylaxSzp(i  ,1  ,1,l) - BylaxSzp(i  ,1  ,1,l),    &
                                BylaxSzp(i  ,1  ,1,l) - BylaxSzp(i  ,1  ,1,l))  
       BySzlaxR(i,1,1,l)  =     BylaxSzm(i  ,1  ,1,l) - FAC  *                    &
                            mcl(BylaxSzp(i  ,1  ,1,l) - BylaxSzp(i  ,1  ,1,l),    &
                                BylaxSzp(i  ,1  ,1,l) - BylaxSzp(i  ,1  ,1,l))
       BzSzlaxL(i,1,1,l)  =     BzlaxSzp(i  ,1  ,1,l) + FAC  *                    &
                            mcl(BzlaxSzp(i  ,1  ,1,l) - BzlaxSzp(i  ,1  ,1,l),    &
                                BzlaxSzp(i  ,1  ,1,l) - BzlaxSzp(i  ,1  ,1,l))  
       BzSzlaxR(i,1,1,l)  =     BzlaxSzm(i  ,1  ,1,l) - FAC  *                    &
                            mcl(BzlaxSzp(i  ,1  ,1,l) - BzlaxSzp(i  ,1  ,1,l),    &
                                BzlaxSzp(i  ,1  ,1,l) - BzlaxSzp(i  ,1  ,1,l))

       BxSzlax(i,1,1,l)   = 0.5d0 * (BxSzlaxL(i,1,1,l)  + BxSzlaxR(i,1,1,l)    )
       BySzlax(i,1,1,l)   = 0.5d0 * (BySzlaxL(i,1,1,l)  + BySzlaxR(i,1,1,l)    )
       BzSzlax(i,1,1,l)   = 0.5d0 * (BzSzlaxL(i,1,1,l)  + BzSzlaxR(i,1,1,l)    )

! -------------------------------------- EGLM--------------------------------------


   end do ! for i

!$OMP END DO

    else if (DIM == 2) then

!$OMP DO 

       do j=-2,jmax+2
          do i=-2,imax+2


          if (const_sp .eq. 1.d0) then
             c_sp = (gamma * (gamma -1.d0) * epsilon(i,j,1))/(1.d0 + gamma * epsilon(i,j,1))
             Sp   = const_sp * abs((sqrt(Vx(i,j,1)**2 + Vy(i,j,1)**2 + Vz(i,j,1)**2) + c_sp))
          else

             Sp = 1.d0

          end if

! Electric Gauss Law flux
!_________________________________________________________________________
       Exlaxpsip(i,j,1,l) =  Exint (i  ,j  ,1  ,l) + Sp * psiint(i  ,j  ,1  ,l)
       Exlaxpsim(i,j,1,l) =  Exint (i  ,j  ,1  ,l) - Sp * psiint(i  ,j  ,1  ,l)
       Eylaxpsip(i,j,1,l) =  Eyint (i  ,j  ,1  ,l) + Sp * psiint(i  ,j  ,1  ,l)
       Eylaxpsim(i,j,1,l) =  Eyint (i  ,j  ,1  ,l) - Sp * psiint(i  ,j  ,1  ,l)
       Ezlaxpsip(i,j,1,l) =  Ezint (i  ,j  ,1  ,l) + Sp * psiint(i  ,j  ,1  ,l)
       Ezlaxpsim(i,j,1,l) =  Ezint (i  ,j  ,1  ,l) - Sp * psiint(i  ,j  ,1  ,l)
! Magnetic Gauss Law flux
!_________________________________________________________________________

       Bxlaxphip(i,j,1,l) =  Bxint (i  ,j  ,1  ,l) + Sp * phiint(i  ,j  ,1  ,l) 
       Bxlaxphim(i,j,1,l) =  Bxint (i  ,j  ,1  ,l) - Sp * phiint(i  ,j  ,1  ,l) 
       Bylaxphip(i,j,1,l) =  Byint (i  ,j  ,1  ,l) + Sp * phiint(i  ,j  ,1  ,l) 
       Bylaxphim(i,j,1,l) =  Byint (i  ,j  ,1  ,l) - Sp * phiint(i  ,j  ,1  ,l) 
       Bzlaxphip(i,j,1,l) =  Bzint (i  ,j  ,1  ,l) + Sp * phiint(i  ,j  ,1  ,l) 
       Bzlaxphim(i,j,1,l) =  Bzint (i  ,j  ,1  ,l) - Sp * phiint(i  ,j  ,1  ,l) 

! Faraday Law flux
!_________________________________________________________________________

! Flujos Bx

       philaxBxp(i,j,1,l) =  phiint(i  ,j  ,1  ,l) + Sp * Bxint (i  ,j  ,1  ,l) 
       philaxBxm(i,j,1,l) =  phiint(i  ,j  ,1  ,l) - Sp * Bxint (i  ,j  ,1  ,l) 
       EzlaxBxp (i,j,1,l) =  Ezint (i  ,j  ,1  ,l) + Sp * Bxint (i  ,j  ,1  ,l) 
       EzlaxBxm (i,j,1,l) =  Ezint (i  ,j  ,1  ,l) - Sp * Bxint (i  ,j  ,1  ,l) 
       EylaxBxp (i,j,1,l) = -Eyint (i  ,j  ,1  ,l) + Sp * Bxint (i  ,j  ,1  ,l)  
       EylaxBxm (i,j,1,l) = -Eyint (i  ,j  ,1  ,l) - Sp * Bxint (i  ,j  ,1  ,l)
! Flujos By

       EzlaxByp (i,j,1,l) = -Ezint (i  ,j  ,1  ,l) + Sp * Byint (i  ,j  ,1  ,l) 
       EzlaxBym (i,j,1,l) = -Ezint (i  ,j  ,1  ,l) - Sp * Byint (i  ,j  ,1  ,l) 
       philaxByp(i,j,1,l) =  phiint(i  ,j  ,1  ,l) + Sp * Byint (i  ,j  ,1  ,l) 
       philaxBym(i,j,1,l) =  phiint(i  ,j  ,1  ,l) - Sp * Byint (i  ,j  ,1  ,l) 
       ExlaxByp (i,j,1,l) =  Exint (i  ,j  ,1  ,l) + Sp * Byint (i  ,j  ,1  ,l) 
       ExlaxBym (i,j,1,l) =  Exint (i  ,j  ,1  ,l) - Sp * Byint (i  ,j  ,1  ,l) 


! Flujos Bz

       EylaxBzp (i,j,1,l) =  Eyint (i  ,j  ,1  ,l) + Sp * Bzint (i  ,j  ,1  ,l)   
       EylaxBzm (i,j,1,l) =  Eyint (i  ,j  ,1  ,l) - Sp * Bzint (i  ,j  ,1  ,l)   
       ExlaxBzp (i,j,1,l) = -Exint (i  ,j  ,1  ,l) + Sp * Bzint (i  ,j  ,1  ,l) 
       ExlaxBzm (i,j,1,l) = -Exint (i  ,j  ,1  ,l) - Sp * Bzint (i  ,j  ,1  ,l) 
       philaxBzp(i,j,1,l) =  phiint(i  ,j  ,1  ,l) + Sp * Bzint (i  ,j  ,1  ,l)  
       philaxBzm(i,j,1,l) =  phiint(i  ,j  ,1  ,l) - Sp * Bzint (i  ,j  ,1  ,l)  


! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxExp(i,j,1,l) =  psiint(i  ,j  ,1  ,l) + Sp * Exint (i  ,j  ,1  ,l) 
       psilaxExm(i,j,1,l) =  psiint(i  ,j  ,1  ,l) - Sp * Exint (i  ,j  ,1  ,l) 
       BzlaxExp (i,j,1,l) = -Bzint (i  ,j  ,1  ,l) + Sp * Exint (i  ,j  ,1  ,l) 
       BzlaxExm (i,j,1,l) = -Bzint (i  ,j  ,1  ,l) - Sp * Exint (i  ,j  ,1  ,l) 
       BylaxExp (i,j,1,l) =  Byint (i  ,j  ,1  ,l) + Sp * Exint (i  ,j  ,1  ,l) 
       BylaxExm (i,j,1,l) =  Byint (i  ,j  ,1  ,l) - Sp * Exint (i  ,j  ,1  ,l) 

! Flujos Ey

       BzlaxEyp (i,j,1,l) =  Bzint (i  ,j  ,1  ,l) + Sp * Eyint (i  ,j  ,1  ,l) 
       BzlaxEym (i,j,1,l) =  Bzint (i  ,j  ,1  ,l) - Sp * Eyint (i  ,j  ,1  ,l) 
       psilaxEyp(i,j,1,l) =  psiint(i  ,j  ,1  ,l) + Sp * Eyint (i  ,j  ,1  ,l) 
       psilaxEym(i,j,1,l) =  psiint(i  ,j  ,1  ,l) - Sp * Eyint (i  ,j  ,1  ,l) 
       BxlaxEyp (i,j,1,l) = -Bxint (i  ,j  ,1  ,l) + Sp * Eyint (i  ,j  ,1  ,l) 
       BxlaxEym (i,j,1,l) = -Bxint (i  ,j  ,1  ,l) - Sp * Eyint (i  ,j  ,1  ,l) 

! Flujos Ez

       BylaxEzp (i,j,1,l) = -Byint (i  ,j  ,1  ,l) + Sp * Ezint (i  ,j  ,1  ,l) 
       BylaxEzm (i,j,1,l) = -Byint (i  ,j  ,1  ,l) - Sp * Ezint (i  ,j  ,1  ,l)  
       BxlaxEzp (i,j,1,l) =  Bxint (i  ,j  ,1  ,l) + Sp * Ezint (i  ,j  ,1  ,l) 
       BxlaxEzm (i,j,1,l) =  Bxint (i  ,j  ,1  ,l) - Sp * Ezint (i  ,j  ,1  ,l) 
       psilaxEzp(i,j,1,l) =  psiint(i  ,j  ,1  ,l) + Sp * Ezint (i  ,j  ,1  ,l) 
       psilaxEzm(i,j,1,l) =  psiint(i  ,j  ,1  ,l) - Sp * Ezint (i  ,j  ,1  ,l) 



! Conserved current flux
!_________________________________________________________________________

       Jxlaxp(i,j,1,l)    =  Jxint(i  ,j  ,1  ,l)  + Sp * qint (i  ,j  ,1  ,l) 
       Jxlaxm(i,j,1,l)    =  Jxint(i  ,j  ,1  ,l)  - Sp * qint (i  ,j  ,1  ,l) 
       Jylaxp(i,j,1,l)    =  Jyint(i  ,j  ,1  ,l)  + Sp * qint (i  ,j  ,1  ,l) 
       Jylaxm(i,j,1,l)    =  Jyint(i  ,j  ,1  ,l)  - Sp * qint (i  ,j  ,1  ,l) 
       Jzlaxp(i,j,1,l)    =  Jzint(i  ,j  ,1  ,l)  + Sp * qint (i  ,j  ,1  ,l) 
       Jzlaxm(i,j,1,l)    =  Jzint(i  ,j  ,1  ,l)  - Sp * qint (i  ,j  ,1  ,l) 

! Conserved Mass flux
!_________________________________________________________________________

       FDxlaxp(i,j,1,l)   =  FDxint(i  ,j  ,1  ,l) + Sp * DDint (i  ,j  ,1  ,l)  
       FDxlaxm(i,j,1,l)   =  FDxint(i  ,j  ,1  ,l) - Sp * DDint (i  ,j  ,1  ,l)  
       FDylaxp(i,j,1,l)   =  FDyint(i  ,j  ,1  ,l) + Sp * DDint (i  ,j  ,1  ,l) 
       FDylaxm(i,j,1,l)   =  FDyint(i  ,j  ,1  ,l) - Sp * DDint (i  ,j  ,1  ,l) 
       FDzlaxp(i,j,1,l)   =  FDzint(i  ,j  ,1  ,l) + Sp * DDint (i  ,j  ,1  ,l) 
       FDzlaxm(i,j,1,l)   =  FDzint(i  ,j  ,1  ,l) - Sp * DDint (i  ,j  ,1  ,l) 

! Conserved Energy flux
!_________________________________________________________________________

       Ftauxlaxp(i,j,1,l) =  Ftauxint(i  ,j  ,1  ,l) + Sp * tauint  (i  ,j  ,1  ,l) 
       Ftauxlaxm(i,j,1,l) =  Ftauxint(i  ,j  ,1  ,l) - Sp * tauint  (i  ,j  ,1  ,l) 
       Ftauylaxp(i,j,1,l) =  Ftauyint(i  ,j  ,1  ,l) + Sp * tauint  (i  ,j  ,1  ,l) 
       Ftauylaxm(i,j,1,l) =  Ftauyint(i  ,j  ,1  ,l) - Sp * tauint  (i  ,j  ,1  ,l) 
       Ftauzlaxp(i,j,1,l) =  Ftauzint(i  ,j  ,1  ,l) + Sp * tauint  (i  ,j  ,1  ,l) 
       Ftauzlaxm(i,j,1,l) =  Ftauzint(i  ,j  ,1  ,l) - Sp * tauint  (i  ,j  ,1  ,l) 
! -------------------------------------- EGLM--------------------------------------
       psilaxtaup(i,j,1,l) = psiint  (i  ,j  ,1  ,l) + Sp * tauint  (i  ,j  ,1  ,l) 
       psilaxtaum(i,j,1,l) = psiint  (i  ,j  ,1  ,l) - Sp * tauint  (i  ,j  ,1  ,l) 
       philaxtaup(i,j,1,l) = phiint  (i  ,j  ,1  ,l) + Sp * tauint  (i  ,j  ,1  ,l)
       philaxtaum(i,j,1,l) = phiint  (i  ,j  ,1  ,l) - Sp * tauint  (i  ,j  ,1  ,l)  
! -------------------------------------- EGLM--------------------------------------

! Conserved Momentum flux tensor
!_________________________________________________________________________

       FSxxlaxp(i,j,1,l)  =  FSxxint(i  ,j  ,1  ,l) + Sp * Sxint (i  ,j  ,1  ,l) 
       FSxxlaxm(i,j,1,l)  =  FSxxint(i  ,j  ,1  ,l) - Sp * Sxint (i  ,j  ,1  ,l) 
       FSxylaxp(i,j,1,l)  =  FSxyint(i  ,j  ,1  ,l) + Sp * Syint (i  ,j  ,1  ,l)
       FSxylaxm(i,j,1,l)  =  FSxyint(i  ,j  ,1  ,l) - Sp * Syint (i  ,j  ,1  ,l)
       FSxzlaxp(i,j,1,l)  =  FSxzint(i  ,j  ,1  ,l) + Sp * Szint (i  ,j  ,1  ,l) 
       FSxzlaxm(i,j,1,l)  =  FSxzint(i  ,j  ,1  ,l) - Sp * Szint (i  ,j  ,1  ,l) 

       FSyxlaxp(i,j,1,l)  =  FSyxint(i  ,j  ,1  ,l) + Sp * Sxint (i  ,j  ,1  ,l) 
       FSyxlaxm(i,j,1,l)  =  FSyxint(i  ,j  ,1  ,l) - Sp * Sxint (i  ,j  ,1  ,l) 
       FSyylaxp(i,j,1,l)  =  FSyyint(i  ,j  ,1  ,l) + Sp * Syint (i  ,j  ,1  ,l) 
       FSyylaxm(i,j,1,l)  =  FSyyint(i  ,j  ,1  ,l) - Sp * Syint (i  ,j  ,1  ,l) 
       FSyzlaxp(i,j,1,l)  =  FSyzint(i  ,j  ,1  ,l) + Sp * Szint (i  ,j  ,1  ,l) 
       FSyzlaxm(i,j,1,l)  =  FSyzint(i  ,j  ,1  ,l) - Sp * Szint (i  ,j  ,1  ,l) 

       FSzxlaxp(i,j,1,l)  =  FSzxint(i  ,j  ,1  ,l) + Sp * Sxint (i  ,j  ,1  ,l) 
       FSzxlaxm(i,j,1,l)  =  FSzxint(i  ,j  ,1  ,l) - Sp * Sxint (i  ,j  ,1  ,l) 
       FSzylaxp(i,j,1,l)  =  FSzyint(i  ,j  ,1  ,l) + Sp * Syint (i  ,j  ,1  ,l) 
       FSzylaxm(i,j,1,l)  =  FSzyint(i  ,j  ,1  ,l) - Sp * Syint (i  ,j  ,1  ,l) 
       FSzzlaxp(i,j,1,l)  =  FSzzint(i  ,j  ,1  ,l) + Sp * Szint (i  ,j  ,1  ,l)   
       FSzzlaxm(i,j,1,l)  =  FSzzint(i  ,j  ,1  ,l) - Sp * Szint (i  ,j  ,1  ,l)   

! -------------------------------------- EGLM--------------------------------------

       BxlaxSxp(i,j,1,l)   =  Bxint  (i  ,j  ,1  ,l) + Sp * Sxint (i  ,j  ,1  ,l)
       BxlaxSxm(i,j,1,l)   =  Bxint  (i  ,j  ,1  ,l) - Sp * Sxint (i  ,j  ,1  ,l)
       BylaxSxp(i,j,1,l)   =  Byint  (i  ,j  ,1  ,l) + Sp * Sxint (i  ,j  ,1  ,l)
       BylaxSxm(i,j,1,l)   =  Byint  (i  ,j  ,1  ,l) - Sp * Sxint (i  ,j  ,1  ,l)
       BzlaxSxp(i,j,1,l)   =  Bzint  (i  ,j  ,1  ,l) + Sp * Sxint (i  ,j  ,1  ,l)
       BzlaxSxm(i,j,1,l)   =  Bzint  (i  ,j  ,1  ,l) - Sp * Sxint (i  ,j  ,1  ,l)

       BxlaxSyp(i,j,1,l)   =  Bxint  (i  ,j  ,1  ,l) + Sp * Syint (i  ,j  ,1  ,l)
       BxlaxSym(i,j,1,l)   =  Bxint  (i  ,j  ,1  ,l) - Sp * Syint (i  ,j  ,1  ,l)
       BylaxSyp(i,j,1,l)   =  Byint  (i  ,j  ,1  ,l) + Sp * Syint (i  ,j  ,1  ,l)
       BylaxSym(i,j,1,l)   =  Byint  (i  ,j  ,1  ,l) - Sp * Syint (i  ,j  ,1  ,l)
       BzlaxSyp(i,j,1,l)   =  Bzint  (i  ,j  ,1  ,l) + Sp * Syint (i  ,j  ,1  ,l)
       BzlaxSym(i,j,1,l)   =  Bzint  (i  ,j  ,1  ,l) - Sp * Syint (i  ,j  ,1  ,l)

       BxlaxSzp(i,j,1,l)   =  Bxint  (i  ,j  ,1  ,l) + Sp * Szint (i  ,j  ,1  ,l)
       BxlaxSzm(i,j,1,l)   =  Bxint  (i  ,j  ,1  ,l) - Sp * Szint (i  ,j  ,1  ,l)
       BylaxSzp(i,j,1,l)   =  Byint  (i  ,j  ,1  ,l) + Sp * Szint (i  ,j  ,1  ,l)
       BylaxSzm(i,j,1,l)   =  Byint  (i  ,j  ,1  ,l) - Sp * Szint (i  ,j  ,1  ,l)
       BzlaxSzp(i,j,1,l)   =  Bzint  (i  ,j  ,1  ,l) + Sp * Szint (i  ,j  ,1  ,l)
       BzlaxSzm(i,j,1,l)   =  Bzint  (i  ,j  ,1  ,l) - Sp * Szint (i  ,j  ,1  ,l)



! -------------------------------------- EGLM--------------------------------------


     end do ! for i
   end do ! for j

!$OMP END DO

!$OMP DO 

   do j=-1,jmax     
      do i=-1,imax


! Electric Gauss Law flux
!_________________________________________________________________________

       ExlaxpsiL(i,j,1,l) =     Exlaxpsip(i  ,j  ,1,l) + FAC  *                    &
                            mcl(Exlaxpsip(i+1,j  ,1,l) - Exlaxpsip(i  ,j  ,1,l) ,  &
                                Exlaxpsip(i  ,j  ,1,l) - Exlaxpsip(i-1,j  ,1,l))
       ExlaxpsiR(i,j,1,l) =     Exlaxpsim(i+1,j  ,1,l) - FAC  *                    &
                            mcl(Exlaxpsip(i+2,j  ,1,l) - Exlaxpsip(i+1,j  ,1,l) ,  &
                                Exlaxpsip(i+1,j  ,1,l) - Exlaxpsip(i  ,j  ,1,l))

       EylaxpsiL(i,j,1,l) =     Eylaxpsip(i  ,j  ,1,l) + FAC  *                    &
                            mcl(Eylaxpsip(i  ,j+1,1,l) - Eylaxpsip(i  ,j  ,1,l) ,  &
                                Eylaxpsip(i  ,j  ,1,l) - Eylaxpsip(i  ,j-1,1,l))
       EylaxpsiR(i,j,1,l) =     Eylaxpsim(i  ,j+1,1,l) - FAC  *                    &
                            mcl(Eylaxpsip(i  ,j+2,1,l) - Eylaxpsip(i  ,j+1,1,l) ,  &
                                Eylaxpsip(i  ,j+1,1,l) - Eylaxpsip(i  ,j  ,1,l))

       EzlaxpsiL(i,j,1,l) =     Ezlaxpsip(i  ,j  ,1,l) + FAC  *                    &
                            mcl(Ezlaxpsip(i  ,j  ,1,l) - Ezlaxpsip(i  ,j  ,1,l) ,  &
                                Ezlaxpsip(i  ,j  ,1,l) - Ezlaxpsip(i  ,j  ,1,l))
       EzlaxpsiR(i,j,1,l) =     Ezlaxpsim(i  ,j  ,1,l) - FAC  *                    &
                            mcl(Ezlaxpsip(i  ,j  ,1,l) - Ezlaxpsip(i  ,j  ,1,l) ,  &
                                Ezlaxpsip(i  ,j  ,1,l) - Ezlaxpsip(i  ,j  ,1,l))

       Exlaxpsi(i,j,1,l)  =   0.5d0 * (ExlaxpsiL(i,j,1,l) + ExlaxpsiR(i,j,1,l))  
       Eylaxpsi(i,j,1,l)  =   0.5d0 * (EylaxpsiL(i,j,1,l) + EylaxpsiR(i,j,1,l))
       Ezlaxpsi(i,j,1,l)  =   0.5d0 * (EzlaxpsiL(i,j,1,l) + EzlaxpsiR(i,j,1,l))

! Magnetic Gauss Law flux
!_________________________________________________________________________

       BxlaxphiL(i,j,1,l) =     Bxlaxphip(i  ,j  ,1,l) + FAC  *                    &
                            mcl(Bxlaxphip(i+1,j  ,1,l) - Bxlaxphip(i  ,j  ,1,l) ,  &
                                Bxlaxphip(i  ,j  ,1,l) - Bxlaxphip(i-1,j  ,1,l))
       BxlaxphiR(i,j,1,l) =     Bxlaxphim(i+1,j  ,1,l) - FAC  *                    &
                            mcl(Bxlaxphip(i+2,j  ,1,l) - Bxlaxphip(i+1,j  ,1,l) ,  &
                                Bxlaxphip(i+1,j  ,1,l) - Bxlaxphip(i  ,j  ,1,l))

       BylaxphiL(i,j,1,l) =     Bylaxphip(i  ,j  ,1,l) + FAC  *                    &
                            mcl(Bylaxphip(i  ,j+1,1,l) - Bylaxphip(i  ,j  ,1,l) ,  &
                                Bylaxphip(i  ,j  ,1,l) - Bylaxphip(i  ,j-1,1,l))
       BylaxphiR(i,j,1,l) =     Bylaxphim(i  ,j+1,1,l) - FAC  *                    &
                            mcl(Bylaxphip(i  ,j+2,1,l) - Bylaxphip(i  ,j+1,1,l) ,  &
                                Bylaxphip(i  ,j+1,1,l) - Bylaxphip(i  ,j  ,1,l))

       BzlaxphiL(i,j,1,l) =     Bzlaxphip(i  ,j  ,1,l) + FAC  *                    &
                            mcl(Bzlaxphip(i  ,j  ,1,l) - Bzlaxphip(i  ,j  ,1,l) ,  &
                                Bzlaxphip(i  ,j  ,1,l) - Bzlaxphip(i  ,j  ,1,l))
       BzlaxphiR(i,j,1,l) =     Bzlaxphim(i  ,j  ,1,l) - FAC  *                    &
                            mcl(Bzlaxphip(i  ,j  ,1,l) - Bzlaxphip(i  ,j  ,1,l) ,  &
                                Bzlaxphip(i  ,j  ,1,l) - Bzlaxphip(i  ,j  ,1,l))

       Bxlaxphi(i,j,1,l)  =   0.5d0 * (BxlaxphiL(i,j,1,l) + BxlaxphiR(i,j,1,l))
       Bylaxphi(i,j,1,l)  =   0.5d0 * (BylaxphiL(i,j,1,l) + BylaxphiR(i,j,1,l))
       Bzlaxphi(i,j,1,l)  =   0.5d0 * (BzlaxphiL(i,j,1,l) + BzlaxphiR(i,j,1,l))

! Faraday Law flux
!_________________________________________________________________________

! Flujos Bx

       philaxBxL(i,j,1,l) =     philaxBxp(i  ,j  ,1,l) + FAC  *                   &
                            mcl(philaxBxp(i+1,j  ,1,l) - philaxBxp(i  ,j  ,1,l),  &
                                philaxBxp(i  ,j  ,1,l) - philaxBxp(i-1,j  ,1,l))
       philaxBxR(i,j,1,l) =     philaxBxm(i+1,j  ,1,l) - FAC  *                   &
                            mcl(philaxBxp(i+2,j  ,1,l) - philaxBxp(i+1,j  ,1,l),  &
                                philaxBxp(i+1,j  ,1,l) - philaxBxp(i  ,j  ,1,l)) 

       EzlaxBxL(i,j,1,l)  =     EzlaxBxp(i  ,j  ,1,l)  + FAC  *                   &
                            mcl(EzlaxBxp(i  ,j+1,1,l) - EzlaxBxp(i  ,j  ,1,l),    &
                                EzlaxBxp(i  ,j  ,1,l) - EzlaxBxp(i  ,j-1,1,l))
       EzlaxBxR(i,j,1,l)  =     EzlaxBxm(i  ,j+1,1,l) - FAC  *                    &
                            mcl(EzlaxBxp(i  ,j+2,1,l) - EzlaxBxp(i  ,j+1,1,l),    &
                                EzlaxBxp(i  ,j+1,1,l) - EzlaxBxp(i  ,j  ,1,l))  

       EylaxBxL(i,j,1,l)  =     EylaxBxp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(EylaxBxp(i  ,j  ,1,l) - EylaxBxp(i  ,j  ,1,l),    &
                                EylaxBxp(i  ,j  ,1,l) - EylaxBxp(i  ,j  ,1,l))
       EylaxBxR(i,j,1,l)  =     EylaxBxm(i  ,j  ,1,l) - FAC  *                    &
                            mcl(EylaxBxp(i  ,j  ,1,l) - EylaxBxp(i  ,j  ,1,l),    &
                                EylaxBxp(i  ,j  ,1,l) - EylaxBxp(i  ,j  ,1,l))   


       philaxBx(i,j,1,l)  =   0.5d0 * (philaxBxL(i,j,1,l) + philaxBxR(i,j,1,l))
       EzlaxBx (i,j,1,l)  =   0.5d0 * (EzlaxBxL (i,j,1,l) + EzlaxBxR (i,j,1,l))
       EylaxBx (i,j,1,l)  = - 0.5d0 * (EylaxBxL (i,j,1,l) + EylaxBxR (i,j,1,l))   



! Flujos By

       EzlaxByL(i,j,1,l)  =     EzlaxByp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(EzlaxByp(i+1,j  ,1,l) - EzlaxByp(i  ,j  ,1,l),    &
                                EzlaxByp(i  ,j  ,1,l) - EzlaxByp(i-1,j  ,1,l))
       EzlaxByR(i,j,1,l)  =     EzlaxBym(i+1,j  ,1,l) - FAC  *                    &
                            mcl(EzlaxByp(i+2,j  ,1,l) - EzlaxByp(i+1,j  ,1,l),    &
                                EzlaxByp(i+1,j  ,1,l) - EzlaxByp(i  ,j  ,1,l))

       philaxByL(i,j,1,l) =     philaxByp(i  ,j  ,1,l) + FAC  *                   &
                            mcl(philaxByp(i  ,j+1,1,l) - philaxByp(i  ,j  ,1,l),  &
                                philaxByp(i  ,j  ,1,l) - philaxByp(i  ,j-1,1,l))
       philaxByR(i,j,1,l) =     philaxBym(i  ,j+1,1,l) - FAC  *                   &
                            mcl(philaxByp(i  ,j+2,1,l) - philaxByp(i  ,j+1,1,l),  &
                                philaxByp(i  ,j+1,1,l) - philaxByp(i  ,j  ,1,l))

       ExlaxByL(i,j,1,l)  =     ExlaxByp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(ExlaxByp(i  ,j  ,1,l) - ExlaxByp(i  ,j  ,1,l),    &
                                ExlaxByp(i  ,j  ,1,l) - ExlaxByp(i  ,j  ,1,l))
       ExlaxByR(i,j,1,l)  =     ExlaxBym(i  ,j  ,1,l) - FAC  *                    &
                            mcl(ExlaxByp(i  ,j  ,1,l) - ExlaxByp(i  ,j  ,1,l),    &
                                ExlaxByp(i  ,j  ,1,l) - ExlaxByp(i  ,j  ,1,l))

       EzlaxBy (i,j,1,l)  = - 0.5d0 * (EzlaxByL (i,j,1,l) + EzlaxByR (i,j,1,l))
       philaxBy(i,j,1,l)  =   0.5d0 * (philaxByL(i,j,1,l) + philaxByR(i,j,1,l)) 
       ExlaxBy (i,j,1,l)  =   0.5d0 * (ExlaxByL (i,j,1,l) + ExlaxByR (i,j,1,l))  

! Flujos Bz

       EylaxBzL(i,j,1,l)  =     EylaxBzp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(EylaxBzp(i+1,j  ,1,l) - EylaxBzp(i  ,j  ,1,l),    &
                                EylaxBzp(i  ,j  ,1,l) - EylaxBzp(i-1,j  ,1,l))  
       EylaxBzR(i,j,1,l)  =     EylaxBzm(i+1,j  ,1,l) - FAC  *                    &
                            mcl(EylaxBzp(i+2,j  ,1,l) - EylaxBzp(i+1,j  ,1,l),    &
                                EylaxBzp(i+1,j  ,1,l) - EylaxBzp(i  ,j  ,1,l))

       ExlaxBzL(i,j,1,l)  =     ExlaxBzp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(ExlaxBzp(i  ,j+1,1,l) - ExlaxBzp(i  ,j  ,1,l),    &
                                ExlaxBzp(i  ,j  ,1,l) - ExlaxBzp(i  ,j-1,1,l))
       ExlaxBzR(i,j,1,l)  =     ExlaxBzm(i  ,j+1,1,l) - FAC  *                    &
                            mcl(ExlaxBzp(i  ,j+2,1,l) - ExlaxBzp(i  ,j+1,1,l),    &
                                ExlaxBzp(i  ,j+1,1,l) - ExlaxBzp(i  ,j  ,1,l))

       philaxBzL(i,j,1,l) =     philaxBzp(i  ,j  ,1,l) + FAC  *                   &
                            mcl(philaxBzp(i  ,j  ,1,l) - philaxBzp(i  ,j  ,1,l),  &
                                philaxBzp(i  ,j  ,1,l) - philaxBzp(i  ,j  ,1,l))
       philaxBzR(i,j,1,l) =     philaxBzm(i  ,j  ,1,l) - FAC  *                   &
                            mcl(philaxBzp(i  ,j  ,1,l) - philaxBzp(i  ,j  ,1,l),  &
                                philaxBzp(i  ,j  ,1,l) - philaxBzp(i  ,j  ,1,l))

       EylaxBz (i,j,1,l)  =   0.5d0 * (EylaxBzL (i,j,1,l)  + EylaxBzR (i,j,1,l))
       ExlaxBz (i,j,1,l)  = - 0.5d0 * (ExlaxBzL (i,j,1,l)  + ExlaxBzR (i,j,1,l))
       philaxBz(i,j,1,l)  =   0.5d0 * (philaxBzL(i,j,1,l)  + philaxBzR(i,j,1,l))


! Ampere Law flux
!_________________________________________________________________________

! Flujos Ex

       psilaxExL(i,j,1,l) =     psilaxExp(i  ,j  ,1,l) + FAC  *                   &
                            mcl(psilaxExp(i+1,j  ,1,l) - psilaxExp(i  ,j  ,1,l),  &
                                psilaxExp(i  ,j  ,1,l) - psilaxExp(i-1,j  ,1,l))  
       psilaxExR(i,j,1,l) =     psilaxExm(i+1,j  ,1,l) - FAC  *                   &
                            mcl(psilaxExp(i+2,j  ,1,l) - psilaxExp(i+1,j  ,1,l),  &
                                psilaxExp(i+1,j  ,1,l) - psilaxExp(i  ,j  ,1,l))

       BzlaxExL(i,j,1,l)  =     BzlaxExp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BzlaxExp(i  ,j+1,1,l) - BzlaxExp(i  ,j  ,1,l),    &
                                BzlaxExp(i  ,j  ,1,l) - BzlaxExp(i  ,j-1,1,l))
       BzlaxExR(i,j,1,l)  =     BzlaxExm(i  ,j+1,1,l) - FAC  *                    &
                            mcl(BzlaxExp(i  ,j+2,1,l) - BzlaxExp(i  ,j+1,1,l),    &
                                BzlaxExp(i  ,j+1,1,l) - BzlaxExp(i  ,j  ,1,l))

       BylaxExL(i,j,1,l)  =     BylaxExp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BylaxExp(i  ,j  ,1,l) - BylaxExp(i  ,j  ,1,l),    &
                                BylaxExp(i  ,j  ,1,l) - BylaxExp(i  ,j  ,1,l))
       BylaxExR(i,j,1,l)  =     BylaxExm(i  ,j  ,1,l) - FAC  *                    &
                            mcl(BylaxExp(i  ,j  ,1,l) - BylaxExp(i  ,j  ,1,l),    &
                                BylaxExp(i  ,j  ,1,l) - BylaxExp(i  ,j  ,1,l))

       psilaxEx(i,j,1,l)  =   0.5d0 * (psilaxExL(i,j,1,l) + psilaxExR(i,j,1,l))
       BzlaxEx(i,j,1,l)   = - 0.5d0 * (BzlaxExL (i,j,1,l) + BzlaxExR (i,j,1,l))  
       BylaxEx(i,j,1,l)   =   0.5d0 * (BylaxExL (i,j,1,l) + BylaxExR (i,j,1,l))  

! Flujos Ey

       BzlaxEyL(i,j,1,l)  =     BzlaxEyp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BzlaxEyp(i+1,j  ,1,l) - BzlaxEyp(i  ,j  ,1,l),    &
                                BzlaxEyp(i  ,j  ,1,l) - BzlaxEyp(i-1,j  ,1,l))  
       BzlaxEyR(i,j,1,l)  =     BzlaxEym(i+1,j  ,1,l) - FAC  *                    &
                            mcl(BzlaxEyp(i+2,j  ,1,l) - BzlaxEyp(i+1,j  ,1,l),    &
                                BzlaxEyp(i+1,j  ,1,l) - BzlaxEyp(i  ,j  ,1,l))

       psilaxEyL(i,j,1,l) =     psilaxEyp(i  ,j  ,1,l) + FAC  *                   &
                            mcl(psilaxEyp(i  ,j+1,1,l) - psilaxEyp(i  ,j  ,1,l),  &
                                psilaxEyp(i  ,j  ,1,l) - psilaxEyp(i  ,j-1,1,l))
       psilaxEyR(i,j,1,l) =     psilaxEym(i  ,j+1,1,l) - FAC  *                   &
                            mcl(psilaxEyp(i  ,j+2,1,l) - psilaxEyp(i  ,j+1,1,l),  &
                                psilaxEyp(i  ,j+1,1,l) - psilaxEyp(i  ,j  ,1,l))

       BxlaxEyL(i,j,1,l)  =     BxlaxEyp(i  ,j  ,1,l) + FAC  *                   &
                            mcl(BxlaxEyp(i  ,j  ,1,l) - BxlaxEyp(i  ,j  ,1,l),   &
                                BxlaxEyp(i  ,j  ,1,l) - BxlaxEyp(i  ,j  ,1,l))
       BxlaxEyR(i,j,1,l)  =     BxlaxEym(i  ,j  ,1,l) - FAC  *                   &
                            mcl(BxlaxEyp(i  ,j  ,1,l) - BxlaxEyp(i  ,j  ,1,l),   &
                                BxlaxEyp(i  ,j  ,1,l) - BxlaxEyp(i  ,j  ,1,l))

       BzlaxEy (i,j,1,l)  =   0.5d0 * (BzlaxEyL (i,j,1,l) + BzlaxEyR (i,j,1,l))
       psilaxEy(i,j,1,l)  =   0.5d0 * (psilaxEyL(i,j,1,l) + psilaxEyR(i,j,1,l))
       BxlaxEy (i,j,1,l)  = - 0.5d0 * (BxlaxEyL (i,j,1,l) + BxlaxEyR (i,j,1,l))

! Flujos Ez

       BylaxEzL(i,j,1,l)  =     BylaxEzp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BylaxEzp(i+1,j  ,1,l) - BylaxEzp(i  ,j  ,1,l),    &
                                BylaxEzp(i  ,j  ,1,l) - BylaxEzp(i-1,j  ,1,l)) 
       BylaxEzR(i,j,1,l)  =     BylaxEzm(i+1,j  ,1,l) -  FAC  *                   &
                            mcl(BylaxEzp(i+2,j  ,1,l) - BylaxEzp(i+1,j  ,1,l),    &
                                BylaxEzp(i+1,j  ,1,l) - BylaxEzp(i  ,j  ,1,l))

       BxlaxEzL(i,j,1,l)  =     BxlaxEzp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BxlaxEzp(i  ,j+1,1,l) - BxlaxEzp(i  ,j  ,1,l),    &
                                BxlaxEzp(i  ,j  ,1,l) - BxlaxEzp(i  ,j-1,1,l))
       BxlaxEzR(i,j,1,l)  =     BxlaxEzm(i  ,j+1,1,l) - FAC  *                    &
                            mcl(BxlaxEzp(i  ,j+2,1,l) - BxlaxEzp(i  ,j+1,1,l),    &
                                BxlaxEzp(i  ,j+1,1,l) - BxlaxEzp(i  ,j  ,1,l))

       psilaxEzL(i,j,1,l) =     psilaxEzp(i  ,j  ,1,l) + FAC  *                   &
                            mcl(psilaxEzp(i  ,j  ,1,l) - psilaxEzp(i  ,j  ,1,l),  &
                                psilaxEzp(i  ,j  ,1,l) - psilaxEzp(i  ,j  ,1,l))
       psilaxEzR(i,j,1,l) =     psilaxEzm(i  ,j  ,1,l) - FAC  *                   &
                            mcl(psilaxEzp(i  ,j  ,1,l) - psilaxEzp(i  ,j  ,1,l),  &
                                psilaxEzp(i  ,j  ,1,l) - psilaxEzp(i  ,j  ,1,l))

       BylaxEz (i,j,1,l)  = - 0.5d0 * (BylaxEzL (i,j,1,l) + BylaxEzR (i,j,1,l))
       BxlaxEz (i,j,1,l)  =   0.5d0 * (BxlaxEzL (i,j,1,l) + BxlaxEzR (i,j,1,l))
       psilaxEz(i,j,1,l)  =   0.5d0 * (psilaxEzL(i,j,1,l) + psilaxEzR(i,j,1,l))



! Conserved current flux
!_________________________________________________________________________

       JxlaxL(i,j,1,l)    =     Jxlaxp(i  ,j  ,1,l) + FAC  *                      &
                            mcl(Jxlaxp(i+1,j  ,1,l) - Jxlaxp(i  ,j  ,1,l),        &
                                Jxlaxp(i  ,j  ,1,l) - Jxlaxp(i-1,j  ,1,l))  
       JxlaxR(i,j,1,l)    =     Jxlaxm(i+1,j  ,1,l) -  FAC  *                     &
                            mcl(Jxlaxp(i+2,j  ,1,l) - Jxlaxp(i+1,j  ,1,l),        &
                                Jxlaxp(i+1,j  ,1,l) - Jxlaxp(i  ,j  ,1,l))

       JylaxL(i,j,1,l)    =     Jylaxp(i  ,j  ,1,l) +  FAC  *                     &
                            mcl(Jylaxp(i  ,j+1,1,l) - Jylaxp(i  ,j  ,1,l),        &
                                Jylaxp(i  ,j  ,1,l) - Jylaxp(i  ,j-1,1,l))
       JylaxR(i,j,1,l)    =     Jylaxm(i  ,j+1,1,l) -  FAC  *                     &
                            mcl(Jylaxp(i  ,j+2,1,l) - Jylaxp(i  ,j+1,1,l),        &
                                Jylaxp(i  ,j+1,1,l) - Jylaxp(i  ,j  ,1,l))

       JzlaxL(i,j,1,l)    =     Jzlaxp(i  ,j  ,1,l) +  FAC  *                     &
                            mcl(Jzlaxp(i  ,j  ,1,l) - Jzlaxp(i  ,j  ,1,l),        &
                                Jzlaxp(i  ,j  ,1,l) - Jzlaxp(i  ,j  ,1,l))
       JzlaxR(i,j,1,l)    =     Jzlaxm(i  ,j  ,1,l) -  FAC  *                     &
                            mcl(Jzlaxp(i  ,j  ,1,l) - Jzlaxp(i  ,j  ,1,l),        &
                                Jzlaxp(i  ,j  ,1,l) - Jzlaxp(i  ,j  ,1,l))

       Jxlax(i,j,1,l)     =   0.5d0 * (JxlaxL(i,j,1,l) + JxlaxR(i,j,1,l))
       Jylax(i,j,1,l)     =   0.5d0 * (JylaxL(i,j,1,l) + JylaxR(i,j,1,l))
       Jzlax(i,j,1,l)     =   0.5d0 * (JzlaxL(i,j,1,l) + JzlaxR(i,j,1,l))

! Conserved Mass flux
!_________________________________________________________________________

       FDxlaxL(i,j,1,l)   =     FDxlaxp(i  ,j  ,1,l) + FAC  *                     &
                            mcl(FDxlaxp(i+1,j  ,1,l) - FDxlaxp(i  ,j  ,1,l),      &
                                FDxlaxp(i  ,j  ,1,l) - FDxlaxp(i-1,j  ,1,l))   
       FDxlaxR(i,j,1,l)   =     FDxlaxm(i+1,j  ,1,l) -  FAC  *                    &
                            mcl(FDxlaxp(i+2,j  ,1,l) - FDxlaxp(i+1,j  ,1,l),      &
                                FDxlaxp(i+1,j  ,1,l) - FDxlaxp(i  ,j  ,1,l))    

       FDylaxL(i,j,1,l)   =     FDylaxp(i  ,j  ,1,l) +  FAC  *                    &
                            mcl(FDylaxp(i  ,j+1,1,l) - FDylaxp(i  ,j  ,1,l),      &
                                FDylaxp(i  ,j  ,1,l) - FDylaxp(i  ,j-1,1,l))
       FDylaxR(i,j,1,l)   =     FDylaxm(i  ,j+1,1,l) -  FAC  *                    &
                            mcl(FDylaxp(i  ,j+2,1,l) - FDylaxp(i  ,j+1,1,l),      &
                                FDylaxp(i  ,j+1,1,l) - FDylaxp(i  ,j  ,1,l))

       FDzlaxL(i,j,1,l)   =     FDzlaxp(i  ,j  ,1,l) +  FAC  *                    &
                            mcl(FDzlaxp(i  ,j  ,1,l) - FDzlaxp(i  ,j  ,1,l),      &
                                FDzlaxp(i  ,j  ,1,l) - FDzlaxp(i  ,j  ,1,l))
       FDzlaxR(i,j,1,l)   =     FDzlaxm(i  ,j  ,1,l) -  FAC  *                    &
                            mcl(FDzlaxp(i  ,j  ,1,l) - FDzlaxp(i  ,j  ,1,l),      &
                                FDzlaxp(i  ,j  ,1,l) - FDzlaxp(i  ,j  ,1,l))

       FDxlax(i,j,1,l)    =   0.5d0 * (FDxlaxL(i,j,1,l) + FDxlaxR(i,j,1,l))
       FDylax(i,j,1,l)    =   0.5d0 * (FDylaxL(i,j,1,l) + FDylaxR(i,j,1,l))
       FDzlax(i,j,1,l)    =   0.5d0 * (FDzlaxL(i,j,1,l) + FDzlaxR(i,j,1,l))

! Conserved Energy flux
!_________________________________________________________________________

       FtauxlaxL(i,j,1,l) =     Ftauxlaxp(i  ,j  ,1,l) + FAC *                    &
                            mcl(Ftauxlaxp(i+1,j  ,1,l) - Ftauxlaxp(i  ,j  ,1,l),  &
                                Ftauxlaxp(i  ,j  ,1,l) - Ftauxlaxp(i-1,j  ,1,l)) 
       FtauxlaxR(i,j,1,l) =     Ftauxlaxm(i+1,j  ,1,l) - FAC *                    &
                            mcl(Ftauxlaxp(i+2,j  ,1,l) - Ftauxlaxp(i+1,j  ,1,l),  &
                                Ftauxlaxp(i+1,j  ,1,l) - Ftauxlaxp(i  ,j  ,1,l)) 

       FtauylaxL(i,j,1,l) =     Ftauylaxp(i  ,j  ,1,l) +  FAC  *                  &
                            mcl(Ftauylaxp(i  ,j+1,1,l) - Ftauylaxp(i  ,j  ,1,l),  &
                                Ftauylaxp(i  ,j  ,1,l) - Ftauylaxp(i  ,j-1,1,l))
       FtauylaxR(i,j,1,l) =     Ftauylaxm(i  ,j+1,1,l) -  FAC  *                  &
                            mcl(Ftauylaxp(i  ,j+2,1,l) - Ftauylaxp(i  ,j+1,1,l),  &
                                Ftauylaxp(i  ,j+1,1,l) - Ftauylaxp(i  ,j  ,1,l))

       FtauzlaxL(i,j,1,l) =     Ftauzlaxp(i  ,j  ,1,l) +  FAC  *                  &
                            mcl(Ftauzlaxp(i  ,j  ,1,l) - Ftauzlaxp(i  ,j  ,1,l),  &
                                Ftauzlaxp(i  ,j  ,1,l) - Ftauzlaxp(i  ,j  ,1,l))
       FtauzlaxR(i,j,1,l) =     Ftauzlaxm(i  ,j  ,1,l) -  FAC  *                  &
                            mcl(Ftauzlaxp(i  ,j  ,1,l) - Ftauzlaxp(i  ,j  ,1,l),  &
                                Ftauzlaxp(i  ,j  ,1,l) - Ftauzlaxp(i  ,j  ,1,l))


       Ftauxlax(i,j,1,l)  =   0.5d0 * (FtauxlaxL(i,j,1,l) + FtauxlaxR(i,j,1,l))
       Ftauylax(i,j,1,l)  =   0.5d0 * (FtauylaxL(i,j,1,l) + FtauylaxR(i,j,1,l))
       Ftauzlax(i,j,1,l)  =   0.5d0 * (FtauzlaxL(i,j,1,l) + FtauzlaxR(i,j,1,l))

! -------------------------------------- EGLM--------------------------------------

       psilaxtauxL(i,j,1,l) =     psilaxtaup(i  ,j  ,1,l) +  FAC  *                   &
                              mcl(psilaxtaup(i+1,j  ,1,l) - psilaxtaup(i  ,j  ,1,l),  &
                                  psilaxtaup(i  ,j  ,1,l) - psilaxtaup(i-1,j  ,1,l))
       psilaxtauxR(i,j,1,l) =     psilaxtaum(i+1,j  ,1,l) -  FAC  *                   &
                              mcl(psilaxtaup(i+2,j  ,1,l) - psilaxtaup(i+1,j  ,1,l),  &
                                  psilaxtaup(i+1,j  ,1,l) - psilaxtaup(i  ,j  ,1,l))

       psilaxtauyL(i,j,1,l) =     psilaxtaup(i  ,j  ,1,l) +  FAC  *                   &
                              mcl(psilaxtaup(i  ,j+1,1,l) - psilaxtaup(i  ,j  ,1,l),  &
                                  psilaxtaup(i  ,j  ,1,l) - psilaxtaup(i  ,j-1,1,l))
       psilaxtauyR(i,j,1,l) =     psilaxtaum(i  ,j+1,1,l) -  FAC  *                   &
                              mcl(psilaxtaup(i  ,j+2,1,l) - psilaxtaup(i  ,j+1,1,l),  &
                                  psilaxtaup(i  ,j+1,1,l) - psilaxtaup(i  ,j  ,1,l))

       psilaxtauzL(i,j,1,l) =     psilaxtaup(i  ,j  ,1,l) +  FAC  *                   &
                              mcl(psilaxtaup(i  ,j  ,1,l) - psilaxtaup(i  ,j  ,1,l),  &
                                  psilaxtaup(i  ,j  ,1,l) - psilaxtaup(i  ,j  ,1,l))
       psilaxtauzR(i,j,1,l) =     psilaxtaum(i  ,j  ,1,l) -  FAC  *                   &
                              mcl(psilaxtaup(i  ,j  ,1,l) - psilaxtaup(i  ,j  ,1,l),  &
                                  psilaxtaup(i  ,j  ,1,l) - psilaxtaup(i  ,j  ,1,l))

       psitauxlax(i,j,1,l)  =   0.5d0 * (psilaxtauxL(i,j,1,l) + psilaxtauxR(i,j,1,l))
       psitauylax(i,j,1,l)  =   0.5d0 * (psilaxtauyL(i,j,1,l) + psilaxtauyR(i,j,1,l))
       psitauzlax(i,j,1,l)  =   0.5d0 * (psilaxtauzL(i,j,1,l) + psilaxtauzR(i,j,1,l))


       philaxtauxL(i,j,1,l) =     philaxtaup(i  ,j  ,1,l) +  FAC  *                   &
                              mcl(philaxtaup(i+1,j  ,1,l) - philaxtaup(i  ,j  ,1,l),  &
                                  philaxtaup(i  ,j  ,1,l) - philaxtaup(i-1,j  ,1,l))
       philaxtauxR(i,j,1,l) =     philaxtaum(i+1,j  ,1,l) -  FAC  *                   &
                              mcl(philaxtaup(i+2,j  ,1,l) - philaxtaup(i+1,j  ,1,l),  &
                                  philaxtaup(i+1,j  ,1,l) - philaxtaup(i  ,j  ,1,l))

       philaxtauyL(i,j,1,l) =     philaxtaup(i  ,j  ,1,l) +  FAC  *                   &
                              mcl(philaxtaup(i  ,j+1,1,l) - philaxtaup(i  ,j  ,1,l),  &
                                  philaxtaup(i  ,j  ,1,l) - philaxtaup(i  ,j-1,1,l))
       philaxtauyR(i,j,1,l) =     philaxtaum(i  ,j+1,1,l) -  FAC  *                   &
                              mcl(philaxtaup(i  ,j+2,1,l) - philaxtaup(i  ,j+1,1,l),  &
                                  philaxtaup(i  ,j+1,1,l) - philaxtaup(i  ,j  ,1,l))

       philaxtauzL(i,j,1,l) =     philaxtaup(i  ,j  ,1,l) +  FAC  *                   &
                              mcl(philaxtaup(i  ,j  ,1,l) - philaxtaup(i  ,j  ,1,l),  &
                                  philaxtaup(i  ,j  ,1,l) - philaxtaup(i  ,j  ,1,l))
       philaxtauzR(i,j,1,l) =     philaxtaum(i  ,j  ,1,l) -  FAC  *                   &
                              mcl(philaxtaup(i  ,j  ,1,l) - philaxtaup(i  ,j  ,1,l),  &
                                  philaxtaup(i  ,j  ,1,l) - philaxtaup(i  ,j  ,1,l))


       phitauxlax(i,j,1,l)  =   0.5d0 * (philaxtauxL(i,j,1,l) + philaxtauxR(i,j,1,l))
       phitauylax(i,j,1,l)  =   0.5d0 * (philaxtauyL(i,j,1,l) + philaxtauyR(i,j,1,l))
       phitauzlax(i,j,1,l)  =   0.5d0 * (philaxtauzL(i,j,1,l) + philaxtauzR(i,j,1,l))

! -------------------------------------- EGLM--------------------------------------


! Conserved Momentum flux tensor
!_________________________________________________________________________

       FSxxlaxL(i,j,1,l)  =     FSxxlaxp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(FSxxlaxp(i+1,j  ,1,l) - FSxxlaxp(i  ,j  ,1,l),    &
                                FSxxlaxp(i  ,j  ,1,l) - FSxxlaxp(i-1,j  ,1,l))  
       FSxxlaxR(i,j,1,l)  =     FSxxlaxm(i+1,j  ,1,l) - FAC  *                    &
                            mcl(FSxxlaxp(i+2,j  ,1,l) - FSxxlaxp(i+1,j  ,1,l),    &
                                FSxxlaxp(i+1,j  ,1,l) - FSxxlaxp(i  ,j  ,1,l))
       FSxylaxL(i,j,1,l)  =     FSxylaxp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(FSxylaxp(i+1,j  ,1,l) - FSxylaxp(i  ,j  ,1,l),    &
                                FSxylaxp(i  ,j  ,1,l) - FSxylaxp(i-1,j  ,1,l))  
       FSxylaxR(i,j,1,l)  =     FSxylaxm(i+1,j  ,1,l) - FAC  *                    &
                            mcl(FSxylaxp(i+2,j  ,1,l) - FSxylaxp(i+1,j  ,1,l),    &
                                FSxylaxp(i+1,j  ,1,l) - FSxylaxp(i  ,j  ,1,l))
       FSxzlaxL(i,j,1,l)  =     FSxzlaxp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(FSxzlaxp(i+1,j  ,1,l) - FSxzlaxp(i  ,j  ,1,l),    &
                                FSxzlaxp(i  ,j  ,1,l) - FSxzlaxp(i-1,j  ,1,l))  
       FSxzlaxR(i,j,1,l)  =     FSxzlaxm(i+1,j  ,1,l) - FAC  *                    &
                            mcl(FSxzlaxp(i+2,j  ,1,l) - FSxzlaxp(i+1,j  ,1,l),    &
                                FSxzlaxp(i+1,j  ,1,l) - FSxzlaxp(i  ,j  ,1,l))


       FSxxlax(i,j,1,l)   =   0.5d0 * (FSxxlaxL(i,j,1,l)  + FSxxlaxR(i,j,1,l)  ) 
       FSxylax(i,j,1,l)   =   0.5d0 * (FSxylaxL(i,j,1,l)  + FSxylaxR(i,j,1,l)  )
       FSxzlax(i,j,1,l)   =   0.5d0 * (FSxzlaxL(i,j,1,l)  + FSxzlaxR(i,j,1,l)  )


       FSyxlaxL(i,j,1,l)  =     FSyxlaxp(i  ,j  ,1,l) +  FAC  *                   &
                            mcl(FSyxlaxp(i  ,j+1,1,l) - FSyxlaxp(i  ,j  ,1,l),    &
                                FSyxlaxp(i  ,j  ,1,l) - FSyxlaxp(i  ,j-1,1,l))
       FSyxlaxR(i,j,1,l)  =     FSyxlaxm(i  ,j+1,1,l) -  FAC  *                   &
                            mcl(FSyxlaxp(i  ,j+2,1,l) - FSyxlaxp(i  ,j+1,1,l),    &
                                FSyxlaxp(i  ,j+1,1,l) - FSyxlaxp(i  ,j  ,1,l))
       FSyylaxL(i,j,1,l)  =     FSyylaxp(i  ,j  ,1,l) +  FAC  *                   &
                            mcl(FSyylaxp(i  ,j+1,1,l) - FSyylaxp(i  ,j  ,1,l),    &
                                FSyylaxp(i  ,j  ,1,l) - FSyylaxp(i  ,j-1,1,l))
       FSyylaxR(i,j,1,l)  =     FSyylaxm(i  ,j+1,1,l) -  FAC  *                   &
                            mcl(FSyylaxp(i  ,j+2,1,l) - FSyylaxp(i  ,j+1,1,l),    &
                                FSyylaxp(i  ,j+1,1,l) - FSyylaxp(i  ,j  ,1,l))
       FSyzlaxL(i,j,1,l)  =     FSyzlaxp(i  ,j  ,1,l) +  FAC  *                   &
                            mcl(FSyzlaxp(i  ,j+1,1,l) - FSyzlaxp(i  ,j  ,1,l),    &
                                FSyzlaxp(i  ,j  ,1,l) - FSyzlaxp(i  ,j-1,1,l))
       FSyzlaxR(i,j,1,l)  =     FSyzlaxm(i  ,j+1,1,l) -  FAC  *                   &
                            mcl(FSyzlaxp(i  ,j+2,1,l) - FSyzlaxp(i  ,j+1,1,l),    &
                                FSyzlaxp(i  ,j+1,1,l) - FSyzlaxp(i  ,j  ,1,l))
 
       FSyxlax(i,j,1,l)   =   0.5d0 * (FSyxlaxL(i,j,1,l)  + FSyxlaxR(i,j,1,l)  )
       FSyylax(i,j,1,l)   =   0.5d0 * (FSyylaxL(i,j,1,l)  + FSyylaxR(i,j,1,l)  )
       FSyzlax(i,j,1,l)   =   0.5d0 * (FSyzlaxL(i,j,1,l)  + FSyzlaxR(i,j,1,l)  )


       FSzxlaxL(i,j,1,l)  =     FSzxlaxp(i  ,j  ,1,l) +  FAC  *                   &
                            mcl(FSzxlaxp(i  ,j  ,1,l) - FSzxlaxp(i  ,j  ,1,l),    &
                                FSzxlaxp(i  ,j  ,1,l) - FSzxlaxp(i  ,j  ,1,l))
       FSzxlaxR(i,j,1,l)  =     FSzxlaxm(i  ,j  ,1,l) -  FAC  *                   &
                            mcl(FSzxlaxp(i  ,j  ,1,l) - FSzxlaxp(i  ,j  ,1,l),    &
                                FSzxlaxp(i  ,j  ,1,l) - FSzxlaxp(i  ,j  ,1,l))
       FSzylaxL(i,j,1,l)  =     FSzylaxp(i  ,j  ,1,l) +  FAC  *                   &
                            mcl(FSzylaxp(i  ,j  ,1,l) - FSzylaxp(i  ,j  ,1,l),    &
                                FSzylaxp(i  ,j  ,1,l) - FSzylaxp(i  ,j  ,1,l))
       FSzylaxR(i,j,1,l)  =     FSzylaxm(i  ,j  ,1,l) -  FAC  *                   &
                            mcl(FSzylaxp(i  ,j  ,1,l) - FSzylaxp(i  ,j  ,1,l),    &
                                FSzylaxp(i  ,j  ,1,l) - FSzylaxp(i  ,j  ,1,l))
       FSzzlaxL(i,j,1,l)  =     FSzzlaxp(i  ,j  ,1,l) +  FAC  *                   &
                            mcl(FSzzlaxp(i  ,j  ,1,l) - FSzzlaxp(i  ,j  ,1,l),    &
                                FSzzlaxp(i  ,j  ,1,l) - FSzzlaxp(i  ,j  ,1,l))
       FSzzlaxR(i,j,1,l)  =     FSzzlaxm(i  ,j  ,1,l) -  FAC  *                   &
                            mcl(FSzzlaxp(i  ,j  ,1,l) - FSzzlaxp(i  ,j  ,1,l),    &
                                FSzzlaxp(i  ,j  ,1,l) - FSzzlaxp(i  ,j  ,1,l))

       FSzxlax(i,j,1,l)  =   0.5d0 * (FSzxlaxL(i,j,1,l)  + FSzxlaxR(i,j,1,l)  )
       FSzylax(i,j,1,l)  =   0.5d0 * (FSzylaxL(i,j,1,l)  + FSzylaxR(i,j,1,l)  )
       FSzzlax(i,j,1,l)  =   0.5d0 * (FSzzlaxL(i,j,1,l)  + FSzzlaxR(i,j,1,l)  )

! -------------------------------------- EGLM--------------------------------------

       BxSxlaxL(i,j,1,l)  =     BxlaxSxp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BxlaxSxp(i+1,j  ,1,l) - BxlaxSxp(i  ,j  ,1,l),    &
                                BxlaxSxp(i  ,j  ,1,l) - BxlaxSxp(i-1,j  ,1,l))  
       BxSxlaxR(i,j,1,l)  =     BxlaxSxm(i+1,j  ,1,l) - FAC  *                    &
                            mcl(BxlaxSxp(i+2,j  ,1,l) - BxlaxSxp(i+1,j  ,1,l),    &
                                BxlaxSxp(i+1,j  ,1,l) - BxlaxSxp(i  ,j  ,1,l))
       BySxlaxL(i,j,1,l)  =     BylaxSxp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BylaxSxp(i  ,j+1,1,l) - BylaxSxp(i  ,j  ,1,l),    &
                                BylaxSxp(i  ,j  ,1,l) - BylaxSxp(i  ,j-1,1,l))  
       BySxlaxR(i,j,1,l)  =     BylaxSxm(i  ,j+1,1,l) - FAC  *                    &
                            mcl(BylaxSxp(i  ,j+2,1,l) - BylaxSxp(i  ,j+1,1,l),    &
                                BylaxSxp(i  ,j+1,1,l) - BylaxSxp(i  ,j  ,1,l))
       BzSxlaxL(i,j,1,l)  =     BzlaxSxp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BzlaxSxp(i  ,j  ,1,l) - BzlaxSxp(i  ,j  ,1,l),    &
                                BzlaxSxp(i  ,j  ,1,l) - BzlaxSxp(i  ,j  ,1,l))  
       BzSxlaxR(i,j,1,l)  =     BzlaxSxm(i  ,j  ,1,l) - FAC  *                    &
                            mcl(BzlaxSxp(i  ,j  ,1,l) - BzlaxSxp(i  ,j  ,1,l),    &
                                BzlaxSxp(i  ,j  ,1,l) - BzlaxSxp(i  ,j  ,1,l))

       BxSxlax(i,j,1,l)   = 0.5d0 * (BxSxlaxL(i,j,1,l)  + BxSxlaxR(i,j,1,l)    )
       BySxlax(i,j,1,l)   = 0.5d0 * (BySxlaxL(i,j,1,l)  + BySxlaxR(i,j,1,l)    )
       BzSxlax(i,j,1,l)   = 0.5d0 * (BzSxlaxL(i,j,1,l)  + BzSxlaxR(i,j,1,l)    )


       BxSylaxL(i,j,1,l)  =     BxlaxSyp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BxlaxSyp(i+1,j  ,1,l) - BxlaxSyp(i  ,j  ,1,l),    &
                                BxlaxSyp(i  ,j  ,1,l) - BxlaxSyp(i-1,j  ,1,l))  
       BxSylaxR(i,j,1,l)  =     BxlaxSym(i+1,j  ,1,l) - FAC  *                    &
                            mcl(BxlaxSyp(i+2,j  ,1,l) - BxlaxSyp(i+1,j  ,1,l),    &
                                BxlaxSyp(i+1,j  ,1,l) - BxlaxSyp(i  ,j  ,1,l))
       BySylaxL(i,j,1,l)  =     BylaxSyp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BylaxSyp(i  ,j+1,1,l) - BylaxSyp(i  ,j  ,1,l),    &
                                BylaxSyp(i  ,j  ,1,l) - BylaxSyp(i  ,j-1,1,l))  
       BySylaxR(i,j,1,l)  =     BylaxSym(i  ,j+1,1,l) - FAC  *                    &
                            mcl(BylaxSyp(i  ,j+2,1,l) - BylaxSyp(i  ,j+1,1,l),    &
                                BylaxSyp(i  ,j+1,1,l) - BylaxSyp(i  ,j  ,1,l))
       BzSylaxL(i,j,1,l)  =     BzlaxSyp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BzlaxSyp(i  ,j  ,1,l) - BzlaxSyp(i  ,j  ,1,l),    &
                                BzlaxSyp(i  ,j  ,1,l) - BzlaxSyp(i  ,j  ,1,l))  
       BzSylaxR(i,j,1,l)  =     BzlaxSym(i  ,j  ,1,l) - FAC  *                    &
                            mcl(BzlaxSyp(i  ,j  ,1,l) - BzlaxSyp(i  ,j  ,1,l),    &
                                BzlaxSyp(i  ,j  ,1,l) - BzlaxSyp(i  ,j  ,1,l))

       BxSylax(i,j,1,l)   = 0.5d0 * (BxSylaxL(i,j,1,l)  + BxSylaxR(i,j,1,l)    )
       BySylax(i,j,1,l)   = 0.5d0 * (BySylaxL(i,j,1,l)  + BySylaxR(i,j,1,l)    )
       BzSylax(i,j,1,l)   = 0.5d0 * (BzSylaxL(i,j,1,l)  + BzSylaxR(i,j,1,l)    )


       BxSzlaxL(i,j,1,l)  =     BxlaxSzp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BxlaxSzp(i+1,j  ,1,l) - BxlaxSzp(i  ,j  ,1,l),    &
                                BxlaxSzp(i  ,j  ,1,l) - BxlaxSzp(i-1,j  ,1,l))  
       BxSzlaxR(i,j,1,l)  =     BxlaxSzm(i+1,j  ,1,l) - FAC  *                    &
                            mcl(BxlaxSzp(i+2,j  ,1,l) - BxlaxSzp(i+1,j  ,1,l),    &
                                BxlaxSzp(i+1,j  ,1,l) - BxlaxSzp(i  ,j  ,1,l))
       BySzlaxL(i,j,1,l)  =     BylaxSzp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BylaxSzp(i  ,j+1,1,l) - BylaxSzp(i  ,j  ,1,l),    &
                                BylaxSzp(i  ,j  ,1,l) - BylaxSzp(i  ,j-1,1,l))  
       BySzlaxR(i,j,1,l)  =     BylaxSzm(i  ,j+1,1,l) - FAC  *                    &
                            mcl(BylaxSzp(i  ,j+2,1,l) - BylaxSzp(i  ,j+1,1,l),    &
                                BylaxSzp(i  ,j+1,1,l) - BylaxSzp(i  ,j  ,1,l))
       BzSzlaxL(i,j,1,l)  =     BzlaxSzp(i  ,j  ,1,l) + FAC  *                    &
                            mcl(BzlaxSzp(i  ,j  ,1,l) - BzlaxSzp(i  ,j  ,1,l),    &
                                BzlaxSzp(i  ,j  ,1,l) - BzlaxSzp(i  ,j  ,1,l))  
       BzSzlaxR(i,j,1,l)  =     BzlaxSzm(i  ,j  ,1,l) - FAC  *                    &
                            mcl(BzlaxSzp(i  ,j  ,1,l) - BzlaxSzp(i  ,j  ,1,l),    &
                                BzlaxSzp(i  ,j  ,1,l) - BzlaxSzp(i  ,j  ,1,l))

       BxSzlax(i,j,1,l)   = 0.5d0 * (BxSzlaxL(i,j,1,l)  + BxSzlaxR(i,j,1,l)    )
       BySzlax(i,j,1,l)   = 0.5d0 * (BySzlaxL(i,j,1,l)  + BySzlaxR(i,j,1,l)    )
       BzSzlax(i,j,1,l)   = 0.5d0 * (BzSzlaxL(i,j,1,l)  + BzSzlaxR(i,j,1,l)    )

! -------------------------------------- EGLM--------------------------------------

     end do ! for i
   end do ! for j

!$OMP END DO

   else

      write(*,*) "STOP: subroutine laxflow"
      write(*,*) "This dimension is not implemented yet"
      stop

   end if

!$OMP END PARALLEL

 end subroutine laxflow
