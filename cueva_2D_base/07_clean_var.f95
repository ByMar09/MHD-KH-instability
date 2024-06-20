  !     *******************************************************************
  !     Subroutine which Clean Sums and Intermediate Variables 
  !     *******************************************************************

  subroutine cleanvarint 


    use scalar
    use parameters
    use threevectors
    use fourvectors

    implicit none

!$OMP PARALLEL
!$OMP DO 

    do l=1,lmax
       do j=0+fix,jmax
          do i=0,imax

                 !     ------------------------------------
                 !     Starts Intermediate Sums
                 !     ------------------------------------
                     
     psiintsum(i,j,1) = 0.d0
     phiintsum(i,j,1) = 0.d0
     Bxintsum(i,j,1)  = 0.d0
     Byintsum(i,j,1)  = 0.d0
     Bzintsum(i,j,1)  = 0.d0
     Exastsum(i,j,1)  = 0.d0
     Eyastsum(i,j,1)  = 0.d0
     Ezastsum(i,j,1)  = 0.d0
     qintsum(i,j,1)   = 0.d0
     Dintsum(i,j,1)   = 0.d0
     tauintsum(i,j,1) = 0.d0
     Sxintsum(i,j,1)  = 0.d0
     Syintsum(i,j,1)  = 0.d0
     Szintsum(i,j,1)  = 0.d0

              !     ------------------------------------
              !     Starts Final Sums
              !     ------------------------------------

     psifinsum(i,j,1) = 0.d0
     phifinsum(i,j,1) = 0.d0
     Bxfinsum(i,j,1)  = 0.d0
     Byfinsum(i,j,1)  = 0.d0
     Bzfinsum(i,j,1)  = 0.d0
     Exfinsum(i,j,1)  = 0.d0
     Eyfinsum(i,j,1)  = 0.d0
     Ezfinsum(i,j,1)  = 0.d0
     qfinsum(i,j,1)   = 0.d0
     Dfinsum(i,j,1)   = 0.d0
     taufinsum(i,j,1) = 0.d0
     Sxfinsum(i,j,1)  = 0.d0
     Syfinsum(i,j,1)  = 0.d0
     Szfinsum(i,j,1)  = 0.d0

  !     Let''s be paranoid: init every AUXILIAR variable to zero

     psiint(i,j,1,l)  = 0.d0
     phiint(i,j,1,l)  = 0.d0
     Exast (i,j,1,l)  = 0.d0
     Eyast (i,j,1,l)  = 0.d0
     Ezast (i,j,1,l)  = 0.d0
     Exint(i,j,1,l)   = 0.d0
     Eyint(i,j,1,l)   = 0.d0
     Ezint(i,j,1,l)   = 0.d0
     Bxint(i,j,1,l)   = 0.d0
     Byint(i,j,1,l)   = 0.d0
     Bzint(i,j,1,l)   = 0.d0
     qint(i,j,1,l)    = 0.d0
     DDint(i,j,1,l)   = 0.d0
     tauint(i,j,1,l)  = 0.d0
     Sxint(i,j,1,l)   = 0.d0
     Syint(i,j,1,l)   = 0.d0
     Szint(i,j,1,l)   = 0.d0


     Exlaxpsip(i,j,1,l) = 0.d0 
     Exlaxpsim(i,j,1,l) = 0.d0
     Eylaxpsip(i,j,1,l) = 0.d0
     Eylaxpsim(i,j,1,l) = 0.d0
     Ezlaxpsip(i,j,1,l) = 0.d0
     Ezlaxpsim(i,j,1,l) = 0.d0

     Bxlaxphip(i,j,1,l) = 0.d0
     Bxlaxphim(i,j,1,l) = 0.d0
     Bylaxphip(i,j,1,l) = 0.d0
     Bylaxphim(i,j,1,l) = 0.d0
     Bzlaxphip(i,j,1,l) = 0.d0
     Bzlaxphim(i,j,1,l) = 0.d0

     philaxBxp(i,j,1,l) = 0.d0
     philaxBxm(i,j,1,l) = 0.d0
     EzlaxBxp (i,j,1,l) = 0.d0
     EzlaxBxm (i,j,1,l) = 0.d0
     EylaxBxp (i,j,1,l) = 0.d0
     EylaxBxm (i,j,1,l) = 0.d0

     EzlaxByp (i,j,1,l) = 0.d0
     EzlaxBym (i,j,1,l) = 0.d0
     philaxByp(i,j,1,l) = 0.d0
     philaxBym(i,j,1,l) = 0.d0
     ExlaxByp (i,j,1,l) = 0.d0
     ExlaxBym (i,j,1,l) = 0.d0

     EylaxBzp (i,j,1,l) = 0.d0
     EylaxBzm (i,j,1,l) = 0.d0
     ExlaxBzp (i,j,1,l) = 0.d0
     ExlaxBzm (i,j,1,l) = 0.d0
     philaxBzp(i,j,1,l) = 0.d0
     philaxBzm(i,j,1,l) = 0.d0

     psilaxExp(i,j,1,l) = 0.d0
     psilaxExm(i,j,1,l) = 0.d0
     BzlaxExp (i,j,1,l) = 0.d0
     BzlaxExm (i,j,1,l) = 0.d0
     BylaxExp (i,j,1,l) = 0.d0
     BylaxExm (i,j,1,l) = 0.d0

     BzlaxEyp (i,j,1,l) = 0.d0
     BzlaxEym (i,j,1,l) = 0.d0
     psilaxEyp(i,j,1,l) = 0.d0
     psilaxEym(i,j,1,l) = 0.d0
     BxlaxEyp (i,j,1,l) = 0.d0
     BxlaxEym (i,j,1,l) = 0.d0

     BylaxEzp (i,j,1,l) = 0.d0
     BylaxEzm (i,j,1,l) = 0.d0
     BxlaxEzp (i,j,1,l) = 0.d0
     BxlaxEzm (i,j,1,l) = 0.d0
     psilaxEzp(i,j,1,l) = 0.d0
     psilaxEzm(i,j,1,l) = 0.d0

     Jxlaxp(i,j,1,l)    = 0.d0
     Jxlaxm(i,j,1,l)    = 0.d0
     Jylaxp(i,j,1,l)    = 0.d0
     Jylaxm(i,j,1,l)    = 0.d0
     Jzlaxp(i,j,1,l)    = 0.d0
     Jzlaxm(i,j,1,l)    = 0.d0

     FDxlaxp(i,j,1,l)   = 0.d0
     FDxlaxm(i,j,1,l)   = 0.d0
     FDylaxp(i,j,1,l)   = 0.d0
     FDylaxm(i,j,1,l)   = 0.d0
     FDzlaxp(i,j,1,l)   = 0.d0
     FDzlaxm(i,j,1,l)   = 0.d0

     Ftauxlaxp(i,j,1,l) = 0.d0
     Ftauxlaxm(i,j,1,l) = 0.d0
     Ftauylaxp(i,j,1,l) = 0.d0
     Ftauylaxm(i,j,1,l) = 0.d0
     Ftauzlaxp(i,j,1,l) = 0.d0
     Ftauzlaxm(i,j,1,l) = 0.d0

     FSxxlaxp(i,j,1,l)  = 0.d0
     FSxxlaxm(i,j,1,l)  = 0.d0
     FSxylaxp(i,j,1,l)  = 0.d0
     FSxylaxm(i,j,1,l)  = 0.d0
     FSxzlaxp(i,j,1,l)  = 0.d0
     FSxzlaxm(i,j,1,l)  = 0.d0

     FSyxlaxp(i,j,1,l)  = 0.d0
     FSyxlaxm(i,j,1,l)  = 0.d0
     FSyylaxp(i,j,1,l)  = 0.d0
     FSyylaxm(i,j,1,l)  = 0.d0
     FSyzlaxp(i,j,1,l)  = 0.d0
     FSyzlaxm(i,j,1,l)  = 0.d0

     FSzxlaxp(i,j,1,l)  = 0.d0
     FSzxlaxm(i,j,1,l)  = 0.d0
     FSzylaxp(i,j,1,l)  = 0.d0
     FSzylaxm(i,j,1,l)  = 0.d0
     FSzzlaxp(i,j,1,l)  = 0.d0
     FSzzlaxm(i,j,1,l)  = 0.d0

     Exlaxpsi(i,j,1,l) =  0.d0
     Eylaxpsi(i,j,1,l) =  0.d0
     Ezlaxpsi(i,j,1,l) =  0.d0

     Bxlaxphi(i,j,1,l) =  0.d0
     Bylaxphi(i,j,1,l) =  0.d0
     Bzlaxphi(i,j,1,l) =  0.d0

     philaxBx(i,j,1,l) =  0.d0
     EzlaxBx (i,j,1,l) =  0.d0
     EylaxBx (i,j,1,l) =  0.d0

     EzlaxBy (i,j,1,l) =  0.d0
     philaxBy(i,j,1,l) =  0.d0
     ExlaxBy (i,j,1,l) =  0.d0

     EylaxBz (i,j,1,l) =  0.d0
     ExlaxBz (i,j,1,l) =  0.d0
     philaxBz(i,j,1,l) =  0.d0

     psilaxEx(i,j,1,l) =  0.d0
     BzlaxEx (i,j,1,l) =  0.d0
     BylaxEx (i,j,1,l) =  0.d0

     BzlaxEy (i,j,1,l) =  0.d0
     psilaxEy(i,j,1,l) =  0.d0
     BxlaxEy (i,j,1,l) =  0.d0

     BylaxEz (i,j,1,l) =  0.d0
     BxlaxEz (i,j,1,l) =  0.d0
     psilaxEz(i,j,1,l) =  0.d0

     Jxlax   (i,j,1,l) =  0.d0
     Jylax   (i,j,1,l) =  0.d0
     Jzlax   (i,j,1,l) =  0.d0

     FDxlax  (i,j,1,l) =  0.d0
     FDylax  (i,j,1,l) =  0.d0
     FDzlax  (i,j,1,l) =  0.d0

     Ftauxlax(i,j,1,l) =  0.d0
     Ftauylax(i,j,1,l) =  0.d0
     Ftauzlax(i,j,1,l) =  0.d0

     FSxxlax (i,j,1,l) =  0.d0
     FSxylax (i,j,1,l) =  0.d0
     FSxzlax (i,j,1,l) =  0.d0
     FSyxlax (i,j,1,l) =  0.d0
     FSyylax (i,j,1,l) =  0.d0
     FSyzlax (i,j,1,l) =  0.d0
     FSzxlax (i,j,1,l) =  0.d0
     FSzylax (i,j,1,l) =  0.d0
     FSzzlax (i,j,1,l) =  0.d0

     ExlaxpsiL(i,j,1,l) =  0.d0
     ExlaxpsiR(i,j,1,l) =  0.d0
     EylaxpsiL(i,j,1,l) =  0.d0
     EylaxpsiR(i,j,1,l) =  0.d0
     EzlaxpsiL(i,j,1,l) =  0.d0
     EzlaxpsiR(i,j,1,l) =  0.d0

     BxlaxphiL(i,j,1,l) =  0.d0
     BxlaxphiR(i,j,1,l) =  0.d0
     BylaxphiL(i,j,1,l) =  0.d0
     BylaxphiR(i,j,1,l) =  0.d0
     BzlaxphiL(i,j,1,l) =  0.d0
     BzlaxphiR(i,j,1,l) =  0.d0


     philaxBxL(i,j,1,l) =  0.d0 
     philaxBxR(i,j,1,l) =  0.d0 
     EzlaxBxL(i,j,1,l)  =  0.d0 
     EzlaxBxR(i,j,1,l)  =  0.d0 
     EylaxBxL(i,j,1,l)  =  0.d0 
     EylaxBxR(i,j,1,l)  =  0.d0 


     EzlaxByL(i,j,1,l)  =  0.d0 
     EzlaxByR(i,j,1,l)  =  0.d0 
     philaxByL(i,j,1,l) =  0.d0 
     philaxByR(i,j,1,l) =  0.d0 
     ExlaxByL(i,j,1,l)  =  0.d0 
     ExlaxByR(i,j,1,l)  =  0.d0 

     EylaxBzL(i,j,1,l)  =  0.d0   
     EylaxBzR(i,j,1,l)  =  0.d0 
     ExlaxBzL(i,j,1,l)  =  0.d0 
     ExlaxBzR(i,j,1,l)  =  0.d0 
     philaxBzL(i,j,1,l) =  0.d0 
     philaxBzR(i,j,1,l) =  0.d0 
    
     psilaxExL(i,j,1,l) =  0.d0
     psilaxExR(i,j,1,l) =  0.d0
     BzlaxExL(i,j,1,l)  =  0.d0
     BzlaxExR(i,j,1,l)  =  0.d0
     BylaxExL(i,j,1,l)  =  0.d0
     BylaxExR(i,j,1,l)  =  0.d0

     BzlaxEyL(i,j,1,l)  =  0.d0 
     BzlaxEyR(i,j,1,l)  =  0.d0 
     psilaxEyL(i,j,1,l) =  0.d0 
     psilaxEyR(i,j,1,l) =  0.d0 
     BxlaxEyL(i,j,1,l)  =  0.d0 
     BxlaxEyR(i,j,1,l)  =  0.d0 
    
     BylaxEzL(i,j,1,l)  =  0.d0 
     BylaxEzR(i,j,1,l)  =  0.d0 
     BxlaxEzL(i,j,1,l)  =  0.d0 
     BxlaxEzR(i,j,1,l)  =  0.d0 
     psilaxEzL(i,j,1,l) =  0.d0 
     psilaxEzR(i,j,1,l) =  0.d0
    
     JxlaxL(i,j,1,l)    =  0.d0
     JxlaxR(i,j,1,l)    =  0.d0
     JylaxL(i,j,1,l)    =  0.d0
     JylaxR(i,j,1,l)    =  0.d0
     JzlaxL(i,j,1,l)    =  0.d0
     JzlaxR(i,j,1,l)    =  0.d0

     FDxlaxL(i,j,1,l)   =  0.d0 
     FDxlaxR(i,j,1,l)   =  0.d0 
     FDylaxL(i,j,1,l)   =  0.d0
     FDylaxR(i,j,1,l)   =  0.d0 
     FDzlaxL(i,j,1,l)   =  0.d0 
     FDzlaxR(i,j,1,l)   =  0.d0 

     FtauxlaxL(i,j,1,l) =  0.d0
     FtauxlaxR(i,j,1,l) =  0.d0
     FtauylaxL(i,j,1,l) =  0.d0
     FtauylaxR(i,j,1,l) =  0.d0
     FtauzlaxL(i,j,1,l) =  0.d0
     FtauzlaxR(i,j,1,l) =  0.d0

     FSxxlaxL(i,j,1,l)  =  0.d0
     FSxxlaxR(i,j,1,l)  =  0.d0
     FSxylaxL(i,j,1,l)  =  0.d0
     FSxylaxR(i,j,1,l)  =  0.d0
     FSxzlaxL(i,j,1,l)  =  0.d0
     FSxzlaxR(i,j,1,l)  =  0.d0

     FSyxlaxL(i,j,1,l)  =  0.d0
     FSyxlaxR(i,j,1,l)  =  0.d0
     FSyylaxL(i,j,1,l)  =  0.d0
     FSyylaxR(i,j,1,l)  =  0.d0
     FSyzlaxL(i,j,1,l)  =  0.d0
     FSyzlaxR(i,j,1,l)  =  0.d0

     FSzxlaxL(i,j,1,l)  =  0.d0 
     FSzxlaxR(i,j,1,l)  =  0.d0
     FSzylaxL(i,j,1,l)  =  0.d0 
     FSzylaxR(i,j,1,l)  =  0.d0 
     FSzzlaxL(i,j,1,l)  =  0.d0 
     FSzzlaxR(i,j,1,l)  =  0.d0 


     E2int(i,j,1,l)    =  0.d0  
     B2int(i,j,1,l)    =  0.d0   
     Edotv             =  0.d0 
     vrotB_x           =  0.d0  
     vrotB_y           =  0.d0 
     vrotB_z           =  0.d0  

     Jxint(i,j,1,l)    =  0.d0  
     Jyint(i,j,1,l)    =  0.d0  
     Jzint(i,j,1,l)    =  0.d0  
     
     FDxint(i,j,1,l)   =  0.d0 
     FDyint(i,j,1,l)   =  0.d0 
     FDzint(i,j,1,l)   =  0.d0 
     
     Ftauxint(i,j,1,l) =  0.d0  
     Ftauyint(i,j,1,l) =  0.d0  
     Ftauzint(i,j,1,l) =  0.d0  

     FSxxint(i,j,1,l)  =  0.d0  
     FSxyint(i,j,1,l)  =  0.d0 
     FSxzint(i,j,1,l)  =  0.d0 
     FSyxint(i,j,1,l)  =  0.d0  
     FSyyint(i,j,1,l)  =  0.d0 
     FSyzint(i,j,1,l)  =  0.d0  
     FSzxint(i,j,1,l)  =  0.d0  
     FSzyint(i,j,1,l)  =  0.d0  
     FSzzint(i,j,1,l)  =  0.d0  

     psiflux(i,j,1,l)  =  0.d0  
     phiflux(i,j,1,l)  =  0.d0  
     
     Bxflux(i,j,1,l)   =  0.d0  
     Byflux(i,j,1,l)   =  0.d0  
     Bzflux(i,j,1,l)   =  0.d0  

     Exastflux(i,j,1,l)=  0.d0  
     Exastsour(i,j,1,l)=  0.d0  
     Eyastflux(i,j,1,l)=  0.d0  
     Eyastsour(i,j,1,l)=  0.d0  
     Ezastflux(i,j,1,l)=  0.d0  
     Ezastsour(i,j,1,l)=  0.d0  
     
     qflux(i,j,1,l)    =  0.d0  
     
     Dflux(i,j,1,l)    =  0.d0  
     
     tauflux(i,j,1,l)  =  0.d0  
     taueglm(i,j,1,l)  =  0.d0  
     
     Sxflux(i,j,1,l)   =  0.d0  
     Syflux(i,j,1,l)   =  0.d0  
     Szflux(i,j,1,l)   =  0.d0  
     Sxeglm(i,j,1,l)   =  0.d0  
     Syeglm(i,j,1,l)   =  0.d0  
     Szeglm(i,j,1,l)   =  0.d0  


              end do !i
        end do !j
     end do !l

!$OMP END DO
!$OMP END PARALLEL

  end subroutine cleanvarint
