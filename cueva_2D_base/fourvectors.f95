

    module fourvectors
  
    use parameters, only: imax, jmax, kmax, lmax
    implicit none 
!    save

  !     Intermediate Primitive Variables
  !     ------------  

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Vxint, Vyint, Vzint, Wint
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: pint, rhoint, enthpyint, epsilonint

   !     Intermediate Augmented Maxwell Fields
  !     ------------  

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psiint, phiint

  !      Intermediate Vectorial Fields
  !     ------------  
  !     Electric (Exast,Eyast,Ezast), Magnetic (Bxint,Byint,Bzint),  
  !     Current density (Jxint,Jyint,Jzint), Mass flux (FDint = \rho
  !      W V), 
  !     Energy flux (Ftauint), Momentum (Sxint,Syint,Szint)

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Exast, Eyast, Ezast
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Bxint, Byint, Bzint
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Jxint, Jyint, Jzint    
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FDxint, FDyint, FDzint
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Ftauxint, Ftauyint, Ftauzint
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Sxint, Syint, Szint

  !     Intermediate  Momentum  Flux Tensor 
  !     ------------  
  !          | FSxxint FSxyint FSxzint |
  !    FS =  | FSyxint FSyyint FSyzint |
  !          | FSzxint FSzyint FSzzint |

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSxxint, FSxyint, FSxzint 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSyxint, FSyyint, FSyzint
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSzxint, FSzyint, FSzzint

  !     Intermediate Scalar Fields
  !     ------------
  !     Density (DDint= \rho W),  Energy density (\tauint),  
  !     Electric field magnitude (Eint^2),  Magnetic (Bint^2), Charge
  !      (qint)

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: DDint, qint, tauint
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: E2int, B2int

  !     Sums of explicit intermediate values
  !     ------------
  !     ------------
  !     IMPLICIT AUXILIAR INTERMEDIATE VALUES, "STIFF TERM" 
  !     ------------

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Exint, Eyint, Ezint

  !     ------------
  !     SIMPLER LAX-FRIEDRICH FLUXS
  !     ------------

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Exlaxpsi,Eylaxpsi,Ezlaxpsi 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Exlaxpsip,Eylaxpsip,Ezlaxpsip
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Exlaxpsim,Eylaxpsim,Ezlaxpsim
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: ExlaxpsiL,EylaxpsiL,EzlaxpsiL  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: ExlaxpsiR,EylaxpsiR,EzlaxpsiR  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Bxlaxphi,Bylaxphi,Bzlaxphi
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Bxlaxphip,Bylaxphip,Bzlaxphip
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Bxlaxphim,Bylaxphim,Bzlaxphim
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxlaxphiL,BylaxphiL,BzlaxphiL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxlaxphiR,BylaxphiR,BzlaxphiR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: philaxBx,EzlaxBx,EylaxBx
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: philaxBxp,EzlaxBxp,EylaxBxp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: philaxBxm,EzlaxBxm,EylaxBxm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: philaxBxL,EzlaxBxL,EylaxBxL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: philaxBxR,EzlaxBxR,EylaxBxR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EzlaxBy,philaxBy,ExlaxBy
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EzlaxByp,philaxByp,ExlaxByp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EzlaxBym,philaxBym,ExlaxBym
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EzlaxByL,philaxByL,ExlaxByL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EzlaxByR,philaxByR,ExlaxByR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EylaxBz,ExlaxBz,philaxBz
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EylaxBzp,ExlaxBzp,philaxBzp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EylaxBzm,ExlaxBzm,philaxBzm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EylaxBzL,ExlaxBzL,philaxBzL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: EylaxBzR,ExlaxBzR,philaxBzR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psilaxEx,BzlaxEx,BylaxEx
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psilaxExp,BzlaxExp,BylaxExp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psilaxExm,BzlaxExm,BylaxExm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psilaxExL,BzlaxExL,BylaxExL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psilaxExR,BzlaxExR,BylaxExR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BzlaxEy,psilaxEy,BxlaxEy
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BzlaxEyp,psilaxEyp,BxlaxEyp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BzlaxEym,psilaxEym,BxlaxEym
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BzlaxEyL,psilaxEyL,BxlaxEyL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BzlaxEyR,psilaxEyR,BxlaxEyR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BylaxEz,BxlaxEz,psilaxEz
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BylaxEzp,BxlaxEzp,psilaxEzp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BylaxEzm,BxlaxEzm,psilaxEzm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BylaxEzL,BxlaxEzL,psilaxEzL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BylaxEzR,BxlaxEzR,psilaxEzR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Jxlax,Jylax,Jzlax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Jxlaxp,Jylaxp,Jzlaxp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Jxlaxm,Jylaxm,Jzlaxm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: JxlaxL,JylaxL,JzlaxL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: JxlaxR,JylaxR,JzlaxR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FDxlax,FDylax,FDzlax 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FDxlaxp,FDylaxp,FDzlaxp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FDxlaxm,FDylaxm,FDzlaxm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FDxlaxL,FDylaxL,FDzlaxL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FDxlaxR,FDylaxR,FDzlaxR    
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Ftauxlax, Ftauylax, Ftauzlax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Ftauxlaxp, Ftauylaxp, Ftauzlaxp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Ftauxlaxm, Ftauylaxm, Ftauzlaxm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FtauxlaxL, FtauylaxL, FtauzlaxL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FtauxlaxR, FtauylaxR, FtauzlaxR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSxxlax, FSxylax, FSxzlax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSxxlaxp, FSxylaxp, FSxzlaxp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSxxlaxm, FSxylaxm, FSxzlaxm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSxxlaxL, FSxylaxL, FSxzlaxL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSxxlaxR, FSxylaxR, FSxzlaxR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSyxlax, FSyylax, FSyzlax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSyxlaxp, FSyylaxp, FSyzlaxp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSyxlaxm, FSyylaxm, FSyzlaxm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSyxlaxL, FSyylaxL, FSyzlaxL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSyxlaxR, FSyylaxR, FSyzlaxR
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSzxlax, FSzylax, FSzzlax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSzxlaxp, FSzylaxp, FSzzlaxp
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSzxlaxm, FSzylaxm, FSzzlaxm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSzxlaxL, FSzylaxL, FSzzlaxL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: FSzxlaxR, FSzylaxR, FSzzlaxR

  !     ------------
  !     FLUX AND SOURCES
  !     ------------

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psiflux,phiflux
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Bxflux, Byflux, Bzflux
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Exastflux, Eyastflux, Ezastflux
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Exastsour, Eyastsour, Ezastsour
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: qflux,  Dflux, tauflux
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Sxflux, Syflux, Szflux


  !     ------------
  !     Extended Generalized Lagrange Multiplier (EGLM)
  !     ------------

    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: taueglm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: Sxeglm, Syeglm, Szeglm
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psilaxtaup, psilaxtaum
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: philaxtaup, philaxtaum
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxlaxSxp,BxlaxSxm,BylaxSxp,BylaxSxm,BzlaxSxp,BzlaxSxm 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxlaxSyp,BxlaxSym,BylaxSyp,BylaxSym,BzlaxSyp,BzlaxSym 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxlaxSzp,BxlaxSzm,BylaxSzp,BylaxSzm,BzlaxSzp,BzlaxSzm 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psilaxtauxL,psilaxtauxR,psilaxtauyL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psilaxtauyR,psilaxtauzL,psilaxtauzR 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: psitauxlax, psitauylax, psitauzlax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: philaxtauxL,philaxtauxR,philaxtauyL
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: philaxtauyR,philaxtauzL,philaxtauzR 
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: phitauxlax, phitauylax, phitauzlax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxSxlaxL,BxSxlaxR,BySxlaxL,BySxlaxR,BzSxlaxL,BzSxlaxR  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxSxlax,BySxlax,BzSxlax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxSylaxL,BxSylaxR,BySylaxL,BySylaxR,BzSylaxL,BzSylaxR  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxSylax,BySylax,BzSylax
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxSzlaxL,BxSzlaxR,BySzlaxL,BySzlaxR,BzSzlaxL,BzSzlaxR  
    doubleprecision,dimension(-6:imax+6,-6:jmax+6,1:kmax,1:lmax) :: BxSzlax,BySzlax,BzSzlax



  end module fourvectors
