

  !     *******************************************************************
  !     Subroutine which Reconstruct variables
  !     *******************************************************************
  subroutine recons_var_primitive

    use scalar
    use parameters
    use threevectors 
    use funciones, only: limiter
	
    implicit none

!$OMP   PARALLEL PRIVATE(psim1,phim1,Bxm1,Bym1,Bzm1,Exm1,Eym1,Ezm1,Vxm1,Vym1,Vzm1,qm1,pm1,rhom1),      &
!$OMP & PRIVATE(psim,phim,Bxm,Bym,Bzm,Exm,Eym,Ezm,Vxm,Vym,Vzm,qm,pm,rhom),                             &
!$OMP & PRIVATE(psip1,phip1,Bxp1,Byp1,Bzp1,Exp1,Eyp1,Ezp1,Vxp1,Vyp1,Vzp1,qp1,pp1,rhop1),               &
!$OMP & PRIVATE(psip2,phip2,Bxp2,Byp2,Bzp2,Exp2,Eyp2,Ezp2,Vxp2,Vyp2,Vzp2,qp2,pp2,rhop2),               &
!$OMP & PRIVATE(sppsi1,spphi1,spBx1,spBy1,spBz1,spEx1,spEy1,spEz1,spVx1,spVy1,spVz1,spq1,spp1,sprho1), &
!$OMP & PRIVATE(smpsi1,smphi1,smBx1,smBy1,smBz1,smEx1,smEy1,smEz1,smVx1,smVy1,smVz1,smq1,smp1,smrho1), &
!$OMP & PRIVATE(sppsi2,spphi2,spBx2,spBy2,spBz2,spEx2,spEy2,spEz2,spVx2,spVy2,spVz2,spq2,spp2,sprho2), &
!$OMP & PRIVATE(smpsi2,smphi2,smBx2,smBy2,smBz2,smEx2,smEy2,smEz2,smVx2,smVy2,smVz2,smq2,smp2,smrho2), &
!$OMP & PRIVATE(psipp,phipp,Bxpp,Bypp,Bzpp,Expp,Eypp,Ezpp,Vxpp,Vypp,Vzpp,qpp,ppp,rhopp),               &
!$OMP & PRIVATE(psimm,phimm,Bxmm,Bymm,Bzmm,Exmm,Eymm,Ezmm,Vxmm,Vymm,Vzmm,qmm,pmm,rhomm)
   
    
     if (DIM == 1) then

        if (recons_method  == 0) then

!$OMP DO
           
         do i=-2,imax

            psim1 = psi(i-1,1,1)
            phim1 = phi(i-1,1,1)
            Bxm1  = Bx (i-1,1,1)
            Bym1  = By (i-1,1,1)
            Bzm1  = Bz (i-1,1,1)
            Exm1  = Ex (i-1,1,1)
            Eym1  = Ey (i-1,1,1)
            Ezm1  = Ez (i-1,1,1)
            Vxm1  = Vx (i-1,1,1)
            Vym1  = Vy (i-1,1,1)
            Vzm1  = Vz (i-1,1,1)
            qm1   = q  (i-1,1,1)
            pm1   = p  (i-1,1,1)
            rhom1 = rho(i-1,1,1)

            psim  = psi(i  ,1,1)
            phim  = phi(i  ,1,1)
            Bxm   = Bx (i  ,1,1)
            Bym   = By (i  ,1,1)
            Bzm   = Bz (i  ,1,1)
            Exm   = Ex (i  ,1,1)
            Eym   = Ey (i  ,1,1)
            Ezm   = Ez (i  ,1,1)
            Vxm   = Vx (i  ,1,1)
            Vym   = Vy (i  ,1,1)
            Vzm   = Vz (i  ,1,1)
            qm    = q  (i  ,1,1)
            pm    = p  (i  ,1,1)
            rhom  = rho(i  ,1,1)

            psip1 = psi(i+1,1,1)
            phip1 = phi(i+1,1,1)
            Bxp1  = Bx (i+1,1,1)
            Byp1  = By (i+1,1,1)
            Bzp1  = Bz (i+1,1,1)
            Exp1  = Ex (i+1,1,1)
            Eyp1  = Ey (i+1,1,1)
            Ezp1  = Ez (i+1,1,1)
            Vxp1  = Vx (i+1,1,1)
            Vyp1  = Vy (i+1,1,1)
            Vzp1  = Vz (i+1,1,1)
            qp1   = q  (i+1,1,1)
            pp1   = p  (i+1,1,1)
            rhop1 = rho(i+1,1,1)

            psip2 = psi(i+2,1,1)
            phip2 = phi(i+2,1,1)
            Bxp2  = Bx (i+2,1,1)
            Byp2  = By (i+2,1,1)
            Bzp2  = Bz (i+2,1,1)
            Exp2  = Ex (i+2,1,1)
            Eyp2  = Ey (i+2,1,1)
            Ezp2  = Ez (i+2,1,1)
            Vxp2  = Vx (i+2,1,1)
            Vyp2  = Vy (i+2,1,1)
            Vzp2  = Vz (i+2,1,1)
            qp2   = q  (i+2,1,1)
            pp2   = p  (i+2,1,1)
            rhop2 = rho(i+2,1,1)

            sppsi1  = psip1 - psim
            spphi1  = phip1 - phim
            spBx1   = Bxp1  - Bxm
            spBy1   = Byp1  - Bym
            spBz1   = Bzp1  - Bzm
            spEx1   = Exp1  - Exm
            spEy1   = Eyp1  - Eym
            spEz1   = Ezp1  - Ezm
            spVx1   = Vxp1  - Vxm
            spVy1   = Vyp1  - Vym
            spVz1   = Vzp1  - Vzm
            spq1    = qp1   - qm
            spp1    = pp1   - pm
            sprho1  = rhop1 - rhom


            smpsi1  = psim - psim1
            smphi1  = phim - phim1
            smBx1   = Bxm  - Bxm1
            smBy1   = Bym  - Bym1
            smBz1   = Bzm  - Bzm1
            smEx1   = Exm  - Exm1
            smEy1   = Eym  - Eym1
            smEz1   = Ezm  - Ezm1
            smVx1   = Vxm  - Vxm1
            smVy1   = Vym  - Vym1
            smVz1   = Vzm  - Vzm1
            smq1    = qm   - qm1
            smp1    = pm   - pm1
            smrho1  = rhom - rhom1

            sppsi2  = psip2 - psip1
            spphi2  = phip2 - phip1
            spBx2   = Bxp2  - Bxp1
            spBy2   = Byp2  - Byp1
            spBz2   = Bzp2  - Bzp1
            spEx2   = Exp2  - Exp1
            spEy2   = Eyp2  - Eyp1
            spEz2   = Ezp2  - Ezp1
            spVx2   = Vxp2  - Vxp1
            spVy2   = Vyp2  - Vyp1
            spVz2   = Vzp2  - Vzp1
            spq2    = qp2   - qp1
            spp2    = pp2   - pp1
            sprho2  = rhop2 - rhop1

            smpsi2  = psip1 - psim
            smphi2  = phip1 - phim
            smBx2   = Bxp1  - Bxm
            smBy2   = Byp1  - Bym
            smBz2   = Bzp1  - Bzm
            smEx2   = Exp1  - Exm
            smEy2   = Eyp1  - Eym
            smEz2   = Ezp1  - Ezm
            smVx2   = Vxp1  - Vxm
            smVy2   = Vyp1  - Vym
            smVz2   = Vzp1  - Vzm
            smq2    = qp1   - qm
            smp2    = pp1   - pm
            smrho2  = rhop1 - rhom


            psipp  = psim + 0.5d0 * limiter(sppsi1,smpsi1)
            phipp  = phim + 0.5d0 * limiter(spphi1,smphi1)
            Bxpp   = Bxm  + 0.5d0 * limiter(spBx1 ,smBx1 )
            Bypp   = Bym  + 0.5d0 * limiter(spBy1 ,smBy1 )
            Bzpp   = Bzm  + 0.5d0 * limiter(spBz1 ,smBz1 )
            Expp   = Exm  + 0.5d0 * limiter(spEx1 ,smEx1 )
            Eypp   = Eym  + 0.5d0 * limiter(spEy1 ,smEy1 )
            Ezpp   = Ezm  + 0.5d0 * limiter(spEz1 ,smEz1 )
            Vxpp   = Vxm  + 0.5d0 * limiter(spVx1 ,smVx1 )
            Vypp   = Vym  + 0.5d0 * limiter(spVy1 ,smVy1 )
            Vzpp   = Vzm  + 0.5d0 * limiter(spVz1 ,smVz1 )
            qpp    = qm   + 0.5d0 * limiter(spq1  ,smq1  )
            ppp    = pm   + 0.5d0 * limiter(spp1  ,smp1  ) 
            rhopp  = rhom + 0.5d0 * limiter(sprho1,smrho1)

            psimm  = psip1 - 0.5d0 * limiter(sppsi2,smpsi2)
            phimm  = phip1 - 0.5d0 * limiter(spphi2,smphi2)
            Bxmm   = Bxp1  - 0.5d0 * limiter(spBx2 ,smBx2 )
            Bymm   = Byp1  - 0.5d0 * limiter(spBy2 ,smBy2 )
            Bzmm   = Bzp1  - 0.5d0 * limiter(spBz2 ,smBz2 )
            Exmm   = Exp1  - 0.5d0 * limiter(spEx2 ,smEx2 )
            Eymm   = Eyp1  - 0.5d0 * limiter(spEy2 ,smEy2 )
            Ezmm   = Ezp1  - 0.5d0 * limiter(spEz2 ,smEz2 )
            Vxmm   = Vxp1  - 0.5d0 * limiter(spVx2 ,smVx2 )
            Vymm   = Vyp1  - 0.5d0 * limiter(spVy2 ,smVy2 )
            Vzmm   = Vzp1  - 0.5d0 * limiter(spVz2 ,smVz2 )
            qmm    = qp1   - 0.5d0 * limiter(spq2  ,smq2  )
            pmm    = pp1   - 0.5d0 * limiter(spp2  ,smp2  ) 
            rhomm  = rhop1 - 0.5d0 * limiter(sprho2,smrho2)

            psiL(i  ,1,1)  = psipp 
            phiL(i  ,1,1)  = phipp 
            BxL (i  ,1,1)  = Bxpp  
            ByL (i  ,1,1)  = Bypp  
            BzL (i  ,1,1)  = Bzpp  
            ExL (i  ,1,1)  = Expp  
            EyL (i  ,1,1)  = Eypp  
            EzL (i  ,1,1)  = Ezpp  
            VxL (i  ,1,1)  = Vxpp  
            VyL (i  ,1,1)  = Vypp  
            VzL (i  ,1,1)  = Vzpp  
            qL  (i  ,1,1)  = qpp   
            pL  (i  ,1,1)  = ppp   
            rhoL(i  ,1,1)  = rhopp 

            psiR(i+1,1,1)  = psimm 
            phiR(i+1,1,1)  = phimm 
            BxR (i+1,1,1)  = Bxmm  
            ByR (i+1,1,1)  = Bymm  
            BzR (i+1,1,1)  = Bzmm  
            ExR (i+1,1,1)  = Exmm  
            EyR (i+1,1,1)  = Eymm  
            EzR (i+1,1,1)  = Ezmm  
            VxR (i+1,1,1)  = Vxmm  
            VyR (i+1,1,1)  = Vymm  
            VzR (i+1,1,1)  = Vzmm  
            qR  (i+1,1,1)  = qmm   
            pR  (i+1,1,1)  = pmm   
            rhoR(i+1,1,1)  = rhomm 

         end do
!$OMP END DO
 
!$OMP DO 
         do i = 0, imax

           psi(i,1,1) = 0.5d0 * ( psiR(i,1,1) + psiL(i,1,1) )
           phi(i,1,1) = 0.5d0 * ( phiR(i,1,1) + phiL(i,1,1) )
           Bx(i,1,1)  = 0.5d0 * ( BxR(i,1,1)  + BxL (i,1,1) )
           By(i,1,1)  = 0.5d0 * ( ByR(i,1,1)  + ByL (i,1,1) )
           Bz(i,1,1)  = 0.5d0 * ( BzR(i,1,1)  + BzL (i,1,1) )
           Ex(i,1,1)  = 0.5d0 * ( ExR(i,1,1)  + ExL (i,1,1) )
           Ey(i,1,1)  = 0.5d0 * ( EyR(i,1,1)  + EyL (i,1,1) )
           Ez(i,1,1)  = 0.5d0 * ( EzR(i,1,1)  + EzL (i,1,1) )
           Vx(i,1,1)  = 0.5d0 * ( VxR(i,1,1)  + VxL (i,1,1) )
           Vy(i,1,1)  = 0.5d0 * ( VyR(i,1,1)  + VyL (i,1,1) )
           Vz(i,1,1)  = 0.5d0 * ( VzR(i,1,1)  + VzL (i,1,1) )
           q(i,1,1)   = 0.5d0 * ( qR(i,1,1)   + qL  (i,1,1) )
           p(i,1,1)   = 0.5d0 * ( pR(i,1,1)   + pL  (i,1,1) )
           rho(i,1,1) = 0.5d0 * ( rhoR(i,1,1) + rhoL(i,1,1) )

         end do

!$OMP END DO

             else if (recons_method  == 1) then 

!!$!//////////////////////////// PSI ///////////////////////////////
!!$
!!$         varname = "not__"
!!$
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = psi(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            psiL(i,1,1) = up(i) 
!!$            psiR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// PSI ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) =phi(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            phiL(i,1,1) = up(i) 
!!$            phiR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// Bx   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) =Bx(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            BxL(i,1,1) = up(i) 
!!$            BxR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// By   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = By(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            ByL(i,1,1) = up(i) 
!!$            ByR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// Bz   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = Bz(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            BzL(i,1,1) = up(i) 
!!$            BzR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// Ex   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = Ex(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            ExL(i,1,1) = up(i) 
!!$            ExR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// Ey   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = Ey(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            EyL(i,1,1) = up(i) 
!!$            EyR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// Ez   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = Ez(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            EzL(i,1,1) = up(i) 
!!$            EzR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// Vx   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = Vx(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            VxL(i,1,1) = up(i) 
!!$            VxR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// Vy   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = Vy(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            VyL(i,1,1) = up(i) 
!!$            VyR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// Vz   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = Vz(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            VzL(i,1,1) = up(i) 
!!$            VzR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// q   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = q(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            qL(i,1,1) = up(i) 
!!$            qR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// p   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = p(i,1,1)
!!$
!!$         end do
!!$
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            pL(i,1,1) = up(i) 
!!$            pR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$!//////////////////////////// rho   ///////////////////////////////
!!$
!!$         do i=-4,imax+4
!!$
!!$            u(i) = rho(i,1,1)
!!$
!!$         end do
!!$
!!$         varname = "rho__"
!!$         call MP5
!!$
!!$         do i=-1,imax+1
!!$
!!$            rhoL(i,1,1) = up(i) 
!!$            rhoR(i,1,1) = um(i) 
!!$
!!$         end do
!!$
!!$         varname = "not__"
!!$
!!$!//////////////////////////// rho   ///////////////////////////////

!!$!$OMP DO 
!!$
!!$ 
!!$         do i = 0, imax
!!$
!!$           psi(i,1,1) = 0.5d0 * ( psiR(i,1,1) + psiL(i,1,1) ) 
!!$           phi(i,1,1) = 0.5d0 * ( phiR(i,1,1) + phiL(i,1,1) ) 
!!$           Bx(i,1,1)  = 0.5d0 * ( BxR (i,1,1) + BxL (i,1,1) ) 
!!$           By(i,1,1)  = 0.5d0 * ( ByR (i,1,1) + ByL (i,1,1) ) 
!!$           Bz(i,1,1)  = 0.5d0 * ( BzR (i,1,1) + BzL (i,1,1) ) 
!!$           Ex(i,1,1)  = 0.5d0 * ( ExR (i,1,1) + ExL (i,1,1) ) 
!!$           Ey(i,1,1)  = 0.5d0 * ( EyR (i,1,1) + EyL (i,1,1) ) 
!!$           Ez(i,1,1)  = 0.5d0 * ( EzR (i,1,1) + EzL (i,1,1) ) 
!!$           Vx(i,1,1)  = 0.5d0 * ( VxR (i,1,1) + VxL (i,1,1) ) 
!!$           Vy(i,1,1)  = 0.5d0 * ( VyR (i,1,1) + VyL (i,1,1) ) 
!!$           Vz(i,1,1)  = 0.5d0 * ( VzR (i,1,1) + VzL (i,1,1) ) 
!!$           q(i,1,1)   = 0.5d0 * ( qR  (i,1,1) + qL  (i,1,1) ) 
!!$           p(i,1,1)   = 0.5d0 * ( pR  (i,1,1) + pL  (i,1,1) ) 
!!$           rho(i,1,1) = 0.5d0 * ( rhoR(i,1,1) + rhoL(i,1,1) ) 
!!$
!!$         end do
!!$
!!$!$OMP END DO

             else

                write(*,*)  "STOP: subroutine recons_var_primitive"
                write(*,*)  " this recons_method parameter is not valid"
                stop

             end if




          else if (DIM == 2) then


!$OMP DO              

         do j=-2,jmax
            do i=-2,imax


            psim1 = psi(i-1,j-1,1)
            phim1 = phi(i-1,j-1,1)
            Bxm1  = Bx (i-1,j-1,1)
            Bym1  = By (i-1,j-1,1)
            Bzm1  = Bz (i-1,j-1,1)
            Exm1  = Ex (i-1,j-1,1)
            Eym1  = Ey (i-1,j-1,1)
            Ezm1  = Ez (i-1,j-1,1)
            Vxm1  = Vx (i-1,j-1,1)
            Vym1  = Vy (i-1,j-1,1)
            Vzm1  = Vz (i-1,j-1,1)
            qm1   = q  (i-1,j-1,1)
            pm1   = p  (i-1,j-1,1)
            rhom1 = rho(i-1,j-1,1)

            psim  = psi(i  ,j  ,1)
            phim  = phi(i  ,j  ,1)
            Bxm   = Bx (i  ,j  ,1)
            Bym   = By (i  ,j  ,1)
            Bzm   = Bz (i  ,j  ,1)
            Exm   = Ex (i  ,j  ,1)
            Eym   = Ey (i  ,j  ,1)
            Ezm   = Ez (i  ,j  ,1)
            Vxm   = Vx (i  ,j  ,1)
            Vym   = Vy (i  ,j  ,1)
            Vzm   = Vz (i  ,j  ,1)
            qm    = q  (i  ,j  ,1)
            pm    = p  (i  ,j  ,1)
            rhom  = rho(i  ,j  ,1)

            psip1 = psi(i+1,j+1,1)
            phip1 = phi(i+1,j+1,1)
            Bxp1  = Bx (i+1,j+1,1)
            Byp1  = By (i+1,j+1,1)
            Bzp1  = Bz (i+1,j+1,1)
            Exp1  = Ex (i+1,j+1,1)
            Eyp1  = Ey (i+1,j+1,1)
            Ezp1  = Ez (i+1,j+1,1)
            Vxp1  = Vx (i+1,j+1,1)
            Vyp1  = Vy (i+1,j+1,1)
            Vzp1  = Vz (i+1,j+1,1)
            qp1   = q  (i+1,j+1,1)
            pp1   = p  (i+1,j+1,1)
            rhop1 = rho(i+1,j+1,1)

            psip2 = psi(i+2,j+2,1)
            phip2 = phi(i+2,j+2,1)
            Bxp2  = Bx (i+2,j+2,1)
            Byp2  = By (i+2,j+2,1)
            Bzp2  = Bz (i+2,j+2,1)
            Exp2  = Ex (i+2,j+2,1)
            Eyp2  = Ey (i+2,j+2,1)
            Ezp2  = Ez (i+2,j+2,1)
            Vxp2  = Vx (i+2,j+2,1)
            Vyp2  = Vy (i+2,j+2,1)
            Vzp2  = Vz (i+2,j+2,1)
            qp2   = q  (i+2,j+2,1)
            pp2   = p  (i+2,j+2,1)
            rhop2 = rho(i+2,j+2,1)

            sppsi1  = psip1 - psim
            spphi1  = phip1 - phim
            spBx1   = Bxp1  - Bxm
            spBy1   = Byp1  - Bym
            spBz1   = Bzp1  - Bzm
            spEx1   = Exp1  - Exm
            spEy1   = Eyp1  - Eym
            spEz1   = Ezp1  - Ezm
            spVx1   = Vxp1  - Vxm
            spVy1   = Vyp1  - Vym
            spVz1   = Vzp1  - Vzm
            spq1    = qp1   - qm
            spp1    = pp1   - pm
            sprho1  = rhop1 - rhom


            smpsi1  = psim - psim1
            smphi1  = phim - phim1
            smBx1   = Bxm  - Bxm1
            smBy1   = Bym  - Bym1
            smBz1   = Bzm  - Bzm1
            smEx1   = Exm  - Exm1
            smEy1   = Eym  - Eym1
            smEz1   = Ezm  - Ezm1
            smVx1   = Vxm  - Vxm1
            smVy1   = Vym  - Vym1
            smVz1   = Vzm  - Vzm1
            smq1    = qm   - qm1
            smp1    = pm   - pm1
            smrho1  = rhom - rhom1

            sppsi2  = psip2 - psip1
            spphi2  = phip2 - phip1
            spBx2   = Bxp2  - Bxp1
            spBy2   = Byp2  - Byp1
            spBz2   = Bzp2  - Bzp1
            spEx2   = Exp2  - Exp1
            spEy2   = Eyp2  - Eyp1
            spEz2   = Ezp2  - Ezp1
            spVx2   = Vxp2  - Vxp1
            spVy2   = Vyp2  - Vyp1
            spVz2   = Vzp2  - Vzp1
            spq2    = qp2   - qp1
            spp2    = pp2   - pp1
            sprho2  = rhop2 - rhop1

            smpsi2  = psip1 - psim
            smphi2  = phip1 - phim
            smBx2   = Bxp1  - Bxm
            smBy2   = Byp1  - Bym
            smBz2   = Bzp1  - Bzm
            smEx2   = Exp1  - Exm
            smEy2   = Eyp1  - Eym
            smEz2   = Ezp1  - Ezm
            smVx2   = Vxp1  - Vxm
            smVy2   = Vyp1  - Vym
            smVz2   = Vzp1  - Vzm
            smq2    = qp1   - qm
            smp2    = pp1   - pm
            smrho2  = rhop1 - rhom


            psipp  = psim + 0.5d0 * limiter(sppsi1,smpsi1)
            phipp  = phim + 0.5d0 * limiter(spphi1,smphi1)
            Bxpp   = Bxm  + 0.5d0 * limiter(spBx1 ,smBx1 )
            Bypp   = Bym  + 0.5d0 * limiter(spBy1 ,smBy1 )
            Bzpp   = Bzm  + 0.5d0 * limiter(spBz1 ,smBz1 )
            Expp   = Exm  + 0.5d0 * limiter(spEx1 ,smEx1 )
            Eypp   = Eym  + 0.5d0 * limiter(spEy1 ,smEy1 )
            Ezpp   = Ezm  + 0.5d0 * limiter(spEz1 ,smEz1 )
            Vxpp   = Vxm  + 0.5d0 * limiter(spVx1 ,smVx1 )
            Vypp   = Vym  + 0.5d0 * limiter(spVy1 ,smVy1 )
            Vzpp   = Vzm  + 0.5d0 * limiter(spVz1 ,smVz1 )
            qpp    = qm   + 0.5d0 * limiter(spq1  ,smq1  )
            ppp    = pm   + 0.5d0 * limiter(spp1  ,smp1  ) 
            rhopp  = rhom + 0.5d0 * limiter(sprho1,smrho1)

            psimm  = psip1 - 0.5d0 * limiter(sppsi2,smpsi2)
            phimm  = phip1 - 0.5d0 * limiter(spphi2,smphi2)
            Bxmm   = Bxp1  - 0.5d0 * limiter(spBx2 ,smBx2 )
            Bymm   = Byp1  - 0.5d0 * limiter(spBy2 ,smBy2 )
            Bzmm   = Bzp1  - 0.5d0 * limiter(spBz2 ,smBz2 )
            Exmm   = Exp1  - 0.5d0 * limiter(spEx2 ,smEx2 )
            Eymm   = Eyp1  - 0.5d0 * limiter(spEy2 ,smEy2 )
            Ezmm   = Ezp1  - 0.5d0 * limiter(spEz2 ,smEz2 )
            Vxmm   = Vxp1  - 0.5d0 * limiter(spVx2 ,smVx2 )
            Vymm   = Vyp1  - 0.5d0 * limiter(spVy2 ,smVy2 )
            Vzmm   = Vzp1  - 0.5d0 * limiter(spVz2 ,smVz2 )
            qmm    = qp1   - 0.5d0 * limiter(spq2  ,smq2  )
            pmm    = pp1   - 0.5d0 * limiter(spp2  ,smp2  ) 
            rhomm  = rhop1 - 0.5d0 * limiter(sprho2,smrho2)

            psiL(i  ,j  ,1)  = psipp 
            phiL(i  ,j  ,1)  = phipp 
            BxL (i  ,j  ,1)  = Bxpp  
            ByL (i  ,j  ,1)  = Bypp  
            BzL (i  ,j  ,1)  = Bzpp  
            ExL (i  ,j  ,1)  = Expp  
            EyL (i  ,j  ,1)  = Eypp  
            EzL (i  ,j  ,1)  = Ezpp  
            VxL (i  ,j  ,1)  = Vxpp  
            VyL (i  ,j  ,1)  = Vypp  
            VzL (i  ,j  ,1)  = Vzpp  
            qL  (i  ,j  ,1)  = qpp   
            pL  (i  ,j  ,1)  = ppp   
            rhoL(i  ,j  ,1)  = rhopp 

            psiR(i+1,j+1,1)  = psimm 
            phiR(i+1,j+1,1)  = phimm 
            BxR (i+1,j+1,1)  = Bxmm  
            ByR (i+1,j+1,1)  = Bymm  
            BzR (i+1,j+1,1)  = Bzmm  
            ExR (i+1,j+1,1)  = Exmm  
            EyR (i+1,j+1,1)  = Eymm  
            EzR (i+1,j+1,1)  = Ezmm  
            VxR (i+1,j+1,1)  = Vxmm  
            VyR (i+1,j+1,1)  = Vymm  
            VzR (i+1,j+1,1)  = Vzmm  
            qR  (i+1,j+1,1)  = qmm   
            pR  (i+1,j+1,1)  = pmm   
            rhoR(i+1,j+1,1)  = rhomm 


         end do
      end do

!$OMP END DO      


!$OMP DO 

         do j=0,jmax
            do i=0,imax

           psi(i,j,1) = 0.5d0 * ( psiR(i,j,1) + psiL(i,j,1) )
           phi(i,j,1) = 0.5d0 * ( phiR(i,j,1) + phiL(i,j,1) )
           Bx(i,j,1)  = 0.5d0 * ( BxR (i,j,1) + BxL (i,j,1) )
           By(i,j,1)  = 0.5d0 * ( ByR (i,j,1) + ByL (i,j,1) )
           Bz(i,j,1)  = 0.5d0 * ( BzR (i,j,1) + BzL (i,j,1) )
           Ex(i,j,1)  = 0.5d0 * ( ExR (i,j,1) + ExL (i,j,1) )
           Ey(i,j,1)  = 0.5d0 * ( EyR (i,j,1) + EyL (i,j,1) )
           Ez(i,j,1)  = 0.5d0 * ( EzR (i,j,1) + EzL (i,j,1) )
           Vx(i,j,1)  = 0.5d0 * ( VxR (i,j,1) + VxL (i,j,1) )
           Vy(i,j,1)  = 0.5d0 * ( VyR (i,j,1) + VyL (i,j,1) )
           Vz(i,j,1)  = 0.5d0 * ( VzR (i,j,1) + VzL (i,j,1) )
           q(i,j,1)   = 0.5d0 * ( qR  (i,j,1) + qL  (i,j,1) )
           p(i,j,1)   = 0.5d0 * ( pR  (i,j,1) + pL  (i,j,1) )
           rho(i,j,1) = 0.5d0 * ( rhoR(i,j,1) + rhoL(i,j,1) )

        end do
     end do

!$OMP END DO
     
    else

       write(*,*)  "STOP: subroutine recons_var_primitive"
       write(*,*)  " this dimension is not implemented yet"
       stop

    end if

!$OMP END PARALLEL    

  end subroutine recons_var_primitive

