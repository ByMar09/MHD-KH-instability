
      subroutine MP9(imp5)  !Tomek: 02/07/17

        use parameters, only: imax, jmax, kmax, mm, nn, DIM
        use scalar    , only: u, um, up, du, d2, vorp, vorm, vmpp, vmpm, coordenate
        use funciones
!*
!*     This surbroutine implements the MP5 method descrived in Suresh & Huynh (1997; JCP, 136, 83-99; SH97)
!*

!*     u: cell avergages 
!*     up: interface value at i+1/2
!*     um: interface value at i-1/2
!*
      implicit none

!*
!*     The value of alpha should limit the maximum CFL to be used: CFL <= 1/(1+alpha).
!*     Depending on the test, sometimes this limit can be overshooted. However, in tests which may develop
!*     very small presures/densities, it is a good idea to reduce alpha to alpha=2.
!*

!Tomek      

      doubleprecision, parameter :: b1 = 1.d0/2520.d0, b2 = 4.d0/3.d0, epsm = 1.d-10, alpha_mp = 2.d0 !   alpha_mp = 1.d0/3.d0 ! original alpha_mp =2.d0 ! MP9
      
!Tomek
      
!*
!*     Local variables
!*

!Tomek
            integer                                    :: im4, im3, im2, im1, ip1, ip2, ip3, ip4, i0
!Tomek
      
      integer                                    :: imp, jmp, jmp_max
      integer, intent(in)                        :: imp5
!      doubleprecision, dimension(-4:imax+4,1:68) :: du, d2, vorp, vorm, vmpp, vmpm
      doubleprecision                            :: dm4jph, dm4jmh, vul, vav, vmd, vlc, vmin, vmax
!*      character*5                           :: varname

!*
!*     Arguments
!*
      
      !*      doubleprecision, dimension(-4:imax+4)  :: u, um, up

      mm=0
      nn=imax

      if (DIM == 1) then

         jmp_max = 43

      else if (DIM == 2) then

         jmp_max = 68

      end if

      if (coordenate .eq. 1) then

         nn=imax

      else if (coordenate .eq. 2) then

         nn=jmax

      endif



      do imp=mm-5,nn+5
         do jmp = 1, jmp_max

         du(imp,jmp,imp5) = u(imp  ,jmp,imp5) - u(imp-1,jmp,imp5)
         d2(imp,jmp,imp5) = u(imp-1,jmp,imp5) - 2d0 * u(imp,jmp,imp5) + u(imp+1,jmp,imp5)

         !         print*, du(imp,jmp,imp), d2(imp,jmp,imp5), i
      end do
   end do


   do imp=mm-2,nn+2    !Tomek: I use a dirty trick to evade a ghostzone /stencil problem, should be fine
         do jmp = 1, jmp_max


!Tomek:  this is just a dirty trick not to think too much about the ghost zones            
!!$  ip4 = min ( imp + 4 ,  nn+2 )
!!$  ip3 = min ( imp + 3 ,  nn+2 )
!!$  ip2 = min ( imp + 2 ,  nn+2 )
!!$  ip1 = min ( imp + 1 ,  nn+2 )
!!$  i0 = imp
!!$  im1 = max ( imp - 1 ,  mm-2 )
!!$  im2 = max ( imp - 2 ,  mm-2 )
!!$  im3 = max ( imp - 3 ,  mm-2 )
!!$  im4 = max ( imp - 4 ,  mm-2 )


  ip4 =  imp + 4 
  ip3 =  imp + 3 
  ip2 =  imp + 2 
  ip1 =  imp + 1 
  i0 =   imp
  im1 =  imp - 1 
  im2 =  imp - 2 
  im3 =  imp - 3 
  im4 =  imp - 4 
  
!Tomek


!Tomek: this is valid for  MP9            
!*        Original (unlimited) "left" interface value at i+1/2 (SH97; Eq. 2.31b)
            vorp(imp,jmp,imp5) = b1*( 4.d0 * u(im4,jmp,imp5) - 41.d0   * u(im3,jmp,imp5) + 199.d0  * u(im2,jmp,imp5)  &
                                   -641.d0 * u(im1,jmp,imp5) + 1879.d0 * u(i0 ,jmp,imp5) + 1375.d0 * u(ip1,jmp,imp5)  &
                                   -305.d0 * u(ip2,jmp,imp5) + 55.d0   * u(ip3,jmp,imp5) - 5.d0    * u(ip4,jmp,imp5) )

!*        Original (unlimited) "right" interface value at i-1/2 (mirror of Eq. 2.31b)
            vorm(imp,jmp,imp5) = b1*( 4.d0 * u(ip4,jmp,imp5) - 41.d0   * u(ip3,jmp,imp5) + 199.d0  * u(ip2,jmp,imp5)  &
                                   -641.d0 * u(ip1,jmp,imp5) + 1879.d0 * u(i0 ,jmp,imp5) + 1375.d0 * u(im1,jmp,imp5)  &
                                   -305.d0 * u(im2,jmp,imp5) + 55.d0   * u(im3,jmp,imp5) - 5.d0    * u(im4,jmp,imp5) )
!Tomek: this is valid for  MP9            

            
!*        Monotonicity preserving median value (SH97; Eq. 2.12)
         vmpp(imp,jmp,imp5) = u(imp,jmp,imp5) + dmm( du(imp+1,jmp,imp5),  alpha_mp*du(imp  ,jmp,imp5) )
         vmpm(imp,jmp,imp5) = u(imp,jmp,imp5) + dmm(-du(imp  ,jmp,imp5), -alpha_mp*du(imp+1,jmp,imp5) )

      end do
   end do

   do imp=mm-2,nn+2
      do jmp = 1, jmp_max

         up(imp,jmp,imp5) = vorp(imp,jmp,imp5)
         um(imp,jmp,imp5) = vorm(imp,jmp,imp5)

      end do
   end do

!*
!*     Compute U^L_{i+1/2}
!*

   do imp=mm-2,nn+2
      do jmp = 1, jmp_max

         if ((vorp(imp,jmp,imp5)-u(imp,jmp,imp5)) * (vorp(imp,jmp,imp5)-vmpp(imp,jmp,imp5)) .le. epsm ) then
            up(imp,jmp,imp5) = vorp(imp,jmp,imp5)
         else
!            djm1 = d2(imp-2,jmp,imp5)      ! v(j-2) - 2d0*v(j-1) + v(j)
!            dj   = d2(imp-1,jmp,imp5)      ! v(j-1) - 2d0*v(j) + v(j+1)
!            djp1 = d2(imp  ,jmp,imp5)      ! v(j) - 2d0 *v(j+1) + v(j+2)
            dm4jph = DM4(4d0*d2(imp,jmp,imp5)-d2(imp+1,jmp,imp5), 4d0*d2(imp+1,jmp,imp5)-d2(imp,jmp,imp5) &
                                                                , d2(imp,jmp,imp5), d2(imp+1,jmp,imp5))
            dm4jmh = DM4(4d0*d2(imp,jmp,imp5)-d2(imp-1,jmp,imp5), 4d0*d2(imp-1,jmp,imp5)-d2(imp,jmp,imp5) &
                                                                , d2(imp,jmp,imp5), d2(imp-1,jmp,imp5))
!*           [SH97; Eq. 2.8]
            vul = u(imp,jmp,imp5) + alpha_mp * du(imp,jmp,imp5)
!*           [SH97; Eq. 2.16]
            vav = .5d0 * ( u(imp,jmp,imp5) + u(imp+1,jmp,imp5) )
!*           [SH97; Eq. 2.28]
            vmd = vav - .5d0*dm4jph
!*           [SH97; Eq. 2.29]
            vlc = u(imp,jmp,imp5) + .5d0 * du(imp,jmp,imp5) + b2*dm4jmh
!*           [SH97; Eq. 2.24a]
            vmin = dmax1(dmin1(u(imp,jmp,imp5),u(imp+1,jmp,imp5),vmd), dmin1(u(imp,jmp,imp5),vul,vlc))
!*           [SH97; Eq. 2.24b]
            vmax = dmin1(dmax1(u(imp,jmp,imp5),u(imp+1,jmp,imp5),vmd), dmax1(u(imp,jmp,imp5),vul,vlc))
!*           [SH97; Eq. 2.26 (note that the median in defined as in Eq. (2.5): median(x,y,z) = x + minmod(y-x,z-x)]
            up(imp,jmp,imp5) = vorp(imp,jmp,imp5) + DMM(vmin-vorp(imp,jmp,imp5),vmax-vorp(imp,jmp,imp5))
!            up(imp,jmp,imp5) = vorp(imp,jmp,imp5) + DMM(vmpp(imp,jmp,imp5)-vorp(imp,jmp,imp5),u(imp,jmp,imp5)-vorp(imp,jmp,imp5))

         endif

      end do
   end do

!*
!*     Compute U^R_{i-1/2}
!*

   do imp=mm-2,nn+2
      do jmp = 1, jmp_max
            
         if ((vorm(imp,jmp,imp5)-u(imp,jmp,imp5)) * (vorm(imp,jmp,imp5)-vmpm(imp,jmp,imp5)) .le. epsm ) then
            um(imp,jmp,imp5) = vorm(imp,jmp,imp5)
         else
!            djm1 = d2(imp-2,jmp,imp5)      ! v(j-2) - 2d0*v(j-1) + v(j)
!            dj   = d2(imp-1,jmp,imp5)      ! v(j-1) - 2d0*v(j) + v(j+1)
!            djp1 = d2(imp  ,jmp,imp5)      ! v(j) - 2d0 *v(j+1) + v(j+2)
            dm4jph = DM4(4d0*d2(imp,jmp,imp5)-d2(imp+1,jmp,imp5), 4d0*d2(imp+1,jmp,imp5)-d2(imp,jmp,imp5) &
                                                                , d2(imp,jmp,imp5), d2(imp+1,jmp,imp5))
            dm4jmh = DM4(4d0*d2(imp,jmp,imp5)-d2(imp-1,jmp,imp5), 4d0*d2(imp-1,jmp,imp5)-d2(imp,jmp,imp5) &
                                                                , d2(imp,jmp,imp5), d2(imp-1,jmp,imp5))
            vul = u(imp,jmp,imp5) - alpha_mp * du(imp+1,jmp,imp5) !-> -alpha_mp * du(imp+1,jmp,imp5)
            vav = .5d0 * ( u(imp,jmp,imp5) + u(imp-1,jmp,imp5) )
            vmd = vav - .5d0*dm4jmh  !-> dm4jmh
            vlc = u(imp,jmp,imp5) - .5d0 * du(imp+1,jmp,imp5) + b2*dm4jph
            vmin = dmax1(dmin1(u(imp,jmp,imp5),u(imp-1,jmp,imp5),vmd), dmin1(u(imp,jmp,imp5),vul,vlc))
            vmax = dmin1(dmax1(u(imp,jmp,imp5),u(imp-1,jmp,imp5),vmd), dmax1(u(imp,jmp,imp5),vul,vlc))
            um(imp,jmp,imp5) = vorm(imp,jmp,imp5) + DMM(vmin-vorm(imp,jmp,imp5),vmax-vorm(imp,jmp,imp5))
!            um(imp,jmp,imp5) = vorm(imp,jmp,imp5) + DMM(vmpm(imp,jmp,imp5)-vorm(imp,jmp,imp5),u(imp,jmp,imp5)-vorm(imp,jmp,imp5))

         endif

      end do
   end do

!*
!*    if (varname .eq.'rho__') then
!*
!*     Limit the ocurrence of negative densities
!*
!   verify if  p_p(24) < 0 or rho_p(25) < 0  x-direction
   
         do imp=mm-2,nn+2
               
            if (up(imp,24,imp5).lt.0d0) then
!**               print*,'UP<0:',i,up(imp,24),u(imp,24)
               up(imp,24,imp5) = u(imp,24,imp5) + .50d0*DMM(du(imp,24,imp5),du(imp+1,24,imp5))
               um(imp,24,imp5) = u(imp,24,imp5)
!**               print*,'UP:::',i,up(imp,24),u(imp,24)            
            endif

            if (up(imp,25,imp5).lt.0d0) then
!**               print*,'UP<0:',i,up(imp,25),u(imp,25)
               up(imp,25,imp5) = u(imp,25,imp5) + .50d0*DMM(du(imp,25,imp5),du(imp+1,25,imp5))
               um(imp,25,imp5) = u(imp,25,imp5)
!**               print*,'UP:::',i,up(imp,25),u(imp,25)            
            endif
            
         end do

         
!   verify if  p_m(24) < 0 or rho_m(25) < 0  x-direction
         
      do imp=mm-2,nn+2
            
            if (um(imp,24,imp5).lt.0d0) then
!**               print*,'UM<0:',i,um(imp,24),u(imp,24)
               um(imp,24,imp5) = u(imp,24,imp5) - .50d0*DMM(du(imp,24,imp5),du(imp+1,24,imp5))
               up(imp,24,imp5) = u(imp,24,imp5)
!**               print*,'UM:::',i,um(imp,24),u(imp,24)
            endif

            if (um(imp,25,imp5).lt.0d0) then
!**               print*,'UM<0:',i,um(imp,25),u(imp,25)
!               up(imp,25) = u(imp,25)
!               um(imp,25) = u(imp,25)
               um(imp,25,imp5) = u(imp,25,imp5) - .50d0*DMM(du(imp,25,imp5),du(imp+1,25,imp5))
               up(imp,25,imp5) = u(imp,25,imp5)
!**               print*,'UM:::',i,um(imp,25),u(imp,25)
            endif
            
       end do

!*___________________________________________

      
       if (DIM ==2) then

          
!   verify if  p_p(61) < 0 or rho_p(62) < 0  y-direction

          do imp=mm-2,nn+2
               
            if (up(imp,61,imp5).lt.0d0) then
!**               print*,'UP<0:',i,up(imp,61,imp5),u(imp,61,imp5)
               up(imp,61,imp5) = u(imp,61,imp5) + .50d0*DMM(du(imp,61,imp5),du(imp+1,61,imp5))
               um(imp,61,imp5) = u(imp,61,imp5)
!**               print*,'UP:::',i,up(imp,61,imp5),u(imp,61,imp5)            
            endif

            if (up(imp,62,imp5).lt.0d0) then
!**               print*,'UP<0:',i,up(imp,62,imp5),u(imp,62,imp5)
               up(imp,62,imp5) = u(imp,62,imp5) + .50d0*DMM(du(imp,62,imp5),du(imp+1,62,imp5))
               um(imp,62,imp5) = u(imp,62,imp5)
!**               print*,'UP:::',i,up(imp,62,imp5),u(imp,62,imp5)            
            endif

         end do

!   verify if  p_m(61) < 0 or rho_m(62) < 0  y-direction
         
      do imp=mm-2,nn+2
            
            if (um(imp,61,imp5).lt.0d0) then
!**               print*,'UM<0:',i,um(imp,61,imp5),u(imp,61,imp5)
               um(imp,61,imp5) = u(imp,61,imp5) - .50d0*DMM(du(imp,61,imp5),du(imp+1,61,imp5))
               up(imp,61,imp5) = u(imp,61,imp5)
!**               print*,'UM:::',i,um(imp,61,imp5),u(imp,61,imp5)
            endif

            if (um(imp,62,imp5).lt.0d0) then
!**               print*,'UM<0:',i,um(imp,62,imp5),u(imp,62,imp5)
               um(imp,62,imp5) = u(imp,62,imp5) - .50d0*DMM(du(imp,62,imp5),du(imp+1,62,imp5))
               up(imp,62,imp5) = u(imp,62,imp5)
!**               print*,'UM:::',i,um(imp,62,imp5),u(imp,62,imp5)
            endif
            
       end do

    end if
!*___________________________________________         

!*
!*      endif

!*
!*     return
!*
     return
!*
!*     end of MP5_1D
!*

   end subroutine MP9
