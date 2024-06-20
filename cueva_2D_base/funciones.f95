  module funciones

  CONTAINS

  !     *******************************************************************
  !     Error Function
  !     *******************************************************************

    doubleprecision function errorf(x)    
!
!     Compute the error funtion using the Mclaurin series of it
!
      implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     include "constants.h"
      integer                    :: r                           !nmax
      doubleprecision            :: x, sum, n_factor, diff, tmp,aerr
      doubleprecision,parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                

      if (abs(x).lt.3.75d0) then 
         sum = x
         n_factor = 1.d0
         diff = 1.d0
         r = 0
         do while (diff.gt.1.d-7)
            r = r + 1
            n_factor = n_factor * r
            tmp = (-1)**r * x**(2*r+1) / ( (2*r + 1.d0) * n_factor )
            diff = abs( tmp / sum )
            sum = sum + tmp
         end do
         errorf = 2.d0 / sqrt(pi) * sum
      else
!
!     for x>~3.75, the Maclaurin expansion converges very slowly, but for large 
!     values of x, there is a good asymptotic expansion 
!     (see Wikipedia: http://es.wikipedia.org/wiki/Función_error )
!
         aerr = - 8.d0 * ( pi - 3.d0 ) / ( 3.d0 * pi *( pi - 4.d0 ) )
         errorf = sqrt( 1.d0 - exp(-x**2 * (4.d0/pi + aerr*x**2) / (1.d0+aerr*x**2))) * sign(1.d0, x)
         
      endif

      return

    end function errorf




  !     *******************************************************************
  !     Flux reconstruction
  !     *******************************************************************
          
      doubleprecision function  mcl(fa, fb)  

        use scalar

        implicit none
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      doubleprecision            :: fa, fb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (slop == 1) then 

         ! MinMod slop Limiter

         mcl = 0.5d0 * ( sign(1.d0,fa) + sign(1.d0,fb) ) * min(abs(fa),abs(fb))

      else if (slop == 2) then

       !  Monotoniced Central Limiter MCL    

         mcl = 0.5d0 * ( sign(1.d0,fa) + sign(1.d0,fb) ) * min(2.d0*abs(fa),2.d0*abs(fb),0.5d0*abs(fa+fb))

      else if (slop == 3) then

       ! Superbee Limiter

         mcl = 0.5d0 * ( sign(1.d0,fa) + sign(1.d0,fb) ) * min(2.d0*abs(fa),2.d0*abs(fb),max(abs(fa),abs(fb)))

      else

         write(*,*) "This parameter for slop limiter is not correctly"
         write(*,*) "slop =", slop
         stop

      end if

      return

    end function mcl




  !     *******************************************************************
  !     Primitive variables reconstruction
  !     *******************************************************************
          
      doubleprecision function  limiter(ga, gb)  

        use scalar

      implicit none
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      doubleprecision            :: ga, gb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      if (limit == 1) then 

! MinMod slop Limiter

         limiter = 0.5d0*( sign(1.d0,ga) + sign(1.d0,gb) ) * min(abs(ga),abs(gb))

      else if (limit == 2) then

!  Monotoniced Central Limiter MCL    

         limiter = 0.5d0 * ( sign(1.d0,ga) + sign(1.d0,gb) ) * min(2.d0*abs(ga),2.d0*abs(gb),0.5d0*abs(ga+gb))

      else if (limit == 3) then

! Superbee Limiter

         limiter = 0.5d0 * ( sign(1.d0,ga) + sign(1.d0,gb) ) * min(2.d0*abs(ga),2.d0*abs(gb),max(abs(ga),abs(gb)))

      else

         write(*,*) "This parameter for slop limiter is not correctly"
         write(*,*) "slop =", slop
         stop
      
      end if

      return

    end function limiter

  !     *******************************************************************
  !     MP9 and MP5 functions
  !     *******************************************************************

      DOUBLEPRECISION FUNCTION DMM(X,Y)
!*
      implicit none
!*      
      doubleprecision x, y
!*
!*     Minmod function for 2 arguments [SH97; Eq. 2.4]
!*
      DMM = .5d0 * ( dsign(1d0,x) + dsign(1d0,y) ) * dmin1( dabs(x), dabs(y) )
!*
!*     return
!*
      return
!*
      end FUNCTION DMM


     DOUBLEPRECISION FUNCTION DM4(W,XX,Y,Z)
!*
      implicit none
!*      
      doubleprecision w, xx, y, z      
!*
!*     Minmod function for 4 arguments [SH97; Eqs. 2.6a - 2.6b] 
!*
      DM4 = .125d0 * ( dsign(1d0,w) + dsign(1d0,xx) )     &
             * dabs( ( dsign(1d0,w) + dsign(1d0,y ) )   * &
                     ( dsign(1d0,w) + dsign(1d0,z ) ) ) * &
            dmin1(dabs(w), dabs(xx), dabs(y), dabs(z) )
!*
!*
!*     return
!*
      return
!*
!*     end of DM4
!*
      end FUNCTION DM4


      doubleprecision function get_hllc_vel(a,b,c,rt)

        
!!$/*@@
!!$   @function  get_hllc_vel
!!$   @date      Tuesday 18 July 2006
!!$   @author    Miguel Angel Aloy Toras
!!$   @desc
!!$   This function obtains the root of the quadratic equation
!!$
!!$       a*v^2 - b*v + c = 0
!!$
!!$   in the physical range [-1,1] usign a Newton-Rapson iterative method.
!!$  @enddesc
!!$@@*/

        
      implicit none

!     Coeffcients of the quadratic equation a*v^2 - b*v + c = 0
      DOUBLEPRECISION a, b, c, rt

      DOUBLEPRECISION v

      INTEGER, PARAMETER :: maxiter = 100
      INTEGER               iter

      DOUBLEPRECISION, PARAMETER:: minerr = 1.d-12
      DOUBLEPRECISION, PARAMETER:: diferr = 1.d-3
      DOUBLEPRECISION fv, dfdv, oldv, err
     
!*
!*     Initial seed (as if a=0) but limited to the range [-1,1]
!*
      oldv = max(-1.d0, min( 1.d0, rt ))
!c      print*,'Initial seed', oldv, a,b,c
      iter = 0
!*
!*     Start Newton-Rapson loop
!*

333   fv = a*v*v - b*v + c
      dfdv = 2.d0*a*v - b

      v = max(-1.d0, min( 1.d0, oldv - fv/dfdv ))
      err = abs( oldv - v )
!c      print'(i3,1x,5e12.5)',iter,oldv,v,fv,dfdv,err
      oldv = v
      iter = iter+1
!*
!*     Maximum number of iterations reached?
!*
      if (iter.ge.maxiter) then
         print'(a,i3,1x,a)',                                            &
              '[GET_HLLC_VEL]: Non convergence in velocity after',iter, &
              'iterations'
         stop
      endif
!*
!*     Error larger than tolerance?
!*
      if (err.gt.minerr) goto 333

      get_hllc_vel = max(-1.d0 + diferr, min( 1.d0 - diferr, v )) ! v
!*
!*
!*
      return
!*
!*     end of get_hllc_vel
!*
    end function get_hllc_vel   

    logical function isnan(a)

      doubleprecision :: a

      if (a.ne.a) then

         isnan = .true.

      else

         isnan = .false.

      end if

      return

    end function isnan


  end module funciones
