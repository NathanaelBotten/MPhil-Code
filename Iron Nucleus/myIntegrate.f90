MODULE myIntegrate

   USE nrtype
   USE quadpack
   USE adapt_quad
   USE DCUHRE90

   implicit none

!
!  (PREVIOUSLY) Local variables
!
   real(DP)             :: zresult, errEstimate
   integer              :: neval,  ifail, last
!
   real(DP), parameter :: Infty = HUGE(1.0_DP)

!  Overload all of the integration routines below into one nice simple to use function: integral
!    Each function takes a different number of arguments so the compiler can tell which to use

   INTERFACE integral
      MODULE PROCEDURE integralND, integral, integralqxgs, integralToInfty, integralOf, integralBreakPts, integralSin
   END INTERFACE


CONTAINS


  function integralND(func, a, b, absErr, relErr)
    implicit none
    real(DP)             :: integralND
    INTERFACE  
       FUNCTION func(x)
         USE nrtype
         IMPLICIT NONE
         REAL(dp), intent(in) :: x(:)
         REAL(dp)             :: func
       END FUNCTION func
    END INTERFACE
    real(DP), intent(in) :: a(:), b(:), absErr, relErr
    !
    !  Local variables
    !
    real(DP)                :: result(1), errEstimateArray(1)
    integer                 :: ndim
    integer                 :: mincls = 0, maxcls = 5000000, key = 0, restart = 0
    integer, parameter      :: nw = 500000
    real(dp), dimension(nw) :: work
    !
    !  Get integral dimension from the limits 
    !
    ndim = SIZE(a)
    if (ndim < 2 .or. ndim > 15) then  ! Routine can only handle limited dimensions
       write(*,*) ' Error in integralND: invalid ndim = ', ndim
    else
!!$       !JDC
!!$       print *,"a limits = ",a
!!$       print *,"b limits = ",b
!!$       !JDC
       CALL DCUHRE(ndim, a, b, mincls, maxcls, func, abserr, relerr,        &
            &      key, nw, restart, result, errEstimateArray, neval, ifail, work)
       integralND  = result(1)
       errEstimate = errEstimateArray(1)
       if ( ifail /= 0 ) then
          write(*,*) ' Warning from DCUHRE: the error code is ', ifail
       end if
    end if
    return
  end function integralND




   function integral(f, a, b, absErr, relErr)
      implicit none
      real(DP)             :: integral
      INTERFACE
         FUNCTION f(x)
            USE nrtype
            IMPLICIT NONE
            REAL(DP)             :: f
            REAL(DP), intent(in) :: x
         END FUNCTION f
      END INTERFACE
      real(DP), intent(in) :: a, b, absErr, relErr
      real(DP)             :: bound=0.0_DP

!!$!qxgs:
!!$      integer              :: last
!!$      integer    :: limit=10000
      integer              :: inf

!
!  Determine if the limits include infinity and call qagi if nessary
!
      if (a == -Infty) then
         if (b == Infty) then
            inf=2
      CALL qagi(f, bound, inf, absErr, relErr, zresult, errEstimate, neval, ifail)
            if ( ifail /= 0 ) then
               write(*,*) ' Warning from qagi: the error code is ', ifail
            end if
         else
            inf = -1
            bound = b
      CALL qagi(f, bound, inf, absErr, relErr, zresult, errEstimate, neval, ifail)
            if ( ifail /= 0 ) then
               write(*,*) ' Warning from qagi: the error code is ', ifail
            end if
         end if
      else 
         if (b == Infty) then
            inf = 1
            bound = a
      CALL qagi(f, bound, inf, absErr, relErr, zresult, errEstimate, neval, ifail)
            if ( ifail /= 0 ) then
               write(*,*) ' Warning from qagi: the error code is ', ifail
            end if
         else
      CALL qags(f, a, b, absErr, relErr, zresult, errEstimate, neval, ifail)
            if ( ifail /= 0 ) then
!outputerrors               write(*,*) ' Warning from qags (myIntegrate): the error code is ', ifail
!               STOP  !! MADEACHANGEBIG!
            end if
!!$            if ( errEstimate .gt. 5.0e-12_dp ) then
!               write(*,'(A,F25.20)')" error estimate from qags is ",errEstimate
!!$            end if
         end if
      end if
      integral = zresult
      return
   end function integral

   function integralqxgs(f, a, b, absErr, relErr, limit, failcode, totalerr)
      implicit none
      real(DP)             :: integralqxgs
      INTERFACE
         FUNCTION f(x)
            USE nrtype
            IMPLICIT NONE
            REAL(DP)             :: f
            REAL(DP), intent(in) :: x
         END FUNCTION f
      END INTERFACE
      real(DP), intent(in) :: a, b, absErr, relErr
      real(DP)             :: bound=0.0_DP

!qxgs:
!      integer              :: last
      integer              :: inf
      integer              :: limit
      integer              :: failcode
      real(dp) :: totalerr
      CALL qxgs(f, a, b, absErr, relErr, zresult, errEstimate, ifail, limit, last)
!!$            if ( ifail /= 0 ) then
!!$               write(*,'(a,i2,a,i4,a)') ' Warning from qxgs (myIntegrate): the error code is ', ifail,". ",last," subintervals used"
!!$            end if
!            print *,"errEst = ",errEstimate
!            if ( errEstimate .gt. 0.0005 ) then
!               write(*,'(A,F20.16)')" error estimate from qxgs is ",errEstimate
 !           end if
      integralqxgs = zresult
      failcode = ifail
      totalerr = errEstimate
      return
   end function integralqxgs

   function integralToInfty(f, bound, absErr, relErr)
      implicit none
      real(DP)             :: integralToInfty
      INTERFACE
         FUNCTION f(x)
            USE nrtype
            IMPLICIT NONE
            REAL(DP)             :: f
            REAL(DP), intent(in) :: x
         END FUNCTION f
      END INTERFACE
      real(DP), intent(in) :: bound, absErr, relErr

      integer, parameter   :: inf=1

      CALL qagi(f, bound, inf, absErr, relErr, zresult, errEstimate, neval, ifail)
      integralToInfty = zresult
      if ( ifail /= 0 ) then
         write(*,*) ' Warning from qagi: the error code is ', ifail
      end if
      return
   end function integralToInfty

   function integralOf(f, absErr, relErr)
      implicit none
      real(DP)             :: integralOf
      INTERFACE
         FUNCTION f(x)
            IMPLICIT NONE
            INTEGER,  parameter  :: DP = KIND(1.0d0)
            REAL(DP)             :: f
            REAL(DP), intent(in) :: x
         END FUNCTION f
      END INTERFACE
      real(DP), intent(in) :: absErr, relErr
      real(DP)             :: bound=0.0d0
      integer, parameter   :: inf=2

      CALL qagi(f, bound, inf, absErr, relErr, zresult, errEstimate, neval, ifail)
      integralOf = zresult
      if ( ifail /= 0 ) then
         write(*,*) ' Warning from qagi: the error code is ', ifail
      end if
      return
   end function integralOf

   function integralBreakPts(f, a, b, absErr, relErr, nBreakPts, BreakPts)
      implicit none
      real(DP)             :: integralBreakPts
      INTERFACE
         FUNCTION f(x)
            USE nrtype
            IMPLICIT NONE
            REAL(DP)             :: f
            REAL(DP), intent(in) :: x
         END FUNCTION f
      END INTERFACE
      real(DP), intent(in)       :: a, b, absErr, relErr
      integer,  intent(in)       :: nBreakPts
      real(DP), intent(in), dimension(nBreakPts) :: BreakPts

      real(DP), dimension(nBreakPts+2) :: BreakPtsP2

      BreakPtsP2(1:nBreakPts) = BreakPts(1:nBreakPts)
      CALL qagp(f, a, b, nBreakPts+2, BreakPtsP2, absErr, relErr, zresult, errEstimate, neval, ifail)
      integralBreakPts = zresult
      if ( ifail /= 0 ) then
         write(*,*) ' Warning from qagp: the error code is ', ifail
      end if
      return
   end function integralBreakPts

   function integralSin(f, a, b, omega, absErr, relErr)

     implicit none
     real(dp)  :: IntegralSin
     INTERFACE 
        FUNCTION f(x)
          use nrtype
          implicit none
          real(DP) :: f
          real(DP), intent(in) :: x
        end FUNCTION f
     end INTERFACE
     
     real(dp), intent(in) :: a, b,absErr, relErr, omega
     integer :: sin = 2
     
     call qawo(f, a, b, omega, sin, absErr, relErr, zresult, errEstimate, neval, ifail)
     integralSin = zresult
     if (ifail /= 0 ) then
        write(*,*) ' Warning from qawo: the error code is ', ifail
     end if
     return
   end function integralSin
   
end MODULE myIntegrate



