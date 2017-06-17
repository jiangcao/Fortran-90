module integrate_poly
	! module contains different numerical integration schemes based on quadrature rule.

use types, only : dp, rfn
use constants, only : pi

implicit none

private
   real(dp), PARAMETER :: EPS=3.0d-15       	!EPS is the relative precision

public Trapezoid_multi, qgss2d, Trapezoid_fc, trapezoid

contains

pure function Trapezoid(x, y) result(r)
  !! Calculates the integral of an array y with respect to x using the trapezoid
  !! approximation. Note that the mesh spacing of x does not have to be uniform.
  real(dp), intent(in)  :: x(:)         !! Variable x
  real(dp), intent(in)  :: y(size(x))   !! Function y(x)
  real(dp)              :: r            !! Integral ∫y(x)·dx
  ! Integrate using the trapezoidal rule
  associate(n => size(x))
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
  end associate
end function

function Trapezoid_fc(x, fx) result(r)
  !! Calculates the integral of an array y with respect to x using the trapezoid
  !! approximation. Note that the mesh spacing of x does not have to be uniform.
  real(dp), intent(in)  :: x(:)         !! Variable x
  procedure(rfn)        :: fx
  real(dp)              :: y(size(x))   !! Function values y(x)
  real(dp)              :: r            !! Integral ∫y(x)·dx
  integer :: i
  ! Integrate using the trapezoidal rule
  associate(n => size(x))
  	do i = 1,n 
  		y(i) = fx(x(i)) 
  	enddo
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
  end associate
end function


pure function Trapezoid_multi(x, y, ynd) result(r)
  !! Calculates the integral of an array y with respect to x using the trapezoid
  !! approximation. Note that the mesh spacing of x does not have to be uniform.
  real(dp), intent(in)  :: x(:)             !! Variable x
  real(dp), intent(in)  :: y(size(x),ynd)   !! Function y(x)
  integer,  intent(in)  :: ynd              !! Dimension of function y
  real(dp)              :: r(ynd)           !! Integral ∫y(x)·dx
  integer :: i
  ! Integrate using the trapezoidal rule
  associate(n => size(x))
  	do i = 1, ynd 
    	r(i) = sum((y(1+1:n-0,i) + y(1+0:n-1,i))*(x(1+1:n-0) - x(1+0:n-1)))/2
    enddo
  end associate
end function

!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************
SUBROUTINE  gauleg(ngp, xabsc, weig)

    implicit none
    INTEGER  i, j, m
    REAL(dp)  p1, p2, p3, pp, z, z1
    INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
    REAL(dp), INTENT(OUT) :: xabsc(ngp), weig(ngp)


    m = (ngp + 1) / 2
    !* Roots are symmetric in the interval - so only need to find half of them  */

    do i = 1, m				! Loop over the desired roots */

        z = cos( pi * (i-0.25d0) / (ngp+0.5d0) )
        !*   Starting with the above approximation to the ith root,
        !*          we enter the main loop of refinement by NEWTON'S method   */
        100     	p1 = 1.0d0
        p2 = 0.0d0
        !*  Loop up the recurrence relation to get the Legendre
        !*  polynomial evaluated at z                 */

        do j = 1, ngp
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        enddo

        !* p1 is now the desired Legendre polynomial. We next compute pp,
        !* its derivative, by a standard relation involving also p2, the
        !* polynomial of one lower order.      */
        pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        z1 = z
        z = z1 - p1/pp             ! Newton's Method  */

        if (dabs(z-z1) .gt. EPS) GOTO  100

        xabsc(i) =  - z                         ! Roots will be bewteen -1.0 & 1.0 */
        xabsc(ngp+1-i) =  + z                   ! and symmetric about the origin  */
        weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp)     ! Compute the weight and its       */
        weig(ngp+1-i) = weig(i)                 ! symmetric counterpart         */

    end do     ! i loop

End subroutine gauleg



!********************************************************************************
!*     Returns the SINGLE integral of the function (of ONE VARIABLE) "func"
!* between x1 and x2 by N-point Gauss-Legendre integration. The function
!* is evaluated exactly N times at interior points in the range of
!* integration.       */
!********************************************************************************
   recursive function qgauss(func, x1, x2, ngp) RESULT(intgrl)
      implicit none
    REAL(dp)  intgrl, x1, x2, func
      REAL(dp)  xm, xl
      INTEGER j
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dp) :: xabsc(ngp), weig(ngp)

      call gauleg(ngp, xabsc, weig)

      intgrl = 0.0d0
     xm = 0.5 * (x2 + x1)
    xl = 0.5 * (x2 - x1)
      do j = 1, ngp
         intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j) )
      END do

    intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   END function qgauss


   recursive function qgaussmat1(func, x1, x2, ngp, frow, fcol) RESULT(intgrl)
      implicit none
      REAL(dp), INTENT(IN) :: x1, x2
      INTEGER :: frow, fcol
    REAL(dp) :: intgrl(frow, fcol), tmpm(frow,fcol)
      REAL(dp) ::  func
      REAL(dp) ::  xm, xl, arg
      INTEGER j

      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dp) :: xabsc(ngp), weig(ngp)

      call gauleg(ngp, xabsc, weig)

      intgrl(:,:) = 0.0d0
      tmpm(:,:) = 0.0d0
     xm = 0.5 * (x2 + x1)
    xl = 0.5 * (x2 - x1)
      do j = 1, ngp
         arg =  xm + xl*xabsc(j)
      PRINT *, 'szhgd ds'
          PRINT *,arg
         tmpm = func(arg)
         intgrl = intgrl + weig(j) * tmpm
      END do

    intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   END function qgaussmat1





!**********************************************************************************
!*           Generic 2D Gaussian Quadraure Routines                               *
!**********************************************************************************
!* Use this function to calculate 2D integral by Gaussian Quadrature
!* Must supply boundary conditions x1,x2, y1(x), y2(x)
!* and the function to be integrated over f(x,y)
   RECURSIVE function qgss2d(origfn, xx1, xx2, yf1, yf2, ngp) RESULT(inth)
       implicit none                               ! returns integral in inth
       REAL(dp) :: inth, xx1, xx2
       REAL(dp) :: xm, xl, xtmp
       INTEGER :: j
       INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
       REAL(dp) :: xabsc(ngp), weig(ngp)

       interface
           function origfn(xp,yp) RESULT(vfun2d)     ! Original Function's Interface
               import
               implicit none
               REAL(dp) ::  vfun2d, xp, yp
           end function origfn
           function yf1(x) RESULT(vy1)
               import
               implicit none
               REAL(dp) :: vy1, x
           end  function yf1
           function yf2(x) RESULT(vy2)
               import
               implicit none
               REAL(dp) :: vy2, x
           end  function yf2
       end interface

       call gauleg(ngp, xabsc, weig)

       inth = 0.0d0
       xm = 0.5 * (xx2 + xx1)
       xl = 0.5 * (xx2 - xx1)

       do j = 1, ngp
           xtmp = xm + xl*xabsc(j)                ! Gauss-Legendre Abcissas
           inth = inth + weig(j) * qgssgy()
       END do

       inth = inth * xl;    !Scale the answer to the range of integration  */

   CONTAINS
       RECURSIVE function qgssgy() RESULT(intg)
           implicit none                                ! returns integral in intg
           REAL(dp) :: intg
           REAL(dp) :: ym, yl, ytmp                ! all undeclared variables are
           INTEGER j                                    !   COOMON with HOST

           intg = 0.0d0
           ym = 0.5 * (yf2(xtmp) + yf1(xtmp))
           yl = 0.5 * (yf2(xtmp) - yf1(xtmp))

           do j = 1, ngp
               ytmp = ym + yl*xabsc(j)                ! Gauss-Legendre Abcissas
               intg = intg + weig(j) * origfn(xtmp,ytmp)
           END do           
           intg = intg * yl;    !Scale the answer to the range of integration  */
       END function qgssgy

   END FUNCTION  qgss2d



end module integrate_poly
