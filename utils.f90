module utils

    use types , only : dp
    
    implicit none

    private

    public seqi, seq, polar2cart
    public cbind, rbind
    public cbind_3c
    public swap
    
    contains

    subroutine swap(a,b)
    implicit none
        real(dp), intent(inout) :: a,b
        real(dp) :: c
        c =a
        a =b
        b =c
    end subroutine swap


    real(dp) function seqi(i,xmin,xmax,n)
        integer, intent(in)  :: n,i
        real(dp), intent(in) :: xmin,xmax
        seqi = xmin + (xmax - xmin ) / dble(n) * dble(i)
    end function seqi

    function seq(xmin,xmax,n) result(s)
        integer, intent(in)  :: n
        real(dp), intent(in) :: xmin,xmax
        real(dp) :: s(n) 
        integer  :: i
        do i =1 , n
            s(i) = xmin + (xmax - xmin ) / dble(n) * dble(i)
        enddo
    end function seq

    ! Function returns the Cartesian coordinates from the Polar coordinates
    ! polar(theta, phi, r)
    function polar2cart(polar) result(cart)
        real(dp), intent(in) :: polar(3)    
        real(dp) :: cart(3)
        real(dp) :: theta, phi, r
        theta = polar(1)
        phi   = polar(2)
        r     = polar(3)
        cart(1) = r * sin(theta) * cos(phi)
        cart(2) = r * sin(theta) * sin(phi)
        cart(3) = r * cos(theta)
    end function polar2cart

    function cbind_3c(x,y,z) result(m)
        real(dp), intent(in) :: x(:),y(:),z(:)
        real(dp) :: m(size(x),3)
        m(:,1) = x(:)
        m(:,2) = y(:)
        m(:,3) = z(:)
    end function cbind_3c


    function cbind(x,y) result(m)
        real(dp), intent(in) :: x(:),y(:)
        real(dp) :: m(size(x),2)
        m(:,1) = x(:)
        m(:,2) = y(:)
    end function cbind

    function rbind(x,y) result(m)
        real(dp), intent(in) :: x(:),y(:)
        real(dp) :: m(2,size(x))
        m(1,:) = x(:)
        m(2,:) = y(:)
    end function rbind



end module utils
