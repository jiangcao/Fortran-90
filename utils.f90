module utils

    use types , only : dp
    
    implicit none

    private

    public seq, seqi
    public cbind, rbind
    
    contains

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
