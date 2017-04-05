module utils

    use types , only : dp
    
    implicit none

    private

    public seqi
    
    contains

    real(dp) function seqi(i,xmin,xmax,n)
        integer, intent(in)  :: n,i
        real(dp), intent(in) :: xmin,xmax
        seqi = xmin + (xmax - xmin ) / dble(n) * dble(i)
    end function seqi


end module utils
