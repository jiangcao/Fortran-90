module integrate
! Module to calculate numerically the integral of a function


    implicit none

    private

    public :: int_1D
    public :: int_test_fc

    contains

        ! Function returns the numerical integration of a function `f(x)` from x(1) to x(2)
        ! Parameters:
        ! - `f` the function to integrate
        ! - `x`   the integration interval
        ! - `by` the integration step
        real(8) function int_1D(f, x, by) 
            real(8), intent(in) :: x(2), by
            real(8) :: dx
            integer :: i, n
            interface
                real(8) function func(x)
                    real(8), intent(in) :: x 
                end function func
            end interface
            procedure(func)               :: f 
            int_1D = 0d0
            n = abs(x(2) - x(1))/by
            if (x(2) > x(1)) then
                dx = abs(by)
            else
                dx = -abs(by)
            endif
            do i = 1, n 
                !print *, i
                int_1D = int_1D + dx * f(x(1) + dx*dble(i))
            enddo
        end function int_1D

        real(8) function int_test_fc(x)
            real(8) , intent(in) :: x 
            int_test_fc = log(x) 
        endfunction int_test_fc

end module integrate
