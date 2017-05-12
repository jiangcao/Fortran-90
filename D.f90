module D
! Module to calculate the derivative numerically    
    use types, only : dp

    implicit none

    private

    real(dp), parameter :: abs_tol = 1e-5_dp     ! tolerance 
    real(dp), parameter :: rel_tol = 1e-3_dp     ! tolerance 
    integer, parameter  :: MaxIter = 100

    public :: D_1D
    public :: D_test_f

    contains


        ! Function returns the numerical derivative of an 1D function `f(x)` on a point `x`
        ! by using directly the definition of derivative. 
        ! Parameters:
        ! - `f` the 1D function
        ! - `x` the value around which to evaluate the derivative
        ! - `info` is -1 if the limit does not converge
        !
        ! @note
        ! This simple method may converge slowly, so many evaluations of `f(x)` may be needed, 
        ! which can be a bottleneck of the performance if the function is difficult to evaluate.
        function D_1D(f,x,info) result(df)
            implicit none
            interface
                real(dp) function func(x)
                    import
                    real(dp), intent(in) :: x 
                end function func
            end interface
            real(dp)                       :: Df
            real(dp),intent(in)            :: x
            integer,intent(out), optional :: info
            procedure(func)               :: f 
            real(dp) :: dx
            real(dp) :: error
            real(dp) :: dfxnew, dfx
            integer :: i                      !! counter of the number of iterations
            integer :: ns                     !! counter of the number of consecutive successes
            integer, parameter :: nsneed = 0  !! number of consecutive sucesses needed
            logical, parameter :: debug = .false.
            dx = 0.1_dp * x
            dfx   =  f(x+dx) / dx - f(x) / dx 
            df = dfx
            error = (Huge(1d0))
            i  = 0
            ns = 0
            if (debug) print *, "<<<< Begin D_1D >>>>"
            if (debug) print *, "x=", x, "f(x)=", f(x)
            do while (( i < MaxIter )  .and. (ns <= nsneed)  )
                dx = dx*0.5_dp 
                dfxnew =  ( f(x+dx) / dx ) - ( f(x)  / dx )
                if (dfx==0) then
                    error = abs(dfxnew - dfx)/abs(dfxnew)
                else
                    error = abs(dfxnew - dfx)/abs(dfx)
                endif
                dfx    = dfxnew
                i = i+1
                ns = ns +1
                if ( error > rel_tol ) then 
                    ns = 0 ! error larger than tolerance, reset the counter
                else
                    df = dfxnew
                endif
                if ( error == 1.0_dp) exit 
                if (debug) print *,ns,error,dfx,df
            enddo
            if ((ns < nsneed) .and. (present(info))) then
                info = -1
            endif
            if (ns < nsneed) then
                print *,"WARNING!!! D_1D NOT CONVERGE"
            else
                if (debug) print *, "D_1D, i=", i
            endif
            if (debug) print *, "<<<< End D_1D >>>>"
        end function D_1D 


        real(dp) function D_test_f(x)
            real(dp) , intent(in) :: x
            real(dp) :: alpha=0.5d0
            D_test_f = dsqrt(x * (1.0d0 + alpha*x))
        end function D_test_f
end module D
