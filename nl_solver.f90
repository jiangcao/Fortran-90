module nl_solver
! Module to solve the non-linear equations
    use types, only: dp
    use D, only : D_1D

    implicit none
    private

    real(dp), parameter :: abs_tol = 1e-10
    real(dp), parameter :: rel_tol = 1e-10

    public  nlsolver_bisec_1D, nlsolver_newton_1D
    public  nlsolver_linsearch_1D
    public  nlsolver_df_test, nlsolver_f_test


    contains
        

        !!!!!!!!!     1D      !!!!!!!!!


        ! 1D Newton-Raphson method to find the root of function `f(x)`
        ! Parameters:
        ! - `f` is the function
        ! - `df` is the derivative of `f`
        ! - `x0` is the initial guess value of x
        function nlsolver_newton_1D(f,df,x0) result(x)
            real(dp) :: x
            real(dp), intent(in) :: x0
            interface
                real(dp) function func(x)
                    import
                    real(dp), intent(in) :: x 
                end function func
            end interface
            procedure(func) :: f, df
            real(dp) :: err, xn
            integer :: i, ns
            integer, parameter :: NMAX = 1000
            integer, parameter :: NSN  = 1
            logical, parameter :: debug = .true.
            real(dp) :: dfx, fx, dx
            if (debug) print *, "<<<< Begin Newton  >>>>"
            err = HUGE(1d0)
            i = 0
            ns= 0
            xn = x0
            if (debug) open(unit = 22, file = "x_inNewton.dat")
            do while ((i<=NMAX) .and. (ns<=NSN))
                print *, xn
                fx = f(xn)
                dfx= df(xn)
                print *, "dfx=",dfx
                if (debug) write (22,*) xn
                if (dfx .ne. 0.0_dp) then 
                    dx = - fx/dfx
                else
                    dx = 1.0_dp
                endif
                print *, dx
                i = i+1
                xn = xn + dx
                err = abs(fx)
                if (debug) print *,i,err
                if (err<abs_tol) then
                    ns = ns + 1
                else
                    ns = 0
                endif
            enddo
            if (debug) close(22)
            if (err > abs_tol) then
                x = 0d0
                print *, "Newton NOT Converged "
                print *, "after ", i, " iterations, error= ", err
                call abort
            else
                x = xn
            endif
            if (debug) print *, "<<<< End Newton  >>>>"
        end function nlsolver_newton_1D


        logical function nlsolver_linsearch_1D(f,xmin,xmax, y, by,x)
            real(dp), intent(in)  :: xmin,xmax,by
            real(dp), intent(out) :: x(2)
            real(dp), intent(in)  :: y
            interface
                real(dp) function func(x)
                    import
                    real(dp), intent(in) :: x 
                end function func
            end interface
            procedure(func) :: f 
            real(dp) :: x1 
            logical, parameter :: debug=.false.
            x1=xmin
            nlsolver_linsearch_1D = .false.
            do while ((x1 < xmax ) .and. (.not. nlsolver_linsearch_1D))
                if ( ((f(x1)-y)*(f(x1+by)-y) <= 0d0) .or.(abs(f(x1)-y)<abs_tol) .or.(abs(f(x1+by)-y)<abs_tol) ) then 
                    nlsolver_linsearch_1D = .true.
                    x = (/x1, x1+by/)
                    if (debug) print *, x1,x1+by
                endif
                x1 = x1+by
            enddo
        end function nlsolver_linsearch_1D




        ! Bisection method to find the root of a non-linear 1D function f(x)=y in a defined interval
        ! Parameters:
        ! - `f` is the 1D function 
        ! - 'x(1)' and 'x(2)' are the bondaries of the bisection
        ! - `info` is an optional parameter for indicating the error 
        function nlsolver_bisec_1D(f,xrange,info, y,abstol) result(x0)
            real(dp), intent(in) :: xrange(2) 
            real(dp), intent(in) :: y
            real(dp), intent(in), optional :: abstol
            interface
                real(dp) function func(x)
                    import
                    real(dp), intent(in) :: x 
                end function func
            end interface
            procedure(func) :: f 
            integer, intent(out), optional :: info
            real(dp) :: x0,x(2)
            integer :: j,i
            logical, parameter :: debug= .true.
            logical :: fin
            integer, parameter  :: NMAX= 1000
            x = xrange
            j = 0
            fin = .false.
            do while ((j<nmax) .and. (.not. fin))
                do i = 1,2
                    if( present(abstol) ) then
                        if (abs(f(x(i))-y) < abstol) then 
                            x0 = x(i)
                            if (present(info)) info = 0
                            fin=.true.
                        endif
                    else
                        if (abs(f(x(i))-y) < abs_tol) then 
                            x0 = x(i)
                            if (present(info)) info = 0
                            fin=.true.
                        endif
                    endif
                enddo
                if ((.not. fin) .and. ((f(x(1))-y)*(f(x(2))-y) > 0d0)) then 
                    if (present(info)) info = -1
                    print *, "<<<< IN nlsolver_bisec_1D >>>>"
                    print *, "level = ",j 
                    print *, "x   = ", x 
                    print *, "f-y = ", f(x(1))-y,f(x(2))-y 
                    call abort
                else
                    if ((f(sum(x)/2d0)-y) * (f(x(1))-y) <= 0d0) then
                        x = (/ x(1), sum(x)/2d0 /)                   
                    else
                        x = (/ sum(x)/2d0, x(2) /)
                    endif
                endif       
                j = j+1
            enddo
            if ((.not. fin) .and. present(info)) info = -2
            if (.not. fin) then 
                print *, "<<<< IN nlsolver_bisec_1D >>>>"
                print *, "level = ",j 
                print *, "x   = ", x 
                print *, "f-y = ", f(x(1))-y,f(x(2))-y 
                call abort
            endif
        end function nlsolver_bisec_1D


real(dp) function nlsolver_f_test(x)
    real(dp), intent(in) :: x
    nlsolver_f_test = x**2-x - 1d0
end function nlsolver_f_test

real(dp) function nlsolver_df_test(x)
    real(dp), intent(in) :: x
    nlsolver_df_test = D_1D(nlsolver_f_test,x) 
end function nlsolver_df_test

end module nl_solver
