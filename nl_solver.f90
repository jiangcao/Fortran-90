module nl_solver
! Module to solve the non-linear equations
    use progressbar
    use types, only: rfn, rndfn, dp
    use D, only : D_1D
    use mpi, only : idnode

    implicit none
    private

    real(dp), parameter :: abs_tol = 1e-10
    real(dp), parameter :: rel_tol = 1e-10

    public  nlsolver_bisec_1D,nlsolver_bisec_1Dg, nlsolver_newton_1D
    public  nlsolver_linsearch_1D
    public  nlsolver_linsearch_1Dg


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

        logical function nlsolver_linsearch_1Dg(f,x0,ivar,xmin,xmax,y,by,x)
            procedure(rndfn) :: f
            real(dp), intent(in) :: x0(:)
            integer, intent(in) :: ivar
            real(dp), intent(in)  :: xmin,by
            real(dp), intent(in), optional :: xmax
            real(dp), intent(out) :: x(2)
            real(dp), intent(in)  :: y

            nlsolver_linsearch_1Dg = nlsolver_linsearch_1D(f1, xmin, xmax, y, by, x)
            
        contains
            
            real(dp) function f1(x)
                real(dp),intent(in) :: x
                real(dp) :: xp(size(x0))
                xp = x0
                xp(ivar) = x
                f1 = f(xp)
            end function f1
        end function nlsolver_linsearch_1Dg


        logical function nlsolver_linsearch_1D(f,xmin,xmax, y, by,x)
            real(dp), intent(in)  :: xmin,by
            real(dp), intent(in), optional :: xmax
            real(dp), intent(out) :: x(2)
            real(dp), intent(in)  :: y
            interface
                real(dp) function func(x)
                    import
                    real(dp), intent(in) :: x 
                end function func
            end interface
            procedure(func) :: f 
            real(dp) :: x1 , nlxmax
            logical, parameter :: debug=.false.
            x1=xmin
            if (present(xmax)) then 
                nlxmax = xmax
            else
                nlxmax = HUGE(1.0_dp)
            endif
            nlsolver_linsearch_1D = .false.
            do while ((x1 < nlxmax ) .and. (.not. nlsolver_linsearch_1D))
                if ( ((f(x1)-y)*(f(x1+by)-y) <= 0d0) .or.(abs(f(x1)-y)<abs_tol) .or.(abs(f(x1+by)-y)<abs_tol) ) then 
                    nlsolver_linsearch_1D = .true.
                    x = (/x1, x1+by/)
                    if (debug) print *, x1,x1+by
                endif
                x1 = x1+by
            enddo
        end function nlsolver_linsearch_1D


        function nlsolver_bisec_1Dg(f, x, ivar, xrange, y, abstol, bar,info) result(x0)
            real(dp), intent(in) :: xrange(2) 
            real(dp), intent(in) :: x(:)
            integer, intent(in) :: ivar
            integer, intent(out), optional :: info
            real(dp), intent(in) :: y
            real(dp), intent(in), optional :: abstol
            procedure(rndfn) :: f 
            logical, intent(in), optional  :: bar
            real(dp) :: x0

            x0 = nlsolver_bisec_1D(f1,xrange=xrange, y=y,abstol=abstol,bar=bar, info=info)

        contains
            real(dp) function f1(xi)
                real(dp),intent(in) :: xi
                real(dp) :: xp(size(x))
                xp = x
                xp(ivar) = xi
                f1 = f(xp)
            end function f1
        end function nlsolver_bisec_1Dg


        ! Bisection method to find the root of a non-linear 1D function f(x)=y in a defined interval
        ! Parameters:
        ! - `f` is the 1D function 
        ! - 'x(1)' and 'x(2)' are the bondaries of the bisection
        ! - `info` is an optional parameter for indicating the error 
        function nlsolver_bisec_1D(f,xrange,info, y,abstol,bar) result(x0)
            real(dp), intent(in) :: xrange(2) 
            real(dp), intent(in) :: y
            real(dp), intent(in), optional :: abstol
            procedure(rfn) :: f 
            integer, intent(out), optional :: info
            logical, intent(in), optional  :: bar
            real(dp) :: x0,x(2)
            real(dp) :: fm,fx(2)
            real(dp) :: labstol
            integer :: j,i
            logical, parameter :: debug= .true.
            logical :: fin
            integer, parameter  :: NMAX= 1000
            x = xrange
            j = 0
            fin = .false.
            labstol = merge(abstol, abs_tol, present(abstol))
            if (present(bar) .and. bar) print *,"<<<< bisection >>>>" 
            if (present(bar) .and. bar) call set_progress_bar()
            fx(1) = f(x(1))
            fx(2) = f(x(2))
            do while ((j<nmax) .and. (.not. fin))
                do i = 1,2
                    if (abs(fx(i)-y) < labstol) then 
                        x0 = x(i)
                        if (present(info)) info = 0
                        fin=.true.
                    endif
                enddo
                if (.not. fin) then
                    if ((fx(1)-y)*(fx(2)-y) > 0d0) then 
                        if (present(info)) info = -1
                        print *, "<<<< IN nlsolver_bisec_1D: xrange too small or too big >>>>"
                        print *, "level = ",j 
                        print *, "x   = ", x 
                        print *, "f-y = ", fx(1)-y,fx(2)-y 
                        print *, idnode
                        call abort
                    else
                        fm = f(sum(x)/2d0)
                        if ((fm-y) * (fx(1)-y) <= 0d0) then
                            x = (/ x(1), sum(x)/2d0 /)                   
                            fx(2) = fm
                        else
                            fx(1) = fm
                            x = (/ sum(x)/2d0, x(2) /)
                        endif
                    endif       
                    j = j+1
                endif
                if (present(bar) .and. bar) call progress_bar(min(1.0_dp,labstol/minval(abs(fx-y))))
            enddo
            if ((.not. fin) .and. present(info)) info = -2
            if (.not. fin) then 
                print *, "<<<< IN nlsolver_bisec_1D: still not fully converge >>>>"
                print *, "x   = ", x 
                print *, "f-y = ", fx(1)-y,fx(2)-y 
                x0 = sum(x)/2.0_dp
            endif
        end function nlsolver_bisec_1D

    end module nl_solver
