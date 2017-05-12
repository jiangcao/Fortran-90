module integrate
! Module to calculate numerically the integral of a function

    use types, only : dp, rfn, r2fn , r3fn
    use progressbar

    implicit none

    private

    public :: int_1D
    public :: int_2D
    public :: int_3D

    contains

        ! Function returns the numerical integration of a function `f(x)` from x(1) to x(2)
        ! Parameters:
        ! - `f` the function to integrate
        ! - `x`   the integration interval
        ! - `by` the integration step
        real(dp) function int_1D(f, x, by, saveto,bar) 
            real(dp), intent(in) :: x(2), by
            real(dp) :: dx, fx, xi
            logical, intent(in), optional :: bar
            logical :: lbar
            integer :: i,np,npdone, n
            procedure(rfn)               :: f 
            character(len=*), intent(in), optional :: saveto
            integer :: u
            lbar = merge( bar, .false., present(bar))
            if (lbar)     print *, " <<< 1D integral >>>"
            if (lbar)     call set_progress_bar() 
            if (present(saveto)) then 
                open(newunit = u, file = trim(saveto))
            endif
            int_1D = 0d0
            n = abs(x(2) - x(1))/by
            if (x(2) > x(1)) then
                dx = abs(by)
            else
                dx = -abs(by)
            endif
            do i = 1, n 
                xi = x(1) + dx*dble(i)
                fx = f( xi )
                int_1D = int_1D + dx * fx 
                if( present(saveto) ) then 
                    write(u, *) xi, fx 
                endif
                if (lbar)       call progress_bar( dble(i)/dble(n) )
            enddo
            if (present(saveto)) close(u)
        end function int_1D

        real(dp) function int_2D(f, x, y, xby, yby, bar)
            real(dp), intent(in) :: x(2), y(2), xby, yby
            procedure(r2fn) :: f
            real(dp) :: tmp_x
            integer :: np,npdone
            logical, intent(in), optional :: bar
            logical :: lbar
            lbar = merge( bar, .false., present(bar))
            npdone=0
            if (lbar)    print *, " <<< 2D integral >>>"
            if (lbar)    call set_progress_bar() 
            np = floor(abs(x(2)-x(1))/xby)*floor(abs(y(2)-y(1))/yby)
            int_2D = int_1D(f=f1, x=x, by=xby)

            contains
                real(dp) function f1(inx)
                    real(dp) , intent(in) :: inx
                    tmp_x = inx
                    f1 = int_1D(f=int_f1, x=y, by=yby)
                end function f1
                
                real(dp) function int_f1(iny)
                    real(dp), intent(in) :: iny
                    int_f1 = f(tmp_x,iny)
                    npdone = npdone+1
                    if (lbar)  call progress_bar( dble(npdone)/dble(np) )
                end function int_f1 
        end function int_2D

        real(dp) function int_3D(f, x, y, z, xby, yby, zby,bar)
            real(dp), intent(in) :: x(2), y(2), z(2), xby, yby, zby
            procedure(r3fn) :: f
            real(dp) ::tmp_x, tmp_y
            integer :: np,npdone
            logical, intent(in), optional :: bar
            logical :: lbar
            lbar = merge( bar, .false., present(bar))
            npdone=0
            if (lbar)  print *, " <<< 3D integral >>>"
            if (lbar)  call set_progress_bar() 
            np = floor(abs(x(2)-x(1))/xby)*floor(abs(y(2)-y(1))/yby)*floor(abs(z(2)-z(1))/zby)
            int_3D = int_2D(f=f1, x=x,y=y,xby=xby,yby=yby,bar=.false.)            
        contains
            real(dp) function f1(inx, iny)
                real(dp), intent(in) :: inx, iny
                tmp_x = inx
                tmp_y = iny
                f1 = int_1D(f=int_f1, x=z, by=zby,bar=.false.)
            end function f1

            real(dp) function int_f1(inz)
                real(dp), intent(in) :: inz
                int_f1 = f(tmp_x, tmp_y, inz)
                npdone = npdone+1
                if (lbar)            call progress_bar( dble(npdone)/dble(np) )
            end function int_f1
        end function int_3D

end module integrate
