module integrate
! Module to calculate numerically the integral of a function

    use types, only : dp, rfn, r2fn , r3fn, rf3d,r3f3d,r2f3d,r2f6d,rf6d
    use progressbar
    use mpi, only : idnode, mxnode 
    use utils, only : int2str

    implicit none

    private

    public :: int_1D,int_1Dg
    public :: int_2D,int_2Dg
    public :: int_3D!,int_3Dg

    contains

        ! Function returns the numerical integration of a function `f(x)` from x(1) to x(2)
        ! Parameters:
        ! - `f` the function to integrate
        ! - `x`   the integration interval
        ! - `by` the integration step
        real(dp) function int_1D(f, x, by, saveto,bar,mpi) 

            include "mpif.h"

            real(dp), intent(in) :: x(2), by
            real(dp) :: dx, fx, xi
            real(dp) :: linteg 
            logical, intent(in), optional :: bar
            logical :: lbar
            logical :: lmpi
            logical, intent(in), optional :: mpi
            integer :: i,np,npdone, n
            procedure(rfn)               :: f 
            character(len=*), intent(in), optional :: saveto
            integer :: ierr,u
            integer :: lidn, lnnd
            lbar = merge( bar, .false., present(bar))
            lmpi = merge( mpi, .false., present(mpi))
            if (lbar)     print *, " <<< 1D integral >>>"
            if (lbar)     call set_progress_bar() 
            if (present(saveto)) then 
                open(newunit = u, file = trim(saveto))
            endif
            int_1D = 0d0
            linteg = 0d0
            n = abs(x(2) - x(1))/by
            if (x(2) > x(1)) then
                dx = abs(by)
            else
                dx = -abs(by)
            endif
            lidn = 0
            lnnd = 1
            if (lmpi) then
                lidn = idnode
                lnnd = mxnode
            endif

            do i = 1+lidn, n, lnnd 
                xi = x(1) + dx*dble(i)
                fx = f( xi )
                if (lmpi) then
                    linteg = linteg + dx * fx
                else
                    int_1D = int_1D + dx * fx 
                endif
                if( present(saveto) .and. (.not. lmpi) ) then 
                    write(u, *) xi, fx 
                endif
                if ((lbar).and.(.not. lmpi))              call progress_bar( dble(i)/dble(n) )
                if ((lbar).and.(lmpi).and.(lidn==0))       call progress_bar( dble(i)*lnnd/dble(n) )
            enddo

            if (lmpi) then
                ! mpi summation
                call MPI_REDUCE(linteg,int_1D,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
                if (ierr .ne. MPI_SUCCESS) then
                    print *, "mpi failed",ierr
                    call abort
                endif
            endif

            if (present(saveto)) close(u)
        end function int_1D

 !       real(dp) function int_2D(f, x, y, xby, yby,saveto, bar)
 !           character(len=*), intent(in), optional :: saveto
 !           real(dp), intent(in) :: x(2), y(2), xby, yby
 !           procedure(r2fn) :: f
 !           real(dp) :: tmp_x
 !           integer ::u, np,npdone
 !           logical, intent(in), optional :: bar
 !           logical :: lbar
 !           if (present(saveto)) then 
 !               open(newunit = u, file = trim(saveto))
 !           endif
 !           lbar = merge( bar, .false., present(bar))
 !           npdone=0
 !           if (lbar)    print *, " <<< 2D integral >>>"
 !           if (lbar)    call set_progress_bar() 
 !           np = floor(abs(x(2)-x(1))/xby)*floor(abs(y(2)-y(1))/yby)
 !           int_2D = int_1D(f=f1, x=x, by=xby)
 !           if (present(saveto)) then 
 !               close(u)
 !           endif

 !           contains
 !               real(dp) function f1(inx)
 !                   real(dp) , intent(in) :: inx
 !                   tmp_x = inx
 !                   f1 = int_1D(f=int_f1, x=y, by=yby)
 !               end function f1
 !               
 !               real(dp) function int_f1(iny)
 !                   real(dp), intent(in) :: iny
 !                   int_f1 = f(tmp_x,iny)
 !                   if (present(saveto)) then 
 !                       write(u,*) tmp_x, iny, int_f1
 !                   endif
 !                   npdone = npdone+1
 !                   if (lbar)  call progress_bar( dble(npdone)/dble(np) )
 !               end function int_f1 
 !       end function int_2D

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
            end function int_f1
        end function int_3D


        function int_1Dg(f, x, by, saveto,bar,mpi) result(intf) 
            implicit none

            include "mpif.h"

            real(dp), intent(in) :: x(2), by
            real(dp) :: dx, fx(6), xi
            real(dp) :: intf(6) 
            logical, intent(in), optional :: bar
            logical :: lbar
            integer :: i,np,npdone, n
            procedure(rf6d)               :: f 
            character(len=*), intent(in), optional :: saveto
            integer :: u
            real(dp) :: linteg(6)
            logical :: lmpi
            logical, intent(in), optional :: mpi
            integer :: lidn, lnnd, ierr
            lbar = merge( bar, .false., present(bar))
            lmpi = merge( mpi, .false., present(mpi))
            if (lbar)     print *, " <<< 1D integral >>>"
            if (lbar)     call set_progress_bar() 
            if (present(saveto)) then 
                open(newunit = u, file = trim(saveto))
            endif
            intf = 0d0
            linteg = 0d0
            n = abs(x(2) - x(1))/by
            if (x(2) > x(1)) then
                dx = abs(by)
            else
                dx = -abs(by)
            endif
            lidn = 0
            lnnd = 1
            if (lmpi) then
                lidn = idnode
                lnnd = mxnode
            endif

            do i = 1+lidn, n, lnnd 
                xi = x(1) + dx*dble(i)
                fx = f( xi )
                if (lmpi) then
                    linteg = linteg + dx * fx
                else
                    intf = intf + dx * fx 
                endif
                if( present(saveto) .and. (.not. lmpi) ) then 
                    write(u, *) xi, fx 
                endif
                if ((lbar).and.(.not. lmpi))              call progress_bar( dble(i)/dble(n) )
                if ((lbar).and.(lmpi).and.(lidn==0))       call progress_bar( dble(i+lnnd)/dble(n) )
            enddo

            if (lmpi) then
                ! mpi summation
                call MPI_REDUCE(linteg,intf,6,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
                if (ierr .ne. MPI_SUCCESS) then
                    print *, "mpi failed",ierr
                    call abort
                endif
                call MPI_BCAST(intf,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
                if (ierr .ne. MPI_SUCCESS) then
                    print *, "mpi failed",ierr
                    call abort
                endif
            endif

            if (present(saveto)) close(u)
        end function int_1Dg


        ! function int_3Dg(f, x, y, z, xby, yby, zby,bar) result (intf)
        !     real(dp), intent(in) :: x(2), y(2), z(2), xby, yby, zby
        !     procedure(r3f3d) :: f
        !     real(dp) ::tmp_x, tmp_y
        !     real(dp) :: intf(3)
        !     integer :: np,npdone
        !     logical, intent(in), optional :: bar
        !     logical :: lbar
        !     lbar = merge( bar, .false., present(bar))
        !     npdone=0
        !     if (lbar)  print *, " <<< 3D integral >>>"
        !     if (lbar)  call set_progress_bar() 
        !     np = floor(abs(x(2)-x(1))/xby)*floor(abs(y(2)-y(1))/yby)*floor(abs(z(2)-z(1))/zby)
        !     intf = int_2Dg(f=f1, x=x,y=y,xby=xby,yby=yby,bar=.false.)            
        ! contains
        !     function f1(inx, iny) result(of1)
        !         real(dp), intent(in) :: inx, iny
        !         real(dp) :: of1(3)
        !         tmp_x = inx
        !         tmp_y = iny
        !         of1 = int_1Dg(f=int_f1, x=z, by=zby,bar=.false.)
        !     end function f1

        !     function int_f1(inz) result(intf1)
        !         real(dp) :: intf1(3)
        !         real(dp), intent(in) :: inz
        !         intf1 = f(tmp_x, tmp_y, inz)
        !         npdone = npdone+1
        !         if (lbar)            call progress_bar( dble(npdone)/dble(np) )
        !     end function int_f1
        ! end function int_3Dg

        function int_2D(f, x, y, xby, yby,saveto, bar, mpi) result(intf)
            implicit none

            include "mpif.h"

            character(len=*), intent(in), optional :: saveto
            real(dp), intent(in) :: x(2), y(2), xby, yby
            procedure(r2fn) :: f
            real(dp) ::  fx
            real(dp) :: linteg,intf , xi, yi, dx, dy            
            integer ::u, n, i, nx,ny
            logical, intent(in), optional :: bar
            logical :: lbar
            logical :: lmpi
            logical, intent(in), optional :: mpi
            integer :: lidn, lnnd, ierr

            lmpi = merge( mpi, .false., present(mpi))
            if (present(saveto)) then 
                open(newunit = u, file = trim(saveto)//int2str(idnode)//'.dat')
            endif
            lbar = merge( bar, .false., present(bar))
            if ((lbar).and.(idnode==0))    print *, " <<< 2D integral >>>"
            if (lbar)    call set_progress_bar() 
  
            intf = 0d0
            linteg = 0d0
            nx = floor(abs(x(2) - x(1))/xby)
            ny = floor(abs(y(2) - y(1))/yby)
            if (x(2) > x(1)) then
                dx = abs(xby)
            else
                dx = -abs(xby)
            endif
            if (y(2) > y(1)) then
                dy = abs(yby)
            else
                dy = -abs(yby)
            endif
            lidn = 0
            lnnd = 1
            if (lmpi) then
                lidn = idnode
                lnnd = mxnode
            endif
            n = nx*ny

            do i = 1+lidn, n, lnnd 
                xi = x(1) + dx*dble(((i-1)/ny)+1)
                yi = y(1) + dy*dble(mod(i-1,ny)+1)                
                fx = f( xi , yi )
                if (lmpi) then
                    linteg = linteg + dx * dy * fx
                else
                    intf = intf + dx * dy * fx 
                endif
                if( present(saveto) ) then 
                    write(u, *) xi, yi, fx 
                endif
                if ((lbar).and.(.not. lmpi))              call progress_bar( dble(i)/dble(n) )
                if ((lbar).and.(lmpi).and.(lidn==0))       call progress_bar( dble(i+lnnd)/dble(n) )
            enddo

            if (lmpi) then
                ! mpi summation
                call MPI_REDUCE(linteg,intf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
                if (ierr .ne. MPI_SUCCESS) then
                    print *, "mpi failed",ierr
                    call abort
                endif
                call MPI_BCAST(intf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
                if (ierr .ne. MPI_SUCCESS) then
                    print *, "mpi failed",ierr
                    call abort
                endif
            endif

            if (present(saveto)) close(u)

        end function int_2D


        function int_2Dg(f, x, y, xby, yby,saveto, bar, mpi) result(intf)
            implicit none

            include "mpif.h"

            character(len=*), intent(in), optional :: saveto
            real(dp), intent(in) :: x(2), y(2), xby, yby
            procedure(r2f6d) :: f
            real(dp) ::  fx(6)
            real(dp) :: linteg(6),intf(6) , xi, yi, dx, dy            
            integer ::u, n, i, nx,ny
            logical, intent(in), optional :: bar
            logical :: lbar
            logical :: lmpi
            logical, intent(in), optional :: mpi
            integer :: lidn, lnnd, ierr

            lmpi = merge( mpi, .false., present(mpi))
            if (present(saveto)) then 
                open(newunit = u, file = trim(saveto)//int2str(idnode)//'.dat')
            endif
            lbar = merge( bar, .false., present(bar))
            if ((lbar).and.(idnode==0))    print *, " <<< 2D integral >>>"
            if (lbar)    call set_progress_bar() 
  
            intf = 0d0
            linteg = 0d0
            nx = floor(abs(x(2) - x(1))/xby)
            ny = floor(abs(y(2) - y(1))/yby)
            if (x(2) > x(1)) then
                dx = abs(xby)
            else
                dx = -abs(xby)
            endif
            if (y(2) > y(1)) then
                dy = abs(yby)
            else
                dy = -abs(yby)
            endif
            lidn = 0
            lnnd = 1
            if (lmpi) then
                lidn = idnode
                lnnd = mxnode
            endif
            n = nx*ny

            do i = 1+lidn, n, lnnd 
                xi = x(1) + dx*dble(((i-1)/ny)+1)
                yi = y(1) + dy*dble(mod(i-1,ny)+1)                
                fx = f( xi , yi )
                if (lmpi) then
                    linteg = linteg + dx * dy * fx
                else
                    intf = intf + dx * dy * fx 
                endif
                if( present(saveto) ) then 
                    write(u, *) xi, yi, fx 
                endif
                if ((lbar).and.(.not. lmpi))              call progress_bar( dble(i)/dble(n) )
                if ((lbar).and.(lmpi).and.(lidn==0))       call progress_bar( dble(i+lnnd)/dble(n) )
            enddo

            if (lmpi) then
                ! mpi summation
                call MPI_REDUCE(linteg,intf,6,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
                if (ierr .ne. MPI_SUCCESS) then
                    print *, "mpi failed",ierr
                    call abort
                endif
                call MPI_BCAST(intf,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
                if (ierr .ne. MPI_SUCCESS) then
                    print *, "mpi failed",ierr
                    call abort
                endif
            endif

            if (present(saveto)) close(u)



        end function int_2Dg


end module integrate
