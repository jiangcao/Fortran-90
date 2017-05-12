module interpol1D
! module to interpolate a function from a data set

    use types, only : dp
    use output

    implicit none

    private

    real(dp),save, allocatable :: dataset(:,:)
    real(dp),save :: xmin, xmax, ymin, ymax 
    integer  :: np=0

    public interpol_set, find_insert, find_neighbor
    public interpol_lin_gy, interpol_lin_fx
    public interpol_freedata
    public interpol_getNP, interpol_getrange,interpol_writedata
    public interpol_adddatapoint
    
    contains

        subroutine interpol_adddatapoint(x,y,point)
            real(dp), optional, intent(in) :: x,y,point(2)
            real(dp), allocatable :: tmpds(:,:)
            if (present(x) .and. present(y)) then
                np = np+1
                allocate(tmpds(np,2))
                if (np>1) then
                    tmpds(1:np-1,:) = dataset
                    deallocate(dataset)
                endif
                tmpds(np,:) = (/x,y/)
                call move_alloc(tmpds, dataset)
            else if (present(point)) then
                np = np+1
                allocate(tmpds(np,2))
                if (np>1) then
                    tmpds(1:np-1,:) = dataset
                    deallocate(dataset)
                endif
                tmpds(np,:) = point 
                call move_alloc(tmpds, dataset)
            else
                print *,"error, missing data"
                call abort
            endif
            xmin = minval(dataset(:,1))
            xmax = maxval(dataset(:,1))
            ymin = minval(dataset(:,2))
            ymax = maxval(dataset(:,2))
        end subroutine interpol_adddatapoint


        function interpol_getNP() result(getnp)
            integer :: getnp
            getnp = np
        end function interpol_getNP

        function interpol_getrange() result(reg)
            real(dp) :: reg(4)
            reg = (/xmin, xmax, ymin, ymax/)
        end function interpol_getrange


        subroutine interpol_writedata(fn)
            implicit none
            character(len=*) , intent(in) :: fn
            call output_2c(fn, dataset)
        end subroutine interpol_writedata

        ! linear interpolation from the two neighbor points 
        function interpol_lin_gy(y) result(x)
            implicit none
            real(dp), intent(in) :: y
            real(dp) :: x
            real(dp) :: p1(2),p2(2),dl(2),p(2)
            integer  :: ind(2)
            if ( (y<ymin) .or. (y>ymax) ) then 
                x = 0.0_dp 
                print *, "y out of range"
                call abort
            else 
                ind = find_neighbor (y, n=1, lst=dataset(:,2))
                p1 = dataset(ind(1),:)
                p2 = dataset(ind(2),:)
                dl = p2 - p1
                p  = dl / dl(2) * (y-p1(2)) + p1
                x  = p(1)
            endif
        end function interpol_lin_gy

        function interpol_lin_fx(x) result(y)
            implicit none
            real(dp) , intent(in) :: x
            real(dp) :: y
            real(dp) :: p1(2),p2(2),dl(2),p(2)
            integer  :: ind(2)
            if ( (x<xmin-TINY(0.0_dp)) .or. (x>xmax+TINY(0.0_dp)) ) then 
                y = 0.0_dp
                print *, "x out of range"
                call abort
            else
                ind = find_neighbor (x, n=1, lst=dataset(:,1))
                p1 = dataset(ind(1),:)
                p2 = dataset(ind(2),:)
                dl = p2 - p1
                p  = dl / dl(1) * (x-p1(1)) + p1
                y  = p(2)
            endif
        end function interpol_lin_fx


        ! function returns the indices of the neighbor points in the data list 
        function find_neighbor(x,n,lst) result(ind)
            implicit none
            real(dp) , intent(in) :: x, lst(:) 
            integer  , intent(in) :: n
            integer :: ind(2*n) , nx 
            real(dp) :: xnp(2*n) ! neighbor points value
            integer  :: npi(2*n) ! neighbor points indices
            integer  :: i, m
            nx = size(lst)
            npi = 0
            xnp(1:n) = -HUGE(1.0_dp)
            xnp(n+1:2*n) =  HUGE(1.0_dp)
            do i = 1, nx
                if ( (lst(i) <= x ) .and. (lst(i) > xnp(1)) ) then 
                    m = find_insert(lst(i), lst=xnp(1:n)) 
                    xnp(1:m-1) = xnp(2:m)
                    xnp(m) = lst(i)
                    npi(1:m-1) = npi(2:m)
                    npi(m) = i
                endif
                if ( (lst(i) > x) .and. (lst(i) < xnp(2*n))) then 
                    m = find_insert(lst(i), lst=xnp(n+1:2*n)) 
                    m = m+n
                    xnp(m+2:2*n) = xnp(m+1:n*2-1)
                    xnp(m+1) = lst(i)
                    npi(m+2:2*n) = npi(m+1:n*2-1)
                    npi(m+1) = i
                endif
            enddo
            ind = npi
        end function find_neighbor
 

        function find_insert(x, lst) result(ind)
            implicit none
            real(dp), intent(in) :: x , lst(:)
            integer :: ind
            integer :: i, nx
            nx = size(lst)
            i = 1
            do while ( ( lst(i) < x ) .and. (i<= nx) )
                i = i+1
            enddo
            ind = i-1
        end function find_insert


            



        subroutine interpol_set(set_np,set_dataset,set_x,set_y)
            implicit none
            integer, intent(in) ::set_np
            real(dp), intent(in), optional :: set_dataset(:,:), set_x(:), set_y(:) 
            allocate ( dataset(set_np, 2) )
            np = set_np
            dataset = 0.0_dp
            if (present(set_dataset)) dataset = set_dataset
            if (present(set_x)) dataset(:,1) = set_x(:)
            if (present(set_y)) dataset(:,2) = set_y(:)
            xmin = minval(dataset(:,1)) 
            xmax = maxval(dataset(:,1)) 
            ymin = minval(dataset(:,2)) 
            ymax = maxval(dataset(:,2)) 
        end subroutine interpol_set


        subroutine interpol_freedata()
            deallocate(dataset)
        end subroutine interpol_freedata
end module interpol1D
