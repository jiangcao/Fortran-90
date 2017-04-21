module interpol1D
! module to interpolate a function from a data set

    use types, only : dp

    implicit none

    private

    real(dp),save, allocatable :: dataset(:,:)
    real(dp),save :: xmin, xmax, ymin, ymax 
    integer  :: np

    public interpol_set, find_insert, find_neighbor
    public interpol_lin_gy, interpol_lin_fx
    public interpol_freedata
    
    contains

        ! linear interpolation from the two neighbor points 
        function interpol_lin_gy(y) result(x)
            implicit none
            real(dp), intent(in) :: y
            real(dp) :: x
            real(dp) :: p1(2),p2(2),dl(2),p(2)
            integer  :: ind(2)
            if ( (y<ymin) .or. (y>ymax) ) then 
                print *, "y out of range"
                call abort 
            endif 
            ind = find_neighbor (y, n=1, lst=dataset(:,2))
            p1 = dataset(ind(1),:)
            p2 = dataset(ind(2),:)
            dl = p2 - p1
            p  = dl / dl(2) * (y-p1(2)) + p1
            x  = p(1)
        end function interpol_lin_gy

        function interpol_lin_fx(x) result(y)
            implicit none
            real(dp) , intent(in) :: x
            real(dp) :: y
            real(dp) :: p1(2),p2(2),dl(2),p(2)
            integer  :: ind(2)
            if ( (x<xmin) .or. (x>xmax) ) then 
                print *, "x out of range"
                call abort 
            endif 
            ind = find_neighbor (x, n=1, lst=dataset(:,1))
            p1 = dataset(ind(1),:)
            p2 = dataset(ind(2),:)
            dl = p2 - p1
            p  = dl / dl(1) * (x-p1(1)) + p1
            y  = p(2)
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
