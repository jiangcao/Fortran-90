module input

    use types, only : dp

    implicit none

    private

    public input_2c
    public input_nc

    contains

        subroutine input_2c(fn, dataset)
            implicit none
            character(len=*), intent(in) :: fn
            real(dp), allocatable, intent(out) :: dataset(:,:)
            call input_nc(fn, 2, dataset)
        end subroutine input_2c


        subroutine input_nc(fn, ncol, dataset)
            implicit none
            character(len=*), intent(in) :: fn
            integer, intent(in) :: ncol
            real(dp), allocatable, intent(out) :: dataset(:,:)
            integer :: u, i, n, err
            character(len=1000) :: line
            open(newunit = u, file = trim(fn))
            n = 0
            ! How many lines in the files
            err = 0
            do while (err >= 0)
                read(u, '(A)', iostat=err)  line
                n = n +1                
            enddo
            n = n-1
            rewind(u)
            allocate(dataset(n, ncol))
            do i = 1,n
                read(u, *) dataset(i,:)
            enddo
            close(u)
        end subroutine input_nc



end module input
