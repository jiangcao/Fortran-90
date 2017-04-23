module output

    use types, only : dp

    implicit none

    private

    public output_2c, show_2c

    contains

        subroutine output_2c(fn, dataset)
            implicit none
            character(len=*), intent(in) :: fn
            real(dp), intent(in) :: dataset(:,:)
            integer :: u, i
            open(newunit = u, file = trim(fn))
            do i = 1, size(dataset,1)
                write(u, *) dataset(i,:)
            enddo
            close(u)
        end subroutine output_2c

        subroutine show_2c(dataset)
            implicit none
            real(dp), intent(in) :: dataset(:,:)
            integer :: i
            print *, "======"
            do i = 1, size(dataset,1)
                print *, dataset(i,:)
            enddo
            print *, "======"
        end subroutine show_2c


end module output
