program main

    use input
    use output
    use utils
    use types, only : dp
    use interpol1D

    implicit none

    real(dp) , allocatable :: d(:,:), x(:), y(:)
    integer :: i


    call input_2c(fn="xy.dat", dataset = d)

    call show_2c(d)

    call interpol_set(size(d,1),set_x=d(:,1), set_y=d(:,2))

    allocate(x(20))
    allocate(y(20))

    x = seq(0.2_dp,0.7_dp,20)
    do i = 1, 20
        y(i) = interpol_lin_fx(x(i))
    enddo

    call output_2c(fn="xy-lin.dat", dataset=cbind(x,y))

    deallocate(d)

    call interpol_freedata


end program main
