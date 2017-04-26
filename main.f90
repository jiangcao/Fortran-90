program main

    use types, only : dp
    use utils, only : seq, cbind
    use interpol1D 
    use output
    integer :: u
    real(dp) :: ys(8),xs(8),t
    real(dp) :: x(20),y(20)
    integer  :: j,i
    xs = seq(0.0_dp,1.0_dp,8) 
    ys = xs**2
    call output_2c(fn="xy.dat",dataset=cbind(xs,ys))


    print *, "=== random x ==="
    do i = 1,size(xs)
        print *, i, " | ", xs(i)
    enddo

    call srand(23154)
    do i = 1,size(xs)
        do j = 1,size(xs)
            if (rand()>0.5_dp) then 
                t = xs(i)
                xs(i) = xs(j)
                xs(j) = t
            endif
        enddo
    enddo
    ys = xs**2

    print *, "==="

    print *, "x=1.6"
    print *, find_neighbor(1.6_dp, n=3, lst=xs)
    print *, xs(find_neighbor(1.6_dp, n=3, lst=xs))

    print *,"=== interpolate ==="
    call interpol_set(8, set_x=xs, set_y=ys)

    do i = 1,20
        y(i) = i/40.0_dp+0.2_dp
        x(i) = interpol_lin_gy(y(i))         
    enddo

    call output_2c(fn="xy-in.dat", dataset=cbind(x,y))

    call interpol_freedata

end program main
