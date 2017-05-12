module progressbar

    use types , only : dp


    private 

    real(dp), save :: perc

    public set_progress_bar, progress_bar

contains


    subroutine set_progress_bar()
        perc = 0.0_dp
        write(*,"(A)", advance='no') "    0%"
    end subroutine set_progress_bar

    subroutine progress_bar(progress)
        real(dp), intent(in) :: progress
        integer :: i, n
        CHARACTER :: CR = CHAR(8)    ! carriage backspace character
        n = floor((progress - perc)*100.0_dp) 
        if (n >= 1) then 
            perc = floor(progress*1e2_dp)/1e2_dp
            write(*,"(4A)", advance='no') CR,CR,CR,CR
            do i = 1, n 
                write(*,"(A)", advance='no') "*"
            enddo
            write(*,"(I3,A)", advance='no') floor(progress*1e2_dp),"%"
        endif
    end subroutine progress_bar


end module
