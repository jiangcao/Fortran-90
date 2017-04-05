module types
implicit none
private
public dp, hp, rfn
integer, parameter :: dp=kind(0.d0), &          ! double precision
                      hp=selected_real_kind(15) ! high precision
interface
    real(dp) function rfn(x)
        import
        real(dp), intent(in) :: x 
    end function rfn
end interface

end module
