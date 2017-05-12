module types
implicit none
private
public dp, hp, rfn, r2fn, r3fn
integer, parameter :: dp=kind(0.d0), &          ! double precision
                      hp=selected_real_kind(15) ! high precision
interface
    real(dp) function rfn(x)
        import
        real(dp), intent(in) :: x 
    end function rfn
    real(dp) function r2fn(x,y)
        import
        real(dp), intent(in) :: x, y
    end function r2fn
    real(dp) function r3fn(x,y,z)
        import
        real(dp), intent(in) :: x, y, z
    end function r3fn

end interface

end module
