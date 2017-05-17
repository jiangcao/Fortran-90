module types
implicit none
private
public dp, hp
public rfn, r2fn, r3fn, rf3d, r2f3d,  r3f3d, rndfn, r3dfn
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
    real(dp) function rndfn(x)
        import
        real(dp), intent(in) :: x(:)
    end function rndfn
    function rf3d(x) result(fx)
        import
        real(dp), intent(in) :: x
        real(dp) :: fx(3)
    end function rf3d
    function r3f3d(x,y,z) result(fx)
        import
        real(dp), intent(in) :: x, y, z
        real(dp) :: fx(3)
    end function r3f3d
    function r2f3d(x,y) result(fx)
        import
        real(dp), intent(in) :: x, y
        real(dp) :: fx(3)
    end function r2f3d
    real(dp) function r3dfn(x)
        import
        real(dp), intent(in) :: x(3)
    end function r3dfn




end interface

end module
