module io

    implicit none

    contains


!!****f* input_output/free_lun
! COPYRIGHT
!     (c) 2011 by Alessandro Cresti, IMEP-LAHC UMR5130, France
! NAME
!     free_lun - it provides a free unit number
! SYNOPSIS
!     lun = free_lun(lun_min,lun_max)
! FUNCTION
!     It gives the first free unit number for input/output operations.
!     The range of the units can be optionally chosen.
!     If no range is set, the unit is chosen between 10 and 10000.
! INPUT
!     lun_min     : minimum allowed value for the unit : integer, optional (default 10)
!     lun_max     : maximum allowed value for the unit : integer, optional (default 10000)
! OUTPUT
!     lun         : number of the free unit            : integer
! INTERNAL
!     opened      : true if the unit is not free       : logical
! HISTORY
!     16-06-2011  :  creation
! AUTHOR
!     Alessandro Cresti
! SOURCE
  function io_free_lun(lun_min,lun_max,lun_out) result (lun) !! free logic unit for files
    implicit none
    integer, intent(in), optional  :: lun_min,lun_max
    integer, intent(out), optional :: lun_out
    integer :: lun,min_lun,max_lun
    logical :: opened
    min_lun=10
    max_lun=10000
    if(present(lun_min).and.present(lun_max)) then
       min_lun=min(lun_min,lun_max)
       max_lun=max(lun_min,lun_max)
    elseif(present(lun_min)) then
       min_lun=lun_min
    elseif(present(lun_max)) then
       max_lun=lun_max
    end if
    do lun=min_lun,max_lun
       inquire(unit=lun,opened=opened)
       if (.not.opened) then 
           if (present(lun_out)) lun_out = lun
           return
       endif
    end do
    lun=-10
  end function io_free_lun
!!***


end module io
