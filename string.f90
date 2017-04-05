module string
!-----------------------------------------------------------------------------------------
!! This module implements operations related to Strings.
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
implicit none
integer ,  parameter :: string_length=1000
type type_string 
	character(len=string_length) :: s 
end type type_string

interface split
 	module Procedure split_string
end interface split

interface assignment(=)
	module Procedure assign_string,assign_stringToChars
end interface assignment(=)

interface get 
	module Procedure get_string
end interface get 


contains

!-----------------------------------------------------------------------------------------
pure function get_string(ss,i) result(s)
!! Procedure returns the i-th string in an array of string, and removes the useless spaces.
implicit none 
type(type_string),intent(in) :: ss(:)
integer, intent(in) :: i 
character(len=string_length) :: s
s = adjustl(ss(i)%s)
!-----------------------------------------------------------------------------------------
end function get_string



!-----------------------------------------------------------------------------------------
recursive subroutine split_string(s,subs,sep) 
!! Procedure splits a string (s) into substrings (subs) based on the seperator character 
!! (sep), and returns an array of the substrings. 
!!
!! @note Each substring is stored in a fixed-length string of length=1000, `trim` is 
!! necessary to remove the useless spaces.
implicit none 
character(len=*), intent(in) :: s 
character(len=1), intent(in), optional :: sep 
type(type_string), allocatable :: subs(:)
type(type_string), allocatable :: subs_sub(:)
character(len=1) :: sep_op
integer :: i
if (allocated(subs)) deallocate(subs)
sep_op = merge(sep,' ',present(sep))
i = index(s,sep_op)
if (i == len(s)) then 
	allocate(subs(1))
	subs(1)%s = s(1:len(s)-1)
else
	call split_string(s(i+1:),subs_sub,sep_op)
	allocate(subs(size(subs_sub)+1))
	FORALL (i=1:size(subs_sub)) subs(i+1) = subs_sub(i)
	subs(2:) = subs_sub(:)
	subs(1)%s = s(:i-1)
	deallocate(subs_sub)
endif 
!-----------------------------------------------------------------------------------------
end subroutine split_string

!-----------------------------------------------------------------------------------------
pure subroutine assign_string(a,b)
type(type_string) ,intent(in) :: b 
type(type_string) ,intent(out):: a 
a%s = b%s 
!-----------------------------------------------------------------------------------------
end subroutine assign_string

!-----------------------------------------------------------------------------------------
pure subroutine assign_stringToChars(a,b)
type(type_string) ,intent(in) :: b 
character(len=*),intent(inout):: a 
a = b%s 
!-----------------------------------------------------------------------------------------
end subroutine assign_stringToChars

end module string