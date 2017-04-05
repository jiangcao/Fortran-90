module Table
!-----------------------------------------------------------------------------------------
!! This module defines the structure [[type_Table]] which contains the 
!! physical information of the different materials. The main data structure is  
!! a 2-D array, in which the column represents different properties and row represents 
!! materials. The values can be accessed by giving the name of material and of property 
!! via the ``get`` interface to [[getProperty]] and [[getAllProperties]] functions.
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
use string
!-----------------------------------------------------------------------------------------
implicit none

private
public :: type_Table
public :: free,ReadTxt,get

type type_Table
	private
	real(8), allocatable :: prop(:,:) 
	!! The Properties table, a 2-D array
	character(len=100),allocatable :: propName(:)
	!! The properties' name
	character(len=100),allocatable :: matName(:)
	integer :: Nmat 
	integer :: Nprop
end type type_Table

interface free 
	module procedure freePropertiesTable
end interface free

interface ReadTxt 
	module procedure ReadPropertiesFromTxt
end interface ReadTxt

interface get 
	module Procedure getProperty, getAllProperties
end interface get

contains

!-----------------------------------------------------------------------------------------
pure function getAllProperties(tab,material) result(r)
!! Function returns an array containing all the properties of a material in a table
!! (tab).
implicit none 
type(type_Table), intent(in) :: tab 
character(len=*),intent(in) :: material
real(8) :: r(tab%Nprop)
integer :: i
do i=1,tab%Nmat 
	if (trim(tab%matName(i)) == trim(material)) exit 
enddo
if (i>tab%Nmat) then 
	r = HUGE(1.0d0)
	return 
endif
r = tab%prop(:,i)
!-----------------------------------------------------------------------------------------
end function getAllProperties



!-----------------------------------------------------------------------------------------
pure function getProperty(tab,material,property) result(r)
!! Function returns the value of property of material in the table, if the 
!! value is not found returns HUGE(1.0D0).
implicit none
type(type_Table), intent(in) :: tab 
character(len=*),intent(in)  :: material,property 
real(8) :: r
integer :: i,j
do i=1,tab%Nmat 
	if (trim(tab%matName(i)) == trim(material)) exit 
enddo
do j=1,tab%Nprop
	if (trim(tab%propName(j))== trim(property)) exit 
enddo
if ((i>tab%Nmat) .or. (j>tab%Nprop)) then 
	r = HUGE(1.0d0)
	return 
endif
r = tab%prop(j,i)
!-----------------------------------------------------------------------------------------
end function getProperty

!-----------------------------------------------------------------------------------------
subroutine ReadPropertiesFromTxt(fname,tab)
!! Procedure reads in the properties from a text file with the following form
!!
!!```
!!         |  property1    |    property2    |    property3  |  ...
!! --------+---------------+-----------------+---------------+-----------
!!  mat1   |         [1]   |                 |               |  ...  
!!  mat2   |               |           [2]   |               |  ...  
!! ======================================================================
!![1] 
!![2]
!!```
!!
!! References are placed in the end of the table, and seperated with `=====`. The procedure
!! will just ignore all the lines after the '=' 
!!  
implicit none
type(type_Table), intent(out) :: tab 
!! Read-in properties table
character(len=*),intent(in)             :: fname
!! File name string
integer,parameter :: handle = 46523
character(len=5000) :: line
character(len=20)   :: numstr
type(type_string),allocatable :: subs(:)
integer :: i, err,nl,j,k
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
open(unit=handle,file=fname)
nl = 0
! How many lines in the file
do  
	read(handle,'(A)',IOSTAT=err) line
	! End of the File
	if (err < 0) EXIT
	! End of the table
	if (line(1:1) == "=") EXIT
	nl = nl+1
enddo
print '(A,I3)','Number of material=',nl-2
! go back to begining of file
rewind(handle)
! read 1st line and decide how many properties
read(handle,'(A)') line
call split(trim(adjustl(line)),subs,'|')
tab%Nprop = size(subs)-1
allocate(tab%propName(tab%Nprop))
FORALL (i=2:size(subs)) tab%propName(i-1) = adjustl(get(subs,i))
print '(A10,100A12)','',(tab%propName(i) ,i=1,tab%Nprop)
read(handle,'(A)') line
print '(A)',Repeat('=',tab%Nprop*12+5) 
! ==============
tab%Nmat = nl-2
allocate(tab%prop(tab%Nprop,tab%Nmat))
allocate(tab%matName(tab%Nmat))
do i = 1,tab%Nmat 
	read(handle,'(A)',IOSTAT=err) line
	call split(trim(adjustl(line)),subs,'|')
	tab%matName(i) = adjustl(get(subs,1))
	! read the values
	do j = 1,tab%Nprop
		numstr = get(subs,j+1)
		k = index(numstr,'[')
        if (k .ne. 0) numstr = trim(numstr(1:k-1))
		read(numstr,'(E)') tab%prop(j,i)
	enddo
	print '(A5,100E12.3)',tab%matName(i),tab%prop(:,i)
enddo
close(handle)
if (allocated(subs)) deallocate(subs, stat=err)
if (err /= 0) print *, "subs: Deallocation request denied"
!-----------------------------------------------------------------------------------------
end subroutine ReadPropertiesFromTxt


	
!-----------------------------------------------------------------------------------------
subroutine freePropertiesTable(tab)
!! Procedure frees the memory allocated to the properties table (tab)
implicit none
type(type_Table), intent(inout) :: tab
integer :: err
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
if (allocated(tab%prop)) deallocate(tab%prop, stat=err)
if (err /= 0) print *, "tab%prop: Deallocation request denied"
if (allocated(tab%propName)) deallocate(tab%propName, stat=err)
if (err /= 0) print *, "tab%propName: Deallocation request denied"
if (allocated(tab%matName)) deallocate(tab%matName, stat=err)
if (err /= 0) print *, "tab%matName: Deallocation request denied"
!-----------------------------------------------------------------------------------------
end subroutine freePropertiesTable

end module Table
