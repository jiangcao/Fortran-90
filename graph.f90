module graph 
!-----------------------------------------------------------------------------------------
!! Library cuts a graph into slices which have only connections with the left
!! and right neighbor slices. [[slice]] is the main subroutine to call.
!-----------------------------------------------------------------------------------------
use types, only : dp
implicit none

private 
public :: AddEdge,ReadTxt,SaveTxt,SaveTxt2col,testSlicing,slice
public :: getPointsInAllSlices,getPointsInSlice

interface slice 
	module procedure slice_1contact, slice_2contacts	
end interface 



interface SaveTxt
	module procedure SaveTxtGraph,SaveTxtAllConnectedPoints
end interface

interface ReadTxt
	module procedure ReadGraphFromText
end interface ReadTxt

contains



!-----------------------------------------------------------------------------------------
function getPointsInAllSlices(S) result(v)
!! Function returns the points in all the slices
implicit none
integer,CONTIGUOUS , intent(in) :: S(:,:)
!! Slices information
integer :: v(sum(S(1,:))-size(S,2))
integer :: i,n
n=0
do i = 1, size(S,2)
	v(n+1:n+S(1,i)-1) = S(2:S(1,i),i)
	n = n+S(1,i)-1
enddo
!-----------------------------------------------------------------------------------------
end function getPointsInAllSlices


!-----------------------------------------------------------------------------------------
function getPointsInSlice(S,i) result(v)
!! Function returns the points in a slice number `i`
implicit none
integer,CONTIGUOUS , intent(in) :: S(:,:)              !! Slices information
integer, intent(in) :: i                               !! Slice number
integer :: v(S(1,i)-1)
v(:) = S(2:S(1,i),i)
!-----------------------------------------------------------------------------------------
end function getPointsInSlice

!-----------------------------------------------------------------------------------------
subroutine SaveSlicesTxt(handle,S,X,Y,Z)
!! Procedure saves the slice information into a text file
implicit none
integer, intent(in) :: handle                           !! file unit number
integer,CONTIGUOUS , intent(in) :: S(:,:)               !! Slices information
real(dp),CONTIGUOUS , intent(in) :: X(:)
real(dp),CONTIGUOUS , intent(in) :: Y(:)
real(dp),CONTIGUOUS , intent(in) :: Z(:)
integer :: i,j
write (handle,*) '     X         Y         Z       Slice# '
do i = 1,size(S,2)
	do j = 2,S(1,i) 
		write(handle,'(3E15.5,I10)') X(S(j,i)),Y(S(j,i)),Z(S(j,i)),i 
	enddo
enddo
!-----------------------------------------------------------------------------------------
end subroutine SaveSlicesTxt


!-----------------------------------------------------------------------------------------
subroutine AddEdge(g,ij)
!! Subroutine update the graph data (g) by adding an edge from node-i to node-j.
!! The memory space for (g) is already allocated before calling this procedure
implicit none	
integer,CONTIGUOUS , intent(out) :: g(:,:)  
integer, intent(in) :: ij(2)
! Find if this Edge is already present in the graph, avoiding duplicates
if (ANY(g(2:g(1,ij(1)) ,ij(1)) == ij(2))) then 
	return 
else 
	g(1,ij)=g(1,ij)+1
	g(g(1,ij(2)),ij(2)) = ij(1)
	g(g(1,ij(1)),ij(1)) = ij(2)
endif
!-----------------------------------------------------------------------------------------
end subroutine AddEdge


!-----------------------------------------------------------------------------------------
subroutine ReadGraphFromText(fname,g)
!! Subroutine for reading in the graph data from an ASCII text file. 
!! 
!! @note The graph connectivity table is stored in a 2D integer array `g(:,:)` which
!! is allocated inside this subroutine, so remember to deallocate it outside. The 2nd 
!! index of `g(:,:)` refers to the point ID and 'g(1,:)-1' is the number of connections 
!! of the point. The text file should contain 2 columns of integers which are the point 
!! IDs of 2 connected points in the graph.
!!
implicit none 
integer, allocatable, intent(out)   :: g(:,:)        !! Graph connectivity table. 
character(len=*), intent(in)        :: fname         !! Unit number of the input text file
integer, allocatable                :: gn(:)         !! Temporary array for sizing the Graph Table
integer                             :: i,j,k         !! Looping variables
integer                             :: NL            !! Number of lines in the text file
integer                             :: IO            !! IO state during reading
integer                             :: M             !! Number of points in the graph 
integer, parameter :: handle = 675
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
open(unit = handle, file=fname)
M = -HUGE(1)
NL = 0
! Read through the file first to know the number of points in the graph
do 
	read(handle,*,IOSTAT = IO) i,j 
	if (IO < 0) exit 
	if (max(i,j)>M) M=max(i,j)
	NL = NL+1
enddo
allocate(gn(M))
gn = 0
! Read again the file to know about the size of the Graph Table
rewind handle
do k =1,NL
	read(handle,*) i,j 
	gn((/i,j/)) = gn((/i,j/)) +1
enddo
allocate(g(maxval(gn)+1,M))
g(:,:) = 0
g(1,:) = 1
deallocate(gn)
! Read last time the file for the Graph Table
rewind handle
do k=1,NL
	read(handle,*,IOSTAT = IO) i,j 
	call AddEdge(g, (/i,j/))
enddo 
close(handle)
!-----------------------------------------------------------------------------------------
endsubroutine ReadGraphFromText



!-----------------------------------------------------------------------------------------
subroutine SaveTxtGraph(handle,g)
!! Subroutine for saving the graph data into an ASCII text file. 
!! 
!! @note The 1st column is the point ID, and the following columns are the connecting points' ID.
!!
implicit none 
integer,CONTIGUOUS , intent(in)           :: g(:,:)  !! Graph connectivity table. 
integer, intent(in), optional :: handle		         !! Unit number of the input text file
character(len=1024)           :: str                 !! String for FORMAT writing
integer 			          :: i                   !! Looping variable
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
write (str, "(I10)") size(g,1)+1
if (present(handle)) then 
	do i=1,size(g,2)
	write(handle, '('//trim(adjustl(str))//'i10)') i, g(2:g(1,i),i)
	enddo
else 
	do i=1,size(g,2)
	print '('//trim(adjustl(str))//'i10)', i, g(2:g(1,i),i)
	enddo
endif
!-----------------------------------------------------------------------------------------
end subroutine SaveTxtGraph





!-----------------------------------------------------------------------------------------
subroutine SaveTxt2col(handle,g)
!! Subroutine for saving the graph data into an ASCII text file with 2 columns format 
!! 
!! @note The 1st and 2nd columns are point IDs. Each connection is represented by 1 line.
!!
implicit none 
integer,CONTIGUOUS , intent(in)              :: g(:,:)  !! Graph connectivity table. 
integer, intent(in), optional    :: handle              !! Unit number of the input text file
integer                          :: i,j                 !! Looping variable
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
do i=1,size(g,2)
	do j=2,g(1,i)
		if (i .le. g(j,i)) then 
			if (present(handle)) then 
				write(handle, '(2i10)') i, g(j,i)
			else 
				write(*, '(2i10)') i, g(j,i)
			endif 
		endif
	enddo
enddo
!-----------------------------------------------------------------------------------------
end subroutine SaveTxt2col


!-----------------------------------------------------------------------------------------
subroutine SaveTxtAllConnectedPoints(handle,g,X,Y,Z)
implicit none 
integer,CONTIGUOUS , intent(in)  :: g(:,:)        !! Graph connectivity table. 
integer, intent(in), optional    :: handle        !! Unit number of the input text file
real(dp),CONTIGUOUS , intent(in) :: X(:)
real(dp),CONTIGUOUS , intent(in) :: Y(:)
real(dp),CONTIGUOUS , intent(in) :: Z(:)
integer                          :: i,j           
!! Looping variable
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
if (present(handle)) then 
	write (handle,*) '     X         Y         Z        Point#'
else 
	print *, '     X         Y         Z         Point#'
endif
do i=1,size(g,2)
	do j=2,g(1,i)
		if (i .le. g(j,i)) then 
			if (present(handle)) then 
				write(handle, '(3E10.2,i10)') X(i),Y(i),Z(i), i
				write(handle, '(3E10.2,i10)') X(g(j,i)),Y(g(j,i)),Z(g(j,i)), g(j,i)
			else 
				write(*, '(3E10.2,i10)') X(i),Y(i),Z(i), i
				write(*, '(3E10.2,i10)') X(g(j,i)),Y(g(j,i)),Z(g(j,i)), g(j,i)
			endif 
		endif
	enddo
enddo
!-----------------------------------------------------------------------------------------
end subroutine SaveTxtAllConnectedPoints


!-----------------------------------------------------------------------------------------
function dist(g,E) result(D)
!! Function returns the distance of all the points to an edge of the graph.
!! Distance means the nomber of steps needed to arrive at this point starting from the edge
!! 
implicit none
integer,CONTIGUOUS , intent(in):: g(:,:)                !! Graph connectivity table. 
integer,CONTIGUOUS , intent(in):: E(:)                  !! Edge points' IDs
integer	:: D(size(g,2))                                 !! Distance table
integer :: ST(0:size(g,2)*2)                            !! Circular Stack for BFS algorithm
integer :: LST                                          !! Length of the Stack
integer :: SLP                                          !! Lower Stack Pointer
integer :: SHP                                          !! Higher Stack Pointer
integer :: PID                                          !! A point ID
integer :: PID2	                                        !! A point ID
integer :: i                                            !! Looping integer
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
D(:) = HUGE(1)
D(E(:)) = 0 
SLP=1
SHP=size(E)+1
ST = 0
ST(SLP:SHP-1) = E(:)
LST = size(ST)
do while (SHP .ne. SLP)
	PID = ST(SLP)								! Take 1 point out from the Stack
	do i = 2, g(1,PID)
		PID2 = g(i,PID)							! Look at its connections
		if (D(PID2) .eq. HUGE(1)) then 
			D(PID2) = D(PID) +1
			ST(SHP) = PID2						! Put the point in Stack for further moves
			SHP = MOD(SHP+1,LST)
		endif		
	enddo
	SLP = MOD(SLP+1,LST)
enddo
!-----------------------------------------------------------------------------------------
end function dist	



!-----------------------------------------------------------------------------------------
function part(DL,DR) result(P)
!! Function returns the partition of a graph depending on the distances to the edges 
!! obtained from [[edge]]
!!
implicit none
integer,CONTIGUOUS ,intent(in)				:: DL(:) !! Distance to the left edge obtained from [[edge]]
integer,CONTIGUOUS ,intent(in)				:: DR(:) !! Distance to the right edge obtained from [[edge]]
integer 						:: P(size(DL))	     !! Partition, 1 for left part, and 2 for right part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
P = 2
where(DL < DR) P=1
!-----------------------------------------------------------------------------------------
end function part



!-----------------------------------------------------------------------------------------
function PtsInParts(P,PID) result(PT)
!! Function returns a list of points' IDs in a partition obtained from [[part]] 
!!
implicit none
integer,CONTIGUOUS ,intent(in)	:: P(:)					!! Partition from [[part]] 
integer,intent(in)        		:: PID   				!! Part Number 
integer 				        :: PT(COUNT(P==PID))	!! List of Points in the part #PID
integer 						:: i,j
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
j = 1
do i = 1,size(P)
	if (P(i) == PID) then 
		PT(j) = i 
		j = j+1
	endif 
enddo
!-----------------------------------------------------------------------------------------
end function PtsInParts




!-----------------------------------------------------------------------------------------
function edge(g,P) result(E)
!! Function returns the edges of a partition of the graph. 
!! The partition is obtained from function [[part]]
!!
implicit none
integer,CONTIGUOUS , intent(in)	:: g(:,:)				  !! Graph connectivity table. 
integer,CONTIGUOUS , intent(in)	:: P(:)	    			  !! Partition Information from [[part]] 
integer 					:: E(size(P)+2,maxval(P))     !! Edge of each part
integer 					:: i,j   				   	  !! Loop integers
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
E = 0
E(1,:) = 1
do i = 1,size(g,2)							! Loop over all the points
	do j = 2,g(1,i)							! Look up its cunnecting points
		if (P(g(j,i)) .ne. P(i)) then   	! They belong to different parts
			E(1,P(i)) = E(1,P(i)) + 1
			E(E(1,P(i)),P(i)) = i 			! Add this point into the Edge list of part-i
			exit
		endif 
	enddo
enddo
!-----------------------------------------------------------------------------------------
end function edge





!-----------------------------------------------------------------------------------------
function subgraph(g,PT) result(sg)
!! Function returns the sub-graph of a graph by providing a list of points in the sub-graph
!! @note The sub-graph has a different point index than the original graph, since it contains 
!! only the selected points. 
!! To find the original point index, one has to refer to array PT. For exemple, point `1` 
!! in the sub-graph corresponds the point `PT(1)` in original graph.
!!
implicit none 
integer,CONTIGUOUS , intent(in) :: g(:,:)	            !! Graph connectivity table. 
integer,CONTIGUOUS , intent(in) :: PT(:)	            !! List of points in the sub-graph
integer 			      :: sg(size(g,1)+1,size(PT))	!! Sub-Graph connectivity table
integer 			      :: i,j,k 					    !! Looping integers
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
sg = 0
sg(1,:) = 1
do i = 1, size(PT)
	do j = 2, g(1,PT(i))
		if (ANY(PT == g(j,PT(i)))) then ! point PT(i) connects to a point belonging to PT
			sg(1,i) = sg(1,i) + 1       ! append the connection list of sub-graph
			sg(sg(1,i),i) = index(PT,g(j,PT(i)))  ! Add the connected point, with new index in sub-graph
		endif 
	enddo 
enddo
!-----------------------------------------------------------------------------------------
end function subgraph



!-----------------------------------------------------------------------------------------
pure function index(A,y) result(x)
!! Function returns the index of one value y in a list A
!!
implicit none
integer, intent(in)			:: A(:), y
integer 					:: x,i			
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
do i=1,size(A)			
	if (y == A(i)) exit 
enddo 
x = i
!-----------------------------------------------------------------------------------------
end function index



!-----------------------------------------------------------------------------------------
function indexes(A,y) result(x)
!! Function returns the index of an array of values y in the list A by calling [[index]]
!!
implicit none
integer,CONTIGUOUS , intent(in)	:: A(:), y(:)
integer 					:: x(size(y)),i			
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
FORALL (i=1:size(y)) x(i) = index(A,y(i))
!-----------------------------------------------------------------------------------------
end function indexes



!-----------------------------------------------------------------------------------------
function resize(A) result(B)
!! Function returns a more compact table list by resizing the array
!!
implicit none
integer,CONTIGUOUS , intent(in)	:: A(:,:)	   !! Table with too many 0, to be compressed 
integer 					:: B(maxval(A(1,:)),size(A,2))	!! Output resized table
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
B(:,:) = A(1:size(B,1),:)
!-----------------------------------------------------------------------------------------
end function resize



!-----------------------------------------------------------------------------------------
subroutine slice_1contact(g,E,S)
!! Procedure returns the slices of the graph going from 1 edge. BFS algorithm is used to 
!! assign a distance to the edge for each point, then the slice number = distance.
!!
!! @note S(1,:) indicates the end of the slice, so S(1,:)-1 is the number of points in a slice
implicit none
integer,CONTIGUOUS , intent(in)	:: g(:,:)		!! Graph connectivity table. 
integer,CONTIGUOUS , intent(in) :: E(:) 		!! Left Edge points' IDs
integer, allocatable,intent(out):: S(:,:)       !! Output the Slices, 2nd index is the Slice number 
integer   						:: D(size(g,2))	!! Distance table of the points
integer   						:: i, NP
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
D = dist(g,E) +1
! Add 1 in order to have distance=1 for the edge points
NP = -1
do i=1,maxval(D)
	if (COUNT(D == i) > NP) NP = COUNT(D == i) 
enddo
allocate(S(NP+1,maxval(D)))
S =1
do i=1,size(g, dim=2)
	S(1,D(i)) = S(1,D(i)) +1
	S(S(1,D(i)), D(i)) = i 
enddo
!-----------------------------------------------------------------------------------------
end subroutine slice_1contact


!-----------------------------------------------------------------------------------------
recursive subroutine slice_2contacts(g,E1,E2,NMAX,S)
!! Function returns the slices of the graph. The problem is solved in a divide-and-conquer  
!! manner. The algorithm is following.
!!
!! 1. Compute the distance of points to left and right edges, by calling [[dist]]
!! 2. Divide the graph into 2 parts based on if the right distance is larger
!! 3. Find the connections between 2 parts, and define new edges
!! 4. ``Slice`` seperately the 2 parts (recursive)
!! 5. Combine the results
!! 
!! Stop condition: the 2 edges touch each other, or the number of points in the part is 
!! small enough.
!!
!! @note S(1,:) indicates the end of the slice, so S(1,:)-1 is the number of points in a slice.
!!
implicit none
integer,CONTIGUOUS , intent(in) :: g(:,:)		                      !! Graph connectivity table. 
integer,CONTIGUOUS , intent(in) :: E1(:) 		                      !! Left Edge points' IDs
integer,CONTIGUOUS , intent(in) :: E2(:) 		                      !! RightEdge points' IDs
integer, intent(in)				:: NMAX			                      !! Maximum number of points in a single slice
integer, allocatable,intent(out):: S(:,:)                             !! Output the Slices, 2nd index is the Slice number 
integer   						:: P(size(g,2))	                      !! Partition table 
integer   						:: E(size(g,2),2)                     !! Edge table 
integer   						:: NP(2)		                      !! Number of points in parts
integer, allocatable			:: S1(:,:)		                      !! Slices in the left part  
integer, allocatable			:: S2(:,:)		                      !! Slices in the right part
integer, allocatable			:: PT(:)		                      !! Point list in a part
integer                         :: i, err
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
do i= 1,size(E1)
	if (ANY(E2(:) == E1(i))) then       		
	! If the 2 edges contact, than we need to put the whole part into one slice
		allocate(S(size(g,2)+1,1))
		S(1,1) = size(g,2)+1
		S(2:,1)= (/ (i,i=1,size(g,2)) /) 
		return
	endif 
enddo
! Split the Graph into 2 parts
P = part(dist(g,E1),dist(g,E2))
FORALL (i=1:2) NP(i) = count(P==i)
! Determinate the edges
E = edge(g,P)
! Slice each part
! If a part contains less points than NMAX, than all the points in this part goes into one slice 
allocate(PT(NP(1)), stat=err)
if (err /= 0) print *, "PT: Allocation request denied"
PT = PtsInParts(P,1)
if ((NP(1) .gt. NMAX).and.(NP(1) .gt. size(E1))) then
	CALL slice(subgraph(g,PT),indexes(PT,E1),indexes(PT,E(2:E(1,1),1)),NMAX,S1)
	! Recover the original point ID
	FORALL (i=2:size(S1,1))	
		where(S1(i,:) .ne. -1)  S1(i,:) = PT(S1(i,:))					
		! -1 to indicate no number
	end forall
else 
	allocate(S1(NP(1)+1,1))
	S1(1,1) = NP(1)+1
	S1(2:,1)= PT
endif
if (allocated(PT)) deallocate(PT, stat=err)
if (err /= 0) print *, "PT: Deallocation request denied"	
allocate(PT(NP(2)), stat=err)
if (err /= 0) print *, "PT: Allocation request denied"
PT = PtsInParts(P,2)
if ((NP(2) .gt. NMAX).and.(NP(2) .gt. size(E2))) then
	CALL slice(subgraph(g,PT),indexes(PT,E(2:E(1,2),2)),indexes(PT,E2),NMAX,S2)	
	! Recover the original point ID		
	FORALL (i=2:size(S2,1)) 
		where(S2(i,:) .ne. -1) S2(i,:) = PT(S2(i,:))					
		! -1 to indicate no number
	end forall
else 
	allocate(S2(NP(2)+1,1))
	S2(1,1)	= NP(2)+1
	S2(2:,1)= PT
endif
if (allocated(PT)) deallocate(PT, stat=err)
if (err /= 0) print *, "PT: Deallocation request denied"	
! Merge the results of 2 parts
allocate(S(max(size(S1,1),size(S2,1)),size(S1,2)+size(S2,2)))
S = -1 
! -1 to indicate no number
S(1:size(S1,1),1:size(S1,2))	= S1 
S(1:size(S2,1),1+size(S1,2):) 	= S2
if (allocated(S1)) deallocate(S1, stat=err)
if (err /= 0) print *, "S1: Deallocation request denied"	
if (allocated(S2)) deallocate(S2, stat=err)
if (err /= 0) print *, "S2: Deallocation request denied"	
!-----------------------------------------------------------------------------------------
end subroutine slice_2contacts




!-----------------------------------------------------------------------------------------
function testSlicing(g,S) result(b)
!! Function tests if a slicing from subroutine [[slice]]  of the graph is consistent by  
!! looking at the neighbors of all the points in each slice. 
!!
implicit none
integer,CONTIGUOUS , intent(in)	:: g(:,:)				!! Graph connectivity table. 
integer							:: b 					!! Test result
integer,CONTIGUOUS , intent(in)	:: S(:,:)               !! Slices, 2nd index is the Slice number 
integer 						:: L(size(g,2))			!! Slice number of points
integer 						:: i, k,j,err			
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
L = -1
L = SliceNum(S)
IF (ANY (L == -1) ) then 
	b = -1 							     	! Not all points belong to a slice 
	return 
endif 
do i=1,size(S,2)							! Loop over the slices
	do j=2,S(1,i)							! Loop over the points in one slice
		do k=2,g(1,S(j,i))					! Loop over all its neighbor points
			if (ABS(L(g(k,S(j,i))) - i) .gt. 1) then 
				b = i 		        		! It connects to a points too far
!				print *,L(g(k,S(j,i))),g(k,S(j,i)), i,S(j,i)
				return 
			endif 
		enddo 
	enddo 
enddo
b = 0
!-----------------------------------------------------------------------------------------
end Function testSlicing




!-----------------------------------------------------------------------------------------
function SliceNum(S) result(N)
!! Function returns an array of indexes of slice to which the points belong 
!!
implicit none
integer,CONTIGUOUS , intent(in)	:: S(:,:)                !! Slices, 2nd index is the Slice number 
integer 						:: N(1:maxval(S(2:,:)))	 !! Slice number of points
integer 						:: i,j
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
do i=1,size(S, 2)							! Loop over the slices
	do j=2,S(1,i)							! Loop over the points in one slice
		N(S(j,i)) = i
	enddo 
enddo
!-----------------------------------------------------------------------------------------
end Function SliceNum




end module graph
