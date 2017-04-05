module matrix_c
!-----------------------------------------------------------------------------------------
!! Complex Matrix Library
!! A 2D complex array, element to form a list/table of complex matrices
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
use types, only : dp
!-----------------------------------------------------------------------------------------
implicit none

type type_matrix_complex
    complex(dp),allocatable :: m(:,:)
    !! complex matrix
    integer :: n(2)
    !! matrix size
end type type_matrix_complex


interface operator(.m.)
	module procedure array_times_array_simple
end interface

interface operator(**)
	module procedure array_power, &
					 array_transpose 
end interface

interface operator(.md.)
	module procedure array_times_array_dagger
end interface

interface operator(.dm.)
	module procedure array_dagger_times_array
end interface

interface sizeof 
	module procedure array_size,matrix_size, matrix_list_size, matrix_size_dim, matrix_list_size2
end interface sizeof

interface ReadTxt
	
end interface ReadTxt

interface SaveTxt 
	module procedure matrix_list_print,array_print
end interface SaveTxt

interface show
	module procedure array_print_on_screen
end interface show

interface alloc 
	module procedure matrix_list_allocElem2,matrix_list_allocElem
	module procedure matrix_alloc, matrix_alloc2, array_alloc, array_alloc2 
end interface

interface free 
	module procedure matrix_free
end interface

interface eye 
	module procedure array_eye 
end interface

interface eig
	module procedure  array_eigen
end interface

interface diag 
	module procedure array_to_diag
end interface

interface inv
	module procedure array_inverse, array_inverse2	
end interface

interface trace_c
	module procedure array_trace,matrix_trace
end interface

contains

!  =====  Allocation/Deallocation  =====
!-----------------------------------------------------------------------------------------
pure subroutine matrix_alloc(M,n,nn,source)	
implicit none 
type(type_matrix_complex), intent(out) :: M 
integer, intent(in)                    :: n
integer, intent(in), optional          :: nn
complex(dp), intent(in), optional       :: source(:,:)
if (present(nn)) then
	call matrix_alloc2(M, (/n,nn/), source=source)
else 
	call matrix_alloc2(M, (/n,n/), source=source)
endif
!-----------------------------------------------------------------------------------------
end subroutine matrix_alloc

!-----------------------------------------------------------------------------------------
pure subroutine matrix_alloc2(M,n,source)	
implicit none 
type(type_matrix_complex), intent(out) :: M 
integer, intent(in)                    :: n(2)
complex(dp), intent(in), optional       :: source(1:n(1),1:n(2))
if (.not. allocated(M%m)) then
	allocate(M%m(n(1),n(2)))
else 
	if ((M%n(1) == n(1)).and.(M%n(2) == n(2))) then 
	else
		deallocate(M%m)
		allocate(M%m(n(1),n(2)))
	endif
endif
if (present(source)) then
    M%m(:,:) = source(:,:)
else 
    M%m(:,:) = dcmplx(0.0d0,0.0d0)
endif
M%n  = n
!-----------------------------------------------------------------------------------------
end subroutine matrix_alloc2

!-----------------------------------------------------------------------------------------
pure subroutine array_alloc(M,n,nn)	
implicit none 
complex(dp), intent(out), allocatable   :: M (:,:)
integer, intent(in)                    :: n
integer, intent(in),optional           :: nn
if (present(nn)) then 
	call array_alloc2(M, (/n,nn/))
else 
	call array_alloc2(M, (/n,n/))
endif
!-----------------------------------------------------------------------------------------
end subroutine array_alloc

!-----------------------------------------------------------------------------------------
pure subroutine array_alloc2(M,n)	
!! This function allocates a 2D complex array. If the array is already allocated, this function
!! will resize the array to the new size. The allocated array is initiated to zero.
implicit none 
complex(dp), intent(out), allocatable   :: M (:,:)
integer, intent(in)                    :: n(2)
if (.not. allocated(M)) then
	allocate(M(n(1),n(2)))
else
	if ((size(M,1) == n(1)).and.(size(M,2) == n(2))) then 
	else 
		deallocate(M)
		allocate(M(n(1),n(2)))
	endif
endif
M  = dcmplx(0.0d0,0.0d0)
!-----------------------------------------------------------------------------------------
end subroutine array_alloc2

!-----------------------------------------------------------------------------------------
pure function array_eye(n) result (R)
implicit none
integer, intent(in) :: n
complex(dp)::R(n,n)
INTEGER   :: ii 
R = dcmplx(0.0d0,0.0d0)
forall (ii=1:n) R(ii,ii) = dcmplx(1.0d0,0.0d0)
!-----------------------------------------------------------------------------------------
end function array_eye

!-----------------------------------------------------------------------------------------
pure subroutine matrix_list_allocElem(this,nx,nm,nn,source) 
	implicit none
	integer, intent(in) 						:: nx
	integer, intent(in) 						:: nm(1:nx)
	integer, intent(in),optional    			:: nn(1:nx)
	type(type_matrix_complex),intent(out) 	    :: this (1:nx)
	complex(dp), intent(in), optional            :: source(:,:,:) !! the source data to put into the matrices
	integer :: ii 
	do ii=1,nx
		if (present(nn)) then 
			call matrix_alloc2(this(ii), (/nm(ii),nn(ii)/), source = source(:,:,ii))
		else 
			call matrix_alloc2(this(ii), (/nm(ii),nm(ii)/), source = source(:,:,ii))			
		endif
	enddo
!-----------------------------------------------------------------------------------------	
end subroutine matrix_list_allocElem


!-----------------------------------------------------------------------------------------
pure subroutine matrix_list_allocElem2(this,nx,n,source) 
	implicit none
	integer, intent(in) 						:: nx , n(2,1:nx)
	type(type_matrix_complex),intent(out) 	    :: this (1:nx)
	complex(dp), intent(in), optional            :: source(:,:,:) !! the source data to put into the matrices
	integer :: ii 
	do ii=1,nx
		call matrix_alloc2(this(ii), n(1:2,ii), source = source(:,:,ii))		
	enddo
!-----------------------------------------------------------------------------------------	
end subroutine matrix_list_allocElem2

!-----------------------------------------------------------------------------------------	
elemental subroutine matrix_free(this)
implicit none
type(type_matrix_complex),intent(out) :: this
 	if (allocated(this%m)) deallocate(this%m)
!-----------------------------------------------------------------------------------------	
end subroutine matrix_free


!  =====  size, print  =====
!-----------------------------------------------------------------------------------------	
pure function matrix_list_size(list,dim) result(nm)
	implicit none 
	type(type_matrix_complex), intent(in) :: list(:) 
	integer ,intent(in)					  :: dim
	INTEGER 							  :: nm(size(list))
	integer :: ii
	forall (ii=1:size(list)) nm(ii) = list(ii)%n(dim)
!-----------------------------------------------------------------------------------------	
end function matrix_list_size

!-----------------------------------------------------------------------------------------	
pure function matrix_list_size2(list) result(nm)
	implicit none 
	type(type_matrix_complex), intent(in) :: list(:) 
	INTEGER 							  :: nm(2,size(list))
	integer :: ii
	forall (ii=1:size(list)) nm(:,ii) = list(ii)%n(:)
!-----------------------------------------------------------------------------------------		
end function matrix_list_size2

!-----------------------------------------------------------------------------------------		
pure function array_size(this) result(s)
	implicit none
	complex(dp),intent(in) :: this(:,:)
	integer :: s(2), ii
	FORALL (ii=1:2) s(ii) = size(this,dim=ii)
!-----------------------------------------------------------------------------------------			
end function array_size

!-----------------------------------------------------------------------------------------		
pure function matrix_size(this) result(s)
	implicit none
	type(type_matrix_complex),intent(in) :: this 
	integer :: s(2)
	s(:) = this%n
!-----------------------------------------------------------------------------------------			
end function matrix_size

!-----------------------------------------------------------------------------------------			
pure function matrix_size_dim(this,dim) result(s)
	implicit none
	type(type_matrix_complex), intent(in) :: this 
	integer, intent(in) :: dim
	integer :: s
	s = this%n(dim)
!-----------------------------------------------------------------------------------------				
end function matrix_size_dim

!-----------------------------------------------------------------------------------------				
subroutine matrix_list_print(handle,this)
	implicit none
	type(type_matrix_complex), intent(in) 	:: this(:)
	integer, intent(in), optional 			:: handle
	integer :: ii,xx,yy
	if (present(handle)) then
		write(handle, '(A)') "Format XY"
		write(handle, '(1(i8))') size(this)
		write(handle, '(3(i8),es15.4,es15.4)') (((ii,xx,yy,this(ii)%m(xx,yy),xx=1,size(this(ii)%m,1)),yy=1,size(this(ii)%m,2)),ii=1,size(this))
	else
		print '(3(i8),es15.4,es15.4)',(((ii,xx,yy,this(ii)%m(xx,yy),xx=1,size(this(ii)%m,1)),yy=1,size(this(ii)%m,2)),ii=1,size(this))
	endif
!-----------------------------------------------------------------------------------------					
end subroutine matrix_list_print

!-----------------------------------------------------------------------------------------				
subroutine array_print(handle,A)
	implicit none
	complex(dp), intent(in) 					:: A(:,:)
	integer, intent(in)          			:: handle
	integer :: xx,yy
		write(handle, '(A)') "Format XY"
		write(handle, '(2(i8))') size(A,1),size(A,2)
		write(handle, '(2(i8),es15.4,es15.4)') ((xx,yy,A(xx,yy),xx=1,size(A,1)),yy=1,size(A,2))
		write(handle, '(A)') "END"
!-----------------------------------------------------------------------------------------						
end subroutine array_print

!-----------------------------------------------------------------------------------------				
subroutine matrix_read(handle,A)
	implicit none
	type(type_matrix_complex),intent(out)   :: A
	integer, intent(in)          			:: handle
	integer    :: xx,yy
	real(dp)    :: re, im
	character(len=100) :: s
		read(handle,'(100A)') s
		if (trim(s) == "Format XY" ) then 
			read(handle,*) xx,yy 
			call matrix_alloc(A,xx,yy)
			read(handle,'(100A)') s 
			do while (trim(s) /= "END")				
				read(s, *) xx,yy,re,im 
				A%m(xx,yy) = dcmplx(re,im)
				read(handle,'(100A)') s 				
			enddo			
 		else

		endif		
!-----------------------------------------------------------------------------------------						
end subroutine matrix_read

subroutine array_print_on_screen(A)
!-----------------------------------------------------------------------------------------						
	implicit none
	complex(dp), intent(in) 					:: A(:,:)
	integer :: xx,yy
    do xx = 1,size(A,1)
		print '(10(A,es8.1,",",es8.1,")"))',("(",A(xx,yy),yy=1,size(A,2))
    enddo
!-----------------------------------------------------------------------------------------						
end subroutine array_print_on_screen


! subroutine blocs_matrix_print(handle,A)
! 	implicit none 
! 	type(type_matrix_complex), intent(in)  :: A (:,:)
! 	integer, intent(in)                    :: handle
! 	integer :: sizes(2,0:size(A, dim=1),0:size(A, dim=2))
! 	integer :: i,j,n(size(A, dim=2)),m(size(A, dim=1)),x,y 
! 	sizes = 0
! 	n=0
! 	m=0
! 	sizes(:,1:,1:) = blocs_matrix_size(A)
! 	do i=1,size(A, dim=2)		
! 		sizes(2,0,i) = maxval(sizes(2,1:,i))
! 	enddo 
! 	do j=1,size(A, dim=1)		
! 		sizes(1,j,0) = maxval(sizes(1,j,1:))
! 	enddo 
! 	do i=1,size(A, dim=2)
! 		n(i) = sum(sizes(2,0,1:i-1))
! 	enddo		
! 	do j=1,size(A, dim=1)
! 		m(j) = sum(sizes(1,1:j-1,0))
! 	enddo
! 	do i=1,size(A, dim=2)
! 		do j=1,size(A, dim=1)
! 			if (allocated(A(j,i)%m)) then
! 				do x=1,size(A(j,i)%m, dim=2)
! 					do y=1,size(A(j,i)%m, dim=1)
! 						write(handle,'(2I8,2E15.5)') y+m(j),x+n(i),A(j,i)%m(y,x)
! 					enddo 
! 				enddo 
! 			endif
! 		enddo 
! 	enddo
! end subroutine blocs_matrix_print

!-----------------------------------------------------------------------------------------	
pure function array_testHermitian(M) result(b)
implicit none 
complex(dp), intent(in) :: M(:,:)
logical :: b 
integer :: i,j 
real(dp),parameter :: TOL = 1.0D-10
b = .true.
do i=1,size(M,2)
	do j=1,i 
		if(abs(M(i,j) - conjg(M(j,i))) .gt. TOL) then 
			b = .false.
			return 
		endif 
	enddo
enddo
!-----------------------------------------------------------------------------------------	
end function array_testHermitian

! pure function blocs_matrix_testHermitian(M) result(b)
! 	implicit none 
! 	type(type_matrix_complex) , intent(in) :: M(:,:)
! 	logical :: b 
! 	integer :: i,j,x,y 
! 	b = .true.
! 	do i = 1,size(M, dim=2)
! 		do j = 1,size(M, dim=1)
! 			if (allocated(M(j,i)%m)) then 
! 				do x = 1,size(M(j,i)%m,2)
! 					do y = 1,size(M(j,i)%m,1)
! 						if (.not. allocated(M(i,j)%m)) then 
! 							b = .false.
! 							return 
! 						endif 
! 						if (M(i,j)%m(x,y) /= conjg(M(j,i)%m(y,x))) then 
! 							b = .false.
! 							return 
! 						endif 
! 					enddo 
! 				enddo 
! 			endif 
! 		enddo 
! 	enddo
! end function blocs_matrix_testHermitian




!  =====  Multiplications and others  =====
!-----------------------------------------------------------------------------------------	
pure function array_times_array_dagger(A,B) result(C)
implicit none 
	complex(dp)			, intent(in)			:: A(:,:),B(:,:)
	complex(dp)									:: C(size(A,1),size(B,1))
	C = array_times_array(A,B,trA=.false.,trB=.true.,cjA=.false.,cjB=.true.)
!-----------------------------------------------------------------------------------------	
end function array_times_array_dagger

!-----------------------------------------------------------------------------------------	
pure function array_dagger_times_array(A,B) result(C)
implicit none 
	complex(dp)			, intent(in)			:: A(:,:),B(:,:)
	complex(dp)									:: C(size(A,2),size(B,2))
	C = array_times_array(A,B,trA=.true.,trB=.false.,cjA=.true.,cjB=.false.)
!-----------------------------------------------------------------------------------------		
end function array_dagger_times_array

!-----------------------------------------------------------------------------------------	
pure function array_times_array(A,B,trA,trB,cjA,cjB) result(C)
	implicit none 
	complex(dp)			, intent(in)			:: A(:,:),B(:,:)
	LOGICAL(KIND=4)		, intent(in)			:: trA,trB
	complex(dp)									:: C(size(A,merge(2,1,trA)),size(B,merge(1,2,trB)))	
	LOGICAL(KIND=4)		, intent(in), optional	:: cjA,cjB
	integer :: lda,ldb,k,m,kb,n
	character :: ctrA,ctrB
	interface 
		pure subroutine ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
			COMPLEX(dp),intent(in) 			:: ALPHA,BETA
    	    INTEGER,intent(in) 				:: K,LDA,LDB,LDC,M,N
     		CHARACTER,intent(in) 			:: TRANSA,TRANSB
     		COMPLEX(dp),intent(in) 			:: A(lda,*),B(ldb,*)
     		COMPLEX(dp),intent(inout) 		:: C(ldc,*)
		end subroutine ZGEMM
	end interface
	lda = size(A,1)
	ldb = size(B,1)
	if (.not. trA) then 
		k = size(A,2)
		m = size(A,1)
		ctrA = 'n'
	else 		
		k = size(A,1)
		m = size(A,2)
		if (present(cjA) .and. cjA) then 
			ctrA = 'c'
		else 
			ctrA = 't'
		endif
	endif
	if (.not. trB)then
		kb = size(B,1)
		n = size(B,2)
		ctrB = 'n'
	else
		kb = size(B,2)
		n = size(B,1)
		if (present(cjB) .and. cjB) then 
			ctrB = 'c'
		else 
			ctrB = 't'
		endif		
	endif
	! if (k /= kb) then 
	! 	print *,'ERROR in array_times_array:','k=',k,'kb=',kb
	! 	call abort
	! endif
	call zgemm(ctrA,ctrB,m,n,k,dcmplx(1.0d0,0.0d0),A,lda,B,ldb,dcmplx(0.0d0,0.0d0),C,m)
!-----------------------------------------------------------------------------------------		
end function array_times_array

!-----------------------------------------------------------------------------------------	
pure function array_times_array_simple(A,B) result(C)
	implicit none 
	complex(dp),intent(in)			        :: A(:,:),B(:,:) 
	complex(dp)								:: C(size(A,1),size(B,2))
	C = array_times_array(A,B,.false.,.false.)
!-----------------------------------------------------------------------------------------		
end function array_times_array_simple

!-----------------------------------------------------------------------------------------	
function array_power(A,n) result (C)
	implicit none
	complex(dp), intent(in) 				  :: A(:,:)
	integer , intent(in)				  :: n 
	complex(dp)				              :: B(size(A,1),size(A,1))
	complex(dp)							  :: C(size(A,1),size(A,1))
	integer :: ii
	if (n > 0) then		
		B = A
		do ii = 2, n
			B = B .m. A
		enddo
		C = B
	elseif (n == 0) then 
		C = array_eye(size(A, dim=1))
	elseif (n == -1) then
		C = array_inverse(A)
	else
		C = array_inverse(A)
		B = C
		do ii = 2, -n
			B = B .m. C
		enddo
		C = B
	endif
!-----------------------------------------------------------------------------------------		
end function array_power

!-----------------------------------------------------------------------------------------	
pure function array_transpose(A,t) result (C)
	implicit none
	complex(dp)	, intent(in) 			  :: A(:,:)
	character , intent(in)				  :: t
	complex(dp)							  :: C(size(A,2),size(A,1))
	if ((t=='t').or.(t=='T')) then 
		C = Transpose(A)
	elseif((t=='c').or.(t=='C')) then 
		C = Transpose(Conjg(A))
	endif
!-----------------------------------------------------------------------------------------		
end function array_transpose

!-----------------------------------------------------------------------------------------	
function array_eigen(A,B,eigvec,itype,uplo) result(eig)
	implicit none 
	complex(dp), intent(in)				:: A(:,:)
	complex(dp), intent(in), optional	:: B(:,:) 
	real(dp)								:: eig(size(A,1))	
	complex(dp), intent(inout),optional 	:: eigvec(size(A,1),size(A,2)) 
	integer , intent(in), optional 		:: itype
	CHARACTER,intent(in), optional 		:: uplo
	integer 		:: LDA, N, LDB, lwork, INFO, itypeop,ii
	CHARACTER 		:: jobz,uploop
	real(dp)		  	:: W(size(A,1)), RWORK(3*size(A,2))
	complex(dp)    	:: work( 1 + 4*size(A,2) + size(A,2)**2), C(size(A,1),size(A,2))
	C(:,:) = A(:,:)
	if (present(eigvec)) then 
		jobz = 'V'
	else
		jobz = 'N'
	endif
	uploop 	= merge(uplo,'U',present(uplo))
	itypeop = merge(itype, 1, present(itype))	
	N 	= size(A, dim=2)
	LDA = size(A, dim=1)
	LWORK = size(WORK)
	if (present(B)) then 
		LDB = size(B, dim=1)
		call zhegv(itypeop,jobz,uploop,N,C,LDA,B,LDB,eig,WORK,LWORK,RWORK,INFO)
		if (INFO .ne. 0) then 
			print *,'@array_eigen ZHEGV fails with INFO=',INFO 
		call abort
	endif
	else
		LDB = LDA
		call zheev(jobz,uploop,N,C,LDA,eig,WORK,LWORK,RWORK,INFO)
		if (INFO .ne. 0) then 
			print *,'@array_eigen ZHEEV fails with INFO=',INFO 
		call abort
	endif
	endif
	if (present(eigvec)) eigvec = C 
!-----------------------------------------------------------------------------------------		
end function array_eigen

!-----------------------------------------------------------------------------------------	
pure function array_to_diag(A) result(diag)
implicit none
complex(dp), intent(in) 		:: A(:,:)
complex(dp)					:: diag(size(A,1))
integer :: ii 
forall (ii=1:size(A,1)) diag(ii) = A(ii,ii)
!-----------------------------------------------------------------------------------------	
end function array_to_diag

!-----------------------------------------------------------------------------------------	
function array_inverse2(A, UPLO) result(C) ! for Hermitian matrix
implicit none
complex(dp), intent(in)		:: A(:,:)
complex(dp)					:: C(size(A, dim=1),size(A, dim=1))
CHARACTER, intent(in) 		:: UPLO 
integer 	:: info,lda,lwork,n,nnz
integer 	:: ipiv(size(A,1))
complex(dp), allocatable 	:: work(:,:)
 n = size(A,1)
 if (n /= size(A,2)) then
 	print *,'@array_inverse, size not square',n,size(A,2)
 	call abort
 endif
 C(:,:) = A(:,:)
 allocate(work(n*n,n*n))
 LDA = size(A,2)
 call zhetrf(UPLO, n, C, LDA, ipiv, WORK, size(WORK), info)
 if (info .ne. 0) then 
 	print *,'@array_inverse2 ZHETRF fails with INFO=', info
 	call abort
 endif
 call zhetri(UPLO, n, C, LDA, ipiv, WORK, info)
 if (info .ne. 0) then 
 	print *,'@array_inverse2 ZHETRI fails with INFO=', info
 	call abort
 endif
!-----------------------------------------------------------------------------------------	
end function array_inverse2

!-----------------------------------------------------------------------------------------	
function array_inverse(A) result(C) ! for General matrix
implicit none
complex(dp), intent(in)		:: A(:,:)
complex(dp)					:: C(size(A, dim=1),size(A, dim=1))
integer 	:: info,lda,lwork,n,nnz
integer 	:: ipiv(size(A,1))
complex(dp), allocatable 	:: work(:,:)
 n = size(A,1)
 if (n /= size(A,2)) then
 	print *,'@array_inverse, size not square',n,size(A,2)
 	call abort
 endif
 C(:,:) = A(:,:)
 allocate(work(n*n,n*n))
 call zgetrf(n,n,C,n,ipiv,info)
 if (info .ne. 0) then 
 	print *,'@array_inverse ZGETRF fails with INFO=', info
 	call abort
 endif
 call zgetri(n,C,n,ipiv,work,n*n,info)
 if (info .ne. 0) then 
 	print *,'@array_inverse ZGETRI fails with INFO=', info
 	call abort
 endif
!-----------------------------------------------------------------------------------------	 
end function array_inverse

!-----------------------------------------------------------------------------------------	
pure function array_trace(A) result (tr)
implicit none
complex(dp), intent(in)		:: A(:,:)
complex(dp)					:: tr 
integer :: ii
tr = sum((/(A(ii,ii),ii=1,size(A, 1))/))
!-----------------------------------------------------------------------------------------	
end function array_trace

!-----------------------------------------------------------------------------------------	
elemental function matrix_trace(M) result (tr)
implicit none
type(type_matrix_complex), intent(in) :: M
complex(dp)					:: tr 
integer :: ii
tr = sum((/(M%m(ii,ii),ii=1,size(M%m, 1))/))
!-----------------------------------------------------------------------------------------	
end function matrix_trace

!-----------------------------------------------------------------------------------------	
subroutine matrix_copy(matrices, tab) 
	implicit none
	type(type_matrix_complex), intent(in) :: matrices(:) 
	complex(dp), intent(out)               :: tab(:,:,:)
!-----------------------------------------------------------------------------------------	
integer :: i 
forall (i=1:size(matrices)) tab(1:matrices(i)%n(1),1:matrices(i)%n(2),i) = matrices(i)%m(:,:)
!-----------------------------------------------------------------------------------------	
end subroutine matrix_copy	

end module matrix_c
