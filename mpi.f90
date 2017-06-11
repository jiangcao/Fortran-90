module mpi

    implicit none
    
    private
    INTEGER, save  :: idnode, mxnode

    public initmpi, exitmpi , idnode, mxnode


    contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initialize openmpi 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE initmpi()

 implicit none
 include "mpif.h"

 INTEGER ierr

 call MPI_init(ierr)
 idnode = 0
 call MPI_COMM_RANK(MPI_COMM_WORLD, idnode ,ierr)
 mxnode = 1
 call MPI_COMM_SIZE(MPI_COMM_WORLD, mxnode, ierr)

 return

END SUBROUTINE initmpi
! Exit openmpi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE exitmpi

 implicit none
 include "mpif.h"
 INTEGER ierr

 call MPI_FINALIZE(ierr)
 return

 END SUBROUTINE exitmpi

end module mpi
