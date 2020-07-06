module mpitranspose

  implicit none

contains

  subroutine mpi_transpose(data, rank, dims, mpicomm, mpierr, workarray)


    !
    ! This routine implements an optimised MPI tranpose, while
    ! leaving the fortran array the same size and dimensions.
    !
    ! The motivation for this routine is to allow Fortran to
    ! use FFTW, which can leave arrays in transpose form.
    !
    ! The user needs to ensure that the array sizes for the
    ! dimensions to be transposed are of the correct size
    ! (ie. multiples of ncpu's to be used).
    !

    include 'mpif.h'

    integer, intent(in) :: mpicomm, rank
    integer, intent(inout) :: mpierr
    integer, dimension(rank), intent(inout) :: dims
    complex(8), dimension(product(dims)), intent(inout) :: data
    complex(8), dimension(size(data)), optional, target, intent(inout) :: workarray

    complex(8), dimension(:), pointer :: work
    integer, dimension(3) :: edims
    integer :: mpisize, mpirank, chunksize, i, j, k, l, m, comm

    
    !
    ! Sanity check
    !
    if (rank < 2) then
       mpierr = -1
       return
    end if

    !
    ! Setup MPI
    !
    call MPI_Comm_dup(mpicomm, comm, mpierr)
    if (mpierr /= 0) return

    call MPI_Comm_rank(comm, mpirank, mpierr)
    if (mpierr /= 0) return

    call MPI_Comm_size(comm, mpisize, mpierr)
    if (mpierr /= 0) return

    !
    ! Sanity check
    !
    if (modulo(dims(rank-1), mpisize) /= 0) then
       mpierr = -2
       return
    end if

    !
    ! Setup work array
    !
    if (present(workarray)) then
       work => workarray
    else
       allocate(work(size(data)))
    end if

    !
    ! Setup effective dimensions
    !
    if (rank == 2) then
       edims(1) = 1
       edims(2) = dims(1)
       edims(3) = dims(2)
    else
       edims(1) = product(dims(:rank-2))
       edims(2) = dims(rank-1)
       edims(3) = dims(rank)
    end if

    !
    ! On processor transpose
    !
    l = 0
    do k = 1, edims(3)
       do j = 1, edims(2)
          do i = 1, edims(1)
             l = l + 1
             work(i+(k-1)*edims(1)+(j-1)*edims(1)*edims(3)) = data(l)
          end do
       end do
    end do

    !
    ! Communication
    !
    chunksize = size(data)/mpisize
    call mpi_alltoall(work, chunksize, MPI_DOUBLE_COMPLEX, data, chunksize, MPI_DOUBLE_COMPLEX, comm, mpierr)
    if (mpierr /= 0) return

    !
    ! On processor block transpose
    !
    l = 0
    do k = 1, edims(2)/mpisize
       do m = 1, mpisize
          do j = 1, edims(3)
             do i = 1, edims(1)
                l = l + 1
                work(l) = data(i+(j-1)*edims(1)+(m-1)*chunksize+(k-1)*edims(1)*edims(3))
             end do
          end do
       end do
    end do

    !
    ! Prepare to leave subroutine
    !
    data = work
    dims(rank-1) = edims(3)*mpisize
    dims(rank) = edims(2)/mpisize

    !
    ! Cleanup work array
    !
    if (present(workarray)) then
       nullify(work)
    else
       deallocate(work)
    end if

    !
    ! Cleanup MPI
    !
    call MPI_Comm_free(comm, mpierr)
    if (mpierr /= 0) return

  end subroutine mpi_transpose

end module mpitranspose
