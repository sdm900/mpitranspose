subroutine output(rank, dims, data, tdata)

  integer, intent(in) :: rank
  integer, dimension(rank), intent(in) :: dims
  real(8), dimension(product(dims)), intent(in) :: data, tdata

  integer :: i, stride
  character(80) :: tbuf
  
  stride = dims(1)

  do i = 1, size(data), stride
     write(tbuf, '(i)') stride
     tbuf = '('//trim(adjustl(tbuf))//'f7.1)'
     write(*, tbuf, advance="NO") data(i:i+stride-1)
     write(*, '(a)', advance="NO") "  ==  "
     write(*, tbuf) tdata(i:i+stride-1)
  end do

end subroutine output

#ifndef NDIM
#define NDIM 2
#endif

#if NDIM == 2
#define N1 6
#define N2 2
#define N N1,N2
#define TN N2,N1
#define PRODN N1*N2
#define RN1 N1
#define RN2 N2*mpisize
#define RN RN1,RN2
#define TRN RN2,RN1
#define PRODRN RN1*RN2
#define RDIMS :,:
#define DO(v1_,l1_,u1_,v2_,l2_,u2_,v3_,l3_,u3_,v4_,l4_,u4_) do v3_=l3_,u3_; do v4_=l4_,u4_
#define ENDDO end do; end do
#define VAR(v_,i1_,i2_,i3_,i4_) v_(i1_,i2_)
#define TVAR(v_,i1_,i2_,i3_,i4_) v_(i2_,i1_)
#endif

#if NDIM == 3
#define N1 6
#define N2 6
#define N3 2
#define N N1,N2,N3
#define TN N1, N3, N2
#define PRODN N1*N2*N3
#define RN1 N1
#define RN2 N2
#define RN3 N3*mpisize
#define RN RN1,RN2,RN3
#define TRN RN1,RN3,RN2
#define PRODRN RN1*RN2*RN3
#define RDIMS :,:,:
#define DO(v1_,l1_,u1_,v2_,l2_,u2_,v3_,l3_,u3_,v4_,l4_,u4_) do v2_=l2_,u2_; do v3_=l3_,u3_; do v4_=l4_,u4_
#define ENDDO end do; end do; end do
#define VAR(v_,i1_,i2_,i3_,i4_) v_(i1_,i2_,i3_)
#define TVAR(v_,i1_,i2_,i3_,i4_) v_(i1_,i3_,i2_)
#endif

#if NDIM == 4
#define N1 7
#define N2 5
#define N3 6
#define N4 2
#define N N1,N2,N3,N4
#define TN N1,N2,N3,N4
#define PRODN N1*N2*N3*N4
#define RN1 N1
#define RN2 N2
#define RN3 N3
#define RN4 N4*mpisize
#define RN RN1,RN2,RN3,RN4
#define TRN RN1,RN2,RN4,RN3
#define PRODRN RN1*RN2*RN3*RN4
#define RDIMS :,:,:,:
#define DO(v1_,l1_,u1_,v2_,l2_,u2_,v3_,l3_,u3_,v4_,l4_,u4_) do v1_=l1_,u1_; do v2_=l2_,u2_; do v3_=l3_,u3_; do v4_=l4_,u4_
#define ENDDO end do; end do; end do; end do
#define VAR(v_,i1_,i2_,i3_,i4_) v_(i1_,i2_,i3_,i4_)
#define TVAR(v_,i1_,i2_,i3_,i4_) v_(i1_,i2_,i4_,i3_)
#endif

program test

  use mpitranspose

  implicit none

  include 'mpif.h'

  integer, dimension(NDIM) :: dims
  complex(8), dimension(N) :: a, work
  complex(8), dimension(RDIMS), allocatable :: t, ra, rt
  integer :: i, j, k, l, ll, mpirank, mpisize, mpierr
  integer(8) :: plan

  dims = (/ N /)

  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierr)

  ll = mpirank*PRODN
  DO(l,1,N4,k,1,N3,j,1,N2,i,1,N1)
     ll = ll + 1
     VAR(a, i, j, k, l) = dble(ll)
  ENDDO

  !
  ! Do a dodgy single cpu transpose
  !
  allocate(ra(RN), rt(TRN), t(TN))
  call mpi_gather(a, size(a), MPI_DOUBLE_COMPLEX, ra, size(a), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, mpierr)

  DO(l,1,RN4,k,1,RN3,j,1,RN2,i,1,RN1)
     TVAR(rt, i, j, k, l) = VAR(ra, i, j, k, l)
  ENDDO

  call mpi_scatter(rt, size(t), MPI_DOUBLE_COMPLEX, t, size(t), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, mpierr)
  

!!$  !
!!$  ! Call my transpose
!!$  !
!!$  call mpi_transpose(a, NDIM, dims, MPI_COMM_WORLD, mpierr)
!!$  
!!$
!!$  if (mpirank == 0)  write(*,*) "==== ", dims, " ===="
!!$
!!$  do i = 0, mpisize-1
!!$     if (mpirank == i) then
!!$        call output(NDIM, dims, real(t), real(a))
!!$     end if
!!$
!!$     call mpi_barrier(MPI_COMM_WORLD, mpierr)
!!$  end do

  !
  ! Call my transpose
  !
  call mpi_transpose(a, NDIM, dims, MPI_COMM_WORLD, mpierr, work)

  if (mpirank == 0)  write(*,*) "==== ", dims, " ===="

  do i = 0, mpisize-1
     if (mpirank == i) then
        call output(NDIM, dims, real(t), real(a))
     end if

     call mpi_barrier(MPI_COMM_WORLD, mpierr)
  end do

  call MPI_Finalize(mpierr)

end program test
