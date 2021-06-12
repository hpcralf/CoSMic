! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
! Made F conformant by Walt Brainerd
module qsort_c_module
  implicit none
  public :: QsortC
  public :: BinarySearch
  private :: Partition

contains
! this code may not suitable for large problem for the reason of memory limitation.
  recursive subroutine QsortC(A)
    integer, intent(in out), dimension(:) :: A
    integer :: iq
    if(size(A) > 1) then
      call Partition(A, iq)
      call QsortC(A(:iq-1))
      call QsortC(A(iq:))
    endif
  end subroutine QsortC

  subroutine Partition(A, marker)
    integer, intent(in out), dimension(:) :: A
    integer, intent(out) :: marker
    integer :: i, j
    integer :: temp
    integer :: x      ! pivot point
    x = A(1)
    i= 0
    j= size(A) + 1
    do
      j = j-1
      do
        if (A(j) <= x) exit
        j = j-1
      end do
      i = i+1
      do
        if (A(i) >= x) exit
        i = i+1
      end do
      if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
      elseif (i == j) then
        marker = i+1
        return
      else
        marker = i
        return
      endif
    end do
  end subroutine Partition

  integer function BinarySearch(elem,array)
    implicit none
    integer, dimension(:)            :: array
    integer                          :: elem,size_in
    integer                          :: L,R,M
    size_in = size(array)
    L = 1
    R = size_in
    M = (L+R)/2
    if ((elem < array(L) .OR. (elem > array(R))))then
        BinarySearch = 0
        return
    end if
!     if (elem == array(L)) then
!         BinarySearch = L
!         return
!     end if
!     if (elem == array(R)) then
!         BinarySearch = R
!         return
!     end if
    
    do while (L <= R)
        if (elem > array(M)) then
            L = M + 1
            M = (L+R)/2
        else if (elem < array(M)) then
            R = M - 1
            M = (L+R)/2
        else if (elem .EQ. array(M)) then
            BinarySearch = M
            return
        endif
    enddo
    BinarySearch = 0
    return
  end function
end module qsort_c_module
