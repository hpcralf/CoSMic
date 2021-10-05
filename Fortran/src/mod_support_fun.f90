!###############################################################################
!###############################################################################
!#      ___      __         _      
!#     / __\___ / _\  /\/\ (_) ___ 
!#    / /  / _ \\ \  /    \| |/ __|
!#   / /__| (_) |\ \/ /\/\ \ | (__ 
!#   \____/\___/\__/\/    \/_|\___|
!#
!#  COVID-19 Spatial Microsimulation  ---  For Germany  ########################
!###############################################################################
!#
!# Authors:      Qifeng Pan
!#               Ralf Schneider
!#               Christian Dudel
!#               Matthias Rosenbaum-Feldbruegge
!#               Sebastian Kluesener
!#
!# Contact:      ralf.schneider@hlrs.de
!#               qifeng.pan@hlrs.de
!#
!#==============================================================================
!#
!# Copyright (C) 2021 High Performance Computing Center Stuttgart (HLRS),
!#                    Federal Institute for Population Research (BIB),
!#                    The Max Planck Institute for Demographic Research (MPIDR)
!# 
!# This program is free software: you can redistribute it and/or modify
!# it under the terms of the GNU General Public License as published by
!# the Free Software Foundation, either version 3 of the License, or
!# (at your option) any later version.
!# 
!# This program is distributed in the hope that it will be useful,
!# but WITHOUT ANY WARRANTY; without even the implied warranty of
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!# GNU General Public License for more details.
!# 
!# You should have received a copy of the GNU General Public License
!# along with this program.  If not, see <http://www.gnu.org/licenses/>.
!#
!###############################################################################
!
! Module with support functions
!
!###############################################################################
module support_fun

  use qsort_c_module
  use quicksort_nr
  
  implicit none
  
interface get_unique
  module procedure get_unique
end interface

interface sum_bygroup
  module procedure sum_bygroup_mod
  module procedure sum_bygroup_array
end interface

interface get_index
   module procedure get_index_single_char
   module procedure get_index_single_integer
   module procedure get_index_mul_char
!    module procedure get_index_mul_char_and
   module procedure get_index_mul_integer
end interface

interface get_index_not
   module procedure get_index_not_single_integer
   module procedure get_index_not_mul_integer
   module procedure get_index_not_single_real
end interface

interface generate_seq
  module procedure seq_date
  module procedure seq_integer
!   module procedure seq_char
!   module procedure seq_real
end interface

interface sample
  module procedure sample_i1
  module procedure sample_i4
end interface sample

interface condition_and
  module procedure condition_and_integer
  module procedure condition_and_real

end interface
! interface get_index_invers
!   module procedure get_index_single_int_inv()
! end interface


contains

  function get_unique(array_in) result(array_out)

    implicit none
    
    integer             :: size_in,j,index
    integer,dimension(:):: array_in
    integer,dimension(size(array_in))   :: temp_in,idx_in
    integer,allocatable :: array_out(:)
    
    if (size(array_in) == 0)then
       allocate(array_out(0))
       return
    end if
    
    if (size(array_in) == 1)then
       array_out = (/1/)
       return
    end if
    

    index = 1
    size_in = size(array_in)
    ! call QsortC(temp_in)
    call SORTRX(size(array_in),array_in,idx_in)
    temp_in = array_in(idx_in)
    !detemine dimension of array_out
    do j=2,size_in
       if (temp_in(j)/=temp_in(j-1)) then
          index = index + 1
       endif
    enddo
    !allocate array and reset index number
    allocate(array_out(index))
    index = 1
    array_out(index) = temp_in(index)
    do j=2,size_in
       if (temp_in(j)/=temp_in(j-1)) then
          index = index + 1
          array_out(index) = temp_in(j)
       endif
    enddo
    
  end function get_unique
  
function get_index_single_char(array_in,tar)result(array_index)
implicit none
integer size_in,index,i
character(len=*), dimension(:)                  :: array_in
character(len=*)                                :: tar
integer,allocatable,dimension(:)                :: array_index
!detemine dimension of array_out
index=0
size_in = size(array_in)
do i=1,size_in
    if (trim(tar) == trim(array_in(i))) then
        index = index + 1
    endif
enddo
!allocate array and reset index number
allocate(array_index(index))
index=0
do i=1,size_in
    if (trim(tar) == trim(array_in(i))) then
        index = index + 1
        array_index(index) = i
    endif
enddo
end function get_index_single_char

function get_index_single_integer(array_in,tar) result(array_index)
implicit none
integer size_in,index,i
integer, dimension(:)                        :: array_in
integer                                      :: tar
integer,allocatable,dimension(:)             :: array_index
!detemine dimension of array_out
index=0
size_in = size(array_in)
do i=1,size_in
    if (tar == array_in(i)) then
        index = index + 1
    endif
enddo
allocate(array_index(index))
!allocate array and reset index number
index=0
do i=1,size_in
    if (tar == array_in(i)) then
        index = index + 1
        array_index(index) = i
    endif
enddo
end function get_index_single_integer

! function get_index_single_integer(array_in,tar) result(array_index)
!     use omp_lib
!     implicit none
!     integer,dimension(:)               :: array_in
!     integer                            :: tar
!     integer,allocatable                :: array_index(:)
!     integer,dimension(size(array_in))  :: array_temp
!     integer,allocatable                :: index(:)
!     integer :: len, from, middle, to_e, nthreads, thread, chunk, chunk2, i,from_temp,to_temp
!     len      = size(array_in)
!     !$OMP parallel
!     nthreads = omp_get_num_threads()
!     !$OMP end parallel
!     chunk    = len / nthreads  
!     allocate(index(nthreads))
! 
!     !$OMP parallel do default(none)          &
!     !$OMP firstprivate(chunk, len, nthreads) &
!     !$OMP shared(array_temp,array_in,index,tar) private(from, to_e, thread)
!     do thread = 0, nthreads-1
!         
!        from = thread           * chunk + 1
!        to_e   = min((thread + 1) * chunk, len)
!        call apply_search(array_in(from:to_e), array_temp(from:to_e),tar,index(thread+1))
! !        order(from:to) = order(from:to) + from - 1
!     end do
!     !$OMP end parallel do
! 
!     allocate(array_index(sum(index)))
!     do thread =0,nthreads-1
!         from = sum(index(1:thread))+1
!         to_e   = sum(index(1:thread+1))
!         from_temp = thread           * chunk + 1
!         to_temp   = from_temp+index(thread+1) -1
!         array_index(from:to_e) = array_temp(from_temp:to_temp)
!     end do
! 
! end function

subroutine apply_search(array_in,array_temp,tar,index)
    integer,dimension(:)            :: array_in,array_temp
    integer                         :: tar,index,i
    index = 0
    do i = 1,size(array_in)
        if (array_in(i)==tar)then
            index=index+1
            array_temp(index) = i
        endif
    end do
end subroutine

! function get_index_single_integer(array_in,tar) result(array_index)
! use qsort_c_module
! implicit none
! integer size_in,index,i
! integer, dimension(:)                        :: array_in
! integer                                      :: tar,LOC,minus,plus
! integer,allocatable,dimension(:)             :: array_index,seq
! integer,dimension(size(array_in))            :: temp_in,idx
! size_in = size(array_in)
! call SORTRX(size(array_in),array_in,idx)
! temp_in = array_in(idx)
! 
! LOC     = BinarySearch(tar,temp_in)
! minus   = 0
! plus    = 0
! do while (((LOC-minus) >= 0) .and. (temp_in(LOC-minus) == tar))
!     minus = minus + 1
! end do
! 
! do while((LOC+plus)<=size_in .and. temp_in(LOC+plus) == tar )
!     plus  = plus + 1 
! end do
! 
! allocate(seq(plus+minus-1))
! do i = 1,size(seq)
!     seq(i) = LOC-minus+i
! end do
! 
! array_index = idx(seq)
! 
! end function get_index_single_integer


! function get_index_mul_integer(array_in,array_tar) result(array_index)
! use qsort_c_module
! implicit none
! integer size_in,index,i,j,total
! integer, dimension(:)                          :: array_in
! integer, dimension(:)                          :: array_tar
! integer, allocatable, dimension(:)             :: array_index
! integer                                        :: LOC,size1,size2
! integer, dimension(size(array_in))             :: idx_in,temp_in
! integer, dimension(size(array_tar))            :: idx_tar,temp_tar
! 
! size1 = size(array_in)
! size2 = size(array_tar)
! 
! call SORTRX(size1,array_in,idx_in)
! call SORTRX(size2,array_tar,idx_tar)
! 
! temp_in  = array_in(idx_in)
! temp_tar = array_tar(idx_tar)
! if (size2 < 10000)then
!     i = 1
!     j = 1
!     total = 0
!     do while((i <= size(temp_in)) .and. (j <= size(temp_tar)))
! 
!         if (temp_in(i) < temp_tar(j))then
!             i = i + 1 
!         else if (temp_in(i) > temp_tar(j))then
!             j = j + 1
!         else
!             total = total + 1
!             i = i + 1
!         endif
!     end do
!     
!     allocate(array_index(total))
!     i = 1
!     j = 1
!     total = 0
!     do while((i <= size(temp_in)) .and. (j <= size(temp_tar)))
!         
!         if (temp_in(i) < temp_tar(j))then
!             i = i + 1 
!         elseif (temp_in(i) > temp_tar(j))then
!             j = j + 1
!         else
!             total = total + 1
!             array_index(total) = idx_in(i)
!             i = i + 1
!         endif
!     end do
! else
!     i = 1
!     j = 1
!     total = 0
!     do while((i <= size(temp_tar)) .and. (j <= size(temp_in)))
! 
!         if (temp_tar(i) < temp_in(j))then
!             i = i + 1 
!         else if (temp_tar(i) > temp_in(j))then
!             j = j + 1
!         else
!             total = total + 1
!             i = i + 1
!         endif
!     end do
!     
!     allocate(array_index(total))
!     i = 1
!     j = 1
!     total = 0
!     do while((i <= size(temp_tar)) .and. (j <= size(temp_in)))
!         
!         if (temp_tar(i) < temp_in(j))then
!             i = i + 1 
!         elseif (temp_tar(i) > temp_in(j))then
!             j = j + 1
!         else
!             total = total + 1
!             array_index(total) = idx_in(j)
!             i = i + 1
!         endif
!     end do
! end if
! end function  

! function get_index_mul_integer(array_in,array_tar) result(array_index)
! use qsort_c_module
! implicit none
! integer size_in,index,i,size_tar
! integer, dimension(:)                          :: array_in
! integer, dimension(:)                          :: array_tar
! integer, allocatable, dimension(:)             :: array_index
! integer                                        :: LOC
! integer, dimension(size(array_in))             :: idx_in,temp_in
! integer, dimension(size(array_tar))            :: idx_tar,temp_tar
! ! sort for array_in in order apply binarysearch
! 
! 
!  
! !detemine dimension of array_out
! index = 0
! size_in = size(array_in)
! if (size(array_tar)<1000)then
! call SORTRX(size(array_tar),array_tar,idx_tar)
! temp_tar = array_tar(idx_tar)
!     do i = 1,size_in
! 
!         LOC = BinarySearch(array_in(i),temp_tar)
!         if (LOC /= 0) then
!             index = index + 1
!         endif
!     enddo
!     !allocate array and reset index number
!     allocate(array_index(index))
!     index = 0
!     do i = 1,size_in
!         LOC = BinarySearch(array_in(i),temp_tar)
!         if (LOC /= 0) then
!             index = index + 1
!             array_index(index) = idx_tar(i)
!         endif
!     enddo
! else    
! call SORTRX(size(array_in),array_in,idx_in)
! temp_in = array_in(idx_in)
!     do i = 1,size(array_tar)
! 
!         LOC = BinarySearch(array_tar(i),temp_in)
!         if (LOC /= 0) then
!             index = index + 1
!         endif
!     enddo
!     !allocate array and reset index number
!     allocate(array_index(index))
!     index = 0
!     do i = 1,size_in
!         LOC = BinarySearch(array_in(i),temp_tar)
!         if (LOC /= 0) then
!             index = index + 1
!             array_index(index) = idx_in(i)
!         endif
!     enddo
! end if
! end function  

function get_index_mul_integer(array_in,array_tar) result(array_index)
implicit none
integer size_in,index,i,j,size_tar
integer, dimension(:)                          :: array_in
integer, dimension(:)                          :: array_tar
integer, allocatable, dimension(:)             :: array_index
integer                                        :: tar,temp_int
temp_int = 0
do i = 1,size(array_tar)
    tar = array_tar(i)
    do j = 1,size(array_in)
        if (array_in(j) == tar)then
            temp_int = temp_int + 1
        end if
    end do
end do
allocate(array_index(temp_int))
temp_int = 0

do i = 1,size(array_tar)
    tar = array_tar(i)
    do j = 1,size(array_in)
        if (array_in(j) == tar)then
            temp_int = temp_int + 1
            array_index(temp_int) = j
        end if
    end do
end do
end function

function sum_byindex(array_in,index_in) result(output)
use qsort_c_module
implicit none
integer,dimension(:)                :: array_in,index_in
integer,dimension(size(index_in))   :: output,temp,idx_d
integer                             :: i,LOC
call SORTRX(size(index_in),index_in,idx_d)
temp = index_in(idx_d)
output = 0
do i = 1,size(array_in)
    LOC = binarysearch(array_in(i),temp)
    if (LOC /= 0)then
        output(idx_d(LOC)) = output(idx_d(LOC)) + 1
    end if
end do
return
end function


function get_index_mul_char(array_in,array_tar) result(array_index)
integer size_in,index,i,j,size_tar
character(len=*), dimension(:)                    :: array_in
character(len=*), dimension(:)                    :: array_tar
integer, allocatable, dimension(:)                :: array_index

index = 0
size_in = size(array_in)
size_tar = size(array_tar)
do i = 1,size_in
    do j = 1, size_tar
        if (trim(array_tar(j)) == trim(array_in(i)))then
            index = index + 1
        endif
    end do
end do
allocate(array_index(index))

index = 0
do i = 1,size_in
    do j = 1,size_tar
        if (trim(array_tar(j)) == trim(array_in(i)))then
            index = index + 1
            array_index(index) = i
        endif
    end do
end do

end function

! function get_index_mul_char_and(array_in,array_tar) result(array_index)
! end function

function get_index_single_int_inv(array_in,tar) result(array_index)
implicit none
integer  size_in,index,i
real,dimension(:)                              :: array_in
real                                           :: tar
integer,allocatable                            :: array_index(:)
index=0
size_in = size(array_in)
do i=1,size_in
    if (tar /= array_in(i)) then
        index = index + 1
    endif
enddo
!allocate array and reset index number
allocate(array_index(index))
index=0
do i=1,size_in
    if (tar /= array_in(i)) then
        index = index + 1
        array_index(index) = i
    endif
enddo
end function

integer function get_file_N(iFileUnit)
implicit none
integer, Intent(IN) :: iFileUnit
character*1         :: cDummy
write(*,*)"operating on unit:",iFileUnit
get_file_N = 0
rewind(iFileUnit)

DO
read(iFileUnit,*,End= 999,Err = 999) cDummy
get_file_N = get_file_N + 1
end do
999 rewind(iFileUnit)
return
end function get_file_N

real function invert_char(Input_char)
implicit none
character(len=30)        Input_char
character*28        temp_char_in
real                temp_real
integer             temp_int,i
character           temp
i = 1
temp = Input_char(1)
write(temp,"(i2)")temp_int
temp_real = real(temp_int)
invert_char = temp_real
i = i + 1
do
    if (Input_char(i) == '.') then
        i = i + 1
        continue
    else if (Input_char(i) == ' ') then
        exit
    else
        temp = Input_char(i)
        write(temp,"(i2)")temp_int
        temp_real = real(temp_int)
        invert_char = invert_char + temp_real/(10**(i-2))
        i = i + 1
    end if
enddo
return
end function invert_char

function smoothing_change(x,steps,types)result(output)
implicit none
real,dimension(:)           :: x
real,dimension(size(x))     :: output
integer                     :: steps
character(len = *)          :: types

integer,dimension(size(x))  :: changex_index_temp
real,dimension(size(x)-1)   :: changex
integer,allocatable         :: changex_index(:)
real                        :: logval(steps),linval(steps),linval_temp(steps+2),start,final,res_temp(steps)
real                        :: inc_log, inc_lin
integer                     :: size_xnchange,i,size_x,count
size_x = size(x)
inc_log = (5 - (-5))/real(steps - 1)
inc_lin = (1 - 0)/real(steps + 2 -1)
output  = x
do i = 1,steps-1
    logval(i) = -5 + i*inc_log
end do

do i = 1,steps+2
    linval_temp(i) = 1 - (i-1) * inc_lin
end do
linval = linval_temp(2:size(linval_temp)-1)
changex = x(2:size_x)-x(1:size_x-1)
! changex_index = get_index_not(changex,0.0)
count = 0;

do i = 1,size(changex)
    if(abs(changex(i))>0.0)then
        count = count+1
        changex_index_temp(count) = i
    endif
end do
changex_index = changex_index_temp(1:count)

size_xnchange = size(changex_index)
    if (trim(types) == "logistic" .and. size_xnchange > 0) then
        do i=1,size(changex_index)
            start  = x(changex_index(i))
            final  = x(changex_index(i)+steps+1)
            
            logval = 1 - 1/( 1 + exp(-logval))
            logval = start * logval + final * (1 - logval)
            output(changex_index(i)+1:changex_index(i) + steps) = logval
        enddo
    endif
    
    if (trim(types) == "linear" .and. size_xnchange > 0) then
        do i=1,size(changex_index)
            start  = x(changex_index(i))
            final  = x(changex_index(i)+steps+1)
            
            res_temp = start * linval + final * (1 - linval)
            output(changex_index(i)+1:changex_index(i) + steps) = res_temp
        end do
    end if

end function
            
function seq_date(start_date,end_date)result(output)
implicit none
character *10             :: start_date,end_date
character *10,allocatable :: output(:)
integer                   :: time_start,time_end,size_out,i

time_start = Date2Unixtime(start_date)
time_end   = Date2Unixtime(end_date)

size_out = (time_end-time_start)/86400 + 1

allocate(output(size_out))

do i =1, size_out
    output(i) = Unixtime2Date(time_start + (i-1)*86400)
end do

end function

function seq_integer(start_n,end_n,advan)result(output)
implicit none
integer start_n,end_n,n,i,advan
integer,allocatable :: output(:)

n = (end_n-start_n)/advan

allocate(output(n+1))
output(1)  =  start_n
do i = 1, n
    output(i+1) = start_n + i * advan
end do

end function


! function seq_char(start_n,end_n,advan)result(output)
! implicit none
! integer start_n,end_n,n,i
! character,allocatable :: output(:)
! 
! n = (end_n-start_n)/advan
! 
! allocate(output(n+1))
! output(1)  =  start_n
! do i = 1, n
!     write(output(i+1),("I2"))start_n + i * advan
! end do
! 
! end function

function get_start_date(date_array)result(res)
    implicit none
    integer i,flag,minimal_time,time_loop
    character(len = *),dimension(:)        :: date_array
    character*10                           :: res
    flag = 1
    minimal_time  = Date2Unixtime(date_array(1))
    do i =2, size(date_array)
        time_loop = Date2Unixtime(date_array(i))
        if (time_loop < minimal_time) then
            flag = i
            minimal_time = time_loop
        endif
    end do
    res = date_array(flag)
end function

function rep(tar,times)result(output)
implicit none
integer,dimension(:)    :: tar
integer,dimension(:)    :: times
integer,allocatable     :: output(:)

integer                 :: i,size_out,cur_pos
size_out = sum(times)
allocate(output(size_out))
cur_pos = 0
do i = 1, size(times)
    output(cur_pos+1:cur_pos+times(i)) = tar(i)
    cur_pos = cur_pos + times(i)
end do

end function

function sample_i4(input_array,sample_N)result(sample_seq)
implicit none
integer,dimension(:)        :: input_array
integer,dimension(size(input_array)) :: temp_input
integer                     :: sample_N
integer,allocatable         :: sample_seq(:)
! integer,allocatable         :: output_array(:)

integer                     :: n,temp,pos,i
real                        :: ran

temp_input = input_array

do i = 1,size(temp_input)
     temp = temp_input(i)
     call random_number(ran)
     pos = int(ran*size(temp_input))
     do while (pos == 0)
        call random_number(ran)
        pos = int(ran*size(temp_input))
     end do
     temp_input(i) = temp_input(pos)
     temp_input(pos) = temp
end do
sample_seq = temp_input(1:sample_N)
! output_array = input_array(sample_N+1:size(input_array))
end function sample_i4

function sample_i1(input_array,sample_N)result(sample_seq)
implicit none
integer(kind=1),dimension(:)        :: input_array
integer,dimension(size(input_array)) :: temp_input
integer                     :: sample_N
integer,allocatable         :: sample_seq(:)
! integer,allocatable         :: output_array(:)

integer                     :: n,temp,pos,i
real                        :: ran

temp_input = input_array

do i = 1,size(temp_input)
     temp = temp_input(i)
     call random_number(ran)
     pos = int(ran*size(temp_input))
     do while (pos == 0)
        call random_number(ran)
        pos = int(ran*size(temp_input))
     end do
     temp_input(i) = temp_input(pos)
     temp_input(pos) = temp
end do
sample_seq = temp_input(1:sample_N)
! output_array = input_array(sample_N+1:size(input_array))
end function sample_i1

function sum_bygroup_mod(t1,dist_id,counties,mod_type)result(output)
implicit none
integer,dimension(:)                :: t1,dist_id,counties
integer,allocatable                 :: temp_index(:)
integer,dimension(size(counties))   :: output
character(len=*)                    :: mod_type
! allocate(output(size(counties)))
if (trim(mod_type)=="ill") then
    temp_index = get_index(t1,(/1,2,3/))
    output     = sum_byindex(dist_id(temp_index),counties)  
end if

if (trim(mod_type)=="healthy") then
    temp_index = get_index(t1,0)
    output     = sum_byindex(dist_id(temp_index),counties)   
    return
end if

if (trim(mod_type)=="inf_ncon") then
    temp_index = get_index(t1,1)
    output     = sum_byindex(dist_id(temp_index),counties)     
    return
end if

if (trim(mod_type)=="inf_con") then
    temp_index = get_index(t1,2)
    output     = sum_byindex(dist_id(temp_index),counties)    
    return
end if

if (trim(mod_type)=="ill_con") then
    temp_index = get_index(t1,3)
    output     = sum_byindex(dist_id(temp_index),counties)    
    return
end if

if (trim(mod_type)=="dead") then
    temp_index = get_index(t1,6)
    output     = sum_byindex(dist_id(temp_index),counties)     
    return
end if

end function
    

function get_index_not_single_integer(array_in,tar) result(output)   
implicit none
integer,dimension(:)        :: array_in
integer,dimension(size(array_in)) :: output_temp
integer                     :: tar,i,count
integer,allocatable         :: temp_index(:),output(:)
count = 0
do i = 1,size(array_in)
    if (array_in(i) /= tar)then
        count = count + 1
        output_temp(count) = i
    end if
end do
output = output_temp(1:count)
end function

function get_index_not_mul_integer(array_in,tar) result(output)
implicit none
integer,dimension(:)              :: array_in,tar
integer,allocatable               :: temp_index(:),output(:)
integer,dimension(size(array_in)) :: seq_index
integer                           :: i

temp_index = get_index(array_in,tar)
seq_index  = (/(i,i = 1,size(seq_index))/)
output     = array_minus(seq_index,temp_index)

end function

function get_index_not_single_real(array_in,tar) result(output)
real,dimension(:)        :: array_in
real,dimension(size(array_in)) :: array_temp
real                     :: tar
integer,allocatable      :: temp_index(:)
integer                  :: i,count_i
integer,allocatable      :: output(:)
count_i = 0
do i = 1,size(array_in)
    if ((array_in(i) - tar) > 0.00000001)then
        count_i = count_i + 1
    end if
end do
allocate(output(count_i))
count_i = 0
do i = 1,size(array_in)
    if ((array_in(i) - tar) > 0.00000001)then
        count_i = count_i + 1
        output(count_i) = i
    end if
end do
return


end function

function array_minus(array_in,minus)result(output)
implicit none
integer,dimension(:)        :: array_in
integer,dimension(size(array_in))   :: array_temp
integer,allocatable         :: minus(:),output(:)
integer                     :: i

! allocate(output(size(array_in)-size(minus)))
array_temp = array_in

array_temp(minus) = 0
output = pack(array_temp,array_temp.ne.0)



end function

function sum_bygroup_array(dist_id)result(output)
implicit none
integer,dimension(:)            :: dist_id
integer,allocatable             :: output(:),temp(:),counties(:)
integer                         :: i
counties = get_unique(dist_id)
allocate(output(size(counties)))
do i = 1,size(counties)
    temp = get_index(dist_id,counties(i))
    output(i) = size(temp)
end do

end function

subroutine sum_bygroup_distID(cases,output,dist_id)
implicit none
integer,dimension(:)            :: cases
integer,allocatable             :: dist_id(:),output(:),temp(:)
integer                         :: i

dist_id = get_unique(cases)
allocate(output(size(dist_id)))
output = 0
do i = 1,size(dist_id)
    temp = get_index(cases,dist_id(i))
    output(i) = size(temp)
end do

end subroutine

subroutine shift_cases(sick_new,mod_inf,dist_id,counties)
implicit none
integer,dimension(:)        :: sick_new,mod_inf,counties,dist_id
integer,allocatable         :: to_modifiy_index(:),seq_index(:),temp_index(:)
integer                     :: i,county


do i = 1,size(mod_inf)
    county = counties(i)
    temp_index = get_index(dist_id,county)
    if (mod_inf(i) < 0) then
        to_modifiy_index = get_index(sick_new(temp_index),1)
        if (size(to_modifiy_index) > abs(mod_inf(i))) then
            seq_index = sample(to_modifiy_index,abs(mod_inf(i)))
            sick_new(seq_index) = 0
        else
            sick_new(to_modifiy_index) = 0
        end if
    end if
    
    if (mod_inf(i) > 0) then
        to_modifiy_index = get_index(sick_new(temp_index),0)
        if (size(to_modifiy_index) > 0)then
            if (size(to_modifiy_index) > mod_inf(i)) then
                seq_index = sample(to_modifiy_index,abs(mod_inf(i)))
                sick_new(seq_index) = 1
            else 
                sick_new(to_modifiy_index) = 1
            endif
        endif
    endif
enddo
end subroutine

function find_and(array1,array2)result(output)
use qsort_c_module
implicit none
integer,dimension(:)    :: array1(:),array2(:)
integer,dimension(size(array1)) :: idx1,temp1
integer,dimension(size(array2)) :: idx2,temp2
integer,allocatable     :: output(:)
integer                 :: i,j,total,size1,size2

size1 = size(array1)
size2 = size(array2)

call SORTRX(size1,array1,idx1)
call SORTRX(size2,array2,idx2)

temp1 = array1(idx1)
temp2 = array2(idx2)
i = 1
j = 1
total = 0
do while((i <= size(temp1)) .and. (j <= size(temp2)))

    if (temp1(i) < temp2(j))then
        i = i + 1 
    else if (temp1(i) > temp2(j))then
        j = j + 1
    else
        total = total + 1
        i = i + 1
    endif
end do
  
allocate(output(total))
i = 1
j = 1
total = 0
do while((i <= size(temp1)) .and. (j <= size(temp2)))
    
    if (temp1(i) < temp2(j))then
        i = i + 1 
    elseif (temp1(i) > temp2(j))then
        j = j + 1
    else
        total = total + 1
        output(total) = temp1(i)
        i = i + 1
    endif
end do

end function

function condition_and_integer(array1,cond1,mod1,array2,cond2,mod2)result(output)
implicit none
integer,dimension(:)        :: array1,array2
integer                     :: cond1,cond2
character*1                 :: mod1,mod2
integer,allocatable         :: output(:),temp1(:),temp2(:)

if (trim(mod1) == "e") then
    temp1 = get_index(array1,cond1)
else if(trim(mod1) == "l") then
    temp1 = find_lessthan(array1,cond1)
else if(trim(mod1) == "g") then
    temp1 = find_greatthan(array1,cond1)
else 
    temp1 = get_index_not_single_integer(array1,cond1)
end if
! print *,temp1
if (trim(mod2) == "e") then   
    temp2 = get_index(array2,cond2)
else if(trim(mod2) == "l") then
    temp2 = find_lessthan(array2,cond2)
else if(trim(mod2) == "g") then   
    temp2 = find_greatthan(array2,cond2)
else
    temp2 =  get_index_not_single_integer(array2,cond2)
end if
! print *,temp2

output = find_and(temp1,temp2)



end function





!version 2, use logical expersion

! function condition_and_integer(array1,cond1,mod1,array2,cond2,mod2)result(output)
! implicit none
! integer,dimension(:)                    :: array1,array2
! integer                                 :: cond1,cond2
! character(len=*)                        :: mod1,mod2
! logical,dimension(size(array1))         :: temp1,temp
! logical,dimension(size(array2))         :: temp2
! integer,allocatable                     :: output(:)
! integer                                 :: temp_int,i
! 
! if (mod1 == "e") then
!     temp1 = (array1 == cond1)
! else if(mod1 == "le") then
!     temp1 = (array1 <= cond1)
! else if (mod1 == "l")then
!     temp1 = (array1 < cond1)
! else if(mod1 == "ge") then
!     temp1 = (array1 >= cond1)
! else if(mod1 == "g")then
!     temp1 = (array1 > cond1)
! else 
!     temp1 = (array1 /= cond1)
! end if
! ! print *,temp1
! if (mod2 == "e") then
!     temp2 = (array2 == cond2)
! else if(mod2 == "le") then
!     temp2 = (array2 <= cond2)
! else if (mod2 == "l")then
!     temp2 = (array2 < cond2)
! else if(mod2 == "ge") then
!     temp2 = (array2 >= cond2)
! else if(mod2 == "g")then
!     temp2 = (array2 > cond2)
! else 
!     temp2 = (array2 /= cond2)
! end if
! 
! temp = temp1.and.temp2
! temp_int = 0
! do i = 1,size(temp)
!     if (temp(i))then
!         temp_int = temp_int + 1
!     end if
! end do
! allocate(output(temp_int))
! temp_int = 0
! do i = 1,size(temp)
!     if (temp(i))then
!         temp_int = temp_int + 1
!         output(temp_int) = i
!     end if
! end do
! 
! return
! 
! end function

function condition_and_real(array1,cond1,mod1,array2,cond2,mod2)result(output)
implicit none
integer,dimension(:)        :: array1,array2
integer                     :: cond1
real                        :: cond2
character(len=*)            :: mod1,mod2
integer,allocatable         :: output(:),temp1(:),temp2(:)

if (mod1 == "e") then
    temp1 = get_index(array1,cond1)
else if(mod1 == "l") then
    temp1 = find_lessthan(array1,cond1)
else if(mod1 == "g") then
    temp1 = find_greatthan(array1,cond1)
else 
    temp1 = get_index_not_single_integer(array1,cond1)
end if

if (mod2 == "e") then   
    temp2 = get_index(array2,nint(cond2))
else if(mod2 == "l") then
    temp2 = find_lessthan(array2,nint(cond2))
else if(mod2 == "g") then   
    temp2 = find_greatthan(array2,nint(cond2))
else
    temp2 =  get_index_not_single_integer(array2,nint(cond2))
end if

output = find_and(temp1,temp2)

end function


function find_lessthan(array_in,cond)result(output)
implicit none
integer,dimension(:)        :: array_in
integer                     :: cond,i,total
integer,allocatable         :: output(:)

total = 0
do i = 1,size(array_in)
    if (array_in(i) < cond) then
        total = total + 1
    end if
end do

allocate(output(total))
total = 0
do i = 1,size(array_in)
    if (array_in(i) < cond) then
        total  = total + 1
        output(total) = i
    end if
end do
return
end function

function find_greatthan(array_in,cond)result(output)
implicit none
integer,dimension(:)        :: array_in
integer                     :: cond,i,total
integer,allocatable         :: output(:)

total = 0
do i = 1,size(array_in)
    if (array_in(i) >= cond) then
        total = total + 1
    end if
end do

allocate(output(total))
total = 0
do i = 1,size(array_in)
    if (array_in(i) >= cond) then
        total  = total + 1
        output(total) = i
    end if
end do

end function

function find_max_date(date_in)result(maxdate)
character(len = *),dimension(:)       :: date_in
integer,dimension(size(date_in)):: time_in
integer                         :: maxdate
integer                         :: i
do i = 1,size(date_in)
    time_in(i) = Date2Unixtime(date_in(i))
end do

maxdate = maxval(time_in)
return
end function

function Date2Unixtime(Date) result(Unixtime)

implicit none
Integer                     :: Unixtime
Character(len=*), Intent(in):: Date
Character*4 :: Cyear
Character*2 :: Cmonth,Cday
Integer     :: year,month,day
Integer     :: iy,im,a
real        :: JD

Cyear  = Date(1:4)
Cmonth = Date(6:7)
Cday   = Date(9:10)

!print*,year,month,day
! read(Cyear,*)year
read(Cyear,'(I4)')year
! read(Cmonth,*)month
read(Cmonth,'(I2)')month
! read(Cday,*)day
read(Cday,'(I2)')day

if (month>2) then
        iy = year
        im = month
else
        iy = year - 1
        im = month + 12
end if

a = INT(iy/100)
a = 2 - a +INT(a/4)
JD = INT(365.25*(iy + 4716)) + INT(30.60001*(im + 1)) + day + a -1524.5
Unixtime = INT((JD-2440587.5)*86400)
return
end function

function Unixtime2Date(Unixtime)result(Date)
implicit none
Integer, Intent(in) :: Unixtime
Character*10        :: Date
Character*4         :: Cyear
Character*2         :: Cmonth,Cday
integer             :: year,month,day,udays,mday
real                :: rday



udays = int(Unixtime/86400)
mday   = udays + 40587
year  = 1858 + int( (mday + 321.51) / 365.25)
rday  = aint( mod(mday + 262.25, 365.25) ) + 0.5
month = 1 + int(mod(rday / 30.6 + 2.0, 12.0) )
day   = 1 + int(mod(rday,30.6))


write(Cyear,'(I4)')year
if (month > 9)then
    write(Cmonth,'(I2)')month
else
    write(Cmonth,'(I2)')month
    Cmonth = "0"//trim(adjustl(Cmonth))
end if

if (day > 9)then
    write(Cday,'(I2)')day
else
    write(Cday,'(I2)')day
    Cday   = "0"//trim(adjustl(Cday))
end if
Date = trim(Cyear)//"-"//trim(Cmonth)//"-"//trim(Cday)
return
end function

function add_date(Date,days)result(res)
implicit none
character*10   Date
character*10   res
integer        days,time

time = Date2Unixtime(Date)
time = time + 86400*days
res = Unixtime2Date(time)
return
end function


end module support_fun
