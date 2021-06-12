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
! Module for data preparation
!
!###############################################################################
module data_preprocessing

  use global_constants
  use global_types
  use support_fun

  implicit none
  
contains

  !=============================================================================
  subroutine data_prepro(iol,counties_integer, &
       sp)

    Type(static_parameters), intent(in) :: sp

    type(iols)              :: iol
    integer,allocatable     :: counties_integer(:)
    integer,allocatable     :: all_counties(:)
    integer,allocatable     :: state_code_index(:)
    integer,allocatable     :: counties_index(:)
    integer,allocatable     :: temp_index(:)
    real,allocatable        :: connect_total(:,:)
    real,allocatable        :: connect_work(:,:)

    integer                 :: i

    !---------------------------------------------------------------------------
    if (import_R0_matrix) then
       !read file here needs to be implemented
    end if

    !list of all counties codes
    all_counties = get_unique(iol%pop_distid)

    if (sp%restrict) then
       !if(all(counties == "none"))
       ! don't know how to implement it wisely
       ! just leave blank here, since it wouldn't
       ! affect the code
       state_code_index = get_index_mul_char(iol%states_name,sp%sim_regions)
       !             print *,state_code_index
       counties_index   = get_index_mul_integer(all_counties/1000,iol%states_code(state_code_index))
       !             allocate(counties_integer(size(counties_index)))
       counties_integer = all_counties(counties_index)

       temp_index = get_index(iol%pop_distid,counties_integer)
       iol%pop_distid = iol%pop_distid(temp_index)
       iol%pop_date   = iol%pop_date(temp_index)
       iol%pop_sex    = iol%pop_sex(temp_index)
       iol%pop_age    = iol%pop_age(temp_index)
       iol%pop_total  = iol%pop_total(temp_index)

    else

       write(*,*)"restrict = .FALSE. is is not yet implemented."
       write(*,*)"Program halted !!!"
       stop

    end if


    ! connectitvity data
    if (sp%restrict) then
       temp_index = get_index_mul_integer(iol%connect_total_distid,counties_integer)
       !     allocate(connect_total(size(temp_index),size(temp_index)))
       connect_total = iol%connect_total(temp_index,temp_index)

       !     deallocate(temp_index)
       temp_index = get_index_mul_integer(iol%connect_total_distid,counties_integer)
       !     allocate(connect_work(size(temp_index),size(temp_index)))
       connect_work = iol%connect_total(temp_index,temp_index)

       ! rescale: rows sum to 1
       if (size(connect_total,dim=1) == size(connect_work,dim=1)) then
          do i = 1, size(connect_total,dim=1)
             connect_total(:,i)  =  connect_total(:,i)/sum(connect_total(:,i))
             connect_work(:,i)   =  connect_work(:,i)/sum(connect_work(:,i))
          end do
       else 
          do i = 1, size(connect_total,dim=1)
             connect_total(:,i)  =  connect_total(:,i)/sum(connect_total(:,i))
          end do
          do i = 1, size(connect_work,dim=1)
             connect_work(:,i)   =  connect_work(:,i)/sum(connect_work(:,i))
          end do
       end if

       ! should i tranpose the matrix also in fortran?
       ! data processing for iol is done

    end if
  end subroutine data_prepro

end module data_preprocessing
