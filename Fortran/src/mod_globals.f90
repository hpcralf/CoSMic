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

!===============================================================================
!> \file mod_globals.f90
!> Modules with global constants, types and variables.
!>
!> \author Ralf Schneider
!> \date 20.08.2021

!===============================================================================
!> Module with globally used constants
Module global_constants

  Implicit None

  Integer,parameter   :: ROOT    = 0 ! The root rank

  character(len=*), parameter :: fmt_file_missing = &
       '("Operation on file:",/,A,/, "failed.",/, "Program halted !!!")'
  
  character(len=*), parameter :: fmt_timing = &
       '(A,T40,"finished in",F10.3,"sec.")'

  character(len=*),Dimension( 9),parameter :: mandatory_params = &
       ["trans_pr     ", "pop_data     ", "inf_cases    ", "dead_cases   ",&
       "connect_total", "connect_work ", "states       ", "counties     ",&
       "R0_effects   "]
 
End Module global_constants

!===============================================================================
!> Module with globally used derived data-types
module global_types

  use precision
  use global_constants

  implicit none

  public
  save

  type TColumnData
 
     character(len=:), Dimension(:), Allocatable :: cd_c
     Integer(Kind=ik), Dimension(:), Allocatable :: cd_i
     Real(Kind=rk)   , Dimension(:), Allocatable :: cd_r
          
  End type TColumnData

  type TTableData
 
     character(len=:), Dimension(:), allocatable :: head        !x

     Integer(kind=4) , Dimension(2)              :: data_size   !x

     Character       , Dimension(:), allocatable :: col_types   !x
     Integer(kind=4), Dimension(:), allocatable  :: col_lengths !x
     
     character(len=:), Dimension(:), allocatable :: rownames    !x
     
     Type(TcolumnData),Dimension(:)  , allocatable :: data      !
          
  End type TTableData
  
  !> Type to hold all baseline input data --------------------------------------
  type iols

     Type(TTableData)                       :: R0_effect
     Type(TTableData)                       :: counties
     Type(TTableData)                       :: seed
     Type(TTableData)                       :: death
     Type(TTableData)                       :: trans_pr
     Type(TTableData)                       :: pop
     Type(TTableData)                       :: connect_total
     Type(TTableData)                       :: connect_work
     Type(TTableData)                       :: states

     Type(TTableData)                       :: obsicu_state
     Type(TTableData)                       :: obsicu_nuts2

  end type iols

  type opt_sim_switchs
     Logical                                :: sim_reuse  ! Reuse the previous sim data
     Logical                                :: sim_reload ! Reloading checkpoint 
     Logical                                :: sim_write  ! Writing restart data
     Integer                                :: read_week  ! Reloading checkpoint of the given week
     Integer                                :: write_week ! Writing restart data for the given week
  end type 

end module global_types

!===============================================================================
!> Module with globally used variables
module global_vars

  Use ISO_FORTRAN_ENV
  
  use precision
  use global_constants
  use global_types

  implicit none

  integer :: time_n
  integer :: un_lf  = OUTPUT_UNIT
  
End module global_vars
