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
! Module containing the CoSMic model loop
!
!###############################################################################
module list_variable
  
  use global_constants
  use global_types

    implicit none

    
    type,public :: seeds
        integer,allocatable             :: dist_id(:)
        integer,allocatable             :: cases(:)
        character*10,allocatable        :: date(:)

    end type seeds
    
!!$    interface   set_value
!!$        module procedure  set_value_seed
!!$        module procedure  set_value_iol
!!$    end interface
contains
    function init_seed(num)
        implicit none
        type(seeds) :: init_seed   
        integer     :: num
        allocate(init_seed%dist_id(num))
        allocate(init_seed%cases(num))
        allocate(init_seed%date(num))
    end function
    
!!$    function set_value_seed(seed,index_list)
!!$        implicit none
!!$        type(seeds)  :: set_value_seed
!!$        type(seeds)  :: seed
!!$        integer      :: index_list(:)
!!$!         integer,dimension(:)      :: value_id
!!$!         integer,dimension(:)      :: value_cases
!!$!         character*10,dimension(:) :: value_date
!!$        
!!$        set_value_seed%dist_id = seed%dist_id(index_list)
!!$        set_value_seed%cases   = seed%cases(index_list)
!!$        set_value_seed%date    = seed%date(index_list)
!!$    end function
    
!!$    function set_value_iol(iol,index_list)
!!$        implicit none
!!$        type(iols)   :: iol
!!$        type(seeds)  :: set_value_iol
!!$        integer      :: index_list(:)
!!$        set_value_iol%dist_id  = iol%seed_distid(index_list)
!!$        set_value_iol%date     = iol%seed_date(index_list)
!!$        set_value_iol%cases    = iol%seed_cases(index_list)
!!$    end function
    
end module list_variable
