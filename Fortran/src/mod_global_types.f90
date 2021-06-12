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
! Module with globally used derived data-types
!
!###############################################################################
module global_types

  use precision
  use global_constants

  implicit none

  public
  save

  !! !!!!!!variables for iol which is a list in R code
  type iols
     ! variables for trans_pr
     character*15,allocatable,dimension(:) :: transpr_age_gr
     character*15,allocatable,dimension(:) :: transpr_sex
     real,allocatable,dimension(:)         :: transpr_surv_ill
     real,allocatable,dimension(:)         :: transpr_icu_risk
     real,allocatable,dimension(:)         :: transpr_surv_icu
     character*15,dimension(5)             :: titel
     ! variables for pop
     integer,allocatable,dimension(:)      :: pop_distid
     character(10),allocatable,dimension(:):: pop_date
     character(1),allocatable,dimension(:) :: pop_sex
     character(2),allocatable,dimension(:) :: pop_age
     integer,allocatable,dimension(:)      :: pop_total
     ! variable for seed
     integer,allocatable,dimension(:)      :: seed_distid
     character*10,allocatable,dimension(:) :: seed_date
     integer,allocatable,dimension(:)      :: seed_cases
     character*15,dimension(3)             :: seed_titel
     !variables for seed_dea
     integer,allocatable,dimension(:)      :: death_distid
     character*10,allocatable,dimension(:) :: death_date
     integer,allocatable,dimension(:)      :: death_cases
     !variable for connect_total
     integer,dimension(401)                :: connect_total_distid
     character*15,dimension(401)           :: connect_total_name
     real,dimension(401,401)               :: connect_total
     character*15                          :: connect_titel
     !variable for connect_work
     integer,dimension(401)                :: connect_work_distid
     character*15,dimension(401)           :: connect_work_name
     real,dimension(401,401)               :: connect_work
     !variable for states
     integer,allocatable,dimension(:)      :: states_code
     integer,allocatable,dimension(:)      :: states_inhabitant
     character*15,allocatable,dimension(:) :: states_shortcut
     character*30,allocatable,dimension(:) :: states_name
     character*15,dimension(4)             :: state_titel
     !variable for counties
     integer,allocatable,dimension(:)      :: counties_dist_id
     character*50,allocatable,dimension(:) :: counties_name
     real,allocatable,dimension(:)         :: counties_area
     integer,allocatable,dimension(:)      :: counties_inhabitants
  end type iols

  !variable declaration for pspace
  type ps_scalar
     real                                  :: param
     character*10                          :: var_type
     character*10                          :: name
  end type ps_scalar

  type ps_vector
     real,allocatable                      :: param(:,:)
     character*20,allocatable              :: param_char(:,:)                        
     character*10                          :: var_type
  end type ps_vector

  type pspaces
     type(ps_scalar),dimension(8)          :: Ps_scalar_list
     type(ps_vector)                       :: ROeffect_ps
  end type pspaces

  !-------------------------------------------------------------------------------
  !> Type to hold execution parameters
  type exec_parameters

     character(len=mcl) :: exec_procedure
     character(len=mcl) :: data_dir
     character(len=mcl) :: output_dir
     character(len=mcl) :: model_version
     character(len=mcl) :: export_name

  end type exec_parameters

  !-------------------------------------------------------------------------------
  !> Type to hold static parameters
  type static_parameters

     character(len=mcl), Dimension(:), Allocatable :: sim_regions

     character(len=mcl)                            :: country
     character(len=10)                             :: seed_date
     Logical                                       :: restrict
     integer                                       :: lhc_samples

     ! File names of input files -----------------------------
     character(len=mcl)                            :: trans_pr=""     
     character(len=mcl)                            :: pop_data=""     
     character(len=mcl)                            :: inf_cases=""    
     character(len=mcl)                            :: dead_cases=""   
     character(len=mcl)                            :: connect_total=""
     character(len=mcl)                            :: connect_work="" 
     character(len=mcl)                            :: states=""       
     character(len=mcl)                            :: counties=""     
     character(len=mcl)                            :: R0_matrix_inp=""
     
  end type static_parameters

  integer                       :: time_n
  logical                       :: import_R0_matrix

end module global_types
