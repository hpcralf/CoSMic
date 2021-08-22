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
!###############################################################################
Program CoSMic

  !! Reworked Modules ------------------------------------------------
  use param_tree
  use timer
  
  !! Legacy Modules --------------------------------------------------
  Use global_constants
  Use global_types
  Use cosmic_io
  Use data_preprocessing
  Use kernel
  
  Implicit None

  Type(iols)                    :: iol
  Type(pspaces)                 :: pspace
  Integer                       :: iter
  Integer,Allocatable           :: counties_integer(:)

  Type(static_parameters)       :: sp
  Type(exec_parameters)         :: ep

  call start_timer("CoSMic")
  call print_cosmic_head()
  
  !=============================================================================
  ! Load Parameters ============================================================
  Call start_timer("Read parameters")
  Call read_param_file("exec_parameters.dat")
  Call read_param_file("static_parameters.dat")

  call write_param_tree(6)

  Call load_exec_parameters(ep,"exec_parameters.dat")
  Call log_exec_parameters(ep)

  !=============================================================================
  ! Load Static Parameters =====================================================
  Call load_static_parameters   (sp,"static_parameters.dat")
  Call log_static_parameters(sp)
  Call end_timer("Read parameters")
stop
  !=============================================================================
  ! Load Inout data ============================================================
  Call start_Timer("Load data")
  Call loaddata(iol, pspace, ep, sp)
  Call end_Timer("Load data")

  !=============================================================================
  ! Prepare input data =========================================================
  Call start_Timer("Prepare data")
  Call data_prepro(iol,counties_integer,sp)
  Call end_Timer("Prepare data")
  
  !=============================================================================
  ! Execute model ==============================================================
  Call start_Timer("Exec. Simulation") 
  Call COVID19_Spatial_Microsimulation_for_Germany(sp,iol,pspace,counties_integer)
  Call end_Timer("Exec. Simulation")

  call end_timer("CoSMic")

  call write_timelist()
  
End Program CoSMic
