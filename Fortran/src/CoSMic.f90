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

  Use global_constants
  Use global_types
  Use global_vars  
  
  !! Legacy Modules --------------------------------------------------
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

  !** Open new logfile ---------------------------------------------------------
  Open(newunit=un_lf, file="CoSMic.log", action="write", status="replace")

  !** Set monitoring unit of param tree library to CoSMic logfile unit. -------
  pt_umon = un_lf
  
  call start_timer("CoSMic")
  call print_cosmic_head()

  !=============================================================================
  ! Load Parameters ============================================================
  Call start_timer("Read parameters")

  Call read_param_file("exec_parameters.dat")
  Call read_param_file("static_parameters.dat")

  call write_param_tree(un_lf)

  Call end_timer("Read parameters")

  !=============================================================================
  ! Load Inout data ============================================================
  Call start_Timer("Load data")
  Call loaddata(iol, ep, sp)
  Call end_Timer("Load data")

!!!!----initalizing pspace ============================================
  pspace%Ps_scalar_list(2)%param                 = 3.5
  pspace%Ps_scalar_list(2)%var_type              = "direct"
  pspace%Ps_scalar_list(2)%name                  = "R0"
  !icu_dur
  pspace%Ps_scalar_list(3)%param                 = 14
  pspace%Ps_scalar_list(3)%var_type              = "direct"
  pspace%Ps_scalar_list(3)%name                  = "icu_dur"
  !mod_surv_ill
  pspace%Ps_scalar_list(4)%param                 = 1
  pspace%Ps_scalar_list(4)%var_type              = "direct"
  pspace%Ps_scalar_list(4)%name                  = "mod_surv_ill"
  !lcokdown_effect
  pspace%Ps_scalar_list(5)%param                 = 0.39
  pspace%Ps_scalar_list(5)%var_type              = "direct"
  pspace%Ps_scalar_list(5)%name                  = "lcokdown_effect"
  !w_int
  pspace%Ps_scalar_list(6)%param                 = 0.9
  pspace%Ps_scalar_list(6)%var_type              = "direct"
  pspace%Ps_scalar_list(6)%name                  = "w_int"
  !w_obs
  pspace%Ps_scalar_list(7)%param                 = 0.0
  pspace%Ps_scalar_list(7)%var_type              = "direct"
  pspace%Ps_scalar_list(7)%name                  = "w_obs"
  !w_obs_by_state
  pspace%Ps_scalar_list(8)%param                 = 0.0
  pspace%Ps_scalar_list(8)%var_type              = "direct"
  pspace%Ps_scalar_list(8)%name                  = "w_obs_by_state"
  !ROeffect_ps

  !=============================================================================
  ! Prepare input data =========================================================
  Call start_Timer("Prepare data")
  Call data_prepro(iol,counties_integer)
  Call end_Timer("Prepare data")

  !=============================================================================
  ! Execute model ==============================================================
  Call start_Timer("Exec. Simulation")
  write(*,*)"COVID19_Spatial_Microsimulation_for_Germany - Start"
  Call COVID19_Spatial_Microsimulation_for_Germany(iol,pspace,counties_integer)
  Call end_Timer("Exec. Simulation")

  call end_timer("CoSMic")

  call write_timelist(unit=un_lf)

End Program CoSMic
