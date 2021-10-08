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

  use OMP_LIB
  
  Implicit None

  Type(iols)                    :: iol

  Integer,Allocatable           :: counties_integer(:)

  Type(static_parameters)       :: sp
  Type(exec_parameters)         :: ep

  Integer(kind=ik)                              :: iter
  Integer(kind=ik)                              :: inf_dur
  Integer(kind=ik)                              :: cont_dur
  Integer(kind=ik)                              :: ill_dur
  Integer(kind=ik)                              :: icu_dur
  Integer(kind=ik), Allocatable, Dimension(:)   :: icu_per_day
  Real(kind=rk)                                 :: less_contagious
  Real(kind=rk)                                 :: R0_force
  Logical                                       :: immune_stop
  Integer(kind=ik), Allocatable, Dimension(:,:) :: R0change
  Logical                                       :: R0delay
  Integer(kind=ik)                              :: R0delay_days
  Character(len=:), Allocatable                 :: R0delay_type
  character(len=:), Allocatable                 :: control_age_sex
  character(len=:), Allocatable                 :: seed_date
  Integer(kind=ik)                              :: seed_before    
  Integer(kind=ik)                              :: sam_size
  Real(kind=rk)                                 :: R0
  
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

  !=============================================================================
  ! Prepare input data =========================================================
  Call start_Timer("Prepare data")
  Call data_prepro(iol,counties_integer)
  Call end_Timer("Prepare data")

  !=============================================================================
  ! Load input data from param tree ============================================
  Call start_Timer("Load pt variables")

  call pt_get("#iter",iter)
  
  call pt_get("#inf_dur" ,inf_dur ) ! 3
  call pt_get("#cont_dur",cont_dur) ! 2
  call pt_get("#ill_dur" ,ill_dur ) ! 8
  call pt_get("#icu_dur" ,icu_dur ) ! 14
  
  call pt_get("#icu_per_day" ,icu_per_day )
  
  call pt_get("#less_contagious" ,less_contagious )
  
  call pt_get("#R0_force",R0_force) ! 0
  
  call pt_get("#immune_stop",immune_stop)! = .True.
  
  call pt_get("#R0change"     ,R0change    )
  call pt_get("#R0delay"      ,R0delay     )
  call pt_get("#R0delay_days" ,R0delay_days)
  call pt_get("#R0delay_type" ,R0delay_type)
  
  call pt_get("#control_age_sex",control_age_sex) ! "age"
  
  call pt_get("#seed_date",seed_date)
  call pt_get("#seed_before",seed_before)
  
  call pt_get("#sam_size",sam_size)
  call pt_get("#R0"      ,R0      )

  Call end_Timer("Load pt variables")
  
  !=============================================================================
  ! Execute model ==============================================================
  Call start_Timer("Exec. Simulation")

  Call COVID19_Spatial_Microsimulation_for_Germany(iol,counties_integer, &
       iter , &
       inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
       less_contagious, R0_force, immune_stop, &
       R0change, R0delay ,R0delay_days, R0delay_type, &
       control_age_sex, seed_date, seed_before, sam_size, R0, &
       iol%R0_effect%data)
  
  Call end_Timer("Exec. Simulation")

  call end_timer("CoSMic")

  call write_timelist(unit=un_lf)

End Program CoSMic
