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
!#               Huan Zhou
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
  Use genetic_algorithm   

  !! Legacy Modules --------------------------------------------------
  Use cosmic_io
  Use kernel

  use OMP_LIB
  use mpi

  Implicit None

  Type(iols)                    :: iol

  Integer(kind=ik)                              :: iter, iter_s, iter_e, n_iter
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

  Integer(kind=ik)                              :: lhc_samples
  Integer(kind=ik)                              :: n_lhc_samples
  Integer(kind=ik)                              :: lhc_samples_s
  Integer(kind=ik)                              :: lhc_samples_e

  Integer(kind=ik)                              :: cur_week
  
  Real(kind=rk),    Allocatable, Dimension(:)     :: R0
  Real(kind=rk),    Allocatable, Dimension(:,:)   :: R0_effects, R0_effects_lhc
  Character(len=:), Allocatable                   :: R0_effects_fn
  Character(len=:), Allocatable                   :: lhc_dir

  Character(len=:), Allocatable                   :: region
  Character(Len=:), allocatable, Dimension(:)     :: region_ids
  Character(Len=:), allocatable, Dimension(:)     :: region_ids_R0e
  integer(kind=ik), allocatable, dimension(:)     :: region_index

  character(len=:), Allocatable                   :: output_dir
  character(len=:), Allocatable                   :: export_name
  Character(len=pt_mcl)                           :: filename   
  Character(len=pt_mcl)                           :: logfile

  !-- MPI Variables -----------------------------------------------------------
  Integer(kind=mpi_ik)                            :: ierr
  Integer(kind=mpi_ik)                            :: rank_mpi, size_mpi

  Integer(kind=pt_ik), Dimension(:), Allocatable :: serial_pt
  Integer(kind=pt_ik)                            :: serial_pt_size

  Integer(kind=pt_ik)                            :: ii, jj
  Integer(kind=ik)                               :: num_counties
  Integer(kind=ik)                               :: num_regions
  !=============================================================================

  Character(len=:), allocatable                  :: exec_type
  Integer(kind=ik),dimension(:), Allocatable     :: full_region_index ! Required by "optimization"
  Integer,dimension(:,:,:), Allocatable          :: ill_ICU_cases_final ! Store the evaluated ICU cases
  Type(sims)                                     :: sim
  Type(opt_sim_switchs)                          :: sim_switch
  Call mpi_init(ierr)

  Call MPI_COMM_RANK(MPI_COMM_WORLD, rank_mpi, ierr)

  Call MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)

  call start_timer("CoSMic")
  
  !** Open new logfile ---------------------------------------------------------
  write(logfile,'(A,I6.6,A)')"CoSMic-proc",rank_mpi,".log"
  Open(newunit=un_lf, file=trim(logfile), action="write", status="replace")

  !** Set monitoring unit of param tree library to CoSMic logfile unit. -------
  pt_umon = un_lf

  !****************************************************************************
  !** Rank 0 -- Process master Init Process and broadcast init params *********
  !****************************************************************************
  If (rank_mpi==0) Then

     call print_cosmic_head()

     !==========================================================================
     ! Load Parameters =========================================================
     Call start_timer("Read parameters")

     Call read_param_file("exec_parameters.dat")
     Call read_param_file("static_parameters.dat")

     call write_param_tree(un_lf)

     Call end_timer("Read parameters")

     !==========================================================================
     ! Load Inout data =========================================================
     Call start_Timer("Load data")
     Call loaddata(iol)
     Call end_Timer("Load data")

     !==========================================================================
     ! Broadcast send parameter tree ===========================================
     Call start_Timer("Broadcast Param Tree")
     
     Call serialize_pt(serial_pt, serial_pt_size)
     
     Call mpi_bcast(serial_pt_size, 1_mpi_ik, MPI_INTEGER4, 0_mpi_ik,&
          MPI_COMM_WORLD, ierr)
     
     Call mpi_bcast(serial_pt, INT(serial_pt_size,mpi_ik), MPI_INTEGER4, 0_mpi_ik,&
          MPI_COMM_WORLD, ierr)
     
     Call end_Timer("Broadcast Param Tree")
     
  Else

     Call start_Timer("Broadcast Param Tree")

     !==========================================================================
     ! Broadcast recieve parameter tree ========================================
     Call mpi_bcast(serial_pt_size, 1_mpi_ik, MPI_INTEGER4, 0_mpi_ik,&
          MPI_COMM_WORLD, ierr)

     !** Serial_root_size == -1 ==> Signal that Domain_Number < size_mpi-1 ****
     If ( serial_pt_size == -1 ) then
        Goto 1001
     End If

     Allocate(serial_pt(serial_pt_size))

     Call mpi_bcast(serial_pt, INT(serial_pt_size,mpi_ik), MPI_INTEGER4, 0_mpi_ik,&
          MPI_COMM_WORLD, ierr)

     Call end_Timer("Broadcast Param Tree")

     Call start_Timer("Deserialize param tree")

     Call deserialize_pt(serial_pt)

     Call end_Timer("Deserialize param tree")

     call write_param_tree(un_lf)
   
  End If

  Call pt_get("#exec.procedure",exec_type)
  Call start_Timer("Broadcast input data")

  Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%counties)
  Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%seed)
  Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%trans_pr)
  Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%death)
  Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%pop)
  Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%connect_work)

  !** Specifically for "optimization" execution procedure 
  If (trim(exec_type) .EQ. "Optimization") Then
     Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%states)
     Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%obsicu_state)
     Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%obsicu_nuts2)
  End If

  call pt_get("#R0_effects",R0_effects_fn)

  If (trim(R0_effects_fn) .NE. "LHC") then
     Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%R0_effect)
  End If
   
  !Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%connect_total)
  !Call mpi_bcast_table(MPI_COMM_WORLD, 0_mpi_ik, rank_mpi, iol%connect_states)
  
  Call end_Timer("Broadcast input data")
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
  
  call pt_get("#region"  ,region  )

  call pt_get("#output.dir" , output_dir)
  call pt_get("#export_name", export_name)

  sim_switch%sim_reuse  = .false.
  sim_switch%sim_reload = .false.
  sim_switch%sim_write  = .false.
  
  If (trim(R0_effects_fn) .NE. "LHC") then
     R0_effects = table_to_real_array(iol%R0_effect)
  End If
  
  Call end_Timer("Load pt variables")

  If (trim(R0_effects_fn) .EQ. "LHC") then
     !=============================================================================
     ! Load input data from param tree ============================================
     Call start_Timer("Load LHC variables")

     call pt_get("#lhc.samples",lhc_samples)

     if ( mod(lhc_samples,size_mpi) .NE. 0 ) then
        n_lhc_samples = Int(lhc_samples / size_mpi ) + 1
     else
        n_lhc_samples = Int(lhc_samples / size_mpi )
     End if
     
     lhc_samples_s = rank_mpi * n_lhc_samples + 1
     lhc_samples_e = min((rank_mpi + 1 ) * n_lhc_samples , lhc_samples)
     
     write(*,*)rank_mpi,lhc_samples_s,lhc_samples_e,n_lhc_samples,lhc_samples

     iter_s = 1
     iter_e = iter

     call pt_get("#current.week",cur_week)
     call pt_get("#lhc.dir",lhc_dir)
     
     write(filename,'(3(A,I0),A)')"lhc_var_",cur_week-2,&
          "_constant_",lhc_samples_s,"_",lhc_samples_e,".csv"
     
     call read_TableData( &
          trim(lhc_dir)//"/"//trim(filename),sep=" ",head=.TRUE., &
          rownames=.TRUE., data=iol%R0_effect &
          )
 
     R0_effects_lhc = table_to_real_array(iol%R0_effect)
    
     if (trim(region) == "states") then

        region_index   = int(get_int_table_column(iol%counties,"dist_id")/1000)
        num_regions    = 16

     Else
        
        region_ids     = get_char_column(iol%counties,'Nuts2')
        region_ids_R0e = iol%R0_effect%head(1:38*(cur_week+1):cur_week+1)(1:4)
        
        num_counties = size(region_ids)
        allocate(region_index(num_counties))
        
        do ii = 1, num_counties
           do jj = 1, size(region_ids_R0e)
              if (trim(region_ids(ii)) == trim(region_ids_R0e(jj)))then
                 region_index(ii) = jj
              endif
           end do
        end do
        
     End if
     
     Call end_Timer("Load LHC variables")

     Call start_Timer("Exec. Simulation")
     Do ii = 1, n_lhc_samples

        R0_effects = reshape(R0_effects_lhc(ii,1:num_regions*(cur_week+1)),(/cur_week+1,num_regions/))

       ill_ICU_cases_final = COVID19_Spatial_Microsimulation_for_Germany(iol, &
             iter_s, iter_e , &
             inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
             less_contagious, R0_force, immune_stop, &
             R0change, R0delay ,R0delay_days, R0delay_type, &
             control_age_sex, seed_date, seed_before, sam_size, R0, &
             R0_effects, region_index, sim, sim_switch, &
             output_dir, export_name, rank_mpi)
        
     End Do
     Call end_Timer("Exec. Simulation")

  Else

     !==========================================================================
     ! Execute single model instance ===========================================
     Call start_Timer("Exec. Simulation")
     
     if (trim(region) == "states") then

        region_index   = int(get_int_table_column(iol%counties,"dist_id")/1000)
        
     Else

        region_ids     = get_char_column(iol%counties,'Nuts2')

        !write(*,*) "region_ids:",region_ids

        num_counties = size(region_ids)
        allocate(region_index(num_counties))
        
        do ii = 1, num_counties
           do jj = 1, size(iol%R0_effect%head)
              if (trim(region_ids(ii)) == trim(iol%R0_effect%head(jj)))then
                 region_index(ii) = jj
              endif
           end do
        end do
        
     end if

     write(un_lf,PTF_M_AI0)"size(region_index) =",size(region_index)

     if (trim(exec_type) .EQ. "Optimization") then ! Optimization

        full_region_index = region_index
        iter_s = 1
        iter_e = iter
        
        call ga(&
             iol, &
             iter_s, iter_e, &
             inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
             less_contagious, R0_force, immune_stop, &
             R0change, R0delay, R0delay_days, R0delay_type, &
             control_age_sex, seed_date, seed_before, sam_size, R0, &
             R0_effects, full_region_index, region_index, output_dir, export_name, rank_mpi, size_mpi)

     Else ! Basic-Param
        
        if ( mod(iter,size_mpi) .NE. 0 ) then
           n_iter = Int(iter / size_mpi ) + 1
        else
           n_iter = Int(iter / size_mpi )
        End if
        
        iter_s = rank_mpi * n_iter + 1
        iter_e = min((rank_mpi + 1 ) * n_iter , iter)
        
        write(un_lf,PTF_M_nAI0)"rank_mpi:",rank_mpi,"iter_s:",iter_s,"iter_e:",iter_e,&
             "n_iter:",n_iter,"iter:",iter

        ill_ICU_cases_final = COVID19_Spatial_Microsimulation_for_Germany(iol, &
             iter_s, iter_e , &
             inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
             less_contagious, R0_force, immune_stop, &
             R0change, R0delay ,R0delay_days, R0delay_type, &
             control_age_sex, seed_date, seed_before, sam_size, R0, &
             R0_effects, region_index, sim, sim_switch, output_dir, export_name, rank_mpi)
     endif
     
     Call end_Timer("Exec. Simulation")

  End If
  
1001 Continue
  
  call end_timer("CoSMic")

  call write_timelist(unit=un_lf)

  Call MPI_FINALIZE(ierr)

End Program CoSMic
