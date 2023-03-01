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
! Module containing the CoSMic model loop and time loop subroutines
!
!###############################################################################
Module kernel

  use timer

  use param_tree

  use precision
  Use global_constants
  use global_types
  use global_vars
  
  Use list_variable
  Use support_fun
  use qsort_c_module
  use urandom
  use CoSMic_IO
  
  use OMP_LIB
  
  Implicit None

  Type icu_risk_lists ! in order to have a similar structure to R
     Character*2,Allocatable         :: age(:)
     Integer,Allocatable             :: agei(:)
     Character*1,Allocatable         :: sex(:)
     Real(kind=rk),Allocatable       :: risk(:)
     Integer,Allocatable             :: dur(:)
  End Type icu_risk_lists

  Type sims
     Integer,Allocatable             :: dist_id(:)
     Integer,Allocatable             :: dist_id_rn(:)
     Integer,Allocatable             :: sample_index(:)
     Character,Allocatable           :: sex(:)
     Integer(kind=1),Allocatable     :: age(:)
     Integer(kind=1),Allocatable     :: t1(:)
     Integer(kind=1),Allocatable     :: t2(:)
     Integer,Allocatable             :: d(:)
  End Type sims

  !> states an individual can reach
  Integer, Parameter ::   missing    = -1
  Integer, Parameter ::   healthy    =  0
  Integer, Parameter ::   inf_noncon =  1
  Integer, Parameter ::   inf_contag =  2
  Integer, Parameter ::   ill_contag =  3
  Integer, Parameter ::   ill_ICU    =  4
  Integer, Parameter ::   immune     =  5
  Integer, Parameter ::   dead       =  6

  Integer, Parameter ::   min_state  = -1
  Integer, Parameter ::   max_state  =  6

  ! Declare the interface for POSIX fsync function
  interface
     function fsync (fd) bind(c,name="fsync")
       use iso_c_binding, only: c_int
       integer(c_int), value :: fd
       integer(c_int) :: fsync
     end function fsync
  end interface
  
Contains
  
  function COVID19_Spatial_Microsimulation_for_Germany( &
       iol, &
       iter_s, iter_e, &
       inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
       less_contagious, R0_force, immune_stop, &
       R0change, R0delay ,R0delay_days, R0delay_type, &
       control_age_sex, seed_date, seed_before, sam_size, R0, &
       R0_effects, region_index, sim_backup, sim_switch, output_dir, export_name, proc ) Result(ill_ICU_cases_final)

    !===========================================================================
    ! Declaration
    !===========================================================================

    Type(iols), Target                           , intent(in)   :: iol
    Integer(kind=ik)                             , intent(in)   :: iter_s, iter_e
    Integer(kind=ik)                             , intent(in)   :: inf_dur
    Integer(kind=ik)                             , intent(in)   :: cont_dur
    Integer(kind=ik)                             , intent(in)   :: ill_dur
    Integer(kind=ik)                             , intent(in)   :: icu_dur
    Integer(kind=ik), Allocatable, Dimension(:)  , intent(in)   :: icu_per_day
    Real(kind=rk)                                , intent(in)   :: less_contagious
    Real(kind=rk)                                , intent(in)   :: R0_force
    Logical                                      , intent(in)   :: immune_stop
    Integer(kind=ik), Allocatable, Dimension(:,:), intent(in)   :: R0change
    Logical                                      , intent(in)   :: R0delay
    Integer(kind=ik)                             , intent(in)   :: R0delay_days
    Character(len=*)                             , intent(in)   :: R0delay_type
    character(len=*)                             , intent(in)   :: control_age_sex
    character(len=*)                             , intent(in)   :: seed_date
    Integer(kind=ik)                             , intent(in)   :: seed_before    
    Integer(kind=ik)                             , intent(in)   :: sam_size
    Real(kind=rk),    Allocatable, Dimension(:)  , intent(in)   :: R0
    Real(kind=rk),    Allocatable, Dimension(:,:), intent(in)   :: R0_effects
    integer(kind=ik), allocatable, dimension(:)  , intent(in)   :: region_index
    Type(sims)                                   , intent(inout):: sim_backup
    Type(opt_sim_switchs)                        , intent(in)   :: sim_switch
    character(len=:), Allocatable                , intent(in)   :: output_dir
    character(len=:), Allocatable                , intent(in)   :: export_name
    Integer(kind=mpi_ik)                         , intent(in)   :: proc
    
    !--------------------------------------------------------------------------

    Integer(kind=ik)   :: i, j, index, temp_int,icounty,county,it_ss,status
    character(len=:), Allocatable   :: seed_date_mod
    Integer(kind=ik), Allocatable , Dimension(:) :: counties_index
    Integer(kind=ik)                             :: n_iter
    
    Character(Len=10)             :: seed_before_char,seed_temp
    Character(Len=10),Allocatable :: seed_seq(:),seed_inf_cont_seq(:),seed_inf_ncont_seq(:)
    Character(Len=10),Allocatable :: seed_d_seq(:)

    Integer(kind=ik)              :: days

    Integer,Allocatable :: temp(:)

    Type(sims)                      :: sim

    Integer,Allocatable             :: rownumbers(:)

    Type(seeds)                     :: seed_ill,seed_inf_cont,seed_inf_ncont,seed_death
    Integer,Allocatable             :: seed_ill_dur(:),seed_inf_cont_dur(:),seed_inf_ncont_dur(:)
    Integer,Allocatable             :: il_d(:),inf_c_d(:),inf_nc_d(:)
    Integer,Allocatable             :: rownumbers_ill(:),rownumbers_cont(:),rownumbers_ncont(:)
    Integer,Allocatable             :: rownumbers_left(:),rownumbers_dea(:)
    Integer,Allocatable             :: gettime(:),ur(:)
    Real(kind=rk),Allocatable                :: getchange(:)
    Integer                         :: inf_ill,inf_cont,inf_ncont,inf_dth

    Real(kind=rk)                            :: R0_daily
    Real(kind=rk),Allocatable                :: R0matrix(:,:),connect(:,:)
    Integer,Allocatable             :: healthy_cases_final(:,:,:)
    Integer,Allocatable             :: ill_ICU_cases_final(:,:,:)
    Integer,Allocatable             :: immune_cases_final(:,:,:)
    Integer,Allocatable             :: inf_noncon_cases_final(:,:,:)
    Integer,Allocatable             :: inf_contag_cases_final(:,:,:)
    Integer,Allocatable             :: dead_cases_final(:,:,:)
    Integer,Allocatable             :: ill_contag_cases_final(:,:,:)

    Integer                         :: timestep
    
    Integer,Allocatable             :: tmp_d_new(:),tmp_count(:)

    Character*10                    :: temp_date

    Integer                         :: max_date,n_change
    character(len=10)               :: l_seed_date
    Real(kind=rk)                            :: iter_pass_handle(6)

    Character(len=8)                :: proc_char
    
    Type(tTimer),pointer            :: timer

    Integer(kind=ik), Allocatable, Dimension(:)     :: dist_id_cref 

    Integer(kind=ik)                                :: pop_size
    Integer(kind=ik)                                :: num_counties
    Integer(kind=ik)                                :: ii,jj,kk

    Real(Kind=rk)   , Allocatable, Dimension(:,:)   :: surv_ill_pas
    Real(Kind=rk)   , Allocatable, Dimension(:,:,:) :: ICU_risk_pasd
    Real(Kind=rk)   , Allocatable, Dimension(:,:)   :: surv_icu_pas

    Character(Len=:), allocatable, Dimension(:)     :: inf_seed_date
    Integer(kind=ik), pointer    , Dimension(:)     :: inf_seed_distid
    Integer(kind=ik), pointer    , Dimension(:)     :: inf_seed_cases

    Character(Len=:), allocatable, Dimension(:)     :: dea_seed_date
    Integer(kind=ik), pointer    , Dimension(:)     :: dea_seed_distid
    Integer(kind=ik), pointer    , Dimension(:)     :: dea_seed_cases

    Integer(kind=ik), pointer    , Dimension(:)     :: pop_total    
    Integer(kind=ik), pointer    , Dimension(:)     :: pop_distid
    Character(Len=:), allocatable, Dimension(:)     :: pop_sex
    Integer(kind=ik), pointer    , Dimension(:)     :: pop_age
    
    Character(Len=:), allocatable, Dimension(:)     :: transpr_age_gr
    Character(Len=:), allocatable, Dimension(:)     :: transpr_sex

    Integer(kind=ik), pointer    , Dimension(:)     :: connect_work_distid

    Integer(Kind=ik)                                :: time_n, start_time
    Logical                                         :: cp_write
    Integer(kind=ik)                                :: cp_time
    Integer                                         :: cp_unit
    
    Logical                                         :: cp_reload
    Integer(kind=ik)                                :: cp_reload_time
    Integer                                         :: cp_reload_unit
    character(len=:), Allocatable                   :: cp_dir, cp_file
    
    character(len=8)                                :: cp_time_char
    Logical                                         :: cp_exist, cp_unit_open 
    
    Integer                                         :: std_int
    Integer                                         :: int_storage_size
    Integer(kind=1)                                 :: std_int1

    Integer                                         :: int1_storage_size,ret
    Integer(kind=8)                                 :: wpos
    Character(len=:), allocatable                   :: exec_type
    !===========================================================================
    ! Implementation
    !===========================================================================

    Call pt_get("#exec.procedure",exec_type)
    !** Init proc character to name files --------
    if(trim(exec_type) .eq. "Optimization") then
      write(proc_char,'("_",I7.7)')ROOT
    else
      write(proc_char,'("_",I7.7)')proc
    endif

    !** Init default integer storage size --------
    int_storage_size  = storage_size(std_int)  / FILE_STORAGE_SIZE
    int1_storage_size = storage_size(std_int1) / FILE_STORAGE_SIZE

    !** Get time_n from R0change -----------------
    time_n = Maxval(R0change) + 1

    !** If Checkpoint should be loaded -----------
    if (trim(exec_type) .eq. "Optimization") then 
       cp_reload = sim_switch%sim_reload
       cp_reload_time = sim_switch%read_week
    else
      call pt_get("#cp.reload"      , cp_reload)
      if (cp_reload) then
         call pt_get("#cp.reload.time" , cp_reload_time)
      end if
    endif
    
    !** Shift seed_date --------------------------
    days             = 1
    seed_date_mod        = add_date(seed_date,days)

    !** Generate seed_sequence -------------------
    days             = -1-seed_before
    seed_before_char = add_date(seed_date_mod,days)

    seed_seq         = generate_seq(seed_before_char,seed_date_mod)

    if (PT_DEBUG) then
       write(un_lf,PTF_sep)
       write(un_lf,PTF_M_A)"Seed sequence for ill cases:",seed_seq
    End if

    !** Derive dates of infections for those that are inf_cont, but --
    !** are not yet aware about it (will be registered the next     --
    !** two days)                                                   --
    seed_temp          = add_date(seed_date_mod,cont_dur)
    seed_inf_cont_seq  = generate_seq(add_date(seed_date_mod,1),seed_temp)

    if (PT_DEBUG) then
       write(un_lf,PTF_sep)
       write(un_lf,PTF_M_A)"Seed sequence for infected contagious cases:",seed_inf_cont_seq
    End if

    !** Derive dates of infections for those that are inf_cont, but --
    !** are not yet aware about it (will be registered the next     --
    !** 3-5 days)                                                   --
    seed_inf_ncont_seq = generate_seq(add_date(seed_date_mod,cont_dur+1),add_date(seed_date_mod,inf_dur+cont_dur))

    if (PT_DEBUG) then
       write(un_lf,PTF_sep)
       write(un_lf,PTF_M_A)"Seed sequence for infected non-contagious cases:",seed_inf_ncont_seq
    End if

    !** ------------------------------------------------------------------------
    !** Init ICU risk per age, sex and duration of illness ---------------------
    call init_ICU_risk(control_age_sex, iol, ill_dur, icu_per_day, ICU_risk_pasd)

    !** ------------------------------------------------------------------------
    !** Init survival chance of ill state by age and sex  ----------------------
    call init_surv_ill(control_age_sex, iol, ill_dur, surv_ill_pas)

    !** ------------------------------------------------------------------------
    !** Init survival chance of icu state by age and sex  ----------------------
    call init_surv_icu(control_age_sex, iol, ill_dur, surv_icu_pas)

    counties_index = get_int_table_column(iol%counties,'dist_id')

    num_counties = Size(counties_index)

    n_iter = iter_e - iter_s + 1

    if (cp_reload) then
       start_time = cp_reload_time
    else
       start_time = 1
    end if
      
    Allocate(healthy_cases_final   (num_counties,time_n,n_iter))
    Allocate(inf_noncon_cases_final(num_counties,time_n,n_iter))
    Allocate(inf_contag_cases_final(num_counties,time_n,n_iter))
    Allocate(ill_contag_cases_final(num_counties,time_n,n_iter))
    Allocate(ill_ICU_cases_final   (num_counties,time_n,n_iter))
    Allocate(immune_cases_final    (num_counties,time_n,n_iter))
    Allocate(dead_cases_final      (num_counties,time_n,n_iter))
         
    inf_seed_date = get_char_column(iol%seed,'date')

    max_date = find_max_date(inf_seed_date)

    !** Allocate and setup dist_id cross reference -----------------------------
    allocate(dist_id_cref(minval(counties_index):maxval(counties_index)))
    
    dist_id_cref = -1

    Do ii = 1, num_counties
       dist_id_cref(counties_index(ii)) = ii
    End Do

    l_seed_date = seed_date_mod

    !!=======================================================================
    !! seed infections

    inf_seed_distid => get_int_column_pointer(iol%seed,'distid')
    inf_seed_cases  => get_int_column_pointer(iol%seed,'cases')
    
    temp = get_index(inf_seed_date,seed_seq)
    seed_ill%dist_id  = inf_seed_distid(temp)
    seed_ill%date     = inf_seed_date(temp)
    seed_ill%cases    = inf_seed_cases(temp)

    if (PT_DEBUG) then
       write(un_lf,*)"sum(seed_ill%cases):",sum(seed_ill%cases),size(seed_ill%cases)       
    End if

    temp = get_index(inf_seed_date,seed_inf_cont_seq)

    seed_inf_cont%dist_id  = inf_seed_distid(temp)
    seed_inf_cont%date     = inf_seed_date(temp)
    seed_inf_cont%cases    = inf_seed_cases(temp)

    if (PT_DEBUG) then
       write(un_lf,*)"sum(seed_inf_cont%cases):",sum(seed_inf_cont%cases)
    end if
    
    temp = get_index(inf_seed_date,seed_inf_ncont_seq)
    seed_inf_ncont%dist_id  = inf_seed_distid(temp)
    seed_inf_ncont%date     = inf_seed_date(temp)
    seed_inf_ncont%cases    = inf_seed_cases(temp)

    if (PT_DEBUG) then
       write(un_lf,*)"sum(seed_inf_ncont%cases):",sum(seed_inf_ncont%cases)
    end if
    
    dea_seed_date   =  get_char_column(iol%death,'date')
    dea_seed_distid => get_int_column_pointer(iol%death,"distid")
    dea_seed_cases  => get_int_column_pointer(iol%death,'deaths')

    temp_date  = get_start_date(dea_seed_date)
    
    seed_d_seq = generate_seq(temp_date,add_date(l_seed_date,-1))
    temp       = get_index(dea_seed_date,seed_d_seq)

    seed_death%dist_id = dea_seed_distid(temp)
    seed_death%date    = dea_seed_date(temp)
    seed_death%cases   = dea_seed_cases(temp)

    If (Allocated(seed_ill_dur)) Deallocate(seed_ill_dur)

    Allocate(seed_ill_dur(Size(seed_ill%date)))

    days             = 0

    l_seed_date      = add_date(l_seed_date,days)

    temp_int = Date2Unixtime(l_seed_date)
    ii=1

    Do i = 1,Size(seed_ill_dur)
       seed_ill_dur(i)  = (temp_int - Date2Unixtime(seed_ill%date(i)))/86400 + 1
    End Do

    temp_int = 1

    temp = condition_and(seed_ill_dur,ill_dur+1,"l",seed_ill%cases,temp_int,"g")

    seed_ill%dist_id  = seed_ill%dist_id(temp)
    seed_ill%date     = seed_ill%date(temp)
    seed_ill%cases    = seed_ill%cases(temp)
    seed_ill_dur      = seed_ill_dur(temp)

    if (PT_DEBUG) then
       write(un_lf,*)"sum(seed_ill%cases):",sum(seed_ill%cases)
    end if
    
    if (PT_DEBUG) then
       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"First 20 and last 20 seeds for ill cases."
       Do ii=1, 20
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_ill%dist_id(ii)    , seed_ill%date(ii), &
               seed_ill%cases(ii),seed_ill_dur(ii)
       End Do
       write(un_lf,"('...')")
       Do ii = Size(temp)-19, Size(temp)
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_ill%dist_id(ii)    , seed_ill%date(ii), &
               seed_ill%cases(ii),seed_ill_dur(ii)
       End Do
    End if

    If (.Not.Allocated(seed_inf_cont_dur))Then
       Allocate(seed_inf_cont_dur(Size(seed_inf_cont%date)))
    End If

    days             = -1

    l_seed_date      = add_date(l_seed_date,days)

    temp_int = Date2Unixtime(l_seed_date)

    Do i = 1,Size(seed_inf_cont_dur)
       seed_inf_cont_dur(i)  = (temp_int - Date2Unixtime(seed_inf_cont%date(i)))/86400&
            + cont_dur + 2
    End Do

    if (allocated(temp)) Deallocate(temp)
    Allocate(temp(Size(seed_inf_cont%cases)))

    temp = 0
    temp_int = 0
    Do i = 1,Size(temp)
       If (seed_inf_cont%cases(i)>0) Then
          temp_int = temp_int + 1
          temp(temp_int) = i
       End If
    End Do

    seed_inf_cont%dist_id  = seed_inf_cont%dist_id(temp(1:temp_int))
    seed_inf_cont%date     = seed_inf_cont%date(temp(1:temp_int))
    seed_inf_cont%cases    = seed_inf_cont%cases(temp(1:temp_int))

    if (PT_DEBUG) then
       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"First 20 and last 20 seeds for contagious cases."
       Do ii=1, 20
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_inf_cont%dist_id(ii)    , seed_inf_cont%date(ii), &
               seed_inf_cont%cases(ii),seed_inf_cont_dur(ii)
       End Do
       write(un_lf,"('...')")
       Do ii = Size(seed_inf_cont%dist_id)-19, Size(seed_inf_cont%dist_id)
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_inf_cont%dist_id(ii)    , seed_inf_cont%date(ii), &
               seed_inf_cont%cases(ii),seed_inf_cont_dur(ii)
       End Do
    End if

    Deallocate(temp)

    If (.Not.Allocated(seed_inf_ncont_dur))Then
       Allocate(seed_inf_ncont_dur(Size(seed_inf_ncont%date)))
    Else
       Deallocate(seed_inf_ncont_dur)
       Allocate(seed_inf_ncont_dur(Size(seed_inf_ncont%date)))
    End If
    temp_int = Date2Unixtime(l_seed_date)
    Do i = 1,Size(seed_inf_ncont_dur)
       seed_inf_ncont_dur(i)  = (temp_int - Date2Unixtime(seed_inf_ncont%date(i)))/86400&
            + inf_dur + cont_dur + 2
    End Do

    Allocate(temp(Size(seed_inf_ncont%cases)))
    temp = 0
    temp_int = 0
    Do i = 1,Size(temp)
       If (seed_inf_ncont%cases(i)>0) Then
          temp_int = temp_int + 1
          temp(temp_int) = i
       End If
    End Do

    seed_inf_ncont%dist_id  = seed_inf_ncont%dist_id(temp(1:temp_int))
    seed_inf_ncont%date     = seed_inf_ncont%date(temp(1:temp_int))
    seed_inf_ncont%cases    = seed_inf_ncont%cases(temp(1:temp_int))

    if (PT_DEBUG) then
       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"First 20 and last 20 seeds for non contagious cases."
       Do ii=1, 20
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_inf_ncont%dist_id(ii)    , seed_inf_ncont%date(ii), &
               seed_inf_ncont%cases(ii),seed_inf_ncont_dur(ii)
       End Do
       write(un_lf,"('...')")
       Do ii = Size(seed_inf_ncont%dist_id)-19, Size(seed_inf_ncont%dist_id)
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_inf_ncont%dist_id(ii)    , seed_inf_ncont%date(ii), &
               seed_inf_ncont%cases(ii),seed_inf_ncont_dur(ii)
       End Do

       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"Seeds per county."
       write(un_lf,'(5(A10,1X))')"county","inf_ncont","inf_cont","inf_ill","inf_dth"
    End if

    !! ----------------------------------------------------------------------
    !! Convert from Weekly to daily R0_effects ------------------------------
    R0_daily = R0_force   * R0 / Real(Real(cont_dur)+Real(ill_dur)*less_contagious) + &
         ( 1 - R0_force ) * R0 / Real(cont_dur+ill_dur)
    
    ! this block simplifies the if judgment
    If (.Not.Allocated(R0matrix))Then
       Allocate(R0matrix(num_counties,time_n-1))
    End If
    R0matrix = R0_daily
    n_change  = Size(R0change,dim=2)
    
    Do i = 1,n_change
       
       getchange = R0_effects(i,region_index)
       
       gettime = generate_seq(R0change(1,i),R0change(2,i),1)
       
       Do j = 1,Size(gettime)
          R0matrix(:,gettime(j)) = R0matrix(:,gettime(j)) * getchange
       End Do
    End Do
    
    If (R0delay) Then
       Do i = 1,Size(R0matrix,dim = 1)
          R0matrix(i,:) = smoothing_change(R0matrix(i,:),R0delay_days,R0delay_type)
       End Do
    End If
    
    !! ==============================================================
    !! Init connectivity matrix
    connect =  transpose(table_to_real_array(iol%connect_work))

    Do i = 1,Size(connect,dim=2)
       connect(:,i) = connect(:,i)/Sum(connect(:,i))
    End Do

    connect = transpose(connect)

!!!=============================================================================
!!! Open checkpoint restart file ===============================================
!!!=============================================================================
    if (trim(exec_type) .eq. "Optimization") then 
       cp_write = sim_switch%sim_write
    else
       call pt_get("#cp.write" , cp_write)
    endif
    if ( cp_write ) then
       if (trim(exec_type) .eq. "Optimization") then 
          cp_time  = time_n ! Write restart file at the last timestep
          write(cp_time_char,'("_",I3.3,"week")')sim_switch%write_week ! The restart file is written in weeks
       else
          call pt_get("#cp.time" , cp_time)
          write(cp_time_char,'("_",I7.7)')cp_time
       endif
       cp_file = trim(output_dir)//"/"//"CoSMic"//cp_time_char//proc_char//".restart"

       ! If we already have a restart file at the given timestep -----
       inquire(exist=cp_exist, file=trim(cp_file))

       cp_unit_open = .FALSE.
       
       if (cp_exist) then
          cp_write = .FALSE.
          write(un_lf,'(A)')"Warning - Checkpoint already exists."
          write(un_lf,'(A)')trim(output_dir)//"/"//"CoSMic"//cp_time_char//proc_char//".restart"
          write(un_lf,'(A)')"Reset cp_write = .FALSE."
       else
          Open(newunit=cp_unit, action="write", status="replace", access="stream", &
               file=trim(cp_file))
          cp_unit_open =.TRUE.
       end if
       
    End if

    if ( cp_reload ) then

       call pt_get("#cp.dir"         , cp_dir)

       if (trim(exec_type) .eq. "Optimization") then 
          write(cp_time_char,'("_",I3.3,"week")')cp_reload_time
       else
          write(cp_time_char,'("_",I7.7)')cp_reload_time
       endif
       Open(newunit=cp_reload_unit, action="read", status="old", access="stream", &
            file=trim(cp_dir)//"/"//"CoSMic"//cp_time_char//proc_char//".restart")

       Write(un_lf,PTF_M_A)"Opend checkpoint file:",&
            trim(cp_dir)//"/"//"CoSMic"//cp_time_char//proc_char//".restart"
    End if      
    
!!!=============================================================================
!!! Non-Deterministic Iteration
!!!=============================================================================

    !!==============================================================
    !! Init population
    pop_total  => get_int_column_pointer(iol%pop,'total')
    pop_distid => get_int_column_pointer(iol%pop,'dist_id')
    pop_sex    =  get_char_column(iol%pop,'sex')
    pop_age    => get_int_column_pointer(iol%pop,'age_gr')

    pop_total = Nint(Real(pop_total)/Real(Sum(pop_total)) * sam_size)
    pop_size      = Sum(pop_total)

    !$OMP PARALLEL if((iter_e-iter_s+1) .gt. 1)&
    !$OMP& default(shared) &
    !$OMP& private(sim) &
    !$OMP& firstprivate(tmp_d_new) &
    !$OMP& firstprivate(index,temp_int,i,temp,days,icounty,county,jj,kk,ii,tmp_count) &
    !$OMP& firstprivate(seed_ill,seed_inf_cont,seed_inf_ncont,seed_death,seed_d_seq,temp_date) &
    !$OMP& firstPRIVATE(il_d,inf_ill,inf_c_d,inf_cont,inf_nc_d,inf_ncont,inf_dth) &
    !$OMP& firstPRIVATE(rownumbers_left,rownumbers) &
    !$OMP& firstPRIVATE(rownumbers_ill,rownumbers_cont,rownumbers_ncont,rownumbers_dea) &
    !$OMP& firstprivate(seed_ill_dur,seed_inf_cont_dur,seed_inf_ncont_dur,l_seed_date) 
    !$OMP DO
    Do it_ss = 1, (iter_e-iter_s+1)

       !call start_timer("Init Sim Loop",reset=.FALSE.)
       !call random_seed(put=(/1,2,3,4,5,6,7,8/))  ! Put seed array into PRNG.
       call random_seed(size=ii)               ! Get size of seed array.
       ii = max(ii,OMP_GET_MAX_THREADS())
       ur = urandom_seed(ii)
       !write(un_lf,*)ii,proc,ur
       call random_seed(put=ur)  ! Put seed array into PRNG.
       
       If (.Not.Allocated(sim%dist_id))Then
          Allocate(sim%dist_id(pop_size))
          Allocate(sim%sample_index(pop_size))
          Allocate(sim%sex(pop_size))
          Allocate(sim%age(pop_size))
          Allocate(sim%t1(pop_size))
          Allocate(sim%t2(pop_size))
          Allocate(sim%d(pop_size))
       Endif

       !$OMP critical
        If ((OMP_GET_THREAD_NUM() == 0) .and. (.Not.Allocated(sim_backup%dist_id))) Then
          Allocate(sim_backup%dist_id(pop_size))
          Allocate(sim_backup%sample_index(pop_size))
          Allocate(sim_backup%sex(pop_size))
          Allocate(sim_backup%age(pop_size))
          Allocate(sim_backup%t1(pop_size))
          Allocate(sim_backup%t2(pop_size))
          Allocate(sim_backup%d(pop_size))
       Endif
       !$OMP END critical

     If (.not.(sim_switch%sim_reuse)) then
       index = 0 ! position index

       Do i = 1, Size(pop_total)
          temp_int = pop_total(i)
          sim%dist_id(index+1: index+temp_int) = Int(pop_distid(i),2)
          sim%sex(index+1: index+temp_int)     = pop_sex(i)
          sim%age( index+1: index+temp_int)    = int(pop_age(i),1)
          index                                = index + temp_int
       End Do

       sim%t1 = healthy
       sim%t2 = missing
       sim%d(:)  = 1

     write(un_lf,PTF_SEP)
     Write(un_lf,PTF_M_AI0)"Size of population is", pop_size

     !** Reshuffle population since ordering according to dist id leads to ***
     !** lower infections in older age groups                              ***
     !** To do so we use
       do ii = 1, pop_size
          sim%d(ii) = ii
       End do

       Allocate(tmp_count(pop_size))

       if ( cp_reload ) then
          if ( trim(exec_type) .eq. "Optimization" ) then 
            wpos = int(int_storage_size,8)*int(pop_size,8)*int(ROOT,8)+1_8
          else
            wpos = int(int_storage_size,8)*int(pop_size,8)*int(OMP_GET_THREAD_NUM(),8)+1_8
          endif
          !call start_timer("cp_reload tmp_count",reset=.FALSE.)
          
          !write(un_lf,PTF_M_AI0)"Reading index shuffle at position:",wpos,OMP_GET_THREAD_NUM()
          read(cp_reload_unit,pos=wpos)tmp_count
          !call end_timer("cp_reload tmp_count")
       Else
          tmp_count = sample(sim%d,pop_size)
       End if
       sim%sample_index = tmp_count

!!!=============================================================================
!!! Write index shuffle for restart since this is non deterministic ============
!!!=============================================================================
       if ( cp_write ) then
          if ( (trim(exec_type) .eq. "Optimization") &
           .and. (OMP_GET_THREAD_NUM() .eq. 0) ) then
            wpos = int(int_storage_size,8)*int(pop_size,8)*int(OMP_GET_THREAD_NUM(),8)+1_8
            write(cp_unit,pos=wpos)tmp_count
          else
            !$OMP critical
            !call start_timer("cp_write",reset=.FALSE.)
            wpos = int(int_storage_size,8)*int(pop_size,8)*int(OMP_GET_THREAD_NUM(),8)+1_8
            !write(un_lf,PTF_M_AI0)"Writing index shuffle at position:",wpos
            write(cp_unit,pos=wpos)tmp_count
            flush(cp_unit)
            ret = fsync(fnum(cp_unit))
            if (ret /= 0) then
               write(*,*) "ret = fsync(fnum(cp_unit)) = ",ret
               stop
            End if
            !call end_timer("cp_write")
            !$OMP end critical
          endif 
       End if
       
       sim%dist_id=sim%dist_id(tmp_count)
       sim%sex    =sim%sex    (tmp_count)
       sim%age    =sim%age    (tmp_count)
       sim%d  = 1

       deallocate(tmp_count)

       if ( cp_reload ) then

          !write(un_lf,*)"Reloading checkpoint at time:",cp_reload_time
          !write(un_lf,PTF_M_AI0)"Reading sim%t1 at position:",&
          !     int(int_storage_size,8)  * int(pop_size,8) * int(OMP_GET_NUM_THREADS(),8) + &
          !     int(int1_storage_size,8) * int(pop_size,8) * int(OMP_GET_THREAD_NUM() ,8) + 1_8,&
          !     OMP_GET_THREAD_NUM()
          !call start_timer("cp_reload sim",reset=.FALSE.) 
          if ( trim(exec_type) .eq. "Optimization" ) then
            read(cp_reload_unit,  &
                 pos = int(int_storage_size,8)  * int(pop_size,8) + &
                       int(int1_storage_size,8) * int(pop_size,8) * int(ROOT ,8) + 1_8)sim%t1
            read(cp_reload_unit, &
                 pos = int(int_storage_size,8)  * int(pop_size,8) + &
                       int(int1_storage_size,8) * int(pop_size,8) + &
                       int(int_storage_size,8)  * int(pop_size,8) * int(ROOT ,8) + 1_8)sim%d
          else
            read(cp_reload_unit,  &
                 pos = int(int_storage_size,8)  * int(pop_size,8) * int(OMP_GET_NUM_THREADS(),8) + &
                       int(int1_storage_size,8) * int(pop_size,8) * int(OMP_GET_THREAD_NUM() ,8) + 1_8)sim%t1

            !write(un_lf,PTF_M_AI0)"Reading sim%d at position:",&
            !     int(int_storage_size,8)  * int(pop_size,8) * int(OMP_GET_NUM_THREADS(),8) + &
            !     int(int1_storage_size,8) * int(pop_size,8) * int(OMP_GET_NUM_THREADS(),8) + &
            !     int(int_storage_size,8)  * int(pop_size,8) * int(OMP_GET_THREAD_NUM() ,8) + 1_8,&
            !     OMP_GET_THREAD_NUM()
            
            read(cp_reload_unit, &
                 pos = int(int_storage_size,8)  * int(pop_size,8) * int(OMP_GET_NUM_THREADS(),8) + &
                       int(int1_storage_size,8) * int(pop_size,8) * int(OMP_GET_NUM_THREADS(),8) + &
                       int(int_storage_size,8)  * int(pop_size,8) * int(OMP_GET_THREAD_NUM() ,8) + 1_8)sim%d
          endif
          !call end_timer("cp_reload sim")
       Else
          
          Do icounty = 1,num_counties

             county = counties_index(icounty)

             rownumbers = get_index(sim%dist_id,county)

             temp   = get_index(seed_ill%dist_id,county)
             il_d   = rep(seed_ill_dur(temp),seed_ill%cases(temp))
             inf_ill= Sum(seed_ill%cases(temp))

             temp   = get_index(seed_inf_cont%dist_id,county)
             inf_c_d= rep(seed_inf_cont_dur(temp),seed_inf_cont%cases(temp))

             inf_cont = Sum(seed_inf_cont%cases(temp))

             temp   = get_index(seed_inf_ncont%dist_id,county)
             inf_nc_d = rep(seed_inf_ncont_dur(temp),seed_inf_ncont%cases(temp))
             inf_ncont = Sum(seed_inf_ncont%cases(temp))

             temp   = get_index(seed_death%dist_id,county)
             inf_dth = Sum(seed_death%cases(temp))

             if (PT_DEBUG) then
                write(un_lf,'(11(I11))')county,inf_ncont,inf_cont,inf_ill,inf_dth,&
                     minval(inf_nc_d),maxval(inf_nc_d),minval(inf_c_d),maxval(inf_c_d),&
                     minval(il_d),maxval(il_d)
             End if

             If(Size(rownumbers)<(inf_ill+inf_cont+inf_ncont+inf_dth)) Then
                Print *,"Number of infected and dead is larger than population size"
                Print *,"only ",Size(rownumbers),"number left"
                Print *,"total cases is",(inf_ill+inf_cont+inf_ncont+inf_dth)
                Print *,"seperate cases are:" ,inf_ill,inf_cont,inf_ncont,inf_dth
                Print *,"timestep is",timestep
                Print *,"it_ss is",it_ss
                Print *,"county is",county
                Print *,"counties are",counties_index
                Print *,"icounty is",icounty
                Call Exit(status)
             End If

             rownumbers_left = rownumbers
             If (inf_ill > 0) Then

                rownumbers_ill = sample(rownumbers_left,inf_ill)
                Call QSortC(rownumbers_ill)

                jj = 1
                kk = 1
                Do ii = 1, Size(rownumbers_left)

                   if (rownumbers_left(ii) == rownumbers_ill(jj)) then
                      if (jj < inf_ill) then
                         jj = jj + 1
                      End if
                   Else
                      rownumbers_left(kk) = rownumbers_left(ii)
                      kk = kk + 1
                   End if

                End Do

                sim%t1(rownumbers_ill) = ill_contag
                sim%d(rownumbers_ill)  = il_d
             End If

             If ( inf_cont > 0) Then
                rownumbers_cont = sample(rownumbers_left(1:(Size(rownumbers)-inf_ill)),inf_cont)

                Call QSortC(rownumbers_cont)             

                jj = 1
                kk = 1
                Do ii = 1, Size(rownumbers_left)

                   if (rownumbers_left(ii) == rownumbers_cont(jj)) then
                      if (jj < inf_cont) then
                         jj = jj + 1
                      End if
                   Else
                      rownumbers_left(kk) = rownumbers_left(ii)
                      kk = kk + 1
                   End if

                End Do

                sim%t1(rownumbers_cont)= inf_contag
                sim%d(rownumbers_cont) = inf_c_d
             End If

             If (inf_ncont > 0) Then
                rownumbers_ncont = sample(rownumbers_left(1:(Size(rownumbers)-inf_ill-inf_cont)),inf_ncont)

                Call QSortC(rownumbers_ncont)


                jj = 1
                kk = 1
                Do ii = 1, Size(rownumbers_left)

                   if (rownumbers_left(ii) == rownumbers_ncont(jj)) then
                      if (jj < inf_ncont) then
                         jj = jj + 1
                      End if
                   Else
                      rownumbers_left(kk) = rownumbers_left(ii)
                      kk = kk + 1
                   End if

                End Do

                sim%t1(rownumbers_ncont) = inf_noncon
                sim%d(rownumbers_ncont)  = inf_nc_d
             End If

             If (inf_dth > 0) Then
                rownumbers_dea = sample(rownumbers_left(1:(Size(rownumbers)-inf_ill-inf_cont-inf_ncont)),inf_dth)

                sim%t1(rownumbers_dea) = dead
             End If

          End Do ! do icounty = 1,size(sim_counties)

       end if
      else ! Reuse the passed sim data
        sim%dist_id      = sim_backup%dist_id
        sim%sample_index = sim_backup%sample_index
        sim%sex          = sim_backup%sex 
        sim%age          = sim_backup%age
        sim%t1           = sim_backup%t1  
        sim%d            = sim_backup%d
        if ( cp_write ) then
          if ( (trim(exec_type) .eq. "Optimization") &
           .and. (OMP_GET_THREAD_NUM() .eq. 0) ) then
            wpos = int(int_storage_size,8)*int(pop_size,8)*int(OMP_GET_THREAD_NUM(),8)+1_8
            write(cp_unit,pos=wpos)sim%sample_index
          endif
        End if
      endif
      
       sim%t2 = missing

       !** Set up dist_id renumbered cross_reference ---------------------------
       sim%dist_id_rn = Int(dist_id_cref(sim%dist_id),2)

       If (.Not.Allocated(tmp_d_new))Then
          Allocate(tmp_d_new(Size(sim%d)))
       End If

       ! call end_timer("Init Sim Loop") 

!!!=============================================================================
!!! Simulation Loop ============================================================
!!!=============================================================================

       if (OMP_GET_THREAD_NUM() == 0) then
          call start_timer("Sim Loop",reset=.FALSE.)
       End if

       ! write(*,PTF_M_AI0)"Entering CoSMic_TimeLoop on proc",proc
       
       Call CoSMic_TimeLoop(time_n, pop_size, size(counties_index), counties_index, &
            Real(R0matrix,rk), connect, surv_ill_pas, ICU_risk_pasd, surv_icu_pas, sim ,&
            healthy_cases_final, inf_noncon_cases_final,inf_contag_cases_final, &
            ill_contag_cases_final, ill_ICU_cases_final,  &
            immune_cases_final,dead_cases_final, &
            it_ss, cp_write, cp_time , cp_unit, cp_reload, cp_reload_time)

       !** Return the ultimate sim data to the caller --------------------------
      
       !$OMP critical
       ! when there is only one single thread
       if(OMP_GET_THREAD_NUM() == 0 .and. (it_ss .eq. (iter_e-iter_s+1))) then
        sim_backup%dist_id      = sim%dist_id
        sim_backup%sample_index = sim%sample_index
        sim_backup%sex          = sim%sex 
        sim_backup%age          = sim%age
        sim_backup%t1           = sim%t1
        sim_backup%t2           = sim%t2 
        sim_backup%d            = sim%d
       endif
       !$OMP end critical
     
       if (OMP_GET_THREAD_NUM() == 0) then
          call end_timer("Sim Loop")

          !timer = get_timer("Sim Loop")

          !write(un_lf,'(A)',ADVANCE="NO")"Time per day:"
          !call write_realtime(frac_realtime(diff_realtimes(timer%rt_end,timer%rt_start),time_n),un_lf)
       End if

       days             = +1

       l_seed_date        = add_date(l_seed_date,days)
        
    End Do     ! end do it_ss
    !$OMP END DO
    !$OMP END PARALLEL

    ! Close restart files ------------------------
    if ( cp_reload ) then
       close(cp_reload_unit)
    End if
    
    if ( cp_write .AND. cp_unit_open ) then
       close(cp_unit)
    End if
  !  if ( cp_reload ) then 
  !     close(cp_reload_unit)
  !  End if
    
    !** Disable the writeout operation for the execution with "Optimization" type ------------------
    if(trim(exec_type) .ne. "Optimization") then  
      call start_timer("+- Writeout",reset=.FALSE.)
      
      call write_data(trim(output_dir)//"healthy_cases_"//trim(export_name)//trim(proc_char)//".dat", &
           healthy_cases_final,R0_effects)                                    
      call write_data(trim(output_dir)//"inf_noncon_cases_"//trim(export_name)//trim(proc_char)//".dat", &
           inf_noncon_cases_final,R0_effects)
      call write_data(trim(output_dir)//"inf_contag_cases_"//trim(export_name)//trim(proc_char)//".dat", &
           inf_contag_cases_final,R0_effects)
      call write_data(trim(output_dir)//"ill_contag_cases_"//trim(export_name)//trim(proc_char)//".dat", &
           ill_contag_cases_final,R0_effects)
      call write_data(trim(output_dir)//"ill_ICU_cases_"//trim(export_name)//trim(proc_char)//".dat", &
           ill_ICU_cases_final,R0_effects)
      call write_data(trim(output_dir)//"dead_cases_"//trim(export_name)//trim(proc_char)//".dat", &
           dead_cases_final,R0_effects)
      call write_data(trim(output_dir)//"immune_cases_"//trim(export_name)//trim(proc_char)//".dat", &
           immune_cases_final,R0_effects)
      
      call end_timer("+- Writeout")
    endif

  End Function COVID19_Spatial_Microsimulation_for_Germany

  Subroutine write_data(filename, data, R0_effects)

    character(len=*)                      , intent(in)        :: filename
    Integer, Dimension(:,:,:), Allocatable, intent(in)        :: data
    Real(kind=rk),    Allocatable, Dimension(:,:),intent(in)  :: R0_effects
    
    logical                                            :: exs
    Integer                                            :: un, i,j
 
    !---------------------------------------------------------------------------
    exs = .FALSE.
    
    Inquire(file = trim(filename), exist=exs)
    
    If (exs) then

       Open (newunit=un, file = trim(filename)//".head",&
            position="Append", status ="old")
       write(un,'(A,3(",",I0))',Advance="NO")"int",shape(data)
       write(un,*)(",",R0_effects(1,j),&
            (",",R0_effects(i,j), i=2, size(R0_effects,dim=1)),&
            j=1,size(R0_effects,dim=2))
       close(un)

       Open (newunit=un, file = trim(filename),&
            position="Append", status ="old", access="stream")
    else

       Open (newunit=un, file = trim(filename)//".head",&
            status ="new")
       write(un,'(A)',ADVANCE="NO")"type , dim1 , dim2 , dim3 "
       write(un,'(*(",",A,I3.3))',ADVANCE="NO")("R",i, i=1, size(R0_effects)-1)
       write(un,'(A,I0)')"R",size(R0_effects)

       write(un,'(A,3(",",I0))',Advance="NO")"int",shape(data)
       write(un,*)(",",R0_effects(1,j),&
            (",",R0_effects(i,j), i=2, size(R0_effects,dim=1)),&
            j=1,size(R0_effects,dim=2))
       
       close(un)

       Open (newunit=un, file = trim(filename),&
            access="stream", status="new")
    End If
    
    write(un)data
    
    Close(un)
    
  End Subroutine write_data

  Subroutine write_data_v2(healthy_cases_final,iter_pass_handle,R0Change,counties_index,&
       type_file,filename,proc,iter_s,iter_e)
    Real(kind=rk),Dimension(:)             :: iter_pass_handle(:)
    Real(kind=rk),Dimension(:,:)  :: R0Change
    Real(kind=rk),Allocatable              :: R0change_rep(:,:),R0change_exp(:,:),temp_output(:)
    Integer,Dimension(:,:,:)      :: healthy_cases_final
    Integer,Dimension(:)          :: counties_index
    Integer                       :: iter_s,iter_e, n_iter
    Integer                       :: type_file, proc
    Character(Len=*), intent(in)  :: filename
    
    Integer :: county_size, count 
    Character*15                :: iter_char(6)
    Character*8                 :: proc_char
    Character*5,Allocatable     :: R0change_name(:)
    Character*2                 :: counties(16)
    Integer,Allocatable         :: counties_index_out(:),iter_array(:)

    Character*10,Allocatable    :: Label(:)

    Integer date_time(8),i,j
    Character*10 b(3)
    Character*4 year
    Character*2 day,month
    Character*8 time
    Character*3 temp_char

    Character*15 dir_prefix

    logical :: exs
    
    n_iter      = iter_e-iter_s+1
    county_size = Size(healthy_cases_final,DIM=1)

    iter_char = (/"sam_size       ",&
         "R0             ",&
         "icu_dur        ",&
         "w_int          ",&
         "w.obs          ",&
         "w.obs.by.sate  "/)

    counties  = (/"SH","HH","NI","HB","NW","HE","RP","BW","BY","SL",&
         "BE","BB","MV","SN","ST","TH"/)

    write(proc_char,'("_",I7.7)')proc
    
    Allocate(R0change_rep(Size(R0Change),1))
    Allocate(R0change_exp(n_iter*county_size,Size(R0change_rep)))
    Allocate(temp_output(n_iter*county_size))
    Allocate(R0change_name(Size(R0Change)))
    Allocate(iter_array(county_size*n_iter))
    Allocate(counties_index_out(county_size*n_iter))
    Allocate(Label(Size(healthy_cases_final,dim=2)))

    Call date_and_Time(b(1),b(2),b(3),date_time)
    
    count = 1
    Do i = 1,Size(counties)
       Do j = 1,Size(R0Change,1)
          Write(temp_char,"(I2)")j
          If (j>9)Then
             R0change_name(count) = counties(i)//temp_char
          Else
             R0change_name(count) = counties(i)//Adjustl(temp_char)
          End If
          count = count + 1
       End Do
    End Do
    count = 1
    Do i =1,Size(healthy_cases_final,dim =2)
       Write(Label(i),"(I3)")i
    End Do

    Write(year,"(I4)")date_time(1)
    Write(month,"(I2)")date_time(2)
    Write(day,"(I2)")date_time(3)
    If (date_time(2)<=9)Then
       month = "0"//Adjustl(month)
    Endif

    dir_prefix = "./output/"

    If(date_time(3)<=9)Then
       day = "0"//Adjustl(day)
    End If
    time = year//month//day
    
    exs = .FALSE.
    Inquire(file = Trim(dir_prefix)//time//proc_char//trim(filename), exist=exs)
    If (exs) then
       Open (101, file = Trim(dir_prefix)//time//proc_char//trim(filename),&
            position="Append", status ="old")
    else
       Open (101, file = Trim(dir_prefix)//time//proc_char//trim(filename),&
            status ="new")
       !write the head file
       Write(101,10)Label,iter_char,R0change_name,"iter","x.dist_id"
    End If
    
    R0change_rep = Reshape(R0Change,(/Size(R0change),1/))
    R0change_exp = Transpose(Spread(R0change_rep(:,1),2,n_iter*county_size))

    count = 1
    Do i = 1,n_iter
       Do j = 1,county_size
          Write(101,10)healthy_cases_final(j,:,i),iter_pass_handle,R0change_exp(count,:),&
               (i-1+iter_s),counties_index(j)
          count = count+1
       End Do
    End Do

    Close(101)
    
10  Format(1x,*(g0,","))
    
  End Subroutine write_data_v2

  !!============================================================================
  !> Model Timeloop
  Subroutine CoSMic_TimeLoop(&
       time_n, pop_size, n_counties, counties, &
       R0matrix, tmp_connect, surv_ill_pas, ICU_risk_pasd, surv_icu_pas, &
       sim, &
       healthy_cases, inf_noncon_cases, inf_contag_cases, ill_contag_cases,&
       ill_ICU_cases , immune_cases, dead_cases , &
       it_ss, cp_write, cp_time , cp_unit, cp_reload, cp_reload_time)
    
    Integer(Kind=ik)                       , Intent(In) :: time_n
    Integer(Kind=ik)                       , Intent(In) :: pop_size
    Integer(Kind=ik)                       , Intent(In) :: n_counties
    Integer(Kind=ik), Dimension(n_counties), Intent(In) :: counties
    
    Real(Kind=rk)   , Dimension(n_counties,2:time_n)  , Intent(In) ::  R0matrix
    Real(kind=rk)   , Allocatable, Dimension(:,:)     , Intent(In) ::  tmp_connect
  !  Real(Kind=rk)   , Dimension(n_counties,n_counties), Intent(In) ::  connect
    Real(Kind=rk)   , Allocatable, Dimension(:,:)     , Intent(In) ::  surv_ill_pas
    Real(Kind=rk)   , Allocatable, Dimension(:,:,:)   , Intent(In) ::  ICU_risk_pasd
    Real(Kind=rk)   , Allocatable, Dimension(:,:)     , Intent(In) ::  surv_icu_pas
    
    Type(sims)                                        , Intent(InOut) :: sim

    Integer         , Allocatable, Dimension(:,:,:), Intent(inout) :: healthy_cases
    Integer         , Allocatable, Dimension(:,:,:), Intent(inout) :: inf_noncon_cases
    Integer         , Allocatable, Dimension(:,:,:), Intent(inout) :: inf_contag_cases
    Integer         , Allocatable, Dimension(:,:,:), Intent(inout) :: ill_contag_cases
    Integer         , Allocatable, Dimension(:,:,:), Intent(inout) :: ill_ICU_cases
    Integer         , Allocatable, Dimension(:,:,:), Intent(inout) :: immune_cases
    Integer         , Allocatable, Dimension(:,:,:), Intent(inout) :: dead_cases

    Integer(Kind=ik)                               , Intent(In)    :: it_ss

    Logical                                        , Intent(In)    :: cp_write
    Integer(Kind=ik)                               , Intent(In)    :: cp_time
    Integer(Kind=ik)                               , Intent(In)    :: cp_unit
    Logical                                        , Intent(In)    :: cp_reload
    Integer(Kind=ik)                               , Intent(In)    :: cp_reload_time
    
    !> Counters -------------------------------------------------
    Integer(Kind=ik)             :: timestep, prng_ii, ii, nn
    
    Integer(Kind=ik)             :: at_risk, new_in_state
    Integer(Kind=ik)             :: start_time
    
    Integer(Kind=ik), Dimension(min_state:max_state) :: state_count_t1
    Integer(Kind=ik), Dimension(min_state:max_state) :: state_count_t2

    Integer(Kind=ik), Allocatable, Dimension(:,:)    :: state_count_pl
    Real(Kind=rk)   , Allocatable, Dimension(:)      :: risk_pl, risk_pi

    Integer(Kind=ik), Allocatable, Dimension(:)      :: d_new, ur

    Real(kind=rk)                                    :: less_contagious
    Real(kind=rk)                                    :: w_int
    Integer(kind=ik)                                 :: inf_dur, cont_dur
    Integer(kind=ik)                                 :: ill_dur, icu_dur

    Integer                                          :: std_int
    Integer                                          :: int_storage_size
    Integer(kind=1)                                  :: std_int1
    Integer                                          :: int1_storage_size

    Real(Kind=rk)   , Dimension(n_counties,n_counties) ::  connect

    Integer                                            :: len, ret

    Character(len=:), allocatable                    :: exec_type
    
    !===========================================================================
    ! Implementation
    !===========================================================================

    len = size(tmp_connect,dim=1)
    connect = 0_rk
    do ii = 1,len
      do nn = 1,len
        connect(ii,nn) = Real(tmp_connect(ii,nn),rk)
      enddo
    enddo

    !** init default integer storage size --------
    int_storage_size  = storage_size(std_int)  / FILE_STORAGE_SIZE
    int1_storage_size = storage_size(std_int1) / FILE_STORAGE_SIZE
    
    !** Allocate and init new counter ----------------------
    Allocate(state_count_pl(min_state:max_state,n_counties))
    state_count_pl = 0
    
    Allocate(risk_pl(n_counties))
    risk_pl = -1._rk

    Allocate(risk_pi(pop_size))
    risk_pi = -1._rk
        
    Allocate(d_new(pop_size))
    d_new = sim%d
    
    call pt_get("#less_contagious", less_contagious)
    call pt_get("#w_int"          , w_int          )
    call pt_get("#inf_dur"        , inf_dur        )
    call pt_get("#cont_dur"       , cont_dur       )
    call pt_get("#ill_dur"        , ill_dur        )
    call pt_get("#icu_dur"        , icu_dur        )
    call pt_get("#exec.procedure" ,exec_type       )

    ! call end_timer("+- Prepare data")

    sim%t2 = sim%t1

    !** Population summary per location ------------------------------
    Call summary_2_int(&
         sim%t1, sim%dist_id_rn, pop_size, &
         state_count_pl, min_state, max_state, 1, n_counties)

    !** Record initial state -----------------------------------------
    if (cp_reload .and. ((trim(exec_type) .ne. "Optimization"))) then
       start_time = cp_reload_time
    else
       start_time = 1
    end if
    
    healthy_cases(:,start_time,it_ss)    = state_count_pl(healthy   ,:)
    inf_noncon_cases(:,start_time,it_ss) = state_count_pl(inf_noncon,:)
    inf_contag_cases(:,start_time,it_ss) = state_count_pl(inf_contag,:)
    ill_contag_cases(:,start_time,it_ss) = state_count_pl(ill_contag,:)
    ill_ICU_cases(:,start_time,it_ss)    = state_count_pl(ill_ICU   ,:)
    dead_cases(:,start_time,it_ss)       = state_count_pl(dead      ,:)
    immune_cases(:,start_time,it_ss)     = state_count_pl(immune    ,:)

    !write(un_lf,PTF_M_AI0)"sum(state_count_pl) =",sum(state_count_pl)
    !write(*,*)"it_ss=",it_ss
    
    Do timestep = start_time+1, time_n

       !call start_timer("+- From healthy to infected",reset=.FALSE.)

!!$       if ( mod(timestep,7) == 0 ) then
!!$          write(*,*)"timestep =",timestep," re-init PRNG"
!!$          call random_seed(size=prng_ii)               ! Get size of seed array.
!!$          prng_ii = max(prng_ii,OMP_GET_MAX_THREADS())
!!$          ur = urandom_seed(prng_ii)
!!$          write(un_lf,*)omp_get_thread_num(),ur
!!$          call random_seed(put=ur)  ! Put seed array into PRNG.
!!$       end if

       
       !** Population summary --------------------------------------------------
       Call summary_1_int(sim%t1, pop_size, state_count_t1, min_state, max_state)
      
       !** DEBUG --- Population Summary ----------------------------------------
       if (PT_DEBUG) then
          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary @ start of step:",timestep
          write(un_lf,'(8(I10))')-1,0,1,2,3,4,5,6
          write(un_lf,'(8(I10))')state_count_t1

          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary per location @ start of step:",timestep
          write(un_lf,'(A6,8(I10))')"Loc",-1,0,1,2,3,4,5,6
          Do ii = 1, n_counties
             write(un_lf,'(11(I10))')ii,state_count_pl(:,ii)
          End Do
       End if
       !** DEBUG ---------------------------------------------------------------

       !** Is there something to calculate ? -----------------------------------
       If ( state_count_t1(dead) + state_count_t1(immune) == pop_size) Then
          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_A)"All are immune or dead."
          write(un_lf,PTF_SEP)
          exit
       End If

       !** Number of people who are at risk in population ---
       at_risk = state_count_t1(healthy)

       Do ii = 1, n_counties
          risk_pl(ii) = Real(state_count_pl(inf_contag,ii),rk) +  Real(state_count_pl(ill_contag,ii),rk) * less_contagious
       End Do
    !   risk_pl = state_count_pl(inf_contag,:) +  state_count_pl(ill_contag,:) * less_contagious

       risk_pl = risk_pl * R0matrix(:,timestep)

       risk_pl = w_int * risk_pl + (1._rk - w_int) * Matmul(risk_pl, connect)


       Do ii = 1, n_counties
          risk_pl(ii) = risk_pl(ii) / sum(state_count_pl([healthy,immune, &
                                                          inf_noncon,inf_contag, &
                                                          ill_contag              ],ii))
       End Do

       risk_pl = w_int * risk_pl + (1._rk - w_int) * Matmul(risk_pl, connect)
       
       !** DEBUG --- Risk per location --------------------------------------
       if (PT_DEBUG) then
          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Number of people at risk :",at_risk
          write(un_lf,PTF_M_AI0)"Risk per location @ timestep:",timestep
          write(un_lf,'(*(f0.9,1X))')risk_pl
          write(un_lf,PTF_M_AI0)"size(risk_pi) = ",size(risk_pi)
          write(un_lf,PTF_M_AI0)"at_risk       = ",at_risk
          write(un_lf,PTF_M_AI0)"pop_size      = ",pop_size
       End if
       !** DEBUG ------------------------------------------------------------

       !** Draw risks for all individuals at risk --------------------
       Call random_Number(risk_pi(1:at_risk))
       nn = 1
       Do ii = 1, pop_size

          if ((sim%t1(ii) == healthy)) then
             !** Individual is healthy ----------------------------------
             
             nn = nn + 1

             if (risk_pi(nn) <= risk_pl(sim%dist_id_rn(ii))) then
                !** Individual becomes infected ----------------------
                sim%t2(ii) = inf_noncon
                d_new(ii)  = 1
             Else
                !** Individual stays healthy -------------------------
                d_new(ii) = d_new(ii)  + 1
                
             End if
             
          Else if ((sim%t1(ii) == inf_noncon)) then
             !** Individual is already infected ----------------------
             if (sim%d(ii) >= inf_dur) then
                sim%t2(ii) = inf_contag
                d_new(ii)  = 1
             Else
                d_new(ii) = d_new(ii)  + 1
             End if
             
          Else if ((sim%t1(ii) == inf_contag)) then
             !** Individual is already infected and contagious -------
             if (sim%d(ii) >= cont_dur) then
                sim%t2(ii) = ill_contag
                d_new(ii)  = 1
             Else
                d_new(ii) = d_new(ii) + 1
             End if
          end if
                          
       End Do
       
       !** Population summary --------------------------------------------------
       Call summary_1_int(sim%t2, pop_size, state_count_t2, min_state, max_state)
       
       !** DEBUG --- Population Summary -------------------------------------
       if (PT_DEBUG) then
          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary after infection @ timestep ",timestep
          write(un_lf,'(8(I10))')-1,0,1,2,3,4,5,6
          write(un_lf,'(8(I10))')state_count_t2
       End if
       !** DEBUG ------------------------------------------------------------
       
       ! call end_timer("+- From healthy to infected")
       
       !call start_timer("+- From infected to ill and icu",reset=.FALSE.)

       at_risk = state_count_t1(ill_contag)

       if (PT_DEBUG) then
          write(un_lf,PTF_M_AI0)"At risk to die after infection:",at_risk
       End if
       
       new_in_state = 0
       
       !** Draw risks for all individuals at risk --------------------
       Call random_Number(risk_pi(1:at_risk))
       nn = 1
       Do ii = 1, pop_size

          if (sim%t1(ii) == ill_contag) then
             !** Individual was already ill_contag -------------------
             nn = nn + 1

             if (sim%sex(ii) == "m") then
                if (risk_pi(nn) >= surv_ill_pas(sim%age(ii),1)) then
                   !** Individual dies -------------------------------
                   sim%t2(ii) = dead
                   d_new(ii)  = 1
                   new_in_state = new_in_state + 1
                End if
             else
                if (risk_pi(nn) >= surv_ill_pas(sim%age(ii),2)) then
                   !** Individual dies -------------------------------
                   sim%t2(ii) = dead
                   d_new(ii)  = 1
                   new_in_state = new_in_state + 1
                End if
             End if
          End if
       ENd Do

       !** Population summary --------------------------------------------------
       Call summary_1_int(sim%t2, pop_size, state_count_t2, min_state, max_state)
       
       !** DEBUG --- Population Summary -------------------------------------
       if (PT_DEBUG) then
          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary after dying from infection @ timestep ",timestep
          write(un_lf,'(8(I10))')-1,0,1,2,3,4,5,6
          write(un_lf,'(8(I10))')state_count_t2
       End if
       !** DEBUG ------------------------------------------------------------

       at_risk = state_count_t1(ill_contag)-new_in_state

       if (PT_DEBUG) then
          write(un_lf,PTF_M_AI0)"At risk to move to icu:",at_risk
       End if
       
       new_in_state = 0
       
       !** Draw risks for all individuals at risk --------------------
       Call random_Number(risk_pi(1:at_risk))
       nn = 0
       Do ii = 1, pop_size

          if ((sim%t1(ii) == ill_contag) .AND. (sim%t2(ii) /= dead)) then
             !** Individual is ill_contag ----------------------------
             nn = nn + 1

             if ((ICU_risk_pasd(sim%age(ii),1,sim%d(ii))  >0.) .or. &
                 (ICU_risk_pasd(sim%age(ii),2,sim%d(ii))  >0.)       ) &
                 new_in_state = new_in_state + 1
             
             if (sim%sex(ii) == "m") then
                if (risk_pi(nn) <= ICU_risk_pasd(sim%age(ii),1,sim%d(ii))) then
                   !** Individual moves to icu -----------------------
                   sim%t2(ii) = ill_ICU
                   d_new(ii)  = 1
                Else
                   if (sim%d(ii) >= ill_dur) then
                      !** Individual becomes immune ------------------------
                      sim%t2(ii) = immune
                      d_new(ii)  = 1
                   Else
                      !** Individual stays ill_contagious ------------------
                      d_new(ii) = d_new(ii)  + 1
                   End if
                End if
             else
                if (risk_pi(nn) <= ICU_risk_pasd(sim%age(ii),2,sim%d(ii))) then
                   !** Individual moves to icu -----------------------
                   sim%t2(ii) = ill_ICU
                   d_new(ii)  = 1
                Else
                   if (sim%d(ii) >= ill_dur) then
                      !** Individual becomes immune ------------------------
                      sim%t2(ii) = immune
                      d_new(ii)  = 1
                   Else
                      !** Individual stays ill_contagious ------------------
                      d_new(ii) = d_new(ii)  + 1
                   End if
                End if
             End if
          End if
         
       ENd Do

       if (PT_DEBUG) then
          write(un_lf,PTF_M_AI0)"Really at risk        :",new_in_state
          write(un_lf,PTF_M_AI0)"risks evaluated       :",nn
       End if
       
       !** Population summary --------------------------------------------------
       Call summary_1_int(sim%t2, pop_size, state_count_t2, min_state, max_state)
       
       !** DEBUG --- Population Summary -------------------------------------
       if (PT_DEBUG) then
          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary after moving to ICU @ timestep ",timestep
          write(un_lf,'(8(I10))')-1,0,1,2,3,4,5,6
          write(un_lf,'(8(I10))')state_count_t2
       End if
       !** DEBUG ------------------------------------------------------------

       ! call end_timer("+- From infected to ill and icu")
       
       at_risk = state_count_t1(ill_ICU)

       if (PT_DEBUG) then
          write(un_lf,PTF_M_AI0)"At risk to die in icu:",at_risk
       End if
       
       !** Draw risks for all individuals at risk --------------------
       Call random_Number(risk_pi(1:at_risk))
       nn = 0
       Do ii = 1, pop_size

          if (sim%t1(ii) == ill_ICU) then
             !** Individual was already ill_ICU ----------------------
             nn = nn + 1

             if (sim%sex(ii) == "m") then
                if (risk_pi(nn) >= surv_icu_pas(sim%age(ii),1)) then
                   !** Individual dies -------------------------------------
                   sim%t2(ii) = dead
                   d_new(ii)  = 1
                Else
                   if (sim%d(ii) >= icu_dur) then
                      !** Individual becomes immune ------------------------
                      sim%t2(ii) = immune
                      d_new(ii)  = 1
                   Else
                      !** Individual stays in ICU --------------------------
                      d_new(ii) = d_new(ii)  + 1
                   End if
                End if
             else
                if (risk_pi(nn) >= surv_icu_pas(sim%age(ii),2)) then
                   !** Individual dies -------------------------------
                   sim%t2(ii) = dead
                   d_new(ii)  = 1
                Else
                   if (sim%d(ii) >= icu_dur) then
                      !** Individual becomes immune ------------------------
                      sim%t2(ii) = immune
                      d_new(ii)  = 1
                   Else
                      !** Individual stays in ICU --------------------------
                      d_new(ii) = d_new(ii)  + 1
                   End if
                End if
             End if
          End if

       End Do

       if (PT_DEBUG) then
          write(un_lf,PTF_M_AI0)"risks evaluated      :",nn
       End if
       
       !** Population summary --------------------------------------------------
       Call summary_1_int(sim%t2, pop_size, state_count_t2, min_state, max_state)
       
       !** DEBUG --- Population Summary -------------------------------------
       if (PT_DEBUG) then
          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary after timestep ",timestep
          write(un_lf,'(8(I10))')-1,0,1,2,3,4,5,6
          write(un_lf,'(8(I10))')state_count_t2
       End if
       !write(*,'(9(I10))')timestep,state_count_t2
       !** DEBUG ------------------------------------------------------------

       !** Population summary per location -------------------------------------
       Call summary_2_int(&
            sim%t2, sim%dist_id_rn, pop_size, &
            state_count_pl, min_state, max_state, 1, n_counties)

       healthy_cases(:,timestep,it_ss)    = state_count_pl(healthy   ,:)
       inf_noncon_cases(:,timestep,it_ss) = state_count_pl(inf_noncon,:)
       inf_contag_cases(:,timestep,it_ss) = state_count_pl(inf_contag,:)
       ill_contag_cases(:,timestep,it_ss) = state_count_pl(ill_contag,:)
       ill_ICU_cases(:,timestep,it_ss)    = state_count_pl(ill_ICU   ,:)
       immune_cases(:,timestep,it_ss)     = state_count_pl(immune    ,:)
       dead_cases(:,timestep,it_ss)       = state_count_pl(dead      ,:)
       
       sim%t1 = sim%t2
       sim%d  = d_new

       !!=======================================================================
       !! Write restart data per thread ========================================
       !!=======================================================================
       if ( cp_write .AND. (timestep == cp_time) ) then

          if( (trim(exec_type) .eq. "Optimization") .and. &
            ( OMP_GET_THREAD_NUM() .eq. 0) ) then
            write(cp_unit,  &
                 pos = int(int_storage_size,8)  * int(pop_size,8) &
                        + &
                       int(int1_storage_size,8) * int(pop_size,8) * &
                       int(OMP_GET_THREAD_NUM() ,8) + 1_8)sim%t1
            
            write(cp_unit, &
                 pos = int(int_storage_size,8)  * int(pop_size,8) &
                        + &
                       int(int1_storage_size,8) * int(pop_size,8) &
                        + &
                       int(int_storage_size,8)  * int(pop_size,8) * &
                       int(OMP_GET_THREAD_NUM() ,8) + 1_8)sim%d
            
          else
            !$OMP critical
            !write(un_lf,PTF_M_AI0)"Writing sim%t1 at position:",&
            !     int(int_storage_size,8)  * int(pop_size,8) * &
            !     int(OMP_GET_NUM_THREADS(),8) + &
            !    int(int1_storage_size,8) * int(pop_size,8) * &
            !     int(OMP_GET_THREAD_NUM() ,8) + 1_8
            write(cp_unit,  &
                 pos = int(int_storage_size,8)  * int(pop_size,8) * &
                       int(OMP_GET_NUM_THREADS(),8) + &
                       int(int1_storage_size,8) * int(pop_size,8) * &
                       int(OMP_GET_THREAD_NUM() ,8) + 1_8)sim%t1

            Flush(cp_unit)
            ret = fsync(fnum(cp_unit))
            if (ret /= 0) then
               write(*,*) "ret = fsync(fnum(cp_unit)) = ",ret
               stop
            End if
            
            !write(un_lf,PTF_M_AI0)"Writing sim%d at position:",&
            !     int(int_storage_size,8)  * int(pop_size,8) * &
            !    int(OMP_GET_NUM_THREADS(),8) + &
            !     int(int1_storage_size,8) * int(pop_size,8) * &
            !     int(OMP_GET_NUM_THREADS(),8) + &
            !     int(int_storage_size,8)  * int(pop_size,8) * &
            !     int(OMP_GET_THREAD_NUM() ,8) + 1_8
            
            write(cp_unit, &
                 pos = int(int_storage_size,8)  * int(pop_size,8) * &
                       int(OMP_GET_NUM_THREADS(),8) + &
                       int(int1_storage_size,8) * int(pop_size,8) * &
                       int(OMP_GET_NUM_THREADS(),8) + &
                       int(int_storage_size,8)  * int(pop_size,8) * &
                       int(OMP_GET_THREAD_NUM() ,8) + 1_8)sim%d
            Flush(cp_unit)
            ret = fsync(fnum(cp_unit))
            if (ret /= 0) then
               write(*,*) "ret = fsync(fnum(cp_unit)) = ",ret
               stop
            End if
            
            !$OMP end critical 
          endif
       End if

    ENd Do
    
  End Subroutine CoSMic_TimeLoop

  Subroutine summary_1_int(arr,nn,cnt,lb,ub)

    Integer(kind=1) , Dimension(nn)    , Intent(in)  :: arr
    Integer(kind=ik)                   , Intent(in)  :: nn,lb,ub
    
    Integer(kind=ik), Dimension(lb:ub) , Intent(Out) :: cnt

    Integer(kind=ik)                             :: ii
    
    cnt = 0

    Do ii = 1, nn      
       cnt(arr(ii)) = cnt(arr(ii)) + 1
    End Do

  End Subroutine summary_1_int

  Subroutine summary_2_int(arr1,arr2,nn, cnt,lb1,ub1,lb2,ub2)

    Integer(kind=1)  , Dimension(nn)    , Intent(in)  :: arr1
    Integer(kind=ik) , Dimension(nn)    , Intent(in)  :: arr2
    Integer(kind=ik)                    , Intent(in)  :: nn,lb1,ub1,lb2,ub2
    
    Integer(kind=ik), Dimension(lb1:ub1,lb2:ub2) , Intent(Out) :: cnt

    Integer(kind=ik)                             :: ii

    cnt = 0

    Do ii = 1, nn      
       cnt(arr1(ii),arr2(ii)) = cnt(arr1(ii),arr2(ii)) + 1
    End Do

  End Subroutine summary_2_int

  Subroutine init_ICU_risk(&
       control_age_sex, iol, ill_dur, icu_per_day, &
       ICU_risk_pasd)

    character(len=*), intent(in)                 :: control_age_sex
    Type(iols)      , intent(in)                 :: iol
    Integer         , intent(in)                 :: ill_dur
    Integer(Kind=ik), intent(in), Dimension(:)   :: icu_per_day

    Real(Kind=rk), Allocatable, Dimension(:,:,:), Intent(Out) :: ICU_risk_pasd
    
    Character(len=2), Dimension(:), Allocatable :: ch_age
    Character(Len=5), Dimension(:), Allocatable :: ch_sex

    Integer, Dimension(:), Allocatable          :: temp,temp1
    Integer, Dimension(:), Allocatable          :: target_icu_risk_index

    Real(kind=rk)   , Dimension(:), Allocatable          :: icu_risk

    Character(Len=:), allocatable, Dimension(:)     :: transpr_sex
    Character(Len=:), allocatable, Dimension(:)     :: transpr_age_gr
    real(kind=rk)   , pointer, dimension(:)     :: transpr_icu_risk
    
    Type(icu_risk_lists)                        :: icu_risk_list

    integer                                     :: ii
    !** ------------------------------------------------------------------------
    
    If (control_age_sex == "NONE") Then
       ch_age = (/"to"/)
       ch_sex = (/"total"/)
    End If

    If (control_age_sex == "age") Then
       !     ch_age = generate_seq(0,90,5)
       ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
            "40","45","50","55","60",&
            "65","70","75","80","85","90"/)
       ch_sex = (/"total"/)
    End If

    If (control_age_sex == "sex") Then
       ch_age = (/"to"/)
       ch_sex = (/"m","f"/)
    End If

    If (control_age_sex == "age_sex") Then
       !     ch_age = generate_seq(0,90,5)
       ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
            "40","45","50","55","60",&
            "65","70","75","80","85","90"/)
       ch_sex = (/"m","f"/)
    End If

    !** Get columns from table structure ---------------------------------------
    transpr_icu_risk => get_real_column_pointer(iol%trans_pr, "icu_risk")
    transpr_age_gr   = get_char_column(iol%trans_pr, "age_gr")
    transpr_sex      = get_char_column(iol%trans_pr, "sex")
   
    temp  = get_index(transpr_age_gr,ch_age)
    temp1 = get_index(transpr_sex   ,ch_sex)

    target_icu_risk_index = find_and(temp,temp1)

    If (control_age_sex == "age") Then
       If (.Not.Allocated(icu_risk))Then
          Allocate(icu_risk(2*Size(target_icu_risk_index)))
       Endif
       icu_risk(1:Size(target_icu_risk_index)) = transpr_icu_risk(target_icu_risk_index)
       icu_risk(1:Size(target_icu_risk_index)) = &
            1.0 - (1.0 - icu_risk(1:Size(target_icu_risk_index))) ** (1.0/ Real(ill_dur))
       icu_risk(Size(target_icu_risk_index)+1 : 2* Size(target_icu_risk_index)) = &
            icu_risk(1:Size(target_icu_risk_index))
    End If

    If (control_age_sex == "sex") Then
       If (.Not.Allocated(icu_risk))Then
          Allocate(icu_risk(2*19))
       End If
       icu_risk(1:2) = transpr_icu_risk(target_icu_risk_index)
       icu_risk(1:2) = 1.0 - (1.0 - icu_risk(1:2)) ** (1.0/ ill_dur)
       icu_risk(Size(target_icu_risk_index)+1 : 2* Size(target_icu_risk_index)) = icu_risk(2)
       icu_risk(1:Size(target_icu_risk_index)) = icu_risk(1)
    End If

    If (.Not.Allocated(icu_risk_list%age))Then
       Allocate(icu_risk_list%age(Size(icu_risk)*ill_dur))
       Allocate(icu_risk_list%agei(Size(icu_risk)*ill_dur))
       Allocate(icu_risk_list%sex(Size(icu_risk)*ill_dur))
       Allocate(icu_risk_list%risk(Size(icu_risk)*ill_dur))
       Allocate(icu_risk_list%dur(Size(icu_risk)*ill_dur))
    End If

    Do ii = 1, ill_dur
       icu_risk_list%age((ii-1)*Size(icu_risk)+1: ii*Size(icu_risk)) = (/transpr_age_gr(target_icu_risk_index),&
            transpr_age_gr(target_icu_risk_index)/)
       icu_risk_list%sex((ii-1)*Size(icu_risk)+1:(ii-1)*Size(icu_risk)+19) = 'm'
       icu_risk_list%sex((ii-1)*Size(icu_risk)+20:ii*Size(icu_risk)) = 'f'
       icu_risk_list%risk((ii-1)*Size(icu_risk)+1: ii*Size(icu_risk)) = icu_risk * icu_per_day(ii)
       icu_risk_list%dur((ii-1)*Size(icu_risk)+1: ii*Size(icu_risk)) = ii
    End Do
    Do ii = 1,size(icu_risk_list%agei)
       Read(icu_risk_list%age(ii),*)icu_risk_list%agei(ii)
    End Do

    if (PT_DEBUG) then
       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"ICU risk per age group, sex, and duration in ICU."
       Do ii=1, Size(icu_risk)*ill_dur
          write(un_lf,'(A4,I4,A2,I3, F6.3 )') &
               icu_risk_list%age(ii)    , icu_risk_list%agei(ii), icu_risk_list%sex(ii), &
               icu_risk_list%dur(ii), icu_risk_list%risk(ii)
       End Do
    End if
    
    Allocate(ICU_risk_pasd( &
         minval(icu_risk_list%agei):maxval(icu_risk_list%agei),&
         1:2, &
         minval(icu_risk_list%dur):maxval(icu_risk_list%dur)))

    ICU_risk_pasd = 0._rk
    
    Do ii=1, size(icu_risk_list%sex)
       if (icu_risk_list%sex(ii) == "m") then
          ICU_risk_pasd(icu_risk_list%agei(ii),1,icu_risk_list%dur(ii)) = Real(icu_risk_list%risk(ii),rk)
       Else
          ICU_risk_pasd(icu_risk_list%agei(ii),2,icu_risk_list%dur(ii)) = Real(icu_risk_list%risk(ii),rk)
       End if
    End Do
    
  End Subroutine init_ICU_risk

  Subroutine init_surv_ill(&
       control_age_sex, iol, ill_dur, &
       surv_ill_pas)

    character(len=*), intent(in)                 :: control_age_sex
    Type(iols)      , intent(in)                 :: iol
    Integer         , intent(in)                 :: ill_dur

    Real(Kind=rk), Allocatable, Dimension(:,:), Intent(Out) :: surv_ill_pas
    
    Character(len=2), Dimension(:), Allocatable :: ch_age
    Character(Len=5), Dimension(:), Allocatable :: ch_sex

    Integer, Dimension(:), Allocatable          :: temp,temp1
    Integer, Dimension(:), Allocatable          :: target_icu_risk_index

    Real(kind=rk)   , Dimension(:), Allocatable          :: surv_ill
    Character(Len=:), allocatable, Dimension(:) :: transpr_sex
    Character(Len=:), allocatable, Dimension(:) :: transpr_age_gr
    real(kind=rk), pointer, dimension(:)        :: transpr_surv_ill
    
    Type(icu_risk_lists)                        :: surv_ill_list

    integer                                     :: ii
    !** ------------------------------------------------------------------------

    
    If (control_age_sex == "NONE") Then
       ch_age = (/"to"/)
       ch_sex = (/"total"/)
    End If
    
    If (control_age_sex == "age") Then
       !     ch_age = generate_seq(0,90,5)
       ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
            "40","45","50","55","60",&
            "65","70","75","80","85","90"/)
       ch_sex = (/"total"/)
    End If

    If (control_age_sex == "sex") Then
       ch_age = (/"to"/)
       ch_sex = (/"m","f"/)
    End If

    If (control_age_sex == "age_sex") Then
       !     ch_age = generate_seq(0,90,5)
       ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
            "40","45","50","55","60",&
            "65","70","75","80","85","90"/)
       ch_sex = (/"m","f"/)
    End If

    !** Get columns from table structure ---------------------------------------
    transpr_surv_ill => get_real_column_pointer(iol%trans_pr, "surv_ill")
    transpr_age_gr   = get_char_column(iol%trans_pr, "age_gr")
    transpr_sex      = get_char_column(iol%trans_pr, "sex")
   
    temp = get_index(transpr_age_gr,ch_age)

    temp1= get_index(transpr_sex,ch_sex)

    target_icu_risk_index = find_and(temp,temp1)

    ! init surv_ill
    If (control_age_sex == "age") Then
       If (.Not.Allocated(surv_ill))Then
          Allocate(surv_ill(2*Size(target_icu_risk_index)))
       End If
       surv_ill = (/transpr_surv_ill(target_icu_risk_index),&
            transpr_surv_ill(target_icu_risk_index)/)
       surv_ill = surv_ill ** (1.0/Real(ill_dur))
    End If

    If (control_age_sex == "sex") Then
       If (.Not.Allocated(surv_ill))Then
          Allocate(surv_ill(2*19))
       Endif
       surv_ill(1:2) = transpr_surv_ill(target_icu_risk_index)
       surv_ill(1:2) = 1.0 - (1.0 - surv_ill(1:2)) ** (1.0/ Real(ill_dur))
       surv_ill(Size(target_icu_risk_index)+1 : 2* Size(target_icu_risk_index)) = surv_ill(2)
       surv_ill(1:Size(target_icu_risk_index)) = surv_ill(1)
    End If
    
    If (.Not.Allocated(surv_ill_list%age))Then
       Allocate(surv_ill_list%age(Size(surv_ill)))
       Allocate(surv_ill_list%agei(Size(surv_ill)))
       Allocate(surv_ill_list%sex(Size(surv_ill)))
       Allocate(surv_ill_list%risk(Size(surv_ill)))
    End If
    surv_ill_list%age = (/transpr_age_gr(target_icu_risk_index),transpr_age_gr(target_icu_risk_index)/)
    surv_ill_list%sex = 'f'
    surv_ill_list%sex(1:19) = 'm'
    surv_ill_list%risk  = surv_ill
    
    Do ii = 1, Size(surv_ill)
       Read(surv_ill_list%age(ii),*)surv_ill_list%agei(ii)
    End Do

    if (PT_DEBUG) then
       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"Chance of survival per age group and sex."
       Do ii=1, 2*19
          write(un_lf,'(A4,I4, A2,F6.3)') &
               surv_ill_list%age(ii)    , surv_ill_list%agei(ii), surv_ill_list%sex(ii), &
               surv_ill_list%risk(ii)
       End Do
    End if
    
    allocate(surv_ill_pas(minval(surv_ill_list%agei):maxval(surv_ill_list%agei),1:2))
    Do ii=1, 2*19
       if (surv_ill_list%sex(ii) == "m") then
          surv_ill_pas(surv_ill_list%agei(ii),1) = surv_ill_list%risk(ii)
       Else
          surv_ill_pas(surv_ill_list%agei(ii),2) = surv_ill_list%risk(ii)
       End if
    End Do

  End Subroutine init_surv_ill

  subroutine init_surv_icu(&
       control_age_sex, iol, ill_dur, &
       surv_icu_pas)
    
    character(len=*), intent(in)                 :: control_age_sex
    Type(iols)      , intent(in)                 :: iol
    Integer         , intent(in)                 :: ill_dur

    Real(Kind=rk), Allocatable, Dimension(:,:), Intent(Out) :: surv_icu_pas
    
    Character(len=2), Dimension(:), Allocatable :: ch_age
    Character(Len=5), Dimension(:), Allocatable :: ch_sex

    Integer, Dimension(:), Allocatable          :: temp,temp1
    Integer, Dimension(:), Allocatable          :: target_icu_risk_index

    Real(kind=rk), Dimension(:), Allocatable    :: surv_icu
    Character(Len=:), allocatable, Dimension(:) :: transpr_sex
    Character(Len=:), allocatable, Dimension(:) :: transpr_age_gr
    real(kind=rk), pointer, dimension(:)        :: transpr_surv_icu
    Type(icu_risk_lists)                        :: surv_icu_list

    integer                                     :: ii
    !** ------------------------------------------------------------------------
    
    If (control_age_sex == "NONE") Then
       ch_age = (/"to"/)
       ch_sex = (/"total"/)
    End If
    
    If (control_age_sex == "age") Then
       !     ch_age = generate_seq(0,90,5)
       ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
            "40","45","50","55","60",&
            "65","70","75","80","85","90"/)
       ch_sex = (/"total"/)
    End If

    If (control_age_sex == "sex") Then
       ch_age = (/"to"/)
       ch_sex = (/"m","f"/)
    End If

    If (control_age_sex == "age_sex") Then
       !     ch_age = generate_seq(0,90,5)
       ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
            "40","45","50","55","60",&
            "65","70","75","80","85","90"/)
       ch_sex = (/"m","f"/)
    End If
    
    !** Get columns from table structure ---------------------------------------
    transpr_surv_icu => get_real_column_pointer(iol%trans_pr, "surv_icu")
    transpr_age_gr   =  get_char_column(iol%trans_pr, "age_gr")
    transpr_sex      =  get_char_column(iol%trans_pr, "sex")
    
    temp = get_index(transpr_age_gr,ch_age)

    temp1= get_index(transpr_sex,ch_sex)

    target_icu_risk_index = find_and(temp,temp1)

    ! init surv_icu
    If (control_age_sex == "age") Then
       If (.Not.Allocated(surv_icu))Then
          Allocate(surv_icu(2*Size(target_icu_risk_index)))
       End If
       surv_icu= (/transpr_surv_icu(target_icu_risk_index),&
            transpr_surv_icu(target_icu_risk_index)/)
       surv_icu = surv_icu ** (1/Real(ill_dur))
    End If

    If (control_age_sex == "sex") Then
       If (.Not.Allocated(surv_icu))Then
          Allocate(surv_icu(2*19))
       End If
       surv_icu(1:2) = transpr_surv_icu(target_icu_risk_index)
       surv_icu(1:2) = 1.0 - (1.0 - surv_icu(1:2)) ** (1.0/ Real(ill_dur))
       surv_icu(Size(target_icu_risk_index)+1 : 2* Size(target_icu_risk_index)) = surv_icu(2)
       surv_icu(1:Size(target_icu_risk_index)) = surv_icu(1)
    End If
    If (.Not.Allocated(surv_icu_list%age))Then
       Allocate(surv_icu_list%age(Size(surv_icu)))
       Allocate(surv_icu_list%agei(Size(surv_icu)))
       Allocate(surv_icu_list%sex(Size(surv_icu)))
       Allocate(surv_icu_list%risk(Size(surv_icu)))
    End If
    surv_icu_list%age = (/transpr_age_gr(target_icu_risk_index),transpr_age_gr(target_icu_risk_index)/)
    surv_icu_list%sex = 'f'
    surv_icu_list%sex(1:19) = 'm'
    surv_icu_list%risk  = surv_icu
    Do ii = 1, size(surv_icu)
       Read(surv_icu_list%age(ii),*)surv_icu_list%agei(ii)
    End Do

    if (PT_DEBUG) then
       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"Chance of survival in ICU per age group and sex."
       Do ii=1, 2*19
          write(un_lf,'(A4,I4,A2,F6.3)') &
               surv_icu_list%age(ii),  surv_icu_list%agei(ii)  , surv_icu_list%sex(ii), &
               surv_icu_list%risk(ii)
       End Do
    End if
    
    allocate(surv_icu_pas(minval(surv_icu_list%agei):maxval(surv_icu_list%agei),1:2))
    Do ii=1, 2*19
       if (surv_icu_list%sex(ii) == "m") then
          surv_icu_pas(surv_icu_list%agei(ii),1) = surv_icu_list%risk(ii)
       Else
          surv_icu_pas(surv_icu_list%agei(ii),2) = surv_icu_list%risk(ii)
       End if
    End Do

  end subroutine init_surv_icu
 
End Module kernel
