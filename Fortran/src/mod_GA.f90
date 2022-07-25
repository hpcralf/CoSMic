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
!# Authors:      Huan Zhou
!#               Ralf Schneider
!#               Christian Dudel
!#               Matthias Rosenbaum-Feldbruegge
!#               Sebastian Kluesener
!#
!# Contact:      ralf.schneider@hlrs.de
!#               huan.zhou@hlrs.de
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
! Module containing the GA and its dependent subroutines
!
!###############################################################################
Module genetic_algorithm

  use param_tree
  use precision
  use urandom
  use global_types
  use precision
  use cosmic_io
  use kernel
  use qsort_c_module
  use quicksort_nr 

  Use cosmic_io

  use mpi

  implicit none

  !** Constants ---------------------------------------
  Real, parameter     :: INV     = 99.0 ! Corresponds to NA for fitness array
  Real, parameter     :: EPS     = EPSILON(1.0) ! Tolerance precision, a magnitude of e-07
  Real, parameter     :: NULVAL  = -1.0 ! Corresponds to NA for icu_nuts2_cases
  Integer,parameter   :: ROOT    = 0 ! The root rank
  Integer,parameter   :: WEEKLEN = 7 ! The length of a week by days

  !** Extend the Five-number summary statistic with a mean value -------
  !** The description of the Five-number summary can be found: ---------
  !** https://en.wikipedia.org/wiki/Five-number_summary ----------------
  type ga_summary
    Real,Allocatable,dimension(:)                 :: maximum
    Real,Allocatable,dimension(:)                 :: mean
    Real,Allocatable,dimension(:)                 :: upper_hinge
    Real,Allocatable,dimension(:)                 :: median
    Real,Allocatable,dimension(:)                 :: lower_hinge
    Real,Allocatable,dimension(:)                 :: minimum
  end type ga_summary

  !** GA results/best solution ------------------------------------------
  type opt_res 
    Real(kind=rk),dimension(:), Allocatable       :: opt_bestsol ! Solution
    Real                                          :: opt_fitnessvalue ! Best fitness
    Integer                                       :: opt_acuiter ! The occured number of generations
  end type opt_res 

Contains

  subroutine ga(&
    iol, &
    iter_s, iter_e, &
    inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
    less_contagious, R0_force, immune_stop, &
    R0change, R0delay ,R0delay_days, R0delay_type, &
    control_age_sex, seed_date, seed_before, sam_size, R0, &
    R0_effects, full_region_index, region_index, output_dir, export_name, rank_mpi, size_mpi)

    !===========================================================================
    ! Declaration
    !===========================================================================

    Type(iols), Target                           , intent(in) :: iol
    Integer(kind=ik)                             , intent(in) :: iter_s, iter_e
    Integer(kind=ik)                             , intent(in) :: inf_dur
    Integer(kind=ik)                             , intent(in) :: cont_dur
    Integer(kind=ik)                             , intent(in) :: ill_dur
    Integer(kind=ik)                             , intent(in) :: icu_dur
    Integer(kind=ik), Allocatable, Dimension(:)  , intent(in) :: icu_per_day
    Real(kind=rk)                                , intent(in) :: less_contagious
    Real(kind=rk)                                , intent(in) :: R0_force
    Logical                                      , intent(in) :: immune_stop
    Integer(kind=ik), Allocatable, Dimension(:,:), intent(in) :: R0change
    Logical                                      , intent(in) :: R0delay
    Integer(kind=ik)                             , intent(in) :: R0delay_days
    Character(len=*)                             , intent(in) :: R0delay_type
    character(len=*)                             , intent(in) :: control_age_sex
    character(len=*)                             , intent(in) :: seed_date
    Integer(kind=ik)                             , intent(in) :: seed_before    
    Integer(kind=ik)                             , intent(in) :: sam_size
    Real(kind=rk)                                , intent(in) :: R0
    Real(kind=rk), Allocatable, Dimension(:,:), intent(inout) :: R0_effects ! As output
    Integer(kind=mpi_ik)                         , intent(in) :: rank_mpi
    Integer(kind=mpi_ik)                         , intent(in) :: size_mpi
    integer(kind=ik), Allocatable, Dimension(:)  , intent(in) :: full_region_index ! A list of all regions (nuts2s or states)
    integer(kind=ik), Allocatable, Dimension(:)  , intent(in) :: region_index ! A list of subregions
    character(len=:), Allocatable                , intent(in) :: output_dir
    character(len=:), Allocatable                , intent(in) :: export_name

    !==========================================================================
    ! General purpose for preprocessing and fitness calculation
    !==========================================================================     
    Real(kind=rk), Allocatable, Dimension(:,:)    :: local_ini_pop ! Store the local population for each rank
    Real(kind=rk), Allocatable, Dimension(:,:)    :: ini_pop, tmp_pop ! ini_pop: store the global population, which are the input to be optimized.
    Real(kind=rk), Allocatable, Dimension(:,:)    :: sorted_pop

    Real, Allocatable, Dimension(:)               :: local_fitness ! Store the local fitness for each rank
    Real, Allocatable, Dimension(:)               :: fitness, tmp_fitness ! fitness: store the calculated fitness 
    Real, Allocatable, Dimension(:)               :: sorted_fitness
    Type(ga_summary)                              :: summary_fitness

    Integer,dimension(:,:,:), Allocatable         :: ill_ICU_cases_final ! Store the evaluated ICU cases
    Integer(kind=ik)                              :: i, j, k, m, pop_ii, gene_ii
    Integer                                       :: startid, endid, max_pos, min_pos
    Real                                          :: tmp_max, run_sum
    Integer                                       :: time_n
    
    Real                                          :: diff_sum, diff_loc_sum, diff_tmp_power
    integer                                       :: diff_len, diff_loc_len
    integer                                       :: half
    Integer,dimension(:,:,:), Allocatable         :: aggreg_sum_by_county
    Real,dimension(:,:), Allocatable              :: avg_by_county ! Record the simulated ICU cases by state/nuts2

    Character(Len=10), allocatable, Dimension(:)  :: obsicu_nuts2_date
    Integer(kind=ik), Allocatable, Dimension(:)   :: obsicu_nuts2_cases
    Integer(kind=ik), Allocatable, Dimension(:)   :: obsicu_counties

    Character(Len=10), allocatable, Dimension(:)  :: obsicu_date
    Integer(kind=ik), allocatable, Dimension(:)   :: obsicu_cases
    Character(Len=15), allocatable, Dimension(:)  :: obsicu_sc   

    Integer                                       :: elapsed_days, num_actudays
    Integer(kind=ik)                              :: iter
    Integer                                       :: R0_effects_dim2size, R0change_dim2size
    Integer                                       :: end_week
    Integer                                       :: tail_pos


    !=================================================
    ! Optimization target "state"
    !=================================================
    Integer                                       :: num_states
    Integer, Allocatable                          :: pos_state(:),uniq_distid(:)
    Character(len=15), Allocatable, Dimension(:)  :: states_shortcut

    ! Record the observed ICU cases by state
    Real,Allocatable,dimension(:,:)               :: seltarget_icu_by_state_icucases
    Character(len=15),dimension(:,:), Allocatable :: seltarget_icu_by_state_icudates

    !=================================================
    ! Optimization target "nuts2":
    !=================================================
    Integer,allocatable,dimension(:)                   :: unique_nuts2s
    Integer                                            :: num_nuts2
    Integer(kind=ik)                                   :: loc_sum
    Integer,Allocatable, dimension(:)                  :: pos_nuts2
    Integer(kind=ik), Allocatable, Dimension(:)        :: counties_dist_id
    ! Record the observed ICU cases by nuts2
    Real,Allocatable,dimension(:,:)                    :: seltarget_icu_by_nuts2_icucases
    Character(len=15),dimension(:,:), Allocatable      :: seltarget_icu_by_nuts2_icudates

    !====================================================================================================
    ! Convulution filter and intersection of the observed and evaluation data according to date
    !====================================================================================================
    Real, Allocatable,dimension(:)                :: tmp_filter
    Character(len=10)                             :: refer_date,refer_date_mod
    Integer                                       :: halve_coff_len, coff_len ! 7
    Integer, Allocatable                          :: obs_lb(:), obs_ub(:), eval_lb(:), eval_ub(:)
    Integer                                       :: inteval_days, status
    Real                                          :: coefficient

    !=================================================
    ! GA - selection->crossover->mutation->elitism
    !=================================================
    Integer, Allocatable, dimension(:)            :: sel_par ! The selected parents
    Type(opt_res)                                 :: ga_res
    Integer, Allocatable, dimension(:)            :: opt_index, opt_week_index
    Integer                                       :: opt_local_size ! The number of local population when there are more than two MPI ranks
    Integer                                       :: nvals ! The total number of parameters to be optimized
     integer, Allocatable, dimension(:)           :: sorted_index
    integer                                       :: idx
    real                                          :: pre_fitv, curr_fitv

    !=================================================
    ! Automation
    !=================================================
    Type(sims)                                    :: sim, sim_backup ! sim_backup exists forau guaranteeing data correctness
    Type(opt_sim_switchs)                         :: sim_switch
    Integer                                       :: start_window, end_window ! [start_window:end_window] spans the calibrated weeks
    Real(kind=rk), Allocatable, Dimension(:,:)    :: tmp_R0_effects ! Temporarily store the R0 values, currently being calibrated
    Integer                                       :: obs_starttime, eval_starttime, diff_days, diff_weeks
    Integer(kind=ik), Allocatable, Dimension(:,:) :: cor_R0change, sim_R0change ! Other than the given R0change

    !=================================================================
    ! Optimization Control Panel / Setting the optimization/automation parameters
    !=================================================================
    Logical                                       :: opt_target_icu !.TRUE. Optimiziation target
    Logical                                       :: opt_target_deaths ! .FALSE. Optimization target
    Logical                                       :: opt_filter ! .TRUE.
    Logical                                       :: sim_reload ! .FALSE. If we need to reload the existing
                                                                ! sim data, indicated by start_week (see below)
    Logical                                       :: auto_calibrate ! .TRUE. Multiple (automation) or one-off calibration
    character(len=mcl)                            :: use_sug_sol ! "NULL"
    character(Len=:),Allocatable                  :: opt_target_region ! "state". Optimized region. Valid values: "state" and "nuts2"
    Integer(kind=ik)                              :: opt_pop_size ! 4 Population size
    Integer(kind=ik)                              :: opt_num_gene ! 2 Maximal number of generations
    Integer(kind=ik)                              :: opt_window_size, opt_shift_size ! opt_window_size: the number of calibrated weeks
                                                                                     ! opt_shift_size : move the window backwards by opt_shift_size
    Integer(kind=ik)                              :: current_week, start_week ! Calibrate until the current_week reaches
    Integer                                       :: opt_elitism ! The number of the chosen elitism
    Real                                          :: opt_pcrossover ! 0.8
    Real                                          :: opt_pmutation ! 0.1

    Character(len=:),Allocatable,dimension(:)     :: opt_names ! The optimized targets
    Real                                          :: opt_lb, opt_ub
    !=======================================================
    ! MPI-related
    !=======================================================
    Integer(kind=mpi_ik)                          :: ierr

    !========================================================================================================== 
    call pt_get("#opt_target_icu",opt_target_icu)
    call pt_get("#opt_target_deaths",opt_target_deaths)
    call pt_get("#opt_filter",opt_filter)
    call pt_get("#opt_pop_size",opt_pop_size)
    call pt_get("#opt_max_iter",opt_num_gene)
    call pt_get("#opt_target_region",opt_target_region)
    call pt_get("#auto_calibrate",auto_calibrate)
    call pt_get("#sim_reload",sim_reload)
    
    opt_local_size = opt_pop_size/size_mpi
    R0change_dim2size   = size(R0change, dim=2)
    if ((opt_local_size .eq. 0) .or. (mod(opt_pop_size,size_mpi) .ne. 0)) then
        print *,"Error: too many ranks or the population is not divisable by rank size!"
        Call exit(status)
    endif
    if (auto_calibrate) then ! Calibrate until the current week reaches
      call pt_get("#opt_window_size",opt_window_size)
      call pt_get("#opt_shift_size",opt_shift_size)
      call pt_get("#current_week",current_week)


      !** Parameter checking -----------------------------------------------------
      if (R0change_dim2size .ne. opt_window_size) then
        print *,"Error: the size of the given R0change is not equal to 'opt_window_size'!"
        Call exit(status)
      endif
      if ((opt_shift_size .le. 0) .or. (opt_window_size .lt. opt_shift_size)) then 
        print *,"Error: the shift size should be a positive value and &
                  the window size should not be smaller than the shift size!"
        Call exit(status)
      endif
    else ! Otherwise, the one-off calibration targets is given by opt_names
      call pt_get("#opt_names",opt_names)
    endif

    opt_pcrossover    = 0.8
    opt_pmutation     = 0.1
    
    !** TODO: shall we obtain the value of opt_lb and opt_ub dynamically (from static_parameters.dat) --
    opt_lb = 0.1
    opt_ub = 1.0
    opt_elitism = max(1,NINT(opt_pop_size*0.05))
    if (rank_mpi .eq. ROOT) then
      call print_ga_header(opt_pop_size, opt_num_gene, opt_elitism, opt_pcrossover, opt_pmutation, &
                    opt_target_region)
    endif
  
    iter = iter_e - iter_s + 1
    coff_len = 7
    R0_effects_dim2size = size(iol%R0_effect%head)

    !** Initializing sim_switch -------
    sim_switch%sim_reuse  = .false.
    sim_switch%sim_write  = .false.
    sim_switch%sim_reload = sim_reload

    select case (trim(opt_target_region))
      !** When the optimized target region is "Nuts2" ------------------
      case("nuts2") 
        counties_dist_id = get_int_table_column(iol%counties,'dist_id')
      
        !** unique_nuts2s: set of affected nuts2; pos_nuts2: see the function "get_unique_nuts2_state" ----
        unique_nuts2s = get_unique_nuts2_state(region_index,pos_nuts2)
        num_nuts2 = size(unique_nuts2s) ! Get the number of affected nuts2

        !** Get the position of the given nuts2 in the full nuts2 list -----
        !** The results are stored back to array "unique_nuts2s" -----------
        do i = 1, num_nuts2
          do j = 1, size(full_region_index)
            if (unique_nuts2s(i) .eq. full_region_index(j)) then
              unique_nuts2s(i) = j
              exit
            endif
          enddo
        enddo
      
        Allocate(obs_lb(num_nuts2))
        Allocate(obs_ub(num_nuts2))
        Allocate(eval_lb(num_nuts2))
        Allocate(eval_ub(num_nuts2))

        obs_lb  = 0
        obs_ub  = 0
        eval_lb = 0
        eval_ub = 0

        obsicu_nuts2_date     = get_char_column(iol%obsicu_nuts2,'date')
        obsicu_nuts2_cases    = get_int_table_column(iol%obsicu_nuts2,'cases')
        obsicu_counties       = get_int_table_column(iol%obsicu_nuts2,'dist_id')
        
        !** TODO: The fitness relavant to the death data ------
        !** The fitness relavant to the icu data --------------
        if (opt_target_icu) then
          !** Precalculate the size of the observed data -----------------------------------------------
          elapsed_days = (Date2Unixtime(trim(obsicu_nuts2_date(size(obsicu_nuts2_date))))& 
                      - (Date2Unixtime(trim(obsicu_nuts2_date(1)))))/86400

          num_actudays = elapsed_days+1
          Allocate(seltarget_icu_by_nuts2_icucases(num_actudays,num_nuts2))
          seltarget_icu_by_nuts2_icucases = NULVAL ! Initialization
          Allocate(seltarget_icu_by_nuts2_icudates(num_actudays,num_nuts2))
          
          do i = 1, num_nuts2
            j = 1
            loc_sum = 0
            seltarget_icu_by_nuts2_icudates(j,i) = trim(obsicu_nuts2_date(j))
            !** Calculate the icu cases of certain nuts2 --------------
            do k = 1, size(obsicu_nuts2_date)
              !** When we reach a different date, ------------------------------------------------------
              !** the loc_sum for the previous date will be assigned to the current nuts2_icucases -----
              if (trim(obsicu_nuts2_date(k)) .ne. trim(seltarget_icu_by_nuts2_icudates(j,i))) then
                seltarget_icu_by_nuts2_icucases(j,i) = Real(loc_sum)
                elapsed_days = (Date2Unixtime(trim(obsicu_nuts2_date(k))) - &
                      Date2Unixtime(trim(seltarget_icu_by_nuts2_icudates(j,i))))/86400

                ! When the two observed dates are more than one day apart
                do m = 1,elapsed_days-1
                  j = j + 1
                  seltarget_icu_by_nuts2_icudates(j,i) = add_date(trim(seltarget_icu_by_nuts2_icudates(j-1,i)),1)
                  seltarget_icu_by_nuts2_icucases(j,i) = NULVAL ! Set to NULL for the unknown date (not observed)
                enddo
                j = j + 1
                seltarget_icu_by_nuts2_icudates(j,i) = trim(obsicu_nuts2_date(k))                    
                loc_sum = 0
                 
              endif
              !** Accumulate only when the given nuts2 is hit ---------------------------
              if (BinarySearch(obsicu_counties(k),&
                    counties_dist_id(unique_nuts2s(i):(unique_nuts2s(i)+pos_nuts2(i+1)-pos_nuts2(i)-1))) .ne. 0) then
                loc_sum = loc_sum + obsicu_nuts2_cases(k)
              endif        
            enddo
            seltarget_icu_by_nuts2_icucases(j,i) = Real(loc_sum)
          
            !** Apply convolution filter ---------------
            if (opt_filter) then
              call conv_filter(seltarget_icu_by_nuts2_icucases(:,i), coff_len, num_actudays)
            endif
            obs_starttime = Date2Unixtime(seltarget_icu_by_nuts2_icudates(1,1))
          enddo
        endif
           
      !** When the optimized target region is "state" ------------------- 
      case("state")

        !!TODO: to unify the pos_state and pos_nuts2, they can have the same names
        ! uniq_distid: the affected states; pos_state: collect the beginning index from which a different state starts
        uniq_distid = get_unique_nuts2_state(region_index, pos_state) 
        num_states = size(uniq_distid) ! Get the number of the affected states

        Allocate(obs_lb(num_states))
        Allocate(obs_ub(num_states))
        Allocate(eval_lb(num_states))
        Allocate(eval_ub(num_states))
        obs_lb  = 0
        obs_ub  = 0
        eval_lb = 0
        eval_ub = 0

        obsicu_date  =  get_char_column(iol%obsicu_state,'date')
        obsicu_cases =  get_int_table_column(iol%obsicu_state,'cases')
        obsicu_sc    =  get_char_column(iol%obsicu_state,'state')

        states_shortcut = get_char_column(iol%states,'Shortcut')

        !** TODO: The fitness relavant to the death data ------
        !** The fitness relavant to the icu data --------------
        if (opt_target_icu) then
          !** Precalculate the size of the observed data --------------------------
          elapsed_days = (Date2Unixtime(trim(obsicu_date(size(obsicu_date))))& 
                        - Date2Unixtime(trim(obsicu_date(1))))/86400
          num_actudays = elapsed_days + 1
          Allocate(seltarget_icu_by_state_icucases(num_actudays,num_states))
          Allocate(seltarget_icu_by_state_icudates(num_actudays,num_states))
          seltarget_icu_by_state_icucases = NULVAL
          !** Presume seltarget_icu_by_state is arranged in an ascending order ----------
          !** Preprocess the observed ICU data according to the affected states ---------
          do i = 1,num_states
            k = 1
            do j = 1, size(obsicu_sc)
              if (trim(states_shortcut(uniq_distid(i))) .eq. trim(obsicu_sc(j))) then
                ! When the first matched state is hit
                if (k .eq. 1) then
                  seltarget_icu_by_state_icudates(k,i) = trim(obsicu_date(j))
                  seltarget_icu_by_state_icucases(k,i) = Real(obsicu_cases(j))
                  k = k + 1
                  cycle
                endif

                elapsed_days = (Date2Unixtime(trim(obsicu_date(j))) - &
                      Date2Unixtime(trim(seltarget_icu_by_state_icudates(k-1,i))))/86400
                ! When the two observed dates are more than one day apart
                do m = 1,elapsed_days-1
                  seltarget_icu_by_state_icudates(k,i) = add_date(trim(seltarget_icu_by_state_icudates(k-1,i)),1)
                  seltarget_icu_by_state_icucases(k,i) = NULVAL ! Set to NULL for the unknown date (not observed)
                  k = k + 1
                enddo
                seltarget_icu_by_state_icudates(k,i) = trim(obsicu_date(j))
                seltarget_icu_by_state_icucases(k,i) = Real(obsicu_cases(j))
                k = k + 1
              endif
            enddo
             
            !** Apply convolution filter ------
            if (opt_filter) then
              call conv_filter(seltarget_icu_by_state_icucases(:,i), coff_len, m)                      
            endif
            obs_starttime = Date2Unixtime(seltarget_icu_by_state_icudates(1,1))
          enddo
        end if
      case default
        write(*,*)  "Otherwise, the given target is not supported yet"
        call exit(status) 
    end select

    !** The checkpoint reloading and wirte are performed in weeks
    if (.not. sim_reload) then
      start_week = 1 ! Start from scratch
    else 
      call pt_get("#start_week",start_week)
      sim_switch%read_week = start_week ! Reload the checkpoint of start_week
    endif
    refer_date = add_date(seed_date,(start_week-1)*WEEKLEN)
    refer_date_mod  = add_date(refer_date,1) ! Used for the calculation of the below data intersection

    if (auto_calibrate) then
      !** Note: the weeks in R0_effects (as output) should start from the seed(first) week up to the current week ----
      deallocate(R0_effects)
      allocate(R0_effects(current_week,R0_effects_dim2size))
      start_window = start_week
      if (current_week .lt. start_window) then
        print *,"Error: the current week should not be smaller than the given start week!"
        Call exit(status)
      endif
      
      eval_starttime = Date2Unixtime(refer_date_mod) 
      diff_days = obs_starttime - eval_starttime
      diff_days = diff_days/86400
      
      !** (diff_days/WEEKLEN) indicates the smallest window size meeting the intersection requirement ----
      if ((diff_days/WEEKLEN) .gt. (opt_window_size - 1)) then
        end_window = diff_days/WEEKLEN + start_window
        if (current_week .lt. end_window) then
          print *,"Error: the current week should be large enough to intersect with the observed data!"
          Call exit(status)
        endif
      else
        end_window = start_window + opt_window_size - 1
        if (current_week .lt. end_window) then
          end_window = current_week
        endif
      endif

      !** diff_weeks indicates the initial window size (the number of calibrated weeks) -----
      diff_weeks = end_window - start_window + 1
      if (diff_weeks .ne. opt_window_size) then
        !** When the initial window size is unequal to the given window size, then set 
        !** the specific cor_R0change, otherwise R0change applies
        Allocate(cor_R0change(size(R0change,dim=1),diff_weeks))
        call ini_R0change(cor_R0change)
        time_n = maxval(cor_R0change)+1
      else
        time_n = maxval(R0change) + 1
      endif

      if ( diff_weeks-(opt_window_size-opt_shift_size) .gt. 0 ) then 
        Allocate(sim_R0change(size(R0change,dim=1),diff_weeks-(opt_window_size-opt_shift_size)))
      else
        Allocate(sim_R0change(size(R0change,dim=1),diff_weeks-1))
      endif
      if (size(sim_R0change,dim=2) .ne. 0) then
        call ini_R0change(sim_R0change)
      endif

      !!** By default the optimization target cover all the nuts2 regions or states, the order is indentical to
      !!** the one shown in the R0effect header. ---
      allocate(tmp_R0_effects(diff_weeks,R0_effects_dim2size))

      nvals = diff_weeks*R0_effects_dim2size
      Allocate(opt_week_index(nvals))
      do i = start_window-start_week+1, end_window-start_week+1
        opt_week_index(i) = start_window+i-1
      enddo
      do i = end_window-start_week+2,nvals
        opt_week_index(i) = opt_week_index(i-diff_weeks)
      enddo
      Allocate(opt_index(nvals))
      do i = 1,nvals
        opt_index(i) = CEILING(i/Real(diff_weeks))
      enddo
    else
      !** The optimization target is given by 'opt_names', if the auto_calibration is disabled ----
      time_n = maxval(R0change) + 1
      !** Preparation for preprocessing the parameters subject to optimization -----------------------------------------------
      !** Collect the indexes where the R0 effect data should be replaced with the input that is subject to optimization -----
      nvals = size(opt_names)
      Allocate(opt_index(nvals))
      opt_index = 0
      do j=1,nvals
        do i=1,size(iol%R0_effect%head)
          if (index(trim(opt_names(j)), trim(iol%R0_effect%head(i))) .ne. 0) then
              opt_index(j) = i
              EXIT
          endif
        enddo
        if (opt_index(j) .eq. 0) then
          print *,"Error: the optimized target region ",trim(opt_names(j))," is not expected!"
          Call exit(status)
        endif
      enddo

      Allocate(opt_week_index(nvals))
      opt_week_index = 0
      tail_pos = len(trim(opt_names(1)))
      do i = 1,nvals
        READ(opt_names(i)(tail_pos:tail_pos), "(I4)") j
        opt_week_index(i) = j
      enddo
      start_window = opt_week_index(1)
      end_window   = opt_week_index(nvals)
      current_week = end_window
      opt_window_size = end_window - start_window + 1
      end_week = start_week + R0change_dim2size - 1

      !** Parameter checking -------------------------------------------------------
      if (start_week .gt. start_window) then
        print *,"Error: the optimized target is ahead of the starting date!"
        Call exit(status)
      endif
      if (end_week .lt. end_window) then
        print *,"Error: R0change is not large enough to cover all the optimized target!"
        Call exit(status)
      endif
      if (end_week .gt. R0_effects_dim2size) then
        print *,"Error: R0change is so large that execeed the upper bound of R0_effects!"
        Call exit(status)
      endif

      allocate(tmp_R0_effects(R0change_dim2size,R0_effects_dim2size))
    endif

    Allocate(local_fitness(opt_local_size))
    !** Each rank initializes the local fitness
    do i = 1,opt_local_size
      local_fitness(i) = INV
    enddo
  
    !** Population: generate beginning population by using a uniform range distribution -----
    Allocate(local_ini_pop(nvals, opt_local_size))

    do j = 1, opt_local_size
      do i = 1, nvals
        !! The random number is between opt_lb and opt_ub.
        local_ini_pop(i,j) = random_uniform(opt_lb,opt_ub)
      enddo
    enddo

    Allocate(ga_res%opt_bestsol(nvals))
    Allocate(ini_pop(nvals, opt_pop_size))
   
    do while(start_window .le. current_week)
      select case (trim(opt_target_region))
        case("nuts2")
        !** Calculate the intersection of the evaluated and observed ICU data by nuts2 and date --
          do i = 1, num_nuts2
            call date_intersection(obs_lb(i), eval_lb(i), obs_ub(i), eval_ub(i), refer_date_mod, &
                                  time_n, num_actudays, seltarget_icu_by_nuts2_icudates(1,i), &
                                  seltarget_icu_by_nuts2_icudates(num_actudays,i))
          enddo
        case("state")
        !** Calculate the intersection of the evaluated and observed ICU data by state and date --
          do i = 1, num_states
            call date_intersection(obs_lb(i), eval_lb(i), obs_ub(i), eval_ub(i), refer_date_mod, &
                                  time_n, num_actudays, seltarget_icu_by_state_icudates(1,i), &
                                  seltarget_icu_by_state_icudates(num_actudays,i))
          enddo
      end select
      if (PT_DEBUG) then
        print *, "Intersection is: "
        print *, "myrank: ",rank_mpi, obs_lb(1),eval_lb(1),obs_ub(1),eval_ub(1)
      endif
     
      !** Gather all the local population from all ranks to the root rank 0 before iteration
      call MPI_GATHER(local_ini_pop, opt_local_size*nvals, MPI_REAL8, ini_pop, opt_local_size*nvals, MPI_REAL8, ROOT, &
                MPI_COMM_WORLD, ierr)
      
      !!! ==============================================
      !!! Core of the genetic algorithm ================
      !!! ==============================================

      do gene_ii=1,opt_num_gene ! Iterate over the generations
        ga_res%opt_acuiter = gene_ii
    
        do pop_ii = 1,opt_local_size ! Iterate over all the population
          !** Proceed only when the corresponding fitness is not given (aka. INV) ----
          if (local_fitness(pop_ii) .ne. INV) then
            cycle
          endif
          !** Preprocess the R0 effect data according to the generated population ----
          do i = 1,nvals
            if (opt_week_index(i) .le. current_week) then
              R0_effects(opt_week_index(i),opt_index(i)) = local_ini_pop(i,pop_ii)
            endif
          enddo

          tmp_R0_effects = 0._rk
          
          if (auto_calibrate) then
            tmp_R0_effects(1:(end_window-start_window+1),:) = R0_effects(start_window:end_window,:)
          else
            tmp_R0_effects = R0_effects(start_week:end_week,:)
          endif
          !** Call the COVID19 spatial simuation and store the simulated ICU data to ill_ICU_cases_final -----
          if ( auto_calibrate .and. ((end_window - start_window + 1) .ne. opt_window_size) ) then
            !** Handle with the corner cases (initial and last)
            ill_ICU_cases_final = COVID19_Spatial_Microsimulation_for_Germany(iol,&
                iter_s, iter_e , &
                inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
                less_contagious, R0_force, immune_stop, &
                cor_R0change, R0delay ,R0delay_days, R0delay_type, &
                control_age_sex, seed_date, seed_before, sam_size, R0, &
                tmp_R0_effects, region_index, sim, sim_switch, output_dir, export_name, rank_mpi)
          else
            ill_ICU_cases_final = COVID19_Spatial_Microsimulation_for_Germany(iol,&
                iter_s, iter_e , &
                inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
                less_contagious, R0_force, immune_stop, &
                R0change, R0delay ,R0delay_days, R0delay_type, &
                control_age_sex, seed_date, seed_before, sam_size, R0, &
                tmp_R0_effects, region_index, sim, sim_switch, output_dir, export_name, rank_mpi)
          endif
          
          if (sim_switch%sim_reuse) then
            sim%t1      = sim_backup%t1
            sim%d       = sim_backup%d
          endif

          !** Fitness calculation ------------------------------------------------------------------
          call fitness_cal(local_fitness(pop_ii), opt_target_region, opt_target_icu, opt_filter, &
                            num_nuts2, num_states, pos_nuts2, pos_state, &
                            time_n, iter, coff_len, ill_ICU_cases_final, rank_mpi, &
                            seltarget_icu_by_state_icucases, seltarget_icu_by_nuts2_icucases, &
                            obs_lb, obs_ub, eval_lb)
          
      
          deallocate(ill_ICU_cases_final)
        enddo

        if (.Not.Allocated(sorted_index)) then
          Allocate(sorted_index(opt_pop_size))
          Allocate(fitness(opt_pop_size))
          Allocate(summary_fitness%maximum(opt_num_gene))
          Allocate(summary_fitness%mean(opt_num_gene))
          Allocate(summary_fitness%upper_hinge(opt_num_gene))
          Allocate(summary_fitness%median(opt_num_gene))
          Allocate(summary_fitness%lower_hinge(opt_num_gene))
          Allocate(summary_fitness%minimum(opt_num_gene))
        endif

        !** All ranks perform allgather operation on local_fitness for the upcoming stop criteria check
        !** This also works as a barrier
        call MPI_ALLGATHER(local_fitness, opt_local_size, MPI_REAL, fitness, opt_local_size, MPI_REAL,&
                MPI_COMM_WORLD, ierr)
        sorted_index = 0
      
        call SORTRX_REAL(opt_pop_size,fitness,sorted_index) ! Sorted_index indicates a sorted fitness in an ascending order

        if(PT_DEBUG) then
          print *,"the sorted_index is: "
          print *,sorted_index
        endif

        !** Add and print summary statistics --------------------------------------------------------
        call sixnum_summary_and_print(fitness, sorted_index, rank_mpi, gene_ii, &
                                     summary_fitness%maximum(gene_ii), & 
                                     summary_fitness%minimum(gene_ii), &
                                     summary_fitness%mean(gene_ii), &
                                     summary_fitness%median(gene_ii), &
                                     summary_fitness%upper_hinge(gene_ii), &
                                     summary_fitness%lower_hinge(gene_ii))

        !! TODO: To record the best fitness value and solution per generation (write them to a output file)

        !!! ==========================================================
        !!! ============ Check the stop criteria =====================
        !!! ==========================================================
        run_sum = 0
        max_pos = sorted_index(opt_pop_size)
        min_pos = sorted_index(1)

        if (gene_ii .eq. 1) then
          tmp_max = fitness(max_pos)
        endif
        if (gene_ii .gt. 1) then
          if (tmp_max .lt. summary_fitness%maximum(gene_ii)) then
            tmp_max = summary_fitness%maximum(gene_ii)
          endif
          do i = 1,gene_ii 
            if (summary_fitness%maximum(i) .ge. (tmp_max - EPS)) then
              run_sum = run_sum + 1
            endif
          enddo
        endif
        !** The genetic algorithm stops when -------------------------------------------------------
        !** 1. reaches the generation limit --------------------------------------------------------
        !** 2. the maximal fitnesses over the passed generations are equal with tolerable error ----
        !** 3. all the values in the fitness are identical -----------------------------------------
        if ((gene_ii .eq. opt_num_gene) .or. (run_sum .ge. opt_num_gene)&
            .or. (fitness(max_pos) .eq. fitness(min_pos))) then
          exit
        endif

        !** Only root rank 0 does the selection/crossover/mutation/elitism on the total population--
        if(rank_mpi .eq. ROOT) then 
          if ( .Not.Allocated(sorted_pop) ) then                
            Allocate(tmp_pop(nvals, opt_pop_size))
            Allocate(sorted_pop(nvals, opt_elitism))
            Allocate(tmp_fitness(opt_pop_size))
            Allocate(sorted_fitness(opt_elitism))
            Allocate(sel_par(opt_pop_size))       
          endif
          
          i = 1
          j = opt_pop_size
          pre_fitv = INV
          sorted_fitness = 0.0
          sorted_pop     = 0._rk
          !** Store the best fitness and population of opt_elitism to the array sorted_fitness and sorted_pop --
          !** Note that the repeated fitness and population is excluded ----------------------------------------
          do while(i .le. opt_elitism)
            if (j .eq. 0) then
              print *,"Too many duplicates lead to insufficient elitism"
              call exit(status)
            endif
            idx = sorted_index(j)
            curr_fitv = fitness(idx)
            if (curr_fitv .ne. pre_fitv) then
              sorted_fitness(i) = curr_fitv
              sorted_pop(:,i) = ini_pop(:,idx)
              i = i + 1
            endif
            j = j - 1
            pre_fitv = curr_fitv
          enddo
          if (PT_DEBUG) then
            print *,"elitism:"
            print *,sorted_fitness
            print *,sorted_pop
          endif

          tmp_fitness = fitness

          !** Select the individual population/chromosome having the highest fitness probability ------
          call selection(tmp_fitness, sel_par)

          if (PT_DEBUG) then
            print *,"The selected chromosomes: "
            print *,sel_par
          endif
    
          do i = 1,opt_pop_size
            tmp_pop(:,i) = ini_pop(:,sel_par(i))
            tmp_fitness(i) = fitness(sel_par(i))
          enddo   
          ini_pop = tmp_pop
          fitness = tmp_fitness
          if (PT_DEBUG) then
            print *,"Before crossover: "
            do i = 1,nvals
              print *,ini_pop(i,:)
            enddo
          endif

          !** The selected parents(pair) do the mating and bread two children to replace them ---------
          call crossover(ini_pop, fitness, opt_pcrossover)
          
          if(PT_DEBUG) then    
            print *,"After crossover: "
            do i = 1,nvals
              print *,ini_pop(i,:)
            enddo
            print *,fitness
          endif
          !** To randomly mutate a parameter/gene ------------------------------------------------------
          if(opt_pmutation .ne. 0.0) then
            call mutation(ini_pop, fitness, opt_pmutation, opt_lb, opt_ub)
          endif
          !** Replace the unfittest populations with the fittest ones, which are stored in sorted_pop ---
          call elitism(fitness, ini_pop, sorted_pop, sorted_fitness)
        endif
    
        !** Perfrom scatter operation before getting into the next generation ------
        call MPI_Scatter(fitness, opt_local_size, MPI_REAL, local_fitness, opt_local_size, MPI_REAL,&
                      ROOT, MPI_COMM_WORLD, ierr) ! Root scatters the updated fitness to all ranks
        call MPI_Scatter(ini_pop, opt_local_size*nvals, MPI_REAL8, local_ini_pop, opt_local_size*nvals, MPI_REAL8,&
                      ROOT, MPI_COMM_WORLD, ierr) ! Root scatters the updated population to all ranks
   
      enddo

      if (rank_mpi .eq. ROOT) then
        ! Store the best solution to "ga_res"
        ga_res%opt_bestsol      = ini_pop(:,max_pos)
        ga_res%opt_fitnessvalue = fitness(max_pos)

        call print_ga_summary(opt_names, ga_res, auto_calibrate, start_window, end_window)
      endif

      !** Root broadcasts the best solution to all other ranks ----
      call MPI_Bcast(ga_res%opt_bestsol, nvals, MPI_REAL8, ROOT, MPI_COMM_WORLD, ierr)
      !** Store the best solution back to R0_effects --------------
      do i = 1,nvals
        if (opt_week_index(i) .le. current_week) then
          R0_effects(opt_week_index(i),opt_index(i)) = ga_res%opt_bestsol(i)
        endif
      enddo
      
      !!! =============================================================================
      !!! Automation: the window moves backwards until the current week reaches =======
      !!! =============================================================================
      if (auto_calibrate) then
        tmp_R0_effects = 0._rk
        tmp_R0_effects(1:(end_window-start_window+1),:) = R0_effects(start_window:end_window,:)

        !** Write restart file before automation ends ----------------------
        if (end_window .ge. current_week) then 
          sim_switch%sim_write  = .true.
          sim_switch%write_week = start_window + size(sim_R0change,dim=2)
        endif

        !** Call the COVID19 spatial simuation once again, its resulting sim data serves as the --
        !** input to the next-round calibration -------------------------------------------------- 
        if (size(sim_R0change,dim=2) .ne. 0) then        
          ill_ICU_cases_final = COVID19_Spatial_Microsimulation_for_Germany(iol,&
              1, 1 , &
              inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
              less_contagious, R0_force, immune_stop, &
              sim_R0change, R0delay ,R0delay_days, R0delay_type, &
              control_age_sex, seed_date, seed_before, sam_size, R0, &
              tmp_R0_effects, region_index, sim, sim_switch, output_dir, export_name, rank_mpi)
          sim_backup%t1           = sim%t1
          sim_backup%d            = sim%d
          sim_switch%sim_reuse    = .true.
          sim_switch%sim_reload   = .false.
        endif
      endif
             
      !** The automation stop criteria reaches ---------------
      if (end_window .ge. current_week) then
        start_window = end_window+1
        cycle
      endif

      local_fitness = INV ! Each rank invalidates its local fitness again

      !** Prepare for the next-round calibration ---------------------------------
      !** , when the initial window size is unequal to the given window size -----
      if((end_window - start_window) .ne. (opt_window_size - 1)) then
        refer_date = add_date(refer_date,WEEKLEN*(diff_weeks-(opt_window_size-opt_shift_size)))
        deallocate(cor_R0change)
        deallocate(sim_R0change)
        Allocate(sim_R0change(size(R0change,dim=1),opt_shift_size))
        call ini_R0change(sim_R0change)
        deallocate(tmp_R0_effects)
        Allocate(tmp_R0_effects(opt_window_size,R0_effects_dim2size))
        nvals = opt_window_size*R0_effects_dim2size
        deallocate(local_ini_pop)
        deallocate(ini_pop)
        Allocate(local_ini_pop(nvals,opt_local_size))
        Allocate(ini_pop(nvals,opt_pop_size))
        if (rank_mpi .eq. ROOT) then
          if (allocated(tmp_pop)) then
            deallocate(tmp_pop)
            deallocate(sorted_pop)
          
            Allocate(tmp_pop(nvals,opt_pop_size))
            Allocate(sorted_pop(nvals,opt_elitism))
          endif
        endif
        deallocate(opt_week_index)
        deallocate(opt_index)
        
        Allocate(opt_week_index(nvals))
        Allocate(opt_index(nvals))
        deallocate(ga_res%opt_bestsol)
        Allocate(ga_res%opt_bestsol(nvals))
        do i = 1,nvals
          opt_index(i) = CEILING(i/Real(opt_window_size))
        enddo
      else
        refer_date = add_date(refer_date,WEEKLEN*opt_shift_size)
      endif

      !** The indexes change as the window moves ---------------------------
      do i = 1, opt_window_size
        opt_week_index(i) = i + end_window - (opt_window_size-opt_shift_size)
      enddo
      do i = opt_window_size+1,nvals
        opt_week_index(i) = opt_week_index(i-opt_window_size)
      enddo
      
      !** Population generation --------------------------
      do j = 1, opt_local_size
        do i = 1, nvals
          ! The random number is between opt_lb and opt_ub.
          local_ini_pop(i,j) = random_uniform(opt_lb,opt_ub)
        enddo
      enddo
      refer_date_mod = add_date(refer_date,1)

      end_window   = end_window + opt_shift_size
      start_window = end_window - opt_window_size + 1
      if (end_window .gt. current_week) then
        end_window = current_week
      endif

      time_n = maxval(R0change) + 1

      !** Prepare for the next-round calibration -----------------------------
      !** , when the lost window size is less than the given window size -----
      if ((end_window - start_window) .lt. (opt_window_size-1)) then
        m = end_window - start_window + 1 ! m is larger than 0
        if (Allocated(cor_R0change)) then
          deallocate(cor_R0change)
        endif
        Allocate(cor_R0change(size(R0change,dim=1),m))
        call ini_R0change(cor_R0change)
        if (Allocated(sim_R0change)) then
          deallocate(sim_R0change)
        endif
        if ( m-(opt_window_size-opt_shift_size) .gt. 0 ) then 
          Allocate(sim_R0change(size(R0change,dim=1),m-(opt_window_size-opt_shift_size)))
        else
          Allocate(sim_R0change(size(R0change,dim=1),m-1))
        endif
        if (size(sim_R0change,dim=2) .ne. 0) then
          call ini_R0change(sim_R0change)
        endif
        time_n = maxval(cor_R0change) + 1
      endif
    enddo
    
    !** Write the updated R0 value to the resulting file ---------------------
    if(rank_mpi .eq. ROOT) then
      write(*,'(A)')"Write the calibrated R0effects to a restart file..."
      call writeback_R0effects(auto_calibrate, opt_target_region, &
                                  opt_names, ga_res, &
                                  R0_effects, iol%R0_effect%head, start_week)
    endif

    !!! ===================================
    !!! Deallocation ======================
    !!! ===================================
    select case (trim(opt_target_region))
      case("nuts2")
        deallocate(seltarget_icu_by_nuts2_icucases)
        deallocate(seltarget_icu_by_nuts2_icudates)
        deallocate(obsicu_nuts2_date)
        deallocate(obsicu_nuts2_cases)
        deallocate(obsicu_counties)
        deallocate(unique_nuts2s)
        deallocate(pos_nuts2)
        deallocate(counties_dist_id)

      case("state")
        deallocate(seltarget_icu_by_state_icucases)
        deallocate(seltarget_icu_by_state_icudates)
        deallocate(obsicu_date)
        deallocate(obsicu_sc)
        deallocate(states_shortcut)
        deallocate(uniq_distid)
        deallocate(pos_state)
      case default
        write(*,*)  "Otherwise, the given target is not supported yet"
        call exit(status) 
    end select
    
    deallocate(obs_ub)
    deallocate(obs_lb)
    deallocate(eval_lb)
    deallocate(eval_ub)
    deallocate(ini_pop)
    deallocate(fitness)
    deallocate(local_fitness)
    deallocate(local_ini_pop)
    deallocate(sorted_index)
  
    if(rank_mpi .eq. ROOT) then
      if (allocated(sorted_pop)) then
        deallocate(sorted_pop)
        deallocate(sorted_fitness)
        deallocate(sel_par)
        deallocate(tmp_fitness)
        deallocate(tmp_pop)
      endif
    endif

    deallocate(summary_fitness%maximum)
    deallocate(summary_fitness%mean)
    deallocate(summary_fitness%upper_hinge)
    deallocate(summary_fitness%median)
    deallocate(summary_fitness%lower_hinge)
    deallocate(summary_fitness%minimum)
 
  end subroutine 

  !! --------------------------------------------------------------------------------------------------------------
  !> Subroutine that initializes the R0 change array with given dimension size
  !>
  subroutine ini_R0change(array_inout)
    Integer(kind=ik),allocatable,dimension(:,:),intent(inout)  :: array_inout
    Integer                         :: i, j
    Integer                         :: array_len_dim1, array_len_dim2

    array_len_dim1 = size(array_inout, dim=1)
    array_len_dim2 = size(array_inout, dim=2)

    array_inout(1,1) = 1
    array_inout(array_len_dim1,1) = WEEKLEN
    do i = 1,array_len_dim1
      do j = 2,array_len_dim2
        array_inout(i,j) = array_inout(i,j-1) + WEEKLEN
      enddo 
    enddo
  end subroutine ini_R0change

  !! ------------------------------------------------------------------------------
  !> Subroutine that applies the convolution filter to the input array.
  !>
  subroutine conv_filter(array_filter, coff_len, array_size)
    Real,dimension(array_size), intent(inout)    :: array_filter
    Integer,intent(in)                           :: coff_len, array_size

    Real                              :: coefficient
    Integer                           :: k, j, halve_coff_len, startid
    Real,dimension(:),Allocatable     :: tmp_filter
    coefficient = 1.0/coff_len
    Allocate(tmp_filter(array_size-coff_len+1))
    halve_coff_len = coff_len / 2
   
    if (mod(coff_len,2) .eq. 0) then
      startid = halve_coff_len
    else
      startid = halve_coff_len + 1
    endif

    k = 1
  
    do j = startid,array_size-halve_coff_len
      tmp_filter(k) = sum(array_filter((j-startid+1):(j+halve_coff_len)))*coefficient
      k = k + 1                  
    enddo
    array_filter(startid:(array_size-halve_coff_len))= tmp_filter
    deallocate(tmp_filter)
  end subroutine conv_filter

  !! --------------------------------------------------------------------------------------------------------------
  !> Subroutine that get a subset of the input array_in (with repeated nuts2/state) storing the unique nuts2/state.
  !> 
  !> The returned arrays are array_pos and array_out.
  !> array_out: store the unique nuts2/state; array_pos: store the position from which a different nuts2/state starts in array_in
  function get_unique_nuts2_state(array_in, array_pos) result(array_out)
    integer                                               :: size_in,j,index
    Integer,Allocatable,dimension(:),intent(inout)        :: array_pos
    Integer(kind=ik),allocatable,dimension(:),intent(in)  :: array_in
    Integer,allocatable,dimension(:)                      :: array_out
    
    if (size(array_in) == 0)then
      allocate(array_out(0))
      return
    end if
    
    index = 1
    size_in = size(array_in)
    !** Detemine dimension of array_out
    do j=2,size_in
      if (array_in(j)/=array_in(j-1)) then
        index = index + 1
      endif
    enddo
    !allocate array and reset index number
    allocate(array_out(index))
    allocate(array_pos(index+1))
    index = 1
    array_out(index) = array_in(index)  
    array_pos(index) = 1
    do j=2,size_in
      if (array_in(j)/=array_in(j-1)) then
        index = index + 1
        array_out(index) = array_in(j)
        array_pos(index) = j
      endif
    enddo
    array_pos(index+1) = size_in+1
  
  end function get_unique_nuts2_state

  !! -------------------------------------------------------------------------------------------
  !> Subroutine that calculate the intersection of the evaluated and observed ICU data by date.
  !>
  !> The intersection is represented by the returned intervals [obs_lb, obs_ub] and [eval_lb, eval_ub].
  subroutine date_intersection(obs_lb, eval_lb, obs_ub, eval_ub, refer_date_mod, time_n, length, s_date, e_date)
    Integer,intent(inout)                                  :: obs_lb, obs_ub, eval_lb, eval_ub
    Integer,intent(in)                                     :: length
    Character(len=*)                                       :: s_date, e_date
    Integer                                                :: inteval_days, status, time_n
    Character*10                                           :: refer_date_mod

    !** Here we assume that the date is monotonically increased by days -----------------------
    obs_lb = 1
    eval_lb = 1

    !TODO: Move the calculation of inteval_days to support module
    inteval_days = Date2Unixtime(trim(s_date))&
                        - Date2Unixtime(refer_date_mod)
    inteval_days = inteval_days/86400

    if (inteval_days .LT. 0) then
      obs_lb = obs_lb - inteval_days
    elseif (inteval_days .GT. 0) then
      eval_lb = eval_lb + inteval_days
    endif

    obs_ub = length
    eval_ub = time_n

    inteval_days = (Date2Unixtime(trim(e_date))& 
        - (Date2Unixtime(refer_date_mod)+86400*(time_n-1)))/86400
    if (inteval_days .LT. 0) then
      eval_ub = eval_ub + inteval_days
    elseif (inteval_days .GT. 0) then
      obs_ub = obs_ub - inteval_days
    endif

    !** Throw an error when the intersection is null ---------------------------
    !!TODO: give the hint that gives the least valid current week.
    if ((eval_lb .GT. eval_ub) .OR. (obs_lb .GT. obs_ub) .OR.&
     ((eval_ub - eval_lb) .NE. (obs_ub - obs_lb))) then
      print *, "Error: intersection is null!"
      print *, "eval_lb: ",eval_lb, " eval_ub: ",eval_ub," obs_lb: ",obs_lb," obs_ub: ",obs_ub
      Call exit(status)    
    endif
  end subroutine date_intersection

  !! --------------------------------------------------------------------------------------------------------------
  !> Subroutine that calculates the fitness of certain solution.
  !> 
  !> ill_ICU_cases_final : the evaluated ICU cases;
  !> (seltarget_icu_by_state_icucases, seltarget_icu_by_nuts2_icucases) : the observed ICU cases
  !> The returned value is val_localfitness.
  subroutine fitness_cal(val_localfitness, opt_target_region, opt_target_icu, opt_filter, &
                            num_nuts2, num_states, pos_nuts2, pos_state, &
                            time_n, iter, coff_len, ill_ICU_cases_final, rank_mpi, &
                            seltarget_icu_by_state_icucases, seltarget_icu_by_nuts2_icucases, &
                            obs_lb, obs_ub, eval_lb)
    Real, intent(out)                                   :: val_localfitness
    character(Len=:), Allocatable, intent(in)           :: opt_target_region
    Logical, intent(in)                                 :: opt_target_icu, opt_filter ! .true.
    Integer, intent(in)                                 :: num_nuts2, num_states
    Integer, Allocatable,dimension(:), intent(in)       :: pos_nuts2, pos_state
    Integer, intent(in)                                 :: time_n, iter, coff_len
    Integer, Allocatable,dimension(:,:,:), intent(in)   :: ill_ICU_cases_final ! Store the evaluated ICU cases
    Integer(kind=mpi_ik), intent(in)                    :: rank_mpi
    Real,Allocatable,dimension(:,:), intent(in)         :: seltarget_icu_by_state_icucases
    Real,Allocatable,dimension(:,:), intent(in)         :: seltarget_icu_by_nuts2_icucases
    Integer, Allocatable,dimension(:), intent(in)       :: obs_lb, obs_ub, eval_lb

    
    Real                                          :: diff_sum, diff_loc_sum, diff_tmp_power
    Integer                                       :: diff_len, diff_loc_len
    Integer, Allocatable,dimension(:,:,:)         :: aggreg_sum_by_county
    Real, Allocatable,dimension(:,:)              :: avg_by_county ! Record the simulated ICU cases by state/nuts2
    Integer                                       :: i, j, k, status
    
    diff_sum = 0.0
    diff_len = 0

    !==================================================!
    !====== Calculate fitness function value ==========!
    !==================================================!

    select case (trim(opt_target_region))
      case("nuts2")  
        !** All the data are represented by nuts2 --------
       
        Allocate(aggreg_sum_by_county(num_nuts2,time_n,iter))
        Allocate(avg_by_county(num_nuts2,time_n))
        
        aggreg_sum_by_county = 0
        avg_by_county = 0.0

        !** TODO: The fitness relavant to the death data ----------
        !** The fitness relavant to the icu data ------------------
        if (opt_target_icu) then
          do j = 1,iter
            do i = 1,time_n 
              do k = 1,num_nuts2
                aggreg_sum_by_county(k,i,j) = sum(ill_ICU_cases_final(pos_nuts2(k):(pos_nuts2(k+1)-1),i,j)) ! Sum over all the counties within one nuts2
              enddo
            enddo
          enddo

          do j = 1,time_n
            do i = 1,num_nuts2
              avg_by_county(i,j) = Real(sum(aggreg_sum_by_county(i,j,:)))/Real(iter) ! Average over the runing iterates
            enddo
          enddo

          if (PT_DEBUG) then
            do j = 1,time_n 
              print *,avg_by_county(:,j)
            enddo
          endif
            !** Calculate the difference between the observed and evaluated/simulated ICU data -----------
          do i = 1,num_nuts2                  
            if (opt_filter) then
              call conv_filter(avg_by_county(i,:), coff_len, time_n)                   
            endif
            diff_loc_sum = 0.0
            diff_loc_len = 0
            j = obs_lb(i)
            k = eval_lb(i)

            do while (j .le. obs_ub(i))
              if (seltarget_icu_by_nuts2_icucases(j,i) .ne. NULVAL) then
                diff_tmp_power = (seltarget_icu_by_nuts2_icucases(j,i) - avg_by_county(i,k))**2
                diff_loc_sum = diff_loc_sum + diff_tmp_power
                diff_loc_len = diff_loc_len + 1
              endif
              j = j + 1
              k = k + 1
            end do
            diff_sum = diff_sum + diff_loc_sum
            diff_len = diff_len + diff_loc_len  
          enddo   
          val_localfitness = -1*sqrt(diff_sum/diff_len)
        end if
      case("state")
        !** All the data are represented by state --------
        if (.not.allocated(aggreg_sum_by_county)) then
          Allocate(aggreg_sum_by_county(num_states,time_n,iter))
          Allocate(avg_by_county(num_states,time_n))
        endif
        aggreg_sum_by_county = 0
        avg_by_county = 0.0

        !** TODO: The fitness relevant to the death data ----------
        !** The fitness relavant to the icu data ------------------
      
        ! TODO: record the number of the affected states
        do j = 1,iter
          do i = 1,time_n 
            do k = 1,num_states
              aggreg_sum_by_county(k,i,j) = sum(ill_ICU_cases_final(pos_state(k):pos_state(k+1)-1,i,j)) ! Sum over all the counties within one state
            enddo
          enddo
        enddo
  
        do j = 1,time_n
          do i = 1,num_states
            avg_by_county(i,j) = Real(sum(aggreg_sum_by_county(i,j,:)))/Real(iter) ! Average over the running iterates
          enddo
        enddo
        if (PT_DEBUG) then
          do j = 1,time_n 
            print *,avg_by_county(:,j)
          enddo
        endif

        !** Calculate the difference between the observed and evaluated/simulated ICU data -----------
        do i = 1,num_states                  
          if (opt_filter) then 
            call conv_filter(avg_by_county(i,:), coff_len, time_n)
          endif

          diff_loc_sum = 0.0
          diff_loc_len = 0
          j = obs_lb(i)
          k = eval_lb(i)

          do while (j .le. obs_ub(i))
            if (seltarget_icu_by_state_icucases(j,i) .ne. NULVAL) then
              diff_tmp_power = (seltarget_icu_by_state_icucases(j,i) - avg_by_county(i,k))**2
              diff_loc_sum = diff_loc_sum + diff_tmp_power
              diff_loc_len = diff_loc_len + 1
            endif
            j = j + 1
            k = k + 1
          end do
          diff_sum = diff_sum + diff_loc_sum
          diff_len = diff_len + diff_loc_len    
        enddo     
        val_localfitness = -1*sqrt(diff_sum/diff_len)
      case default
        write(*,*)  "Otherwise, the given target is not supported yet"
        call exit(status)
    end select
    deallocate(aggreg_sum_by_county)
    deallocate(avg_by_county) 

  end subroutine fitness_cal

  !! ------------------------------------------------------------------------------
  !> Subroutine that selects the fittest parents
  !>
  !> Return the result indexes to sel_par.
  subroutine selection(fitness, sel_par)
    Real, Allocatable, Dimension(:),intent(inout)          :: fitness
    Integer, Dimension(size(fitness)),intent(inout)        :: sel_par

    Real, Dimension(size(fitness))                       :: fscaled, prob
    Real                                                 :: fmin, fave, fmax
    Real                                                 :: delta, a, b
    Integer                                              :: sfactor ! 2

    sfactor = 2
    fmin = minval(fitness)
    if (fmin .lt. 0.0) then
      fitness(:) = fitness(:) - fmin
      fmin = minval(fitness)
    end if
    fmax = maxval(fitness)
    fave = sum(fitness)/size(fitness)

    if (fmin .gt. (sfactor*fave - fmax)/(sfactor - 1)) then
      delta = fmax - fave
      a = (sfactor - 1.0)*fave/delta 
      b = fave * (fmax - sfactor*fave)/delta 
    else
      delta = fave - fmin
      a = fave/delta 
      b = -1*fmin*fave/delta 
    endif
    fscaled = a*fitness(:) + b
    prob = abs(fscaled)/sum(abs(fscaled)) ! Get the fitness probability
    ! Sample the indexes according to the weights indicated in prob
    sel_par = sample_weight(size(fitness), prob)    
  end subroutine selection

  !! ------------------------------------------------------------------------------
  !> Subroutine that mates the selected pairs and breed two children for each
  !>
  !> The results are stored back to ini_pop and fitness.
  subroutine crossover(ini_pop, fitness, pcrossover) 
    Real,Allocatable,dimension(:),intent(inout)                   :: fitness
    Real(kind=rk),Allocatable, dimension(:,:),intent(inout)       :: ini_pop
    
    Real,intent(in)                                            :: pcrossover

    Integer,dimension(:),Allocatable                           :: sample_input, sample_parents
    Integer                                                    :: start_id, end_id, nmating,i,j,len
    Real(kind=rk)                                              :: ran
    Real(kind=rk),dimension(:),Allocatable                     :: lchild, rchild

    nmating = floor(size(fitness)/2.0) !get the floor
    len = size(ini_pop,dim=1)
    allocate(sample_input(nmating*2))
    allocate(sample_parents(nmating*2))
    allocate(lchild(len))
    allocate(rchild(len))

    do i = 1, 2*nmating 
      sample_input(i) = i
    enddo
  
    sample_parents = sample_i4(sample_input,2*nmating) ! Sample the parent pairs
    if (PT_DEBUG) then
      print *,"The crossovered pairs: "
      print *,sample_parents
    endif
    
    do i = 1, nmating 
      call random_number(ran)
      end_id = i * 2 
      start_id = end_id - 1 
      start_id = sample_parents(start_id)
      end_id = sample_parents(end_id)

      !** Also skip crossover when the parents are identical
      if((pcrossover > ran) .and.&
       (fitness(start_id) .ne. fitness(end_id))) then
        call random_number(ran)
        !** Mating ----------
        do j=1,len
          lchild(j) = ran*ini_pop(j,start_id) + (1._rk-ran)*ini_pop(j,end_id)
          rchild(j) = (1._rk-ran)*ini_pop(j,start_id) + ran*ini_pop(j,end_id)

          ini_pop(j,start_id) = lchild(j)
          ini_pop(j,end_id)   = rchild(j)
        enddo

        !** Invalid the relevant fitness ----
        fitness(start_id) = INV 
        fitness(end_id)   = INV 
      endif
    enddo
    deallocate(sample_input)
    deallocate(sample_parents)
    deallocate(lchild)
    deallocate(rchild)
  end subroutine crossover

  !! ------------------------------------------------------------------------------
  !> Subroutine that randomly mutates certain parameter
  !>
  !> The results are stored back to ini_pop and fitness.
  subroutine mutation(ini_pop, fitness, pmutation, lb, ub)
    Real,Allocatable,dimension(:),intent(inout)                :: fitness
    Real(kind=rk),Allocatable,dimension(:,:),intent(inout)     :: ini_pop
    
    Real,intent(in)                                         :: pmutation
    Real, intent(in)                                        :: lb, ub

    Integer,dimension(size(ini_pop,dim=1))                  :: sample_input
    Integer,dimension(1)                                    :: sample_parents
    Integer                                                 :: i
    Real                                                    :: ran
    Real(kind=rk)                                           :: Res

    do i = 1,size(sample_input)
      sample_input(i) = i
    enddo
    do i = 1,size(fitness)
      call random_number(ran)
      if (pmutation > ran) then
        sample_parents = sample_i4(sample_input,size(sample_parents)) ! Randomly generate the mutated gene/parameter
        Res = random_uniform(lb,ub)
        ini_pop(sample_parents(1),i) = Res  
      !  print *,"Mutation: chromosome ",i,",gene: ",sample_parents(1),",res: ", Res
        fitness(i) = INV ! Invalidate the mutated chromosome/population individual
      endif
    enddo
  end subroutine mutation

  !! ------------------------------------------------------------------------------
  !> Subroutine that replaces the unfittest individules with the chosen elitim
  !>
  !> The results are stored back to ini_pop and fitness.
  subroutine elitism(fitness, ini_pop, sorted_pop, sorted_fitness)
    Real(kind=rk),Allocatable,dimension(:,:),intent(inout)        :: ini_pop
    Real,Allocatable,dimension(:),intent(inout)                   :: fitness
    Real(kind=rk),dimension(:,:),Allocatable,intent(in)           :: sorted_pop                                   
    Real,dimension(:),Allocatable,intent(in)                      :: sorted_fitness
    Integer,dimension(size(fitness))                              :: sorted_index
    Integer                                                       :: num

    call SORTRX_REAL(size(fitness),fitness,sorted_index) ! Fitness in an ascending order
    num = size(sorted_fitness)
    !** Replace the unfittest individules with the elitism chosen in the first place ---
    !** The elitism is represented as sorted_fitness and sorted_pop
    fitness(sorted_index(1:num)) = sorted_fitness
    ini_pop(:,sorted_index(1:num)) = sorted_pop
    
    if(PT_DEBUG) then
      print *,fitness 
      do num = 1,size(ini_pop,dim=1)
        print *,ini_pop(num,:)
      enddo
    endif
  end subroutine elitism

  !! ------------------------------------------------------------------------------
  !> Function that returns a sample with the weight in mind
  !>
  !> The weights/probabilities are stored in prop.
  !> The high the weight is, the most possibility it is chosen
  function sample_weight(input, prop) Result(sel_par)
    Integer                  :: input
    Real,dimension(input)    :: prop
    Integer,dimension(input) :: sel_par

    Integer                  :: idx,i 
    Real                     :: sample_rmd, acum_prop

    do i = 1,input
      call random_number(sample_rmd)

      acum_prop = 0.0
      do idx = 1,input
        acum_prop = acum_prop + prop(idx)
        if(sample_rmd .le. acum_prop) then
          exit
        end if
      end do
      sel_par(i) = idx
    end do
  end function

  !! ------------------------------------------------------------------------------
  !> Function that returns a six number summary statistics
  !>
  !> This adds a mean value on top of the existing five number summary
  subroutine sixnum_summary_and_print(fitness, sorted_index, rank_mpi, gene_ii,&
                             max, min, mean, median, upper_hinge, lower_hinge)
    Real,Allocatable,dimension(:),intent(in)         :: fitness
    Integer,Allocatable,dimension(:),intent(in)      :: sorted_index
    Integer(kind=mpi_ik), intent(in)                 :: rank_mpi
    Integer(kind=ik), intent(in)                     :: gene_ii
    Real, intent(out)                                :: max, min, mean, median, upper_hinge, lower_hinge

    Integer                                     :: half, i, j, length

    length = size(fitness)
    max = fitness(sorted_index(length))
    min = fitness(sorted_index(1))
    mean = Sum(fitness)/length

    median = get_median(fitness, sorted_index, length)
    half = length/2.0
    upper_hinge = get_median(fitness, sorted_index((half+1):length), length-half) ! Get the median of the second half
    lower_hinge = get_median(fitness, sorted_index(1:half), half) ! Get the median of the first half

    if (rank_mpi .eq. ROOT) then
      !** Output summary statistics --------------------------------------------------------
      print *,"intermediate fitness output:"
      print *,fitness
      print *,"generation: ",gene_ii
      write (*, '(A9,A18,A20,A15,A20,A12)')"max","mean","upper_hinge","median","lower_hinge","min"
      print *,max,mean,upper_hinge,median,lower_hinge,min
    endif
  end subroutine sixnum_summary_and_print
 
  !! ------------------------------------------------------------------------------
  !> Function that returns the median value of an ordered array
  !>
  function get_median(dataset, sorted_index, len) Result(med)
    Real,Allocatable,dimension(:),intent(in)       :: dataset
    Integer,intent(in)                             :: len
    Integer,dimension(len),intent(in)              :: sorted_index
    Real                                           :: med

    Integer                                        :: length, half

    length = size(sorted_index)
    half = length/2

    if (mod(length,2) .eq. 0) then
      med = (dataset(sorted_index(half))+dataset(sorted_index(half+1)))/2.0
    else
      med = dataset(sorted_index(half+1))
    endif
  end function get_median

  !! ------------------------------------------------------------------------------
  !> Subroutine that prints the value of generic GA parameters
  !> This subroutine is only called once
  !>
  subroutine print_ga_header(opt_pop_size, opt_num_gene, opt_elitism, opt_pcrossover, opt_pmutation, &
                  opt_target_region)  
    character(len=*), intent(in)                              :: opt_target_region ! Optimized region. Valid values: "state" and "nuts2"
    Integer, intent(in)                                       :: opt_pop_size ! 4 Population size
    Integer, intent(in)                                       :: opt_num_gene ! 2 Maximal number of generations
    Integer, intent(in)                                       :: opt_elitism ! The number of the chosen elitism
    Real, intent(in)                                          :: opt_pcrossover ! 0.8
    Real, intent(in)                                          :: opt_pmutation ! 0.1
    Integer                                                   :: i

    write(*,'(A)')"!############ Launch Genetic Algorithm #############!"
    write(*,'(A)')"!---------------- GA settings ---------------!"
    print *,"Type: ","real-valued"
    print *,"Population size: ",opt_pop_size
    print *,"Number of generations: ",opt_num_gene ! maxiter
    print *,"Elitism: ",opt_elitism
    print *,"Crossover probability: ",opt_pcrossover
    print *,"Mutation probability: ",opt_pmutation
    print *,"Optimization target: ",trim(opt_target_region)
    write(*,'(A)')"!--------------------------------------------!"
  end subroutine print_ga_header

  !! ------------------------------------------------------------------------------
  !> Subroutine that prints the GA results containing best solution
  !>
  subroutine print_ga_summary(opt_names, ga_res, auto_cali, start_window, end_window)  
    character(len=:), Allocatable, Dimension(:), intent(in)   :: opt_names ! Names of parameters (by state or Nuts-2) to optimize
    type(opt_res),intent(in)                                  :: ga_res
    Logical,intent(in)                                        :: auto_cali
    Integer, intent(in)                                       :: start_window, end_window
    Integer                                                   :: i

    write(*,'(A)')"!----------------- GA results ----------------!"
    print *,"Iterations: ",ga_res%opt_acuiter
    print *,"Fitness function value: ",ga_res%opt_fitnessvalue
    
    if (.not.auto_cali) then
      print *,"Affected targets:"
      do i = 1,size(opt_names)
        write(*, fmt="(1x,a)", advance="no") opt_names(i)
      enddo
      write(*,*)
    else
      print *,"By default all the targets are affected!"
      write(*, fmt="(1x,a)",advance="no") "Optimize the R0_effects value for weeks:"
      do i = start_window, end_window
        write(*, fmt="(1x,a,I3.3)", advance="no") "W",i
      enddo
    endif
    write(*,*)
    write(*,'(A)')"!---------------------------------------------!"
  end subroutine print_ga_summary

  !! ------------------------------------------------------------------------------
  !> Subroutine that write the updated R0 values to the given file
  !>
  !> There are two different kinds of file:
  !> 1. for auto_calibrate (include prognosis); 2. one-off calibration
  !> 
  Subroutine writeback_R0effects(auto_calibrate, opt_target_region, &
                                  opt_names, ga_res, &
                                  R0_effects, head, start_week)

    Logical,intent(in)                                       :: auto_calibrate
    character(len=*), intent(in)                             :: opt_target_region
    Character(len=:),Allocatable,dimension(:),intent(in)     :: opt_names
    Type(opt_res),intent(in)                                 :: ga_res
    Real(kind=rk),    Allocatable, Dimension(:,:),intent(in) :: R0_effects
    character(len=:), Allocatable, Dimension(:)  ,intent(in) :: head
    Integer                                      ,intent(in) :: start_week
    
    logical                                            :: exs
    Integer                                            :: un, i,j
    character(3)                                       :: c_index
    character(len=40)                                  :: filename
 
    !---------------------------------------------------------------------------
    exs = .FALSE.

    If (auto_calibrate) then
      write(filename,'("R0effects_prog_update.auto_",A,".csv")')trim(opt_target_region)
    Else
      write(filename,'("R0effects_calibrate_update.",A,".csv")')trim(opt_target_region)
    Endif
    
    Inquire(file = trim(filename), exist=exs)
    
    If (exs) then
      Open (newunit=un, file = trim(filename),&
            position="Append", status ="old")
    Else
      Open (newunit=un, file = trim(filename),&
            status ="new")
      If (auto_calibrate) then
        write(un,'(*(3A))')('"',trim(head(i)),'" ',i=1,size(head)) 
      Endif
    End If
    If (auto_calibrate) then
      do i = start_week, size(R0_effects,dim=1)
        write(c_index,'(i3.3)') i
        write(un,'(3A)',ADVANCE="NO")'"',c_index,'"'
        write(un,'(*(F20.17))')(R0_effects(i,j),&
              j=1,size(R0_effects,dim=2))
      enddo
    Else
      write(un,'(*(3A))')('"',trim(opt_names(i)),'" ',i=1,size(opt_names))
      write(un,'(*(F20.17))')(ga_res%opt_bestsol(i),&
              i=1,size(ga_res%opt_bestsol))
    Endif
    Close(un)
  End Subroutine writeback_R0effects
End Module genetic_algorithm
