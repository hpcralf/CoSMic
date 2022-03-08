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
! Module containing the GA-related routines
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

  implicit none

  !** Constants ---------------------------------------
  Real, parameter     :: INV    = 99.0 ! Corresponds to NA for fitness array
  Real, parameter     :: EPS    = EPSILON(1.0) ! Tolerance precision, a magnitude of e-07

  Real, parameter     :: NULVAL = -1.0 ! Corresponds to NA for icu_nuts2_cases
  Integer,parameter   :: nvals  = 35 ! The total number of parameters to be optimized

  !** Extend the fivenum summary statistic with a mean value -------
  type ga_summary
    Real,dimension(:),allocatable                 :: maximum
    Real,dimension(:),allocatable                 :: mean 
    Real,dimension(:),allocatable                 :: upper_hinge
    Real,dimension(:),allocatable                 :: median
    Real,dimension(:),allocatable                 :: lower_hinge
    Real,dimension(:),allocatable                 :: minimum
  end type ga_summary

  !** Optimization Control Panel / setting the optimization --------
!  type opt_parameters     
!    Logical                                       :: opt_target_icu
!    Logical                                       :: opt_target_deaths ! Optimization targets
!    Logical                                       :: opt_filter ! .TRUE.
!    character(len=mcl)                            :: use_sug_sol ! "NULL"
!    character(len=8)                              :: opt_target_region ! Optimized region. Valid values: "state" and "nuts2"
!    character(len=8), Dimension(4)                :: opt_names ! Names of parameters (by state or Nuts-2) to optimize
!    Integer,Dimension(4)                          :: opt_names_dur ! The duration of certain parameter is calculated by weeks
!    Real                                          :: opt_lb ! 0.1 Lower bounds of optimized parameters 
!    Real                                          :: opt_ub ! 1.0 Upper bounds of optimized parameters
!    Integer                                       :: opt_pop_size ! 4 Population size
!    Integer                                       :: opt_num_gene ! 2 Maximal number of generations
!    Integer                                       :: opt_elitism ! The number of the chosen elitism
!    Real                                          :: opt_pcrossover ! 0.8
!    Real                                          :: opt_pmutation ! 0.1
!  end type opt_parameters

  !** GA results/best solution -------------------------------------
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
     R0_effects, full_region_index, region_index, output_dir, export_name, rank_mpi) !type opt_parameters

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
    Real(kind=rk),    Allocatable, Dimension(:,:)             :: R0_effects
    Integer(kind=mpi_ik)                         , intent(in) :: rank_mpi
    integer(kind=ik), Allocatable, Dimension(:)  , intent(in) :: full_region_index
    integer(kind=ik), Allocatable, Dimension(:)  , intent(in) :: region_index
    character(len=:), Allocatable                , intent(in) :: output_dir
    character(len=:), Allocatable                , intent(in) :: export_name

    !==========================================================================
    ! General purpose for preprocessing and fitness calculation
    !==========================================================================     
    Real(kind=rk),Allocatable                     :: ini_pop(:,:), tmp_pop(:,:), sorted_pop(:,:) ! ini_pop: store the population, which are the input to be optimized. 

    Real, Allocatable                             :: test_fitness(:), fitness(:), tmp_fitness(:) ! fitness: store the calculated fitness
    Real, Allocatable, Dimension(:)               :: sorted_fitness
    Type(ga_summary)                              :: summary_fitness

    Integer,dimension(:,:,:), Allocatable         :: ill_ICU_cases_final ! Store the evaluated ICU cases
    Integer(kind=ik)                              :: i, j, k, m, pop_ii, gene_ii
    Integer                                       :: startid, endid, max_pos, min_pos
    Real                                          :: tmp_max, run_sum
    Integer                                       :: time_n
      
    Real,dimension(:), Allocatable                :: diff
    Real                                          :: diff_sum, diff_loc_sum, diff_tmp_power
    integer                                       :: diff_len, diff_loc_len
    integer                                       :: half
    Integer,dimension(:,:,:), Allocatable         :: aggreg_sum_by_county
    Real,dimension(:,:), Allocatable              :: avg_by_county ! Record the simulated ICU cases by state/nuts2

    Character(Len=:), allocatable, Dimension(:)   :: obsicu_nuts2_date
    Integer(kind=ik), Allocatable, Dimension(:)   :: obsicu_nuts2_cases
    Integer(kind=ik), Allocatable, Dimension(:)   :: obsicu_counties

    Character(Len=:), allocatable, Dimension(:)   :: obsicu_date
    Integer(kind=ik), allocatable, Dimension(:)   :: obsicu_cases
    Character(Len=:), allocatable, Dimension(:)   :: obsicu_sc   

    Integer                                       :: elapsed_days, num_actudays
    Integer(kind=ik)                              :: iter

    Integer                                       :: R0_effects_dim1size

    !=================================================
    ! Optimization target "state"
    !=================================================
    Integer                                       :: num_states
    Integer, Allocatable                          :: state_sc_index(:),pos_se(:),uniq_distid(:)
    Character*15,Dimension(1)                     :: given_sc
    !  Type(target_obsicu),Allocatable,dimension(:)  :: seltarget_icu_by_state 
    Character*10, allocatable, Dimension(:)       :: states_shortcut

    ! Record the observed ICU cases by state
    Real,dimension(:,:), Allocatable              :: seltarget_icu_by_state_icucases
    Character*15,dimension(:,:), Allocatable      :: seltarget_icu_by_state_icudates

    !=================================================
    ! Optimization target "nuts2":
    !=================================================
   
    Integer,allocatable,dimension(:)                   :: unique_nuts2s
    !  Type(target_nuts2_obsicu),dimension(:),Allocatable :: seltarget_icu_by_nuts2 ! Record the observed ICU cases by nuts2
    Integer                                            :: num_nuts2
    Integer(kind=ik)                                   :: loc_sum
    Integer,Allocatable, dimension(:)                  :: pos_nuts2
    Integer(kind=ik), Allocatable, Dimension(:)        :: counties_dist_id
    ! Record the observed ICU cases by nuts2
    Real,dimension(:,:), Allocatable                   :: seltarget_icu_by_nuts2_icucases
    Character*15,dimension(:,:), Allocatable           :: seltarget_icu_by_nuts2_icudates

    !====================================================================================================
    ! Convulution filter and intersection of the observed and evaluation data according to date
    !====================================================================================================
    Real, Allocatable,dimension(:)                :: tmp_filter
    Character*10                                  :: seed_date_mod
    Integer                                       :: halve_coff_len, coff_len ! 7
    Integer, Allocatable                          :: obs_lb(:), obs_ub(:), eval_lb(:), eval_ub(:)
    Integer                                       :: inteval_days, status
    Real                                          :: coefficient

    !=================================================
    ! GA
    !=================================================
    Integer, Allocatable,dimension(:)             :: sel_par ! The selected parents
    Type(opt_res)                                 :: ga_res

    !=================================================
    ! Sorting
    !=================================================
    integer,Allocatable, dimension(:)             :: sorted_index
    integer                                       :: idx
    real                                          :: pre_fitv, curr_fitv

    !=======================================================
    ! Optimization Control Panel / setting the optimization
    !=======================================================
    Logical                                       :: opt_target_icu !.TRUE. Optimiziation target
    Logical                                       :: opt_target_deaths ! .FALSE. Optimization target
    Logical                                       :: opt_filter ! .TRUE.
    character(len=mcl)                            :: use_sug_sol ! "NULL"
  !  character(Len=:),allocatable                  :: opt_region ! "state". Optimized region. Valid values: "state" and "nuts2"
    Integer(kind=ik)                              :: opt_pop_size ! 4 Population size
    Integer(kind=ik)                              :: opt_num_gene ! 2 Maximal number of generations
    Integer                                       :: opt_elitism ! The number of the chosen elitism
    Real                                          :: opt_pcrossover ! 0.8
    Real                                          :: opt_pmutation ! 0.1
    Integer                                       :: tail_pos

    Character*10, dimension(nvals)                :: opt_names
    Real                                          :: opt_lb, opt_ub
    Character*10                                  :: opt_target_region
    Integer, dimension(nvals)                     :: opt_index

    !** Temporarily static initialization of opt ----- TODO: should be written as input data
    call pt_get("#opt_target_icu",opt_target_icu)
    call pt_get("#opt_target_deaths",opt_target_deaths)
    call pt_get("#opt_pop_size",opt_pop_size)
    call pt_get("#opt_max_iter",opt_num_gene)
    !call pt_get("#opt_target_region",opt_region)

    opt_pcrossover    = 0.8
    opt_pmutation     = 0.1
    opt_target_region = "nuts2"

  !  opt_names = (/"SH1","SH2","SH3","SH4","SH5","SH6","HH1","HH2","HH3","HH4","HH5","HH6","NI1","NI2",&
  !    "NI3","NI4","NI5","NI6","HB1","HB2","HB3","HB4","HB5","HB6"/)

    opt_names = (/"def01","def02","def03","def04","def05","de601","de602","de603","de604","de605",&
        "de911","de912","de913","de914","de915","de921","de922",&
        "de923","de924","de925","de931","de932","de933","de934","de935",&
        "de941","de942","de943","de944","de945","de501","de502","de503","de504","de505"/)
    opt_lb = 0.1
    opt_ub = 1.0
   
    opt_elitism = max(1,NINT(opt_pop_size*0.05))

    time_n = maxval(R0change) + 1       
    seed_date_mod  = add_date(seed_date,1)
    iter = iter_e - iter_s + 1

    coff_len = 7
    R0_effects_dim1size = size(R0_effects,dim=1)
    
    select case (trim(opt_target_region))
      !** When the optimized target region is "Nuts2" ------------------
      case("nuts2") 
        counties_dist_id = get_int_table_column(iol%counties,'dist_id')
     
      !  if (.not. allocated(pos_nuts2)) then
      !    Allocate(pos_nuts2(num_nuts2+1))
      !  endif
       
      !    pos_nuts2 = 0
      !    j = 1
      !    pos_nuts2(j) = 1
      !    id_st = region_index(1)
      !    do i = 2,size(region_index)
      !        if(region_index(i) .NE. id_st) then
      !            j = j + 1
      !            pos_nuts2(j) = i
      !            id_st = region_index(i)
      !        endif
      !    enddo
      !   pos_nuts2(num_nuts2+1) = size(region_index) + 1
      
        ! unique_nuts2s: set of affected nuts2; pos_nuts2: see the function "get_unique_nuts2"
        unique_nuts2s = get_unique_nuts2(region_index,pos_nuts2)
        num_nuts2 = size(unique_nuts2s) ! Get the number of affected nuts2

        ! Get the position of the given nuts2 in the full nuts2 list
        ! The results are stored back to array "unique_nuts2s"
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
            !** Calculate the intersection of the evaluated and observed ICU data by date -------------
            call date_intersection(obs_lb(i), eval_lb(i), obs_ub(i), eval_ub(i), seed_date_mod, &
                                      time_n, num_actudays, seltarget_icu_by_nuts2_icudates(1,1), &
                                      seltarget_icu_by_nuts2_icudates(num_actudays,1))
            print *, "nuts2: ", obs_lb(i),eval_lb(i),obs_ub(i),eval_ub(i)
           
          enddo
        endif
           
      !** When the optimized target region is "state" ------------------- 
      case("state")
        uniq_distid = get_unique(region_index) ! Get the affected states
        num_states = size(uniq_distid) ! Get the number of the affected states
        !** Collect the beginning index from which a different state starts ---------------------
        Allocate(pos_se(num_states+1))
        pos_se = 0
        j = 1
        pos_se(j) = 1
        startid = region_index(1)
        do i = 2,size(region_index)
            if(region_index(i) .NE. startid) then
                j = j + 1
                pos_se(j) = i
                startid = region_index(i)
            endif
        enddo
        pos_se(num_states+1) = size(region_index) + 1

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
        !** The fitness relavant to the icu data ----------------------------------
        if (opt_target_icu) then
          !** Precalculate the size of the observed data --------------------------
          elapsed_days = (Date2Unixtime(trim(obsicu_date(size(obsicu_date))))& 
                        - Date2Unixtime(trim(obsicu_date(1))))/86400
          num_actudays = elapsed_days + 1
          Allocate(seltarget_icu_by_state_icucases(num_actudays,num_states))
          Allocate(seltarget_icu_by_state_icudates(num_actudays,num_states))
          seltarget_icu_by_state_icucases = NULVAL
          !** Presume seltarget_icu_by_state is arranged in an ascending order ----------
          !** Preprocess the observed ICU data according to the affected states --------------
          do i = 1,num_states
            given_sc(1) = states_shortcut(uniq_distid(i))
            k = 1
            do j = 1, size(obsicu_sc)
              if (trim(given_sc(1)) .eq. trim(obsicu_sc(j))) then
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
           
            !** Calculate the intersection of the evaluated and observed ICU data by date -------------
            call date_intersection(obs_lb(i), eval_lb(i), obs_ub(i), eval_ub(i), seed_date_mod, &
                                    time_n, num_actudays, seltarget_icu_by_state_icudates(1,1), &
                                    seltarget_icu_by_state_icudates(num_actudays,1))
            if (PT_DEBUG) then
              print *, "Intersection is: "
              print *, "state: ", obs_lb(i),eval_lb(i),obs_ub(i),eval_ub(i)
            endif
          enddo
            
        end if
      case default
        write(*,*)  "Otherwise, the given target is not supported yet"
        call exit(status) 
    end select
      
 
    !** Preparation for preprocessing the parameters subject to optimization -----------------------------------------------
    !** Collect the indexes where the R0 effect data should be replaced with the input that is subject to optimization -----
    opt_index = 0
   
    do j=1,nvals
      do i=1,size(iol%R0_effect%head)
          if (index(trim(opt_names(j)), trim(iol%R0_effect%head(i))) .ne. 0) then
              opt_index(j) = i
              EXIT
          endif
      enddo
      if (opt_index(j) .eq. 0) then
        print *,"Error: the optimized target region is not expected!"
        Call exit(status)
      endif
    enddo
     
    Allocate(fitness(opt_pop_size))
    fitness = 0
    do i = 1,opt_pop_size
      fitness(i) = INV
    enddo
  
    !** Generate beginning population by using a uniform range distribution -----
    Allocate(ini_pop(nvals, opt_pop_size))
    ini_pop = 0
    do j = 1, opt_pop_size
      do i = 1, nvals
        ini_pop(i,j) = random_uniform(opt_lb,opt_ub)
      enddo
    enddo
   
    !** Core of the genetic algorithm -------------------
    do gene_ii=1,opt_num_gene ! Iterate over the generations
      ga_res%opt_acuiter = gene_ii
  
      do pop_ii = 1,opt_pop_size ! Iterate over all the population
        !** Proceed only when the corresponding fitness is not given (aka. INV) ----
        if (fitness(pop_ii) .ne. INV) then
          cycle
        endif
        !** Preprocess the R0 effect data according to the generated population ----
        tail_pos = len(trim(opt_names(1)))
        startid = 1
        do i = 1,nvals
          READ(opt_names(i)(tail_pos:tail_pos), "(I4)") j
          ! Parameter checking
          if (j .gt. R0_effects_dim1size) then
            print *,"Error: the given optimized target is beyond the expected range of R0_effects!"
            Call exit(status)
          endif
          R0_effects(j,opt_index(i)) = ini_pop(startid,pop_ii)
          startid = startid + 1
        enddo
           
        !** Call the COVID19 spatial simuation and store the simulated ICU data to ill_ICU_cases_final -----          
        ill_ICU_cases_final = COVID19_Spatial_Microsimulation_for_Germany(iol,&
            iter_s, iter_e , &
            inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
            less_contagious, R0_force, immune_stop, &
            R0change, R0delay ,R0delay_days, R0delay_type, &
            control_age_sex, seed_date, seed_before, sam_size, R0, &
            R0_effects, region_index, output_dir, export_name, rank_mpi)
        
        diff_sum = 0.0
        diff_len = 0

        !==================================================!
        !====== Calculate fitness function value ==========!
        !==================================================!

        select case (trim(opt_target_region))
          case("nuts2")  
            !** All the data are represented by nuts2 --------
            if (.not.allocated(aggreg_sum_by_county)) then
              Allocate(aggreg_sum_by_county(num_nuts2,time_n,iter))
              Allocate(avg_by_county(num_nuts2,time_n))
            endif
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

            !  if (PT_DEBUG) then
                  do j = 1,time_n 
                    print *,avg_by_county(:,j)
                  enddo
            !  endif
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
                  if(seltarget_icu_by_nuts2_icucases(j,i) .ne. NULVAL) then
                    diff_tmp_power = (seltarget_icu_by_nuts2_icucases(j,i) - avg_by_county(i,k))**2
                    diff_loc_sum = diff_loc_sum + diff_tmp_power
                    diff_loc_len = diff_loc_len + 1
                  endif
                  j = j + 1
                  k = k + 1
                end do
                diff_sum = diff_sum + diff_loc_sum
                diff_len = diff_len + diff_loc_len
              !  allocate(diff(obs_ub(i)-obs_lb(i)+1))
              !  diff = seltarget_icu_by_nuts2_icucases(obs_lb(i):obs_ub(i),i)&
              !   - avg_by_county(i,eval_lb(i):eval_ub(i))
              !  diff = diff*diff
              !  diff_sum = diff_sum + sum(diff)
              !  diff_len = diff_len + obs_ub(i) - obs_lb(i) + 1
              !  deallocate(diff)   
              enddo    
              fitness(pop_ii) = -1*sqrt(diff_sum/diff_len)
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
                  aggreg_sum_by_county(k,i,j) = sum(ill_ICU_cases_final(pos_se(k):pos_se(k+1)-1,i,j)) ! Sum over all the counties within one state
                enddo
              enddo
            enddo
      
            do j = 1,time_n
              do i = 1,num_states
                avg_by_county(i,j) = Real(sum(aggreg_sum_by_county(i,j,:)))/Real(iter) ! Average over the running iterates
              enddo
            enddo
          !  if (PT_DEBUG) then
                do j = 1,time_n 
                    print *,avg_by_county(:,j)
                enddo
          !  endif

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
                if(seltarget_icu_by_state_icucases(j,i) .ne. NULVAL) then
                  diff_tmp_power = (seltarget_icu_by_state_icucases(j,i) - avg_by_county(i,k))**2
                  diff_loc_sum = diff_loc_sum + diff_tmp_power
                  diff_loc_len = diff_loc_len + 1
                endif
                j = j + 1
                k = k + 1
              end do
              diff_sum = diff_sum + diff_loc_sum
              diff_len = diff_len + diff_loc_len

          !    diff = seltarget_icu_by_state_icucases(obs_lb(i):obs_ub(i),i)&
          !     - avg_by_county(i,eval_lb(i):eval_ub(i))
          !    diff = diff*diff
          !    diff_sum = diff_sum + sum(diff)
          !    diff_len = diff_len + obs_ub(i) - obs_lb(i) + 1
          !    deallocate(diff)      
            enddo    
            fitness(pop_ii) = -1*sqrt(diff_sum/diff_len)
          case default
            write(*,*)  "Otherwise, the given target is not supported yet"
            call exit(status)
        end select 
        deallocate(ill_ICU_cases_final)
      enddo


      if ( .Not.Allocated(tmp_pop) ) then                
        Allocate(tmp_pop(nvals, opt_pop_size))
        Allocate(sorted_pop(nvals, opt_elitism))

        Allocate(tmp_fitness(opt_pop_size))
        Allocate(sorted_fitness(opt_elitism))

        Allocate(sel_par(opt_pop_size))


        Allocate(sorted_index(opt_pop_size))

        Allocate(summary_fitness%maximum(opt_num_gene))
        Allocate(summary_fitness%mean(opt_num_gene))
        Allocate(summary_fitness%upper_hinge(opt_num_gene))
        Allocate(summary_fitness%median(opt_num_gene))
        Allocate(summary_fitness%lower_hinge(opt_num_gene))
        Allocate(summary_fitness%minimum(opt_num_gene))

        print *,"allocated"
      endif

      print *,"intermediate fitness output:"
      print *,fitness
      sorted_index = 0
    
      call SORTRX_REAL(opt_pop_size,fitness,sorted_index) ! Sorted_index indicates a sorted fitness in an ascending order

      !** Add summary statistics --------------------------------------------------------
      call sixnum_summary(fitness, sorted_index, summary_fitness%maximum(gene_ii), & 
                                   summary_fitness%minimum(gene_ii), &
                                   summary_fitness%mean(gene_ii), &
                                   summary_fitness%median(gene_ii), &
                                   summary_fitness%upper_hinge(gene_ii), &
                                   summary_fitness%lower_hinge(gene_ii))
      print *,"generation: ",gene_ii
      print *,"  max                ","mean               ",&
                "upper_hinge         ","median       ",&
                "lower_hinge           ","min"
      print *,summary_fitness%maximum(gene_ii)," ",summary_fitness%mean(gene_ii),&
                        " ",summary_fitness%upper_hinge(gene_ii),&
                        " ", summary_fitness%median(gene_ii), " ", &
                        summary_fitness%lower_hinge(gene_ii), " ", summary_fitness%minimum(gene_ii) 

        
      !! TODO: To record the best fitness value and solution per generation (write them to a output file)

      ! ========================================================== !
      ! ============ Check the stop criteria ===================== !
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
      !** 3. all the values in the fitness are equal ---------------------------------------------
      if ((gene_ii .eq. opt_num_gene) .or. (run_sum .ge. opt_num_gene)&
          .or. (fitness(max_pos) .eq. fitness(min_pos))) then
        exit
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
          print *,"too many duplicates lead to insufficient elitism"
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
      print *,sel_par
      do i = 1,opt_pop_size
        tmp_pop(:,i) = ini_pop(:,sel_par(i))
        tmp_fitness(i) = fitness(sel_par(i))
      enddo   
      ini_pop = tmp_pop
      fitness = tmp_fitness

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
    enddo

    ! Store the best solution to "ga_res"
    ga_res%opt_bestsol      = ini_pop(:,max_pos)
    ga_res%opt_fitnessvalue = fitness(max_pos)

    call print_ga_summary(opt_pop_size, opt_num_gene, opt_elitism, opt_pcrossover, opt_pmutation, &
              opt_target_region, opt_names, ga_res)
      
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
        deallocate(pos_se)
      case default
        write(*,*)  "Otherwise, the given target is not supported yet"
        call exit(status) 
    end select

    deallocate(aggreg_sum_by_county)
    deallocate(avg_by_county)
    
    deallocate(obs_ub)
    deallocate(obs_lb)
    deallocate(eval_lb)
    deallocate(eval_ub)
    deallocate(ini_pop)
    deallocate(fitness)
  
    deallocate(tmp_pop)
    deallocate(tmp_fitness)
    deallocate(sorted_pop)
    deallocate(sorted_fitness)
    deallocate(sorted_index)
    deallocate(sel_par)

    deallocate(summary_fitness%maximum)
    deallocate(summary_fitness%mean)
    deallocate(summary_fitness%upper_hinge)
    deallocate(summary_fitness%median)
    deallocate(summary_fitness%lower_hinge)
    deallocate(summary_fitness%minimum)
  end subroutine ga

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
  !> Subroutine that get a subset of the input array_in (with repeated nuts2) storing the unique nuts2.
  !> 
  !> The returned arrays are array_pos and array_out.
  !> array_out: store the unique nuts2; array_pos: store the position from which a different nuts2 starts in array_in
  function get_unique_nuts2(array_in, array_pos) result(array_out)
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
    !detemine dimension of array_out
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
  
  end function get_unique_nuts2

  !! -------------------------------------------------------------------------------------------
  !> Subroutine that Calculate the intersection of the evaluated and observed ICU data by date.
  !>
  !> The intersection is represented by the returned intervals [obs_lb, obs_ub] and [eval_lb, eval_ub].
  subroutine date_intersection(obs_lb, eval_lb, obs_ub, eval_ub, seed_date_mod, time_n, length, s_date, e_date)
    Integer,intent(inout)                                  :: obs_lb, obs_ub, eval_lb, eval_ub
    Integer,intent(in)                                     :: length
    Character(len=*)                                       :: s_date, e_date
    Integer                                                :: inteval_days, status, time_n
    Character*10                                           :: seed_date_mod

    !** Here we assume that the date is monotonically increased by days -----------------------
    obs_lb = 1
    eval_lb = 1

    print *,s_date 
    print *,e_date

    !TODO: Move the calculation of inteval_days to support module
    inteval_days = Date2Unixtime(trim(s_date))&
                        - Date2Unixtime(seed_date_mod)
    inteval_days = inteval_days/86400

    if (inteval_days .LT. 0) then
      obs_lb = obs_lb - inteval_days
    elseif (inteval_days .GT. 0) then
      eval_lb = eval_lb + inteval_days
    endif

    obs_ub = length
    eval_ub = time_n

    inteval_days = (Date2Unixtime(trim(e_date))& 
        - (Date2Unixtime(seed_date_mod)+86400*(time_n-1)))/86400
    if (inteval_days .LT. 0) then
      eval_ub = eval_ub + inteval_days
    elseif (inteval_days .GT. 0) then
      obs_ub = obs_ub - inteval_days
    endif

    print *,eval_lb," ",eval_ub, " ",obs_lb, " ",obs_ub

    !** Throw an error when the intersection is null ---------------------------
    if ((eval_lb .GT. eval_ub) .AND. (obs_lb .GT. obs_ub) .AND.&
     ((eval_ub - eval_lb) .NE. (obs_ub - obs_lb))) then
      print *, "Error: intersection is null!"
      print *, "eval_lb: ",eval_lb, " eval_ub: ",eval_ub," obs_lb: ",obs_lb," obs_ub: ",obs_ub
      Call exit(status)    
    endif
  end subroutine date_intersection

  !! ------------------------------------------------------------------------------
  !> Subroutine that generates random population using a uniform range distribution.
  !>
  !> The random number is between lb and ub.
  subroutine population(ini_pop, lb, ub, popsize, nvals)
    Real(kind=rk),dimension(:,:),Allocatable,intent(inout) :: ini_pop
    Real,intent(in)                               :: lb 
    Real,intent(in)                               :: ub 
    Integer,intent(in)                                     :: popsize 
    Integer,intent(in)                                     :: nvals

    Integer                                                :: j,i

    ini_pop = 0
    do j = 1, popsize
      do i = 1, nvals
      ini_pop(i,j) = random_uniform(lb,ub)
      enddo
    enddo
  end subroutine population

  !! ------------------------------------------------------------------------------
  !> Subroutine that selects the fittest parents
  !>
  !> Return the result indexes to sel_par.
  subroutine selection(fitness, sel_par)
    Real, Allocatable, Dimension(:),intent(inout)         :: fitness
    Integer, Dimension(size(fitness)),intent(inout)       :: sel_par

    Real, Dimension(size(fitness))                        :: fscaled, prob
    Real                                                  :: fmin, fave, fmax
    Real                                                  :: delta, a, b
    Integer                                               :: sfactor ! 2

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
    if (PT_DEBUG) then
      print *,"The selected chromosomes: "
      print *,sel_par
    endif
      
  end subroutine selection

  !! ------------------------------------------------------------------------------
  !> Subroutine that mates the selected pairs and breed two children for each
  !>
  !> The results are stored back to ini_pop and fitness.
  subroutine crossover(ini_pop, fitness, pcrossover) 
    Real(kind=rk),dimension(:,:),Allocatable,intent(inout)    :: ini_pop
    Real,dimension(:),Allocatable,intent(inout)               :: fitness
    Real,intent(in)                                           :: pcrossover

    Integer,dimension(:),Allocatable                          :: sample_input, sample_parents
    Integer                                                   :: start_id, end_id, nmating,i,j,len
    Real(kind=rk)                                             :: ran
    Real(kind=rk),dimension(:),Allocatable                    :: lchild, rchild

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
    Real(kind=rk),dimension(:,:),Allocatable,intent(inout)    :: ini_pop
    Real,dimension(:),Allocatable,intent(inout)               :: fitness
    Real,intent(in)                                           :: pmutation
    Real, intent(in)        :: lb, ub

    Integer,dimension(size(ini_pop,dim=1))                    :: sample_input
    Integer,dimension(1)                                      :: sample_parents
    Integer                                                   :: i
    Real                                                      :: ran
    Real(kind=rk)                                             :: Res

    do i = 1,size(sample_input)
      sample_input(i) = i
    enddo
    do i = 1,size(fitness)
      call random_number(ran)
      if (pmutation > ran) then
        sample_parents = sample_i4(sample_input,size(sample_parents)) ! Randomly generate the mutated gene/parameter
        Res = random_uniform(lb,ub)
        ini_pop(sample_parents(1),i) = Res  
        print *,"Mutation: chromosome ",i,",gene: ",sample_parents(1),",res: ", Res
        fitness(i) = INV ! Invalidate the mutated chromosome/population
      endif
    enddo
  end subroutine mutation

  !! ------------------------------------------------------------------------------
  !> Subroutine that replaces the unfittest individules with the chosen elitim
  !>
  !> The results are stored back to ini_pop and fitness.
  subroutine elitism(fitness, ini_pop, sorted_pop, sorted_fitness)
    Real(kind=rk),dimension(:,:),Allocatable,intent(inout)    :: ini_pop
    Real,dimension(:),Allocatable,intent(inout)               :: fitness
    Real(kind=rk),dimension(:,:),Allocatable,intent(in)       :: sorted_pop                                   
    Real,dimension(:),Allocatable,intent(in)                  :: sorted_fitness
    Integer,dimension(size(fitness))                          :: sorted_index
    Integer                                                   :: num

    call SORTRX_REAL(size(fitness),fitness,sorted_index) ! Fitness in an ascending order
    num = size(sorted_fitness)
    !** Replace the unfittest individules with the elitism chosen in the first place ---
    !** The elitism is represented as sorted_fitness and sorted_pop
    fitness(sorted_index(1:num)) = sorted_fitness
    ini_pop(:,sorted_index(1:num)) = sorted_pop
    
    if(PT_DEBUG) then
      print *,fitness 
      do num = 1,size(ini_pop,1)
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
  subroutine sixnum_summary(fitness, sorted_index, max, min, mean, median, upper_hinge, lower_hinge)
    Real,Allocatable,dimension(:),intent(in)    :: fitness
    Integer,Allocatable,dimension(:),intent(in) :: sorted_index

    Real, intent(out)                           :: max, min, mean, median, upper_hinge, lower_hinge
    Integer                                     :: pop_size, half, i, j

    pop_size = size(fitness)
    max = fitness(sorted_index(pop_size))
    min = fitness(sorted_index(1))
    mean = Sum(fitness)/pop_size

    median = get_median(fitness, sorted_index, pop_size)
    half = pop_size/2.0
    upper_hinge = get_median(fitness, sorted_index((half+1):pop_size), pop_size-half) ! Get the median of the second half
    lower_hinge = get_median(fitness, sorted_index(1:half), half) ! Get the median of the first half
  end subroutine sixnum_summary
 
  !! ------------------------------------------------------------------------------
  !> Function that returns the median value of an ordered array
  !>
  function get_median(dataset, sorted_index, len) Result(med)
    Real,Allocatable,dimension(:),intent(in)     :: dataset
    Integer,intent(in)                           :: len
    Integer,dimension(len),intent(in)            :: sorted_index
    Real                                         :: med
    Integer                                      :: length, half

    length = size(sorted_index)
    half = length/2

    if (mod(length,2) .eq. 0) then
      med = (dataset(sorted_index(half))+dataset(sorted_index(half+1)))/2.0
    else
      med = dataset(sorted_index(half+1))
    endif
  end function get_median


  subroutine print_ga_summary(opt_pop_size, opt_num_gene, opt_elitism, opt_pcrossover, opt_pmutation, &
                  opt_target_region, opt_names, ga_res)  
    character(len=*), intent(in)                              :: opt_target_region ! Optimized region. Valid values: "state" and "nuts2"
    character(len=*), Dimension(nvals), intent(in)            :: opt_names ! Names of parameters (by state or Nuts-2) to optimize
   
    Integer, intent(in)                                       :: opt_pop_size ! 4 Population size
    Integer, intent(in)                                       :: opt_num_gene ! 2 Maximal number of generations
    Integer, intent(in)                                       :: opt_elitism ! The number of the chosen elitism
    Real, intent(in)                                          :: opt_pcrossover ! 0.8
    Real, intent(in)                                          :: opt_pmutation ! 0.1
    type(opt_res),intent(in)                                  :: ga_res
    Integer                                                   :: i

    write(*,'(A)')"!----------------------- Genetic Algorithm ----------------------------!"
    write(*,'(A)')"!------ GA settings ------!"
    print *,"Type: ","real-valued"
    print *,"Population size: ",opt_pop_size
    print *,"Number of generations: ",opt_num_gene ! maxiter
    print *,"Elitism: ",opt_elitism
    print *,"Crossover probability: ",opt_pcrossover
    print *,"Mutation probability: ",opt_pmutation
    print *,"Optimization target: ",trim(opt_target_region)
    write(*,'(A)')"!------ GA results -------!"
    print *,"Iterations: ",ga_res%opt_acuiter
    print *,"Fitness function value: ",ga_res%opt_fitnessvalue
    print *,"Solution:"
    print *,ga_res%opt_bestsol
    print *,"Affected nuts2:"
    do i = 1, size(opt_names)
      print *,opt_names(i)
    enddo
  end subroutine print_ga_summary
End Module genetic_algorithm
