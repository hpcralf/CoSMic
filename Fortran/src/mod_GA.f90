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

    use precision
    use urandom
    use global_types
    use precision
    use cosmic_io
    use kernel
    use qsort_c_module
    use quicksort_nr 

    Use data_preprocessing
    Use cosmic_io

    implicit none

    !** Constants ---------------------------------------

    Real, parameter     :: INV = 99.0 ! Corresponds to NA for fitness array 
    Real, parameter     :: EPS = EPSILON(1.0) ! Tolerance precision, a magnitude of e-07

    !** Extend the fivenum summary statistic with a mean value -------
    type ga_summary
        Real                                          :: max
        Real                                          :: mean 
        Real                                          :: upper_hinge
        Real                                          :: median
        Real                                          :: lower_hinge
        Real                                          :: min
    end type ga_summary

    !** Optimization Control Panel / setting the optimization --------
    type opt_parameters			
        Logical                                       :: opt_target_icu
        Logical                                       :: opt_target_deaths ! Optimization targets
        Logical                                       :: opt_filter  = .TRUE.
        character(len=mcl)							  :: use_sug_sol = "NULL"
        character(len=8)							  :: opt_target_region = "nuts2" ! Optimized region. Valid values: "state" and "nuts2"
        character(len=8), Dimension(:), Allocatable   :: opt_names ! Names of parameters (by Nuts-2) to optimize
        Integer,Dimension(:), Allocatable             :: opt_names_dur ! The duration of certain parameter is calculated by weeks
        Real                        				  :: opt_lb = 0.1 ! Lower bounds of optimized parameters 
        Real                                          :: opt_ub = 1.0 ! Upper bounds of optimized parameters
        Integer 									  :: opt_pop_size = 8 ! Population size
        Integer    								      :: opt_num_gene = 4 ! Maximal number of generations
        Integer                                       :: opt_elitism ! The number of the chosen elitism
        Real                                          :: opt_pcrossover = 0.8
        Real                                          :: opt_pmutation = 0.1
    end type opt_parameters
 
    !** GA results/best solution -------------------------------------
    type opt_res 
        Real(kind=rk),dimension(:), Allocatable       :: opt_bestsol ! Solution
        Real                                          :: opt_fitnessvalue ! Best fitness
        Integer                                       :: opt_acuiter ! The occured number of generations
    end type opt_res 

    !** Observed icu cases by state ----------------------------------
    type target_obsicu
        Real,dimension(:), Allocatable                :: icucases
        Character*10,dimension(:), Allocatable        :: icudates
    end type target_obsicu

    !** Observed icu cases by nuts2 ----------------------------------
    type target_nuts2_obsicu
        Real,dimension(:), Allocatable                :: nuts2_icucases
        Character*10,dimension(:), Allocatable        :: nuts2_icudates
    end type target_nuts2_obsicu

Contains

    subroutine ga(&
       iol, counties_index, &
       iter , &
       inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
       less_contagious, R0_force, immune_stop, &
       R0change, R0delay ,R0delay_days, R0delay_type, &
       control_age_sex, seed_date, seed_before, sam_size, R0, &
       R0_effects) !type opt_parameters

        !===========================================================================
        ! Declaration
        !===========================================================================

        Type(iols)                                               :: iol
        Integer(kind=ik)             , Dimension(:),intent(in)   :: counties_index
        Integer(kind=ik)                           ,intent(in)   :: iter
        Integer(kind=ik)                           ,intent(in)   :: inf_dur
        Integer(kind=ik)                           ,intent(in)   :: cont_dur
        Integer(kind=ik)                           ,intent(in)   :: ill_dur
        Integer(kind=ik)                           ,intent(in)   :: icu_dur
        Integer(kind=ik), Allocatable, Dimension(:),intent(in)   :: icu_per_day
        Real(kind=rk)                              ,intent(in)   :: less_contagious
        Real(kind=rk)                              ,intent(in)   :: R0_force
        Logical                                    ,intent(in)   :: immune_stop
        Integer(kind=ik), Allocatable, Dimension(:,:) ,intent(in):: R0change
        Logical                                       ,intent(in):: R0delay
        Integer(kind=ik)                             ,intent(in) :: R0delay_days
        Character(len=:), Allocatable                 ,intent(in):: R0delay_type
        character(len=:), Allocatable                 ,intent(in):: control_age_sex
        character(len=:), Allocatable                 ,intent(in):: seed_date
        Integer(kind=ik)                              ,intent(in):: seed_before    
        Integer(kind=ik)                              ,intent(in):: sam_size
        Real(kind=rk)                                 ,intent(in):: R0
        Real(kind=rk),    Allocatable, Dimension(:,:) ,intent(in):: R0_effects

        !==========================================================================
        ! General purpose for preprocessing and fitness calculation
        !==========================================================================
        Integer                                       :: nvals ! The total number of parameters to be optimized
        Type(opt_parameters)                          :: opt
        
        Real(kind=rk),Dimension(:,:), Allocatable     :: ini_pop, tmp_pop, sorted_pop ! ini_pop: store the population, which are the input to be optimized. 

        Real, Allocatable, Dimension(:)               :: fitness, tmp_fitness ! fitness: store the calculated fitness
        Real, Allocatable, Dimension(:)               :: sorted_fitness
        Type(ga_summary), Allocatable, Dimension(:)   :: summary_fitness

        Integer,dimension(:,:,:), Allocatable         :: ill_ICU_cases_final ! Store the evaluated ICU cases
        Integer, dimension(:), Allocatable            :: opt_index
        Integer                                       :: i, j, k, pop_ii, gene_ii
        Integer                                       :: startid, endid, max_pos, min_pos
        Real                                          :: tmp_max, run_sum
        Integer                                       :: time_n
        
        Real,dimension(:), Allocatable                :: diff
        Real                                          :: diff_sum = 0.0
        integer                                       :: diff_len = 0
        Integer,dimension(:,:,:), Allocatable         :: aggreg_sum_by_county
        Real,dimension(:,:), Allocatable              :: avg_by_county ! Record the simulated ICU cases by state/nuts2
        

        !=================================================
        ! Target region "state"
        !=================================================
        Integer                                       :: num_states
        Integer,dimension(:), Allocatable             :: state_sc_index,pos_se,uniq_distid
        Character*15,Dimension(1)                     :: given_sc
        Type(target_obsicu),dimension(:),Allocatable  :: seltarget_icu_by_state ! Record the observed ICU cases by state

        !=================================================
        ! Target region "nuts2"
        !=================================================
        character*10,allocatable,dimension(:)              :: nuts2_set, unique_nuts2s
        Type(target_nuts2_obsicu),dimension(:),Allocatable :: seltarget_icu_by_nuts2 ! Record the observed ICU cases by nuts2
        Integer                                            :: num_nuts2, elapsed_days, num_actudays, loc_sum
        Integer,dimension(:),Allocatable                   :: distid_index, pos_nuts2

        !====================================================================================================
        ! Convulution filter and intersection of the observed and evaluation data according to date
        !====================================================================================================
        Real,dimension(:), Allocatable                :: tmp_filter
        Character*10                                  :: seed_date_mod
        Integer                                       :: halve_coff_len, coff_len=7
        Integer,dimension(:), Allocatable             :: obs_lb, obs_ub, eval_lb, eval_ub
        Integer                                       :: inteval_days, status
        Real                                          :: coefficient

        !=================================================
        ! GA
        !=================================================
        Integer,dimension(:), Allocatable             :: sel_par ! The selected parents
        Type(opt_res)                                 :: ga_res

        !=================================================
        ! Sorting
        !=================================================
        integer,dimension(:),Allocatable              :: sorted_index
        integer                                       :: idx
        real                                          :: pre_fitv, curr_fitv


        ! TODO: parameters checking
        ! ====================================================================
        ! ====================================================================


        !** Temporarily static initialization of opt ----- TODO: should be written as input data
        opt%opt_target_icu    = .TRUE.
        opt%opt_target_deaths = .FALSE.
        opt%opt_names = (/ "def0", "de60", "de91", "de92", "de93", "de94", "de50" /)
        opt%opt_names_dur = (/ 6, 6, 6, 6, 6, 6, 6 /)
        nvals = sum(opt%opt_names_dur)
        
        opt%opt_elitism = max(1,NINT(opt%opt_pop_size*0.05))
    
        Allocate(ini_pop(nvals, opt%opt_pop_size))
        Allocate(fitness(opt%opt_pop_size))
        Allocate(tmp_pop(nvals, opt%opt_pop_size))
        Allocate(tmp_fitness(opt%opt_pop_size))
        Allocate(opt_index(size(opt%opt_names)))

        Allocate(sorted_pop(nvals, opt%opt_elitism))
        Allocate(sorted_fitness(opt%opt_elitism))
        Allocate(sorted_index(opt%opt_pop_size))
        Allocate(summary_fitness(opt%opt_num_gene))

        time_n = maxval(R0change) + 1       
       
        Allocate(sel_par(opt%opt_pop_size))

        seed_date_mod        = add_date(seed_date,1)

        select case (trim(opt%opt_target_region))
            !** When the optimized target region is "Nuts2" ------------------
            case("nuts2")
                distid_index = get_index_mul_integer(iol%counties_dist_id,counties_index)
                Allocate(nuts2_set(size(distid_index)))
                nuts2_set = iol%counties_nuts2(distid_index)
                ! unique_nuts2s: set of affected nuts2; pos_nuts2: see the function "get_unique_nuts2"
                unique_nuts2s = get_unique_nuts2(nuts2_set,pos_nuts2)
                num_nuts2 = size(unique_nuts2s) ! Get the number of affected nuts2
                deallocate(distid_index)
                deallocate(nuts2_set)
                deallocate(unique_nuts2s)

                Allocate(obs_lb(num_nuts2))
                Allocate(obs_ub(num_nuts2))
                Allocate(eval_lb(num_nuts2))
                Allocate(eval_ub(num_nuts2))
                !** TODO: The fitness relavant to the death data ------
                !** The fitness relavant to the icu data --------------
                if (opt%opt_target_icu) then
                    Allocate(seltarget_icu_by_nuts2(num_nuts2))

                    !** Predict the size of the observed data, which won't exceed the time_n -----------------------------------------------
                    !** In this case, the ultimate filter results are probably slightly different from the original one (refer to R) -------
                    elapsed_days = (Date2Unixtime(iol%obsicu_nuts2_date(size(iol%obsicu_nuts2_date)))& 
                                - (Date2Unixtime(iol%obsicu_nuts2_date(1))+86400*(time_n-1)))/86400
                    if (elapsed_days .lt. time_n) then
                        num_actudays = elapsed_days+1
                    else
                        num_actudays = time_n
                    endif

                    do i = 1, num_nuts2
                        allocate(seltarget_icu_by_nuts2(i)%nuts2_icucases(num_actudays))
                        allocate(seltarget_icu_by_nuts2(i)%nuts2_icudates(num_actudays))
                    enddo
                   
                    do i = 1, num_nuts2
                        j = 1
                        loc_sum = 0
                        seltarget_icu_by_nuts2(i)%nuts2_icudates(j) = iol%obsicu_nuts2_date(j)
                        !** Calculate the icu cases of certain nuts2 during the limited time period (num_actudays) --------------
                        do k = 1, size(iol%obsicu_nuts2_date)
                            !** When we reach a different date, ------------------------------------------------------
                            !** the loc_sum for the previous date will be assigned to the current nuts2_icucases -----
                            if (trim(iol%obsicu_nuts2_date(k)) .ne. trim(seltarget_icu_by_nuts2(i)%nuts2_icudates(j))) then
                                seltarget_icu_by_nuts2(i)%nuts2_icucases(j) = Real(loc_sum)
                                loc_sum = 0
                                j = j + 1
                                if (j .gt. num_actudays) then 
                                    exit
                                else
                                    seltarget_icu_by_nuts2(i)%nuts2_icudates(j) = iol%obsicu_nuts2_date(k)
                                endif
                            endif
                            !** Accumulate only when the given nuts2 is hit ---------------------------
                            if (BinarySearch(iol%obsicu_counties(k),counties_index(pos_nuts2(i):pos_nuts2(i+1)-1)) .ne. 0) then
                                loc_sum = loc_sum + iol%obsicu_nuts2_cases(k)
                            endif        
                        enddo
                        !** Apply convolution filter ---------------
                        if (opt%opt_filter) then
                            call conv_filter(seltarget_icu_by_nuts2(i)%nuts2_icucases, coff_len, num_actudays)
                        endif
                        !** Calculate the intersection of the evaluated and observed ICU data by date -------------
                        call date_intersection(obs_lb(i), eval_lb(i), obs_ub(i), eval_ub(i), seed_date_mod,&
                                                                             time_n, seltarget_icu_by_nuts2(i)%nuts2_icudates)
                    !    print *, "nuts2: ", obs_lb(i),eval_lb(i),obs_ub(i),eval_ub(i)
                    enddo
                end if
            !** When the optimized target region is "state" ------------------- 
            case("state")
                !** Collect the beginning index from which a different state starts ---------------------
                uniq_distid = get_unique(counties_index/1000) ! Get the affected states
                num_states = size(uniq_distid) ! Get the number of the affected states
                Allocate(pos_se(num_states+1))
                j = 1
                pos_se(j) = 1
                startid = counties_index(1)
                do i = 2,size(counties_index)
                    if(counties_index(i)/1000 .NE. startid/1000) then
                        j = j + 1
                        pos_se(j) = i
                        startid = counties_index(i)
                    endif
                enddo
                pos_se(num_states+1) = size(counties_index) + 1

                Allocate(obs_lb(num_states))
                Allocate(obs_ub(num_states))
                Allocate(eval_lb(num_states))
                Allocate(eval_ub(num_states))
                !** TODO: The fitness relavant to the death data ------
                !** The fitness relavant to the icu data --------------
                if (opt%opt_target_icu) then
                    Allocate(seltarget_icu_by_state(num_states))
                    do i = 1,num_states
                        given_sc(1) = iol%states_shortcut(uniq_distid(i))
                        state_sc_index = get_index_mul_char(iol%obsicu_states_shortcut,given_sc)

                        !** Presume seltarget_icu_by_state is arranged in an ascending order ----------
                        !** Preprocess the observed ICU data according to the affected states --------------
                        seltarget_icu_by_state(i)%icucases = Real(iol%obsicu_cases(state_sc_index))
                        seltarget_icu_by_state(i)%icudates = iol%obsicu_date(state_sc_index)

                        !** Apply convolution filter ------
                        if (opt%opt_filter) then
                            call conv_filter(seltarget_icu_by_state(i)%icucases, coff_len, size(state_sc_index))                       
                        endif
                        deallocate(state_sc_index)

                        !** Calculate the intersection of the evaluated and observed ICU data by date -------------
                        call date_intersection(obs_lb(i), eval_lb(i), obs_ub(i), eval_ub(i), seed_date_mod,&
                                                                             time_n, seltarget_icu_by_state(i)%icudates)
                    !    print *, "state: ", obs_lb(i),eval_lb(i),obs_ub(i),eval_ub(i)
                    enddo
                end if
            case default
                write(*,*)  "Otherwise, the given target is not supported yet" 
        end select

        !** Preparation for preprocessing the parameters subject to optimization -----------------------------------------------
        !** Collect the indexes where the R0 effect data should be replaced with the input that is subject to optimization -----
        do j=1,size(opt%opt_names)
            do i=1,size(iol%R0_effect%head)
                if (index(trim(iol%R0_effect%head(i)), trim(opt%opt_names(j))) .NE. 0) then
                    opt_index(j) = i
                    EXIT
                endif
            enddo
        enddo

        !** Generate beginning population --------------------
        call population(ini_pop, opt%opt_lb, opt%opt_ub, opt%opt_pop_size, nvals)
        fitness(:) = INV ! initialize fitness with "NA"
                
        !** Core of the genetic algorithm -------------------
        do gene_ii=1,opt%opt_num_gene ! Iterate over the generations
            ga_res%opt_acuiter = gene_ii
            do pop_ii = 1,opt%opt_pop_size ! Iterate over all the population
                !** Proceed only when the corresponding fitness is not given (aka. INV) ----
                if (fitness(pop_ii) .ne. INV) then
                    cycle
                endif
                !** Preprocess the R0 effect data according to the generated population ----
                startid = 1
                do i=1,size(opt%opt_names)
                    endid = startid + opt%opt_names_dur(i) - 1
                    iol%R0_effect%data(1:opt%opt_names_dur(i),opt_index(i)) = ini_pop(startid:endid,pop_ii)
                    startid = endid + 1
                enddo
                !** Call the COVID19 spatial simuation and store the simulated ICU data to ill_ICU_cases_final -----
                Call COVID19_Spatial_Microsimulation_for_Germany(iol,counties_index, &
                    iter , &
                    inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
                    less_contagious, R0_force, immune_stop, &
                    R0change, R0delay ,R0delay_days, R0delay_type, &
                    control_age_sex, seed_date, seed_before, sam_size, R0, &
                    iol%R0_effect%data,ill_ICU_cases_final)

                diff_sum = 0.0
                diff_len = 0

                !==================================================!
                !====== Calculate fitness function value ==========!
                !==================================================!

                select case (trim(opt%opt_target_region))
                    case("nuts2")
                        !** All the data are represented by nuts2 --------
                        Allocate(aggreg_sum_by_county(num_nuts2,time_n,iter))
                        Allocate(avg_by_county(num_nuts2,time_n))

                        !** TODO: The fitness relavant to the death data ----------
                        !** The fitness relavant to the icu data ------------------
                        if (opt%opt_target_icu) then
                            do j = 1,iter
                                do i = 1,time_n 
                                    do k = 1,num_nuts2
                                        aggreg_sum_by_county(k,i,j) = sum(ill_ICU_cases_final(pos_nuts2(k):pos_nuts2(k+1)-1,i,j)) ! Sum over all the counties within one nuts2
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
                                if (opt%opt_filter) then
                                    call conv_filter(avg_by_county(i,:), coff_len, time_n)                   
                                endif
                                allocate(diff(obs_ub(i)-obs_lb(i)+1))
                                diff = seltarget_icu_by_nuts2(i)%nuts2_icucases(obs_lb(i):obs_ub(i))&
                                 - avg_by_county(i,eval_lb(i):eval_ub(i))
                                diff = diff*diff
                                diff_sum = diff_sum + sum(diff)
                                diff_len = diff_len + obs_ub(i) - obs_lb(i) + 1
                                deallocate(diff)   
                            enddo    
                            fitness(pop_ii) = -1*sqrt(diff_sum/diff_len)
                        end if

                    case("state")
                        !** All the data are represented by state --------
                        Allocate(aggreg_sum_by_county(num_states,time_n,iter))
                        Allocate(avg_by_county(num_states,time_n))

                        !** TODO: The fitness relavant to the death data ----------
                        !** The fitness relavant to the icu data ------------------
                        if (opt%opt_target_icu) then
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
                                    avg_by_county(i,j) = Real(sum(aggreg_sum_by_county(i,j,:)))/Real(iter) ! Average over the runing iterates
                                enddo
                            enddo
                            if (PT_DEBUG) then
                                do j = 1,time_n 
                                    print *,avg_by_county(:,j)
                                enddo
                            endif
                            
                            !** Calculate the difference between the observed and evaluated/simulated ICU data -----------
                            do i = 1,num_states                  
                                if (opt%opt_filter) then 
                                    call conv_filter(avg_by_county(i,:), coff_len, time_n)
                                endif
                                diff = seltarget_icu_by_state(i)%icucases(obs_lb(i):obs_ub(i))&
                                 - avg_by_county(i,eval_lb(i):eval_ub(i))
                                diff = diff*diff
                                diff_sum = diff_sum + sum(diff)
                                diff_len = diff_len + obs_ub(i) - obs_lb(i) + 1
                                deallocate(diff)      
                            enddo    
                            fitness(pop_ii) = -1*sqrt(diff_sum/diff_len)
                        end if
                    case default
                        write(*,*)  "Otherwise, the given target is not supported yet"
                end select    
                if(opt%opt_target_icu) then                          
                    deallocate(aggreg_sum_by_county)
                    deallocate(avg_by_county)
                end if 
            !    deallocate(ill_ICU_cases_final)
            enddo
        
            print *,"intermediate fitness output:"
            print *,fitness
          
            call SORTRX_REAL(opt%opt_pop_size,fitness,sorted_index) ! Sorted_index indicates a sorted fitness in an ascending order
       
            !** Add summary statistics --------------------------------------------------------
            summary_fitness(gene_ii) = sixnum_summary(fitness, sorted_index)
            write(*,'(A12,I2)')"generation: ",gene_ii
            write(*,'(A13,A15,A20,A15,A20,A15)')"max","mean","upper_hinge","median","lower_hinge","min"
            print *,summary_fitness(gene_ii)

            !! TODO: To record the best fitness value and solution per generation (write them to a output file)

            ! ========================================================== !
            ! ============ Check the stop criteria ===================== !
            run_sum = 0
            max_pos = sorted_index(opt%opt_pop_size)
            min_pos = sorted_index(1)
            
            if (gene_ii .gt. 1) then
                tmp_max = maxval(summary_fitness(1:gene_ii)%max)
                do i = 1,gene_ii 
                    if (summary_fitness(i)%max .ge. (tmp_max - EPS)) then
                        run_sum = run_sum + 1
                    endif
                enddo
            endif
            !** The genetic algorithm stops when -------------------------------------------------------
            !** 1. reaches the generation limit --------------------------------------------------------
            !** 2. the maximal fitnesses over the passed generations are equal with tolerable error ----
            !** 3. all the values in the fitness are equal ---------------------------------------------
            if ((gene_ii .eq. opt%opt_num_gene) .or. (run_sum .ge. opt%opt_num_gene)&
                .or. (fitness(max_pos) .eq. fitness(min_pos))) then
                exit
            endif

            i = 1
            j = opt%opt_pop_size
            pre_fitv = INV
            !** Store the best fitness and population of opt_elitism to the array sorted_fitness and sorted_pop --
            !** Note that the repeated fitness and population is excluded ----------------------------------------
            do while(i .le. opt%opt_elitism)
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
            do i = 1,opt%opt_pop_size
                tmp_pop(:,i) = ini_pop(:,sel_par(i))
                tmp_fitness(i) = fitness(sel_par(i))
            enddo   
            ini_pop = tmp_pop
            fitness = tmp_fitness

            !** The selected parents(pair) do the mating and bread two children to replace them ---------
            call crossover(ini_pop, fitness, opt%opt_pcrossover)
            
            if(PT_DEBUG) then
                do i = 1,nvals
                    print *,ini_pop(i,:)
                enddo
                print *,fitness
            endif

            !** To randomly mutate a parameter/gene ------------------------------------------------------
            if(opt%opt_pmutation .ne. 0.0) then
                call mutation(ini_pop, fitness, opt%opt_pmutation, opt%opt_lb, opt%opt_ub)
            endif
            !** Replace the unfittest populations with the fittest ones, which are stored in sorted_pop ---
            call elitism(fitness, ini_pop, sorted_pop, sorted_fitness)
        enddo

        ga_res%opt_bestsol = ini_pop(:,max_pos)
        ga_res%opt_fitnessvalue = fitness(max_pos)

        call print_ga_summary(opt, ga_res)

        !** Deallocating memory --------------------------
        select case (trim(opt%opt_target_region))   
            case("nuts2")
                deallocate(pos_nuts2)
                if (opt%opt_target_icu) then
                    do i = 1,num_nuts2
                        deallocate(seltarget_icu_by_nuts2(i)%nuts2_icucases)
                        deallocate(seltarget_icu_by_nuts2(i)%nuts2_icudates)
                    enddo
                    deallocate(seltarget_icu_by_nuts2)
                endif

            case("state")
                deallocate(uniq_distid)
                deallocate(pos_se)
                if (opt%opt_target_icu) then
                    do i = 1,num_states
                        deallocate(seltarget_icu_by_state(i)%icucases)
                        deallocate(seltarget_icu_by_state(i)%icudates)
                    enddo
                    deallocate(seltarget_icu_by_state)
                endif
            case default
                write(*,*)  "Otherwise, the given target is not supported yet"
        end select
    
        deallocate(obs_ub)
        deallocate(obs_lb)
        deallocate(eval_lb)
        deallocate(eval_ub)
        deallocate(ini_pop)
        deallocate(fitness)
        
        deallocate(opt_index)
        deallocate(tmp_pop)
        deallocate(tmp_fitness)
        deallocate(sorted_pop)
        deallocate(sorted_fitness)
        deallocate(sorted_index)
        deallocate(summary_fitness) 
        deallocate(sel_par)
        deallocate(opt%opt_names)
        deallocate(opt%opt_names_dur)
        deallocate(ga_res%opt_bestsol)
    end subroutine ga

    !! ------------------------------------------------------------------------------
    !> Subroutine that apply the convolution filter to the input array.
    !>
    subroutine conv_filter(array_filter, coff_len, array_size)
        Real,dimension(array_size),intent(inout)        :: array_filter
        Integer,intent(in)                              :: coff_len, array_size

        Real                              :: coefficient
        Integer                           :: k, j, halve_coff_len
        Real,dimension(:),Allocatable     :: tmp_filter
 
        coefficient = 1.0/coff_len
        Allocate(tmp_filter(array_size-coff_len+1))
        halve_coff_len = coff_len / 2

        k = 1
        do j = halve_coff_len+1,array_size-halve_coff_len
            tmp_filter(k) = sum(array_filter((j-halve_coff_len):(j+halve_coff_len)))*coefficient
            k = k + 1                  
        enddo
        array_filter((halve_coff_len+1):(array_size-halve_coff_len))= tmp_filter

        deallocate(tmp_filter)
    end subroutine conv_filter


    !! --------------------------------------------------------------------------------------------------------------
    !> Subroutine that get a subset of the input array_in (with repeated nuts2) storing the unique nuts2.
    !> 
    !> The returned arrays are array_pos and array_out.
    !> array_out: store the unique nuts2; array_pos: store the position from which a different nuts2 starts in array_in
    function get_unique_nuts2(array_in, array_pos) result(array_out)
    
        integer                                          :: size_in,j,index
        Integer,Allocatable,dimension(:),intent(inout)   :: array_pos
        character*10,allocatable,dimension(:),intent(in) :: array_in
        character*10,allocatable,dimension(:)            :: array_out
        
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
    subroutine date_intersection(obs_lb, eval_lb, obs_ub, eval_ub, seed_date_mod, time_n, dateset)

        Integer,intent(inout)                                  :: obs_lb, obs_ub, eval_lb, eval_ub
        Character*10,dimension(:), Allocatable, intent(in)     :: dateset
        Integer                                                :: inteval_days, status, time_n
        Character*10                                           :: seed_date_mod

        !** Here we assume that the date is monotonically increased by days -----------------------
        obs_lb = 1
        eval_lb = 1

        !TODO: Move the calculation of inteval_days to support module
        inteval_days = Date2Unixtime(dateset(1))&
                            - Date2Unixtime(seed_date_mod)
        inteval_days = inteval_days/86400

        if (inteval_days .LT. 0) then
            obs_lb = obs_lb - inteval_days
        elseif (inteval_days .GT. 0) then
            eval_lb = eval_lb + inteval_days
        endif

        obs_ub = size(dateset)
        eval_ub = time_n

        inteval_days = (Date2Unixtime(dateset(obs_ub))& 
            - (Date2Unixtime(seed_date_mod)+86400*(time_n-1)))/86400
        if (inteval_days .LT. 0) then
            eval_ub = eval_ub + inteval_days
        elseif (inteval_days .GT. 0) then
            obs_ub = obs_ub - inteval_days
        endif

        !** Throw an error when the intersection is null ---------------------------
        if ((eval_lb .GT. eval_ub) .AND. (obs_lb .GT. obs_ub) .AND.&
         ((eval_ub - eval_lb) .NE. (obs_ub - obs_lb))) then
            print *, "intersection is null"
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
        Real,intent(in)                                        :: lb 
        Real,intent(in)                                        :: ub 
        Integer,intent(in)                                     :: popsize 
        Integer,intent(in)                                     :: nvals

        Integer                                                :: j

        do j = 1, popsize
            ini_pop(:,j) = random_uniform(nvals,lb,ub)
        enddo
    end subroutine population

    !! ------------------------------------------------------------------------------
    !> Subroutine that selects the fittest parents
    !>
    !> Return the result indexes to sel_par.
    subroutine selection(fitness, sel_par)
        Real, Allocatable, Dimension(:),intent(inout)         :: fitness
        Integer, Dimension(:),intent(out)                     :: sel_par

        Real, Dimension(size(fitness))                        :: fscaled, prob
        Real                                                  :: fmin, fave, fmax
        Real                                                  :: delta, a, b
        Integer                                               :: sfactor = 2

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

    !    print *,sel_par
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
    !    print *,sample_parents
        
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
        Real,intent(in)                                           :: lb, ub

        Integer,dimension(size(ini_pop,dim=1))                    :: sample_input
        Integer,dimension(1)                                      :: sample_parents
        Integer                                                   :: i
        Real                                                      :: ran
        Real(kind=rk),dimension(1)                                :: Res

        do i = 1,size(sample_input)
            sample_input(i) = i
        enddo
        do i = 1,size(fitness)
            call random_number(ran)
            if (pmutation > ran) then
                sample_parents = sample_i4(sample_input,size(sample_parents)) ! Randomly generate the mutated gene/parameter
                Res = random_uniform(1,lb,ub)
                ini_pop(sample_parents(1),i) = Res(1)  
            !    print *,"the number",i," sample_parents: ",sample_parents(1)," res: ", Res(1)
                fitness(i) = INV ! Invalid the mutated chromosome/population
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
    function sixnum_summary(fitness, sorted_index) Result(summary)
        Real,Allocatable,dimension(:),intent(in)    :: fitness
        Integer,Allocatable,dimension(:),intent(in) :: sorted_index

        type(ga_summary)                            :: summary
        Integer                                     :: pop_size, half, i, j

        pop_size = size(fitness)
        summary%max = fitness(sorted_index(pop_size))
        summary%min = fitness(sorted_index(1))
        summary%mean = Sum(fitness)/pop_size

        summary%median = get_median(fitness, sorted_index, pop_size)
        half = pop_size/2.0
        summary%upper_hinge = get_median(fitness, sorted_index((half+1):pop_size), pop_size-half) ! Get the median of the second half
        summary%lower_hinge = get_median(fitness, sorted_index(1:half), half) ! Get the median of the first half
    end function sixnum_summary
   
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

    !! ------------------------------------------------------------------------------
    !> Function that prints the GA settings and results
    !>
    !> The formats still lack
    subroutine print_ga_summary(opt, ga_res)
        type(opt_parameters),intent(in)             :: opt
        type(opt_res),intent(in)                    :: ga_res

        write(*,'(A)')"!----------------------- Genetic Algorithm ----------------------------!"
        write(*,'(A)')"!------ GA settings ------!"
        print *,"Type: ","real-valued"
        print *,"Population size: ",opt%opt_pop_size
        print *,"Number of generations: ",opt%opt_num_gene ! maxiter
        print *,"Elitism: ",opt%opt_elitism
        print *,"Crossover probability: ",opt%opt_pcrossover
        print *,"Mutation probability: ",opt%opt_pmutation
        print *,"Optimization target: ",opt%opt_target_region
        write(*,'(A)')"!------ GA results -------!"
        print *,"Iterations: ",ga_res%opt_acuiter
        print *,"Fitness function value: ",ga_res%opt_fitnessvalue
        print *,"Solution:"
        print *,ga_res%opt_bestsol
        print *,"Affected nuts2: ",opt%opt_names
        print *,"Affected duration for each ",opt%opt_names_dur
    end subroutine print_ga_summary

  
End Module genetic_algorithm
