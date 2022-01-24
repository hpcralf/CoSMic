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

    Real, parameter :: INV = 99.0 ! Corresponds to NA for fitness array 
    Real, parameter :: EPS = EPSILON(1.0) ! Tolerance precision, a magnitude of e-07

    !** Optimization Control Panel ----------------------
    type opt_parameters			
        Logical                                       :: opt_target_icu
        Logical                                       :: opt_target_deaths ! Optimization targets
        Logical                                       :: opt_filter  = .TRUE.
        character(len=mcl)							  :: use_sug_sol = "NULL"
        character(len=8)							  :: opt_target_region = "state" ! Optimized region. Valid values: "country", "state" and "nuts2"
        character(len=8), Dimension(:), Allocatable   :: opt_names ! Names of parameters (by Nuts-2) to optimize
        Integer,Dimension(:), Allocatable             :: opt_names_dur ! The duration of certain parameter is calculated by weeks
        Real                        				  :: opt_lb = 0.1 ! Lower bounds of optimized parameters 
        Real                                          :: opt_ub = 1.0 ! Upper bounds of optimized parameters
        Integer 									  :: opt_pop_size = 8 ! Population size
        Integer    								      :: opt_num_gene = 4 ! Maximal generation
        Integer                                       :: opt_elitism ! The number of the chosen eltism
    end type opt_parameters

    !** Observed icu cases by state -------------------------------
    type target_obsicu
        Real,dimension(:), Allocatable                :: icucases
        Character*10,dimension(:), Allocatable        :: icudates
    end type target_obsicu

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
        Real                                          :: pcrossover = 0.8
        Real                                          :: pmutation  = 0.1
        Integer                                       :: nvals ! The total number of parameters to be optimized
        Type(opt_parameters)                          :: opt

        Real(kind=rk),Dimension(:,:), Allocatable     :: ini_pop, tmp_pop, sorted_pop ! ini_pop: store the population, which are the input to be optimized. 

        Real, Allocatable, Dimension(:)               :: fitness, tmp_fitness ! fitness: store the calculated fitness
        Real, Allocatable, Dimension(:)               :: sorted_fitness, summary_fitness

        Integer,dimension(:,:,:), Allocatable         :: ill_ICU_cases_final ! Store the evaluated ICU cases
        Integer,dimension(:), Allocatable             :: uniq_distid
        Integer, dimension(:), Allocatable            :: opt_index
        Integer                                       :: i, j, k, pop_ii, gene_ii
        Integer                                       :: startid, endid, max_pos, min_pos
        Integer                                       :: run_stop ! For the convergence, by default equals to opt_num_gene
        Real                                          :: tmp_max, run_sum
        Integer                                       :: time_n, num_states
        Character*15,Dimension(1)                     :: given_sc
        Type(target_obsicu),dimension(:),Allocatable  :: seltarget_icu_by_state ! Record the observed ICU cases by state
        Real,dimension(:), Allocatable                :: diff
        Real                                          :: diff_sum = 0.0
        integer                                       :: diff_len = 0
        Integer,dimension(:,:,:), Allocatable         :: aggreg_sum_by_county
        Integer,dimension(:), Allocatable             :: pos_se
        Real,dimension(:,:), Allocatable              :: avg_by_county ! Record the simulated ICU cases by state
        Integer,dimension(:), Allocatable             :: state_sc_index

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
        ! Selection
        !=================================================
        Integer,dimension(:), Allocatable             :: sel_par ! The selected parents

        !=================================================
        ! Sorting
        !=================================================
        integer,dimension(:),Allocatable              :: sorted_index
        integer                                       :: idx
        real                                          :: pre_fitv, curr_fitv


        ! TODO: parameters checking
        ! ====================================================================
        ! ====================================================================


        !** Initialization of opt ----- TODO: should be written as input data
        opt%opt_target_icu    = .TRUE.
        opt%opt_target_deaths = .FALSE.
        opt%opt_names = (/ "def0", "de60", "de91", "de92", "de93", "de94", "de50" /)
        opt%opt_names_dur = (/ 6, 6, 6, 6, 6 ,6 ,6 /)
        nvals = sum(opt%opt_names_dur)
        
        opt%opt_elitism = max(1,NINT(opt%opt_pop_size*0.05))

        
        Allocate(ini_pop(nvals, opt%opt_pop_size))
        Allocate(fitness(opt%opt_pop_size))
        Allocate(tmp_pop(nvals, opt%opt_pop_size))
        Allocate(tmp_fitness(opt%opt_pop_size))
        Allocate(opt_index(size(opt%opt_names)))

        allocate(sorted_pop(nvals, opt%opt_elitism))
        allocate(sorted_fitness(opt%opt_elitism))
        allocate(sorted_index(opt%opt_pop_size))

        allocate(summary_fitness(opt%opt_num_gene))  !tempaporily store the maximum of fitness among all the population for each generation

        time_n = maxval(R0change) + 1

        uniq_distid = get_unique(counties_index/1000) ! Get the affected states
        num_states = size(uniq_distid) ! Get the number of the affected states
        Allocate(pos_se(num_states+1))
        
        Allocate(seltarget_icu_by_state(num_states))
        Allocate(obs_lb(num_states))
        Allocate(obs_ub(num_states))
        Allocate(eval_lb(num_states))
        Allocate(eval_ub(num_states))

        allocate(sel_par(opt%opt_pop_size))

        !** Collect the beginning index from which a different state starts ---------------------
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


        !** When the optimized target region is "state" -------------------
        if (trim(opt%opt_target_region) == "state") then
            !** The fitness relavant to the death data ------------

            !** The fitness relavant to the icu data --------------
            if (opt%opt_target_icu) then
                do i = 1,num_states
                    given_sc(1) = iol%states_shortcut(uniq_distid(i))
                    state_sc_index = get_index_mul_char(iol%obsicu_states_shortcut,given_sc)

                    !** Presume seltarget_icudates_by_state is arranged in an ascending order ----------
                    !** Preprocess the observed ICU data according to the affected states --------------
                    seltarget_icu_by_state(i)%icucases = Real(iol%obsicu_cases(state_sc_index))
                    seltarget_icu_by_state(i)%icudates = iol%obsicu_date(state_sc_index)

                    !** Apply convolution filter ------
                    if (opt%opt_filter) then
                        coefficient = 1.0/coff_len
                        Allocate(tmp_filter(size(state_sc_index)-coff_len+1))
                        halve_coff_len = coff_len / 2

                        k = 1
                        do j = halve_coff_len+1,size(state_sc_index)-halve_coff_len
                            tmp_filter(k) = sum(seltarget_icu_by_state(i)%icucases((j-halve_coff_len):(j+halve_coff_len)))&
                            *coefficient
                            k = k + 1                  
                        enddo
                        seltarget_icu_by_state(i)%icucases((halve_coff_len+1):(size(state_sc_index)-halve_coff_len)) = tmp_filter

                        deallocate(tmp_filter)
                    endif
                    deallocate(state_sc_index)

                    !** Calculate the intersection of the evaluated and observed ICU data by date -------------
                    !** Here we assume that the date is monotonically increased by days -----------------------
                    seed_date_mod        = add_date(seed_date,1)
                    obs_lb(i) = 1
                    eval_lb(i) = 1

                    !TODO: move the calculation of inteval_days to support module
                    inteval_days = Date2Unixtime(seltarget_icu_by_state(i)%icudates(1))&
                                        - Date2Unixtime(seed_date_mod)
                    inteval_days = inteval_days/86400
                    
                    if (inteval_days .LT. 0) then
                        obs_lb(i) = obs_lb(i) - inteval_days
                    elseif (inteval_days .GT. 0) then
                        eval_lb(i) = eval_lb(i) + inteval_days
                    endif

                    obs_ub(i) = size(seltarget_icu_by_state(i)%icudates)
                    eval_ub(i) = time_n

                    inteval_days = (Date2Unixtime(seltarget_icu_by_state(i)%icudates(obs_ub(i)))& 
                        - (Date2Unixtime(seed_date_mod)+86400*(time_n-1)))/86400
                    if (inteval_days .LT. 0) then
                        eval_ub(i) = eval_ub(i) + inteval_days
                    elseif (inteval_days .GT. 0) then
                        obs_ub(i) = obs_ub(i) - inteval_days
                    endif

                    !** Throw an error when the intersection is null ---------------------------
                    if ((eval_lb(i) .GT. eval_ub(i)) .AND. (obs_lb(i) .GT. obs_ub(i)) .AND.&
                     ((eval_ub(i) - eval_lb(i)) .NE. (obs_ub(i) - obs_lb(i)))) then
                        print *, "intersection is null"
                        print *, "eval_lb: ",eval_lb(i), " eval_ub: ",eval_ub(i)," obs_lb: ",obs_lb(i)," obs_ub: ",obs_ub(i)
                        Call exit(status)    
                    endif
                enddo
            end if 
        end if

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
            print *,"gene_ii: ",gene_ii
            do pop_ii = 1,opt%opt_pop_size ! Iterate over all the population
                print *,"pop_ii: ",pop_ii
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
                if (trim(opt%opt_target_region) == "state") then
                    !** All the data are represented by state --------
                    Allocate(aggreg_sum_by_county(num_states,time_n,iter))
                    Allocate(avg_by_county(num_states,time_n))

                    !** TODO: The fitness relavant to the death data ----------

                    !** the fitness relavant to the icu data ------------------
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
                                avg_by_county(i,j) = sum(aggreg_sum_by_county(i,j,:))/Real(iter) ! Average over the runing iterates
                            enddo
                        enddo

                        !** Calculate the difference between the observed and evaluated/simulated ICU data -----------
                        do i = 1,num_states                  
                            if (opt%opt_filter) then                   
                                ! ====================================================================================!
                                ! =========== TODO: the following code is repeated and should be modularied ==========!
                                Allocate(tmp_filter(size(avg_by_county(i,:))-coff_len+1))
                                k = 1
                               
                                do j = halve_coff_len+1,size(avg_by_county(i,:))-halve_coff_len
                                    tmp_filter(k) = sum(avg_by_county(i,(j-halve_coff_len):(j+halve_coff_len)))*coefficient
                                    k = k + 1                  
                                enddo
                                
                                avg_by_county(i,(halve_coff_len+1):(size(avg_by_county(i,:))-halve_coff_len)) = tmp_filter
                                deallocate(tmp_filter)
                                
                            endif
                            diff = seltarget_icu_by_state(i)%icucases(obs_lb(i):obs_ub(i)) - avg_by_county(i,eval_lb(i):eval_ub(i))
                            diff = diff*diff
                            diff_sum = diff_sum + sum(diff)
                            diff_len = diff_len + obs_ub(i) - obs_lb(i) + 1
                            deallocate(diff)
                        enddo
                        fitness(pop_ii) = -1*sqrt(diff_sum/diff_len)
                         
                    end if
                    deallocate(aggreg_sum_by_county)
                    deallocate(avg_by_county)
                end if
                deallocate(ill_ICU_cases_final)
            enddo

            do i = 1,opt%opt_pop_size
                print *,i,": ",fitness(i)
            enddo

            ! ========================================================== !
            ! ============ Check the stop criteria ===================== !
            call SORTRX_REAL(opt%opt_pop_size,fitness,sorted_index) ! Sorted_index indicates a sorted fitness in an ascending order
            
            run_sum = 0
            max_pos = sorted_index(opt%opt_pop_size)
            min_pos = sorted_index(1)
            summary_fitness(gene_ii) = fitness(max_pos)
            if (gene_ii > 1) then
                tmp_max = maxval(summary_fitness(1:gene_ii))
                do i = 1,gene_ii 
                    if (summary_fitness(i) .ge. (tmp_max - EPS)) then
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

            print *,"elitism:"
            print *,sorted_fitness
            print *,sorted_pop

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
            call crossover(ini_pop, fitness, pcrossover)
            do i = 1,nvals
                print *,ini_pop(i,:)
            enddo

            print *,fitness

            !** To randomly mutate a parameter/gene ------------------------------------------------------
            if(pmutation .ne. 0.0) then
                call mutation(ini_pop, fitness, pmutation, opt%opt_lb, opt%opt_ub)

                do i = 1,nvals
                    print *,ini_pop(i,:)
                enddo
                print *,fitness
            endif
            !** Replace the unfittest populations with the fittest ones, which are stored in sorted_pop ---
            call elitism(fitness, ini_pop, sorted_pop, sorted_fitness)
        enddo

        deallocate(pos_se)
        deallocate(seltarget_icu_by_state)
        deallocate(obs_ub)
        deallocate(obs_lb)
        deallocate(eval_lb)
        deallocate(eval_ub)
        deallocate(ini_pop)
        deallocate(fitness)
        deallocate(uniq_distid)
        deallocate(opt_index)
        deallocate(tmp_pop)
        deallocate(tmp_fitness)
        deallocate(sorted_pop)
        deallocate(sorted_fitness)
        deallocate(sorted_index)
        deallocate(summary_fitness) 
        deallocate(sel_par)
        !===== TODO: if we need to deallocate the member seltarget_icu_by_state(i)%icucases ====
    end subroutine ga

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

        print *,"comare: ",(sfactor*fave - fmax)/(sfactor - 1)
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
        prob = abs(fscaled)/sum(abs(fscaled))
        print *,"a: ",a," b: ",b," prob: ",prob
        sel_par = sample_weight(size(fitness), prob)

        print *,sel_par
    end subroutine selection

    !! ------------------------------------------------------------------------------
    !> Subroutine that mate the selected pair and breed two children
    !>
    !> The results are stored back to ini_pop and fitness.
    subroutine crossover(ini_pop, fitness, pcrossover) 
        Real(kind=rk),dimension(:,:),Allocatable,intent(inout)    :: ini_pop
        Real,dimension(:),Allocatable,intent(inout)               :: fitness
        Real,intent(in)                                           :: pcrossover

        Integer,dimension(:),Allocatable                          :: sample_input, sample_parents
        Integer                                                   :: start_id, end_id, nmating,i,j
        Real(kind=rk)                                             :: ran
        Real(kind=rk),dimension(size(ini_pop,dim=1))              :: lchild, rchild

        nmating = size(fitness)/2 !get the floor
        allocate(sample_input(nmating*2))
        allocate(sample_parents(nmating*2))

        print *,size(ini_pop,dim=1)
        do i = 1,size(lchild)
            print *,ini_pop(i,:)
        enddo
        print *,fitness

        do i = 1, 2*nmating 
            sample_input(i) = i
        enddo
      
        sample_parents = sample_i4(sample_input,2*nmating)
        print *,sample_parents
        
        do i = 1, nmating 
            call random_number(ran)
            print *,"crossover ran: ",ran
            end_id = i * 2 
            start_id = end_id - 1 
            start_id = sample_parents(start_id)
            end_id = sample_parents(end_id)

            if((pcrossover > ran) .and.&
             (fitness(start_id) .ne. fitness(end_id))) then
                call random_number(ran)
                print *,"ran: ",ran
                lchild = ran*ini_pop(:,start_id) + (1.0-ran)*ini_pop(:,end_id)
                rchild = (1.0-ran)*ini_pop(:,start_id) + ran*ini_pop(:,end_id)

                ini_pop(:,start_id) = lchild
                ini_pop(:,end_id)   = rchild

                fitness(start_id) = INV 
                fitness(end_id)   = INV 
            endif
        enddo
    end subroutine crossover

    !! ------------------------------------------------------------------------------
    !> Subroutine that randomly mutate certain parameter
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
        print *,"mutation"
        do i = 1,size(fitness)
            call random_number(ran)
            if (pmutation > ran) then
                sample_parents = sample_i4(sample_input,size(sample_parents))
                Res = random_uniform(1,lb,ub)
                ini_pop(sample_parents(1),i) = Res(1)
                print *,"the number",i," sample_parents: ",sample_parents(1)," res: ", Res(1)
                fitness(i) = INV
            endif
        enddo

    end subroutine mutation

    !! ------------------------------------------------------------------------------
    !> Subroutine that replace the unfittest individules with the chosen elitim
    !>
    !> The results are stored back to ini_pop and fitness.
    subroutine elitism(fitness, ini_pop, sorted_pop, sorted_fitness)
        Real(kind=rk),dimension(:,:),Allocatable,intent(inout)    :: ini_pop
        Real,dimension(:),Allocatable,intent(inout)               :: fitness
        Real(kind=rk),dimension(:,:),Allocatable,intent(in)       :: sorted_pop                                   
        Real,dimension(:),Allocatable,intent(in)                  :: sorted_fitness
        Integer,dimension(size(fitness))                          :: sorted_index
        Integer                                                   :: num

        print *,"elitism"

        call SORTRX_REAL(size(fitness),fitness,sorted_index)
        num = size(sorted_fitness)
        fitness(sorted_index(1:num)) = sorted_fitness
        ini_pop(:,sorted_index(1:num)) = sorted_pop
        
        print *,fitness 
        do num = 1,size(ini_pop,1)
           print *,ini_pop(num,:)
        enddo


    end subroutine elitism

    !! ------------------------------------------------------------------------------
    !> Function that return a sample with the weight in mind
    !>
    !> The weights/probabilities are stored in prop.
    function sample_weight(input, prop)Result(sel_par)
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


    subroutine loadobsdata(iol)

        Type(iols),intent(INOUT)   :: iol

        Integer                    :: i,j,k,it_ss
        Integer                    :: index, un_in
       

        !=================================================================
        !Character(len=:), allocatable :: data_dir
        Character(len=:), allocatable :: filename
        !=================================================================

        !call pt_get("#data_dir",data_dir)

        ! Read data from file: observed icu cases by state ------------------
        !call pt_get("#observ_icu_state",filename)

        filename = "./ICU_cases_Germany_09_11_bundesland_new.csv"
        call open_and_index(trim(filename),un_in,index)

        ! allocation for the readin variables ------------------
        Allocate(iol%obsicu_date(index-1))
        Allocate(iol%obsicu_cases(index-1))
        Allocate(iol%obsicu_states_shortcut(index-1))

        ! read the first line(character) -----------------------
        Read(un_in,*) iol%titel

        Do i = 1, index-1
           Read(un_in,*) iol%obsicu_date(i),iol%obsicu_cases(i),&
                iol%obsicu_states_shortcut(i)
        End Do
        Close(un_in)

    !    write(*,*)"the observed icu cases by state",iol%obsicu_date,iol%obsicu_cases,iol%obsicu_states_shortcut

    end subroutine loadobsdata

End Module genetic_algorithm
