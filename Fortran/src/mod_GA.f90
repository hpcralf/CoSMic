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

    Use data_preprocessing
    Use cosmic_io

    implicit none

    type opt_parameters
    					
        Logical                                       :: opt_target_icu
        Logical                                       :: opt_target_deaths
        Logical                                       :: opt_filter  = .TRUE.
        character(len=mcl)							  :: use_sug_sol = "NULL"
        character(len=8)							  :: opt_target_region = "state"
        character(len=8), Dimension(:), Allocatable   :: opt_names 
        Integer,Dimension(:), Allocatable             :: opt_names_dur
        Real                        				  :: opt_lb = 0.1
        Real                                          :: opt_ub = 1.0
        Integer 									  :: opt_pop_size = 8
        Integer    								      :: opt_num_gene = 2

    end type opt_parameters

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
        character(len=:), Allocatable             ,intent(in)    :: seed_date
        Integer(kind=ik)                          ,intent(in)    :: seed_before    
        Integer(kind=ik)                          ,intent(in)    :: sam_size
        Real(kind=rk)                              ,intent(in)   :: R0
        Real(kind=rk),    Allocatable, Dimension(:,:) ,intent(in):: R0_effects

        !---------------------------------------------------------------------------
        Real                                          :: pcrossover = 0.8
        Real                                          :: pmutation  = 0.1
        Integer                                       :: nvals
        Type(opt_parameters)                          :: opt

        Real, Allocatable, Dimension(:)               :: fitness

        Real(kind=rk),Dimension(:,:), Allocatable     :: ini_pop ! in line with the R0_effect
        Integer,dimension(:,:,:), Allocatable         :: ill_ICU_cases_final
        Integer,dimension(:), Allocatable             :: uniq_distid
        Integer, dimension(:), Allocatable            :: opt_index
        Integer                                       :: i, j, k, pop_ii
        Integer                                       :: startid, endid

        !--------------------------------------------------------------------------
        Integer,dimension(:,:,:), Allocatable         :: aggreg_sum_by_county
        Integer,dimension(:), Allocatable             :: pos_se
        Real,dimension(:,:), Allocatable              :: avg_by_county
        Integer,dimension(:), Allocatable             :: state_sc_index
        Type(target_obsicu),dimension(:),Allocatable  :: seltarget_icu_by_state
        Real,dimension(:), Allocatable                :: tmp_filter
        
        Integer                                       :: time_n, num_states
        Character(len=:),Dimension(:),allocatable     :: sim_regions
        Character(len=8),Dimension(1)                 :: given_sc

        Character*10                                  :: seed_date_mod
        Real,dimension(:), Allocatable                :: diff
        Real                                          :: diff_sum = 0.0
        integer                                       :: diff_len = 0
        Integer                                       :: halve_coff_len, coff_len=7
        Integer,dimension(:), Allocatable             :: obs_lb, obs_ub, eval_lb, eval_ub
        Integer                                       :: inteval_days, status
        Real                                          :: coefficient


        ! TODO: parameters checking
        ! ====================================================================
        ! ====================================================================


        ! Initialization of opt. TODO: should be written as input data
        opt%opt_target_icu    = .TRUE.
        opt%opt_target_deaths = .FALSE.
        opt%opt_names = (/ "def0", "de60", "de91", "de92", "de93", "de94", "de50" /)
        opt%opt_names_dur = (/ 6, 6, 6, 6, 6 ,6 ,6 /)
        nvals = sum(opt%opt_names_dur)

        ! Generate beginning population
        Allocate(ini_pop(nvals, opt%opt_pop_size))
        Allocate(fitness(opt%opt_pop_size))
        Allocate(opt_index(size(opt%opt_names)))

        time_n = maxval(R0change) + 1
        uniq_distid = get_unique(counties_index/1000)
        num_states = size(uniq_distid)
        Allocate(pos_se(num_states+1))
        
        Allocate(seltarget_icu_by_state(num_states))
        Allocate(obs_lb(num_states))
        Allocate(obs_ub(num_states))
        Allocate(eval_lb(num_states))
        Allocate(eval_ub(num_states))

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


         if (trim(opt%opt_target_region) == "state") then
            ! ===== the fitness relavant to the death data ======

            ! ===== the fitness relavant to the icu data ========
            if (opt%opt_target_icu) then
                do i = 1,num_states
                    given_sc(1) = iol%states_shortcut(uniq_distid(i))
                    state_sc_index = get_index_mul_char(iol%obsicu_states_shortcut,given_sc)

                    ! presume seltarget_icudates_by_state is arranged in an ascending order
                    seltarget_icu_by_state(i)%icucases = Real(iol%obsicu_cases(state_sc_index))
                    seltarget_icu_by_state(i)%icudates = iol%obsicu_date(state_sc_index)
                    !deallocate(state_sc_index)

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

                    !============= intersection ===================

                    seed_date_mod        = add_date(seed_date,1)
                    obs_lb(i) = 1
                    eval_lb(i) = 1

                    !TODO: move the calculation of inteval_days to support module
                    !==== here we assume that the date is increased by days
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

                    if ((eval_lb(i) .GT. eval_ub(i)) .AND. (obs_lb(i) .GT. obs_ub(i)) .AND.&
                     ((eval_ub(i) - eval_lb(i)) .NE. (obs_ub(i) - obs_lb(i)))) then
                        print *, "intersection is null"
                        print *, "eval_lb: ",eval_lb(i), " eval_ub: ",eval_ub(i)," obs_lb: ",obs_lb(i)," obs_ub: ",obs_ub(i)
                        Call exit(status)    
                    endif
                enddo

            end if 
        end if

        !============== preprocessing: update the data subject to optimization =================================
        do j=1,size(opt%opt_names)
            do i=1,size(iol%R0_effect%head)
                if (index(trim(iol%R0_effect%head(i)), trim(opt%opt_names(j))) .NE. 0) then
                    opt_index(j) = i
                    EXIT
                endif
            enddo
        enddo

        call population(ini_pop, opt%opt_lb, opt%opt_ub, opt%opt_pop_size, nvals)
        
        do pop_ii = 1,opt%opt_pop_size
            startid = 1
            do i=1,size(opt%opt_names)
                endid = startid + opt%opt_names_dur(i) - 1
                iol%R0_effect%data(1:opt%opt_names_dur(i),opt_index(i)) = ini_pop(startid:endid,pop_ii)
                startid = endid + 1
            enddo

            Call COVID19_Spatial_Microsimulation_for_Germany(iol,counties_index, &
                iter , &
                inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
                less_contagious, R0_force, immune_stop, &
                R0change, R0delay ,R0delay_days, R0delay_type, &
                control_age_sex, seed_date, seed_before, sam_size, R0, &
                iol%R0_effect%data,ill_ICU_cases_final)


            diff_sum = 0.0
            diff_len = 0

            !==================================================
            !====== calculate fitness function value ==========
            if (trim(opt%opt_target_region) == "state") then
                ! ===== the fitness relavant to the death data ======
                Allocate(aggreg_sum_by_county(num_states,time_n,iter))
                Allocate(avg_by_county(num_states,time_n)) !TODO: note that it is column-majored

                ! ===== the fitness relavant to the icu data ========
                if (opt%opt_target_icu) then
                ! TODO: record the number of the affected states
                    do j = 1,iter
                        do i = 1,time_n 
                            do k = 1,num_states
                                aggreg_sum_by_county(k,i,j) = sum(ill_ICU_cases_final(pos_se(k):pos_se(k+1)-1,i,j))
                            enddo
                        enddo
                    enddo

                    do j = 1,time_n
                        do i = 1,num_states
                            avg_by_county(i,j) = sum(aggreg_sum_by_county(i,j,:))/Real(iter)
                        enddo
                    enddo

                    !do j = 1,time_n 
                    !    print *,avg_by_county(:,j)
                    !enddo


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
        !===== TODO: if we need to deallocate the member seltarget_icu_by_state(i)%icucases ====
    end subroutine ga

    subroutine population(ini_pop, lb, ub, popsize, nvals)
        Real(kind=rk),dimension(:,:),Allocatable :: ini_pop
        Real                                     :: lb 
        Real                                     :: ub 
        Integer                                  :: popsize 
        Integer                                  :: nvals

        Integer                                  :: j

        do j = 1, popsize
            ini_pop(:,j) = random_uniform(nvals,lb,ub)
        enddo


    end subroutine population

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
