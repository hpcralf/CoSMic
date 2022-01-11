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
    use param_tree
    use cosmic_io
    use kernel

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
        Integer 									  :: opt_pop_size = 4
        Integer    								      :: opt_num_gene = 2

    end type opt_parameters

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

        Type(iols)                                                :: iol
        Integer(kind=ik)             , Dimension(:)  , intent(in) :: counties_index
        Integer(kind=ik)                             , intent(in) :: iter
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
        Real(kind=rk),    Allocatable, Dimension(:,:), intent(in) :: R0_effects

        !---------------------------------------------------------------------------
        Real                                          :: pcrossover = 0.8
        Real                                          :: pmutation  = 0.1
        Integer                                       :: nvals
        Type(opt_parameters)                          :: opt

        Real, Allocatable, Dimension(:)               :: fitness

        Real,Dimension(:,:), Allocatable              :: ini_pop
        Integer,dimension(:,:,:), Allocatable         :: ill_ICU_cases_final
        Integer,dimension(:), Allocatable             :: uniq_distid
        Integer, dimension(7)                         :: opt_index
        Integer                                       :: i, j, k
        Integer                                       :: startid, endid

        !--------------------------------------------------------------------------
        Integer,dimension(:,:,:), Allocatable         :: aggreg_sum_by_county
        Integer,dimension(:), Allocatable             :: pos_se
        Real,dimension(:,:), Allocatable              :: avg_by_county
        Integer,dimension(:), Allocatable             :: state_sc_index
        Real,dimension(:), Allocatable                :: seltarget_icucases_by_state
        Real,dimension(:), Allocatable                :: tmp_filter
        Character*10,dimension(:), Allocatable        :: seltarget_icudates_by_state
        Integer                                       :: time_n, num_states
        Character(len=:),Dimension(:),allocatable     :: sim_regions
        Character(len=8),Dimension(1)                 :: given_sc

        Character*10                                  :: seed_date_mod
        Real,dimension(:), Allocatable                :: diff
        Real                                          :: diff_sum = 0.0
        integer                                       :: diff_len = 0
        Integer                                       :: halve_coff_len, coff_len=7
        Integer                                       :: obs_lb, obs_ub, eval_lb, eval_ub, inteval_days, status
        Real                                          :: coefficient


        ! TODO: parameters checking
        ! ====================================================================
        ! ====================================================================


        ! Initialization of opt. TODO: should be written as input data
        opt%opt_target_icu    = .TRUE.
        opt%opt_target_deaths = .FALSE.
    !    Allocate(opt%opt_names(7))
        Allocate(opt%opt_names_dur(7))
        opt%opt_names = (/ "def0", "de60", "de91", "de92", "de93", "de94", "de50" /)
        opt%opt_names_dur = (/ 6, 6, 6, 6, 6 ,6 ,6 /)

        nvals = sum(opt%opt_names_dur)

        ! Generate beginning population
        Allocate(ini_pop(nvals, opt%opt_pop_size))
        Allocate(fitness(opt%opt_pop_size))

        ! Load the observed icu cases data
        !call loadobsdata(iol) 
        !print *,iol%obsicu_cases
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

        startid = 1
        do i=1,size(opt%opt_names)
            endid = startid + opt%opt_names_dur(i) - 1
            iol%R0_effect%data(1:opt%opt_names_dur(i),opt_index(i)) = ini_pop(startid:endid,1)
            startid = endid + 1
        enddo

        
        !do j = 1, nvals
        !    print *,iol%R0_effect%data(:,j)
        !    print *,j
        !enddo

        Call COVID19_Spatial_Microsimulation_for_Germany(iol,counties_index, &
            iter , &
            inf_dur, cont_dur, ill_dur, icu_dur, icu_per_day, &
            less_contagious, R0_force, immune_stop, &
            R0change, R0delay ,R0delay_days, R0delay_type, &
            control_age_sex, seed_date, seed_before, sam_size, R0, &
            iol%R0_effect%data, ill_ICU_cases_final)

        !==================================================
        !====== calculate fitness function value ==========
        if (trim(opt%opt_target_region) == "state") then
            ! ===== the fitness relavant to the death data ======

            ! ===== the fitness relavant to the icu data ========
            if (opt%opt_target_icu) then
            ! TODO: record the number of the affected states
                time_n = maxval(R0change) + 1
                uniq_distid = get_unique(counties_index/1000)
                num_states = size(uniq_distid)
                Allocate(pos_se(num_states+1))
                Allocate(aggreg_sum_by_county(num_states,time_n,iter))
                Allocate(avg_by_county(num_states,time_n)) !TODO: note that it is column-majored

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

                do i = 1,num_states
                    given_sc(1) = iol%states_shortcut(uniq_distid(i))
                    state_sc_index = get_index_mul_char(iol%obsicu_states_shortcut,given_sc)
                    ! presume seltarget_icudates_by_state is arranged in an ascending order
                    seltarget_icucases_by_state = iol%obsicu_cases(state_sc_index)
                    seltarget_icudates_by_state = iol%obsicu_date(state_sc_index)
                    deallocate(state_sc_index)

                    if (opt%opt_filter) then
                        coefficient = 1.0/coff_len
                        Allocate(tmp_filter(size(state_sc_index)-coff_len+1))
                        halve_coff_len = coff_len / 2

                        k = 1
                        do j = halve_coff_len+1,size(state_sc_index)-halve_coff_len
                            tmp_filter(k) = sum(seltarget_icucases_by_state((j-halve_coff_len):(j+halve_coff_len)))*coefficient
                            k = k + 1                  
                        enddo
                        seltarget_icucases_by_state((halve_coff_len+1):(size(state_sc_index)-halve_coff_len)) = tmp_filter

                    !    print *, seltarget_icucases_by_state
                        deallocate(tmp_filter)
                         
                        ! ====================================================================================!
                        ! =========== TODO: the following code is repeated and should be modularied ==========!
                        Allocate(tmp_filter(size(avg_by_county(i,:))-coff_len+1))
                        k = 1
                        do j = halve_coff_len+1,size(avg_by_county(i,:))-halve_coff_len
                            tmp_filter(k) = sum(avg_by_county(i,(j-halve_coff_len):(j+halve_coff_len)))*coefficient
                            k = k + 1                  
                        enddo
                        ! size(avg_by_county(i,:)) equals to time_n
                        avg_by_county(i,(halve_coff_len+1):(size(avg_by_county(i,:))-halve_coff_len)) = tmp_filter
                        !    print *, avg_by_county(i,:)
                        deallocate(tmp_filter)
                    endif

                    !============= intersection ===================
                    seed_date_mod        = add_date(seed_date,1)
                    obs_lb = 1
                    eval_lb = 1

                    !TODO: move the calculation of inteval_days to support module
                    !==== here we assume that the date is increased by days
                    inteval_days = Date2Unixtime(seltarget_icudates_by_state(1))&
                                        - Date2Unixtime(seed_date_mod)
                    inteval_days = inteval_days/86400
                    
                    if (inteval_days .LT. 0) then
                        obs_lb = obs_lb - inteval_days
                    elseif (inteval_days .GT. 0) then
                        eval_lb = eval_lb + inteval_days
                    endif

                    obs_ub = size(seltarget_icudates_by_state)
                    eval_ub = time_n

                    inteval_days = (Date2Unixtime(seltarget_icudates_by_state(obs_ub))& 
                        - (Date2Unixtime(seed_date_mod)+86400*(time_n-1)))/86400
                    if (inteval_days .LT. 0) then
                        eval_ub = eval_ub + inteval_days
                    elseif (inteval_days .GT. 0) then
                        obs_ub = obs_ub - inteval_days
                    endif
                    
                    if ((eval_lb .le. eval_ub) .and. (obs_lb .le. obs_ub) .and. ((eval_ub - eval_lb) .eq. (obs_ub-obs_lb))) then
                        !Allocate(diff(obs_ub-obs_lb+1))
                        diff = seltarget_icucases_by_state(obs_lb:obs_ub) - avg_by_county(i,eval_lb:eval_ub)
                        diff = diff*diff
                        diff_sum = diff_sum + sum(diff)
                        diff_len = diff_len + obs_ub - obs_lb + 1

                        deallocate(diff)
                    else 
                        print *, "intersection is null"
                        print *, "eval_lb: ",eval_lb, " eval_ub: ",eval_ub," obs_lb: ",obs_lb," obs_ub: ",obs_ub
                        Call exit(status)
                    end if
                    deallocate(seltarget_icucases_by_state)
                    deallocate(seltarget_icudates_by_state)
                enddo
                fitness(1) = -1*sqrt(diff_sum/diff_len)
                print *,fitness(1)
                
                deallocate(pos_se)
                deallocate(aggreg_sum_by_county)
                deallocate(avg_by_county)
                
              
            end if
        end if
    end subroutine ga

    subroutine population(ini_pop, lb, ub, popsize, nvals)
        Real,dimension(:,:)             :: ini_pop
        Real                            :: lb 
        Real                            :: ub 
        Integer                         :: popsize 
        Integer                         :: nvals

        Integer                         :: j

        do j = 1, popsize
            ini_pop(:,j) = random_uniform(nvals,lb,ub)
        end do


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
