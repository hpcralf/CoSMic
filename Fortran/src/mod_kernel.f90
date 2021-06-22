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
! Module containing the CoSMic model loop
!
!###############################################################################
Module kernel

  Use global_constants
  use global_types
  Use list_variable
  Use support_fun
    
  Implicit None

Contains
  
  Subroutine COVID19_Spatial_Microsimulation_for_Germany( &
       sp, &
       iol, pspace, iter, counties_index)

    Include 'mpif.h'

    Type(static_parameters) :: sp

    Type(pspaces)       :: pspace
    Type(iols)          :: iol
    Integer             :: i, j, k, index, temp_int,icounty,county,it_ss,iter,status
    Character*1         :: mod1,mod2
    Integer,Dimension(:):: counties_index
    Logical             :: seed_in_inner_loop,seed_mpi_in_inner_loop
!!!!!-----1.states of the model ------
    Integer             :: healthy,inf_noncon,inf_contag,ill_contag,ill_ICU,immune,dead,missing

!!!!!-----2.derive the initial population ------
    Real                :: sam_prop_ps(16)=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/) ! scaling coeff for population
    Character*15        :: sim_pop

!!!!!-----3.Set seed infections in population ------
    Integer             :: ini_infected
    Character*5         :: seed_infections
    Integer             :: seed_before

!!!!!-----4. define disease characteristics -----
    Real                :: R0_force
    !?? pspace
    Integer             :: inf_dur,cont_dur,ill_dur,icu_dur
    Real                :: icu_per_day(8)=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,8.0/)
    Real                :: less_contagious
    Logical             :: immune_stop

!!!!!-----5. define reductions in social contacts -----
    Integer             :: R0change(2,20)
    Character*5         :: R0county(20)
    Logical             :: R0delay
    Integer             :: R0delay_days
    Character*6         :: R0delay_type
    Logical             :: endogenous_lockdown
    Real                :: lockdown_connect
    Integer             :: lockdown_threshold
    Integer             :: lockdown_days

!!!!!-----7.  Define whether transition probabilities should differ by age and sex
    Character*10        :: control_age_sex
    Character*10        :: seed_before_char,seed_temp
    Character*10,Allocatable :: seed_seq(:),seed_inf_cont_seq(:),seed_inf_ncont_seq(:)
    Character*10,Allocatable :: seed_d_seq(:)
    Integer             :: days

    Integer             :: n_direct,n_directv,n_directl,n_dist
    Integer             :: size_lhc
    Real,Allocatable    :: lhc(:,:)
    Character*10,Allocatable :: lhc_name(:)
!!!!!-----8. variables for do loop -------
    Real,Allocatable    :: icu_risk(:),surv_ill(:),surv_icu(:)
    Integer,Allocatable :: temp1(:),temp(:),targert_icu_risk_index(:)
    Character*3,Allocatable :: temp_char_mul(:)
    Character               :: temp_char_sig
    Type icu_risk_lists ! in order to have a similar structure to R
       Character*2,Allocatable         :: age(:)
       Character*1,Allocatable         :: sex(:)
       Real,Allocatable                :: risk(:)
       Integer,Allocatable             :: dur(:)
    End Type icu_risk_lists

    Type(icu_risk_lists)            :: icu_risk_list,surv_ill_list,surv_icu_list

    Type sims
       Integer,Allocatable             :: dist_id(:)
       Character,Allocatable           :: sex(:)
       Character*2,Allocatable         :: age(:)
       Integer,Allocatable             :: t1(:)
       Integer,Allocatable             :: t2(:)
       Integer,Allocatable             :: d(:)
    End Type sims

    Type(sims)                      :: sim,tmp
    Integer,Allocatable             :: tmp_dnew(:)
    Integer,Allocatable             :: sim_counties(:),rownumbers(:)

    Type(seeds)                     :: seed_ill,seed_inf_cont,seed_inf_ncont,target_inf
    Integer,Allocatable             :: seed_ill_dur(:),seed_inf_cont_dur(:),seed_inf_ncont_dur(:)
    Integer,Allocatable             :: il_d(:),inf_c_d(:),inf_nc_d(:)
    Integer,Allocatable             :: rownumbers_ill(:),rownumbers_cont(:),rownumbers_ncont(:),rownumbers_dea(:),&
         rownumbers_left(:)
    Integer,Allocatable             :: gettime(:)
    Real,Allocatable                :: getchange(:)
    Integer                         :: inf_ill,inf_cont,inf_ncont,inf_dth

    Real                            :: R0_daily
    Real,Allocatable                :: R0matrix(:,:),connect(:,:)
    Integer,Allocatable             :: healthy_cases_final(:,:,:),&
         ill_ICU_cases_final(:,:,:),immune_cases_final(:,:,:),&
         inf_noncon_cases_final(:,:,:),inf_contag_cases_final(:,:,:),&
         dead_cases_final(:,:,:),ill_contag_cases_final(:,:,:)
    Integer,Allocatable             :: inf_cases(:,:),icu_cases(:,:),healthy_cases(:,:),inf_noncon_cases(:,:),&
         inf_contag_cases(:,:),ill_contag_cases(:,:)
    Integer,Allocatable             :: ill_ICU_cases(:,:),immune_cases(:,:),dead_cases(:,:),dead_cases_bICU(:,:),&
         mod_inf_cases(:,:),org_noncon_cases(:,:)
    Integer,Allocatable             :: lockdowns(:),start_value_tot(:)
    Integer                         :: timestep

    Integer,Allocatable             :: tmp_index(:),tmp_d_new(:),tmp_count(:)
    Integer,Allocatable             :: susceptible(:),contagious_dist_id(:),contagious_index(:),denominator(:),&
         revers_proj(:),final_count(:),dist_id_temp(:),ill_index(:),&
         ill_dist_id(:)
    Real,Allocatable                :: contagious(:)
    Integer                         :: at_risk
    Integer                         :: initial_sick
    Real                            :: n_contagious,between_weight,within_weight
    Real,Allocatable                :: exp_infect(:)
    Integer,Allocatable             :: check_days(:)
    Logical                         :: flag_lockdown
    Real,Allocatable                :: risk(:),prob(:),runif(:),prop_target(:)
    Character*10                    :: target_date,temp_date

    Integer,Allocatable             :: sick(:),case_count(:),state_id(:),prop_inf_cases(:),&
         prop_target_inf(:),mod_inf(:)
    Character*6,Allocatable         :: age_sex(:),surv_ill_label(:),surv_icu_label(:)
    Character*6,Allocatable         :: age_sex_dur(:),icu_risk_label(:)
    Character*2,Allocatable         :: temp_character(:)
    Real,Allocatable                :: surv_ill_i(:),die(:),icu_risk_i(:),icu(:),surv_icu_i(:)
    Integer,Allocatable             :: die_count(:),icu_count(:),die_icu_count(:)
    Character*2,Allocatable         :: ch_age(:)
    Character*5,Allocatable         :: ch_sex(:)
    Integer                         :: max_date,n_change
    Character*10                    :: temp_mod, seed_date
    Integer                         :: ierror,size_of_process,my_rank
    Integer                         :: index_final(7),block_size
    Integer,Allocatable             :: req(:)

    Real                            :: iter_pass_handle(6)
    ! should import some reliable romdon seed generation code here
    !seed_base = ??

    index_final = 1
    Call MPI_INIT(ierror)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD,size_of_process,ierror)
    Call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierror)

    if ( size_of_process < 2 ) then
       write(*,*)"Sorry this application requires at least 2 MPI-Processes"
       write(*,*)"Programm halted !!!"
       goto 1000
    end if
    
    Allocate(req(size_of_process))

    seed_in_inner_loop = .False.
    seed_mpi_in_inner_loop = .True.

!!!!!-----1.states of the model ------

    !------state that the simulated individuals can reach
    healthy    = 0
    inf_noncon = 1
    inf_contag = 2
    ill_contag = 3
    ill_ICU    = 4
    immune     = 5
    dead       = 6

    missing     = -99

!!!!!-----2.derive the initial population ------
    sim_pop = "proportional"

!!!!!-----3.Set seed infections in population ------
    ini_infected = 10
    seed_infections = "data"
    seed_before = 7

!!!!!-----4.define disease characteristics------

    inf_dur  = 3
    cont_dur = 2
    ill_dur  = 8
    icu_dur  = 14
    ! icu_per_day = (/0,0,0,0,0,0,0,8/)
    If (Size(icu_per_day) /= ill_dur) Then
       Print *,'Length icu_per_day not equal to ill_dur'
    Endif
    If ((Sum(icu_per_day)/Size(icu_per_day)) /= 1) Then
       Print *,'Mean icu per day not equal to 1'
    End If

    less_contagious = 0.7

    R0_force = 0
    immune_stop = .True.

!!!!!-----5. define reductions in social contacts -----
    ! R0change = reshape((/1,7,8,13,14,19,20,25,26,32,33,39,40,46,47,53,54,60,61,67,68,74,75,81,82,88,&
    !                     89,95,96,102,103,109,110,116,117,123,124,130,131,137,138,144,145,151,152,158,&
    !                     159,164/),shape(R0change))
    R0change = Reshape((/1,7,8,13,14,19,20,25,26,32,33,39,40,46,47,53,54,60,61,67,68,74,75,81,82,88,&
         89,95,96,102,103,109,110,116,117,123,124,130,131,137/),Shape(R0change))
    R0county            = "ALL"
    R0delay             = .True.
    R0delay_days        = 5
    R0delay_type        = "linear"
    endogenous_lockdown = .False.
    lockdown_connect    = 0.5
    lockdown_threshold  = 100
    lockdown_days       = 10

!!!!!-----7.  Define whether transition probabilities should differ by age and sex
    control_age_sex     = "age"
    days             = 1
    seed_date        = add_date(sp%seed_date,days)
    days             = -1-seed_before
    seed_before_char = add_date(seed_date,days)
    seed_seq         = generate_seq(seed_before_char,seed_date)

    !Derive dates of infections for those that are inf_cont,
    !but are not yet aware about it (will be registered the
    !next two days)
    seed_temp        = add_date(seed_date,cont_dur)
    seed_inf_cont_seq  =  generate_seq(add_date(seed_date,1),seed_temp)
    !Derive dates of infections for those that are inf_cont,
    !but are not yet aware about it (will be registered the
    !next 3-5 days)
    seed_inf_ncont_seq = generate_seq(add_date(seed_date,cont_dur+1),add_date(seed_date,inf_dur+cont_dur))

    If (Maxval(R0change)> time_n) Then
       time_n = Maxval(R0change) + 1
    End If

    ! this part should done by the code, but here is set manually for simplicity
    n_direct  = 8
    n_directv = 0
    n_directl = 1
    n_dist    = 0

    If (n_dist > 0) Then
       ! code for randomLHS
       ! since it would not affect the code, it can be apllied later
    Else
       size_lhc = 0                    ! the first position is reserved for sam_size
    End If

    If (n_direct > 0) Then
       size_lhc = size_lhc + n_direct
    End If

    If (n_directl > 0) Then
       size_lhc = size_lhc + Size(pspace%ROeffect_ps%param)
    End If

    Allocate(lhc(size_lhc,iter))

    ! lhc(1,:)            = pspace%sam_size%param
    Do i = 1,n_direct
       lhc(i,:) = pspace%Ps_scalar_list(i)%param
    End Do


    Do i = 1,iter
       lhc(n_direct+1:Size(lhc,dim=1),i) = Reshape(pspace%ROeffect_ps%param,Shape(lhc(n_direct+1:Size(lhc,dim=1),1)))
    End Do
    ! print *, "after reshape is",reshape(pspace%ROeffect_ps%param,shape(lhc(n_direct+1:size(lhc,dim=1),1)))
    Allocate(lhc_name(size_lhc))

    ! do i = 1, 9
    !     lhc(i,:) = pspace%Ps_scalar_list(i)%param
    ! end do

    ! temp_int = 10
    !skip the following line,since lhc_name is no longer used.
    ! do i =1,size(iol%states_shortcut)
    !     do j = 1,size(pspace%ROeffect_ps%param,dim= 2)
    !         lhc_name(temp_int) = iol%states_shortcut // pspace%ROeffect_ps%param(1,j)
    !         temp_int = temp_int + 1
    !     end do
    ! end do

    If (control_age_sex == "NONE") Then
       ch_age = (/"total"/)
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
       ch_age = (/"total"/)
       ch_sex = (/"m","f"/)
    End If

    If (control_age_sex == "age_sex") Then
       !     ch_age = generate_seq(0,90,5)
       ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
            "40","45","50","55","60",&
            "65","70","75","80","85","90"/)
       ch_sex = (/"m","f"/)
    End If


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! the following part was in the loop at R code
!!! in order to aviod to raise problem of memery
!!! allocation and redundent calculation, put it
!!! outside the loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    temp = get_index(iol%transpr_age_gr,ch_age)
    temp1= get_index(iol%transpr_sex,ch_sex)
    targert_icu_risk_index = find_and(temp,temp1)

    If (control_age_sex == "age") Then
       If (.Not.Allocated(icu_risk))Then
          Allocate(icu_risk(2*Size(targert_icu_risk_index)))
       Endif
       icu_risk(1:Size(targert_icu_risk_index)) = iol%transpr_icu_risk(targert_icu_risk_index)
       icu_risk(1:Size(targert_icu_risk_index)) = 1.0 - (1.0 - icu_risk(1:Size(targert_icu_risk_index))) ** (1.0/ Real(ill_dur))
       icu_risk(Size(targert_icu_risk_index)+1 : 2* Size(targert_icu_risk_index)) = icu_risk(1:Size(targert_icu_risk_index))
    End If

    If (control_age_sex == "sex") Then
       If (.Not.Allocated(icu_risk))Then
          Allocate(icu_risk(2*19))
       End If
       icu_risk(1:2) = iol%transpr_icu_risk(targert_icu_risk_index)
       icu_risk(1:2) = 1.0 - (1.0 - icu_risk(1:2)) ** (1.0/ ill_dur)
       icu_risk(Size(targert_icu_risk_index)+1 : 2* Size(targert_icu_risk_index)) = icu_risk(2)
       icu_risk(1:Size(targert_icu_risk_index)) = icu_risk(1)
    End If
    !     skip line 960 to 940 in R code
    If (.Not.Allocated(icu_risk_list%age))Then
       Allocate(icu_risk_list%age(Size(icu_risk)*ill_dur))
       Allocate(icu_risk_list%sex(Size(icu_risk)*ill_dur))
       Allocate(icu_risk_list%risk(Size(icu_risk)*ill_dur))
       Allocate(icu_risk_list%dur(Size(icu_risk)*ill_dur))
    End If

    Do i = 1, ill_dur
       icu_risk_list%age((i-1)*Size(icu_risk)+1: i*Size(icu_risk)) = (/iol%transpr_age_gr(targert_icu_risk_index),&
            iol%transpr_age_gr(targert_icu_risk_index)/)
       icu_risk_list%sex((i-1)*Size(icu_risk)+1:(i-1)*Size(icu_risk)+19) = 'm'
       icu_risk_list%sex((i-1)*Size(icu_risk)+20:i*Size(icu_risk)) = 'f'
       icu_risk_list%risk((i-1)*Size(icu_risk)+1: i*Size(icu_risk)) = icu_risk * icu_per_day(i)
       icu_risk_list%dur((i-1)*Size(icu_risk)+1: i*Size(icu_risk)) = i
    End Do

    ! init surv_ill
    If (control_age_sex == "age") Then
       If (.Not.Allocated(surv_ill))Then
          Allocate(surv_ill(2*Size(targert_icu_risk_index)))
       End If
       surv_ill = (/iol%transpr_surv_ill(targert_icu_risk_index),&
            iol%transpr_surv_ill(targert_icu_risk_index)/)
       surv_ill = surv_ill ** (1.0/Real(ill_dur))
    End If

    If (control_age_sex == "sex") Then
       If (.Not.Allocated(surv_ill))Then
          Allocate(surv_ill(2*19))
       Endif
       surv_ill(1:2) = iol%transpr_surv_ill(targert_icu_risk_index)
       surv_ill(1:2) = 1.0 - (1.0 - surv_ill(1:2)) ** (1.0/ Real(ill_dur))
       surv_ill(Size(targert_icu_risk_index)+1 : 2* Size(targert_icu_risk_index)) = surv_ill(2)
       surv_ill(1:Size(targert_icu_risk_index)) = surv_ill(1)
    End If
    If (.Not.Allocated(surv_ill_list%age))Then
       Allocate(surv_ill_list%age(Size(surv_ill)))
       Allocate(surv_ill_list%sex(Size(surv_ill)))
       Allocate(surv_ill_list%risk(Size(surv_ill)))
    End If
    surv_ill_list%age = (/iol%transpr_age_gr(targert_icu_risk_index),iol%transpr_age_gr(targert_icu_risk_index)/)
    surv_ill_list%sex = 'f'
    surv_ill_list%sex(1:19) = 'm'
    surv_ill_list%risk  = surv_ill

    ! init surv_icu
    If (control_age_sex == "age") Then
       If (.Not.Allocated(surv_icu))Then
          Allocate(surv_icu(2*Size(targert_icu_risk_index)))
       End If
       surv_icu= (/iol%transpr_surv_icu(targert_icu_risk_index),&
            iol%transpr_surv_icu(targert_icu_risk_index)/)
       surv_icu = surv_icu ** (1/Real(ill_dur))
    End If

    If (control_age_sex == "sex") Then
       If (.Not.Allocated(surv_icu))Then
          Allocate(surv_icu(2*19))
       End If
       surv_icu(1:2) = iol%transpr_surv_icu(targert_icu_risk_index)
       surv_icu(1:2) = 1.0 - (1.0 - surv_icu(1:2)) ** (1.0/ Real(ill_dur))
       surv_icu(Size(targert_icu_risk_index)+1 : 2* Size(targert_icu_risk_index)) = surv_icu(2)
       surv_icu(1:Size(targert_icu_risk_index)) = surv_icu(1)
    End If
    If (.Not.Allocated(surv_icu_list%age))Then
       Allocate(surv_icu_list%age(Size(surv_icu)))
       Allocate(surv_icu_list%sex(Size(surv_icu)))
       Allocate(surv_icu_list%risk(Size(surv_icu)))
    End If
    surv_icu_list%age = (/iol%transpr_age_gr(targert_icu_risk_index),iol%transpr_age_gr(targert_icu_risk_index)/)
    surv_icu_list%sex = 'f'
    surv_icu_list%sex(1:19) = 'm'
    surv_icu_list%risk  = surv_icu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! the following part was in the loop at R code
!!! put it outside to avoid memoery problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (my_rank == 0)Then

       Allocate(healthy_cases_final(Size(counties_index),time_n,iter))
       Allocate(inf_noncon_cases_final(Size(counties_index),time_n,iter))
       Allocate(inf_contag_cases_final(Size(counties_index),time_n,iter))
       Allocate(ill_contag_cases_final(Size(counties_index),time_n,iter))
       Allocate(ill_ICU_cases_final(Size(counties_index),time_n,iter))
       Allocate(immune_cases_final(Size(counties_index),time_n,iter))
       Allocate(dead_cases_final(Size(counties_index),time_n,iter))

       Allocate(inf_cases(Size(counties_index),time_n))
       Allocate(icu_cases(Size(counties_index),time_n))

       Allocate(healthy_cases(Size(counties_index),time_n))
       Allocate(inf_noncon_cases(Size(counties_index),time_n))
       Allocate(inf_contag_cases(Size(counties_index),time_n))
       Allocate(ill_contag_cases(Size(counties_index),time_n))
       Allocate(ill_ICU_cases(Size(counties_index),time_n))
       Allocate(immune_cases(Size(counties_index),time_n))
       Allocate(dead_cases(Size(counties_index),time_n))
       Allocate(dead_cases_bICU(Size(counties_index),time_n))
       Allocate(mod_inf_cases(Size(counties_index),time_n))
       Allocate(org_noncon_cases(Size(counties_index),time_n))

    Else
       Allocate(inf_cases(Size(counties_index),time_n))
       Allocate(icu_cases(Size(counties_index),time_n))

       Allocate(healthy_cases(Size(counties_index),time_n))
       Allocate(inf_noncon_cases(Size(counties_index),time_n))
       Allocate(inf_contag_cases(Size(counties_index),time_n))
       Allocate(ill_contag_cases(Size(counties_index),time_n))
       Allocate(ill_ICU_cases(Size(counties_index),time_n))
       Allocate(immune_cases(Size(counties_index),time_n))
       Allocate(dead_cases(Size(counties_index),time_n))
       Allocate(dead_cases_bICU(Size(counties_index),time_n))
       Allocate(mod_inf_cases(Size(counties_index),time_n))
       Allocate(org_noncon_cases(Size(counties_index),time_n))

    Endif
    Allocate(lockdowns(time_n))
    block_size = Size(counties_index) * time_n

    max_date = find_max_date(iol%seed_date)

!!$OMP PARALLEL DO private(sim)
    If (my_rank/=0)Then
       Do it_ss = my_rank,Size(lhc,dim=2),size_of_process-1
          Call random_Seed()
!!!!!!!!!!!
          ! initalizing icu_risk and surv_icu
!!!!!!!!!!!

          !     allocate(icu_risk(size(targert_icu_risk_index)))
          !     icu_risk = iol%tran_pr_icu_risk(targert_icu_risk_index)

          !     icu_risk = 1 - (1 - icu_risk) ** (1/ ill_dur)
!!!!!!!!!!!
          !random Sampling
!!!!!!!!!!!

          !!comment this block out only for testing!
          If (sim_pop == "proportional") Then
             iol%pop_total = Nint(Real(iol%pop_total)/Real(Sum(iol%pop_total))* Real(lhc(1,it_ss)))
             If (.Not.Allocated(sim%dist_id))Then
                Allocate(sim%dist_id(Sum(iol%pop_total)))
                Allocate(sim%sex(Sum(iol%pop_total)))
                Allocate(sim%age(Sum(iol%pop_total)))
                Allocate(sim%t1(Sum(iol%pop_total)))
                Allocate(sim%t2(Sum(iol%pop_total)))
                Allocate(sim%d(Sum(iol%pop_total)))
             Endif

             index = 0 ! position index
             ! set all male for testing

             Do i = 1, Size(iol%pop_total)
                temp_int = iol%pop_total(i)
                sim%dist_id(index+1: index+temp_int) = iol%pop_distid(i)
                sim%sex(index+1: index+temp_int) = iol%pop_sex(i)
                sim%age(index+1: index+temp_int) = iol%pop_age(i)
                index                            = index + temp_int
             End Do

             sim%t1 = healthy
             sim%t2 = missing
             sim%d(:)  = 1
          End If
          If (my_rank == 1)Then
             Print *,"size of population is",Size(sim%d)
          End If

          If (sim_pop == "random") Then
             ! can be filled later
          End If

!!!!!!!!!!!!!
          !seed infections
!!!!!!!!!!!!!
          If (seed_infections == "random") Then
             ! remains blank
          End If

          If (seed_infections == "county") Then
             ! remains blank
          End If

          If (seed_infections == "data") Then
             ! skip the first two lines in R since it would not affect
             ! the code, code skip : line 1102-1104
             !         sim_counties = get_unique(sim%dist_id)

             temp = get_index(iol%death_distid,counties_index)

             iol%death_distid = iol%death_distid(temp)
             iol%death_date   = iol%death_date(temp)
             iol%death_cases  = iol%death_cases(temp)
             !skip line 1113, 1115

             temp = get_index(iol%seed_date,seed_seq)
             !         print *,size(temp)
             seed_ill%dist_id  = iol%seed_distid(temp)
             seed_ill%date     = iol%seed_date(temp)
             seed_ill%cases    = iol%seed_cases(temp)

             temp = get_index(iol%seed_date,seed_inf_cont_seq)
             seed_inf_cont%dist_id  = iol%seed_distid(temp)
             seed_inf_cont%date     = iol%seed_date(temp)
             seed_inf_cont%cases    = iol%seed_cases(temp)

             temp = get_index(iol%seed_date,seed_inf_ncont_seq)
             seed_inf_ncont%dist_id  = iol%seed_distid(temp)
             seed_inf_ncont%date     = iol%seed_date(temp)
             seed_inf_ncont%cases    = iol%seed_cases(temp)

             temp_date  =  get_start_date(iol%death_date)
             seed_d_seq = generate_seq(temp_date,add_date(seed_date,-1))

             temp       = get_index(iol%death_date,seed_d_seq)
             iol%death_distid = iol%death_distid(temp)
             iol%death_date   = iol%death_date(temp)
             iol%death_cases  = iol%death_cases(temp)

             If (.Not.Allocated(seed_ill_dur))Then
                Allocate(seed_ill_dur(Size(seed_ill%date)))
             End If
             temp_int = Date2Unixtime(seed_date)

             Do i = 1,Size(seed_ill_dur)
                seed_ill_dur(i)  = (temp_int - Date2Unixtime(seed_ill%date(i)))/86400 + 1
             End Do
             temp_int = 0
             !         print *,seed_ill%cases
             !         mod1 = "l"
             !         mod2 = "g"
             temp = condition_and(seed_ill_dur,ill_dur+1,"l",seed_ill%cases,temp_int,"g")
             !          print *,"size of temp",size(temp)
             seed_ill%dist_id  = seed_ill%dist_id(temp)
             seed_ill%date     = seed_ill%date(temp)
             seed_ill%cases    = seed_ill%cases(temp)

             If (.Not.Allocated(seed_inf_cont_dur))Then
                Allocate(seed_inf_cont_dur(Size(seed_inf_cont%date)))
             End If
             temp_int = Date2Unixtime(seed_date)
             Do i = 1,Size(seed_inf_cont_dur)
                seed_inf_cont_dur(i)  = (temp_int - Date2Unixtime(seed_inf_cont%date(i)))/86400&
                     + cont_dur + 1
             End Do

             Deallocate(temp)
             Allocate(temp(Size(seed_inf_cont%cases)))

             temp = 0
             temp_int = 0
             Do i = 1,Size(temp)
                If (seed_inf_cont%cases(i)>0) Then
                   temp_int = temp_int + 1
                   temp(temp_int) = i
                End If
             End Do

             !         seed_inf_cont = set_value(seed_inf_cont,temp(1:temp_int))
             seed_inf_cont%dist_id  = seed_inf_cont%dist_id(temp(1:temp_int))
             seed_inf_cont%date     = seed_inf_cont%date(temp(1:temp_int))
             seed_inf_cont%cases    = seed_inf_cont%cases(temp(1:temp_int))
             Deallocate(temp)

             If (.Not.Allocated(seed_inf_ncont_dur))Then
                Allocate(seed_inf_ncont_dur(Size(seed_inf_ncont%date)))
             Else
                Deallocate(seed_inf_ncont_dur)
                Allocate(seed_inf_ncont_dur(Size(seed_inf_ncont%date)))
             End If
             temp_int = Date2Unixtime(seed_date)
             Do i = 1,Size(seed_inf_ncont_dur)
                seed_inf_ncont_dur(i)  = (temp_int - Date2Unixtime(seed_inf_ncont%date(i)))/86400&
                     + inf_dur + cont_dur + 1
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
             !         seed_inf_ncont = set_value(seed_inf_ncont,temp(1:temp_int))
             seed_inf_ncont%dist_id  = seed_inf_ncont%dist_id(temp(1:temp_int))
             seed_inf_ncont%date     = seed_inf_ncont%date(temp(1:temp_int))
             seed_inf_ncont%cases    = seed_inf_ncont%cases(temp(1:temp_int))
             !       skip line 1163,1166,1169,1172 since the mechanism of aggregate is not clear, and it seems not changing anything at all
             !
             !       do scaling
             ! skip scaling ,since they are all 1
             !         seed_ill%cases = seed_ill%cases * sam_prop_ps(seed_ill%dist_id/1000)
             !
             !         seed_inf_cont%cases = seed_inf_cont%cases * sam_prop_ps(seed_ill%dist_id/1000)
             !
             !         seed_inf_ncont%cases = seed_inf_ncont%cases * sam_prop_ps(seed_ill%dist_id/1000)

             Do icounty = 1,Size(counties_index)

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

                temp   = get_index(iol%death_distid,county)
                inf_dth = Sum(iol%death_cases(temp))

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
                   rownumbers_left = rownumbers_left(inf_ill+1:Size(rownumbers))
                   sim%t1(rownumbers_ill) = ill_contag
                   sim%d(rownumbers_ill)  = il_d
                End If

                If ( inf_cont > 0) Then
                   rownumbers_cont = sample(rownumbers_left,inf_cont)
                   rownumbers_left = rownumbers_left(inf_cont+1:Size(rownumbers_left))
                   sim%t1(rownumbers_cont)= inf_contag
                   sim%d(rownumbers_cont) = inf_c_d
                End If

                If (inf_ncont > 0) Then
                   rownumbers_ncont = sample(rownumbers_left,inf_ncont)
                   rownumbers_left = rownumbers_left(inf_ncont+1:Size(rownumbers_left))
                   sim%t1(rownumbers_ncont) = inf_noncon
                   sim%d(rownumbers_ncont)  = inf_nc_d
                End If

                If (inf_dth > 0) Then
                   rownumbers_dea = sample(rownumbers_left,inf_dth)
                   sim%t1(rownumbers_dea) = dead
                End If
             End Do ! do icounty = 1,size(sim_counties)
          End If ! if (seed_infections == "data")

          If (import_R0_matrix) Then
             !         temp = get_index(iol%R0_raw)
             !         R0_raw1 =
          End If
          R0_daily = R0_force *lhc(2,it_ss)/Real(Real(cont_dur)+Real(ill_dur)*less_contagious) + &
               (1-R0_force)*lhc(2,it_ss)/Real(cont_dur+ill_dur)
          ! this block simplifies the if judgment
          If (.Not.Allocated(R0matrix))Then
             Allocate(R0matrix(Size(counties_index),time_n-1))
          End If
          R0matrix = R0_daily
          n_change  = Size(R0change,dim=2)

          Do i = 1,n_change

             If(R0county(i) == "ALL") Then
                !use counties
                ! give up using character to find the number, but use the index directly
                ! if is "all" use full counties_index
                temp  = (i-1) * Size(pspace%ROeffect_ps%param,dim= 1) + 9 + ((counties_index/1000)-1)

                getchange = lhc(temp,it_ss)

                gettime = generate_seq(R0change(1,i),R0change(2,i),1)
                Do j = 1,Size(gettime)
                   R0matrix(:,gettime(j)) = R0matrix(:,gettime(j)) * getchange
                End Do
             End If
          End Do

          If (R0delay) Then
             Do i = 1,Size(R0matrix,dim = 1)
                R0matrix(i,:) = smoothing_change(R0matrix(i,:),R0delay_days,R0delay_type)
             End Do
          End If

          !     print *,"R0delay_days",R0delay_days
          !     print *,"R0delay_type",R0delay_type



          start_value_tot        = sum_bygroup(sim%t1,sim%dist_id,counties_index,"ill")
          inf_cases(:,1)         = start_value_tot
          icu_cases(:,1)         = 0
          temp_mod = "healthy   "
          healthy_cases(:,1)     = sum_bygroup(sim%t1,sim%dist_id,counties_index,temp_mod)
          temp_mod = "inf_ncon  "
          inf_noncon_cases(:,1)  = sum_bygroup(sim%t1,sim%dist_id,counties_index,temp_mod)
          temp_mod = "inf_con   "
          inf_contag_cases(:,1)  = sum_bygroup(sim%t1,sim%dist_id,counties_index,temp_mod)
          temp_mod = "ill_con   "
          ill_contag_cases(:,1)  = sum_bygroup(sim%t1,sim%dist_id,counties_index,temp_mod)
          temp_mod = "dead      "
          dead_cases(:,1)        = sum_bygroup(sim%t1,sim%dist_id,counties_index,temp_mod)
          ill_ICU_cases(:,1)     = 0
          immune_cases(:,1)      = 0
          dead_cases_bICU(:,1)   = 0
          mod_inf_cases(:,1)     = 0
          org_noncon_cases(:,1)  = 0
          sim%t2                 = sim%t1
          If (.Not.Allocated(tmp_d_new))Then
             Allocate(tmp_d_new(Size(sim%d)))
          End If

          Do timestep = 2,time_n
             !skip using prev_state and up_state
             tmp = sim

             !         print *,"now in the ",timestep, " th loop"

             tmp_index = get_index(tmp%t1,(/inf_contag,inf_noncon,ill_ICU,ill_contag,healthy/))

             tmp_d_new = missing

             If (Size(tmp_index) == 0)Then
                ! leave it blank and deal with it later
             End If
             tmp_count = get_index(tmp%t1,healthy)
             !         susceptible = sum_byindex(tmp%dist_id(tmp_count),counties_index)

             at_risk = Size(tmp_count)

             contagious_index   = get_index(tmp%t1,inf_contag)
             ill_index          = get_index(tmp%t1,ill_contag)

             contagious_dist_id = get_unique(tmp%dist_id(contagious_index))
             !         if(my_rank == 1)then
             !              print *,size(contagious_dist_id)
             !         endif
             !         print *,"ill case is",size(ill_index)
             ill_dist_id        = get_unique(tmp%dist_id(ill_index))
             If (.Not.Allocated(contagious))Then
                Allocate(contagious(Size(counties_index)))
             End If
             !         do i = 1,size(counties_index)
             !             temp1 = get_index(tmp%dist_id(contagious_index),counties_index(i))
             !             contagious(i) = real(size(temp1))
             !             temp1 = get_index(tmp%dist_id(ill_index),counties_index(i))
             !             contagious(i) = contagious(i) + real(size(temp1)) * less_contagious
             !         end do
             contagious = Real(sum_byindex(tmp%dist_id(contagious_index),counties_index))
             !         if(my_rank == 1)then
             !              print *,size(contagious)
             !         endif
             contagious = contagious + Real(sum_byindex(tmp%dist_id(ill_index),counties_index)) * less_contagious
             !         print *,"ill cases are: ",real(size(temp1))
             !         print *,"nint of nint(sum(contagious)) is",nint(sum(contagious))
             !         print *,"sum(susceptible) is ",sum(susceptible)
             n_contagious = Sum(contagious)
             temp = get_index(tmp%t1,inf_noncon)
             !         if(my_rank == 1)then
             !             print *,at_risk
             !         endif
             If (at_risk > 0 .And. n_contagious > 0) Then
                exp_infect = contagious * R0matrix(:,timestep-1)
                !             if(my_rank == 1)then
                !                 print *,at_risk
                !             endif
                !             print *,"contagious is",contagious
                If (endogenous_lockdown) Then
                   temp = generate_seq(timestep-1,timestep-lockdown_days,-1)
                   check_days = Pack(temp,temp>0)
                   temp = Sum(inf_cases(:,check_days),dim = 1)
                   flag_lockdown  = .False.
                   temp = Pack(temp,temp>lockdown_threshold)
                   If (Size(temp)>0)Then
                      flag_lockdown = .True.
                   End If
                   If (flag_lockdown)Then
                      exp_infect = exp_infect *lhc(5,it_ss)! times lock_downeffect, which is the 5th element of the array
                      lockdowns(timestep) = 1
                   Else
                      lockdowns(timestep) = 0
                   End If

                Else
                   lockdowns(timestep) = 0
                End If ! end if (endogenous_lockdown)
                temp = get_index(iol%connect_work_distid,counties_index)

                connect = iol%connect_work(temp,temp)
                Do i = 1,Size(contagious_dist_id)
                   connect(:,i) = connect(:,i)/Sum(connect(:,i))
                End Do

                !             connect = transpose(connect)
                between_weight = 1 - lhc(6,it_ss) ! 1 - lhc(it_ss,"w_int")

                If (endogenous_lockdown)Then
                   between_weight = between_weight * lockdown_connect
                End If

                within_weight = 1 - between_weight
                exp_infect = within_weight * exp_infect + between_weight * Matmul (exp_infect, connect)!matmul (connect,exp_infect)
                !             do i = 1,size(contagious_dist_id)
                !                 exp_infect(i) = within_weight * exp_infect(i) + between_weight * sum(exp_infect*connect(:,i))
                !             end do
                If (.Not.Allocated(denominator))Then
                   Allocate(denominator(Size(counties_index)))
                End If
                If (.Not.immune_stop) Then
                   !                 denominator = sum_bygroup(tmp%dist_id(temp))
                   !                 do i = 1,size(counties_index)
                   !                     temp1 = get_index(tmp%dist_id(temp),counties_index(i))
                   !                     denominator(i) = size(temp1)
                   !                 end do
                   denominator = sum_byindex(tmp%dist_id(tmp_count),counties_index)
                Else
                   temp = get_index(tmp%t1,(/healthy,immune,inf_noncon,inf_contag,ill_contag/))
                   !                 do i = 1,size(counties_index)
                   !                     temp1 = get_index(tmp%dist_id(temp),counties_index(i))
                   !                     denominator(i) = size(temp1)
                   !                 end do
                   denominator = sum_byindex(tmp%dist_id(temp),counties_index)
                End If

                risk = exp_infect/Real(denominator)
                Where (risk > 1000)
                   risk = 0
                End Where
                Where (risk > 1)
                   risk = 1
                End Where

                risk = lhc(6,it_ss)*risk + (1-lhc(6,it_ss))* Matmul(risk,connect)!matmul (connect,risk)


                temp = get_index(counties_index,tmp%dist_id(tmp_count))
                prob = risk(temp)

                If (.Not.Allocated(runif))Then
                   Allocate(runif(at_risk))
                   Allocate(sick(at_risk))
                Else
                   Deallocate(runif)
                   Deallocate(sick)
                   Allocate(runif(at_risk))
                   Allocate(sick(at_risk))
                End If
                Call random_Number(runif)
                !check if probs in runif is higher than
                temp_int = 0
                sick = 0
                Do i = 1,Size(runif)
                   If (runif(i) <= prob(i)) Then
                      temp_int = temp_int + 1
                      sick(i)  = 1
                   End If
                End Do

                If (Allocated(revers_proj))Then
                   Deallocate(revers_proj)
                   Allocate(revers_proj(temp_int))
                Else
                   Allocate(revers_proj(temp_int))
                Endif
                temp_int = 0

                Do i = 1,Size(runif)
                   If (runif(i) <= prob(i)) Then
                      temp_int = temp_int + 1
                      revers_proj(temp_int) = tmp_count(i) ! record index of infected person in healthy group
                   End If
                End Do

                initial_sick = Size(sick)

                If (.Not.Allocated(final_count))Then
                   Allocate(final_count(Size(counties_index)))
                   !                 deallocate(dist_id_temp)
                Else
                   Deallocate(final_count)
                   Allocate(final_count(Size(counties_index)))
                End If
                final_count = 0
                !             call sum_bygroup_distID(tmp%dist_id(revers_proj),final_count,dist_id_temp)
                !             do i =1,size(counties_index)
                !                 temp = get_index(tmp%dist_id(revers_proj),counties_index(i))
                !                 final_count(i) = size(temp)
                !             enddo
                final_count = sum_byindex(tmp%dist_id(revers_proj),counties_index)
                inf_cases(:,timestep) = final_count
                org_noncon_cases(:,timestep) = final_count

                target_date = add_date(seed_date,timestep-1+7)

                If (lhc(7,it_ss) > 0 .And. Date2Unixtime(target_date) < max_date .And. initial_sick > 0)Then
                   !                print *,"i'm in"
                   temp = get_index(iol%seed_date,target_date)
                   target_inf%cases   = iol%seed_cases(temp)
                   target_inf%dist_id = iol%seed_distid(temp)
                   target_inf%date    = iol%seed_date(temp)
                   If (.Not.Allocated(final_count))Then
                      Allocate(final_count(Size(counties_index)))
                   Else
                      Deallocate(final_count)
                      Allocate(final_count(Size(counties_index)))
                   End If
                   final_count = 0
                   !                 do i =1,size(counties_index)
                   !                     temp = get_index(target_inf%dist_id,counties_index(i))
                   !                     final_count(i) = size(temp)
                   !                 end do
                   final_count = sum_byindex(target_inf%dist_id,counties_index)
                   If (Sum(final_count) > 0) Then
                      If(lhc(8,it_ss) == 1) Then
                         state_id = get_unique(counties_index/1000)
                         If (.Not.Allocated(mod_inf))Then
                            Allocate(mod_inf(Size(counties_index)))
                         End If
                         mod_inf = 0
                         Do i = 1, Size(state_id)
                            temp = get_index(counties_index/1000,state_id(i))
                            prop_inf_cases = inf_cases(temp,timestep)
                            temp_int       = Sum(prop_inf_cases)
                            If (temp_int > 0) Then
                               prop_inf_cases = prop_inf_cases / temp_int
                            Else
                               prop_inf_cases = 0
                            End If

                            prop_target_inf = final_count(temp)
                            temp_int = Sum(prop_target_inf)
                            If (temp_int > 0) Then
                               prop_target_inf = prop_target_inf / temp_int
                            Else
                               prop_target_inf = 0
                            End If
                            prop_target = lhc(7,it_ss) * Real(prop_target_inf) + (1 - lhc(7,it_ss))*Real(prop_inf_cases)

                            mod_inf(temp) = Nint(prop_target* Real(temp_int)) - inf_cases(temp,timestep)
                         End Do

                      Else
                         prop_inf_cases = inf_cases(:,timestep) / initial_sick
                         prop_target = lhc(7,it_ss) * Real(final_count)/Real(Sum(final_count)) + &
                              (1 - lhc(7,it_ss))*Real(prop_inf_cases)
                         ! line 1675
                         mod_inf = Nint(prop_target * Real(initial_sick))-inf_cases(:,timestep)

                      End If

                      Call shift_cases(sick,mod_inf,tmp%dist_id(tmp_count),counties_index)
!!!!!!!!!!!!!!!!!!!!!!!!
!!!can we write the following code in a more elgant
!!!way?
!!!!!!!!!!!!!!!!!!!!!!!!

                      temp = get_index(sick,1)
                      ! do a revers projection
                      temp = tmp_count(temp)
                      If (.Not.Allocated(final_count))Then
                         Allocate(final_count(Size(counties_index)))
                      End If
                      If (Allocated(dist_id_temp))Then
                         Deallocate(dist_id_temp)
                      End If
                      If (Allocated(case_count))Then
                         Deallocate(case_count)
                      End If

                      !                     call sum_bygroup_distID(tmp%dist_id(temp),case_count,dist_id_temp)
                      final_count = 0
                      !                     do i =1,size(counties_index)
                      !                         temp1 = get_index(tmp%dist_id(temp),counties_index(i))
                      !                         final_count(i) = size(temp1)
                      !                     enddo
                      final_count = sum_byindex(tmp%dist_id(temp),counties_index)
                      inf_cases(:,timestep) = final_count

                      mod_inf_cases(:,timestep) = mod_inf
                   Else
                      Print *,"no cases shifted"
                   End If

                Else
                   mod_inf_cases(:,timestep) = 0
                Endif

                !skip 1742 and 1743 since they are not actually do
                !anythin?

                tmp%t2(tmp_count) = sick

                temp = condition_and(tmp%t1,healthy,"e",tmp%t2,inf_noncon,"e")


                tmp_d_new(temp) = 1
                !             tmp%d(temp) = 1

                temp = condition_and(tmp%t1,healthy,"e",tmp%t2,healthy,"e")
                tmp_d_new(temp) = tmp%d(temp) + 1
                !             tmp%d(temp) =    tmp%d(temp) + 1
             End If

             If (at_risk > 0 .And. Nint(n_contagious) == 0) Then
                !             temp = get_index(tmp%t1,healthy)
                tmp%t2(tmp_count) = healthy
                !             if (.not.allocated(tmp_d_new))then
                !                 allocate(tmp_d_new(size(tmp%d)))
                !             else
                !                 deallocate(tmp_d_new)
                !                 allocate(tmp_d_new(size(tmp%d)))
                !             endif

                temp = condition_and(tmp%t1,healthy,"e",tmp%t2,healthy,"e")

                tmp_d_new(temp) = tmp%d(temp) + 1
                !             tmp%d(temp)       =     tmp%d(temp) + 1
                inf_cases(:,timestep) = 0

             End If

             If (at_risk == 0) Then
                inf_cases(:,timestep) = 0
             End If


             !If day limit not reach: Stays the same
             temp = condition_and(tmp%t1,inf_noncon,"e",tmp%d,inf_dur,"l")
             tmp%t2(temp) = inf_noncon
             !         print *,"number of infnoncon",size(temp)
             !Update day(current count plus 1)
             temp = condition_and(tmp%t1,inf_noncon,"e",tmp%t2,inf_noncon,"e")
             !         tmp%d(temp)= tmp%d(temp) + 1
             tmp_d_new(temp) = tmp%d(temp) + 1
             !         print *,"number of infnoncon update",size(temp)
             ! if day limit reached: move to contagious
             temp = condition_and(tmp%t1,inf_noncon,"e",tmp%d,inf_dur-1,"g")
             tmp%t2(temp) = inf_contag
             !         print *,"number of infcontag",size(temp)
             !Update days (Reset counter to 1 if moved)
             temp = condition_and(tmp%t1,inf_noncon,"e",tmp%t2,inf_contag,"e")
             !         tmp%d(temp) = 1
             tmp_d_new(temp) = 1
             !         print *,"number of infcontag update",size(temp)
!!!state infected,contagious
             !if day limit noe reach: stay the same
             temp = condition_and(tmp%t1,inf_contag,"e",tmp%d,cont_dur,"l")
             tmp%t2(temp) = inf_contag
             !         print *,"number of infcontag stay",size(temp)
             !update days
             temp = condition_and(tmp%t1,inf_contag,"e",tmp%t2,inf_contag,"e")
             !         tmp%d(temp) = tmp%d(temp) + 1
             tmp_d_new(temp) = tmp%d(temp) + 1
             !         print *,"number of infcontag stay update",size(temp)
             ! if day limit reached: move to contagious
             temp = condition_and(tmp%t1,inf_contag,"e",tmp%d,cont_dur-1,"g")
             !         print *,sim%d
             tmp%t2(temp) = ill_contag
             !         print *,"number of illcontag",size(temp)
             !update days(reset to 1 if becoming ill/contagious)
             temp = condition_and(tmp%t1,inf_contag,"e",tmp%t2,ill_contag,"e")
             !         tmp%d(temp) = 1
             tmp_d_new(temp) = 1
             !         print *,"number of illcontag update",size(temp)

!!!state:ill,contagious
             !dying population at risk
             temp = get_index(tmp%t1,ill_contag)
             at_risk = Size(temp)

             If (at_risk > 0) Then
                age_sex = tmp%age(temp) // tmp%sex(temp)
                surv_ill_label = surv_ill_list%age //surv_ill_list%sex
                temp = get_index(surv_ill_label,age_sex)
                surv_ill_i = surv_ill_list%risk(temp)
                If (.Not.Allocated(die)) Then
                   Allocate(die(Size(temp)))
                   Allocate(die_count(Size(temp)))
                Else
                   Deallocate(die)
                   Deallocate(die_count)
                   Allocate(die(Size(temp)))
                   Allocate(die_count(Size(temp)))
                End If
                die_count = ill_contag
                Call random_Number(die)
                Do i = 1,Size(die)
                   If (surv_ill_i(i) <= die(i))Then
                      die_count(i) = dead
                   Endif
                End Do
                temp = get_index(tmp%t1,ill_contag)

                tmp%t2(temp) = die_count

             Endif

             ! Save dead before moving to ICU
             temp1 = get_index(tmp%t2,dead)
             !         if (allocated(dist_id_temp))then
             !             deallocate(dist_id_temp)
             !         end if
             !         if (allocated(final_count))then
             !             deallocate(final_count)
             !         end if
             !         call sum_bygroup_distID(tmp%dist_id(temp1),final_count,dist_id_temp)
             final_count = 0
             !         do i = 1,size(counties_index)
             !             temp = get_index(tmp%dist_id(temp1),counties_index(i))
             !             final_count(i) = size(temp)
             !         end do

             final_count = sum_byindex(tmp%dist_id(temp1),counties_index)
             dead_cases_bICU(:,timestep) = final_count


             !moving to icu: population at risk
             temp = condition_and(tmp%t1,ill_contag,"e",tmp%t2,dead,"n")

             at_risk = Size(temp)



             If(at_risk > 0) Then
                If (.Not.Allocated(temp_character))Then
                   Allocate(temp_character(Size(temp)))
                Else
                   Deallocate(temp_character)
                   Allocate(temp_character(Size(temp)))
                End If
                If (.Not.Allocated(age_sex_dur))Then
                   Allocate(age_sex_dur(Size(temp)))
                Else
                   Deallocate(age_sex_dur)
                   Allocate(age_sex_dur(Size(temp)))
                End If
                Write(temp_character,'(I2)')tmp%d(temp)
                Do i =1,Size(temp)
                   age_sex_dur(i) = Trim(tmp%age(temp(i))) // Trim(tmp%sex(temp(i)))//Trim(temp_character(i))
                End Do

                If (.Not.Allocated(temp_character))Then
                   Allocate(temp_character(Size(icu_risk_list%dur)))
                Else
                   Deallocate(temp_character)
                   Allocate(temp_character(Size(icu_risk_list%dur)))
                End If
                Write(temp_character,'(I2)')icu_risk_list%dur
                If (.Not.Allocated(icu_risk_label))Then
                   Allocate(icu_risk_label(Size(icu_risk_list%age)))
                End If
                Do i = 1,Size(icu_risk_list%age)
                   icu_risk_label(i) = Trim(icu_risk_list%age(i))//Trim(icu_risk_list%sex(i))&
                        //Trim(temp_character(i))
                End Do

                temp1 = get_index(icu_risk_label,age_sex_dur)


                icu_risk_i = icu_risk_list%risk(temp1)

                If (.Not.Allocated(icu))Then
                   Allocate(icu(Size(temp)))
                   Allocate(icu_count(Size(temp)))
                Else
                   Deallocate(icu)
                   Deallocate(icu_count)
                   Allocate(icu(Size(temp)))
                   Allocate(icu_count(Size(temp)))
                Endif
                icu_count = ill_contag
                Call random_Number(icu)
                Do i = 1,Size(icu_risk_i)
                   If(icu(i) <= icu_risk_i(i))Then
                      icu_count(i) = ill_ICU
                   Endif
                Enddo

                tmp%t2(temp) = icu_count
                temp = condition_and(tmp%t1,ill_contag,"e",tmp%t2,ill_ICU,"e")
                !             tmp%d(temp)  = 1
                tmp_d_new(temp) = 1
             End If

             If (at_risk == 0)Then
                icu_cases(:,timestep) = 0
             Endif

             ! staying ill/contagious
             temp = condition_and(tmp%t1,ill_contag,"e",tmp%d,ill_dur,"l")
             temp1= get_index(tmp%t2,(/healthy,immune,inf_noncon,inf_contag,ill_contag/))
             temp = find_and(temp,temp1)

             tmp%t2(temp) = ill_contag

             ! update days
             temp = condition_and(tmp%t1,ill_contag,"e",tmp%t2,ill_contag,"e")
             !         tmp%d(temp) = tmp%d(temp) + 1
             tmp_d_new(temp) = tmp%d(temp) + 1

             !getting healthy/immune if day limit reached
             temp = condition_and(tmp%t1,ill_contag,"e",tmp%d,ill_dur-1,"g")! -1 because in R is >=
             temp1= get_index(tmp%t2,(/healthy,immune,inf_noncon,inf_contag,ill_contag/))
             temp = find_and(temp,temp1)
             tmp%t2(temp) = immune

             !update days
             temp = condition_and(tmp%t1,ill_contag,"e",tmp%t2,immune,"e")
             !         tmp%d(temp) = 1
             tmp_d_new(temp) = 1

!!!state ICU
             temp = get_index(tmp%t1,ill_ICU)
             at_risk = Size(temp)

             If(at_risk > 0) Then
                age_sex = tmp%age(temp) // tmp%sex(temp)
                surv_icu_label = surv_icu_list%age // surv_icu_list%sex
                temp1 = get_index(surv_icu_label,age_sex)

                surv_icu_i = surv_icu_list%risk(temp1)
                !who dies
                If (.Not.Allocated(die_icu_count))Then
                   Allocate(die_icu_count(at_risk))
                Else
                   Deallocate(die_icu_count)
                   Allocate(die_icu_count(at_risk))
                Endif
                die_icu_count = ill_contag
                Call random_Number(die)
                Do i = 1,Size(surv_icu_i)
                   If (die(i) >= surv_icu_i(i))Then
                      die_icu_count(i) = dead
                   End If
                End Do
                tmp%t2(temp) = die_icu_count
             End If

             !         tmp_count = get_index(tmp%t2,healthy)

             ! stay ill/ICU
             temp = condition_and(tmp%t1,ill_ICU,"e",tmp%d,Int(lhc(3,it_ss)),"l")
             temp1= get_index(tmp%t2,(/healthy,ill_ICU,inf_noncon,inf_contag,ill_contag,immune/))
             temp = find_and(temp,temp1)
             tmp%t2(temp) = ill_ICU

             ! update days
             temp = condition_and(tmp%t1,ill_ICU,"e",tmp%t2,ill_ICU,"e")
             !         tmp%d(temp) = tmp%d(temp) + 1
             tmp_d_new(temp) = tmp%d(temp) + 1

             !getting healthy /ICU if day limited reached
             temp = condition_and(tmp%t1,ill_ICU,"e",tmp%d,Int(lhc(3,it_ss))-1,"g")
             temp1= get_index_not(tmp%t2,dead)
             temp = find_and(temp,temp1)
             tmp%t2(temp) = immune

             !update days
             temp = condition_and(tmp%t1,ill_ICU,"e",tmp%t2,immune,"e")
             !         tmp%d(temp) = 1
             tmp_d_new(temp) = 1

             !check NA set to missing
             !skip this line 1995

             !move to main data frame
             temp = get_index(sim%t1,(/healthy,ill_ICU,inf_noncon,inf_contag,ill_contag/))
             sim%t2(temp) = tmp%t2(temp)
             !         sim%d(temp)  = tmp%d(temp)
             !         sim%d(temp)  = tmp_d_new(temp)
             !         sim%t2 = tmp%t2
             sim%d  = tmp_d_new
             !immune and dead remain the same
             temp = get_index(sim%t1,immune)
             sim%t2(temp) = immune
             temp = get_index(sim%t1,dead)
             sim%t2(temp) = dead




!!!calculate total daily numbers by counties
             ! 0 /healthy
             temp = get_index(sim%t2,healthy)
             final_count = 0
             !         do i = 1,size(counties_index)
             !             temp1 = get_index(tmp%dist_id(temp),counties_index(i))
             !             final_count(i) = size(temp)
             !         end do
             final_count = sum_byindex(tmp%dist_id(temp),counties_index)
             !         call sum_bygroup_distID(sim%dist_id(temp),final_count,dist_id_temp)
             !         temp1 = get_index(counties_index,dist_id_temp)
             healthy_cases(:,timestep) = final_count

             ! 1/infected
             temp = get_index(sim%t2,inf_noncon)
             final_count = 0
             !         do i = 1,size(counties_index)
             !             temp1 = get_index(tmp%dist_id(temp),counties_index(i))
             !             final_count(i) = size(temp1)
             !         end do
             final_count = sum_byindex(tmp%dist_id(temp),counties_index)
             !         call sum_bygroup_distID(sim%dist_id(temp),final_count,dist_id_temp)
             !         temp1 = get_index(counties_index,dist_id_temp)
             inf_noncon_cases(:,timestep) = final_count

             ! 2/infected,contagious
             temp = get_index(sim%t2,inf_contag)
             final_count = 0
             !         do i = 1,size(counties_index)
             !             temp1 = get_index(tmp%dist_id(temp),counties_index(i))
             !             final_count(i) = size(temp1)
             !         end do
             final_count = sum_byindex(tmp%dist_id(temp),counties_index)
             !         call sum_bygroup_distID(sim%dist_id(temp),final_count,dist_id_temp)
             !         temp1 = get_index(counties_index,dist_id_temp)
             inf_contag_cases(:,timestep) = final_count

             ! 3/ill,contagious
             temp = get_index(sim%t2,ill_contag)
             final_count = 0
             !         do i = 1,size(counties_index)
             !             temp1 = get_index(tmp%dist_id(temp),counties_index(i))
             !             final_count(i) = size(temp1)
             !         end do
             final_count = sum_byindex(tmp%dist_id(temp),counties_index)
             !         call sum_bygroup_distID(sim%dist_id(temp),final_count,dist_id_temp)
             !         temp1 = get_index(counties_index,dist_id_temp)
             ill_contag_cases(:,timestep) = final_count

             ! 4/ill,icu
             temp = get_index(sim%t2,ill_ICU)
             final_count = 0
             !         do i = 1,size(counties_index)
             !             temp1 = get_index(tmp%dist_id(temp),counties_index(i))
             !             final_count(i) = size(temp1)
             !         end do
             final_count = sum_byindex(tmp%dist_id(temp),counties_index)
             !         call sum_bygroup_distID(sim%dist_id(temp),final_count,dist_id_temp)
             !         temp1 = get_index(counties_index,dist_id_temp)
             ill_ICU_cases(:,timestep) = final_count

             ! 5/immune
             temp = get_index(sim%t2,immune)
             final_count = 0
             !         do i = 1,size(counties_index)
             !             temp1 = get_index(tmp%dist_id(temp),counties_index(i))
             !             final_count(i) = size(temp1)
             !         end do
             final_count = sum_byindex(tmp%dist_id(temp),counties_index)
             !         call sum_bygroup_distID(sim%dist_id(temp),final_count,dist_id_temp)
             !         temp1 = get_index(counties_index,dist_id_temp)
             immune_cases(:,timestep) = final_count

             ! 6/dead
             temp = get_index(sim%t2,dead)
             final_count = 0
             !         do i = 1,size(counties_index)
             !             temp1 = get_index(tmp%dist_id(temp),counties_index(i))
             !             final_count(i) = size(temp1)
             !         end do
             final_count = sum_byindex(tmp%dist_id(temp),counties_index)
             !         call sum_bygroup_distID(sim%dist_id(temp),final_count,dist_id_temp)
             !         temp1 = get_index(counties_index,dist_id_temp)
             dead_cases(:,timestep) = final_count

             sim%t1 = sim%t2




          End Do ! end do timestep =2,time_n
          !   if (my_rank /= 0)then   
          !        call MPI_Send(healthy_cases(1,1),block_size,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierror)
          !        print *,"ierror in branch",ierror
          !        call MPI_Send(inf_noncon_cases(1,1),block_size,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierror)
          !        call MPI_Send(inf_contag_cases(1,1),block_size,MPI_INTEGER,0,2,MPI_COMM_WORLD,ierror)
          !        call MPI_Send(ill_contag_cases(1,1),block_size,MPI_INTEGER,0,3,MPI_COMM_WORLD,ierror)
          !        call MPI_Send(ill_ICU_cases(1,1),block_size,MPI_INTEGER,0,4,MPI_COMM_WORLD,ierror)
          !        call MPI_Send(immune_cases(1,1),block_size,MPI_INTEGER,0,5,MPI_COMM_WORLD,ierror)
          !        call MPI_Send(dead_cases(1,1),block_size,MPI_INTEGER,0,6,MPI_COMM_WORLD,ierror)
          !    else
          !        healthy_cases_final(:,:,index_final(1)) = healthy_cases
          !        inf_noncon_cases_final(:,:,index_final(2)) = inf_noncon_cases
          !        inf_contag_cases_final(:,:,index_final(3)) = inf_contag_cases
          !        ill_contag_cases_final(:,:,index_final(4)) = ill_contag_cases
          !        ill_ICU_cases_final(:,:,index_final(5)) = ill_ICU_cases
          !        immune_cases_final(:,:,index_final(6)) = immune_cases
          !        dead_cases_final(:,:,index_final(7)) = dead_cases
          !        index_final = index_final + 1
          !        call MPI_Recv(healthy_cases_final(1,1,index_final(1)),block_size,MPI_INTEGER,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,ierror)
          !        index_final(1) = index_final(1) + 1
          !!        healthy_cases_final(:,:,index_final(1)) = healthy_cases
          !      
          !        call MPI_Recv(inf_noncon_cases_final(1,1,index_final(2)),block_size,MPI_INTEGER,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,ierror)
          !        index_final(2) = index_final(2) + 1
          !!        inf_noncon_cases_final(:,:,index_final(2)) = inf_noncon_cases
          !        call MPI_Recv(inf_contag_cases_final(1,1,index_final(3)),block_size,MPI_INTEGER,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,ierror)
          !        index_final(3) = index_final(3) + 1
          !        call MPI_Recv(ill_contag_cases_final(1,1,index_final(4)),block_size,MPI_INTEGER,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,ierror)
          !        index_final(4) = index_final(4) + 1
          !        call MPI_Recv(ill_ICU_cases_final(1,1,index_final(5)),block_size,MPI_INTEGER,MPI_ANY_SOURCE,4,MPI_COMM_WORLD,ierror)
          !        index_final(5) = index_final(5) + 1
          !        call MPI_Recv(immune_cases_final(1,1,index_final(6)),block_size,MPI_INTEGER,MPI_ANY_SOURCE,5,MPI_COMM_WORLD,ierror)
          !        index_final(6) = index_final(6) + 1
          !        call MPI_Recv(dead_cases_final(1,1,index_final(7)),block_size,MPI_INTEGER,MPI_ANY_SOURCE,6,MPI_COMM_WORLD,ierror)
          !        index_final(7) = index_final(7) + 1
          !    endif
          !   if (my_rank /= 0)then
          !        call MPI_Send(healthy_cases(1,1),block_size,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierror)
          !   else
          !        healthy_cases_final(:,:,index_final(1)) = healthy_cases
          !        index_final(1) = index_final(1) + 1
          !        if(index_final(1)<=iter)then
          !         call MPI_Recv(healthy_cases_final(1,1,index_final(1)),block_size,MPI_INTEGER,&
          !                      MPI_ANY_SOURCE,0,MPI_COMM_WORLD,req,ierror)
          !         index_final(1) = index_final(1) + 1
          !        end if
          !   end if
          !   
          !   if (my_rank /= 0)then
          !        call MPI_Send(inf_noncon_cases(1,1),block_size,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierror)
          !   else
          !        inf_noncon_cases_final(:,:,index_final(2)) = inf_noncon_cases
          !        index_final(2) = index_final(2) + 1
          !        if (index_final(2)<=iter)then
          !         call MPI_Recv(inf_noncon_cases_final(1,1,index_final(2)),block_size,MPI_INTEGER,&
          !                      MPI_ANY_SOURCE,1,MPI_COMM_WORLD,req,ierror)
          !        index_final(2) = index_final(2) + 1
          !        endif
          !         print *,"req is",req
          !   end if
          !   
          !   if (my_rank /= 0)then
          !        call MPI_Send(inf_contag_cases(1,1),block_size,MPI_INTEGER,0,2,MPI_COMM_WORLD,ierror)
          !   else
          !        inf_contag_cases_final(:,:,index_final(3)) = inf_contag_cases
          !        index_final(3) = index_final(3) + 1
          !        if (index_final(3)<=iter)then
          !         call MPI_Recv(inf_contag_cases_final(1,1,index_final(3)),block_size,MPI_INTEGER,&
          !                      MPI_ANY_SOURCE,2,MPI_COMM_WORLD,req,ierror)
          ! 
          !         index_final(3) = index_final(3) + 1
          !         end if
          !   end if
          ! 
          !   if (my_rank /= 0)then
          !        call MPI_Send(ill_contag_cases(1,1),block_size,MPI_INTEGER,0,3,MPI_COMM_WORLD,ierror)
          !   else
          !        ill_contag_cases_final(:,:,index_final(4)) = ill_contag_cases
          !        index_final(4) = index_final(4) + 1
          !        if (index_final(4)<=iter)then
          !         call MPI_Recv(ill_contag_cases_final(1,1,index_final(4)),block_size,MPI_INTEGER,&
          !                      MPI_ANY_SOURCE,3,MPI_COMM_WORLD,req,ierror)
          !         index_final(4) = index_final(4) + 1
          !         end if
          !   end if
          !   
          !   if (my_rank /= 0)then
          !        call MPI_Send(ill_ICU_cases(1,1),block_size,MPI_INTEGER,0,4,MPI_COMM_WORLD,ierror)
          !   else
          !        ill_ICU_cases_final(:,:,index_final(5)) = ill_ICU_cases
          !        index_final(5) = index_final(5) + 1
          !        if (index_final(5)<iter)then
          !         call MPI_Recv(ill_ICU_cases_final(1,1,index_final(5)),block_size,MPI_INTEGER,&
          !                      MPI_ANY_SOURCE,4,MPI_COMM_WORLD,req,ierror)  
          !         index_final(5) = index_final(5) + 1
          !         end if
          !   end if  
          !   
          !   if (my_rank /= 0)then
          !        call MPI_Send(immune_cases(1,1),block_size,MPI_INTEGER,0,5,MPI_COMM_WORLD,ierror)
          !   else
          !        immune_cases_final(:,:,index_final(6)) = immune_cases
          !        index_final(6) = index_final(6) + 1
          !        if (index_final(6)<=iter)then
          !          call MPI_Recv(immune_cases_final(1,1,index_final(6)),block_size,MPI_INTEGER,&
          !                      MPI_ANY_SOURCE,5,MPI_COMM_WORLD,req,ierror)
          !        index_final(6) = index_final(6) + 1
          !         end if
          !   end if  
          !   
          !   if (my_rank /= 0)then
          !        call MPI_Send(dead_cases(1,1),block_size,MPI_INTEGER,0,6,MPI_COMM_WORLD,ierror)
          !   else
          !        dead_cases_final(:,:,index_final(7)) = dead_cases
          !        index_final(7) = index_final(7) + 1
          !        if (index_final(7)<=iter)then
          !         call MPI_Recv(dead_cases_final(1,1,index_final(7)),block_size,MPI_INTEGER,&
          !                      MPI_ANY_SOURCE,6,MPI_COMM_WORLD,req,ierror)
          !         index_final(7) = index_final(7) + 1
          !         end if
          !   end if   
          Call MPI_ISend(healthy_cases(1,1),block_size,MPI_INTEGER,0,0,MPI_COMM_WORLD,req,ierror)
          Call MPI_ISend(inf_noncon_cases(1,1),block_size,MPI_INTEGER,0,1,MPI_COMM_WORLD,req,ierror)
          Call MPI_ISend(inf_contag_cases(1,1),block_size,MPI_INTEGER,0,2,MPI_COMM_WORLD,req,ierror)
          Call MPI_ISend(ill_contag_cases(1,1),block_size,MPI_INTEGER,0,3,MPI_COMM_WORLD,req,ierror)
          Call MPI_ISend(ill_ICU_cases(1,1),block_size,MPI_INTEGER,0,4,MPI_COMM_WORLD,req,ierror)
          Call MPI_ISend(immune_cases(1,1),block_size,MPI_INTEGER,0,5,MPI_COMM_WORLD,req,ierror)
          Call MPI_ISend(dead_cases(1,1),block_size,MPI_INTEGER,0,6,MPI_COMM_WORLD,req,ierror)
       End Do     ! end do it_ss
    Else ! the very begining if:my_rank /= 0
       Do k = 1,iter
          Call MPI_Recv(healthy_cases_final(1,1,index_final(1)),block_size,MPI_INTEGER,&
               MPI_ANY_SOURCE,0,MPI_COMM_WORLD,req,ierror)
          index_final(1) = index_final(1) + 1
          Call MPI_Recv(inf_noncon_cases_final(1,1,index_final(2)),block_size,MPI_INTEGER,&
               MPI_ANY_SOURCE,1,MPI_COMM_WORLD,req,ierror)  
          index_final(2) = index_final(2) + 1
          Call MPI_Recv(inf_contag_cases_final(1,1,index_final(3)),block_size,MPI_INTEGER,&
               MPI_ANY_SOURCE,2,MPI_COMM_WORLD,req,ierror)
          index_final(3) = index_final(3) + 1
          Call MPI_Recv(ill_contag_cases_final(1,1,index_final(4)),block_size,MPI_INTEGER,&
               MPI_ANY_SOURCE,3,MPI_COMM_WORLD,req,ierror)
          index_final(4) = index_final(4) + 1
          Call MPI_Recv(ill_ICU_cases_final(1,1,index_final(5)),block_size,MPI_INTEGER,&
               MPI_ANY_SOURCE,4,MPI_COMM_WORLD,req,ierror)
          index_final(5) = index_final(5) + 1
          Call MPI_Recv(immune_cases_final(1,1,index_final(6)),block_size,MPI_INTEGER,&
               MPI_ANY_SOURCE,5,MPI_COMM_WORLD,req,ierror)
          index_final(6) = index_final(6) + 1
          Call MPI_Recv(dead_cases_final(1,1,index_final(7)),block_size,MPI_INTEGER,&
               MPI_ANY_SOURCE,6,MPI_COMM_WORLD,req,ierror)
          index_final(7) = index_final(7) + 1
       End Do
    End If

!!$OMP END PARALLEL DO

    If (my_rank == 0)Then
       iter_pass_handle = (/lhc(1,iter),lhc(2,iter),lhc(3,iter),&
            lhc(6,iter),lhc(7,iter),lhc(8,iter)/)
       !     do k = 1,iter
       !         print *,healthy_cases_final(:,:,k)
       !     end do
       Call write_data_v2(healthy_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,1)
       Call write_data_v2(inf_noncon_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,2)
       Call write_data_v2(inf_contag_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,3)
       Call write_data_v2(ill_contag_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,4)
       Call write_data_v2(ill_ICU_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,5)
       Call write_data_v2(immune_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,6)
       Call write_data_v2(dead_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,7)
    End If

    1000 continue
    Call MPI_Finalize(ierror)

  End Subroutine COVID19_Spatial_Microsimulation_for_Germany


  Subroutine write_data_v2(healthy_cases_final,iter_pass_handle,R0Change,counties_index,type_file)
    Real,Dimension(:)        :: iter_pass_handle(:)
    Real,Dimension(:,:)      :: R0Change
    Real,Allocatable         :: R0change_rep(:,:),R0change_exp(:,:),temp_output(:)
    Integer,Dimension(:,:,:) :: healthy_cases_final
    Integer,Dimension(:)     :: counties_index
    Integer                  :: iter
    Integer                  :: type_file

    Integer :: county_size, count 
    Character*15                :: iter_char(6)
    Character*5,Allocatable     :: R0change_name(:)
    Character*2                 :: counties(16)
    Integer,Allocatable         :: counties_index_out(:),iter_array(:)

    Character*10,Allocatable    :: Label(:)

    Integer date_time(8),i,j,k
    Character*10 b(3)
    Character*4 year
    Character*2 day,month
    Character*8 time
    Character*3 temp_char

    Character*15 dir_prefix

    iter = Size(healthy_cases_final,DIM=3)
    county_size = Size(healthy_cases_final,DIM=1)
    ! print *,"iter is",iter
    ! print "county_size is",county_size
    iter_char = (/"sam_size       ",&
         "R0             ",&
         "icu_dur        ",&
         "w_int          ",&
         "w.obs          ",&
         "w.obs.by.sate  "/)

    counties  = (/"SH","HH","NI","HB","NW","HE","RP","BW","BY","SL",&
         "BE","BB","MV","SN","ST","TH"/)

    Allocate(R0change_rep(Size(R0Change),1))
    Allocate(R0change_exp(iter*county_size,Size(R0change_rep)))
    Allocate(temp_output(iter*county_size))
    Allocate(R0change_name(Size(R0Change)))
    Allocate(iter_array(county_size*iter))
    Allocate(counties_index_out(county_size*iter))
    Allocate(Label(Size(healthy_cases_final,dim=2)))

    Call date_and_Time(b(1),b(2),b(3),date_time)

    count = 1
    Do i = 1,Size(counties)
       Do j = 1,Size(R0Change,2)
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
    If (type_file == 1)Then
       Open (101, file = Trim(dir_prefix)//time//"healthy_cases.csv")
    End If

    If (type_file == 2)Then
       Open (101, file = Trim(dir_prefix)//time//"inf_noncon_cases.csv")
    End If

    If (type_file == 3)Then
       Open (101, file = Trim(dir_prefix)//time//"inf_contag_cases.csv")
    End If

    If (type_file == 4)Then
       Open (101, file = Trim(dir_prefix)//time//"ill_contag_cases.csv")
    End If

    If (type_file == 5)Then
       Open (101, file = Trim(dir_prefix)//time//"ill_ICU_cases.csv")
    End If

    If (type_file == 6)Then
       Open (101, file = Trim(dir_prefix)//time//"immune_cases.csv")
    End If

    If (type_file == 7)Then
       Open (101, file = Trim(dir_prefix)//time//"dead_cases.csv")
    End If

10  Format(1x,*(g0,","))

    !write the head file
    Write(101,10)Label,iter_char,R0change_name,"iter","x.dist_id"

    R0change_rep = Reshape(R0Change,(/Size(R0change),1/))
    R0change_exp = Transpose(Spread(R0change_rep(:,1),2,iter*county_size))

    count = 1
    Do i = 1,iter
       Do j = 1,county_size
          Write(101,10)healthy_cases_final(j,:,i),iter_pass_handle,R0change_exp(count,:),&
               i,counties_index(j)
          count = count+1
       End Do
    End Do

    Close(101)
  End Subroutine write_data_v2

End Module kernel
