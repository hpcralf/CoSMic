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

  use timer

  use param_tree

  use precision
  Use global_constants
  use global_types
  use global_vars
  
  Use list_variable
  Use support_fun
  
  Implicit None

Contains
  
  Subroutine COVID19_Spatial_Microsimulation_for_Germany( &
       iol, pspace, counties_index &
       )

    !===========================================================================
    ! Declaration
    !===========================================================================

    !Include 'mpif.h'

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
    Integer(kind=ik), Allocatable, Dimension(:,:) :: R0change
    Logical             :: R0delay
    Integer             :: R0delay_days
    Character(len=:),Allocatable :: R0delay_type

!!!!!-----7.  Define whether transition probabilities should differ by age and sex
    Character*10        :: control_age_sex
    Character*10        :: seed_before_char,seed_temp
    Character*10,Allocatable :: seed_seq(:),seed_inf_cont_seq(:),seed_inf_ncont_seq(:)
    Character*10,Allocatable :: seed_d_seq(:)
    Integer             :: days

    Integer             :: n_direct,n_directv,n_directl,n_dist
    Integer             :: size_lhc
    Real(kind=rk),Allocatable    :: lhc(:,:)
    Integer(kind=ik)             :: tmp_i8

!!!!!-----8. variables for do loop -------
    Real,Allocatable    :: icu_risk(:),surv_ill(:),surv_icu(:)
    Integer,Allocatable :: temp1(:),temp(:),targert_icu_risk_index(:)

    Integer,Allocatable,Dimension(:) :: temp_s
    
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
       Integer,Allocatable             :: dist_id_rn(:)
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

    Integer,Allocatable             :: start_value_tot(:)

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

    Integer                         :: timestep

    Integer,Allocatable             :: tmp_d_new(:),tmp_count(:)
    Integer,Allocatable             :: susceptible(:),contagious_dist_id(:),contagious_index(:),denominator(:),&
         revers_proj(:),final_count(:),dist_id_temp(:),ill_index(:),&
         ill_dist_id(:)
    Real,Allocatable                :: contagious(:)
    Integer                         :: at_risk
    Integer                         :: initial_sick
    Real                            :: n_contagious,between_weight,within_weight
    Real,Allocatable                :: exp_infect(:)
    Integer,Allocatable             :: check_days(:)

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
    Character*10                    :: temp_mod
    character(len=:),allocatable    :: seed_date
    Integer                         :: ierror,size_of_process,my_rank
    Integer                         :: index_final(7),block_size
    Integer,Allocatable             :: req(:)

    Real                            :: iter_pass_handle(6)

    Type(tTimer)                    :: timer
    Integer, Dimension(8)           :: rt

    integer                         :: tar

    Real(kind=pt_rk),Dimension(:,:),Allocatable :: R0_effects
    Integer(kind=ik),Dimension(:)  ,Allocatable :: dist_id_cref 
    
    Integer(kind=ik), Dimension(0:16)         :: istate_count
    Integer(kind=ik)                         :: pop_size
    Integer(kind=ik)                         :: num_counties
    Integer(kind=ik)                         :: ii
    
    ! should import some reliable romdon seed generation code here
    !seed_base = ??

    !===========================================================================
    ! Implementation
    !===========================================================================
    call pt_get("#iter",iter)

    index_final = 1

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
    call pt_get("#R0change"     ,R0change    )
    call pt_get("#R0delay"      ,R0delay     )
    call pt_get("#R0delay_days" ,R0delay_days)
    call pt_get("#R0delay_type" ,R0delay_type)

    time_n = Maxval(R0change) + 1
    
!!!!!-----7.  Define whether transition probabilities should differ by age and sex
    control_age_sex     = "age"
    days             = 1

    call pt_get("#seed_date",seed_date)
    seed_date        = add_date(seed_date,days)

    days             = -1-seed_before
    seed_before_char = add_date(seed_date,days)
    seed_seq         = generate_seq(seed_before_char,seed_date)
    
    write(un_lf,PTF_sep)
    write(un_lf,PTF_M_A)"Seed sequence for ill cases:",seed_seq
    
    !Derive dates of infections for those that are inf_cont,
    !but are not yet aware about it (will be registered the
    !next two days)
    seed_temp        = add_date(seed_date,cont_dur)
    seed_inf_cont_seq  =  generate_seq(add_date(seed_date,1),seed_temp)

    write(un_lf,PTF_sep)
    write(un_lf,PTF_M_A)"Seed sequence for infected contagious cases:",seed_inf_cont_seq

    !Derive dates of infections for those that are inf_cont,
    !but are not yet aware about it (will be registered the
    !next 3-5 days)
    seed_inf_ncont_seq = generate_seq(add_date(seed_date,cont_dur+1),add_date(seed_date,inf_dur+cont_dur))

    write(un_lf,PTF_sep)
    write(un_lf,PTF_M_A)"Seed sequence for infected non-contagious cases:",seed_inf_ncont_seq

    !! Setup of latin hypercube. ===============================================
    !! This part should done by the code, but here is set manually for 
    !! simplicity. This part should be seperated away into the preprocessing step
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
       size_lhc = size_lhc + size(iol%R0_effect%data)
    End If

    Allocate(lhc(size_lhc,iter))

    ! lhc(1,:)            = pspace%sam_size%param
    Do i = 1,n_direct
       lhc(i,:) = pspace%Ps_scalar_list(i)%param
    End Do

    call pt_get("#sam_size",tmp_i8)
    lhc(1,:) = tmp_i8

    Do i = 1,iter
       lhc(n_direct+1:Size(lhc,dim=1),i) = Reshape(transpose(iol%R0_effect%data),&
            Shape(lhc(n_direct+1:Size(lhc,dim=1),1)))
    End Do
    ! print *, "after reshape is",reshape(pspace%ROeffect_ps%param,shape(lhc(n_direct+1:size(lhc,dim=1),1)))

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

    write(un_lf,PTF_SEP)
    write(un_lf,PTF_M_A)"ICU risk per age group, sex, and duration in ICU."
    Do ii=1, Size(icu_risk)*ill_dur
       write(un_lf,'(A4,A2,F6.3,I3)') &
            icu_risk_list%age(ii)    , icu_risk_list%sex(ii), &
            icu_risk_list%risk(ii), icu_risk_list%dur(ii)
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

    write(un_lf,PTF_SEP)
    write(un_lf,PTF_M_A)"Chance of survival per age group and sex."
    Do ii=1, 2*19
       write(un_lf,'(A4,A2,F6.3)') &
            surv_ill_list%age(ii)    , surv_ill_list%sex(ii), &
            surv_ill_list%risk(ii)
    End Do

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

    write(un_lf,PTF_SEP)
    write(un_lf,PTF_M_A)"Chance of survival in ICU per age group and sex."
    Do ii=1, 2*19
       write(un_lf,'(A4,A2,F6.3)') &
            surv_icu_list%age(ii)    , surv_icu_list%sex(ii), &
            surv_icu_list%risk(ii)
    End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! the following part was in the loop at R code
!!! put it outside to avoid memoery problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    num_counties = Size(counties_index)
    
    Allocate(healthy_cases_final(num_counties,time_n,iter))
    Allocate(inf_noncon_cases_final(num_counties,time_n,iter))
    Allocate(inf_contag_cases_final(num_counties,time_n,iter))
    Allocate(ill_contag_cases_final(num_counties,time_n,iter))
    Allocate(ill_ICU_cases_final(num_counties,time_n,iter))
    Allocate(immune_cases_final(num_counties,time_n,iter))
    Allocate(dead_cases_final(num_counties,time_n,iter))

    Allocate(inf_cases(num_counties,time_n))
    Allocate(icu_cases(num_counties,time_n))

    Allocate(healthy_cases(num_counties,time_n))
    Allocate(inf_noncon_cases(num_counties,time_n))
    Allocate(inf_contag_cases(num_counties,time_n))
    Allocate(ill_contag_cases(num_counties,time_n))
    Allocate(ill_ICU_cases(num_counties,time_n))
    Allocate(immune_cases(num_counties,time_n))
    Allocate(dead_cases(num_counties,time_n))
    Allocate(dead_cases_bICU(num_counties,time_n))
    Allocate(mod_inf_cases(num_counties,time_n))
    Allocate(org_noncon_cases(num_counties,time_n))

    block_size = num_counties * time_n

    max_date = find_max_date(iol%seed_date)

    !** Allocate and setup dist_id cross reference -----------------------------
    allocate(dist_id_cref(minval(counties_index):maxval(counties_index)))
    dist_id_cref = -1
    
    Do ii = 1, num_counties
       dist_id_cref(counties_index(ii)) = ii
    End Do
    
!!!=============================================================================
!!! Iteration over parameter space
!!!=============================================================================
    Do it_ss = 1,  Size(lhc,dim=2)

       call start_timer("Init Sim Loop",reset=.FALSE.)
       
       Call random_Seed()

       !!=======================================================================
       !! Init population                 
       iol%pop_total = Nint(Real(iol%pop_total)/Real(Sum(iol%pop_total))* Real(lhc(1,it_ss)))
       pop_size      = Sum(iol%pop_total)
       
       If (.Not.Allocated(temp_s)) allocate(temp_s(pop_size))
       
       If (.Not.Allocated(sim%dist_id))Then
          Allocate(sim%dist_id(pop_size))
          Allocate(sim%sex(pop_size))
          Allocate(sim%age(pop_size))
          Allocate(sim%t1(pop_size))
          Allocate(sim%t2(pop_size))
          Allocate(sim%d(pop_size))
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

       write(un_lf,PTF_SEP)
       Write(un_lf,PTF_M_AI0)"Size of population is", sum(iol%pop_total)

       !!=======================================================================
       !! seed infections
       temp = get_index(iol%death_distid,counties_index)
       
       iol%death_distid = iol%death_distid(temp)
       iol%death_date   = iol%death_date(temp)
       iol%death_cases  = iol%death_cases(temp)
       !skip line 1113, 1115

       temp = get_index(iol%seed_date,seed_seq)
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

       days             = -1

       seed_date        = add_date(seed_date,days)

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

       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"First 20 and last 20 seeds for ill cases."
       Do ii=1, 20
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_ill%dist_id(ii)    , seed_ill%date(ii), &
               seed_ill%cases(ii),temp(ii)
       End Do
       write(un_lf,"('...')")
       Do ii = Size(temp)-19, Size(temp)
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_ill%dist_id(ii)    , seed_ill%date(ii), &
               seed_ill%cases(ii),temp(ii)
       End Do
       
       If (.Not.Allocated(seed_inf_cont_dur))Then
          Allocate(seed_inf_cont_dur(Size(seed_inf_cont%date)))
       End If
       temp_int = Date2Unixtime(seed_date)
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

       !         seed_inf_cont = set_value(seed_inf_cont,temp(1:temp_int))
       seed_inf_cont%dist_id  = seed_inf_cont%dist_id(temp(1:temp_int))
       seed_inf_cont%date     = seed_inf_cont%date(temp(1:temp_int))
       seed_inf_cont%cases    = seed_inf_cont%cases(temp(1:temp_int))

       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"First 20 and last 20 seeds for contagious cases."
       Do ii=1, 20
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_inf_cont%dist_id(ii)    , seed_inf_cont%date(ii), &
               seed_inf_cont%cases(ii),seed_inf_cont_dur(ii)
       End Do
       write(un_lf,"('...')")
       Do ii = Size(temp)-19, Size(temp)
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_inf_cont%dist_id(ii)    , seed_inf_cont%date(ii), &
               seed_inf_cont%cases(ii),seed_inf_cont_dur(ii)
       End Do
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
       !         seed_inf_ncont = set_value(seed_inf_ncont,temp(1:temp_int))
       seed_inf_ncont%dist_id  = seed_inf_ncont%dist_id(temp(1:temp_int))
       seed_inf_ncont%date     = seed_inf_ncont%date(temp(1:temp_int))
       seed_inf_ncont%cases    = seed_inf_ncont%cases(temp(1:temp_int))

       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"First 20 and last 20 seeds for non contagious cases."
       Do ii=1, 20
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_inf_ncont%dist_id(ii)    , seed_inf_ncont%date(ii), &
               seed_inf_ncont%cases(ii),seed_inf_ncont_dur(ii)
       End Do
       write(un_lf,"('...')")
       Do ii = Size(temp)-19, Size(temp)
          write(un_lf,'(I6,A12,I6,I6)') &
               seed_inf_ncont%dist_id(ii)    , seed_inf_ncont%date(ii), &
               seed_inf_ncont%cases(ii),seed_inf_ncont_dur(ii)
       End Do

       !       skip line 1163,1166,1169,1172 since the mechanism of 
       !       aggregate is not clear, and it seems not changing anything at all
       !       do scaling
       ! skip scaling ,since they are all 1
       !         seed_ill%cases = seed_ill%cases * sam_prop_ps(seed_ill%dist_id/1000)
       !
       !         seed_inf_cont%cases = seed_inf_cont%cases * sam_prop_ps(seed_ill%dist_id/1000)
       !
       !         seed_inf_ncont%cases = seed_inf_ncont%cases * sam_prop_ps(seed_ill%dist_id/1000)

       write(un_lf,PTF_SEP)
       write(un_lf,PTF_M_A)"Seeds per county."
       write(un_lf,'(5(A10,1X))')"county","inf_ncont","inf_cont","inf_ill","inf_dth"
       
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

          temp   = get_index(iol%death_distid,county)
          inf_dth = Sum(iol%death_cases(temp))

          write(un_lf,'(11(I11))')county,inf_ncont,inf_cont,inf_ill,inf_dth,&
               minval(inf_nc_d),maxval(inf_nc_d),minval(inf_c_d),maxval(inf_c_d),&
               minval(il_d),maxval(il_d)
          
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

       !! ----------------------------------------------------------------------
       !! Convert from Weekly to daily R0_effects ------------------------------
       R0_daily = R0_force *lhc(2,it_ss)/Real(Real(cont_dur)+Real(ill_dur)*less_contagious) + &
            (1-R0_force)*lhc(2,it_ss)/Real(cont_dur+ill_dur)

       ! this block simplifies the if judgment
       If (.Not.Allocated(R0matrix))Then
          Allocate(R0matrix(num_counties,time_n-1))
       End If
       R0matrix = R0_daily
       n_change  = Size(R0change,dim=2)

       Do i = 1,n_change

          !use counties
          ! give up using character to find the number, but use the index directly
          temp  = 8 + (i-1) * Size(iol%R0_effect%data,dim= 2) + (counties_index/1000)

          getchange = lhc(temp,it_ss)

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

       !! ----------------------------------------------------------------------
       !! Init result fields
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
       sim%t2                 = missing !sim%t1

       !** Set up dist_id renumbered cross_reference ---------------------------
       sim%dist_id_rn = dist_id_cref(sim%dist_id)

       If (.Not.Allocated(tmp_d_new))Then
          Allocate(tmp_d_new(Size(sim%d)))
       End If

       call end_timer("Init Sim Loop")
    
!!!=============================================================================
!!! Simulation Loop ============================================================
!!!=============================================================================
       call start_timer("Sim Loop",reset=.FALSE.)

       Do timestep = 2,time_n

          call start_timer("+- From healthy to infected",reset=.FALSE.)
          !** Loop local Copy of population ***
          tmp = sim

          !** Init new durations ***
          tmp_d_new = missing

          !** Get number of all who are in an ill or suceptible state and ***
          !** check if anyone is left                                     ***
          istate_count = 0
          Do ii = 1, pop_size
             istate_count(tmp%t1(ii)) = istate_count(tmp%t1(ii)) + 1
             !istate_count(tmp%d(ii)) = istate_count(tmp%d(ii)) + 1
          End Do

          !** DEBUG --- Population Summary -------------------------------------
          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary @ start of step:",timestep-1
          write(un_lf,'(11(I10))')0,1,2,3,4,5,6,7,8,9,10
          write(un_lf,'(11(I10))')istate_count(0:10)
          !** DEBUG --- Population Summary -------------------------------------
          
          If (sum(istate_count((/&
               inf_contag,inf_noncon, &
               ill_ICU,   ill_contag, &
               healthy                &
               /))) == 0)Then
             write(*,*)"All are dead or immune."
             goto 1000
          End If

          !** Get Indicees of all who are suceptible ***
          tmp_count = get_index(tmp%t1,healthy)

          !** Number of people who are at risk ***
          at_risk = Size(tmp_count)

          !** Number of people in states who can infect others ***
          contagious_index   = get_index(tmp%t1,inf_contag)
          ill_index          = get_index(tmp%t1,ill_contag)

          !** Where do people who can infect others live ? ***
          contagious_dist_id = get_unique(tmp%dist_id(contagious_index))
          ill_dist_id        = get_unique(tmp%dist_id(ill_index))

          If (.Not.Allocated(contagious))Then
             Allocate(contagious(num_counties))
          End If

          contagious = Real(sum_byindex(tmp%dist_id(contagious_index),counties_index))
          contagious = contagious + Real(sum_byindex(tmp%dist_id(ill_index),counties_index)) * less_contagious

          n_contagious = Sum(contagious)
          write(un_lf,PTF_M_AI0)"n_contagious:",nint(n_contagious),"at_risk to be infected:",at_risk

          If (at_risk > 0 .And. n_contagious > 0) Then

             exp_infect = contagious * R0matrix(:,timestep-1)
            
             temp = get_index(iol%connect_work_distid,counties_index)

             connect = iol%connect_work(temp,temp)
             Do i = 1,Size(contagious_dist_id)
                connect(:,i) = connect(:,i)/Sum(connect(:,i))
             End Do

             between_weight = 1 - lhc(6,it_ss) ! 1 - lhc(it_ss,"w_int")

             within_weight = 1 - between_weight
             exp_infect = within_weight * exp_infect + between_weight * Matmul (exp_infect, connect)

             If (.Not.Allocated(denominator))Then
                Allocate(denominator(num_counties))
             End If

             
             If (.Not.immune_stop) Then
                denominator = sum_byindex(tmp%dist_id(tmp_count),counties_index)
             Else
                temp = get_index(tmp%t1,(/healthy,immune,inf_noncon,inf_contag,ill_contag/))

                denominator = sum_byindex(tmp%dist_id(temp),counties_index)
             End If

             risk = exp_infect/Real(denominator)
             Where (risk > 1000)
                risk = 0
             End Where
             Where (risk > 1)
                risk = 1
             End Where

             risk = lhc(6,it_ss)*risk + (1-lhc(6,it_ss))* Matmul(risk,connect)
             
             temp = sim%dist_id_rn(tmp_count)
                                     
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

             temp_int = 0
             sick = 0
             Do i = 1,Size(runif)
                If (runif(i) <= prob(i)) Then
                   temp_int = temp_int + 1
                   sick(i)  = 1
                End If
             End Do
             
             write(un_lf,PTF_M_AI0)"# of newly infected:",sum(sick)

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
                Allocate(final_count(num_counties))

             Else
                Deallocate(final_count)
                Allocate(final_count(num_counties))
             End If
             final_count = 0

             final_count = sum_byindex(tmp%dist_id(revers_proj),counties_index)
             inf_cases(:,timestep) = final_count
             org_noncon_cases(:,timestep) = final_count

             target_date = add_date(seed_date,timestep-1+7)

             
!!$             If (lhc(7,it_ss) > 0 .And. Date2Unixtime(target_date) < max_date .And. initial_sick > 0)Then
!!$
!!$                temp = get_index(iol%seed_date,target_date)
!!$                target_inf%cases   = iol%seed_cases(temp)
!!$                target_inf%dist_id = iol%seed_distid(temp)
!!$                target_inf%date    = iol%seed_date(temp)
!!$                If (.Not.Allocated(final_count))Then
!!$                   Allocate(final_count(num_counties))
!!$                Else
!!$                   Deallocate(final_count)
!!$                   Allocate(final_count(num_counties))
!!$                End If
!!$                final_count = 0
!!$                final_count = sum_byindex(target_inf%dist_id,counties_index)
!!$                If (Sum(final_count) > 0) Then
!!$                   If(lhc(8,it_ss) == 1) Then
!!$                      state_id = get_unique(counties_index/1000)
!!$                      If (.Not.Allocated(mod_inf))Then
!!$                         Allocate(mod_inf(num_counties))
!!$                      End If
!!$                      mod_inf = 0
!!$                      Do i = 1, Size(state_id)
!!$                         temp = get_index(counties_index/1000,state_id(i))
!!$                         prop_inf_cases = inf_cases(temp,timestep)
!!$                         temp_int       = Sum(prop_inf_cases)
!!$                         If (temp_int > 0) Then
!!$                            prop_inf_cases = prop_inf_cases / temp_int
!!$                         Else
!!$                            prop_inf_cases = 0
!!$                         End If
!!$
!!$                         prop_target_inf = final_count(temp)
!!$                         temp_int = Sum(prop_target_inf)
!!$                         If (temp_int > 0) Then
!!$                            prop_target_inf = prop_target_inf / temp_int
!!$                         Else
!!$                            prop_target_inf = 0
!!$                         End If
!!$                         prop_target = lhc(7,it_ss) * Real(prop_target_inf) + (1 - lhc(7,it_ss))*Real(prop_inf_cases)
!!$
!!$                         mod_inf(temp) = Nint(prop_target* Real(temp_int)) - inf_cases(temp,timestep)
!!$                      End Do
!!$
!!$                   Else
!!$                      prop_inf_cases = inf_cases(:,timestep) / initial_sick
!!$                      prop_target = lhc(7,it_ss) * Real(final_count)/Real(Sum(final_count)) + &
!!$                           (1 - lhc(7,it_ss))*Real(prop_inf_cases)
!!$                      ! line 1675
!!$                      mod_inf = Nint(prop_target * Real(initial_sick))-inf_cases(:,timestep)
!!$
!!$                   End If
!!$
!!$                   Call shift_cases(sick,mod_inf,tmp%dist_id(tmp_count),counties_index)
!!$!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!can we write the following code in a more elgant
!!$!!!way?
!!$!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$                   temp = get_index(sick,1)
!!$                   ! do a revers projection
!!$                   temp = tmp_count(temp)
!!$                   If (.Not.Allocated(final_count))Then
!!$                      Allocate(final_count(num_counties))
!!$                   End If
!!$                   If (Allocated(dist_id_temp))Then
!!$                      Deallocate(dist_id_temp)
!!$                   End If
!!$                   If (Allocated(case_count))Then
!!$                      Deallocate(case_count)
!!$                   End If
!!$
!!$                   !                     call sum_bygroup_distID(tmp%dist_id(temp),case_count,dist_id_temp)
!!$                   final_count = 0
!!$                   !                     do i =1,size(counties_index)
!!$                   !                         temp1 = get_index(tmp%dist_id(temp),counties_index(i))
!!$                   !                         final_count(i) = size(temp1)
!!$                   !                     enddo
!!$                   final_count = sum_byindex(tmp%dist_id(temp),counties_index)
!!$                   inf_cases(:,timestep) = final_count
!!$
!!$                   mod_inf_cases(:,timestep) = mod_inf
!!$                Else
!!$                   Print *,"no cases shifted"
!!$                End If
!!$
!!$             Else
                mod_inf_cases(:,timestep) = 0
!!$             Endif

             !skip 1742 and 1743 since they are not actually do
             !anything?

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
          
          !! #############################################################
          !! State: Infected, non-contagious #############################
          
          !If day limit not reach: Stays the same
          temp = condition_and(tmp%t1,inf_noncon,"e",tmp%d,inf_dur,"l")
          tmp%t2(temp) = inf_noncon

          write(un_lf,PTF_M_AI0)"# staying in inf_noncon:",size(temp)
          
          !Update day(current count plus 1)
          temp = condition_and(tmp%t1,inf_noncon,"e",tmp%t2,inf_noncon,"e")
          !         tmp%d(temp)= tmp%d(temp) + 1
          tmp_d_new(temp) = tmp%d(temp) + 1

          !** DEBUG --- Population Summary -------------------------------------
          istate_count = 0
          Do ii = 1, pop_size
             istate_count(tmp%t2(ii)) = istate_count(tmp%t2(ii)) + 1
             !istate_count(tmp%d(ii)) = istate_count(tmp%d(ii)) + 1
          End Do

          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary after new inf_noncon",timestep-1
          write(un_lf,'(11(I10))')0,1,2,3,4,5,6,7,8,9,10
          write(un_lf,'(11(I10))')istate_count(0:10)
          !** DEBUG --- Population Summary -------------------------------------

          !         print *,"number of infnoncon update",size(temp)
          ! if day limit reached: move to contagious
          temp = condition_and(tmp%t1,inf_noncon,"e",tmp%d,inf_dur,"g")
          tmp%t2(temp) = inf_contag

          write(un_lf,PTF_M_AI0)"# moving to inf_contag:",size(temp)
          
          !Update days (Reset counter to 1 if moved)
          temp = condition_and(tmp%t1,inf_noncon,"e",tmp%t2,inf_contag,"e")
          !         tmp%d(temp) = 1
          tmp_d_new(temp) = 1
          !         print *,"number of infcontag update",size(temp)
          call end_timer("+- From healthy to infected")

!!!=============================================================================
!!!== state infected,contagious ================================================
          call start_timer("+- From infected to inf_contag",reset=.FALSE.)

          !if day limit not reached: stay the same
          temp = condition_and(tmp%t1,inf_contag,"e",tmp%d,cont_dur,"l")
          tmp%t2(temp) = inf_contag
          !         print *,"number of infcontag stay",size(temp)
          write(un_lf,PTF_M_AI0)"# staying in inf_contag:",size(temp)
          
          !update days
          temp = condition_and(tmp%t1,inf_contag,"e",tmp%t2,inf_contag,"e")
          !         tmp%d(temp) = tmp%d(temp) + 1
          tmp_d_new(temp) = tmp%d(temp) + 1
                   
          !** DEBUG --- Population Summary -------------------------------------
          istate_count = 0
          Do ii = 1, pop_size
             istate_count(tmp%t2(ii)) = istate_count(tmp%t2(ii)) + 1
             !istate_count(tmp%d(ii)) = istate_count(tmp%d(ii)) + 1
          End Do

          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary after new inf_contag",timestep-1
          write(un_lf,'(11(I10))')0,1,2,3,4,5,6,7,8,9,10
          write(un_lf,'(11(I10))')istate_count(0:10)
          !** DEBUG --- Population Summary -------------------------------------

          !         print *,"number of infcontag stay update",size(temp)
          ! if day limit reached: move to contagious
          temp = condition_and(tmp%t1,inf_contag,"e",tmp%d,cont_dur,"g")
          !         print *,sim%d
          tmp%t2(temp) = ill_contag

          write(un_lf,PTF_M_AI0)"# Getting ill:",size(temp)

          !update days(reset to 1 if becoming ill/contagious)
          temp = condition_and(tmp%t1,inf_contag,"e",tmp%t2,ill_contag,"e")
          !         tmp%d(temp) = 1
          tmp_d_new(temp) = 1
          !         print *,"number of illcontag update",size(temp)
          call end_timer("+- From infected to inf_contag")
         
!!!state:ill,contagious
          !dying population at risk
          temp = get_index(tmp%t1,ill_contag)
          at_risk = Size(temp)

          write(un_lf,PTF_M_AI0)"at_risk to die:",at_risk
          
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

          !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !** Get number of all who are in an ill or suceptible state and ***
          !** check if anyone is left                                     ***
!!$          istate_count = 0
!!$          Do ii = 1, pop_size
!!$             istate_count(tmp%t2(ii)) = istate_count(tmp%t2(ii)) + 1
!!$          End Do
!!$          write(un_lf,PTF_SEP)
!!$          write(un_lf,PTF_M_A)"Population Summary after infection evaluation"
!!$          Do ii = 0, 6
!!$             write(un_lf,'(i2,I10)')ii,istate_count(ii)
!!$          End Do
          !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          ! Save dead before moving to ICU
          temp1 = get_index(tmp%t2,dead)

          final_count = sum_byindex(tmp%dist_id(temp1),counties_index)
          dead_cases_bICU(:,timestep) = final_count


          !moving to icu: population at risk
          temp = condition_and(tmp%t1,ill_contag,"e",tmp%t2,dead,"n")

          at_risk = Size(temp)
          
          write(un_lf,PTF_M_AI0)"at_risk to move to ICU:",at_risk
          
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

          write(un_lf,PTF_M_AI0)"# staying in ill_contag:",size(temp)
          
          ! update days
          temp = condition_and(tmp%t1,ill_contag,"e",tmp%t2,ill_contag,"e")
          !         tmp%d(temp) = tmp%d(temp) + 1
          tmp_d_new(temp) = tmp%d(temp) + 1

          !getting healthy/immune if day limit reached
          temp = condition_and(tmp%t1,ill_contag,"e",tmp%d,ill_dur,"g")
          temp1= get_index(tmp%t2,(/healthy,immune,inf_noncon,inf_contag,ill_contag/))
          temp = find_and(temp,temp1)
          tmp%t2(temp) = immune

          write(un_lf,PTF_M_AI0)"# becoming immune:",size(temp)
          
          !update days
          temp = condition_and(tmp%t1,ill_contag,"e",tmp%t2,immune,"e")
          !         tmp%d(temp) = 1
          tmp_d_new(temp) = 1

          !** DEBUG --- Population Summary -------------------------------------
          istate_count = 0
          Do ii = 1, pop_size
             istate_count(tmp%t2(ii)) = istate_count(tmp%t2(ii)) + 1
             !istate_count(tmp%d(ii)) = istate_count(tmp%d(ii)) + 1
          End Do

          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary after new ill_contag",timestep-1
          write(un_lf,'(11(I10))')0,1,2,3,4,5,6,7,8,9,10
          write(un_lf,'(11(I10))')istate_count(0:10)

          !** DEBUG --- Population Summary -------------------------------------
          
!!!state ICU
          temp = get_index(tmp%t1,ill_ICU)
          at_risk = Size(temp)

          write(un_lf,PTF_M_AI0)"at_risk to die in ICU:",at_risk
                    
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
          temp = condition_and(tmp%t1,ill_ICU,"e",tmp%d,Int(lhc(3,it_ss)),"g")
          temp1= get_index_not(tmp%t2,dead)
          temp = find_and(temp,temp1)
          tmp%t2(temp) = immune

          !update days
          temp = condition_and(tmp%t1,ill_ICU,"e",tmp%t2,immune,"e")
          !         tmp%d(temp) = 1
          tmp_d_new(temp) = 1

          
          !** DEBUG --- Population Summary -------------------------------------
          istate_count = 0
          Do ii = 1, pop_size
             istate_count(tmp%t2(ii)) = istate_count(tmp%t2(ii)) + 1
             !istate_count(tmp%d(ii)) = istate_count(tmp%d(ii)) + 1
          End Do

          write(un_lf,PTF_SEP)
          write(un_lf,PTF_M_AI0)"Population Summary after step",timestep-1
          write(un_lf,'(11(I10))')0,1,2,3,4,5,6,7,8,9,10
          write(un_lf,'(11(I10))')istate_count(0:10)
          write(*,'(7(I10))')istate_count(0:6)
          !** DEBUG --- Population Summary -------------------------------------
          
          !check NA set to missing
          !skip this line 1995

          !move to main data frame
          !temp = get_index(sim%t1,(/healthy,ill_ICU,inf_noncon,inf_contag,ill_contag/))
          !sim%t2(temp) = tmp%t2(temp)
          !         sim%d(temp)  = tmp%d(temp)
          !         sim%d(temp)  = tmp_d_new(temp)
          !         sim%t2 = tmp%t2



          
          sim%d  = tmp_d_new
          !immune and dead remain the same
          !temp = get_index(sim%t1,immune)
          !sim%t2(temp) = immune
          !temp = get_index(sim%t1,dead)
          sim%t2 = tmp%t2




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

       call end_timer("Sim Loop")

       timer = get_timer("Sim Loop")
       write(*,'(A)',ADVANCE="NO")"Time per day:"
       call write_realtime(frac_realtime(diff_realtimes(timer%rt_end,timer%rt_start),time_n))

       healthy_cases_final(:,:,it_ss) = healthy_cases
       inf_noncon_cases_final(:,:,it_ss) = inf_noncon_cases
       inf_contag_cases_final(:,:,it_ss) = inf_contag_cases
       ill_contag_cases_final(:,:,it_ss) = ill_contag_cases
       ill_ICU_cases_final(:,:,it_ss) =    ill_ICU_cases
       immune_cases_final(:,:,it_ss) =     immune_cases
       dead_cases_final(:,:,it_ss) =       dead_cases

    End Do     ! end do it_ss

    call start_timer("+- Writeout",reset=.FALSE.)
    
    iter_pass_handle = (/lhc(1,iter),lhc(2,iter),lhc(3,iter),&
         lhc(6,iter),lhc(7,iter),lhc(8,iter)/)

    Call write_data_v2(healthy_cases_final,iter_pass_handle,iol%R0_effect%data,counties_index,1)
    Call write_data_v2(inf_noncon_cases_final,iter_pass_handle,iol%R0_effect%data,counties_index,2)
    Call write_data_v2(inf_contag_cases_final,iter_pass_handle,iol%R0_effect%data,counties_index,3)
    Call write_data_v2(ill_contag_cases_final,iter_pass_handle,iol%R0_effect%data,counties_index,4)
    Call write_data_v2(ill_ICU_cases_final,iter_pass_handle,iol%R0_effect%data,counties_index,5)
    Call write_data_v2(immune_cases_final,iter_pass_handle,iol%R0_effect%data,counties_index,6)
    Call write_data_v2(dead_cases_final,iter_pass_handle,iol%R0_effect%data,counties_index,7)
!!$    End If
    call end_timer("+- Writeout")
1000 continue
!!$    Call MPI_Finalize(ierror)

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
    write(*,*)shape(R0change)
    
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



  !!============================================================================
  !> Model
  Subroutine Sim(time_n)

    Integer(Kind=ik), Intent(In) :: time_n
    
    Integer(Kind=ik)             :: timestep
    
    Do timestep = 2, time_n

       
       
    End Do
    
  End Subroutine Sim
  
End Module kernel
