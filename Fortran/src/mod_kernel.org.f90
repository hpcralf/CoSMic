module mod_kernel

contains

subroutine COVID19_Spatial_Microsimulation_for_Germany(iol,pspace,iter,counties_index)
use list_variable
use MOD_GLOBAL_VAR
use support_fun
!use csv_file
implicit none
include 'mpif.h'
type(pspaces)       :: pspace
type(iols)          :: iol
integer             :: i, j, k, index, temp_int,icounty,county,it_ss,iter,status
character*1         :: mod1,mod2
integer,dimension(:):: counties_index
logical             :: seed_in_inner_loop,seed_mpi_in_inner_loop
!!!!!-----1.states of the model ------
integer             :: healthy,inf_noncon,inf_contag,ill_contag,ill_ICU,immune,dead,missing

!!!!!-----2.derive the initial population ------
real                :: sam_prop_ps(16)=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/) ! scaling coeff for population
character*15        :: sim_pop

!!!!!-----3.Set seed infections in population ------
integer             :: ini_infected
character*5         :: seed_infections
integer             :: seed_before

!!!!!-----4. define disease characteristics -----
real                :: R0_force
!?? pspace
integer             :: inf_dur,cont_dur,ill_dur,icu_dur
real                :: icu_per_day(8)=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,8.0/)
real                :: less_contagious
logical             :: immune_stop

!!!!!-----5. define reductions in social contacts -----
integer             :: R0change(2,20)
character*5         :: R0county(20)
logical             :: R0delay
integer             :: R0delay_days
character*6         :: R0delay_type
logical             :: endogenous_lockdown
real                :: lockdown_connect
integer             :: lockdown_threshold
integer             :: lockdown_days

!!!!!-----7.  Define whether transition probabilities should differ by age and sex
character*10        :: control_age_sex
character*10        :: seed_before_char,seed_temp
character*10,allocatable :: seed_seq(:),seed_inf_cont_seq(:),seed_inf_ncont_seq(:)
character*10,allocatable :: seed_d_seq(:)
integer             :: days

integer             :: n_direct,n_directv,n_directl,n_dist
integer             :: size_lhc
real,allocatable    :: lhc(:,:)
character*10,allocatable :: lhc_name(:)
!!!!!-----8. variables for do loop -------
real,allocatable    :: icu_risk(:),surv_ill(:),surv_icu(:)
integer,allocatable :: temp1(:),temp(:),targert_icu_risk_index(:)
character*3,allocatable :: temp_char_mul(:)
character               :: temp_char_sig
type icu_risk_lists ! in order to have a similar structure to R
character*2,allocatable         :: age(:)
character*1,allocatable         :: sex(:)
real,allocatable                :: risk(:)
integer,allocatable             :: dur(:)
end type

type(icu_risk_lists)            :: icu_risk_list,surv_ill_list,surv_icu_list

type sims
integer,allocatable             :: dist_id(:)
character,allocatable           :: sex(:)
character*2,allocatable         :: age(:)
integer,allocatable             :: t1(:)
integer,allocatable             :: t2(:)
integer,allocatable             :: d(:)
end type

type(sims)                      :: sim,tmp
integer,allocatable             :: tmp_dnew(:)
integer,allocatable             :: sim_counties(:),rownumbers(:)

type(seeds)                     :: seed_ill,seed_inf_cont,seed_inf_ncont,target_inf
integer,allocatable             :: seed_ill_dur(:),seed_inf_cont_dur(:),seed_inf_ncont_dur(:)
integer,allocatable             :: il_d(:),inf_c_d(:),inf_nc_d(:)
integer,allocatable             :: rownumbers_ill(:),rownumbers_cont(:),rownumbers_ncont(:),rownumbers_dea(:),&
                                   rownumbers_left(:)
integer,allocatable             :: gettime(:)
real,allocatable                :: getchange(:)
integer                         :: inf_ill,inf_cont,inf_ncont,inf_dth

real                            :: R0_daily
real,allocatable                :: R0matrix(:,:),connect(:,:)
integer,allocatable             :: healthy_cases_final(:,:,:),&
                                   ill_ICU_cases_final(:,:,:),immune_cases_final(:,:,:),&
                                   inf_noncon_cases_final(:,:,:),inf_contag_cases_final(:,:,:),&
                                   dead_cases_final(:,:,:),ill_contag_cases_final(:,:,:)
integer,allocatable             :: inf_cases(:,:),icu_cases(:,:),healthy_cases(:,:),inf_noncon_cases(:,:),&
                                   inf_contag_cases(:,:),ill_contag_cases(:,:)
integer,allocatable             :: ill_ICU_cases(:,:),immune_cases(:,:),dead_cases(:,:),dead_cases_bICU(:,:),&
                                   mod_inf_cases(:,:),org_noncon_cases(:,:)
integer,allocatable             :: lockdowns(:),start_value_tot(:)
integer                         :: timestep

integer,allocatable             :: tmp_index(:),tmp_d_new(:),tmp_count(:)
integer,allocatable             :: susceptible(:),contagious_dist_id(:),contagious_index(:),denominator(:),&
                                   revers_proj(:),final_count(:),dist_id_temp(:),ill_index(:),&
                                   ill_dist_id(:)
real,allocatable                :: contagious(:)
integer                         :: at_risk
integer                         :: initial_sick
real                            :: n_contagious,between_weight,within_weight
real,allocatable                :: exp_infect(:)
integer,allocatable             :: check_days(:)
logical                         :: flag_lockdown
real,allocatable                :: risk(:),prob(:),runif(:),prop_target(:)
character*10                    :: target_date,temp_date

integer,allocatable             :: sick(:),case_count(:),state_id(:),prop_inf_cases(:),&
                                    prop_target_inf(:),mod_inf(:)
character*6,allocatable         :: age_sex(:),surv_ill_label(:),surv_icu_label(:)
character*6,allocatable         :: age_sex_dur(:),icu_risk_label(:)
character*2,allocatable         :: temp_character(:)
real,allocatable                :: surv_ill_i(:),die(:),icu_risk_i(:),icu(:),surv_icu_i(:)
integer,allocatable             :: die_count(:),icu_count(:),die_icu_count(:)
character*2,allocatable         :: ch_age(:)
character*5,allocatable         :: ch_sex(:)
integer                         :: max_date,n_change
character*10                    :: temp_mod
integer                         :: ierror,size_of_process,my_rank
integer                         :: index_final(7),block_size
integer,allocatable             :: req(:)

real                            :: iter_pass_handle(6)
! should import some reliable romdon seed generation code here
!seed_base = ??

index_final = 1
call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,size_of_process,ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierror)

allocate(req(size_of_process))

seed_in_inner_loop = .FALSE.
seed_mpi_in_inner_loop = .TRUE.

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
if (size(icu_per_day) /= ill_dur) then
    print *,'Length icu_per_day not equal to ill_dur'
endif
if ((sum(icu_per_day)/size(icu_per_day)) /= 1) then
    print *,'Mean icu per day not equal to 1'
end if

less_contagious = 0.7

R0_force = 0
immune_stop = .TRUE.

!!!!!-----5. define reductions in social contacts -----
! R0change = reshape((/1,7,8,13,14,19,20,25,26,32,33,39,40,46,47,53,54,60,61,67,68,74,75,81,82,88,&
!                     89,95,96,102,103,109,110,116,117,123,124,130,131,137,138,144,145,151,152,158,&
!                     159,164/),shape(R0change))
R0change = reshape((/1,7,8,13,14,19,20,25,26,32,33,39,40,46,47,53,54,60,61,67,68,74,75,81,82,88,&
                    89,95,96,102,103,109,110,116,117,123,124,130,131,137/),shape(R0change))
R0county            = "ALL"
R0delay             = .TRUE.
R0delay_days        = 5
R0delay_type        = "linear"
endogenous_lockdown = .FALSE.
lockdown_connect    = 0.5
lockdown_threshold  = 100
lockdown_days       = 10

!!!!!-----7.  Define whether transition probabilities should differ by age and sex
control_age_sex     = "age"
days             = 1
seed_date        = add_date(seed_date,days)
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

if (maxval(R0change)> time_n) then
    time_n = maxval(R0change) + 1
end if

! this part should done by the code, but here is set manually for simplicity
n_direct  = 8
n_directv = 0
n_directl = 1
n_dist    = 0

if (n_dist > 0) then
    ! code for randomLHS
    ! since it would not affect the code, it can be apllied later
else
    size_lhc = 0                    ! the first position is reserved for sam_size
end if

if (n_direct > 0) then
    size_lhc = size_lhc + n_direct
end if

if (n_directl > 0) then
    size_lhc = size_lhc + size(pspace%ROeffect_ps%param)
end if

allocate(lhc(size_lhc,iter))

! lhc(1,:)            = pspace%sam_size%param
do i = 1,n_direct
    lhc(i,:) = pspace%Ps_scalar_list(i)%param
end do


do i = 1,iter
lhc(n_direct+1:size(lhc,dim=1),i) = reshape(pspace%ROeffect_ps%param,shape(lhc(n_direct+1:size(lhc,dim=1),1)))
end do
! print *, "after reshape is",reshape(pspace%ROeffect_ps%param,shape(lhc(n_direct+1:size(lhc,dim=1),1)))
allocate(lhc_name(size_lhc))

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

if (control_age_sex == "NONE") then
    ch_age = (/"total"/)
    ch_sex = (/"total"/)
end if

if (control_age_sex == "age") then
!     ch_age = generate_seq(0,90,5)
    ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
    "40","45","50","55","60",&
    "65","70","75","80","85","90"/)
    ch_sex = (/"total"/)
end if

if (control_age_sex == "sex") then
    ch_age = (/"total"/)
    ch_sex = (/"m","f"/)
end if

if (control_age_sex == "age_sex") then
!     ch_age = generate_seq(0,90,5)
    ch_age = (/"0 ","5 ","10","15","20","25","30","35",&
    "40","45","50","55","60",&
    "65","70","75","80","85","90"/)
    ch_sex = (/"m","f"/)
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! the following part was in the loop at R code
!!! in order to aviod to raise problem of memery
!!! allocation and redundent calculation, put it
!!! outside the loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp = get_index(iol%transpr_age_gr,ch_age)
temp1= get_index(iol%transpr_sex,ch_sex)
targert_icu_risk_index = find_and(temp,temp1)

if (control_age_sex == "age") then
        if (.not.allocated(icu_risk))then
            allocate(icu_risk(2*size(targert_icu_risk_index)))
        endif
        icu_risk(1:size(targert_icu_risk_index)) = iol%transpr_icu_risk(targert_icu_risk_index)
        icu_risk(1:size(targert_icu_risk_index)) = 1.0 - (1.0 - icu_risk(1:size(targert_icu_risk_index))) ** (1.0/ real(ill_dur))
        icu_risk(size(targert_icu_risk_index)+1 : 2* size(targert_icu_risk_index)) = icu_risk(1:size(targert_icu_risk_index))
end if

if (control_age_sex == "sex") then
        if (.not.allocated(icu_risk))then
            allocate(icu_risk(2*19))
        end if
        icu_risk(1:2) = iol%transpr_icu_risk(targert_icu_risk_index)
        icu_risk(1:2) = 1.0 - (1.0 - icu_risk(1:2)) ** (1.0/ ill_dur)
        icu_risk(size(targert_icu_risk_index)+1 : 2* size(targert_icu_risk_index)) = icu_risk(2)
        icu_risk(1:size(targert_icu_risk_index)) = icu_risk(1)
end if
!     skip line 960 to 940 in R code
if (.not.allocated(icu_risk_list%age))then
        allocate(icu_risk_list%age(size(icu_risk)*ill_dur))
        allocate(icu_risk_list%sex(size(icu_risk)*ill_dur))
        allocate(icu_risk_list%risk(size(icu_risk)*ill_dur))
        allocate(icu_risk_list%dur(size(icu_risk)*ill_dur))
end if

do i = 1, ill_dur
        icu_risk_list%age((i-1)*size(icu_risk)+1: i*size(icu_risk)) = (/iol%transpr_age_gr(targert_icu_risk_index),&
                                                                    iol%transpr_age_gr(targert_icu_risk_index)/)
        icu_risk_list%sex((i-1)*size(icu_risk)+1:(i-1)*size(icu_risk)+19) = 'm'
        icu_risk_list%sex((i-1)*size(icu_risk)+20:i*size(icu_risk)) = 'f'
        icu_risk_list%risk((i-1)*size(icu_risk)+1: i*size(icu_risk)) = icu_risk * icu_per_day(i)
        icu_risk_list%dur((i-1)*size(icu_risk)+1: i*size(icu_risk)) = i
end do

    ! init surv_ill
if (control_age_sex == "age") then
        if (.not.allocated(surv_ill))then
            allocate(surv_ill(2*size(targert_icu_risk_index)))
        end if
        surv_ill = (/iol%transpr_surv_ill(targert_icu_risk_index),&
                   iol%transpr_surv_ill(targert_icu_risk_index)/)
        surv_ill = surv_ill ** (1.0/real(ill_dur))
end if

if (control_age_sex == "sex") then
        if (.not.allocated(surv_ill))then
            allocate(surv_ill(2*19))
        endif
        surv_ill(1:2) = iol%transpr_surv_ill(targert_icu_risk_index)
        surv_ill(1:2) = 1.0 - (1.0 - surv_ill(1:2)) ** (1.0/ real(ill_dur))
        surv_ill(size(targert_icu_risk_index)+1 : 2* size(targert_icu_risk_index)) = surv_ill(2)
        surv_ill(1:size(targert_icu_risk_index)) = surv_ill(1)
end if
if (.not.allocated(surv_ill_list%age))then
        allocate(surv_ill_list%age(size(surv_ill)))
        allocate(surv_ill_list%sex(size(surv_ill)))
        allocate(surv_ill_list%risk(size(surv_ill)))
end if
surv_ill_list%age = (/iol%transpr_age_gr(targert_icu_risk_index),iol%transpr_age_gr(targert_icu_risk_index)/)
surv_ill_list%sex = 'f'
surv_ill_list%sex(1:19) = 'm'
surv_ill_list%risk  = surv_ill

    ! init surv_icu
if (control_age_sex == "age") then
        if (.not.allocated(surv_icu))then
            allocate(surv_icu(2*size(targert_icu_risk_index)))
        end if
        surv_icu= (/iol%transpr_surv_icu(targert_icu_risk_index),&
                   iol%transpr_surv_icu(targert_icu_risk_index)/)
        surv_icu = surv_icu ** (1/real(ill_dur))
end if

if (control_age_sex == "sex") then
        if (.not.allocated(surv_icu))then
            allocate(surv_icu(2*19))
        end if
        surv_icu(1:2) = iol%transpr_surv_icu(targert_icu_risk_index)
        surv_icu(1:2) = 1.0 - (1.0 - surv_icu(1:2)) ** (1.0/ real(ill_dur))
        surv_icu(size(targert_icu_risk_index)+1 : 2* size(targert_icu_risk_index)) = surv_icu(2)
        surv_icu(1:size(targert_icu_risk_index)) = surv_icu(1)
end if
if (.not.allocated(surv_icu_list%age))then
        allocate(surv_icu_list%age(size(surv_icu)))
        allocate(surv_icu_list%sex(size(surv_icu)))
        allocate(surv_icu_list%risk(size(surv_icu)))
end if
surv_icu_list%age = (/iol%transpr_age_gr(targert_icu_risk_index),iol%transpr_age_gr(targert_icu_risk_index)/)
surv_icu_list%sex = 'f'
surv_icu_list%sex(1:19) = 'm'
surv_icu_list%risk  = surv_icu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! the following part was in the loop at R code
!!! put it outside to avoid memoery problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (my_rank == 0)then

allocate(healthy_cases_final(size(counties_index),time_n,iter))
allocate(inf_noncon_cases_final(size(counties_index),time_n,iter))
allocate(inf_contag_cases_final(size(counties_index),time_n,iter))
allocate(ill_contag_cases_final(size(counties_index),time_n,iter))
allocate(ill_ICU_cases_final(size(counties_index),time_n,iter))
allocate(immune_cases_final(size(counties_index),time_n,iter))
allocate(dead_cases_final(size(counties_index),time_n,iter))

allocate(inf_cases(size(counties_index),time_n))
allocate(icu_cases(size(counties_index),time_n))

allocate(healthy_cases(size(counties_index),time_n))
allocate(inf_noncon_cases(size(counties_index),time_n))
allocate(inf_contag_cases(size(counties_index),time_n))
allocate(ill_contag_cases(size(counties_index),time_n))
allocate(ill_ICU_cases(size(counties_index),time_n))
allocate(immune_cases(size(counties_index),time_n))
allocate(dead_cases(size(counties_index),time_n))
allocate(dead_cases_bICU(size(counties_index),time_n))
allocate(mod_inf_cases(size(counties_index),time_n))
allocate(org_noncon_cases(size(counties_index),time_n))

else
allocate(inf_cases(size(counties_index),time_n))
allocate(icu_cases(size(counties_index),time_n))

allocate(healthy_cases(size(counties_index),time_n))
allocate(inf_noncon_cases(size(counties_index),time_n))
allocate(inf_contag_cases(size(counties_index),time_n))
allocate(ill_contag_cases(size(counties_index),time_n))
allocate(ill_ICU_cases(size(counties_index),time_n))
allocate(immune_cases(size(counties_index),time_n))
allocate(dead_cases(size(counties_index),time_n))
allocate(dead_cases_bICU(size(counties_index),time_n))
allocate(mod_inf_cases(size(counties_index),time_n))
allocate(org_noncon_cases(size(counties_index),time_n))

endif
allocate(lockdowns(time_n))
block_size = size(counties_index) * time_n

max_date = find_max_date(iol%seed_date)

!!$OMP PARALLEL DO private(sim)
if (my_rank/=0)then
do it_ss = my_rank,size(lhc,dim=2),size_of_process-1
    call random_seed()
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
    if (sim_pop == "proportional") then
        iol%pop_total = nint(real(iol%pop_total)/real(sum(iol%pop_total))* real(lhc(1,it_ss)))
        if (.not.allocated(sim%dist_id))then
            allocate(sim%dist_id(sum(iol%pop_total)))
            allocate(sim%sex(sum(iol%pop_total)))
            allocate(sim%age(sum(iol%pop_total)))
            allocate(sim%t1(sum(iol%pop_total)))
            allocate(sim%t2(sum(iol%pop_total)))
            allocate(sim%d(sum(iol%pop_total)))
        endif

        index = 0 ! position index
        ! set all male for testing

        do i = 1, size(iol%pop_total)
            temp_int = iol%pop_total(i)
            sim%dist_id(index+1: index+temp_int) = iol%pop_distid(i)
            sim%sex(index+1: index+temp_int) = iol%pop_sex(i)
            sim%age(index+1: index+temp_int) = iol%pop_age(i)
            index                            = index + temp_int
        end do

        sim%t1 = healthy
        sim%t2 = missing
        sim%d(:)  = 1
    end if
    if (my_rank == 1)then
        print *,"size of population is",size(sim%d)
    end if

    if (sim_pop == "random") then
        ! can be filled later
    end if

    !!!!!!!!!!!!!
    !seed infections
    !!!!!!!!!!!!!
    if (seed_infections == "random") then
        ! remains blank
    end if

    if (seed_infections == "county") then
        ! remains blank
    end if

    if (seed_infections == "data") then
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

        if (.not.allocated(seed_ill_dur))then
            allocate(seed_ill_dur(size(seed_ill%date)))
        end if
        temp_int = Date2Unixtime(seed_date)

        do i = 1,size(seed_ill_dur)
            seed_ill_dur(i)  = (temp_int - Date2Unixtime(seed_ill%date(i)))/86400 + 1
        end do
        temp_int = 0
!         print *,seed_ill%cases
!         mod1 = "l"
!         mod2 = "g"
        temp = condition_and(seed_ill_dur,ill_dur+1,"l",seed_ill%cases,temp_int,"g")
!          print *,"size of temp",size(temp)
        seed_ill%dist_id  = seed_ill%dist_id(temp)
        seed_ill%date     = seed_ill%date(temp)
        seed_ill%cases    = seed_ill%cases(temp)

        if (.not.allocated(seed_inf_cont_dur))then
            allocate(seed_inf_cont_dur(size(seed_inf_cont%date)))
        end if
        temp_int = Date2Unixtime(seed_date)
        do i = 1,size(seed_inf_cont_dur)
            seed_inf_cont_dur(i)  = (temp_int - Date2Unixtime(seed_inf_cont%date(i)))/86400&
                                    + cont_dur + 1
        end do

        deallocate(temp)
        allocate(temp(size(seed_inf_cont%cases)))

        temp = 0
        temp_int = 0
        do i = 1,size(temp)
            if (seed_inf_cont%cases(i)>0) then
                temp_int = temp_int + 1
                temp(temp_int) = i
            end if
        end do

!         seed_inf_cont = set_value(seed_inf_cont,temp(1:temp_int))
        seed_inf_cont%dist_id  = seed_inf_cont%dist_id(temp(1:temp_int))
        seed_inf_cont%date     = seed_inf_cont%date(temp(1:temp_int))
        seed_inf_cont%cases    = seed_inf_cont%cases(temp(1:temp_int))
        deallocate(temp)

        if (.not.allocated(seed_inf_ncont_dur))then
            allocate(seed_inf_ncont_dur(size(seed_inf_ncont%date)))
        else
            deallocate(seed_inf_ncont_dur)
            allocate(seed_inf_ncont_dur(size(seed_inf_ncont%date)))
        end if
        temp_int = Date2Unixtime(seed_date)
        do i = 1,size(seed_inf_ncont_dur)
            seed_inf_ncont_dur(i)  = (temp_int - Date2Unixtime(seed_inf_ncont%date(i)))/86400&
                                     + inf_dur + cont_dur + 1
        end do
        allocate(temp(size(seed_inf_ncont%cases)))
        temp = 0
        temp_int = 0
        do i = 1,size(temp)
            if (seed_inf_ncont%cases(i)>0) then
                temp_int = temp_int + 1
                temp(temp_int) = i
            end if
        end do
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

        do icounty = 1,size(counties_index)

            county = counties_index(icounty)

            rownumbers = get_index(sim%dist_id,county)


            temp   = get_index(seed_ill%dist_id,county)
            il_d   = rep(seed_ill_dur(temp),seed_ill%cases(temp))
            inf_ill= sum(seed_ill%cases(temp))

            temp   = get_index(seed_inf_cont%dist_id,county)
            inf_c_d= rep(seed_inf_cont_dur(temp),seed_inf_cont%cases(temp))

            inf_cont = sum(seed_inf_cont%cases(temp))

            temp   = get_index(seed_inf_ncont%dist_id,county)
            inf_nc_d = rep(seed_inf_ncont_dur(temp),seed_inf_ncont%cases(temp))
            inf_ncont = sum(seed_inf_ncont%cases(temp))

            temp   = get_index(iol%death_distid,county)
            inf_dth = sum(iol%death_cases(temp))

            if(size(rownumbers)<(inf_ill+inf_cont+inf_ncont+inf_dth)) then
                print *,"Number of infected and dead is larger than population size"
                print *,"only ",size(rownumbers),"number left"
                print *,"total cases is",(inf_ill+inf_cont+inf_ncont+inf_dth)
                print *,"seperate cases are:" ,inf_ill,inf_cont,inf_ncont,inf_dth
                print *,"timestep is",timestep
                print *,"it_ss is",it_ss
                print *,"county is",county
                print *,"counties are",counties_index
                print *,"icounty is",icounty
                call exit(status)
            end if

            rownumbers_left = rownumbers
            if (inf_ill > 0) then

                rownumbers_ill = sample(rownumbers_left,inf_ill)
                rownumbers_left = rownumbers_left(inf_ill+1:size(rownumbers))
                sim%t1(rownumbers_ill) = ill_contag
                sim%d(rownumbers_ill)  = il_d
            end if

            if ( inf_cont > 0) then
                rownumbers_cont = sample(rownumbers_left,inf_cont)
                rownumbers_left = rownumbers_left(inf_cont+1:size(rownumbers_left))
                sim%t1(rownumbers_cont)= inf_contag
                sim%d(rownumbers_cont) = inf_c_d
            end if

            if (inf_ncont > 0) then
                rownumbers_ncont = sample(rownumbers_left,inf_ncont)
                rownumbers_left = rownumbers_left(inf_ncont+1:size(rownumbers_left))
                sim%t1(rownumbers_ncont) = inf_noncon
                sim%d(rownumbers_ncont)  = inf_nc_d
            end if

            if (inf_dth > 0) then
                rownumbers_dea = sample(rownumbers_left,inf_dth)
                sim%t1(rownumbers_dea) = dead
            end if
        end do ! do icounty = 1,size(sim_counties)
    end if ! if (seed_infections == "data")

    if (import_R0_matrix) then
!         temp = get_index(iol%R0_raw)
!         R0_raw1 =
    end if
    R0_daily = R0_force *lhc(2,it_ss)/real(real(cont_dur)+real(ill_dur)*less_contagious) + &
               (1-R0_force)*lhc(2,it_ss)/real(cont_dur+ill_dur)
    ! this block simplifies the if judgment
    if (.not.allocated(R0matrix))then
        allocate(R0matrix(size(counties_index),time_n-1))
    end if
    R0matrix = R0_daily
    n_change  = size(R0change,dim=2)

    do i = 1,n_change
        if(R0county(i) == "ALL") then
            !use counties
            ! give up using character to find the number, but use the index directly
            ! if is "all" use full counties_index
            temp  = (i-1) * size(pspace%ROeffect_ps%param,dim= 1) + 9 + ((counties_index/1000)-1)

            getchange = lhc(temp,it_ss)

            gettime = generate_seq(R0change(1,i),R0change(2,i),1)
            do j = 1,size(gettime)
                R0matrix(:,gettime(j)) = R0matrix(:,gettime(j)) * getchange
            end do
        end if
    end do

    if (R0delay) then
        do i = 1,size(R0matrix,dim = 1)
            R0matrix(i,:) = smoothing_change(R0matrix(i,:),R0delay_days,R0delay_type)
        end do
    end if

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
    if (.not.allocated(tmp_d_new))then
        allocate(tmp_d_new(size(sim%d)))
    end if

    do timestep = 2,time_n
        !skip using prev_state and up_state
        tmp = sim

!         print *,"now in the ",timestep, " th loop"

        tmp_index = get_index(tmp%t1,(/inf_contag,inf_noncon,ill_ICU,ill_contag,healthy/))

        tmp_d_new = missing

        if (size(tmp_index) == 0)then
            ! leave it blank and deal with it later
        end if
        tmp_count = get_index(tmp%t1,healthy)
!         susceptible = sum_byindex(tmp%dist_id(tmp_count),counties_index)

        at_risk = size(tmp_count)

        contagious_index   = get_index(tmp%t1,inf_contag)
        ill_index          = get_index(tmp%t1,ill_contag)

        contagious_dist_id = get_unique(tmp%dist_id(contagious_index))
!         if(my_rank == 1)then
!              print *,size(contagious_dist_id)
!         endif
!         print *,"ill case is",size(ill_index)
        ill_dist_id        = get_unique(tmp%dist_id(ill_index))
        if (.not.allocated(contagious))then
            allocate(contagious(size(counties_index)))
        end if
!         do i = 1,size(counties_index)
!             temp1 = get_index(tmp%dist_id(contagious_index),counties_index(i))
!             contagious(i) = real(size(temp1))
!             temp1 = get_index(tmp%dist_id(ill_index),counties_index(i))
!             contagious(i) = contagious(i) + real(size(temp1)) * less_contagious
!         end do
        contagious = real(sum_byindex(tmp%dist_id(contagious_index),counties_index))
!         if(my_rank == 1)then
!              print *,size(contagious)
!         endif
        contagious = contagious + real(sum_byindex(tmp%dist_id(ill_index),counties_index)) * less_contagious
!         print *,"ill cases are: ",real(size(temp1))
!         print *,"nint of nint(sum(contagious)) is",nint(sum(contagious))
!         print *,"sum(susceptible) is ",sum(susceptible)
         n_contagious = sum(contagious)
         temp = get_index(tmp%t1,inf_noncon)
!         if(my_rank == 1)then
!             print *,at_risk
!         endif
        if (at_risk > 0 .and. n_contagious > 0) then
            exp_infect = contagious * R0matrix(:,timestep-1)
!             if(my_rank == 1)then
!                 print *,at_risk
!             endif
!             print *,"contagious is",contagious
            if (endogenous_lockdown) then
                temp = generate_seq(timestep-1,timestep-lockdown_days,-1)
                check_days = pack(temp,temp>0)
                temp = sum(inf_cases(:,check_days),dim = 1)
                flag_lockdown  = .FALSE.
                temp = pack(temp,temp>lockdown_threshold)
                if (size(temp)>0)then
                    flag_lockdown = .TRUE.
                end if
                if (flag_lockdown)then
                    exp_infect = exp_infect *lhc(5,it_ss)! times lock_downeffect, which is the 5th element of the array
                    lockdowns(timestep) = 1
                else
                    lockdowns(timestep) = 0
                end if

            else
                lockdowns(timestep) = 0
            end if ! end if (endogenous_lockdown)
            temp = get_index(iol%connect_work_distid,counties_index)

            connect = iol%connect_work(temp,temp)
            do i = 1,size(contagious_dist_id)
                connect(:,i) = connect(:,i)/sum(connect(:,i))
            end do

!             connect = transpose(connect)
            between_weight = 1 - lhc(6,it_ss) ! 1 - lhc(it_ss,"w_int")

            if (endogenous_lockdown)then
                between_weight = between_weight * lockdown_connect
            end if

            within_weight = 1 - between_weight
            exp_infect = within_weight * exp_infect + between_weight * matmul (exp_infect, connect)!matmul (connect,exp_infect)
!             do i = 1,size(contagious_dist_id)
!                 exp_infect(i) = within_weight * exp_infect(i) + between_weight * sum(exp_infect*connect(:,i))
!             end do
            if (.not.allocated(denominator))then
                allocate(denominator(size(counties_index)))
            end if
            if (.not.immune_stop) then
!                 denominator = sum_bygroup(tmp%dist_id(temp))
!                 do i = 1,size(counties_index)
!                     temp1 = get_index(tmp%dist_id(temp),counties_index(i))
!                     denominator(i) = size(temp1)
!                 end do
                denominator = sum_byindex(tmp%dist_id(tmp_count),counties_index)
            else
                temp = get_index(tmp%t1,(/healthy,immune,inf_noncon,inf_contag,ill_contag/))
!                 do i = 1,size(counties_index)
!                     temp1 = get_index(tmp%dist_id(temp),counties_index(i))
!                     denominator(i) = size(temp1)
!                 end do
                denominator = sum_byindex(tmp%dist_id(temp),counties_index)
            end if

            risk = exp_infect/real(denominator)
            where (risk > 1000)
                risk = 0
            end where
            where (risk > 1)
                risk = 1
            end where

            risk = lhc(6,it_ss)*risk + (1-lhc(6,it_ss))* matmul(risk,connect)!matmul (connect,risk)


            temp = get_index(counties_index,tmp%dist_id(tmp_count))
            prob = risk(temp)

            if (.not.allocated(runif))then
                allocate(runif(at_risk))
                allocate(sick(at_risk))
            else
                deallocate(runif)
                deallocate(sick)
                allocate(runif(at_risk))
                allocate(sick(at_risk))
            end if
            call random_number(runif)
            !check if probs in runif is higher than
            temp_int = 0
            sick = 0
            do i = 1,size(runif)
                if (runif(i) <= prob(i)) then
                    temp_int = temp_int + 1
                    sick(i)  = 1
                end if
            end do

            if (allocated(revers_proj))then
                deallocate(revers_proj)
                allocate(revers_proj(temp_int))
            else
                allocate(revers_proj(temp_int))
            endif
            temp_int = 0

            do i = 1,size(runif)
                if (runif(i) <= prob(i)) then
                    temp_int = temp_int + 1
                    revers_proj(temp_int) = tmp_count(i) ! record index of infected person in healthy group
                end if
            end do

            initial_sick = size(sick)

            if (.not.allocated(final_count))then
                allocate(final_count(size(counties_index)))
!                 deallocate(dist_id_temp)
            else
                deallocate(final_count)
                allocate(final_count(size(counties_index)))
            end if
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

            if (lhc(7,it_ss) > 0 .and. Date2Unixtime(target_date) < max_date .and. initial_sick > 0)then
!                print *,"i'm in"
                temp = get_index(iol%seed_date,target_date)
                target_inf%cases   = iol%seed_cases(temp)
                target_inf%dist_id = iol%seed_distid(temp)
                target_inf%date    = iol%seed_date(temp)
                if (.not.allocated(final_count))then
                    allocate(final_count(size(counties_index)))
                else
                    deallocate(final_count)
                    allocate(final_count(size(counties_index)))
                end if
                final_count = 0
!                 do i =1,size(counties_index)
!                     temp = get_index(target_inf%dist_id,counties_index(i))
!                     final_count(i) = size(temp)
!                 end do
                final_count = sum_byindex(target_inf%dist_id,counties_index)
                if (sum(final_count) > 0) then
                    if(lhc(8,it_ss) == 1) then
                        state_id = get_unique(counties_index/1000)
                        if (.not.allocated(mod_inf))then
                            allocate(mod_inf(size(counties_index)))
                        end if
                        mod_inf = 0
                        do i = 1, size(state_id)
                            temp = get_index(counties_index/1000,state_id(i))
                            prop_inf_cases = inf_cases(temp,timestep)
                            temp_int       = sum(prop_inf_cases)
                            if (temp_int > 0) then
                                prop_inf_cases = prop_inf_cases / temp_int
                            else
                                prop_inf_cases = 0
                            end if

                            prop_target_inf = final_count(temp)
                            temp_int = sum(prop_target_inf)
                            if (temp_int > 0) then
                                prop_target_inf = prop_target_inf / temp_int
                            else
                                prop_target_inf = 0
                            end if
                            prop_target = lhc(7,it_ss) * real(prop_target_inf) + (1 - lhc(7,it_ss))*real(prop_inf_cases)

                            mod_inf(temp) = nint(prop_target* real(temp_int)) - inf_cases(temp,timestep)
                        end do

                    else
                        prop_inf_cases = inf_cases(:,timestep) / initial_sick
                        prop_target = lhc(7,it_ss) * real(final_count)/real(sum(final_count)) + &
                                      (1 - lhc(7,it_ss))*real(prop_inf_cases)
                        ! line 1675
                        mod_inf = nint(prop_target * real(initial_sick))-inf_cases(:,timestep)

                    end if

                    call shift_cases(sick,mod_inf,tmp%dist_id(tmp_count),counties_index)
                    !!!!!!!!!!!!!!!!!!!!!!!!
                    !!!can we write the following code in a more elgant
                    !!!way?
                    !!!!!!!!!!!!!!!!!!!!!!!!

                    temp = get_index(sick,1)
                    ! do a revers projection
                    temp = tmp_count(temp)
                    if (.not.allocated(final_count))then
                        allocate(final_count(size(counties_index)))
                    end if
                    if (allocated(dist_id_temp))then
                        deallocate(dist_id_temp)
                    end if
                    if (allocated(case_count))then
                        deallocate(case_count)
                    end if

!                     call sum_bygroup_distID(tmp%dist_id(temp),case_count,dist_id_temp)
                    final_count = 0
!                     do i =1,size(counties_index)
!                         temp1 = get_index(tmp%dist_id(temp),counties_index(i))
!                         final_count(i) = size(temp1)
!                     enddo
                    final_count = sum_byindex(tmp%dist_id(temp),counties_index)
                    inf_cases(:,timestep) = final_count

                    mod_inf_cases(:,timestep) = mod_inf
                else
                    print *,"no cases shifted"
                end if

            else
                mod_inf_cases(:,timestep) = 0
            endif

            !skip 1742 and 1743 since they are not actually do
            !anythin?

            tmp%t2(tmp_count) = sick

            temp = condition_and(tmp%t1,healthy,"e",tmp%t2,inf_noncon,"e")


            tmp_d_new(temp) = 1
!             tmp%d(temp) = 1

            temp = condition_and(tmp%t1,healthy,"e",tmp%t2,healthy,"e")
            tmp_d_new(temp) = tmp%d(temp) + 1
!             tmp%d(temp) =    tmp%d(temp) + 1
        end if

        if (at_risk > 0 .and. nint(n_contagious) == 0) then
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

        end if

        if (at_risk == 0) then
            inf_cases(:,timestep) = 0
        end if


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
        at_risk = size(temp)

        if (at_risk > 0) then
            age_sex = tmp%age(temp) // tmp%sex(temp)
            surv_ill_label = surv_ill_list%age //surv_ill_list%sex
            temp = get_index(surv_ill_label,age_sex)
            surv_ill_i = surv_ill_list%risk(temp)
            if (.not.allocated(die)) then
                allocate(die(size(temp)))
                allocate(die_count(size(temp)))
            else
                deallocate(die)
                deallocate(die_count)
                allocate(die(size(temp)))
                allocate(die_count(size(temp)))
            end if
            die_count = ill_contag
            call random_number(die)
            do i = 1,size(die)
                if (surv_ill_i(i) <= die(i))then
                    die_count(i) = dead
                endif
            end do
            temp = get_index(tmp%t1,ill_contag)

            tmp%t2(temp) = die_count

        endif

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

        at_risk = size(temp)


        
        if(at_risk > 0) then
            if (.not.allocated(temp_character))then
                allocate(temp_character(size(temp)))
            else
                deallocate(temp_character)
                allocate(temp_character(size(temp)))
            end if
            if (.not.allocated(age_sex_dur))then
                allocate(age_sex_dur(size(temp)))
            else
                deallocate(age_sex_dur)
                allocate(age_sex_dur(size(temp)))
            end if
            write(temp_character,'(I2)')tmp%d(temp)
            do i =1,size(temp)
                age_sex_dur(i) = trim(tmp%age(temp(i))) // trim(tmp%sex(temp(i)))//trim(temp_character(i))
            end do

            if (.not.allocated(temp_character))then
                allocate(temp_character(size(icu_risk_list%dur)))
            else
                deallocate(temp_character)
                allocate(temp_character(size(icu_risk_list%dur)))
            end if
            write(temp_character,'(I2)')icu_risk_list%dur
            if (.not.allocated(icu_risk_label))then
                allocate(icu_risk_label(size(icu_risk_list%age)))
            end if
            do i = 1,size(icu_risk_list%age)
                icu_risk_label(i) = trim(icu_risk_list%age(i))//trim(icu_risk_list%sex(i))&
                                    //trim(temp_character(i))
            end do

            temp1 = get_index(icu_risk_label,age_sex_dur)


            icu_risk_i = icu_risk_list%risk(temp1)

            if (.not.allocated(icu))then
                allocate(icu(size(temp)))
                allocate(icu_count(size(temp)))
            else
                deallocate(icu)
                deallocate(icu_count)
                allocate(icu(size(temp)))
                allocate(icu_count(size(temp)))
            endif
            icu_count = ill_contag
            call random_number(icu)
            do i = 1,size(icu_risk_i)
                if(icu(i) <= icu_risk_i(i))then
                    icu_count(i) = ill_ICU
                endif
            enddo
            
            tmp%t2(temp) = icu_count
            temp = condition_and(tmp%t1,ill_contag,"e",tmp%t2,ill_ICU,"e")
!             tmp%d(temp)  = 1
            tmp_d_new(temp) = 1
        end if

        if (at_risk == 0)then
            icu_cases(:,timestep) = 0
        endif

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
        at_risk = size(temp)

        if(at_risk > 0) then
            age_sex = tmp%age(temp) // tmp%sex(temp)
            surv_icu_label = surv_icu_list%age // surv_icu_list%sex
            temp1 = get_index(surv_icu_label,age_sex)

            surv_icu_i = surv_icu_list%risk(temp1)
            !who dies
            if (.not.allocated(die_icu_count))then
                allocate(die_icu_count(at_risk))
            else
                deallocate(die_icu_count)
                allocate(die_icu_count(at_risk))
            endif
            die_icu_count = ill_contag
            call random_number(die)
            do i = 1,size(surv_icu_i)
                if (die(i) >= surv_icu_i(i))then
                    die_icu_count(i) = dead
                end if
            end do
            tmp%t2(temp) = die_icu_count
        end if

!         tmp_count = get_index(tmp%t2,healthy)

        ! stay ill/ICU
        temp = condition_and(tmp%t1,ill_ICU,"e",tmp%d,int(lhc(3,it_ss)),"l")
        temp1= get_index(tmp%t2,(/healthy,ill_ICU,inf_noncon,inf_contag,ill_contag,immune/))
        temp = find_and(temp,temp1)
        tmp%t2(temp) = ill_ICU

        ! update days
        temp = condition_and(tmp%t1,ill_ICU,"e",tmp%t2,ill_ICU,"e")
!         tmp%d(temp) = tmp%d(temp) + 1
        tmp_d_new(temp) = tmp%d(temp) + 1

        !getting healthy /ICU if day limited reached
        temp = condition_and(tmp%t1,ill_ICU,"e",tmp%d,int(lhc(3,it_ss))-1,"g")
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




    end do ! end do timestep =2,time_n
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
    call MPI_ISend(healthy_cases(1,1),block_size,MPI_INTEGER,0,0,MPI_COMM_WORLD,req,ierror)
    call MPI_ISend(inf_noncon_cases(1,1),block_size,MPI_INTEGER,0,1,MPI_COMM_WORLD,req,ierror)
    call MPI_ISend(inf_contag_cases(1,1),block_size,MPI_INTEGER,0,2,MPI_COMM_WORLD,req,ierror)
    call MPI_ISend(ill_contag_cases(1,1),block_size,MPI_INTEGER,0,3,MPI_COMM_WORLD,req,ierror)
    call MPI_ISend(ill_ICU_cases(1,1),block_size,MPI_INTEGER,0,4,MPI_COMM_WORLD,req,ierror)
    call MPI_ISend(immune_cases(1,1),block_size,MPI_INTEGER,0,5,MPI_COMM_WORLD,req,ierror)
    call MPI_ISend(dead_cases(1,1),block_size,MPI_INTEGER,0,6,MPI_COMM_WORLD,req,ierror)
end do     ! end do it_ss
else ! the very begining if:my_rank /= 0
    do k = 1,iter
        call MPI_Recv(healthy_cases_final(1,1,index_final(1)),block_size,MPI_INTEGER,&
                    MPI_ANY_SOURCE,0,MPI_COMM_WORLD,req,ierror)
        index_final(1) = index_final(1) + 1
        call MPI_Recv(inf_noncon_cases_final(1,1,index_final(2)),block_size,MPI_INTEGER,&
                    MPI_ANY_SOURCE,1,MPI_COMM_WORLD,req,ierror)  
        index_final(2) = index_final(2) + 1
        call MPI_Recv(inf_contag_cases_final(1,1,index_final(3)),block_size,MPI_INTEGER,&
                    MPI_ANY_SOURCE,2,MPI_COMM_WORLD,req,ierror)
        index_final(3) = index_final(3) + 1
        call MPI_Recv(ill_contag_cases_final(1,1,index_final(4)),block_size,MPI_INTEGER,&
                    MPI_ANY_SOURCE,3,MPI_COMM_WORLD,req,ierror)
        index_final(4) = index_final(4) + 1
        call MPI_Recv(ill_ICU_cases_final(1,1,index_final(5)),block_size,MPI_INTEGER,&
                    MPI_ANY_SOURCE,4,MPI_COMM_WORLD,req,ierror)
        index_final(5) = index_final(5) + 1
        call MPI_Recv(immune_cases_final(1,1,index_final(6)),block_size,MPI_INTEGER,&
                    MPI_ANY_SOURCE,5,MPI_COMM_WORLD,req,ierror)
        index_final(6) = index_final(6) + 1
        call MPI_Recv(dead_cases_final(1,1,index_final(7)),block_size,MPI_INTEGER,&
                    MPI_ANY_SOURCE,6,MPI_COMM_WORLD,req,ierror)
        index_final(7) = index_final(7) + 1
    end do
end if
                  
!!$OMP END PARALLEL DO
call MPI_Finalize(ierror)
!call write_data()


if (my_rank == 0)then
    iter_pass_handle = (/lhc(1,iter),lhc(2,iter),lhc(3,iter),&
                         lhc(6,iter),lhc(7,iter),lhc(8,iter)/)
!     do k = 1,iter
!         print *,healthy_cases_final(:,:,k)
!     end do
   call write_data_v2(healthy_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,1)
   call write_data_v2(inf_noncon_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,2)
   call write_data_v2(inf_contag_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,3)
   call write_data_v2(ill_contag_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,4)
   call write_data_v2(ill_ICU_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,5)
   call write_data_v2(immune_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,6)
   call write_data_v2(dead_cases_final,iter_pass_handle,pspace%ROeffect_ps%param,counties_index,7)
end if
end subroutine


subroutine write_data_v2(healthy_cases_final,iter_pass_handle,R0Change,counties_index,type_file)
real,dimension(:)        :: iter_pass_handle(:)
real,dimension(:,:)      :: R0Change
real,allocatable         :: R0change_rep(:,:),R0change_exp(:,:),temp_output(:)
integer,dimension(:,:,:) :: healthy_cases_final
integer,dimension(:)     :: counties_index
integer                  :: iter
integer                  :: type_file

integer :: county_size, count 
character*15                :: iter_char(6)
character*5,allocatable     :: R0change_name(:)
character*2                 :: counties(16)
integer,allocatable         :: counties_index_out(:),iter_array(:)

character*10,allocatable    :: Label(:)

integer date_time(8),i,j,k
character*10 b(3)
character*4 year
character*2 day,month
character*8 time
character*3 temp_char

character*15 dir_prefix

iter = size(healthy_cases_final,DIM=3)
county_size = size(healthy_cases_final,DIM=1)
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

allocate(R0change_rep(size(R0Change),1))
allocate(R0change_exp(iter*county_size,size(R0change_rep)))
allocate(temp_output(iter*county_size))
allocate(R0change_name(size(R0Change)))
allocate(iter_array(county_size*iter))
allocate(counties_index_out(county_size*iter))
allocate(Label(size(healthy_cases_final,dim=2)))

call date_and_time(b(1),b(2),b(3),date_time)

count = 1
do i = 1,size(counties)
    do j = 1,size(R0Change,2)
        write(temp_char,"(I2)")j
        if (j>9)then
            R0change_name(count) = counties(i)//temp_char
        else
            R0change_name(count) = counties(i)//adjustl(temp_char)
        end if
        count = count + 1
    end do
end do
count = 1
do i =1,size(healthy_cases_final,dim =2)
        write(Label(i),"(I3)")i
end do


write(year,"(I4)")date_time(1)
write(month,"(I2)")date_time(2)
write(day,"(I2)")date_time(3)
if (date_time(2)<=9)then
    month = "0"//adjustl(month)
endif


dir_prefix = "./output/"

if(date_time(3)<=9)then
    day = "0"//adjustl(day)
end if
time = year//month//day
if (type_file == 1)then
    open (101, file = trim(dir_prefix)//time//"healthy_cases.csv")
end if

if (type_file == 2)then
    open (101, file = trim(dir_prefix)//time//"inf_noncon_cases.csv")
end if

if (type_file == 3)then
    open (101, file = trim(dir_prefix)//time//"inf_contag_cases.csv")
end if

if (type_file == 4)then
    open (101, file = trim(dir_prefix)//time//"ill_contag_cases.csv")
end if

if (type_file == 5)then
    open (101, file = trim(dir_prefix)//time//"ill_ICU_cases.csv")
end if

if (type_file == 6)then
    open (101, file = trim(dir_prefix)//time//"immune_cases.csv")
end if

if (type_file == 7)then
    open (101, file = trim(dir_prefix)//time//"dead_cases.csv")
end if

10 format(1x,*(g0,","))

!write the head file
write(101,10)Label,iter_char,R0change_name,"iter","x.dist_id"

R0change_rep = reshape(R0Change,(/size(R0change),1/))
R0change_exp = transpose(spread(R0change_rep(:,1),2,iter*county_size))

count = 1
do i = 1,iter
    do j = 1,county_size
        write(101,10)healthy_cases_final(j,:,i),iter_pass_handle,R0change_exp(count,:),&
                     i,counties_index(j)
        count = count+1
    end do
end do

close(101)
end subroutine

end module
