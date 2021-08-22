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
! Module with Input/Output routines
!
!###############################################################################
module CoSMic_IO

  Use global_constants
  Use global_types
  Use support_fun, Only: get_file_N
  
  implicit none

contains

  !=============================================================================
  ! Subroutine to load execution input parameters
  !=============================================================================
  Subroutine load_exec_parameters(ep,ep_infile)

    Type(exec_parameters), Intent(InOut) :: ep
    Character(Len=*)      , Intent(in)   :: ep_infile

    Integer                                :: un_spi, io_stat, pos, spos, ll, ii
    INteger                                :: no_sim_regions, n_lines
    Character(Len=1024),Allocatable,Dimension(:)  :: lines
    
    !---------------------------------------------------------------------------
    Open(newunit=un_spi, file=ep_infile, action="read", status="old", &
         iostat=io_stat)

    if (io_stat .ne. 0 ) then
       write(*,fmt_file_missing)ep_infile
       stop
    end if
    
    n_lines = get_file_N(un_spi)

    Allocate(lines(n_lines))

    Read(un_spi,'(A)',iostat=io_stat)lines

    Close(un_spi)

    !-------------------------------------------------------

    ii = 1

    Do while (ii < n_lines)

       if (lines(ii)(1:9) == "#data_dir") then !---------------------------------
          ii = ii + 1
          ep%data_dir = lines(ii)
          
       else  !-- Unknown keyword -----------------------------------------------

          
          write(*,*)"Found unknown keyword: ",trim(lines(ii))
          ii = ii + 1
          
          Do while ((lines(ii)(1:1) .NE. "#") .AND. (io_stat==0)) 

             Write(*,*)"Line image in unknown keyword:",trim(lines(ii))
             ii = ii +1 

          End Do
          
       End if

       ii = ii + 1
       
    End Do
    
    close(un_spi)
        
  end Subroutine load_exec_parameters

  !=============================================================================
  ! Subroutine to log static input parameters
  !=============================================================================
  Subroutine log_exec_parameters(ep)

    Type(exec_parameters), intent(in) :: ep

    character(len=*),parameter                    :: lca = '(A,T18,"| ",A)'
    character(len=*),parameter                    :: lci = '(A,T18,"| ",I0)'
    character(len=*),parameter                    :: lcl = '(A,T18,"| ",L)'
    character(len=*),parameter                    :: tab_h = '(17("-"),"+",42("-"))'
    
    write(*,'(60("-"))')    
    write(*,'("--",1X,A)')"execution parameters"
    write(*,tab_h)
    write(*,lca)"data_dir"     , trim(ep%data_dir)
    write(*,tab_h)
    write(*,*)
     
  End Subroutine log_exec_parameters
  
  !=============================================================================
  ! Subroutine to load static input data
  !=============================================================================
  Subroutine load_static_parameters(sp,sp_infile)

    Type(static_parameters), Intent(InOut) :: sp
    Character(Len=*)      , Intent(in)     :: sp_infile

    Integer                                :: un_spi, io_stat, pos, spos, ll,ii, jj
    INteger                                :: no_sim_regions, n_lines
    Character(Len=1024),Allocatable,Dimension(:)  :: lines

    character(len=64), Dimension(4)        :: c64_x_4

    !---------------------------------------------------------------------------
    Open(newunit=un_spi, file=sp_infile, action="read", status="old", &
         iostat=io_stat)

    if (io_stat .ne. 0 ) then
       write(*,fmt_file_missing)sp_infile
       stop
    end if
    
    n_lines = get_file_N(un_spi)

    Allocate(lines(n_lines))

    Read(un_spi,'(A)',iostat=io_stat)lines

    Close(un_spi)

    !Do ii = 1, n_lines
    !   write(*,*)trim(lines(ii))
    !End Do
    !-------------------------------------------------------

    ii = 1
    Do while (ii <= n_lines)

       if (lines(ii)(1:12) == "#sim_regions") then !----------------------------

          Read(lines(ii)(13:len_trim(lines(ii))),*)no_sim_regions

          if (allocated(sp%sim_regions)) then

             write(*,*)"sp%sim_regions already allocated! Did you specify #sim_regions twice?"
             stop

          End if
          
          Allocate(sp%sim_regions(no_sim_regions))

          pos=1
          Do jj = 1, (no_sim_regions-1) / 4 + 1

             c64_x_4 = ""
             Read(lines(ii+jj),*,iostat=io_stat)c64_x_4
             sp%sim_regions(pos:pos+3) = c64_x_4
             pos = pos + 4

          End Do
          
          ii = ii +  (no_sim_regions-1) / 4 + 1

       else if (trim(lines(ii)) == "#country") then !---------------------------
          ii = ii + 1
          sp%country = lines(ii)

       else if (trim(lines(ii)) == "#seed_date") then !-------------------------
          ii = ii + 1
          sp%seed_date = lines(ii)

       else if (trim(lines(ii)) == "#restrict") then !--------------------------
          ii = ii + 1
          Read(lines(ii),*)sp%restrict

       else if (trim(lines(ii)) == "#trans_pr") then !--------------------------
          ii = ii + 1
          sp%trans_pr = lines(ii)
          
       else if (trim(lines(ii)) == "#lhc_samples") then !--------------------------
          ii = ii + 1
          Read(lines(ii),*)sp%lhc_samples
          
       else if (trim(lines(ii)) == "#pop_data") then !--------------------------
          ii = ii + 1
          sp%pop_data = lines(ii)

       else if (trim(lines(ii)) == "#inf_cases") then !-------------------------
          ii = ii + 1
          sp%inf_cases = lines(ii)

       else if (trim(lines(ii)) == "#dead_cases") then !------------------------
          ii = ii + 1
          sp%dead_cases = lines(ii)

       else if (trim(lines(ii)) == "#connect_work") then !----------------------
          ii = ii + 1
          sp%connect_work = lines(ii)

       else if (trim(lines(ii)) == "#connect_total") then !---------------------
          ii = ii + 1
          sp%connect_total = lines(ii)

       else if (trim(lines(ii)) == "#states") then !----------------------------
          ii = ii + 1
          sp%states = lines(ii)

       else if (trim(lines(ii)) == "#counties") then !--------------------------
          ii = ii + 1
          sp%counties = lines(ii)

       else if (trim(lines(ii)) == "#R0_effects") then !---------------------
          ii = ii + 1
          sp%R0_effects = lines(ii)

       else  !-- Unknown keyword -----------------------------------------------

          write(*,*)"Found unknown keyword: ",trim(lines(ii))
          ii = ii + 1
          
          Do while ((lines(ii)(1:1) .NE. "#") .AND. (io_stat==0)) 

             Write(*,*)"Line image in unknown keyword:",trim(lines(ii))
             ii = ii +1 

          End Do

          ii = ii - 1
             
       End if

       ii = ii + 1
       
    End Do

  end Subroutine load_static_parameters

  !=============================================================================
  ! Subroutine to log static input parameters
  !=============================================================================
  Subroutine log_static_parameters(sp)

    Type(static_parameters), intent(in) :: sp

    character(len=*),parameter                    :: lca = '(A,T18,"| ",A)'
    character(len=*),parameter                    :: lci = '(A,T18,"| ",I0)'
    character(len=*),parameter                    :: lcl = '(A,T18,"| ",L)'
    character(len=*),parameter                    :: tab_h = '(17("-"),"+",42("-"))'

    integer                                       :: ii
    
    write(*,'(60("-"))')    
    write(*,'("--",1X,A)')"static parameters"
    write(*,tab_h)
    write(*,lca)"country"      , trim(sp%country      )
    write(*,lca)"seed_date"    , trim(sp%seed_date    )
    write(*,lcl)"restrict"     , sp%restrict     
    write(*,lci)"lhc_samples"  , sp%lhc_samples
    write(*,lca)"trans_pr"     , trim(sp%trans_pr     )
    write(*,lca)"pop_data"     , trim(sp%pop_data     )
    write(*,lca)"inf_cases"    , trim(sp%inf_cases    )
    write(*,lca)"dead_cases"   , trim(sp%dead_cases   )
    write(*,lca)"connect_total", trim(sp%connect_total)
    write(*,lca)"connect_work" , trim(sp%connect_work )
    write(*,lca)"states"       , trim(sp%states       )
    write(*,lca)"counties"     , trim(sp%counties     )
    write(*,lca)"R0_effects"   , trim(sp%R0_effects   )
    write(*,tab_h)
    write(*,lci)"# sim_regions", size(sp%sim_regions  )
    do ii = 1, size(sp%sim_regions)
       write(*,lca)" ",trim(sp%sim_regions(ii))
    End do
    write(*,tab_h)
    write(*,*)
     
  End Subroutine log_static_parameters

  !=============================================================================
  ! Subroutine to load input data from files
  !=============================================================================
  Subroutine loaddata(iol, pspace, ep, sp)

    Type(exec_parameters)  , intent(in) :: ep
    Type(static_parameters), intent(in) :: sp

    Integer                    :: i,j,k,it_ss
    Integer                    :: index, un_in
    Type(iols)                 :: iol
    Type(pspaces)              :: pspace

    !first version, set the array length manually
    !!
    ! read data from file to tran_pr
    Open(newunit=un_in, file=trim(ep%data_dir)//trim(sp%trans_pr),&
         access='sequential',form="formatted",iostat=k, status="old")

    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(ep%data_dir)//trim(sp%trans_pr)
       stop
    end if
    
    index = get_file_N(un_in)
    
    ! allocation for the readin variables
    Allocate(iol%transpr_age_gr(index-1))
    Allocate(iol%transpr_sex(index-1))
    Allocate(iol%transpr_surv_ill(index-1))
    Allocate(iol%transpr_icu_risk(index-1))
    Allocate(iol%transpr_surv_icu(index-1))
    ! read the first line(character)
    Read(un_in,*,iostat= k) iol%titel

    Do i = 1, index-1
       Read(un_in,*,iostat=k) iol%transpr_age_gr(i),iol%transpr_sex(i),&
            iol%transpr_surv_ill(i),&
            iol%transpr_icu_risk(i),iol%transpr_surv_icu(i)
    End Do
    Close(un_in)

    !read data from file to pop
    Open(13,file=Trim(ep%data_dir)//trim(sp%pop_data),&
         access='sequential',form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(ep%data_dir)//trim(sp%pop_data)
       stop
    end if
    
    index = get_file_N(13)
    ! allocation for the readin variables
    Allocate(iol%pop_distid(index-1))
    Allocate(iol%pop_date(index-1))
    Allocate(iol%pop_sex(index-1))
    Allocate(iol%pop_age(index-1))
    Allocate(iol%pop_total(index-1))
    ! read the first line(character)
    Read(13,*,iostat= k) iol%titel

    Do i = 1,index-1
       Read(13,*,iostat=k) iol%pop_distid(i),iol%pop_date(i),iol%pop_sex(i),&
            iol%pop_age(i),iol%pop_total(i)
    Enddo
    Close(13)

    !read data from file to seed
    Open(14,file=Trim(ep%data_dir)//trim(sp%inf_cases),access='sequential',&
         form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(ep%data_dir)//trim(sp%trans_pr)
       stop
    end if

    index = get_file_N(14)
    ! allocation for the readin variables
    Allocate(iol%seed_distid(index-1))
    Allocate(iol%seed_date(index-1))
    Allocate(iol%seed_cases(index-1))
    ! read the first line(character)
    Read(14,*,iostat= k) iol%seed_titel

    Do i = 1,index-1
       Read(14,*,iostat=k) iol%seed_distid(i),iol%seed_date(i),iol%seed_cases(i)
    Enddo
    Close(14)

    !read data from file to seeddeath

    Open(15,file=Trim(ep%data_dir)//trim(sp%dead_cases), &
         access='sequential',form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(ep%data_dir)//trim(sp%trans_pr)
       stop
    end if

    index = get_file_N(15)
    ! allocation for the readin variables
    Allocate(iol%death_distid(index-1))
    Allocate(iol%death_date(index-1))
    Allocate(iol%death_cases(index-1))
    ! read the first line(character)
    Read(15,*,iostat= k) iol%seed_titel

    Do i = 1,Size(iol%death_distid)
       Read(15,*,iostat=k) iol%death_distid(i),iol%death_date(i),iol%death_cases(i)
    Enddo
    Close(15)


    !read data from file to connect_total

    Open(16,file=Trim(ep%data_dir)//trim(sp%connect_total), &
         access='sequential',form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(ep%data_dir)//trim(sp%trans_pr)
       stop
    end if

    Read(16,*,iostat= k) iol%connect_titel,iol%connect_total_name(:)
    Do i = 1,Size(iol%connect_total_distid)
       Read(16,*,iostat=k) iol%connect_total_distid(i),iol%connect_total(:,i)
    Enddo
    Close(16)

    !read data from file to connect_total

    Open(17,file=Trim(ep%data_dir)//trim(sp%connect_work),access='sequential',&
         form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(ep%data_dir)//trim(sp%trans_pr)
       stop
    end if

    Read(17,*,iostat= k) iol%connect_titel,iol%connect_work_name(:)
    Do i = 1,Size(iol%connect_work_distid)
       Read(17,*,iostat=k) iol%connect_work_distid(i),iol%connect_work(:,i)
    Enddo


    Close(17)
    !read data from file to connect_total
    Open(18,file=Trim(ep%data_dir)//trim(sp%states), &
         access='sequential',form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(ep%data_dir)//trim(sp%trans_pr)
       stop
    end if

    index = get_file_N(18)
    ! allocation for the readin variables
    Allocate(iol%states_code(index-1))
    Allocate(iol%states_inhabitant(index-1))
    Allocate(iol%states_shortcut(index-1))
    Allocate(iol%states_name(index-1))
    Read(18,*,iostat= k) iol%state_titel(:)
    Do i = 1,index-1 
       Read(18,*,iostat=k) iol%states_code(i),iol%states_inhabitant(i),iol%states_shortcut(i),iol%states_name(i)
    Enddo
    Close(18)

    Open(19,file=Trim(ep%data_dir)//trim(sp%counties),access='sequential',&
         form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(ep%data_dir)//trim(sp%trans_pr)
       stop
    end if

    index = get_file_N(19)
    ! allocation for the readin variables
    Allocate(iol%counties_dist_id(index-1))
    Allocate(iol%counties_name(index-1))
    Allocate(iol%counties_area(index-1))
    Allocate(iol%counties_inhabitants(index-1))
    Read(19,*,iostat= k) iol%state_titel(:)
    Do i = 1,index-1
       Read(19,*,iostat=k) iol%counties_dist_id(i),iol%counties_name(i),iol%counties_area(i),iol%counties_inhabitants(i)
    Enddo
    Close(19)


!!!!----initalizing pspace ============================================
    !sam_size
    pspace%Ps_scalar_list(1)%param                 = 232000
    pspace%Ps_scalar_list(1)%var_type              = "direct"
    pspace%Ps_scalar_list(1)%name                  = "same_size"
    !R0
    pspace%Ps_scalar_list(2)%param                 = 3.5
    pspace%Ps_scalar_list(2)%var_type              = "direct"
    pspace%Ps_scalar_list(2)%name                  = "R0"
    !icu_dur
    pspace%Ps_scalar_list(3)%param                 = 14
    pspace%Ps_scalar_list(3)%var_type              = "direct"
    pspace%Ps_scalar_list(3)%name                  = "icu_dur"
    !mod_surv_ill
    pspace%Ps_scalar_list(4)%param                 = 1
    pspace%Ps_scalar_list(4)%var_type              = "direct"
    pspace%Ps_scalar_list(4)%name                  = "mod_surv_ill"
    !lcokdown_effect
    pspace%Ps_scalar_list(5)%param                 = 0.39
    pspace%Ps_scalar_list(5)%var_type              = "direct"
    pspace%Ps_scalar_list(5)%name                  = "lcokdown_effect"
    !w_int
    pspace%Ps_scalar_list(6)%param                 = 0.9
    pspace%Ps_scalar_list(6)%var_type              = "direct"
    pspace%Ps_scalar_list(6)%name                  = "w_int"
    !w_obs
    pspace%Ps_scalar_list(7)%param                 = 0.0
    pspace%Ps_scalar_list(7)%var_type              = "direct"
    pspace%Ps_scalar_list(7)%name                  = "w_obs"
    !w_obs_by_state
    pspace%Ps_scalar_list(8)%param                 = 0.0
    pspace%Ps_scalar_list(8)%var_type              = "direct"
    pspace%Ps_scalar_list(8)%name                  = "w_obs_by_state"
    !ROeffect_ps

    Open(20,file=Trim(ep%data_dir)//trim(sp%R0_effects),&
         access='sequential',form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)Trim(ep%data_dir)//trim(sp%R0_effects)
       stop
    end if
    
    index = get_file_N(20)

    Allocate(pspace%ROeffect_ps%param_char(17,index))
    Allocate(pspace%ROeffect_ps%param(16,index-1))
    Read(20,*,iostat = k)
    Do i = 1,index-1
       Read(20,*,iostat = k) pspace%ROeffect_ps%param_char(:,i)
    End Do

    Do i = 1,index-1
       Do j= 1,16
          Read(pspace%ROeffect_ps%param_char(j+1,i),"(f20.4)") pspace%ROeffect_ps%param(j,i)
       End Do
    End Do


    pspace%ROeffect_ps%var_type = "directl"
    ! convert character_param to numeric
    ! some code here, should be here

  End Subroutine loaddata


subroutine print_cosmic_head()
  
  write(*,*)
  write(*,'(A)')"!###############################################################################"
  write(*,'(A)')"!###############################################################################"
  write(*,'(A)')"!#      ___      __         _      "
  write(*,'(A)')"!#     / __\___ / _\  /\/\ (_) ___ "
  write(*,'(A)')"!#    / /  / _ \\ \  /    \| |/ __|"
  write(*,'(A)')"!#   / /__| (_) |\ \/ /\/\ \ | (__ "
  write(*,'(A)')"!#   \____/\___/\__/\/    \/_|\___|"
  write(*,'(A)')"!#"
  write(*,'(A)')"!#  COVID-19 Spatial Microsimulation  ---  For Germany  ########################"
  write(*,'(A)')"!###############################################################################"
  write(*,*)

end subroutine print_cosmic_head

end module CoSMic_IO
