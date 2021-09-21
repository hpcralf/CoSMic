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

  Use param_tree
  Use strings
  
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
  Subroutine loaddata(iol, ep, sp)

    Type(exec_parameters)  , intent(in) :: ep
    Type(static_parameters), intent(in) :: sp

    Integer                    :: i,j,k,it_ss
    Integer                    :: index, un_in
    Type(iols)                 :: iol

    !=================================================================
    Character(len=:), allocatable :: data_dir
    Character(len=:), allocatable :: filename
    !=================================================================

    call pt_get("#data_dir",data_dir)
    
    ! Read data from file: Transition probabilities ------------------
    call pt_get("#trans_pr",filename)

    call open_and_index(Trim(data_dir)//trim(filename),un_in,index)

    ! allocation for the readin variables ------------------
    Allocate(iol%transpr_age_gr(index-1))
    Allocate(iol%transpr_sex(index-1))
    Allocate(iol%transpr_surv_ill(index-1))
    Allocate(iol%transpr_icu_risk(index-1))
    Allocate(iol%transpr_surv_icu(index-1))

    ! read the first line(character) -----------------------
    Read(un_in,*,iostat= k) iol%titel

    Do i = 1, index-1
       Read(un_in,*,iostat=k) iol%transpr_age_gr(i),iol%transpr_sex(i),&
            iol%transpr_surv_ill(i),&
            iol%transpr_icu_risk(i),iol%transpr_surv_icu(i)
    End Do
    Close(un_in)

    !read data from file to pop --------------------------------------
    call pt_get("#pop_data",filename)

    call open_and_index(Trim(data_dir)//trim(filename),un_in,index)

    ! allocation for the readin variables ------------------
    Allocate(iol%pop_distid(index-1))
    Allocate(iol%pop_date(index-1))
    Allocate(iol%pop_sex(index-1))
    Allocate(iol%pop_age(index-1))
    Allocate(iol%pop_total(index-1))

    ! read the first line(character) -----------------------
    Read(un_in,*,iostat= k) iol%titel

    Do i = 1,index-1
       Read(un_in,*,iostat=k) iol%pop_distid(i),iol%pop_date(i),iol%pop_sex(i),&
            iol%pop_age(i),iol%pop_total(i)
    Enddo
    Close(un_in)

    !read data from file to seed ---------------------------
    call pt_get("#inf_cases",filename)

    call open_and_index(Trim(data_dir)//trim(filename),un_in,index)

    ! allocation for the readin variables ------------------
    Allocate(iol%seed_distid(index-1))
    Allocate(iol%seed_date(index-1))
    Allocate(iol%seed_cases(index-1))
    ! read the first line(character)
    Read(un_in,*,iostat= k) iol%seed_titel

    Do i = 1,index-1
       Read(un_in,*,iostat=k) iol%seed_distid(i),iol%seed_date(i),iol%seed_cases(i)
    Enddo
    Close(un_in)

    !read data from file to seeddeath --------------------------------
    call pt_get("#dead_cases",filename)

    call open_and_index(Trim(data_dir)//trim(filename),un_in,index)

    ! allocation for the readin variables ------------------
    Allocate(iol%death_distid(index-1))
    Allocate(iol%death_date(index-1))
    Allocate(iol%death_cases(index-1))
    
    ! read the first line(character) -----------------------
    Read(un_in,*,iostat= k) iol%seed_titel

    Do i = 1,Size(iol%death_distid)
       Read(un_in,*,iostat=k) iol%death_distid(i),iol%death_date(i),iol%death_cases(i)
    Enddo
    Close(un_in)

    !read data from file to connect_total  ---------------------------
    call pt_get("#connect_total",filename)

    call open_and_index(Trim(data_dir)//trim(filename),un_in,index)

    Read(un_in,*,iostat= k) iol%connect_titel,iol%connect_total_name(:)
    Do i = 1,Size(iol%connect_total_distid)
       Read(un_in,*,iostat=k) iol%connect_total_distid(i),iol%connect_total(:,i)
    Enddo
    Close(un_in)

    !read data from file to connect_work -----------------------------
    call pt_get("#connect_work",filename)

    call open_and_index(Trim(data_dir)//trim(filename),un_in,index)

    Read(un_in,*,iostat= k) iol%connect_titel,iol%connect_work_name(:)
    Do i = 1,Size(iol%connect_work_distid)
       Read(un_in,*,iostat=k) iol%connect_work_distid(i),iol%connect_work(:,i)
    Enddo

    Close(un_in)

    !read data from file to states -----------------------------------
    call pt_get("#states",filename)

    call open_and_index(Trim(data_dir)//trim(filename),un_in,index)

    ! allocation for the readin variables ------------------
    Allocate(iol%states_code(index-1))
    Allocate(iol%states_inhabitant(index-1))
    Allocate(iol%states_shortcut(index-1))
    Allocate(iol%states_name(index-1))
    
    Read(un_in,*,iostat= k) iol%state_titel(:)
    Do i = 1,index-1 
       Read(un_in,*,iostat=k) iol%states_code(i),iol%states_inhabitant(i),iol%states_shortcut(i),iol%states_name(i)
    Enddo
    Close(un_in)

    !read data from file to counties ---------------------------------
    call pt_get("#counties",filename)

    call open_and_index(Trim(data_dir)//trim(filename),un_in,index)

    ! allocation for the readin variables ------------------
    Allocate(iol%counties_dist_id(index-1))
    Allocate(iol%counties_name(index-1))
    Allocate(iol%counties_area(index-1))
    Allocate(iol%counties_inhabitants(index-1))
    
    Read(un_in,*,iostat= k) iol%state_titel(:)
    Do i = 1,index-1
       Read(un_in,*,iostat=k) iol%counties_dist_id(i),iol%counties_name(i),&
            iol%counties_area(i),iol%counties_inhabitants(i)
    Enddo
    Close(un_in)

    !Read R0_effects -------------------------------------------------
    call pt_get("#R0_effects",filename)
    
    call read_TableData( &
         trim(data_dir)//trim(filename),sep=" ",head=.TRUE., &
         rownames=.TRUE., data=iol%R0_effect &
         )

  End Subroutine loaddata
  
  !! ===========================================================================
  !> Print the CoSMic application logotted file and return its line count 
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
  
  !! ===========================================================================
  !> Subroutine to open an ASCII formatted file and return its line count 
  subroutine open_and_index(path,un,index)

    Integer         , intent(inout) :: un
    character(len=*), intent(in)    :: path
    integer, intent(out)            :: index

    integer                         :: k

    Open(newunit=un, file=Trim(path), &
         access='sequential',form="formatted",iostat=k, status="old")
    if (k .ne. 0 ) then
       write(*,fmt_file_missing)trim(path)
       stop
    end if

    index = get_file_N(un)

  end subroutine open_and_index

  !! ===========================================================================
  !> Subroutine to read ASCII formatted table data
  Subroutine read_TableData(filename,sep,head,rownames,data)

    Character(len=*), intent(in)           :: filename
        
    Type(TableData) , intent(out)          :: data 

    Character       , intent(in), optional :: sep   
    Logical         , Intent(in), optional :: head, rownames
    
    Character                              :: loc_sep
    Logical                                :: loc_head, loc_rownames

    Integer                                :: io_stat, un, ii
    Integer                                :: no_lines, dim1, dim2
    character(len=2048)                    :: l_head
    character(len=:),Dimension(:),allocatable :: str_arr

    !---------------------------------------------------------------------------
    if (present(sep)) then
       loc_sep = sep
    else
       loc_sep = ","
    End if

    if (present(head)) then
       loc_head = head
    else
       loc_head = .TRUE.
    End if

    if (present(rownames)) then
       loc_rownames = rownames
    else
       loc_rownames = .FALSE.
    End if

    !! Open file -----------------------------------------------------
    Open(newunit=un, file=Trim(filename), access='sequential', &
         form="formatted",iostat=io_stat, status="old")
    
    if (io_stat .ne. 0 ) then
       write(*,fmt_file_missing)trim(filename)
       stop
    end if

    !! Get number of lines in file -------------------------
    no_lines = get_file_N(un)

    if (loc_head) then
       !! If a headline is given: Number of lines - 1 gives dim1 ---------------
       
       dim1 = no_lines - 1

       Read(un,'(A)')l_head

       !! Warning in case maximum line length is almost used up ------
       if (len_trim(l_head) >= (len(l_head)*0.9)) then
          write(*,'("WW read_TableData:",A,I0)') &
               "Length of header line reaches limit of ",len(l_head)
       end if
    
       data%head = strtok(trim(l_head),loc_sep)
    
       dim2 = size(data%head)

    Else
       !! If we have no headline, determine dim2 from Number of seperators in --
       !! first line. ----------------------------------------------------------
       
       dim1 = no_lines

       Read(un,'(A)')l_head

       !! Warning in case maximum line length is almost used up ------
       if (len_trim(l_head) >= (len(l_head)*0.9)) then
          write(*,'("WW read_TableData:",A,I0)') &
               "Length of header line reaches limit of ",len(l_head)
       end if
    
       str_arr = strtok(trim(l_head),loc_sep)
    
       dim2 = size(str_arr)

       Rewind(un)
       
    End if

    Allocate(data%data(dim1,dim2))
    
    if (loc_rownames) then
       Allocate(data%rownames(dim1),mold=l_head(1:16))
       Do ii = 1, dim1
          Read(un,*)data%rownames(ii),data%data(ii,:)
       End Do
    Else
       Read(un,*)data%data
    end if

  End Subroutine read_TableData
  
end module CoSMic_IO
