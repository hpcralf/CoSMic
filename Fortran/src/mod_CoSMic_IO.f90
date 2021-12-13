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

  use mpi
  
  implicit none

contains

  !=============================================================================
  !> Subroutine to load input data from files
  Subroutine loaddata(iol)

    Type(iols), Intent(InOut)  :: iol
    
    Integer                    :: i,j,k,it_ss
    Integer                    :: index, un_in

    !=================================================================
    Character(len=:), allocatable :: filename
    !=================================================================

    ! Read data from file: Transition probabilities ------------------
    call pt_get("#trans_pr",filename)

    call read_TableData( &
         trim(filename),sep=",",head=.TRUE., &
         data=iol%trans_pr &
         )

    !read data from file to pop --------------------------------------
    call pt_get("#pop",filename)

    call read_TableData( &
         trim(filename),sep=",",head=.TRUE., &
         data=iol%pop &
         )
    
    !read data from file to seed ---------------------------
    call pt_get("#seed",filename)

    call read_TableData( &
         trim(filename),sep=",",head=.TRUE., &
         data=iol%seed &
         )

    !read data from file to seeddeath --------------------------------
    call pt_get("#seed_dea",filename)

    call read_TableData( &
         trim(filename),sep=",",head=.TRUE., &
         data=iol%death &
         )
        
    !read data from file to connect_total  ---------------------------
    call pt_get("#connect_total",filename)

    call read_TableData( &
         trim(filename),sep=",",head=.TRUE.,rownames=.TRUE., &
         data=iol%connect_total &
         )

    !read data from file to connect_work -----------------------------
    call pt_get("#connect_work",filename)

    call read_TableData( &
         trim(filename),sep=",",head=.TRUE.,rownames=.TRUE., &
         data=iol%connect_work &
         )
    
    !read data from file to states -----------------------------------
    call pt_get("#states",filename)

    call read_TableData( &
         trim(filename),sep=",",head=.TRUE., &
         data=iol%states &
         )

    !read data from file to counties ---------------------------------
    call pt_get("#counties",filename)

    call read_TableData( &
         trim(filename),sep=",",head=.TRUE., &
         data=iol%counties &
         )

    call open_and_index(trim(filename),un_in,index)

    !Read R0_effects -------------------------------------------------
    call pt_get("#R0_effects",filename)

    If (trim(filename) .NE. "LHC") then
       call read_TableData( &
            trim(filename),sep=" ",head=.TRUE., &
            rownames=.TRUE., data=iol%R0_effect &
            )
    End If

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
  !>
  !> The subroutine reads ASCII formatted table data. Per default it assumes to
  !> read the table into a 2D float array. In case characters are found in the
  !> first line that is supposed to contain data (the second line of the input
  !> file in case head=.TRUE.) It creates a data structure with columns.
  Subroutine read_TableData(filename,sep,head,rownames,data)

    Character(len=*), intent(in)           :: filename

    Type(TTableData) , intent(out)         :: data 

    Character       , intent(in), optional :: sep   
    Logical         , Intent(in), optional :: head, rownames

    Character                              :: loc_sep
    Logical                                :: loc_head, loc_rownames

    Integer                                :: io_stat, un, ii, jj, ls
    Integer                                :: no_lines, dim1, dim2, sep_pos
    character(len=16384)                   :: l_head
    character(len=:),allocatable           :: loc_str_chars

    character(len=:),Dimension(:),allocatable :: str_arr

    Character,Dimension(:), allocatable    :: tmp_col_types
    Integer,Dimension(:), allocatable      :: tmp_col_lengths

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

    !** Reduce seperator from sequence of characters signaling string data -----
    sep_pos = scan(str_chars,loc_sep)

    if (sep_pos > 1) then
       loc_str_chars = str_chars(1:sep_pos-1)//str_chars((sep_pos+1):len_trim(str_chars))
    else
       loc_str_chars = str_chars(2:len_trim(str_chars))
    End if

    !! -------------------------------------------------------------------------
    !! Open file ---------------------------------------------------------------
    Open(newunit=un,       file=Trim(filename), access='sequential', &
         form="formatted", iostat=io_stat,      status="old")

    if (io_stat .ne. 0 ) then
       write(*,fmt_file_missing)trim(filename)
       stop
    end if

    write(*,'(80("-"))')
    write(*,'(A,1X,A)')"Read report for file",trim(filename)
    
    !! -------------------------------------------------------------------------
    !! Determine dimension of table --------------------------------------------

    !! Get number of lines in file -------------------------
    no_lines = get_file_N(un)

    !! If a headline is given: Number of lines - 1 gives dim1 ----
    if (loc_head) then
       dim1 = no_lines - 1
    Else
       dim1 = no_lines
    End if

    !** Read first line in file ---
    Read(un,'(A)')l_head

    !! Warning in case maximum line length is almost used up ------
    if (len_trim(l_head) >= (len(l_head)*0.9)) then
       write(*,'("WW read_TableData:",A,I0)') &
            "Length of header line reaches limit of ",len(l_head)
    end if

    str_arr = strtok(trim(l_head),loc_sep)

    !! In case headline is included ----------------------------------
    if (loc_head) then
       !! Number of cols in head line gives dim2 ---
       dim2 = size(str_arr)

       !! Init head --------------
       data%head = str_arr

       !! Read first data line ---
       Read(un,'(A)')l_head

       !! Remove Quotes --------------------------
       Do ii = 1, dim2
          data%head(ii) = unquote(data%head(ii))
       End Do
       
    Else

       !! In case rownames are included number of columns - 1 gives dim2 ---
       if( loc_rownames ) then
          dim2 = size(str_arr)-1
       Else
          dim2 = size(str_arr)
       End if

       !! Init head ------------------------------
       Allocate(data%head(dim2),mold="123456")
       Do ii = 1, dim2
          write(data%head(ii),'("C",I5.5)')ii
       End Do

    End if

    write(*,'(A,1X,I0,1X,A)',ADVANCE="NO")"Head with size",size(data%head),"is set to      "
    write(*,'(*("|",A))')data%head
    
    !! -------------------------------------------------------------------------
    !! Determine type and in case of character the length of the columns -------
    tmp_col_types   = type_strtok(l_head,loc_sep)
    tmp_col_lengths = len_strtok(l_head,loc_sep)

    if (loc_rownames) then
       tmp_col_types   = tmp_col_types(2:dim2+1)
       tmp_col_lengths = tmp_col_lengths(2:dim2+1)
    End if
    
    data%col_lengths = tmp_col_lengths
    data%col_types   = tmp_col_types

    !! size(data%col_lengths) has to be equal dim2+1.
    if ( size(data%col_lengths) .NE. dim2 ) then
       write(*,'("WW read_TableData: ",A,/,A,L,2(/,A,I0))') &
            "Wrong number of data columns in first data line", &
            "Rownames                            :",loc_rownames, &
            "Number of columns in first data line:",size(data%col_lengths), &
            "Number of columns in header         :",dim2   
    End if

    !! Determine column type and maximum data length --------------
    !! First data line was already analysed so start at 2. --------
    Do ii = 2, dim1

       Read(un,'(A)')l_head

       tmp_col_lengths = len_strtok(l_head,loc_sep)
       tmp_col_types   = type_strtok(l_head,loc_sep)

       !! Warning in case maximum line length is almost used up ------
       if (len_trim(l_head) >= (len(l_head)*0.9)) then
          write(*,'("WW read_TableData:",A,I0)') &
               "Length of header line reaches limit of ",len(l_head)
       end if
    
       !! In case a different number of columns than before was found ---
       if ( (       loc_rownames .AND. (size(data%col_lengths) .NE. size(tmp_col_lengths)-1) ) .OR. &
            ( .NOT. loc_rownames .AND. (size(data%col_lengths) .NE. size(tmp_col_lengths)  ) ) ) then
          write(*,'("WW read_TableData:",A,1X,I0)') &
               "Different number of columns in data line", ii
          write(*,'("WW read_TableData:",A,1X,I0)') &
               "size(data%col_lengths) =",size(data%col_lengths)
          write(*,'("WW read_TableData:",A,1X,I0)') &
               "size(tmp_col_lengths) =",size(tmp_col_lengths)
       End if

       if (loc_rownames) then
          tmp_col_types   = tmp_col_types(2:dim2+1)
          tmp_col_lengths = tmp_col_lengths(2:dim2+1)
       End if

       !! Update column type and maximum length -------------------
       Do jj = 1, dim2
          if (data%col_types(jj) .NE. "c") then
             if (data%col_types(jj) .NE. "r") then
                data%col_types(jj) = tmp_col_types(jj)
             End if
          End if
          data%col_lengths(jj) = max(data%col_lengths(jj),tmp_col_lengths(jj))
       End Do

    End Do

    rewind(un)

!    if (PT_DEBUG) then
       write(*,PTF_M_A)"Column types  :",data%col_types
       write(*,PTF_M_AI0)"Column lengths:",data%col_lengths
!    End if

    !! If rownmaes are included allocate the rowname field ---
    if (loc_rownames) then
       Allocate(data%rownames(dim1),mold=l_head(1:data%col_lengths(1)))
    End if

    !! -------------------------------------------------------------------------
    !! Set data size -----------------------------------------------------------
    data%data_size(1) = dim1
    data%data_size(2) = dim2
    
    !! -------------------------------------------------------------------------
    !! Allocate columns --------------------------------------------------------
    Allocate(data%data(dim2))

    Do ii = 1, dim2

       if (data%col_types(ii) == "i") then
          allocate(data%data(ii)%cd_i(dim1))
       Else if (data%col_types(ii) == "r") then
          allocate(data%data(ii)%cd_r(dim1))
       Else
          allocate(data%data(ii)%cd_c(dim1),mold=l_head(1:data%col_lengths(ii)))
       End if
    End Do

    !! -------------------------------------------------------------------------
    !! Read data ---------------------------------------------------------------
    if (loc_head) then
       read(un,*)l_head
    End if

    if (loc_rownames) then

       Do ii = 1, dim1

          Read(un,'(A)')l_head
          str_arr = strtok(trim(l_head),loc_sep)

          data%rownames(ii) = trim(str_arr(1))
          
          Do jj = 1, dim2
             if (data%col_types(jj) == "i") then
                Read(str_arr(jj+1),*)data%data(jj)%cd_i(ii)
             Else if (data%col_types(jj) == "r") then
                Read(str_arr(jj+1),*)data%data(jj)%cd_r(ii)
             Else
                data%data(jj)%cd_c(ii) = trim(str_arr(jj+1))
             End if
          End Do

       End Do

    Else

       Do ii = 1, dim1
          Read(un,'(A)')l_head
          str_arr = strtok(trim(l_head),loc_sep)
          
          Do jj = 1, dim2
             if (data%col_types(jj) == "i") then
                Read(str_arr(jj),*)data%data(jj)%cd_i(ii)
             Else if (data%col_types(jj) == "r") then
                Read(str_arr(jj),*)data%data(jj)%cd_r(ii)
             Else
                data%data(jj)%cd_c(ii) = trim(str_arr(jj))
             End if
          End Do

       End Do

    End if

    ! In case we have character columns unquote ----------------------
    Do jj = 1, dim2

       if (data%col_types(jj) =="c") then
          Do ii = 1, dim1
             data%data(jj)%cd_c(ii) = unquote(data%data(jj)%cd_c(ii))
          End Do
       End if
       
    End Do
       
  End Subroutine read_TableData

  !! ---------------------------------------------------------------------------
  !> Function to unquote a string from " or '
  Function unquote(str) Result(out_str)

    Character(Len=*), intent(in)  :: str
    Character(Len=:), allocatable :: out_str

    integer                       :: ls
    
    out_str = adjustl(str)
    ls = len_trim(out_str)
    if ( ((out_str(1:1) == "'") .AND. (out_str(ls:ls) == "'")) .OR. &
         ((out_str(1:1) == '"') .AND. (out_str(ls:ls) == '"')) ) then
       out_str = out_str(2:len_trim(out_str)-1)
    End if

  End Function unquote
  
  !! ---------------------------------------------------------------------------
  !> Function to retreve all real columns of kind rk from the table
  Function table_to_real_array(table) Result(array)

    Type(TTableData)                          , intent(in)  :: table
    Real(kind=rk), allocatable, Dimension(:,:)              :: array

    Integer                                                 :: ii, nc

    nc = 0
    Do ii = 1, table%data_size(2)
       if (table%col_types(ii) == "r") then
          nc = nc + 1
       End if
    End Do
    
    If (nc == 0) then
       Write(*,*)"WW table_to_real_array: Found no real columns !"
    End If
    
    Allocate(array(table%data_size(1),nc))

    nc = 1
    Do ii = 1, table%data_size(2)
       if (table%col_types(ii) == "r") then
          array(:,nc) = table%data(ii)%cd_r
          nc = nc + 1
       End if
    End Do
    
  End Function table_to_real_array
  
  !! ---------------------------------------------------------------------------
  !> Function to retrieve an integer column of kind=ik from the table 
  Function get_int_table_column(table,col_name) Result(col)

    Type(TTableData)                          , intent(in)  :: table
    character(len=*)                          , intent(in)  :: col_name
    Integer(kind=ik), allocatable, Dimension(:)             :: col

    Integer                                                 :: ii

    Do ii = 1, table%data_size(2)

       if (trim(table%head(ii)) == trim(col_name)) then
          Allocate(col(table%data_size(1)))
          col = table%data(ii)%cd_i
          exit
       End if
       
    End Do
    
    if (.not. allocated(col)) then
       write(*,*)"WW get_int_column:",trim(col_name)," could not be found."
    End if
    
  End Function get_int_table_column

  !! ---------------------------------------------------------------------------
  !> Function to retrieve a real column of kind=rk from the table 
  Function get_real_table_column(table,col_name) Result(col)

    Type(TTableData)                          , intent(in)  :: table
    character(len=*)                          , intent(in)  :: col_name
    Real(kind=rk), allocatable, Dimension(:)                :: col

    Integer                                                 :: ii

    Do ii = 1, table%data_size(2)

       if (trim(table%head(ii)) == trim(col_name)) then
          Allocate(col(table%data_size(1)))
          col = table%data(ii)%cd_r
          exit
       End if
       
    End Do
    
    if (.not. allocated(col)) then
       write(*,*)"WW get_real_column:",trim(col_name)," could not be found."
    End if
    
  End Function get_real_table_column

  !! ---------------------------------------------------------------------------
  !> Function to retrieve a character column 
  Function get_char_column(table, col_name) Result(col)
    
    Type(TTableData)                           , intent(in)    :: table
    character(len=*)                           , intent(in)    :: col_name
    Character(Len=:), allocatable, Dimension(:)                :: col
    
    Integer                                                    :: ii
    
    Do ii = 1, table%data_size(2)

       if (trim(table%head(ii)) == trim(col_name)) then
          Allocate(col(table%data_size(1)),mold=table%data(ii)%cd_c(1))
          col = table%data(ii)%cd_c

          exit
       End if
       
    End Do
    
    if (.not. allocated(col)) then
       write(*,*)"WW get_char_column:",trim(col_name)," could not be found."
    End if
    
  End Function get_char_column
  
  !! ---------------------------------------------------------------------------
  !> Function to retrieve a pointer to an integer column of kind=ik
  Function get_int_column_pointer(table,col_name) Result(col)

    Type(TTableData) ,Target                  , intent(in)  :: table
    character(len=*)                          , intent(in)  :: col_name
    Integer(kind=ik), pointer, Dimension(:)                 :: col

    Integer                                                 :: ii

    col => null()
    
    Do ii = 1, table%data_size(2)

       if (trim(table%head(ii)) == trim(col_name)) then
          col => table%data(ii)%cd_i
          exit
       End if
       
    End Do

    if (.not. associated(col)) then
       write(*,*)"WW get_int_column_pointer:",trim(col_name)," could not be found."
    End if
    
  End Function get_int_column_pointer

  !! ---------------------------------------------------------------------------
  !> Function to retrieve a pointer to a real column of kind=rk
  Function get_real_column_pointer(table,col_name) Result(col)

    Type(TTableData),Target                   , intent(in)  :: table
    character(len=*)                          , intent(in)  :: col_name
    Real(kind=rk), pointer, Dimension(:)                    :: col

    Integer                                                 :: ii

    col => null()
    
    Do ii = 1, table%data_size(2)

       if (trim(table%head(ii)) == trim(col_name)) then
          col => table%data(ii)%cd_r
          exit
       End if
       
    End Do
    
    if (.not. associated(col)) then
       write(*,*)"WW get_real_column_pointer:",trim(col_name)," could not be found."
    End if
    
  End Function get_real_column_pointer

  !! ---------------------------------------------------------------------------
  !> Subroutine which broadcasts a complete table data structure      
  subroutine mpi_bcast_table(MPI_COMM, root_rank_mpi, rank_mpi, table)

    Integer(kind=mpi_ik), Intent(in) :: MPI_COMM
    Integer(kind=mpi_ik), Intent(in) :: root_rank_mpi
    Integer(kind=mpi_ik), Intent(in) :: rank_mpi
    
    Type(TTableData), Intent(InOut) :: table

    Integer(kind=mpi_ik)               :: ierr, ii
    Integer(kind=mpi_ik),Dimension(4)  :: data_dim

    Character(len=16384)                :: m_char
    
    !! -------------------------------------------------------------------------

    if ( rank_mpi == root_rank_mpi ) then
       data_dim(1:2) = table%data_size
       data_dim(3)   = len(table%head(1))
       if ( allocated(table%rownames) ) then
          data_dim(4)   = len(table%rownames(1))
       Else
          data_dim(4) = 0
       End if
    End if
    
    !** Broadcast table size, head and rownmae string length ---------
    Call mpi_bcast(data_dim, 4_mpi_ik, MPI_INTEGER4, &
         root_rank_mpi, MPI_COMM, ierr)
    
    !** Broadcast head -----------------------------------------------
    if ( rank_mpi .NE. root_rank_mpi ) then
       table%data_size = data_dim(1:2)
       Allocate(table%head(table%data_size(2)),mold=m_char(1:data_dim(3)))
    End if
    
    Call mpi_bcast(table%head, Int(table%data_size(2)*data_dim(3),mpi_ik), MPI_CHAR, &
         root_rank_mpi, MPI_COMM, ierr)

    !** Broadcast rownames -------------------------------------------
    If ( data_dim(4) > 0 ) then
       if ( rank_mpi .NE. root_rank_mpi ) then
          Allocate(table%rownames(table%data_size(1)),mold=m_char(1:data_dim(4)))
       End if

       Call mpi_bcast(table%rownames, Int(table%data_size(1)*data_dim(4),mpi_ik), MPI_CHAR, &
            root_rank_mpi, MPI_COMM, ierr)
    End If

    if ( rank_mpi .NE. root_rank_mpi ) then
       Allocate(table%col_types(table%data_size(2)))
       Allocate(table%col_lengths(table%data_size(2)))
    End if
    
    !** Broadcast col_types ------------------------------------------
    Call mpi_bcast(table%col_types, Int(table%data_size(2),mpi_ik), MPI_CHAR, &
         root_rank_mpi, MPI_COMM, ierr)

    !** Broadcast col_lengths ----------------------------------------
    Call mpi_bcast(table%col_lengths, Int(table%data_size(2),mpi_ik), MPI_INTEGER4, &
         root_rank_mpi, MPI_COMM, ierr)

    !** Allocate and broadcast data ----------------------------------
    if ( rank_mpi .NE. root_rank_mpi ) then
       Allocate(table%data(table%data_size(2)))
    End if

    Do ii = 1, table%data_size(2)

       if (table%col_types(ii) == "i") then
          
          if ( rank_mpi .NE. root_rank_mpi ) then         
             allocate(table%data(ii)%cd_i(table%data_size(1)))
          End if
          Call mpi_bcast(table%data(ii)%cd_i, Int(table%data_size(1),mpi_ik), &
               MPI_INTEGER4,root_rank_mpi, MPI_COMM, ierr)
    
       Else if (table%col_types(ii) == "r") then
    
          if ( rank_mpi .NE. root_rank_mpi ) then         
             allocate(table%data(ii)%cd_r(table%data_size(1)))
          End if
          Call mpi_bcast(table%data(ii)%cd_r, Int(table%data_size(1),mpi_ik), &
               MPI_REAL8, root_rank_mpi, MPI_COMM, ierr)
       Else

          if ( rank_mpi .NE. root_rank_mpi ) then         
             allocate(table%data(ii)%cd_c(table%data_size(1)),mold=m_char(1:table%col_lengths(ii)))
          End if

          Call mpi_bcast(table%data(ii)%cd_c, &
               Int(table%data_size(1)*table%col_lengths(ii),mpi_ik), &
               MPI_CHAR, root_rank_mpi, MPI_COMM, ierr)
       End if
    End Do

  End subroutine mpi_bcast_table
  
end module CoSMic_IO
