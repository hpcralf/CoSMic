!==============================================================================
!> \file mod_strings.f90
!> Module for string handling.
!>
!> \author Ralf Schneider
!> \date 20.08.2021
!>

!===============================================================================
!> String handling routines.
!>
!> The module provides functionalities to handle strings. Like splitting,
!> case conversion etc. .
!>
!> \author Ralf Schneider
!> \date 20.08.2021
!>
Module strings

  implicit none
  
  Character(Len=*), Parameter :: str_chars='!"#$%&'//&
       "'()*+,/:;<=>?@ABCDFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdfghijklmnopqrstuvwxyz{|}~"
  Character(Len=*), Parameter :: real_chars='.eE'
  
contains
  
  !! ---------------------------------------------------------------------------
  !> Function to determine number of lines in ASCII file
  Function get_lines_in_file(un) Result(nn)

    Integer, Intent(in) :: un
    Integer             :: nn

    Character           :: cDummy

    nn = 0

    Rewind(un)

    Do
       Read(un,*,End= 999,Err = 999) cDummy
       nn = nn + 1
    End Do

999 Rewind(un)

  End Function get_lines_in_file

  !! ---------------------------------------------------------------------------
  !> Function to split a string into tokens according to the passed seperator
  Function strtok(string,sep_char) Result(str_arr)
    
    character(len=*)     , Intent(in) :: string
    character            , Intent(in) :: sep_char
    
    integer, Dimension(:),allocatable :: pos_arr
    
    character(len=:),Dimension(:),allocatable :: str_arr

    integer                           :: ii, nn, lpos
    integer                           :: maxlen, str_len, nf, fnsc

    !---------------------------------------------------------------------------

    str_len = len(string)

    !! Determine number of seperators in string ----------------------
    nn = 0
    Do ii = 1, str_len
       if (string(ii:ii) == sep_char) nn = nn + 1 
    End Do

    !! Number of fields = number of seperators + 1 -------------------
    nf = nn + 1

    !! Save position of sep_chars in string ------------------------------
    if (str_len > 0) then

       lpos   = 0
       maxlen = -huge(maxlen)

       Allocate(pos_arr(nf+1))
       pos_arr = 0

       nn = 1

       Do ii = 1, str_len
          if (string(ii:ii) == sep_char) then
             nn = nn + 1
             pos_arr(nn) = ii
             maxlen = max(maxlen,ii-lpos-1)
             lpos = ii
          End if
       End Do

       if (lpos == 0) lpos = -1
       maxlen = max(maxlen,str_len-lpos-1)

       pos_arr(NF+1) = str_len+1

       Allocate(str_arr(nf),mold=string(1:maxlen))

       Do ii = 1, nf

          str_arr(ii) = string((pos_arr(ii)+1):pos_arr(ii+1)-1)

          !! Ajust string to the left --------------------------------
          if (str_arr(ii)(1:1) == " ") then
             fnsc = scan(str_arr(ii)," ")
             if (fnsc < len_trim(str_arr(ii))) then
                str_arr(ii) = str_arr(ii)(fnsc+1:len_trim(str_arr(ii)))
             End if
          End if
       End Do

    else

       Allocate(str_arr(0),mold="")

    End if

  End Function strtok

  !! ---------------------------------------------------------------------------
  !> Function to return the lenghts of the substrings in case a string would
  !> be split into tokens according to the passed seperator
  Function len_strtok(string,sep_char) Result(str_lengths)
    
    character(len=*)     , Intent(in) :: string
    character            , Intent(in) :: sep_char
    
    integer, Dimension(:),allocatable :: pos_arr

    integer, Dimension(:),allocatable :: str_lengths

    integer                           :: ii, nn, lpos
    integer                           :: str_len, nf

    !---------------------------------------------------------------------------

    str_len = len_trim(string)

    !! Determine number of seperators in string ----------------------
    nn = 0
    Do ii = 1, str_len
       if (string(ii:ii) == sep_char) nn = nn + 1 
    End Do

    !! Number of fields = number of seperators + 1 -------------------
    nf = nn + 1

    !! Save position of sep_chars in string ------------------------------
    if (str_len > 0) then

       lpos   = 0

       Allocate(pos_arr(nf+1))
       pos_arr(1) = 0

       nn = 1

       Do ii = 1, str_len
          if (string(ii:ii) == sep_char) then
             nn = nn + 1
             pos_arr(nn) = ii
             lpos = ii
          End if
       End Do

       pos_arr(NF+1) = len_trim(string)+1

       Allocate(str_lengths(nf))

       Do ii = 1, nf
          str_lengths(ii) = pos_arr(ii+1)-pos_arr(ii)-1
       End Do

    else

       Allocate(str_lengths(0))

    End if

  End Function len_strtok

  !! ---------------------------------------------------------------------------
  !> Function to return the types of the substrings in case a string would
  !> be split into tokens according to the passed seperator.
  Function type_strtok(in_string,sep_char) Result(str_types)
    
    character(len=*)     , Intent(in)  :: in_string
    character            , Intent(in)  :: sep_char
    
    integer, Dimension(:),allocatable  :: pos_arr

    Character,Dimension(:),allocatable :: str_types

    integer                            :: ii, nn, lpos
    integer                            :: str_len, nf, sep_pos
    integer                            :: io_stat
    
    Real                               :: tmp_real
    integer                            :: tmp_int
    
    character(len=:),allocatable       :: loc_str_chars

    character(len=:),allocatable       :: string
    
    !---------------------------------------------------------------------------

    string = adjustl(in_string)
    
    str_len = len_trim(string)

    !! Determine number of seperators in string ----------------------
    nn = 0
    Do ii = 1, str_len
       if (string(ii:ii) == sep_char) nn = nn + 1 
    End Do

    !! Number of fields = number of seperators + 1 -------------------
    nf = nn + 1

    !! Save position of sep_chars in string ------------------------------
    if (str_len > 0) then

       lpos   = 0

       Allocate(pos_arr(nf+1))
       pos_arr(1) = 0

       nn = 1

       Do ii = 1, str_len
          if (string(ii:ii) == sep_char) then
             nn = nn + 1
             pos_arr(nn) = ii
             lpos = ii
          End if
       End Do

       pos_arr(NF+1) = len_trim(string)+1

       Allocate(str_types(nf))
       
       !** Reduce seperator from sequence of characters signaling string data -----
       sep_pos = scan(str_chars,sep_char)

       if (sep_pos > 1) then
          loc_str_chars = str_chars(1:sep_pos-1)//str_chars((sep_pos+1):len_trim(str_chars))
       else
          loc_str_chars = str_chars(2:len_trim(str_chars))
       End if
       
       Do ii = 1, nf

          if( scan(string((pos_arr(ii)+1):pos_arr(ii+1)-1), loc_str_chars) > 0 ) then
             !**  Check whether we have chars signaling characters -------------
             str_types(ii) = "c"
          Else
             !** If not check whether we have chars signaling real -------------
             if (scan(string((pos_arr(ii)+1):pos_arr(ii+1)-1),real_chars) > 0 ) then

                !** Check whether real can be read -----------------------------
                Read(string((pos_arr(ii)+1):pos_arr(ii+1)-1),*,iostat=io_stat)tmp_real
                if (io_stat .NE. 0) then
                   str_types(ii) = "c"
                else
                   str_types(ii) = "r"
                end if
             Else
                Read(string((pos_arr(ii)+1):pos_arr(ii+1)-1),*,iostat=io_stat)tmp_int
                if (io_stat .NE. 0) then
                   str_types(ii) = "c"
                else
                   str_types(ii) = "i"
                end if
             End if
          End if
       End Do

    else

       Allocate(str_types(0))

    End if

  End Function type_strtok

  !! ---------------------------------------------------------------------------
  !> Function to convert characters to lower case
  Function ToLowerCase(chr) Result(lc_chr)

    character, intent(in) :: chr
    character             :: lc_chr

    if ( (65 <= ichar(chr)) .AND. (ichar(chr) <= 90) ) then
       lc_chr = achar(ichar(chr)+32)
    Else
       lc_chr = chr
    End if

  End Function ToLowerCase

  !! ---------------------------------------------------------------------------
  !> Function to convert characters to upper case
  Function ToUpperCase(chr) Result(uc_chr)

    character, intent(in) :: chr
    character             :: uc_chr

    if ( (97 <= ichar(chr)) .AND. (ichar(chr) <= 122) ) then
       uc_chr = achar(ichar(chr)-32)
    Else
       uc_chr = chr
    End if

  End Function ToUpperCase
  
End Module strings
