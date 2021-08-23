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
  Function strtok(string,char) Result(str_arr)
    
    character(len=*)     , Intent(in) :: string
    character            , Intent(in) :: char
    
    integer, Dimension(:),allocatable :: pos_arr
    
    character(len=:),Dimension(:),allocatable :: str_arr

    integer                           :: ii, nn, lpos
    integer                           :: maxlen, str_len, nf, fnsc

    !---------------------------------------------------------------------------

    str_len = len(string)

    !! Determine number of seperators in string ----------------------
    nn = 0
    Do ii = 1, str_len
       if (string(ii:ii) == char) nn = nn + 1 
    End Do

    !! Number of fields = number of seperators + 1 -------------------
    nf = nn + 1

    !! Save position of chars in string ------------------------------
    if (str_len > 0) then

       lpos   = 0
       maxlen = -huge(maxlen)

       Allocate(pos_arr(nf+1))
       pos_arr(1) = 0

       nn = 1

       Do ii = 1, str_len
          if (string(ii:ii) == char) then
             nn = nn + 1
             pos_arr(nn) = ii
             maxlen = max(maxlen,ii-lpos-1)
             lpos = ii
          End if
       End Do

       maxlen = max(maxlen,str_len-lpos-1)

       pos_arr(NF+1) = str_len+1

       Allocate(str_arr(nf),mold=string(1:maxlen))

       Do ii = 1, nf

          str_arr(ii) = string((pos_arr(ii)+1):pos_arr(ii+1)-1)

          !! Ajust string to the left --------------------------------
          fnsc = scan(str_arr(ii)," ")
          if (fnsc < len_trim(str_arr(ii))) then
             str_arr(ii) = str_arr(ii)(fnsc+1:len_trim(str_arr(ii)))
          End if

       End Do

    else

       Allocate(str_arr(0),mold="")

    End if

  End Function strtok

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
