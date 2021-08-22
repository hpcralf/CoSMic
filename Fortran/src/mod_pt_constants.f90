!==============================================================================
!> \file mod_pt_constants.f90
!> Type modules for the param tree input-data handling library.
!>
!> \author Ralf Schneider
!> \date 20.08.2021
!>

!===============================================================================
!> Global constants for the param tree input-data handling library
!>
!> \author Ralf Schneider
!> \date 20.08.2021
!>
Module pt_constants

  Implicit None

  !> Keyword line field seperator ----------------------------------------------
  character, parameter :: fsep = ","

  ! Character constants for nice output ----------------------------------------
  Character(Len=*), Parameter :: PTF_E_A    = "('EE ',A)"
  Character(Len=*), Parameter :: PTF_E_AI0  = "('EE ',*(A,1X,I0))"
  Character(Len=*), Parameter :: PTF_E_STOP = &
       "('EE PROGRAM STOPPED ..... ',/,'<',78('='),'>')"

  Character(Len=*), Parameter :: PTF_W_A    = "('WW ',A)"
  Character(Len=*), Parameter :: PTF_W_AI0  = "('WW ',*(A,1X,I0))"

  Character(Len=*), Parameter :: PTF_M_A    = "('MM ',A)"
  Character(Len=*), Parameter :: PTF_M_AI0  = "('MM ',A,1X,I0)"

  Character(Len=*), Parameter :: PTF_TIME   = "('MM ',A,1X,F0.6,' sec')"

  Character(Len=*), Parameter :: PTF_SEP    = "('<',78('='),'>')"

  Character(len=*), Parameter :: PTF_file_missing = &
       '("Operation on file:",/,A,/, "failed.",/, "Program halted !!!")'

End Module pt_constants
