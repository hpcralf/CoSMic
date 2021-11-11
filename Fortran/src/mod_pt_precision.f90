!==============================================================================
!> \file mod_pt_precision.f90
!> Type modules for the param tree input-data handling library.
!>
!> \author Ralf Schneider
!> \date 20.08.2021
!>

!==============================================================================
!> Global precision settings for the param tree input-data handling library
!>
!> The module holds global integer, real kinds and character lengths which
!> determine the precision of the param tree input-data handling library
!>
!> \author Ralf Schneider
!> \date 20.08.2021
!>
Module pt_precision

  Use ISO_FORTRAN_ENV

  Implicit None

  !> Integer Kind
  Integer, Parameter :: pt_ik = 4
  !> Real Kind
  Integer, Parameter :: pt_rk = 8
  !> MPI Integer Kind
  Integer, Parameter :: pt_mpi_ik = 4
  !> Maximum character length used in param_tree library
  Integer, Parameter :: pt_mcl = 512

  Integer, Parameter :: pt_ce = 512/4
  
End Module pt_precision
