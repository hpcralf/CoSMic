!==============================================================================
!> \file mod_pt_types.f90
!> Type modules for the param tree input-data handling library.
!>
!> \author Ralf Schneider
!> \date 20.08.2021
!>

!===============================================================================
!> Global derived data-types for the param tree input-data handling library
!>
!> \author Ralf Schneider
!> \date 20.08.2021
!>
Module pt_types

  Use pt_precision

  Implicit None

  Type pt_leaf

     !> Parameter Name
     Character(Len=pt_mcl)             :: name

     !> Data Type
     !>
     !> Currently the reference is <br>
     !> I: 8 Byte Integer data <BR>
     !> R: 8 Byte Floating point data <BR>
     !> C: 1 Byte Character data <BR>
     Character                                    :: dat_ty

     !> Dimension of data
     Integer(Kind=pt_ik)                          :: dat_dim

     !> Number of data
     Integer(Kind=pt_ik),Allocatable,dimension(:) :: dat_no

     !> Data for 8 Byte integer data
     Integer(Kind=8), Dimension(:), Allocatable :: i8
     !> Data for 8 Byte integer data
     Real(Kind=8)   , Dimension(:), Allocatable :: r8
     !> Data for string data
     Character(len=pt_mcl), Dimension(:), Allocatable :: ch
     !> Data for logical data
     Logical, Dimension(:), Allocatable :: l

  End Type pt_leaf

  Type pt_branch

     !> ASCII description of what the data are
     Character(Len=pt_mcl)                      :: desc 
     !> Number of children of type pt_Branch
     Integer(Kind=pt_ik)                        :: no_branches = 0
     !> Number of children of type pt_Leaf
     Integer(Kind=pt_ik)                        :: no_leaves = 0

     !> Children of type pt_branch
     Type(pt_branch), Allocatable, Dimension(:) :: branches
     !> Children of type pt_leaf
     Type(pt_leaf  ), Allocatable, Dimension(:) :: leaves

  End Type pt_branch

End Module pt_types
