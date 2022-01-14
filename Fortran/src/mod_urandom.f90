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
!!##############################################################################
!!
!! Module containing functions to seed random numbers in shared memory setups.
!! 
!!##############################################################################
module urandom
use precision

implicit none


contains

     !! --------------------------------------------------
     !> https://masuday.github.io/fortran_tutorial/random.html
     !! standard uniform
     subroutine random_stduniform(n,u)

          integer,intent(in)            :: n  
          real,dimension(n),intent(out) :: u


          real,dimension(n)             :: r
          call random_number(r)
          u = 1 - r
     end subroutine random_stduniform

     !! uniform distribution
     function random_uniform(n,a,b) Result(x)
          integer,intent(in)           :: n
          real,intent(in)              :: a,b
          real(kind=rk),dimension(n)   :: x
          real         ,dimension(n)   :: u
          call random_stduniform(n,u)
          x = (b-a)*u + a
     !     print *,size(x)
     end function random_uniform

  !! ---------------------------------------------------------------------------
  !> Function that returns a seed array filled with random values.
  !>
  !> The values are read from `/dev/urandom`. detailed description can be
  !> found here:
  !> https://cyber.dabamos.de/programming/modernfortran/random-numbers.html
  function urandom_seed(n, stat) Result(urandom)
    
    integer, intent(in)               :: n
    integer, intent(out), optional    :: stat
    integer, Dimension(n)             :: urandom
    integer                           :: fu, rc
    
    open(access='stream', action='read', file='/dev/urandom', &
         form='unformatted', iostat=rc, newunit=fu)
    
    if (present(stat)) &
         stat = rc
    
    if (rc == 0) &
         read (fu) urandom
    
    close (fu)
    
  end function urandom_seed

End module urandom
