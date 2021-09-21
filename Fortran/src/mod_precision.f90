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
! Module with precision constants
!
!###############################################################################
Module precision

  Implicit None

  Integer, Parameter :: mcl = 128
  Integer, Parameter :: rk  = 8
  Integer, Parameter :: ik  = 4
  
End Module precision
