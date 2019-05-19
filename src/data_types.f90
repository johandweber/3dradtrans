!===============================================================================
!                           ****MODULE M_data_types****
!===============================================================================
! This defines kind parameters for real and integer variables

module M_data_types
  integer, parameter:: i2b = selected_int_kind(3)
  integer, parameter:: i4b = selected_int_kind(8)
  integer, parameter:: i8b = selected_int_kind(16)
  integer, parameter::  sp = selected_real_kind(6,30)
  integer, parameter::  dp = selected_real_kind(13,300)
  integer, parameter:: ldp = selected_real_kind(17,2000)
  integer, parameter::  qp = selected_real_kind(27,2500)
end module M_data_types

!===============================================================================
!                        ****END MODULE M_data_types
!===============================================================================
