MODULE ModelPrecision

INTEGER, PARAMETER :: sp   = selected_REAL_kind(6, 37)     ! 32-bit
INTEGER, PARAMETER :: dp   = selected_REAL_kind(15, 307)   ! 64-bit
INTEGER, PARAMETER :: qp   = selected_REAL_kind(33, 4931)  ! 128-bit
INTEGER, PARAMETER :: prec = dp                            ! Specify the precision here

END MODULE ModelPrecision
