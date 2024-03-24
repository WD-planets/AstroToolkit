Module aovconst
!
! AovConst.f90 : global definitions for Aov*.f90
      Implicit None
!
      Integer, Parameter :: SP = Kind (1e0)! type of measurements
      Integer, Parameter :: CP = 8 ! type of measurements
      Integer, Parameter :: TIME = Kind (1d0)! type of time & frequency variables
      Real (SP), Parameter :: PI2 = &
         6.283185307179586476925286766559005768394_SP
      Integer, Parameter :: CTMIN = 5 !minimum bin occupancy
      Real (SP), Parameter :: MINVAR = epsilon (1.0_SP)!minimum variance difference
!
End Module aovconst
