module scalar_type

  type :: scalar
#if defined USE_COMPLEX
     complex(8) :: c
#else
     real(8) :: d
#endif
    
  end type scalar
end module scalar_type
  

program main
  use scalar_type
  type(scalar) :: a = 1.0d0
  print *, a
  
end program main
