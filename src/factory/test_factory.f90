program main

  use factory_pattern

  implicit none

  type(CFactory) :: factory
  class(Connection), pointer :: db_connect => null()

  call factory%init("Oracle")
  db_connect => factory%create_connection()   !! Create Oracle DB
  call db_connect%description()

  !! The same factory can be used to create different connections
  call factory%init("MySQL")                  !! Create MySQL DB

  !! 'connect' is a 'class' pointer. So can be used for either Oracle or MySQL
  db_connect => factory%create_connection()
  call db_connect%description()

  call factory%final()        ! Destroy the object

end program main
