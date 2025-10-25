module factory_pattern
  implicit none
  
  type CFactory
     private
     character(len=20) :: factory_type               !! Descriptive name for database
     class(Connection), pointer :: connection_type   !! Which type of database ?
   contains                                            !! Note 'class' not 'type' !
     procedure :: init                               !! Constructor
     procedure :: create_connection                  !! Connect to database
     procedure :: final                              !! Destructor
  end type CFactory

  type, abstract :: Connection
   contains
     procedure(generic_desc), deferred, pass(self) :: description
  end type Connection

  abstract interface
     subroutine generic_desc(self)
       import :: Connection
       class(Connection), intent(in) :: self
     end subroutine generic_desc
  end interface

  ! An Oracle connection
  type, extends(Connection) :: OracleConnection
   contains
     procedure, pass(self) :: description => oracle_desc
  end type OracleConnection
  
  ! A MySQL connection
  type, extends(Connection) :: MySQLConnection
   contains
     procedure, pass(self) :: description => mysql_desc
  end type MySQLConnection

contains
  
  subroutine oracle_desc(self)
    class(OracleConnection), intent(in) :: self
    write(*,'(A)') "You are now connected with Oracle"
  end subroutine oracle_desc

  subroutine mysql_desc(self)
    class(MySQLConnection), intent(in) :: self
    write(*,'(A)')  "You are now connected with MySQL"
  end subroutine mysql_desc



  subroutine init(self, string)
    class(CFactory), intent(inout) :: self
    character(len=*), intent(in) :: string
    self%factory_type = trim(string)
    self%connection_type => null()            !! pointer is nullified
  end subroutine init

  subroutine final(self)
    class(CFactory), intent(inout) :: self
    deallocate(self%connection_type)          !! Free the memory
    nullify(self%connection_type)
  end subroutine final




  function create_connection(self)  result(ptr)
    class(CFactory) :: self
    class(Connection), pointer :: ptr

    if(self%factory_type == "Oracle") then
       if(associated(self%connection_type))   deallocate(self%connection_type)
       allocate(OracleConnection :: self%connection_type)
       ptr => self%connection_type
    elseif(self%factory_type == "MySQL") then
       if(associated(self%connection_type))   deallocate(self%connection_type)
       allocate(MySQLConnection :: self%connection_type)
       ptr => self%connection_type
    end if
  end function create_connection

end module factory_pattern
