module meggermo_io_toml

   use tomlf, only : toml_table, get_value

   implicit none

   private

   public :: get_title

contains

   subroutine get_title(table, title)
      type(toml_table), intent(inout) :: table
      character(len=:), allocatable, intent(out) :: title
      call get_value(table, "title", title)
   end subroutine

end module
