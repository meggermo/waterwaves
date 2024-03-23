module meggermo_io_toml_test

   use tomlf
   use testdrive, only : error_type, unittest_type, new_unittest, check
   use meggermo_io_toml, only : get_title
   implicit none

   private

   public :: test_toml_parsing

contains

   subroutine test_toml_parsing(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("test_title", test_title)]
   end subroutine

   subroutine test_title(error)
      type(error_type), allocatable, intent(out) :: error
      type(toml_table), allocatable :: table
      character(len=:), allocatable :: title
      character(len=:), allocatable :: input

      input = &
      & '# This is a TOML document.' // new_line("a") // &
      & 'title = "TOML Example"' // new_line("a")

      call toml_loads(table, input)

      if (allocated(table)) then
         call get_title(table, title)
         call check(error, title, "TOML Example")
      else
         stop "Table not allocated"
      end if

   end subroutine

end module

program tester

   use, intrinsic :: iso_fortran_env, only : error_unit
   use testdrive, only : run_testsuite
   use meggermo_io_toml_test, only : test_toml_parsing

   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(test_toml_parsing, error_unit, stat)

   if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
