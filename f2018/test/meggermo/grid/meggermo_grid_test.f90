module meggermo_grid_test

   use, intrinsic :: iso_fortran_env, only: real64, output_unit
   use testdrive, only: &
      error_type, &
      unittest_type, &
      new_unittest, &
      check

   use meggermo_grid, only: T_Grid, allocate_grid

   implicit none
   private

   public :: test_meggermo_grid
   integer, parameter :: dp = selected_real_kind(15)

contains

   subroutine test_meggermo_grid(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("test_grid_allocation", test_grid_allocation)]
   end subroutine

   subroutine test_grid_allocation(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Grid) :: grid

      grid = allocate_grid(1)
      call check(error, grid%nr_of_elements(), 1)

      grid = allocate_grid(2)
      call check(error, grid%nr_of_elements(), 2)

   end subroutine

end module
program tester

   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use meggermo_grid_test, only: test_meggermo_grid

   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(test_meggermo_grid, error_unit, stat)

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester

