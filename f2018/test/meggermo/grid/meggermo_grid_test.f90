module meggermo_grid_test

   use, intrinsic :: iso_fortran_env, only: real64, output_unit
   use testdrive, only: &
      error_type, &
      unittest_type, &
      new_unittest, &
      check

   use meggermo_grid, only: T_Grid, allocate_grid, create_linear_grid

   implicit none
   private

   public :: test_meggermo_grid
   integer, parameter :: dp = selected_real_kind(15)

contains

   subroutine test_meggermo_grid(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("test_grid_allocation", test_grid_allocation), &
                  new_unittest('test_grid_creation', test_grid_creation)]
   end subroutine

   subroutine test_grid_allocation(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Grid) :: grid

      grid = allocate_grid(1)
      call check(error, grid%nr_of_elements(), 1)

      grid = allocate_grid(2)
      call check(error, grid%nr_of_elements(), 2)

   end subroutine

   subroutine test_grid_creation(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Grid) :: grid
      integer :: i, ne
      real(dp) :: J
      real(dp) :: x_b(2)
      real(dp) :: x_e(2)

      data x_b/0.0, 0.0/
      data x_e/3.0, 3.0/

      ne = 3
      grid = allocate_grid(ne)
      call create_linear_grid(x_b, x_e, grid)

      J = sqrt(2.0_dp)
      do i = 0, ne + 2
         call check(error, grid%x(1, i), -1.0_dp + i)
         call check(error, grid%x(2, i), -1.0_dp + i)
         call check(error, grid%J(i), J)
         call check(error, grid%n(1, i), 1.0_dp/J)
         call check(error, grid%n(2, i), -1.0_dp/J)
      end do

      call grid%compute_geom()
      do i = 0, ne + 2
         call check(error, grid%x(1, i), -1.0_dp + i)
         call check(error, grid%x(2, i), -1.0_dp + i)
         call check(error, grid%J(i), J)
         call check(error, grid%n(1, i), 1.0_dp/J)
         call check(error, grid%n(2, i), -1.0_dp/J)
      end do

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

