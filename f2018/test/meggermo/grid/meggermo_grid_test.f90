module meggermo_grid_test

   use, intrinsic :: iso_fortran_env, only: real64, output_unit
   use testdrive, only: &
      error_type, &
      unittest_type, &
      new_unittest, &
      check

   use meggermo_grid, only: T_Grid, allocate_grid, create_linear_grid, create_cubic_grid

   implicit none
   private

   public :: test_meggermo_grid
   integer, parameter :: dp = selected_real_kind(15)

contains

   subroutine test_meggermo_grid(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("test_grid_allocation", test_grid_allocation), &
                  new_unittest('test_linear_grid_creation', test_linear_grid_creation), &
                  new_unittest('test_cubic_grid_creation', test_cubic_grid_creation), &
                  new_unittest('test_glo_2_loc', test_glo_2_loc), &
                  new_unittest('test_cubic_compressed_grid', test_cubic_compressed_grid)]
   end subroutine

   subroutine test_grid_allocation(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Grid) :: grid

      grid = allocate_grid(1)
      call check(error, grid%nr_of_elements(), 1)

      grid = allocate_grid(2)
      call check(error, grid%nr_of_elements(), 2)

   end subroutine

   subroutine test_linear_grid_creation(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Grid) :: grid
      integer :: i, ne
      real(dp) :: J
      real(dp) :: x_b(2)
      real(dp) :: x_e(2)
      real(dp) :: d_x(2)

      data x_b/0.0, 0.0/
      data x_e/3.0, 3.0/
      data d_x/1.0, 1.0/

      ne = 3
      grid = allocate_grid(ne)
      call create_linear_grid(x_b, x_e, grid)

      J = sqrt(2.0_dp)
      do i = 1, ne
         call check(error, grid%x(1, i), -1.0_dp + i)
         call check(error, grid%x(2, i), -1.0_dp + i)
         call check(error, grid%J(i), J)
         call check(error, grid%K(1, i), 0.0_dp)
         call check(error, grid%K(2, i), 0.0_dp)
         call check(error, grid%n(1, i), 1.0_dp/J)
         call check(error, grid%n(2, i), -1.0_dp/J)
      end do

      call grid%compute_geom()

      do i = 1, ne
         call check(error, grid%x(1, i), -1.0_dp + i)
         call check(error, grid%x(2, i), -1.0_dp + i)
         call check(error, grid%J(i), J)
         call check(error, grid%K(1, i), 0.0_dp)
         call check(error, grid%K(2, i), 0.0_dp)
         call check(error, grid%n(1, i), 1.0_dp/J)
         call check(error, grid%n(2, i), -1.0_dp/J)
      end do

   end subroutine

   subroutine test_cubic_compressed_grid(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Grid) :: grid
      integer :: i, ne
      real(dp) :: J
      real(dp) :: x_b(2)
      real(dp) :: x_e(2)
      real(dp) :: d_x(2, 2)

      data x_b/-2.0, 0.0/
      data x_e/2.0, 0.0/
      data d_x/ &
         0.9, 1.0, &
         0.9, 1.0/

      ne = 4
      grid = allocate_grid(ne)
      call create_cubic_grid(x_b, x_e, d_x, grid)
      call grid%compute_geom()

   end subroutine

   subroutine test_cubic_grid_creation(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Grid) :: grid
      integer :: i, ne
      real(dp) :: J
      real(dp) :: x_b(2)
      real(dp) :: x_e(2)
      real(dp) :: d_x(2, 2)

      data x_b/0.0, 0.0/
      data x_e/3.0, 3.0/
      data d_x/ &
         1.0, 1.0, &
         1.0, 1.0/

      ne = 3
      J = sqrt(2.0_dp)
      grid = allocate_grid(ne)
      call create_cubic_grid(x_b, x_e, d_x, grid)
      call grid%compute_geom()

      do i = 1, ne
         d_x(:, 1) = abs(grid%x(:, i) + 1.0_dp - i)
         call check(error, d_x(1, 1) .LT. 1.0E-6, .TRUE.)
         call check(error, d_x(2, 1) .LT. 1.0E-6, .TRUE.)
         call check(error, abs(grid%J(i) - J) .LT. 1.0E-6)
         call check(error, grid%n(1, i), 1.0_dp/J)
         call check(error, grid%n(2, i), -1.0_dp/J)
      end do

   end subroutine

   subroutine test_glo_2_loc(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Grid) :: grid
      integer :: i
      integer, parameter :: ne = 3
      real(dp) :: J
      real(dp) :: x_b(2)
      real(dp) :: x_e(2)
      real(dp), allocatable :: f_l(:, :), f_g(:, :)

      data x_b/0.0, 0.0/
      data x_e/3.0, 3.0/

      grid = allocate_grid(ne)
      call create_linear_grid(x_b, x_e, grid)

      allocate (f_g(2, ne + 1))
      allocate (f_l(2, ne + 1))
      f_g(1, :) = 1.0
      f_g(2, :) = -1.0

      call grid%glo_2_loc(f_g, f_l)
      call check(error, f_l(1, 1), -sqrt(2.0_dp))
      call grid%loc_2_glo(f_l, f_g)

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

