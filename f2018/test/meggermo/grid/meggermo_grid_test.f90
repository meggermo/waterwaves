module meggermo_grid_test

   use testdrive, only: &
      error_type, &
      unittest_type, &
      new_unittest, &
      check

   use :: meggermo, only: rk
   use :: meggermo_grid, only: t_grid, t_gridtype
   use :: meggermo_grid_linear, only: t_linegridtype

   implicit none
   private

   public :: test_meggermo_grid

contains

   subroutine test_meggermo_grid(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest('test_cubic_grid_creation', test_cubic_grid_creation), &
         new_unittest('test_glo_2_loc', test_glo_2_loc), &
         new_unittest('test_cubic_compressed_grid', test_cubic_compressed_grid)]
   end subroutine

   subroutine test_cubic_compressed_grid(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_linegridtype) :: grid_type
      type(t_grid) :: grid
      integer :: i, ne
      real(rk) :: J
      real(rk) :: x_b(2)
      real(rk) :: x_e(2)
      real(rk) :: d_x(2, 2)
      real(rk), parameter :: pi = 3.1415926535897932384626433, two_pi = 2.0*pi

      data x_b/   0.0, 0.0/
      data x_e/two_pi, 0.0/

      grid_type%x_b = x_b
      grid_type%x_e = x_e

      ne = 64
      grid = grid_type%initialize(ne)
      call grid%apply_y_function(f)
      call grid%compute_geom()

   contains
      real(rk) function f(x)
         real(rk), intent(in) :: x
         f = sin(x)
      end function

   end subroutine

   subroutine test_cubic_grid_creation(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_linegridtype) :: grid_type
      type(t_grid) :: grid
      integer :: i, ne
      real(rk) :: J
      real(rk) :: x_b(2)
      real(rk) :: x_e(2)
      real(rk) :: d_x(2)

      data x_b/0.0, 0.0/
      data x_e/2.0, 2.0/

      J = sqrt(2.0_rk)
      grid_type%x_b = x_b
      grid_type%x_e = x_e

      ne = 2
      grid = grid_type%initialize(ne)

      do i = 1, ne
         d_x = abs(grid%x(:, i) + 1.0_rk - i)
         call check(error, d_x(1) .LT. 1.0E-6_rk, .TRUE.)
         call check(error, d_x(2) .LT. 1.0E-6_rk, .TRUE.)
         call check(error, abs(grid%J(i) - J) .LT. 1.0E-6)
         call check(error, grid%n(1, i), 1.0_rk/J)
         call check(error, grid%n(2, i), -1.0_rk/J)
      end do

   end subroutine

   subroutine test_glo_2_loc(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_grid) :: grid
      type(t_linegridtype) :: lgt
      integer :: i
      integer, parameter :: ne = 3
      real(rk) :: J
      real(rk) :: x_b(2)
      real(rk) :: x_e(2)
      real(rk), allocatable :: f_l(:, :), f_g(:, :)

      data x_b/0.0, 0.0/
      data x_e/3.0, 3.0/

      lgt%x_b = x_b
      lgt%x_e = x_e
      grid = lgt%initialize(ne)

      allocate (f_g(2, ne + 1))
      allocate (f_l(2, ne + 1))
      f_g(1, :) = 1.0
      f_g(2, :) = -1.0

      call grid%glo_2_loc(f_g, f_l)
      call check(error, f_l(1, 1), -sqrt(2.0_rk))
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

