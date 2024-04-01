module meggermo_interpolation_test

   use, intrinsic :: iso_fortran_env, only: real64
   use testdrive, only: error_type, unittest_type, new_unittest, check

   use meggermo_interpolation, only: &
      spline_weights, &
      overhauser_n1, &
      overhauser_n2, &
      overhauser_n3, &
      overhauser_n4, &
      eval, n_1, n_2, n_3, n_4, n_weights, dn_weights

   implicit none
   private

   public :: test_meggermo_interpolation
   integer, parameter :: dp = selected_real_kind(15)

contains

   subroutine test_meggermo_interpolation(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("test_spine_weights", test_spline_weights), &
                  new_unittest("test_overhauser_weights", test_overhauser_weights), &
                  new_unittest("test_wighted_overhauser", test_weighted_overhauser), &
                  new_unittest("test_wighted_overhauser_derivs", test_weighted_overhauser_derivs)]
   end subroutine

   subroutine test_weighted_overhauser_derivs(error)
      type(error_type), allocatable, intent(out) :: error
      !
      real(dp) :: f(5), g(2), w(4), x, k1, k2

      data f/1.0, 2.0, 4.0, 8.0, 16.0/

      ! L_i+1 = 2 L_i => k_i = 2 / (1 + 2) - 1
      k1 = -1.0/3.0
      k2 = -1.0/3.0

      call n_weights(1.0_dp, k1, k2, w)
      g(2) = dot_product(w, f(1:4))
      call n_weights(-1.0_dp, k1, k2, w)
      g(1) = dot_product(w, f(2:5))
      call check(error, abs(g(1) - g(2)) .LT. 1.0E-6)

      call dn_weights(-1.0_dp, k1, k2, w)
      g(1) = dot_product(w, f(1:4))
      call dn_weights(1.0_dp, k1, k2, w)
      g(2) = dot_product(w, f(1:4))
      write (*, *) 'g = ', g
      call check(error, abs(g(1) - g(2)) .LT. 1.0E-6)

      call dn_weights(-1.0_dp, k1, k2, w)
      g(1) = dot_product(w, f(2:5))
      call dn_weights(1.0_dp, k1, k2, w)
      g(2) = dot_product(w, f(2:5))
      write (*, *) 'g = ', g
      call check(error, abs(g(1) - g(2)) .LT. 1.0E-6)

   end subroutine

   subroutine test_weighted_overhauser(error)
      !
      type(error_type), allocatable, intent(out) :: error
      !
      real(dp) :: s(3)
      real(dp), dimension(3, 4) :: w
      real(dp) :: k

      data s/-1.0_dp, 0.0_dp, 1.0_dp/
      w(:, 1) = n_1(s, 0.0_dp)
      w(:, 2) = n_2(s, 0.0_dp, 0.0_dp)
      w(:, 3) = n_3(s, 0.0_dp, 0.0_dp)
      w(:, 4) = n_4(s, 0.0_dp)

      call check(error, w(1, 1), 0.0_dp)
      call check(error, w(2, 1), -1.0/16.0_dp)
      call check(error, w(3, 1), 0.0_dp)

      call check(error, w(1, 2), 1.0_dp)
      call check(error, w(2, 2), 9.0/16.0_dp)
      call check(error, w(3, 2), 0.0_dp)

      call check(error, w(1, 3), 0.0_dp)
      call check(error, w(2, 3), 9.0/16.0_dp)
      call check(error, w(3, 3), 1.0_dp)

      call check(error, w(1, 4), 0.0_dp)
      call check(error, w(2, 4), -1.0/16.0_dp)
      call check(error, w(3, 4), 0.0_dp)

      k = 1.0/3.0_dp
      w(:, 1) = n_1(s, k)
      w(:, 2) = n_2(s, k, k)
      w(:, 3) = n_3(s, k, k)
      w(:, 4) = n_4(s, k)

      call check(error, sum(w(1, :)), 1.0_dp)
      call check(error, sum(w(2, :)), 1.0_dp)
      call check(error, sum(w(3, :)), 1.0_dp)

      k = 1.0/6.0_dp
      w(:, 1) = n_1(s, k)
      w(:, 2) = n_2(s, k, k)
      w(:, 3) = n_3(s, k, k)
      w(:, 4) = n_4(s, k)

      call check(error, sum(w(1, :)), 1.0_dp)
      call check(error, sum(w(2, :)), 1.0_dp)
      call check(error, sum(w(3, :)), 1.0_dp)

   end subroutine

   subroutine test_overhauser_weights(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: s(3)
      real(dp), dimension(3, 4) :: w

      data s/0.0_dp, 0.5_dp, 1.0_dp/
      w(:, 1) = overhauser_n1(s)
      w(:, 2) = overhauser_n2(s)
      w(:, 3) = overhauser_n3(s)
      w(:, 4) = overhauser_n4(s)

      call check(error, w(1, 1), 0.0_dp)
      call check(error, w(2, 1), -1.0/16.0_dp)
      call check(error, w(3, 1), 0.0_dp)

      call check(error, w(1, 2), 1.0_dp)
      call check(error, w(2, 2), 9.0/16.0_dp)
      call check(error, w(3, 2), 0.0_dp)

      call check(error, w(1, 3), 0.0_dp)
      call check(error, w(2, 3), 9.0/16.0_dp)
      call check(error, w(3, 3), 1.0_dp)

      call check(error, w(1, 4), 0.0_dp)
      call check(error, w(2, 4), -1.0/16.0_dp)
      call check(error, w(3, 4), 0.0_dp)

   end subroutine

   subroutine test_spline_weights(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: s
      real(dp), dimension(4, 2) :: w

      s = 0.0_dp
      call spline_weights(s, w)
      call check(error, w(1, 1), 1.0_dp)
      call check(error, w(2, 1), 0.0_dp)
      call check(error, w(3, 1), 0.0_dp)
      call check(error, w(4, 1), 0.0_dp)
      call check(error, w(1, 2), 0.0_dp)
      call check(error, w(2, 2), 0.0_dp)
      call check(error, w(3, 2), 1.0_dp)
      call check(error, w(4, 2), 0.0_dp)

      s = 0.5_dp
      call spline_weights(s, w)
      call check(error, w(1, 1), 0.5_dp)
      call check(error, w(2, 1), 0.5_dp)
      call check(error, w(3, 1), 0.125_dp)
      call check(error, w(4, 1), -0.125_dp)
      call check(error, w(1, 2), -1.5_dp)
      call check(error, w(2, 2), 1.5_dp)
      call check(error, w(3, 2), -0.25_dp)
      call check(error, w(4, 2), -0.25_dp)

      s = 1.0_dp
      call spline_weights(s, w)
      call check(error, w(1, 1), 0.0_dp)
      call check(error, w(2, 1), 1.0_dp)
      call check(error, w(3, 1), 0.0_dp)
      call check(error, w(4, 1), 0.0_dp)
      call check(error, w(1, 2), 0.0_dp)
      call check(error, w(2, 2), 0.0_dp)
      call check(error, w(3, 2), 0.0_dp)
      call check(error, w(4, 2), 1.0_dp)

   end subroutine

end module
program tester

   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use meggermo_interpolation_test, only: test_meggermo_interpolation

   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(test_meggermo_interpolation, error_unit, stat)

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester

