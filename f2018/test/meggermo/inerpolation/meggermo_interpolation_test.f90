module meggermo_interpolation_test

   use meggermo
   use meggermo_params
   use testdrive, only: error_type, unittest_type, new_unittest, check

   use meggermo_interpolation, only: eval, n_weights, dn_weights

   implicit none
   private

   public :: test_meggermo_interpolation

contains

   subroutine test_meggermo_interpolation(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("test_wighted_overhauser", test_weighted_overhauser), &
                  new_unittest("test_wighted_overhauser_derivs", test_weighted_overhauser_derivs)]
   end subroutine

   subroutine test_weighted_overhauser_derivs(error)
      type(error_type), allocatable, intent(out) :: error
      !
      real(rk) :: f(5), g(2), w(4), x, k1, k2

      data f/1.0, 2.0, 4.0, 8.0, 16.0/

      ! L_i+1 = 2 L_i => k_i = 2 / (1 + 2) - 1
      k1 = -1.0/3.0
      k2 = -1.0/3.0

      call n_weights(1.0_rk, k1, k2, w)
      g(2) = dot_product(w, f(1:4))
      call n_weights(-1.0_rk, k1, k2, w)
      g(1) = dot_product(w, f(2:5))
      call check(error, abs(g(1) - g(2)) .LT. 1.0E-6)

      call dn_weights(-1.0_rk, k1, k2, w)
      g(1) = dot_product(w, f(1:4))
      call dn_weights(1.0_rk, k1, k2, w)
      g(2) = dot_product(w, f(1:4))
      call check(error, abs(g(1) - g(2)) .LT. 1.0E-6)

      call dn_weights(-1.0_rk, k1, k2, w)
      g(1) = dot_product(w, f(2:5))
      call dn_weights(1.0_rk, k1, k2, w)
      g(2) = dot_product(w, f(2:5))
      call check(error, abs(g(1) - g(2)) .LT. 1.0E-6)

   end subroutine

   subroutine test_weighted_overhauser(error)
      !
      type(error_type), allocatable, intent(out) :: error
      !
      real(rk) :: x(3)
      real(rk), dimension(4) :: w, v
      real(rk) :: k

      v = (/0.0, 1.0, 0.0, 0.0/)
      call n_weights(-one, k, k, w)
      call check_w()

      v = (/-1.0, 9.0, 9.0, -1.0/)/16.0_rk
      call n_weights(0.0_rk, k, k, w)
      call check_w()

      v = (/0.0, 0.0, 1.0, 0.0/)
      call n_weights(1.0_rk, k, k, w)
      call check_w()

      k = 1.0/3.0_rk
      call n_weights(-one, k, k, w)
      call check(error, sum(w), one)
      call n_weights(zero, k, k, w)
      call check(error, sum(w), one)
      call n_weights(one, k, k, w)
      call check(error, sum(w), one)

      k = 1.0/6.0_rk
      call n_weights(-one, k, k, w)
      call check(error, sum(w), one)
      call n_weights(zero, k, k, w)
      call check(error, sum(w), one)
      call n_weights(one, k, k, w)
      call check(error, sum(w), one)
   contains
      subroutine check_w
         integer :: i
         do i = 1, 4
            call check(error, w(i), v(i))
         end do
      end subroutine
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

