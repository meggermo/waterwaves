module meggermo_quadrature_test

   use, intrinsic :: iso_fortran_env, only: rk => real64
   use testdrive, only: check, error_type, new_unittest, unittest_type

   use meggermo_integration, only: t_funparams, integrate_trapezoid

   implicit none
   private

   public :: run_unit_tests

   type, extends(t_funparams) :: T_FP
      real(rk) :: a = 1.0_rk
   end type

contains

   real(rk) function f(fp, x)
      class(t_funparams), intent(in) :: fp
      real(rk), intent(in) :: x
      select type (fp)
       class is (T_FP)
         f = fp%a*x
      end select
   end function

   subroutine run_unit_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("test_integrate_trapezoid", test_integrate_trapezoid)]
   end subroutine

   subroutine test_integrate_trapezoid(error)
      type(error_type), allocatable, intent(out) :: error

      real(rk), parameter :: x_b = 0.0_rk, x_e = 1.0_rk
      real(rk) ::result
      type(T_FP) :: params

      call params%initialize(f)
      call params%integrate_trapezoid(x_b, x_e, result)
      call check(error, result, 0.5_rk)

   end subroutine

end module

program tester

   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use meggermo_quadrature_test, only: run_unit_tests

   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(run_unit_tests, error_unit, stat)

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester

