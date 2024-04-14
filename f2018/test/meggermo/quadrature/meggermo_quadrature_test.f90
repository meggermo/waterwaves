module meggermo_quadrature_test

   use :: testdrive, only: check, error_type, new_unittest, unittest_type

   use :: meggermo, only: rk
   use :: meggermo_integration, only: t_quad_params, t_quad_result, do_quad => quad
   use :: meggermo_kernel, only: t_elem_params, t_kernel_ref, allocate_kernel
   use :: meggermo_quadrature, only: t_quad, initialize

   implicit none
   private

   public :: run_unit_tests

contains

   subroutine run_unit_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("test_integrate_trapezoid", test_integrate_trapezoid)]
   end subroutine

   subroutine test_integrate_trapezoid(error)
      type(error_type), allocatable, intent(out) :: error

      type(t_quad_params) :: quad_params
      type(t_elem_params) :: elem_params
      real(kind=rk) :: result(4,2)
      type(t_quad) :: kernels

      data elem_params%x_e &
         / -3.0, -1.0, 1.0, 3.0 &
         ,  0.0,  0.0, 0.0, 0.0 &
         /
      data elem_params%q / 0.0, 1.0 /

      kernels = initialize(quad_params)
      call kernels%integrate(elem_params, result)

      write(*, '(4E18.10)') result

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

