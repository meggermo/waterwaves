module meggermo_kernel_test

   use testdrive, only: &
      error_type, &
      unittest_type, &
      new_unittest, &
      check

   use :: meggermo_params

   use meggermo_kernel, only: &
      t_kernel, t_gkernel, t_hkernel, t_kernelparams, t_elemparams, &
      kernel_params, kernel_function, &
      G_ij, &
      H_ij, &
      integrate_kernels

   implicit none
   private

   public :: test_meggermo_kernel

contains

   subroutine test_meggermo_kernel(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("test_kernel_params", test_kernel_params), &
         new_unittest("test_integration", test_integration)]
   end subroutine

   subroutine test_integration(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ierr
      real(rk) :: t, g, h, g_ans, h_ans, err, a, b, tol
      real(rk) :: g_int(4)
      real(rk) :: h_int(4)
      type(t_elemparams) :: ep
      type(t_kernelparams) :: kp
      type(t_gkernel) :: gk
      type(t_hkernel) :: hk

      integer :: npoints

      data ep%x_e(1, :)/-1.0, 0.0/
      data ep%x_e(2, :)/ 0.0, 0.0/
      data ep%x_e(3, :)/ 1.0, 0.0/
      data ep%x_e(4, :)/ 2.0, 0.0/

      data ep%q/0.0, 0.0/

      tol = 1.0E-12
      npoints = 6

      gk%ep = ep
      hk%ep = ep
      call integrate_kernels(tol, npoints, gk, hk, g_int, h_int)

   end subroutine

   subroutine test_kernel_params(error)
      type(error_type), allocatable, intent(out) :: error
      real(rk) :: t, g, h
      type(t_elemparams) :: ep
      type(t_kernelparams) :: kp

      integer :: i

      data ep%x_e(1, :)/-1.5, 1.5/
      data ep%x_e(2, :)/-0.5, 0.5/
      data ep%x_e(3, :)/0.5, -0.5/
      data ep%x_e(4, :)/1.5, -1.5/

      data ep%q/-0.0, 0.0/

      do i = 0, 50
         t = 0.02*i
         call ep%to_kernel_params(t, kp)
         call check(error, kp%Jac, sqrt(2.0_rk))
         call check(error, dot_product(kp%n, kp%n), 1.0_rk)
         g = G_ij(kp)
         h = H_ij(kp)
      end do

   end subroutine

end module
program tester

   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use meggermo_kernel_test, only: test_meggermo_kernel

   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(test_meggermo_kernel, error_unit, stat)

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester

