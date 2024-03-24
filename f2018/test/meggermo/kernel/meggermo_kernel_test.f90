module meggermo_kernel_test

   use, intrinsic :: iso_fortran_env, only: real64, output_unit
   use testdrive, only: &
      error_type, &
      unittest_type, &
      new_unittest, &
      check

   use meggermo_kernel, only: &
      Kernel, G_Kernel, H_Kernel, KernelParams, ElemParams, &
      kernel_params, kernel_function, &
      G_ij, &
      H_ij, &
      integrate_kernels

   use meggermo_interpolation, only: &
      n1 => overhauser_n1, n2 => overhauser_n2, &
      n3 => overhauser_n3, n4 => overhauser_n4

   implicit none
   private

   public :: test_meggermo_kernel
   integer, parameter :: dp = selected_real_kind(15)

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
      real(dp) :: t, g, h, g_ans, h_ans, err, a, b, tol
      real(dp) :: g_int(4)
      real(dp) :: h_int(4)
      type(ElemParams) :: ep
      type(KernelParams) :: kp
      type(G_Kernel) :: gk
      type(H_Kernel) :: hk

      integer :: i, npoints

      data ep%x_e(1, :)/-1.0, 0.0/
      data ep%x_e(2, :)/0.0, 0.0/
      data ep%x_e(3, :)/1.0, 0.0/
      data ep%x_e(4, :)/2.0, 0.0/

      data ep%q/-0.0, 0.0/

      tol = 1.0E-6
      npoints = 6

      gk%ep = ep
      hk%ep = ep
      call integrate_kernels(tol, npoints, gk, hk, g_int, h_int)
      write (output_unit, '(1X,4E18.10)') g_int, h_int

   end subroutine

   subroutine test_kernel_params(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ierr
      real(dp) :: t, g, h, g_ans, h_ans, err, a, b, tol
      real(dp) :: g_int(4)
      real(dp) :: h_int(4)
      type(ElemParams) :: ep
      type(KernelParams) :: kp
      type(G_Kernel) :: gk
      type(H_Kernel) :: hk

      integer :: i, npoints

      data ep%x_e(1, :)/-1.5, 1.5/
      data ep%x_e(2, :)/-0.5, 0.5/
      data ep%x_e(3, :)/0.5, -0.5/
      data ep%x_e(4, :)/1.5, -1.5/

      data ep%q/-0.0, 0.0/

      tol = 1.0E-6
      npoints = 6

      gk%ep = ep
      hk%ep = ep
      call integrate_kernels(tol, npoints, gk, hk, g_int, h_int)
      write (output_unit, '(1X,4E18.10)') g_int, h_int

      do i = 0, 50
         t = 0.02*i
         call ep%to_kernel_params(t, kp)
         call check(error, kp%Jac, sqrt(2.0_dp))
         call check(error, dot_product(kp%n, kp%n), 1.0_dp)
         g = G_ij(kp)
         h = H_ij(kp)
         write (output_unit, '(1X,10(";",E18.8))') t, &
            n1(t)*g, n2(t)*g, n3(t)*g, n4(t)*g, &
            n1(t)*h, n2(t)*h, n3(t)*h, n4(t)*h
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

