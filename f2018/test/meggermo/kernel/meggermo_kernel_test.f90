module meggermo_kernel_test

   use, intrinsic :: iso_fortran_env, only: real64
   use testdrive, only: error_type, unittest_type, new_unittest, check

   use meggermo_kernel, only: ElemParams, KernelParams, kernel_params, G_ij, H_ij

   implicit none
   private

   public :: test_meggermo_kernel
   integer, parameter :: dp = selected_real_kind(15)

contains

   subroutine test_meggermo_kernel(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("test_kernel_params", test_kernel_params)]
   end subroutine

   subroutine test_kernel_params(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: t
      type(ElemParams) :: ep
      type(KernelParams) :: kp

      integer :: i

      data ep%x_e(1, :)/-1.5, 1.5/
      data ep%x_e(2, :)/-0.5, 0.5/
      data ep%x_e(3, :)/0.5, -0.5/
      data ep%x_e(4, :)/1.5, -1.5/

      data ep%q/0.1, -0.2/

      do i = 0, 50
         t = 0.02*i
         call ep%to_kernel_params(t, kp)
         call check(error, kp%Jac, sqrt(2.0_dp))
         call check(error, dot_product(kp%n, kp%n), 1.0_dp)
         write (*, '(1X,1E16.8,";",1E16.8,";",1E16.8)') t, G_ij(kp), H_IJ(kp)
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

