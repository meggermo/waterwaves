module meggermo_kernel_test

   use :: testdrive, only: error_type, unittest_type, new_unittest, check
   use :: meggermo, only: rk
   use :: meggermo_interpolation
   use :: meggermo_kernel

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
      type(t_kernel_ref) :: g1, g2, g3, g4
      type(t_kernel_ref) :: h1, h2, h3, h4
      type(t_elem_params) :: elem_params
      real(kind=rk) :: x_e(4,2), q(2), K(2), x, dx
      integer :: j, n

      g1 = allocate_kernel('G', 1)
      g2 = allocate_kernel('G', 2)
      g3 = allocate_kernel('G', 3)
      g4 = allocate_kernel('G', 4)
      data elem_params%x_e &
         / -2.0, -1.0, 1.0, 2.0 &
         ,  0.0,  0.0, 0.0, 0.0 &
         /
      data elem_params%q / 0.5, 0.5 /
      g1%elem_params = elem_params
      g2%elem_params = elem_params
      g3%elem_params = elem_params
      g4%elem_params = elem_params
      n = 20
      dx = 2.0 / n
      do j = 0, n
         x = -1.0 + j * dx
         write(*,'(X,5E18.10)') x, g2%eval(x), g3%eval(x), g1%eval(x), g4%eval(x)
      end do

   end subroutine

   subroutine test_kernel_params(error)
      type(error_type), allocatable, intent(out) :: error
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

