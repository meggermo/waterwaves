module meggermo_integration

   use :: meggermo, only: rk

   implicit none
   private

   public T_FunParams, integrate_trapezoid

   type, abstract, private :: T_Root
      private
   end type

   type, extends(T_Root), public :: T_FunParams
      private
      procedure(eval), pointer :: f => null()
      real(rk) :: x_b = 0.0
      real(rk) :: x_e = 1.0
   contains
      private
      procedure, public :: initialize => fp_initialize
      procedure, public :: integrate_trapezoid => integrate_trapezoid
   end type

   interface
      real(rk) function eval(fp, x)
         import T_FunParams, rk
         implicit none
         class(T_FunParams), intent(in) :: fp
         real(rk), intent(in) :: x
      end function
   end interface

contains

   subroutine fp_initialize(fp, f)
      class(T_FunParams), intent(inout) :: fp
      procedure(eval) :: f
      fp%f => f
   end subroutine

   subroutine integrate_trapezoid(fp, x_b, x_e, result)

      class(T_FunParams), intent(inout) :: fp
      real(rk), intent(in) :: x_b, x_e
      real(rk), intent(out) ::result

      integer :: i
      integer, parameter :: steps = 10
      real(rk) :: dx

      result = 0.5*(fp%f(x_b) + fp%f(x_e))
      dx = (x_e - x_b)/steps
      do i = 2, steps
         result = result + fp%f(x_b + dx*(i - 1))
      end do
      result = dx*result

   end subroutine

end module
