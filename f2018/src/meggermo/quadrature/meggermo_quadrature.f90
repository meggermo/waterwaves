module meggermo_quadrature

   use :: meggermo, only: rk

   implicit none
   private

   public T_FunParams, P_Eval, integrate_trapezoid

   type, abstract :: T_FunParams
      procedure(P_Eval), pointer, pass(fp) :: eval => null()
   end type

   abstract interface
      function P_Eval(x, fp) result(f)
         import :: T_FunParams, rk
         class(T_FunParams) :: fp
         real(rk), intent(in) :: x
         real(rk) :: f
      end function P_Eval
   end interface

contains

   subroutine integrate_trapezoid(fp, x_b, x_e, result)

      class(T_FunParams), intent(in) :: fp
      real(rk), intent(in) :: x_b, x_e
      real(rk), intent(out) ::result

      integer :: i
      integer, parameter :: steps = 100
      real(rk) :: dx

      result = (fp%eval(x_b) + fp%eval(x_e))/2.0
      dx = (x_e - x_b)/steps
      do i = 2, steps
         result = result + fp%eval(x_b + dx*(i - 1))
      end do

   end subroutine

end module
