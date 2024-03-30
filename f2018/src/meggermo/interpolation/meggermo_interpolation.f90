module meggermo_interpolation

   use, intrinsic :: iso_fortran_env, only: real64

   implicit none
   private

   public :: &
      spline_weights, &
      overhauser_weights, &
      overhauser, &
      overhauser_n1, d_overhauser_n1, &
      overhauser_n2, d_overhauser_n2, &
      overhauser_n3, d_overhauser_n3, &
      overhauser_n4, d_overhauser_n4

contains


   function overhauser(i, t) result(n)
      integer, intent(in) :: i
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      select case (i)
      case (1)
         n = overhauser_n1(t)
      case (2)
         n = overhauser_n2(t)
      case (3)
         n = overhauser_n3(t)
      case (4)
         n = overhauser_n4(t)
      end select

   end function

   elemental function overhauser_n1(t) result(n)
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      n = -0.5*t*(t - 1.0)**2
   end function
   elemental function d_overhauser_n1(t) result(n)
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      n = -1.5*t*t + 2.0*t - 0.5
   end function

   elemental function overhauser_n2(t) result(n)
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      n = 0.5*(t - 1.0)*(3.0*t**2 - 2.0*t - 2.0)
   end function
   elemental function d_overhauser_n2(t) result(n)
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      n = 4.5*t*t - 5.0*t
   end function

   elemental function overhauser_n3(t) result(n)
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      n = -0.5*t*(3.0*t**2 - 4.0*t - 1.0)
   end function
   elemental function d_overhauser_n3(t) result(n)
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      n = -4.5*t*t + 4.0*t + 0.5
   end function

   elemental function overhauser_n4(t) result(n)
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      n = 0.5*t*t*(t - 1.0)
   end function
   elemental function d_overhauser_n4(t) result(n)
      real(kind=real64), intent(in) :: t
      real(kind=real64) :: n
      n = 1.5*t*t - t
   end function

   subroutine overhauser_weights(t, w)
      real(kind=real64), intent(in) :: t
      real(kind=real64), intent(out) :: w(4)
      w(1) = overhauser_n1(t)
      w(2) = overhauser_n2(t)
      w(3) = overhauser_n3(t)
      w(4) = overhauser_n4(t)
   end subroutine

   subroutine overhauser_weights_1st_deriv(t, w)
      real(kind=real64), intent(in) :: t
      real(kind=real64), intent(out) :: w(4)
      real(kind=real64), save :: N(4, 3)
      data N(1, :)/-0.5, 2.0, -1.5/
      data N(2, :)/0.0, -5.0, 4.5/
      data N(3, :)/0.5, 4.0, -4.5/
      data N(4, :)/0.0, -1.0, 1.5/
      w(1) = 1.0
      w(2) = t
      w(3) = t*t
      w = matmul(N, w(1:3))
   end subroutine

   subroutine spline_weights(s, w)
      real(kind=real64), intent(in) :: s
      real(kind=real64), intent(out), dimension(4, 2) :: w
      w(1, 1) = (s - 1.0)**2*(1.0 + 2.0*s)
      w(2, 1) = s**2*(3.0 - 2.0*s)
      w(3, 1) = (s - 1.0)**2*s
      w(4, 1) = s**2*(s - 1.0)
      w(1, 2) = 6.0*s*(s - 1.0)
      w(2, 2) = -6.0*s*(s - 1.0)
      w(3, 2) = (3.0*s - 1.0)*(s - 1.0)
      w(4, 2) = (3.0*s - 2.0)*s
   end subroutine

end module
