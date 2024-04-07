module meggermo_interpolation_spline

   use meggermo, only: rk

   implicit none
   private

   public :: spline_weights

contains

   subroutine spline_weights(s, w)
      real(kind=rk), intent(in) :: s
      real(kind=rk), intent(out), dimension(4, 2) :: w
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
