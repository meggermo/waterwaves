module meggermo_interpolation

   use :: meggermo, only: rk

   implicit none
   private

   public :: f_interp, eval, n_weights, dn_weights, n_1, n_2, n_3, n_4

   interface
      pure function f_interp(x, k1, k2) result(I)
         import rk
         real(kind=rk), intent(in) :: x, k1, k2
         real(kind=rk) :: I
      end function
   end interface


contains

   real(rk) function eval(x, k1, k2, f)
      real(rk), intent(in) :: x, k1, k2, f(4)
      ! Local variables
      real(rk) :: w(4)
      ! Body
      call n_weights(x, k1, k2, w)
      eval = dot_product(w, f)
   end function

   subroutine n_weights(x, k1, k2, w)
      real(rk), intent(in) :: x, k1, k2
      real(rk), intent(out) :: w(4)
      ! Body
      w(1) = n_1(x, k1, k2)
      w(2) = n_2(x, k1, k2)
      w(3) = n_3(x, k1, k2)
      w(4) = n_4(x, k1, k2)
   end subroutine

   subroutine dn_weights(x, k1, k2, w)
      real(rk), intent(in) :: x, k1, k2
      real(rk), intent(out) :: w(4)
      ! Body
      w(1) = dn_1(x, k1, k2)
      w(2) = dn_2(x, k1, k2)
      w(3) = dn_3(x, k1, k2)
      w(4) = dn_4(x, k2, k2)
   end subroutine

   pure real(rk) function n_1(x, k1, k2)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: k1, k2
      ! Body
      n_1 = -(1.0 - k1)**2/(1 + k1)*(x**3 - x**2 - x + 1.0)/16.0
   end function

   pure real(rk) function dn_1(x, k1, k2)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: k1, k2
      ! Body
      dn_1 = -(1.0 - k1)**2/(1 + k1)*(3.0*x**2 - 2.0*x - 1.0)/8.0
   end function

   pure real(rk) function n_2(x, k1, k2)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: k1, k2
      ! local variables
      real(rk) :: a(4)
      ! Body
      a(4) = 3.0 - k1 + k2 + k1*k2
      a(3) = -1.0 + 3.0*k1 + k2 + k1*k2
      a(2) = -11.0 - 7.0*k1 - k2 - k1*k2
      a(1) = 9.0 + 5.0*k1 - k2 - k1*k2
      n_2 = 1.0/(1 + k1)*(a(4)*x**3 + a(3)*x**2 + a(2)*x + a(1))/16.0
   end function

   pure real(rk) function dn_2(x, k1, k2)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: k1, k2
      ! local variables
      real(rk) :: a(3)
      ! Body
      a(3) = 3.0 - k1 + k2 + k1*k2
      a(2) = -1.0 + 3.0*k1 + k2 + k1*k2
      a(1) = -11.0 - 7.0*k1 - k2 - k1*k2
      dn_2 = 1.0/(1 + k1)*(3.0*a(3)*x**2 + 2.0*a(2)*x + a(1))/8.0
   end function

   pure real(rk) function n_3(x, k1, k2)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: k1, k2
      ! local variables
      real(rk) :: a3, a2, a1, a0
      ! Body
      a3 = 3.0 - k1 + k2 + k1*k2
      a2 = 1.0 + k1 + 3.0*k2 - k1*k2
      a1 = -11.0 + k1 + 7.0*k2 - k1*k2
      a0 = -9.0 - k1 + 5.0*k2 + k1*k2
      n_3 = -1.0/(1 - k2)*(a3*x**3 + a2*x**2 + a1*x + a0)/16.0
   end function

   pure real(rk) function dn_3(x, k1, k2)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: k1, k2
      ! local variables
      real(rk) :: a3, a2, a1
      ! Body
      a3 = 3.0 - k1 + k2 + k1*k2
      a2 = 1.0 + k1 + 3.0*k2 - k1*k2
      a1 = -11.0 + k1 + 7.0*k2 - k1*k2
      dn_3 = -1.0/(1 - k2)*(3.0*a3*x**2 + 2.0*a2*x + a1)/8.0
   end function

   pure real(rk) function n_4(x, k1, k2)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: k1, k2
      ! Body
      n_4 = (1.0 + k2)**2/(1 - k2)*(x**3 + x**2 - x - 1.0)/16.0
   end function

   pure real(rk) function dn_4(x, k1, k2)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: k1, k2
      ! Body
      dn_4 = (1.0 + k2)**2/(1 - k2)*(3.0*x**2 + 2.0*x - 1.0)/8.0
   end function

end module
