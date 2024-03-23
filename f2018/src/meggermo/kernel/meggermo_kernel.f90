module meggermo_kernel

   use, intrinsic :: iso_fortran_env, only: real64

   implicit none
   private

   public ElemParams, KernelParams, kernel_params, G_ij, H_ij

   type KernelParams
      real(kind=real64) :: Jac
      real(kind=real64) :: X
      real(kind=real64) :: p(2)
      real(kind=real64) :: n(2)
   end type

   type ElemParams
      real(kind=real64) :: x_e(4, 2)
      real(kind=real64) :: q(2)
   contains
      procedure :: to_kernel_params => kernel_params
   end type

contains

   subroutine kernel_params(ep, t, kp)

      class(ElemParams), intent(in) :: ep
      real(kind=real64), intent(in) :: t
      type(KernelParams), intent(out) :: kp

      real(kind=real64) :: x_e(4), y_e(4)
      real(kind=real64) :: Axy(2), A, Ab, Ac, Ad
      real(kind=real64) :: Bxy(2), B, Bc, Bd
      real(kind=real64) :: Cxy(2), C, Cd
      real(kind=real64) :: Dxy(2), D

      x_e = ep%x_e(:, 1)
      y_e = ep%x_e(:, 2)

      Axy(1) = -0.5*(x_e(1) - 3.0*x_e(2) + 3.0*x_e(3) - x_e(4))
      Axy(2) = -0.5*(y_e(1) - 3.0*y_e(2) + 3.0*y_e(3) - y_e(4))
      A = dot_product(Axy, AXy)

      Bxy(1) = 0.5*(2.0*x_e(1) - 5.0*x_e(2) + 4.0*x_e(3) - x_e(4))
      Bxy(2) = 0.5*(2.0*y_e(1) - 5.0*y_e(2) + 4.0*y_e(3) - y_e(4))
      B = dot_product(Bxy, Bxy)

      Cxy(1) = -0.5*(x_e(1) - x_e(3))
      Cxy(2) = -0.5*(y_e(1) - y_e(3))
      C = dot_product(Cxy, Cxy)

      Dxy(1) = x_e(2) - ep%q(1)
      Dxy(2) = y_e(2) - ep%q(2)
      D = dot_product(Dxy, Dxy)

      Ab = 2.0*dot_product(Axy, Bxy)
      Ac = 2.0*dot_product(Axy, Cxy)
      Ad = 2.0*dot_product(Axy, Dxy)

      Bc = 2.0*dot_product(Bxy, Cxy)
      Bd = 2.0*dot_product(Bxy, Dxy)

      Cd = 2.0*dot_product(cxy, Dxy)

      kp%Jac = sqrt(t*t*(9.0*A*t*t + 6.0*Ab*t + 3.0*Ac + 4.0*B) + 2.0*Bc + C)
      kp%X = t*(t*(t*(t*(t*(A*t + Ab) + Ac + B) + Ac + Ad) + Bc + C) + Cd) + D

      kp%p(1) = t*(t*(Axy(1)*t + Bxy(1)) + Cxy(1)) + Dxy(1)
      kp%p(2) = t*(t*(Axy(2)*t + Bxy(2)) + Cxy(2)) + Dxy(2)

      kp%n(1) = (t*(3.0*Axy(1)*t + Bxy(2)) + Cxy(2))/kp%Jac
      kp%n(2) = (t*(3.0*Axy(1)*t + Bxy(1)) + Cxy(1))/kp%Jac

   end subroutine

   real(kind=real64) function H_ij(kp)
      type(KernelParams), intent(in) :: kp
      H_ij = kp%Jac*dot_product(kp%p, kp%n)/kp%X
   end function

   real(kind=real64) function G_ij(kp)
      type(KernelParams), intent(in) :: kp
      G_ij = kp%Jac*log(kp%X)
   end function

end module
