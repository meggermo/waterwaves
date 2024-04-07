module meggermo_kernel


   use :: quadrature_module, only:integration_class_1d, quadrature_method

   use :: meggermo, only: rk
   use :: meggermo_interpolation_overhauser, only: overhauser

   implicit none
   private

   public t_elemparams, t_kernelparams, t_kernel, t_gkernel, t_hkernel, &
      kernel_params, G_ij, H_ij, kernel_function, integrate_kernels

   type t_kernelparams
      real(kind=rk) :: Jac
      !! Jacobian at t
      real(kind=rk) :: X
      !! Distance from p to q at t
      real(kind=rk) :: p(2)
      !! Point on curve at t
      real(kind=rk) :: n(2)
      !! Normal at t
   end type

   type t_elemparams
      real(kind=rk) :: x_e(4, 2)
      real(kind=rk) :: q(2)
   contains
      procedure :: to_kernel_params => kernel_params
   end type

   type, abstract, extends(integration_class_1d) :: t_kernel
      type(t_elemparams) :: ep
      integer :: i
   end type
   type, extends(t_kernel) :: t_gkernel
   end type
   type, extends(t_kernel) :: t_hkernel
   end type

contains

   subroutine integrate_kernels(tol, npoints, gk, hk, g_int, h_int)
      real(kind=rk), intent(in) :: tol
      integer, intent(in) :: npoints
      class(t_gkernel), intent(inout) :: gk
      class(t_hkernel), intent(inout) :: hk
      real(kind=rk), intent(out) :: g_int(4)
      real(kind=rk), intent(out) :: h_int(4)

      integer :: i, ierr
      real(kind=rk) :: err
      real(kind=rk), parameter :: a = 0.0
      real(kind=rk), parameter :: b = 1.0
      do i = 1, 4
         gk%i = i
         hk%i = i
         call gk%initialize(kernel_function, a, b, tol, npoints)
         call gk%integrate(g_int(i), ierr, err)
         call hk%initialize(kernel_function, a, b, tol, npoints)
         call hk%integrate(h_int(i), ierr, err)
      end do
   end subroutine

   real(kind=rk) function kernel_function(ker, t) result(I)

      class(integration_class_1d), intent(inout)  :: ker
      real(kind=rk), intent(in) :: t

      type(t_kernelparams) :: kp

      select type (ker)
       class is (t_kernel)
         call kernel_params(ker%ep, t, kp)
         select type (ker)
          class is (t_gkernel)
            I = overhauser(ker%i, t)*G_IJ(kp)
          class is (t_hkernel)
            I = overhauser(ker%i, t)*H_IJ(kp)
         end select
      end select
   end function

   subroutine kernel_params(ep, t, kp)

      class(t_elemparams), intent(in) :: ep
      real(kind=rk), intent(in) :: t
      type(t_kernelparams), intent(out) :: kp

      real(kind=rk) :: x_e(4), y_e(4)
      real(kind=rk) :: Axy(2), A, Ab, Ac, Ad
      real(kind=rk) :: Bxy(2), B, Bc, Bd
      real(kind=rk) :: Cxy(2), C, Cd
      real(kind=rk) :: Dxy(2), D

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

   real(kind=rk) function H_ij(kp)
      type(t_kernelparams), intent(in) :: kp
      H_ij = kp%Jac*dot_product(kp%p, kp%n)/kp%X
   end function

   real(kind=rk) function G_ij(kp)
      type(t_kernelparams), intent(in) :: kp
      G_ij = kp%Jac*log(kp%X)
   end function

end module
