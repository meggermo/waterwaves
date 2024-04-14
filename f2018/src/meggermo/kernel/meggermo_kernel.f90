module meggermo_kernel

   use :: meggermo, only: rk
   use :: meggermo_interpolation, only: f_interp, n_1, n_2, n_3, n_4

   implicit none
   private

   public :: p_g_kernel, p_h_kernel, allocate_kernel, t_kernel_ref, t_elem_params

   character, parameter :: p_g_kernel = 'G'
   character, parameter :: p_h_kernel = 'H'

   type t_kernel_params
      real(kind=rk) :: Jac
      !! Jacobian at t
      real(kind=rk) :: X
      !! Distance from p to q at t
      real(kind=rk) :: p(2)
      !! Point on curve at t
      real(kind=rk) :: q(2)
      !! Field point
      real(kind=rk) :: n(2)
      !! Normal at t
   end type

   type t_elem_params
      real(kind=rk) :: x_e(4, 2)
      real(kind=rk) :: q(2)
      real(kind=rk) :: K(2)
   end type

   type, private, abstract :: t_kernel
   contains
      procedure(f_x), deferred :: n
      procedure(f_eval), deferred :: eval
   end type

   type, private, abstract, extends(t_kernel) :: t_g_kernel
   contains
      procedure :: eval => eval_g
   end type
   type, private, extends(t_g_kernel) :: t_g_k1
   contains
      procedure :: n => f_g_n1
   end type
   type, private, extends(t_g_kernel) :: t_g_k2
   contains
      procedure :: n => f_g_n2
   end type
   type, private, extends(t_g_kernel) :: t_g_k3
   contains
      procedure :: n => f_g_n3
   end type
   type, private, extends(t_g_kernel) :: t_g_k4
   contains
      procedure :: n => f_g_n4
   end type

   type, private, abstract, extends(t_kernel) :: t_h_kernel
   contains
      procedure :: eval => eval_h
   end type
   type, private, extends(t_h_kernel) :: t_h_k1
   contains
      procedure :: n => f_h_n1
   end type
   type, private, extends(t_h_kernel) :: t_h_k2
   contains
      procedure :: n => f_h_n2
   end type
   type, private, extends(t_h_kernel) :: t_h_k3
   contains
      procedure :: n => f_h_n3
   end type
   type, private, extends(t_h_kernel) :: t_h_k4
   contains
      procedure :: n => f_h_n4
   end type

   interface
      function f_eval(kernel, x) result(f)
         import t_kernel,t_kernel_params, rk
         class(t_kernel), intent(in) :: kernel
         class(t_kernel_params), intent(in) :: x
         real(kind=rk) :: f
      end function
      function f_x(kernel, x, k1, k2) result(y)
         import t_kernel, rk
         class(t_kernel), intent(in) :: kernel
         real(kind=rk), intent(in) :: x, k1, k2
         real(kind=rk) :: y
      end function
   end interface

   type :: t_kernel_ref
      type(t_elem_params) :: elem_params
      class(t_kernel), private, allocatable :: kernel
   contains
      procedure :: eval
   end type

contains

   function f_g_n1(kernel, x, k1, k2) result(f)
      class(t_g_k1), intent(in) :: kernel
      real(kind=rk), intent(in) :: x, k1, k2
      real(kind=rk) :: f
      f = n_1(x, k1, k2)
   end function
   function f_g_n2(kernel, x, k1, k2) result(f)
      class(t_g_k2), intent(in) :: kernel
      real(kind=rk), intent(in) :: x, k1, k2
      real(kind=rk) :: f
      f = n_2(x, k1, k2)
   end function
   function f_g_n3(kernel, x, k1, k2) result(f)
      class(t_g_k3), intent(in) :: kernel
      real(kind=rk), intent(in) :: x, k1, k2
      real(kind=rk) :: f
      f = n_3(x, k1, k2)
   end function
   function f_g_n4(kernel, x, k1, k2) result(f)
      class(t_g_k4), intent(in) :: kernel
      real(kind=rk), intent(in) :: x, k1, k2
      real(kind=rk) :: f
      f = n_4(x, k1, k2)
   end function

   function f_h_n1(kernel, x, k1, k2) result(f)
      class(t_h_k1), intent(in) :: kernel
      real(kind=rk), intent(in) :: x, k1, k2
      real(kind=rk) :: f
      f = n_1(x, k1, k2)
   end function
   function f_h_n2(kernel, x, k1, k2) result(f)
      class(t_h_k2), intent(in) :: kernel
      real(kind=rk), intent(in) :: x, k1, k2
      real(kind=rk) :: f
      f = n_2(x, k1, k2)
   end function
   function f_h_n3(kernel, x, k1, k2) result(f)
      class(t_h_k3), intent(in) :: kernel
      real(kind=rk), intent(in) :: x, k1, k2
      real(kind=rk) :: f
      f = n_3(x, k1, k2)
   end function
   function f_h_n4(kernel, x, k1, k2) result(f)
      class(t_h_k4), intent(in) :: kernel
      real(kind=rk), intent(in) :: x, k1, k2
      real(kind=rk) :: f
      f = n_4(x, k1, k2)
   end function

   function allocate_kernel(kernel_type, n) result(kernel_ref)
      character, intent(in) :: kernel_type
      integer, intent(in) :: n
      type(t_kernel_ref) :: kernel_ref
      !
      select case (kernel_type)
       case (p_g_kernel)
         select case(n)
          case (1)
            kernel_ref%kernel = t_g_k1()
          case (2)
            kernel_ref%kernel = t_g_k2()
          case (3)
            kernel_ref%kernel = t_g_k3()
          case (4)
            kernel_ref%kernel = t_g_k4()
         end select
       case (p_h_kernel)
         select case(n)
          case (1)
            kernel_ref%kernel = t_h_k1()
          case (2)
            kernel_ref%kernel = t_h_k2()
          case (3)
            kernel_ref%kernel = t_h_k3()
          case (4)
            kernel_ref%kernel = t_h_k4()
         end select
      end select
   end function


   function eval(kernel_ref, x) result(f)
      class(t_kernel_ref), intent(in) :: kernel_ref
      real(kind=rk), intent(in) :: x
      real(kind=rk) :: f
      ! Local variables
      type(t_elem_params) :: elem_params
      real(kind=rk) :: k(2), kernel_val, interp_val
      type(t_kernel_params) :: kernel_params
      !
      elem_params = kernel_ref%elem_params
      k = elem_params%K
      !
      call compute_kernel_params(elem_params, x, kernel_params)

      kernel_val = kernel_ref%kernel%eval(kernel_params)
      interp_val = kernel_ref%kernel%n(x, k(1), k(2))

      f = interp_val * kernel_val

   end function


   subroutine compute_kernel_params(ep, t, kp)

      use :: meggermo_interpolation, only: n_weights, dn_weights

      class(t_elem_params), intent(in) :: ep
      real(kind=rk), intent(in) :: t
      type(t_kernel_params), intent(out) :: kp

      real(kind=rk) :: dx(2), w(4)

      kp%q = ep%q

      call n_weights(t, ep%K(1), ep%K(2), w)
      kp%p(1) = dot_product(w, ep%x_e(:, 1))
      kp%p(2) = dot_product(w, ep%x_e(:, 2))
      kp%X = norm2(kp%p - ep%q)

      call dn_weights(t, ep%K(1), ep%K(2), w)
      dx(1) = dot_product(w, ep%x_e(:, 1))
      dx(2) = dot_product(w, ep%x_e(:, 2))
      kp%Jac = norm2(dx)
      kp%n(1) = -dx(2) / kp%Jac
      kp%n(2) =  dx(1) / kp%Jac

   end subroutine


   real(kind=rk) pure function eval_g(kernel, x)
      class(t_g_kernel), intent(in) :: kernel
      class(t_kernel_params), intent(in) :: x
      eval_g = x%Jac*log(x%X)
   end function

   real(kind=rk) pure function eval_h(kernel, x)
      class(t_h_kernel), intent(in) :: kernel
      class(t_kernel_params), intent(in) :: x
      eval_h = x%Jac*dot_product(x%p - x%q, x%n)/x%X
   end function

end module
