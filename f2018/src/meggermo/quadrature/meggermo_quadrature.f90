module meggermo_quadrature

   use :: meggermo, only: rk
   use :: meggermo_integration, only: t_quad_params, t_quad_result, do_quad => quad
   use :: meggermo_kernel, only: t_kernel_ref, t_elem_params, allocate_kernel

   implicit none
   private

   public :: t_quad, initialize

   type :: t_quad
      type(t_quad_params) :: quad_params
      type(t_kernel_ref), allocatable, dimension(:, :) :: kernels
   contains
      procedure :: integrate
   end type

contains

   function initialize(quad_params)
      type(t_quad_params), intent(in) :: quad_params
      type(t_quad) :: quad, initialize
      allocate(quad%kernels(4,2))
      quad%kernels(1,1) = allocate_kernel('G', 1)
      quad%kernels(2,1) = allocate_kernel('G', 1)
      quad%kernels(3,1) = allocate_kernel('G', 1)
      quad%kernels(4,1) = allocate_kernel('G', 1)
      quad%kernels(1,2) = allocate_kernel('H', 1)
      quad%kernels(2,2) = allocate_kernel('H', 1)
      quad%kernels(3,2) = allocate_kernel('H', 1)
      quad%kernels(4,2) = allocate_kernel('H', 1)
      initialize = quad
   end function

   subroutine integrate(quad, elem_params, result)
      class(t_quad), intent(in) :: quad
      class(t_elem_params), intent(in) :: elem_params
      real(kind=rk), intent(out) ::  result(4,2)

      integer :: i, j
      class(t_kernel_ref), allocatable :: ref
      type(t_quad_result) :: res
      type(t_quad_params) :: par

      par = quad%quad_params
      do j = 1, 2
         do i = 1, 4
            ref = quad%kernels(i,j)
            ref%elem_params = elem_params
            res = integrate_kernel(ref, par)
            result(i,j) = res%anwser
         end do
      end do

   end subroutine

   function integrate_kernel(kernel_ref, quad_params)
      class(t_kernel_ref), intent(in) :: kernel_ref
      type(t_quad_params), intent(in) :: quad_params
      type(t_quad_result) :: integrate_kernel
      integrate_kernel = do_quad(f, quad_params)
   contains
      function f(x)
         real(kind=rk), intent(in) :: x
         real(kind=rk) :: f
         f = kernel_ref%eval(x)
      end function
   end function


end module
