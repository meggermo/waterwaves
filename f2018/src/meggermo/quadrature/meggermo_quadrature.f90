module meggermo_quadrature

   use quadrature_module, wp => quadrature_wp

   implicit none

   type,extends(integration_class_1d),public :: log_kernel
   end type

   type,extends(integration_class_1d),public :: log_kernel_n2
   end type

   type,extends(integration_class_1d),public :: log_kernel_n3
   end type

   type,extends(integration_class_1d),public :: log_kernel_n4
   end type


   type integration_parameters
      real    :: xmin
      real    :: xmax
      integer :: steps
   end type

   type, abstract :: function_parameters
      procedure(eval), pointer, pass(params) :: feval
   end type

   abstract interface
      real function eval(x, params)
         import function_parameters
         real, intent(in) :: x
         class(function_parameters), intent(in) :: params
      end function
   end interface

   private

   public :: trapezoid

contains

   subroutine trapezoid(params, xmin, xmax, steps, result)

      interface
         real function f(s, ps)
            import :: function_parameters
            real, intent(in) :: s
            class(function_parameters), intent(in) :: ps
         end function f
      end interface

      class(function_parameters), intent(in) :: params
      real, intent(in) :: xmin
      real, intent(in) :: xmax
      integer, intent(in) :: steps
      real, intent(out) :: result

      integer :: i
      real :: x
      real :: dx

      dx = (xmax - xmin)/steps
      result = 0.5*(params%feval(xmin) + params%feval(xmax))
      do i = 2, steps
         x = xmin + (i - 1)*dx
         result = result + params%feval(x)
      end do

   end subroutine

end module
