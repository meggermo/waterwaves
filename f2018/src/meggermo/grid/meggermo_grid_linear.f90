module meggermo_grid_linear

   use :: meggermo, only: rk
   use :: meggermo_grid

   implicit none
   private

   public t_linegridtype

   type, extends(t_gridtype) :: t_linegridtype
      real(kind=rk) :: x_b(2)
      real(kind=rk) :: x_e(2)
      real(kind=rk) :: alpha(2) = (/1.0,1.0/)
   contains
      procedure :: generate
   end type

contains

   subroutine generate(grid_type,grid)
      class(t_linegridtype), intent(in) :: grid_type
      type(t_grid), intent(inout) :: grid
      !
      call create_grid(grid_type%x_b, grid_type%x_e, grid_type%alpha, grid)
   end subroutine


   subroutine create_grid(x_b, x_e, alpha, grid)
      !
      real(rk), intent(in) :: x_b(2)
      real(rk), intent(in) :: x_e(2)
      real(rk), intent(in) :: alpha(2)
      class(t_grid), intent(inout) :: grid
      !
      real(rk) :: t, dt, dx(2), a(3)
      integer :: i, ne
      !
      ne = grid%nr_of_elements()
      dx = x_e - x_b
      dt = 1.0/ne

      a(1) = alpha(1)
      a(2) = 3.0 - 2.0*alpha(1) - alpha(2)
      a(3) = alpha(1) + alpha(2) - 2.0

      do i = 0, ne + 2
         t = cubic(dt*(i - 1))
         grid%x(1, i) = x_b(1) + dx(1)*t
         grid%x(2, i) = x_b(2) + dx(2)*t
      end do

   contains
      function cubic(x) result(f)
         real(rk), intent(in) :: x
         real(rk) :: f
         f = a(1)*x + a(2)*x*x + a(3)*x*x*x
      end function

   end subroutine


end module
