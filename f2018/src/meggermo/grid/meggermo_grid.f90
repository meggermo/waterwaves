module meggermo_grid

   use, intrinsic :: iso_fortran_env, only: rk => real64

   implicit none
   private

   public :: T_Grid, allocate_grid, create_linear_grid

   type T_Grid
      real(rk), allocatable :: x(:, :)
      real(rk), allocatable :: n(:, :)
      real(rk), allocatable :: J(:)
   contains
      procedure :: nr_of_elements => grid_nr_of_elements
      procedure :: element_view => grid_element_view
      procedure :: compute_geom => grid_compute_geom
   end type

contains

   subroutine create_linear_grid(x_b, x_e, grid)
      !
      real(rk), intent(in) :: x_b(2)
      real(rk), intent(in) :: x_e(2)
      class(T_Grid), intent(inout) :: grid
      !
      real(rk) :: dx(2), n(2), J
      integer :: i, ne
      !
      ne = grid%nr_of_elements()
      dx = (x_e - x_b)/ne
      J = sqrt(dot_product(dx, dx))
      n(1) = dx(2)/J
      n(2) = -dx(1)/J
      
      do i = -1, ne + 1
         grid%x(1:2, i + 1) = x_b + i*dx
         grid%J(i + 1) = J
         grid%n(1, i + 1) = n(1)
         grid%n(2, i + 1) = n(2)
      end do
      
   end subroutine

   ! ---------------------------------
   ! Allocation functions
   ! ---------------------------------
   type(T_Grid) function allocate_grid(nr_of_elements)
      !
      integer, intent(in):: nr_of_elements
      !
      real(rk), allocatable :: x(:, :)
      real(rk), allocatable :: J(:)
      real(rk), allocatable :: n(:, :)
      !
      allocate (x(2, 0:nr_of_elements + 2))
      allocate (n(2, 0:nr_of_elements + 2))
      allocate (J(0:nr_of_elements + 2))
      !
      allocate_grid = T_Grid(x, n, J)
   end function

   ! ---------------------------------
   ! Element view functions
   ! ---------------------------------
   type(T_Grid) function grid_element_view(grid, i)
      !
      class(T_Grid), intent(in) :: grid
      integer, intent(in) :: i
      !
      grid_element_view = T_Grid(grid%x(:, i - 1:i + 2), grid%n(:, i - 1:i + 2), grid%J(i - 1:i + 2))
   end function

   integer function grid_nr_of_elements(grid)
      !
      class(T_Grid), intent(inout) :: grid
      !
      grid_nr_of_elements = size(grid%x, 2) - 3
   end function

   subroutine grid_compute_geom(grid)
      !
      class(T_Grid), intent(inout) :: grid
      !
      integer :: i, n
      real(rk) :: x_t(2)
      !

      n = grid%nr_of_elements()

      do i = 1, n + 1
         x_t = 0.5*(grid%x(:, i + 1) - grid%x(:, i - 1))
         grid%J(i) = sqrt(dot_product(x_t, x_t))
         grid%n(1, i) = x_t(2)/grid%J(i)
         grid%n(2, i) = -x_t(1)/grid%J(i)
      end do

      i = 0
      grid%x(:, i) = 2.0*(grid%x(:, i + 1) - grid%x(:, i + 3)) + grid%x(:, i + 4)
      grid%n(:, i) = 2.0*(grid%n(:, i + 1) - grid%n(:, i + 3)) + grid%n(:, i + 4)
      grid%J(i) = 2.0*(grid%J(i + 1) - grid%J(i + 3)) + grid%J(i + 4)

      i = n + 2
      grid%x(:, i) = 2.0*(grid%x(:, i - 1) - grid%x(:, i - 3)) + grid%x(:, i - 4)
      grid%n(:, i) = 2.0*(grid%n(:, i - 1) - grid%n(:, i - 3)) + grid%n(:, i - 4)
      grid%J(i) = 2.0*(grid%J(i - 1) - grid%J(i - 3)) + grid%J(i - 4)

   end subroutine

end module
