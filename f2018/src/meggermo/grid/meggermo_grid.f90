module meggermo_grid

   use, intrinsic :: iso_fortran_env, only: rk => real64, output_unit
   use :: meggermo_interpolation, only:n_weights, dn_weights

   implicit none
   private

   public :: T_Grid, allocate_grid, create_linear_grid, create_cubic_grid

   type T_Grid
      real(rk), allocatable :: x(:, :)
      real(rk), allocatable :: n(:, :)
      real(rk), allocatable :: J(:)
      real(rk), allocatable :: K(:)
   contains
      procedure :: nr_of_elements => grid_nr_of_elements
      procedure :: element_view => grid_element_view
      procedure :: compute_geom => grid_compute_geom
      procedure :: glo_2_loc => grid_g2l
      procedure :: loc_2_glo => grid_l2g
      procedure :: print => grid_print
   end type

contains

   ! ---------------------------------
   ! Allocation functions
   ! ---------------------------------
   type(T_Grid) function allocate_grid(nr_of_elements)
      !
      integer, intent(in):: nr_of_elements
      !
      real(rk), allocatable :: x(:, :)
      real(rk), allocatable :: n(:, :)
      real(rk), allocatable :: J(:)
      real(rk), allocatable :: K(:)
      !
      allocate (x(2, 0:nr_of_elements + 2))
      allocate (n(2, 0:nr_of_elements + 2))
      allocate (J(0:nr_of_elements + 2))
      allocate (K(0:nr_of_elements + 2))
      !
      allocate_grid = T_Grid(x, n, J, K)
   end function

   ! ---------------------------------
   ! Element view functions
   ! ---------------------------------
   type(T_Grid) function grid_element_view(grid, i)
      !
      class(T_Grid), intent(in) :: grid
      integer, intent(in) :: i
      !
      grid_element_view = T_Grid(grid%x(:, i - 1:i + 2), grid%n(:, i - 1:i + 2), grid%J(i - 1:i + 2), grid%K(i - 1:i + 2))
   end function

   integer function grid_nr_of_elements(grid)
      !
      class(T_Grid), intent(inout) :: grid
      !
      grid_nr_of_elements = size(grid%x, 2) - 3
   end function

   ! ---------------------------------
   ! Grid generation routines
   ! ---------------------------------
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

      do i = 0, ne + 2
         grid%x(1:2, i) = x_b + (i - 1)*dx
         grid%J(i) = J
         grid%K(i) = 0.0
         grid%n(1, i) = n(1)
         grid%n(2, i) = n(2)
      end do

   end subroutine

   subroutine create_cubic_grid(x_b, x_e, d_x, grid)
      !
      real(rk), intent(in) :: x_b(2)
      real(rk), intent(in) :: x_e(2)
      real(rk), intent(in) :: d_x(2)
      class(T_Grid), intent(inout) :: grid
      !
      real(rk) :: t, dt, dx(2)
      integer :: i, ne
      !
      ne = grid%nr_of_elements()
      dx = x_e - x_b
      dt = 1.0/ne

      do i = 0, ne + 2
         t = dt*(i - 1)
         grid%x(1:2, i) = x_b + dx*cubic(t)
      end do

   contains
      function cubic(x) result(f)
         real(rk), intent(in) :: x
         real(rk) :: f
         f = d_x(1)*x + (3.0 - 2.0*d_x(1) - d_x(2))*x*x + (d_x(1) + d_x(2) - 2.0)*x*x*x
      end function

   end subroutine

   subroutine grid_compute_geom(grid)
      class(T_Grid), intent(inout) :: grid
      ! Local variables
      integer :: n
      !
      n = grid%nr_of_elements()
      call compute_kappa(grid%x, grid%K)
      call compute_jac_and_normal(grid%x, grid%K, grid%J, grid%n)

   contains

      subroutine compute_kappa(x, kappa)
         real(rk), allocatable, intent(in) :: x(:, :)
         real(rk), allocatable, intent(inout) :: kappa(:)
         ! Local variables
         integer :: i
         real(rk) :: dx1(2), dx2(2), arc_len(2)
         ! Assume arc lenght of virtual points are equal to adjacent
         ! element results in a normalization parameter (kappa) of zero.
         kappa(0) = 0.0
         kappa(n + 2) = 0.0
         ! Compute acr-length approximations and derive  kappa from it
         dx1 = x(:, 1) - x(:, 0)
         arc_len(1) = sqrt(dot_product(dx1, dx1))
         do i = 1, n + 1
            dx2 = x(:, i + 1) - x(:, i)
            arc_len(2) = sqrt(dot_product(dx2, dx2))
            kappa(i) = 2.0*arc_len(1)/sum(arc_len) - 1.0
            dx1 = dx2
            arc_len(1) = arc_len(2)
         end do
      end subroutine

      subroutine compute_jac_and_normal(x, kappa, jac, normal)
         real(rk), allocatable, intent(in) :: x(:, :)
         real(rk), allocatable, intent(in) :: kappa(:)
         real(rk), allocatable, intent(inout) :: jac(:)
         real(rk), allocatable, intent(inout) :: normal(:, :)
         ! Local variables
         integer :: i
         real(rk) :: dw(4), dx(2)

         i = 0
         call dn_weights(-2.0_rk, kappa(i), kappa(i + 1), dw)
         dx(1) = dot_product(dw, x(1, i:i + 3))
         dx(2) = dot_product(dw, x(2, i:i + 3))
         jac(i) = sqrt(dot_product(dx, dx))
         normal(1, i) = dx(2)/jac(i)
         normal(2, i) = -dx(1)/jac(i)

         do i = 1, n
            call dn_weights(-1.0_rk, kappa(i), kappa(i + 1), dw)
            dx(1) = dot_product(dw, x(1, i - 1:i + 2))
            dx(2) = dot_product(dw, x(2, i - 1:i + 2))
            jac(i) = sqrt(dot_product(dx, dx))
            normal(1, i) = dx(2)/jac(i)
            normal(2, i) = -dx(1)/jac(i)
         end do
         i = n + 1
         call dn_weights(1.0_rk, kappa(i), kappa(i + 1), dw)
         dx(1) = dot_product(dw, x(1, i - 2:i + 1))
         dx(2) = dot_product(dw, x(2, i - 2:i + 1))
         jac(i) = sqrt(dot_product(dx, dx))
         normal(1, i) = dx(2)/jac(i)
         normal(2, i) = -dx(1)/jac(i)

         i = n + 2
         call dn_weights(2.0_rk, kappa(i - 1), kappa(i), dw)
         dx(1) = dot_product(dw, x(1, i - 3:i))
         dx(2) = dot_product(dw, x(2, i - 3:i))
         jac(i) = sqrt(dot_product(dx, dx))
         normal(1, i) = dx(2)/jac(i)
         normal(2, i) = -dx(1)/jac(i)
      end subroutine

   end subroutine

   ! ---------------------------------
   ! Coordinate transformations
   ! ---------------------------------
   subroutine grid_g2l(grid, f_g, f_l)
      !
      class(T_Grid), intent(in) :: grid
      real(rk), intent(in), allocatable :: f_g(:, :)
      real(rk), intent(inout), allocatable :: f_l(:, :)
      !
      f_l(1, :) = grid%n(2, :)*f_g(1, :) + grid%n(1, :)*f_g(2, :)
      f_l(2, :) = grid%n(2, :)*f_g(2, :) - grid%n(1, :)*f_g(1, :)
   end subroutine

   subroutine grid_l2g(grid, f_l, f_g)
      !
      class(T_Grid), intent(in) :: grid
      real(rk), intent(in), allocatable :: f_l(:, :)
      real(rk), intent(inout), allocatable :: f_g(:, :)
      !
      f_g(1, :) = grid%n(2, :)*f_l(1, :) - grid%n(1, :)*f_l(2, :)
      f_g(2, :) = grid%n(2, :)*f_l(2, :) + grid%n(1, :)*f_l(1, :)
   end subroutine

   ! ---------------------------------
   ! Printing
   ! ---------------------------------
   subroutine grid_print(grid)
      class(T_Grid), intent(in) :: grid
      write (output_unit, '(2A18)') "X", "Y"
      write (output_unit, '(2E18.8)') grid%x
      write (output_unit, '(2A18)') "Nx", "Ny"
      write (output_unit, '(2E18.8)') grid%n
      write (output_unit, '(A18)') "Jac"
      write (output_unit, '(E18.8)') grid%J
      write (output_unit, '(A18)') "Kappa"
      write (output_unit, '(E18.8)') grid%K
   end subroutine

end module
