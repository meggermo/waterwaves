module meggermo_grid

   use :: meggermo, only:rk
   use :: meggermo_interpolation, only:n_weights, dn_weights

   implicit none
   private

   public :: T_Grid, allocate_grid, create_linear_grid, create_cubic_grid

   type T_Grid
      real(rk), allocatable :: x(:, :)
      real(rk), allocatable :: n(:, :)
      real(rk), allocatable :: J(:)
      real(rk), allocatable :: K(:, :)
   contains
      procedure :: nr_of_elements => grid_nr_of_elements
      procedure :: element_view => grid_element_view
      procedure :: compute_geom => grid_compute_geom
      procedure :: apply_y_function => grid_apply_y_function
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
      real(rk), allocatable :: K(:, :)
      !
      allocate (x(2, 0:nr_of_elements + 2))
      allocate (n(2, 1:nr_of_elements + 1))
      allocate (J(1:nr_of_elements + 1))
      allocate (K(2, 1:nr_of_elements))
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
      grid_element_view = T_Grid(grid%x(:, i - 1:i + 2), grid%n(:, i - 1:i + 2), grid%J(i - 1:i + 2), grid%K(:, i:i))
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

      grid%x(1:2, 0) = x_b - dx
      do i = 1, ne
         grid%x(1:2, i) = x_b + (i - 1)*dx
         grid%J(i) = J
         grid%n(1, i) = n(1)
         grid%n(2, i) = n(2)
         grid%K(:, i) = 0.0
      end do
      grid%x(1:2, ne + 1) = x_b + ne*dx
      grid%J(ne + 1) = J
      grid%n(1, ne + 1) = n(1)
      grid%n(2, ne + 1) = n(2)
      grid%x(1:2, ne + 2) = x_b + (ne + 1)*dx

   end subroutine

   subroutine create_cubic_grid(x_b, x_e, d_x, grid)
      !
      real(rk), intent(in) :: x_b(2)
      real(rk), intent(in) :: x_e(2)
      real(rk), intent(in) :: d_x(2, 2)
      class(T_Grid), intent(inout) :: grid
      !
      real(rk) :: t, dt, dx(2), a(3, 2)
      integer :: i, ne
      !
      ne = grid%nr_of_elements()
      dx = x_e - x_b
      dt = 1.0/ne

      a(1, :) = d_x(:, 1)
      a(2, :) = 3.0 - 2.0*d_x(:, 1) - d_x(:, 2)
      a(3, :) = d_x(:, 1) + d_x(:, 2) - 2.0

      do i = 0, ne + 2
         t = dt*(i - 1)
         grid%x(1, i) = x_b(1) + dx(1)*cubic(t, a(:, 1))
         grid%x(2, i) = x_b(2) + dx(2)*cubic(t, a(:, 2))
      end do

   contains
      function cubic(x, alpha) result(f)
         real(rk), intent(in) :: x
         real(rk), intent(in) :: alpha(3)
         real(rk) :: f
         f = alpha(1)*x + alpha(2)*x*x + alpha(3)*x*x*x
      end function

   end subroutine

   subroutine grid_apply_y_function(grid, f)
      interface
         function f(x)
            import rk
            real(rk), intent(in) :: x
            real(rk) :: f
         end function
      end interface
      class(T_Grid), intent(inout) :: grid
      !
      integer :: i, ne
      ne = grid%nr_of_elements()
      do i = 0, ne + 2
         grid%x(2, i) = f(grid%x(1, i))
      end do
   end subroutine

   subroutine grid_compute_geom(grid)
      class(T_Grid), intent(inout) :: grid
      call compute_kappa(grid%x, grid%K)
      call compute_jac_and_normal(grid%x, grid%K, grid%J, grid%n)

   contains

      subroutine compute_kappa(x, kappa)
         real(rk), allocatable, intent(in) :: x(:, :)
         real(rk), allocatable, intent(inout) :: kappa(:, :)
         ! Local variables
         integer :: ie, ne
         real(rk) :: dx1(2), dx2(2), dx3(2), arc_len(3)
         ! Compute acr-length approximations and derive  kappa from it
         ne = grid%nr_of_elements()
         do ie = 1, ne
            dx1 = x(:, ie) - x(:, ie - 1)
            arc_len(1) = sqrt(dot_product(dx1, dx1))
            dx2 = x(:, ie + 1) - x(:, ie)
            arc_len(2) = sqrt(dot_product(dx2, dx2))
            dx3 = x(:, ie + 2) - x(:, ie + 1)
            arc_len(3) = sqrt(dot_product(dx3, dx3))
            ! Because the ratio of arc lenghts is what we need
            ! it is not so important that we use a crude approximation
            kappa(1, ie) = 2.0*arc_len(1)/sum(arc_len(1:2)) - 1.0
            kappa(2, ie) = 2.0*arc_len(3)/sum(arc_len(2:3)) - 1.0
         end do
      end subroutine

      subroutine compute_jac_and_normal(x, kappa, jac, normal)
         real(rk), allocatable, intent(in) :: x(:, :)
         real(rk), allocatable, intent(in) :: kappa(:, :)
         real(rk), allocatable, intent(inout) :: jac(:)
         real(rk), allocatable, intent(inout) :: normal(:, :)
         ! Local variables
         integer :: ie, ne
         real(rk) :: dw(4), dx(2)

         ne = grid%nr_of_elements()
         do ie = 1, ne
            call dn_weights(-1.0_rk, kappa(1, ie), kappa(2, ie), dw)
            write (*, '(I2,2E14.4,4E14.2)') ie, dx, dw
            dx(1) = dot_product(dw, x(1, ie - 1:ie + 2))
            dx(2) = dot_product(dw, x(2, ie - 1:ie + 2))
            jac(ie) = sqrt(dot_product(dx, dx))
            normal(1, ie) = dx(2)/jac(ie)
            normal(2, ie) = -dx(1)/jac(ie)
         end do
         ie = ne
         ! Use right side interpolation to approximate internal point of last element
         call dn_weights(1.0_rk, kappa(1, ie), kappa(2, ie), dw)
         dx(1) = dot_product(dw, x(1, ie - 1:ie + 2))
         dx(2) = dot_product(dw, x(2, ie - 1:ie + 2))
         jac(ne + 1) = sqrt(dot_product(dx, dx))
         normal(1, ne + 1) = dx(2)/jac(ne + 1)
         normal(2, ne + 1) = -dx(1)/jac(ne + 1)

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
      use, intrinsic :: iso_fortran_env, only: output_unit
      class(T_Grid), intent(in) :: grid
      write (output_unit, '(2A18)') "X", "Y"
      write (output_unit, '(2E18.8)') grid%x
      write (output_unit, '(2A18)') "Nx", "Ny"
      write (output_unit, '(2E18.8)') grid%n
      write (output_unit, '(A18)') "Jac"
      write (output_unit, '(E18.8)') grid%J
      write (output_unit, '(2A18)') "K_1", "K_2"
      write (output_unit, '(2E18.8)') grid%K
   end subroutine

end module
