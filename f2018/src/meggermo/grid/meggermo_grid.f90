module meggermo_grid

   use :: meggermo, only: rk
   use :: meggermo_interpolation, only: n_weights, dn_weights

   implicit none
   private

   public :: t_grid, t_gridtype

   type :: t_gridparams
      logical :: disable_kappa = .FALSE.
   end type

   type :: t_grid
      real(rk), allocatable :: x(:, :)
      real(rk), allocatable :: n(:, :)
      real(rk), allocatable :: J(:)
      real(rk), allocatable :: K(:, :)
      type(t_gridparams) :: grid_params
   contains
      procedure :: nr_of_elements => grid_nr_of_elements
      procedure :: element_view => grid_element_view
      procedure :: compute_geom => grid_compute_geom
      procedure :: apply_y_function => grid_apply_y_function
      procedure :: glo_2_loc => grid_g2l
      procedure :: loc_2_glo => grid_l2g
      procedure :: print => grid_print
   end type

   type, abstract:: t_gridtype
   contains
      procedure(s_generate_grid),deferred :: generate
      procedure :: initialize
   end type
   interface
      subroutine s_generate_grid(grid_type, grid)
         import t_grid, t_gridtype
         class(t_gridtype), intent(in) :: grid_type
         type(t_grid), intent(inout) :: grid
      end subroutine
   end interface


contains

   type(t_grid) function initialize(grid_type, nr_of_elements)
      class(t_gridtype), intent(in) :: grid_type
      integer, intent(in) :: nr_of_elements
      ! local variables
      type(t_grid) :: grid
      !
      grid = allocate_grid(nr_of_elements)
      call grid_type%generate(grid)
      call grid%compute_geom()
      initialize = grid
   end function

   ! ---------------------------------
   ! Allocation functions
   ! ---------------------------------
   type(t_grid) function allocate_grid(nr_of_elements)
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
      allocate_grid = t_grid(x, n, J, K)
   end function

   ! ---------------------------------
   ! Element view functions
   ! ---------------------------------
   type(t_grid) function grid_element_view(grid, i)
      !
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      !
      grid_element_view = t_grid(&
         grid%x(:, i - 1:i + 2), &
         grid%n(:, i - 1:i + 2), &
         grid%J(i - 1:i + 2), &
         grid%K(:, i:i), &
         grid%grid_params)
   end function

   integer function grid_nr_of_elements(grid)
      class(t_grid), intent(inout) :: grid
      !
      grid_nr_of_elements = size(grid%x, 2) - 3
   end function

   subroutine grid_apply_y_function(grid, f)
      interface
         function f(x)
            import rk
            real(rk), intent(in) :: x
            real(rk) :: f
         end function
      end interface
      class(t_grid), intent(inout) :: grid
      !
      integer :: i, ne
      ne = grid%nr_of_elements()
      do i = 0, ne + 2
         grid%x(2, i) = f(grid%x(1, i))
      end do
   end subroutine

   subroutine grid_compute_geom(grid)
      class(t_grid), intent(inout) :: grid
      !
      if(grid%grid_params%disable_kappa) then
         grid%K = 0.0
      else
         call compute_kappa(grid%x, grid%K)
      end if
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
      class(t_grid), intent(in) :: grid
      real(rk), intent(in), allocatable :: f_g(:, :)
      real(rk), intent(inout), allocatable :: f_l(:, :)
      !
      f_l(1, :) = grid%n(2, :)*f_g(1, :) + grid%n(1, :)*f_g(2, :)
      f_l(2, :) = grid%n(2, :)*f_g(2, :) - grid%n(1, :)*f_g(1, :)
   end subroutine

   subroutine grid_l2g(grid, f_l, f_g)
      !
      class(t_grid), intent(in) :: grid
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
      class(t_grid), intent(in) :: grid
      write (output_unit, '(2A18)') "X", "Y"
      write (output_unit, '(2E18.10)') grid%x
      write (output_unit, '(2A18)') "Nx", "Ny"
      write (output_unit, '(2E18.10)') grid%n
      write (output_unit, '(A18)') "Jac"
      write (output_unit, '(E18.10)') grid%J
      write (output_unit, '(2A18)') "K_1", "K_2"
      write (output_unit, '(2E18.10)') grid%K
   end subroutine

end module
