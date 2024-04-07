module meggermo_network

   use :: meggermo_grid, only: t_grid
   use :: meggermo_vars, only: t_vars

   implicit none
   private

   public t_network, t_elementview

   type t_network
      type(t_grid) :: grid
      type(t_vars) :: vars
   contains
      procedure :: element_view => network_element_view
   end type

   type t_elementview
      integer :: element_index
      type(t_grid) :: grid
      type(t_vars) :: vars
   contains
      procedure :: x_element => get_x_element
   end type

contains

   ! ---------------------------------
   ! Allocation functions
   ! ---------------------------------
   type(t_network) function allocate_network(grid, vars)
      !
      class(t_grid), intent(in) :: grid
      class(t_vars), intent(in) :: vars
      !
      allocate_network = t_network(grid, vars)
   end function

   ! ---------------------------------
   ! Element view functions
   ! ---------------------------------
   type(t_elementview) function network_element_view(nw, i)
      !
      class(t_network), intent(in) :: nw
      integer, intent(in) :: i
      !
      network_element_view = t_elementview(i, nw%grid%element_view(i), nw%vars%element_view(i))
   end function

   function get_x_element(ev) result(x_e)
      class(t_elementview), intent(in):: ev
      real(kind=real64), dimension(4, 2) :: x_e
      x_e(:, 1) = ev%grid%x(1, 1:4)
      x_e(:, 2) = ev%grid%x(2, 1:4)
      x_e
   end function

end module
