module meggermo_network

   use, intrinsic :: iso_fortran_env, only: real64
   use :: meggermo_grid, only:T_Grid
   use :: meggermo_vars, only:T_Vars
   use :: meggermo_kernel, only :: ElemParams

   implicit none
   private

   public T_Network, T_ElementView

   type T_Network
      type(T_Grid) :: grid
      type(T_Vars) :: vars
   contains
      procedure :: element_view => network_element_view
   end type
   type T_ElementView
      integer :: element_index
      type(T_Grid) :: grid
      type(T_Vars) :: vars
   contains
      procedure :: x_element => get_x_element
   end type

contains

   ! ---------------------------------
   ! Allocation functions
   ! ---------------------------------
   type(T_Network) function allocate_network(grid, vars)
      !
      class(T_Grid), intent(in) :: grid
      class(T_Vars), intent(in) :: vars
      !
      allocate_network = T_Network(grid, vars)
   end function

   ! ---------------------------------
   ! Element view functions
   ! ---------------------------------
   type(T_ElementView) function network_element_view(nw, i)
      !
      class(T_Network), intent(in) :: nw
      integer, intent(in) :: i
      !
      network_element_view = T_ElementView(i, nw%grid%element_view(i), nw%vars%element_view(i))
   end function

   function get_x_element(ev) result(x_e)
      class(T_ElementView), intent(in):: ev
      real(kind=real64), dimension(4, 2) :: x_e
      x_e(:, 1) = ev%grid%x(1, 1:4)
      x_e(:, 2) = ev%grid%x(2, 1:4)
      x_e
   end function

end module
