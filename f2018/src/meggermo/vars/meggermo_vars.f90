module meggermo_vars

   use :: meggermo, only: rk

   implicit none
   private

   public :: t_vars

   type :: t_vars
      real(kind=rk), allocatable :: phi(:, :)
      real(kind=rk), allocatable :: phn(:, :)
   contains
      procedure :: element_view => vars_element_view
   end type

contains

   ! ---------------------------------
   ! Allocation functions
   ! ---------------------------------
   type(t_vars) function allocate_vars(nr_of_elemtents)
      !
      integer, intent(in):: nr_of_elemtents
      !
      real(kind=rk), allocatable :: phi(:, :)
      real(kind=rk), allocatable :: phn(:, :)
      !
      allocate (phi(2, 0:nr_of_elemtents + 2))
      allocate (phn(2, 0:nr_of_elemtents + 2))
      !
      allocate_vars = t_vars(phi, phn)
   end function

   ! ---------------------------------
   ! Element view functions
   ! ---------------------------------
   type(t_vars) function vars_element_view(vars, i)
      !
      class(t_vars), intent(in) :: vars
      integer, intent(in) :: i
      !
      vars_element_view = t_vars(vars%phi(:, i - 1:i + 2), vars%phn(:, i - 1:i + 2))
   end function

end module
