module meggermo_vars

   use :: meggermo, only: rk

   implicit none
   private

   public T_Vars

   type T_Vars
      real(kind=rk), allocatable :: phi(:, :)
      real(kind=rk), allocatable :: phn(:, :)
   contains
      procedure :: element_view => vars_element_view
   end type

contains

   ! ---------------------------------
   ! Allocation functions
   ! ---------------------------------
   type(T_Vars) function allocate_vars(nr_of_elemtents)
      !
      integer, intent(in):: nr_of_elemtents
      !
      real(kind=rk), allocatable :: phi(:, :)
      real(kind=rk), allocatable :: phn(:, :)
      !
      allocate (phi(2, 0:nr_of_elemtents + 2))
      allocate (phn(2, 0:nr_of_elemtents + 2))
      !
      allocate_vars = T_Vars(phi, phn)
   end function

   ! ---------------------------------
   ! Element view functions
   ! ---------------------------------
   type(T_Vars) function vars_element_view(vars, i)
      !
      class(T_Vars), intent(in) :: vars
      integer, intent(in) :: i
      !
      vars_element_view = T_Vars(vars%phi(:, i - 1:i + 2), vars%phn(:, i - 1:i + 2))
   end function

end module
