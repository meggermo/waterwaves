module meggermo_vars

   use, intrinsic :: iso_fortran_env, only: real64

   implicit none
   private

   public T_Vars

   type T_Vars
      real(kind=real64), allocatable :: phi(:, :)
      real(kind=real64), allocatable :: phn(:, :)
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
      real(kind=real64), allocatable :: phi(:, :)
      real(kind=real64), allocatable :: phn(:, :)
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
