module meggermo_integration

   use :: quadpack, only: dqage
   use :: meggermo, only: rk

   implicit none
   private

   public t_quad_params, t_quad_result, quad


   type t_quad_params
      real(rk) :: a = -1.0
      real(rk) :: b =  1.0
      real(rk) :: eps_abs = 1.0e-6
      real(rk) :: eps_rel = 1.0e-3
      integer  :: key = 1
      integer  :: limit = 8
   end type

   type t_quad_result
      real(rk) :: anwser
      real(rk) :: err_abs
      integer  :: i_er
      integer  :: n_eval
      integer  :: last
      type(t_quad_params) :: quad_params
   end type

contains

   function quad(f, quad_params)
      interface
         function f(x)
            import rk
            real(kind=rk), intent(in) :: x
            real(kind=rk) ::f
         end function
      end interface
      class(t_quad_params), intent(in) :: quad_params
      type(t_quad_result)  :: quad
      !
      type(t_quad_result) :: r
      real(rk), dimension(quad_params%limit) :: alist, blist, rlist, elist
      integer :: iord(quad_params%limit)

      r%quad_params = quad_params

      call dqage(f, &
         quad_params%a, quad_params%b, quad_params%eps_abs, quad_params%eps_rel, quad_params%key, quad_params%limit, &
         r%anwser, r%err_abs, r%n_eval, r%i_er, &
         alist, blist, rlist, elist, &
         iord, r%last)

      quad = r
   end function

end module
