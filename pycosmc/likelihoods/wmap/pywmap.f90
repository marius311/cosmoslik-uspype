
module pywmap

use wmap_likelihood_7yr
use WMAP_OPTIONS
use WMAP_UTIL

logical :: pywmap_initalized

contains

    subroutine WMAPInit()

        call wmap_likelihood_init
        pywmap_initalized = .true.

    end subroutine


    ! To call from Python use WMAPLnLike(cl_tt) where cl_tt[0] is ell=2
    function WMAPLnLike(cl_tt, cl_ee, cl_bb, cl_te, lmax)
      real WMAPLnLike
      real(8), dimension(2:lmax) :: cl_tt,cl_te,cl_ee,cl_bb
      real(8) :: like(num_WMAP)

      WMAPLnLike = 1e30
      if (lmax < ttmax) then
         print *, "WMAPLnLike: Please provide c_ell's up to", ttmax
         return
      end if
      if (.not. pywmap_initalized) then
         call WMAPInit()
      end if
         
      call wmap_likelihood_compute(cl_tt(2:ttmax),cl_te(2:ttmax),cl_ee(2:ttmax),cl_bb(2:ttmax),like)
      
      if (wmap_likelihood_ok) then
         WMAPLnLike = sum(like)
      else
         print *, "WMAPLnLike: Error computing WMAP likelihood."
      endif

    end function

end module
