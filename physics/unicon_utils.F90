!============================================================
! Utility routines for the UNIFIED CONVECTION SCHEME (UNICON)
!============================================================

module unicon_utils

use machine,  only: kind_phys
use wv_saturation2,   only : qsat, findsp
!use cam_abortutils,  only : endrun
!use cam_logfile,     only : iulog

implicit none
private
save

public :: &
   unicon_utils_init, &
   exnf,              &
   conden,            &
   slope,             &
   area_overlap,      &
   envcon_flux,       &
   prod_prep_up,      &
   evap_prep_dn,      &
   progup_thlqt,      &
   progup_uv,         &
   progup_wu2,        &
   compute_dp,        &
   buosort_downdraft, &
   compute_PDF,       &
   compute_epsdelnod, &
   buosorts_UW,       &
   positive_moisture, &
   positive_tracer,   &
   findsp_single

real(kind_phys), parameter :: alpha_max =  2._kind_phys       !  Upper limit of mixing parameter of updraft mass flux PDF [ no unit ]
real(kind_phys), parameter :: nonzero   =  1.e-20_kind_phys   !  Non-zero minimal positive constant [ no unit ]
real(kind_phys), parameter :: tmax_fice =  263.15_kind_phys   ! Temperature where ice starts to be formed [ K ]
real(kind_phys), parameter :: tmin_fice =  233.15_kind_phys   ! Temperature where ice fraction becomes 1  [ K ]

real(kind_phys) :: xlv       !  Latent heat of vaporization
real(kind_phys) :: xlf       !  Latent heat of fusion
real(kind_phys) :: xls       !  Latent heat of sublimation
real(kind_phys) :: cp        !  Specific heat of dry air
real(kind_phys) :: zvir      !  rh2o/rair - 1
real(kind_phys) :: r         !  Gas constant for dry air
real(kind_phys) :: g         !  Gravitational constant
real(kind_phys) :: p00       !  Reference pressure for exner function
real(kind_phys) :: rovcp     !  R/cp

real(kind_phys) :: droprad_liq    !  Effectie droplet radius of detrained liquid [ m ]  
real(kind_phys) :: droprad_ice    !  Effectie droplet radius of detrained    ice [ m ]
real(kind_phys) :: density_liq    !  Density of cloud liquid droplets [ kg/m3 ]  
real(kind_phys) :: density_ice    !  Density of cloud ice    crystals [ kg/m3 ]

integer  :: mclimit         !  If '1' ( '0' ), impose (not impose ) 'ql + qi > criqc' at the top interface
                            !  after precipitation fall-out.


!==================================================================================================
contains
!==================================================================================================
  
subroutine unicon_utils_init(&
   xlv_in, cp_in, xlf_in, zvir_in, r_in,                                 &
   g_in, droprad_liq_in, droprad_ice_in, density_liq_in, density_ice_in, &
   mclimit_in)

   real(kind_phys), intent(in) :: xlv_in     !  Latent heat of vaporization
   real(kind_phys), intent(in) :: xlf_in     !  Latent heat of fusion
   real(kind_phys), intent(in) :: cp_in      !  Specific heat of dry air
   real(kind_phys), intent(in) :: zvir_in    !  rh2o/rair - 1
   real(kind_phys), intent(in) :: r_in       !  Gas constant for dry air
   real(kind_phys), intent(in) :: g_in       !  Gravitational constant
   real(kind_phys), intent(in) :: droprad_liq_in    !  Effectie droplet radius of detrained liquid [ m ]  
   real(kind_phys), intent(in) :: droprad_ice_in    !  Effectie droplet radius of detrained    ice [ m ]
   real(kind_phys), intent(in) :: density_liq_in    !  Density of cloud liquid droplets [ kg/m3 ]  
   real(kind_phys), intent(in) :: density_ice_in    !  Density of cloud ice    crystals [ kg/m3 ]
   integer , intent(in) :: mclimit_in         !  If '1' ( '0' ), impose (not impose ) 'ql + qi > criqc' at the top interface

   xlv   = xlv_in
   xlf   = xlf_in
   xls   = xlv + xlf
   cp    = cp_in
   zvir  = zvir_in
   r     = r_in
   g     = g_in
   p00   = 1.e5_kind_phys
   rovcp = r/cp

   droprad_liq = droprad_liq_in
   droprad_ice = droprad_ice_in
   density_liq = density_liq_in
   density_ice = density_ice_in
   mclimit     = mclimit_in

end subroutine unicon_utils_init

!--------------------------------------------------------------------------------------------------

real(kind_phys) function exnf(pressure)
   real(kind_phys), intent(in)  :: pressure
   exnf = (pressure/p00)**rovcp
end function exnf

!--------------------------------------------------------------------------------------------------

subroutine conden(p,thl,qt,th,qv,ql,qi,rvls,id_check)

   ! --------------------------------------------------------------------- !
   ! Calculate thermodynamic properties from a given set of ( p, thl, qt ) !
   ! Note that this subroutine assumes horizontal homogeneity in the grid. ! 
   ! --------------------------------------------------------------------- !

   real(kind_phys), intent(in)  :: p, thl, qt
   real(kind_phys), intent(out) :: th, qv, ql, qi, rvls
   integer,  intent(out) :: id_check
   real(kind_phys)              :: es, qs, gam
   integer   i 
   real(kind_phys)  qc, t, fice, ficeg, ficeg0, leff, leffg 
   real(kind_phys)  f, fg

   ! ---------------------------------------------------------- !
   ! Main Computation Loop : Find Final Equilibrium Temperature !
   ! ---------------------------------------------------------- !

   id_check =  1
   ficeg0   = -1._kind_phys/(tmax_fice-tmin_fice)
   t        =  thl*exnf(p)
   call qsat(t, p, es, qs)

   if( qs .ge. qt ) then  

      qv = qt
      qc = 0._kind_phys
      ql = 0._kind_phys
      qi = 0._kind_phys
      th = t/exnf(p)
      rvls = qs
      id_check = 0
      return

   else 

      do i = 1, 10

         fice = max( 0._kind_phys, min( 1._kind_phys, (tmax_fice-t)/(tmax_fice-tmin_fice) ) )
          
         if( t .lt. tmin_fice-1._kind_phys ) then
            ficeg = 0._kind_phys
         elseif( t .ge. tmin_fice-1._kind_phys .and. t .lt. tmin_fice+1._kind_phys ) then
            ficeg =  ficeg0 * ( t - tmin_fice + 1._kind_phys )/2._kind_phys
         elseif( t .ge. tmin_fice+1._kind_phys .and. t .lt. tmax_fice-1._kind_phys ) then
            ficeg =  ficeg0
         elseif( t .ge. tmax_fice-1._kind_phys .and. t .lt. tmax_fice+1._kind_phys ) then
            ficeg = -ficeg0 * ( t - tmax_fice - 1._kind_phys )/2._kind_phys
         elseif( t .ge. tmax_fice+1._kind_phys ) then
            ficeg = 0._kind_phys
         endif

         leff   = fice *xls + (1._kind_phys - fice)*xlv
         leffg  = ficeg*(xls-xlv)
         call qsat(t, p, es, qs, gam=gam)
         f      = qt - (cp/leff)*(t-exnf(p)*thl)-qs
         fg     = (cp/leff)*((leffg/leff)*(t-exnf(p)*thl)-1._kind_phys-gam)

         if( abs(fg) .lt. 1.e-12_kind_phys ) then
            t = t + 0.1_kind_phys
         else 
            t = t - f/fg
         endif
         if( abs(f/fg) .lt. 1.e-3_kind_phys ) then
            qc = max(qt - qs,0._kind_phys)
            qv = qt - qc
            ql = qc*(1._kind_phys - fice)
            qi = fice*qc
            th = t/exnf(p)             
            rvls = qs
            id_check = 0
            return
         endif

      enddo

   end if

!   write(iulog,*) 'Warning : Convergence in conden is not achived and final value is used in unicon.F90'
   qc = max(qt - qs,0._kind_phys)
   qv = qt - qc
   ql = qc*(1._kind_phys - fice)
   qi = fice*qc
   th = t/exnf(p)             
   rvls = qs
   id_check = 0

end subroutine conden

!--------------------------------------------------------------------------------------------------

function slope(mkx,field,p0)

   ! ------------------------------------------------------------------ !
   ! Function performing profile reconstruction of conservative scalars !
   ! in each layer. This is identical to profile reconstruction used in !
   ! UW-PBL scheme but from bottom to top layer here.     At the lowest !
   ! layer near to surface, slope is defined using the two lowest layer !
   ! mid-point values. I checked this subroutine and it is correct.     !
   ! ------------------------------------------------------------------ !

   integer,  intent(in) :: mkx
   real(kind_phys), intent(in) :: field(mkx)
   real(kind_phys), intent(in) :: p0(mkx)
    
   real(kind_phys)             :: slope(mkx)
   real(kind_phys)             :: below
   real(kind_phys)             :: above
   integer              :: k

   below = (field(2) - field(1))/(p0(2) - p0(1))
   do k = 2, mkx
      above = (field(k) - field(k-1))/(p0(k) - p0(k-1))
      if (above .gt. 0._kind_phys) then
         slope(k-1) = max(0._kind_phys,min(above,below))
      else 
         slope(k-1) = min(0._kind_phys,max(above,below))
      end if
      below = above
   end do
   slope(mkx) = slope(mkx-1)

   ! Sep.22.2011. Set the slope in the lowest model layer to be zero to reduce the
   !              sensitivity of the diurnal cycle to the surface heating.
   ! slope(1) = 0._kind_phys 

end function slope

!--------------------------------------------------------------------------------------------------

function area_overlap(x1,y1,a1,x2,y2,a2,cn)

   ! ----------------------------------------------------------------- !
   ! Function to compute overlapping area between two disks located at !
   ! (x1,y1) in unit of [m] with fractional area a1, and               !
   ! (x2,y2) in unit of [m] with fractional area a2, both of which has !
   ! the same number concentration of 'cn' [ # / m^2 ].                !
   ! The resulting overlapping area has no unit : fractional area.     ! 
   ! ----------------------------------------------------------------- !

   real(kind_phys)             :: area_overlap
   real(kind_phys), intent(in) :: x1, y1, a1, x2, y2, a2, cn
   real(kind_phys)             :: r1, r2
   real(kind_phys)             :: rmin, rmax, d
   real(kind_phys)             :: arg1, arg2, arg3

   r1 = sqrt(max(0._kind_phys,a1/cn/3.141592_kind_phys))
   r2 = sqrt(max(0._kind_phys,a2/cn/3.141592_kind_phys))
   d = sqrt( (x2-x1)**2._kind_phys + (y2-y1)**2._kind_phys )
   rmin = min( r1, r2 )
   rmax = max( r1, r2 )

   if( rmin .eq. 0._kind_phys .or. rmax .eq. 0._kind_phys ) then
      area_overlap = 0._kind_phys
      return
   else
      if( d .le. ( rmax - rmin ) ) then
         area_overlap = min( a1, a2 )
         return
      elseif( d .ge. ( rmax + rmin ) ) then
         area_overlap = 0._kind_phys
         return
      else
         arg1 = (d**2._kind_phys+rmin**2._kind_phys-rmax**2._kind_phys)/(2._kind_phys*d*rmin)    
         arg2 = (d**2._kind_phys+rmax**2._kind_phys-rmin**2._kind_phys)/(2._kind_phys*d*rmax)    
         arg3 = (-d+rmin+rmax)*(d+rmin-rmax)*(d-rmin+rmax)*(d+rmin+rmax)
         arg1 = max(-1._kind_phys,min(1._kind_phys,arg1))
         arg2 = max(-1._kind_phys,min(1._kind_phys,arg2))
         arg3 = max(0._kind_phys,arg3)

         area_overlap = cn * &
                        (rmin**2._kind_phys*acos(arg1) + rmax**2._kind_phys*acos(arg2) - 0.5_kind_phys*sqrt(arg3))

         ! Apr.25.2012. I checked that below is always satisfied with round-off error.
         ! So, I safely added below safety constraint.
         area_overlap = min( a1, min( a2, area_overlap ) )
         return 
      end if
   end if
end function area_overlap

!--------------------------------------------------------------------------------------------------

subroutine envcon_flux(ki,mkx,umi,dmi,a0,ssa0,ps0,au,ad)

   ! --------------------------------------------------------------------------- !
   ! Compute mean-environmental values of conservative scalar for computation of !
   ! convective fluxes by considering the displacement of flux interface induced !
   ! by convective updraft and downdraft mass fluxes and associated compensating !
   ! downwelling and upwelling.                                                  !
   ! ki  : interface index that is considered                                    !
   ! umi : updraft   mass flux in unit of [Pa] during dt ( umi >= 0 )            !
   ! dmi : downdraft mass flux in unit of [Pa] during dt ( dmi >= 0 )            !
   ! a   : environmental conservative scalar that is considered                  !
   ! Done.                                                                       !
   ! --------------------------------------------------------------------------- ! 

   integer,  intent(in)   :: ki, mkx 
   real(kind_phys), intent(in)   :: umi, dmi
   real(kind_phys), intent(in)   :: a0(mkx), ssa0(mkx)
   real(kind_phys), intent(in)   :: ps0(0:mkx)
   real(kind_phys), intent(out)  :: au, ad
   integer   k, ku, kd 
   real(kind_phys)  um, dm
   real(kind_phys)  dp, dpu, dpd, pbot, ptop, dptop, a_dptop, dpbot, a_dpbot

   ! Impose a limiting on the updraft (um) and downdraft mass flux (dm) such that
   ! it cannot be larger than the available mass above the interface (um) and 
   ! below the displaced interface (dm) by updraft mass flux. Note that ps0(0) is
   ! surface interface while ps0(mkx) is top-most interface. Note umi, dmi > 0.

   um = max( 0._kind_phys, min( umi, ps0(ki) - ps0(mkx) ) )
   dm = max( 0._kind_phys, min( dmi, ps0(0)  - ps0(ki) + um ) ) 

   ! Treatment of updraft
 
   ! if( um .eq. 0._kind_phys ) then
   if( um .lt. 1.e-5_kind_phys ) then ! To avoid dividing by zero ( dpu = 0 ) by round-off error.
      if( ki .eq. mkx ) then
         au = a0(ki)
      else
         au = a0(ki+1) + 0.5_kind_phys * ssa0(ki+1) * ( ps0(ki) - ps0(ki+1) ) 
      endif
      goto 50
   endif

   ku = ki + 1
   do k = ki, mkx
      if( ps0(k) .lt. ( ps0(ki) - um ) ) then  
         ku = k 
         goto 10
      endif
   enddo
10 continue

   au = 0._kind_phys
   dpu = 0._kind_phys
   if( ( ku - 1 ) .ge. ( ki + 1 ) ) then
      do k = ki + 1, ku - 1 
         dp = ps0(k-1) - ps0(k)
         au = au + a0(k) * dp
         dpu = dpu + dp
      enddo
   endif

   ptop = ps0(ki) - um
   dptop = ps0(ku-1) - ptop
   a_dptop = a0(ku) + 0.5_kind_phys * ssa0(ku) * ( ptop - ps0(ku) )
   au = au + a_dptop * dptop
   dpu = dpu + dptop
   ! I checked that dpu = 0 happans when umi is very small, 1.e-15.
   if( dpu .eq. 0._kind_phys ) then
!      write(iulog,*) 'ki, ku, um, umi, dmi = ', ki, ku, um, umi, dmi
!      write(iulog,*) 'ptop, dptop, a_dptop, au, dpu = ', ptop, dptop, a_dptop, au, dpu
!      do k = 1, mkx 
!         write(iulog,*) 'ps0(k), a0(k) =', ps0(k), a0(k)
!      enddo
!      call endrun('UNICON : Zero dpu within envcon_flux')  
       call abort
   endif
   au = au / dpu

50 continue

   ! Treatment of downdraft
    
   ! if( dm .eq. 0._kind_phys ) then
   if( dm .lt. 1.e-5_kind_phys ) then ! To avoid dividing by zero ( dpd = 0 ) by round-off error.
      ad = a0(ku) + ssa0(ku) * ( ptop - 0.5_kind_phys * ( ps0(ku-1) - ps0(ku ) ) )
      return
   endif

   pbot = ps0(ki) - um + dm
   kd = ku
   do k = ku, 1, -1
      if( ps0(k) .ge. pbot ) then
         kd = k + 1
         goto 20
      endif
   enddo
20 continue

   ad = 0._kind_phys
   dpd = 0._kind_phys
   if( ( ku - 1 ) .ge. ( kd + 1 ) ) then
      do k =  kd + 1, ku - 1
         dp = ps0(k-1) - ps0(k)
         ad = ad + a0(k) * dp
         dpd = dpd + dp
      enddo
   endif

   if( pbot .le. ps0(ku-1) ) then
      dpbot = dm
      a_dpbot = a0(ku) + 0.5_kind_phys * ssa0(ku) * ( pbot + ptop - ps0(ku-1) - ps0(ku) )
      ad = ad + a_dpbot * dpbot
      dpd = dpd + dpbot
      ad = ad / dpd   
      return
   else
      dpbot = pbot - ps0(kd)
      a_dpbot = a0(kd) + 0.5_kind_phys * ssa0(kd) * ( pbot + ps0(kd) - ps0(kd-1) - ps0(kd) )
      ad = ad + a_dpbot * dpbot + a_dptop * dptop
      dpd = dpd + dpbot + dptop
      ad = ad / dpd
      return
   endif

end subroutine envcon_flux

!--------------------------------------------------------------------------------------------------

  subroutine prod_prep_up( z_b, z_t, p_b, p_t, exn_t, exn_m, w_b, w_t,                   & 
                           thl_in, qt_in, ql_in, qi_in, tr_in,                           &
                           S_b_ql_in, S_b_qi_in, iprd_prep,                              &
                           ql_b, qi_b, epsb,                                             &
                           thl_m, ssthl_m, thl_b, qt_m, ssqt_m, qt_b,                    & 
                           ncnst, ixcldliq, ixcldice, ixnumliq, ixnumice, ii, kk, lchnk, &
                           flxrain, flxsnow, a_p, a_u, a_pu,                             &
                           caer, criqc, c0_ac,                                           &
                           exql, exqi, extr, S_t_ql, S_t_qi, evpR, evpS, evpRStr )
  ! ------------------------------------------------------------------------------------------------------- ! 
  ! Compute 'exql, exqi >= 0' [kg/kg] and 'S_t_ql, S_t_qi >= 0' [kg/kg/Pa].                                 !
  ! ------------------------------------------------------------------------------------------------------- !

    integer,  intent(in)     :: ncnst, ixcldliq, ixcldice, ixnumliq, ixnumice, ii, kk, lchnk, iprd_prep 
    real(kind_phys), intent(in)     :: z_b, z_t, p_b, p_t, exn_t, exn_m
    real(kind_phys), intent(in)     :: w_b, w_t
    real(kind_phys), intent(in)     :: thl_in, qt_in, ql_in, qi_in, tr_in(ncnst)
    real(kind_phys), intent(in)     :: S_b_ql_in, S_b_qi_in
    real(kind_phys), intent(in)     :: ql_b, qi_b, epsb
    real(kind_phys), intent(in)     :: thl_m, ssthl_m, thl_b, qt_m, ssqt_m, qt_b
    real(kind_phys), intent(in)     :: flxrain, flxsnow
    real(kind_phys), intent(in)     :: a_p, a_u, a_pu
    real(kind_phys), intent(in)     :: caer, criqc, c0_ac 
    real(kind_phys), intent(out)    :: exql, exqi, extr(ncnst)
    real(kind_phys), intent(out)    :: evpR, evpS, evpRStr(ncnst)
    real(kind_phys), intent(out)    :: S_t_ql, S_t_qi
    integer   mt, iter, id_exit, id_check, niter
    real(kind_phys)  lambda
    real(kind_phys)  tmp1, tmp2
    real(kind_phys)  tmp_thl, tmp_qt, tmp_th, tmp_qv, tmp_ql, tmp_qi, tmp_qs
    real(kind_phys)  dp, dz, wm, delta_t
    real(kind_phys)  flxrain_in, flxsnow_in
    real(kind_phys)  ql, qi
    real(kind_phys)  S_b_ql, S_b_qi 
    real(kind_phys)  S_ql, S_qi 
    real(kind_phys)  S_t_ql_pre, S_t_qi_pre
    real(kind_phys)  exql_pre, exqi_pre
    real(kind_phys)  dia_thl, dia_qt
    ! ----------------------- !
    ! Compute basic variables !
    ! ----------------------- !

    niter = 1
    if( iprd_prep .eq. -1 .or. iprd_prep .eq. -5 ) niter = 0   ! Forward Method
    lambda     = 0.5_kind_phys 
    dp         = p_b - p_t              ! [  Pa ] >= 0.
    dz         = z_t - z_b              ! [  z  ] >= 0.
    wm         = 0.5_kind_phys * ( w_b + w_t ) ! [ m/s ] >  0.
    delta_t    = dz / wm                ! [  s  ] >= 0.
    flxrain_in = flxrain / max( nonzero, a_p )
    flxsnow_in = flxsnow / max( nonzero, a_p )

    ! ------------------------------------------------------------------------------------------------------------------- !
    ! Current formulation only contains a simple auto-conversion process as a unique function of in-cloud LWC/IWC.        !
    ! In future, I should add more advanced formula as a function of droplet radius, vertical velocity, and               !
    ! precipitation flux falling into the current layer.                                                                  !
    ! IMPORTANT :                                                                                                         !
    !    (1) Current precipitation formula conserve 'the fraction (f(T)) of in-cumulus liquid / ice as a function         !
    !        of T' as is assumed in the subroutine 'conden'. Thus, without using 'conden', I can construct this iteration !
    !        loop associated with precipitation production (i.e., temperature T does not change during precipitation      !
    !        production process ), which saves computation time a lot.                                                    !
    !    (2) With future double-moment microphysics, this f(T) assumed within 'conden' does not hold after precipitation  !
    !        production. Thus, I should not use 'conden' subroutine with future microphysics, which will also save        !
    !        computation time in future.                                                                                  !
    ! ------------------------------------------------------------------------------------------------------------------- !

    if( iprd_prep .eq. -5 ) then

        if( ( ql_b + qi_b ) .gt. criqc ) then
            S_ql = c0_ac * ( ( ql_b + qi_b ) - criqc ) * ( ql_b / ( ql_b + qi_b ) )
            S_qi = c0_ac * ( ( ql_b + qi_b ) - criqc ) * ( qi_b / ( ql_b + qi_b ) )
        else
            S_ql = 0._kind_phys
            S_qi = 0._kind_phys
        endif      
        dia_qt  =   S_ql + S_qi 
        dia_thl = - ( ( xlv / cp / exn_m ) * S_ql + ( xls / cp / exn_m ) * S_qi )
        call progup_thlqt( epsb, 0._kind_phys, dia_thl, p_b, p_t, thl_m, ssthl_m, thl_b, tmp_thl )
        call progup_thlqt( epsb, 0._kind_phys, dia_qt,  p_b, p_t, qt_m,  ssqt_m,  qt_b,  tmp_qt  )
        call conden( p_t, tmp_thl, tmp_qt, tmp_th, tmp_qv, tmp_ql, tmp_qi, tmp_qs, id_check )
        exql = min( max( ql_in - tmp_ql, 0._kind_phys ), 0.99_kind_phys * ql_in )
        exqi = min( max( qi_in - tmp_qi, 0._kind_phys ), 0.99_kind_phys * qi_in )
        if( mclimit .eq. 1 ) then
            tmp1 = exql + exqi
            tmp2 = min( tmp1, max( ql_in + qi_in - criqc, 0._kind_phys ) ) ! To impose a continuous variation across ql + qi = criqc.
            exql = exql * ( tmp2 / max( tmp1, nonzero ) )
            exqi = exqi * ( tmp2 / max( tmp1, nonzero ) )
        endif
        S_ql = exql / max( dp, nonzero )
        S_qi = exqi / max( dp, nonzero )
        S_t_ql = S_ql
        S_t_qi = S_qi

    else

    ! ------------------------------------------------------------------------------------------------------ !
    ! Compute initially-precipitated updraft state variable 'ql,qi' at the top interface using precipitation !
    ! tendency at the base interface.                                                                        !
    ! ------------------------------------------------------------------------------------------------------ !

    S_ql = S_b_ql_in
    S_qi = S_b_qi_in
    exql = S_ql * dp 
    exqi = S_qi * dp
    exql = min( max( 0._kind_phys, exql ), 0.99_kind_phys * ql_in )
    exqi = min( max( 0._kind_phys, exqi ), 0.99_kind_phys * qi_in )
    if( mclimit .eq. 1 ) then
        tmp1 = exql + exqi
        tmp2 = min( tmp1, max( ql_in + qi_in - criqc, 0._kind_phys ) ) ! To impose a continuous variation across ql_in + qi_in = criqc.
        exql = exql * ( tmp2 / max( tmp1, nonzero ) )
        exqi = exqi * ( tmp2 / max( tmp1, nonzero ) )
    endif
    S_b_ql = exql / max( dp, nonzero )
    S_b_qi = exqi / max( dp, nonzero )
    ql = ql_in - exql
    qi = qi_in - exqi
    if( iprd_prep .eq. 1 ) then
        ql = ql_in
        qi = qi_in
    endif    
    if( ( ql + qi ) .gt. criqc ) then
        S_t_ql = c0_ac * ( ( ql + qi ) - criqc ) * ( ql / ( ql + qi ) ) ! [ kg/kg/Pa ]
        S_t_qi = c0_ac * ( ( ql + qi ) - criqc ) * ( qi / ( ql + qi ) ) ! [ kg/kg/Pa ]
    else
        S_t_ql = 0._kind_phys
        S_t_qi = 0._kind_phys
    endif

    ! ---------------------------------------------------------------------------------------------------------------- !
    ! Perform implicit iteration                                                                                       !
    ! The requires output from the iteration loop :                                                                    !
    !   (1) 'exql, exqi'                                                                                               ! 
    !   (2) 'S_t_ql, S_t_qi'                                                                                           !
    ! where the above (1) and (2) are fully consistent ( exql = 0.5_kind_phys * ( S_b_ql + S_t_ql ) * dp, exqi = ... )        !
    ! regardless of the iteration number. However, the consistency between the 'ql = ql_in - exql' and 'S_t_ql' at the !
    ! top interface can be obtained as iteration is executed many times.                                               !
    ! By setting 'do iter = 1, 1', I can use centered difference instead of forward difference which is a default.     !
    ! ---------------------------------------------------------------------------------------------------------------- !

    id_exit = 0
    do iter = 1, niter
    
       ! ------------------------------------------------------ !
       ! Compute 'raw' precipitation rate at the top interface. !
       ! ------------------------------------------------------ !

       if( ( ql + qi ) .gt. criqc ) then
           S_t_ql = c0_ac * ( ( ql + qi ) - criqc ) * ( ql / ( ql + qi ) ) ! [ kg/kg/Pa ]
           S_t_qi = c0_ac * ( ( ql + qi ) - criqc ) * ( qi / ( ql + qi ) ) ! [ kg/kg/Pa ]
       else
           S_t_ql = 0._kind_phys
           S_t_qi = 0._kind_phys
       endif

       ! ------------------------------------------------------------------------------- !
       ! Impose a limiter on the computed 'raw' precipitation rate at the top interface. !
       ! Use 'ql_in, qi_in' which does not include precipitation fall-out at the top.    !
       ! At the end of this block, 'exql,exqi' is fully consistent with 'S_t_ql,S_t_qi'. ! 
       ! ------------------------------------------------------------------------------- !  

       S_ql = 0.5_kind_phys * ( S_b_ql + S_t_ql )
       S_qi = 0.5_kind_phys * ( S_b_qi + S_t_qi )
       if( iprd_prep .eq. 1 ) then
           S_ql = S_t_ql
           S_qi = S_t_qi
       endif
       exql = S_ql * dp 
       exqi = S_qi * dp 
       exql = min( max( 0._kind_phys, exql ), 0.99_kind_phys * ql_in )
       exqi = min( max( 0._kind_phys, exqi ), 0.99_kind_phys * qi_in )
       if( mclimit .eq. 1 ) then
           tmp1 = exql + exqi
           tmp2 = min( tmp1, max( ql_in + qi_in - criqc, 0._kind_phys ) ) ! To impose a continuous variation across ql_in + qi_in = criqc.
           exql = exql * ( tmp2 / max( tmp1, nonzero ) )
           exqi = exqi * ( tmp2 / max( tmp1, nonzero ) )
       endif
       S_ql = exql / max( dp, nonzero )
       S_qi = exqi / max( dp, nonzero )
       S_t_ql = 2._kind_phys * S_ql - S_b_ql ! IMPORTANT : This must be allowed to be negative.
       S_t_qi = 2._kind_phys * S_qi - S_b_qi ! IMPORTANT : This must be allowed to be negative.  	
       if( iprd_prep .eq. 1 ) then
           S_t_ql = 1._kind_phys * S_ql      ! IMPORTANT : This must be allowed to be negative.
           S_t_qi = 1._kind_phys * S_qi      ! IMPORTANT : This must be allowed to be negative.  	
       endif

       ! ------------------------------------------------------------------------------------------------------------- !
       ! Compute 'implicit' precipitation rate at the top interface by averaging the 'current' precipitation rate with !
       ! the 'previous' precipitation rate computed at the previous iteration loop.                                    !
       ! Since both 'current' and 'previous' precipitation rates satisfies the limiters,                               !
       ! the average 'implicit' precipitation rate also satisfies the limiter automatically.                           !
       ! Within if block, the implicit 'exql,exqi' is fully consistent with 'S_t_ql,S_t_qi'.                           ! 
       ! ------------------------------------------------------------------------------------------------------------- !

       if( iter .gt. 1 ) then
           S_t_ql = lambda * S_t_ql + ( 1._kind_phys - lambda ) * S_t_ql_pre
           S_t_qi = lambda * S_t_qi + ( 1._kind_phys - lambda ) * S_t_qi_pre
           S_ql   = 0.5_kind_phys * ( S_b_ql + S_t_ql )
           S_qi   = 0.5_kind_phys * ( S_b_qi + S_t_qi )
           if( iprd_prep .eq. 1 ) then
               S_ql   = S_t_ql
               S_qi   = S_t_qi
           endif
           exql   = S_ql * dp
           exqi   = S_qi * dp
         ! if( kk .ge. 10 .and. kk .le. 12 ) then
         !     write(6,*)
         !     write(6,*) 'UNICON : Convergence test within the subroutine prod_prep_up'
         !     write(6,*) 'kk, iter, abs( exql + exqi - exql_pre - exqi_pre ) = ', kk, iter, abs( exql + exqi - exql_pre - exqi_pre )
         !     write(6,*) 'kk, iter, exql, exqi, S_t_ql, S_t_qi, dp           = ', kk, iter, exql, exqi, S_t_ql, S_t_qi, dp  
         !     write(6,*)
         ! endif  
           if( abs( exql + exqi - exql_pre - exqi_pre ) .lt. 1.e-6_kind_phys ) then
               id_exit = 1
           endif 
       endif

       S_t_ql_pre = S_t_ql 
       S_t_qi_pre = S_t_qi 
       exql_pre   = exql 
       exqi_pre   = exqi 

       ! ----------------------------------------------------------------------------------------------------------------- !
       ! Update state variable at the top interface.                                                                       ! 
       ! At this stage, 'exql = 0.5_kind_phys * ( S_b_ql + S_t_ql ) * dp' is exactly satisfied.                                   !
       ! However, our 'S_t_ql' becomes inconsistent with the below updated 'ql = ql_in - exql' since 'S_t_ql' was computed !
       ! using 'ql = ql_in - exql(old)' where 'exql(old)' differs from 'exql'.                                             ! 
       ! This inconsistency will be removed as iteration goes on.                                                          !
       ! ----------------------------------------------------------------------------------------------------------------- !

       ql = ql_in - exql
       qi = qi_in - exqi
       if( id_exit .eq. 1 ) goto 10

    enddo
 10 S_t_ql = max( 0._kind_phys, S_t_ql ) ! Reset to non-negative value before sending to output.
    S_t_qi = max( 0._kind_phys, S_t_qi ) ! Reset to non-negative value before sending to output.

    endif ! End of 'iprd_prep = -5' choice

    ! ------------------------------- !
    ! Treatment of In-Cumulus Tracers !
    ! ------------------------------- !

    do mt = 1, ncnst
       if( mt .eq. 1 ) then
           extr(mt) = 0._kind_phys
       elseif( mt .eq. ixcldliq ) then
           extr(mt) = exql
       elseif( mt .eq. ixcldice ) then
           extr(mt) = exqi
       elseif( mt .eq. ixnumliq ) then
           extr(mt) = exql * 3._kind_phys / ( 4._kind_phys * 3.141592_kind_phys * droprad_liq**3 * density_liq )
       elseif( mt .eq. ixnumice ) then
           extr(mt) = exqi * 3._kind_phys / ( 4._kind_phys * 3.141592_kind_phys * droprad_ice**3 * density_ice )
       else
           ! ----------------------------------------------------------------------------------------- !
           ! Wet deposition of aerosols (both interstitial and cloud-borne) within convective updarft. !
           ! Below is a very simple treatment which should be refined in future.                       !
           ! Note that I should use 'qt_in, tr_in' (i.e., input values) within below block.            !
           ! ----------------------------------------------------------------------------------------- !
           extr(mt) = tr_in(mt) * ( ( exql + exqi ) / max( qt_in, nonzero ) )
         ! Nov.26.2013. Following the reviewer's comments, set 'extr(mt) = 0' since current formulation of
         !              extr(mt) only treats auto-conversion not accretion.
         ! extr(mt) = 0._kind_phys 
         ! Nov.29.2013. Following the reviewer's comments, use 'ql_in + qi_in' instead of 'qt_in' 
         !              in computing 'extr(mt)' above.
           extr(mt) = caer * tr_in(mt) * min( 1._kind_phys, ( ( exql + exqi ) / max( ql_in + qi_in, nonzero ) ) )
       endif
    enddo

    ! ------------------------------------------------------------ !
    ! Evaporation.                                                 !
    ! Temporary set it to be zero, but should be refined in future !
    ! ------------------------------------------------------------ !

    evpR = 0._kind_phys
    evpS = 0._kind_phys
    do mt = 1, ncnst
       evpRStr(mt) = 0._kind_phys
    enddo

    return

  end subroutine prod_prep_up

!--------------------------------------------------------------------------------------------------

  subroutine evap_prep_dn( z_b, z_t, p_b, p_t, w_dt, bogtop,                                                     &  
                           th_in, qv_in, ql_in, qi_in, tr_in, qmin,                                              &
                           S_t_qvR_in, S_t_qvS_in, ievp_prep,                                                    & 
                           flxrain_bot_upeesm, flxsnow_bot_upeesm, flxtrrs_bot_upeesm, a_p_msfc,                 &
                           ncnst, ixcldliq, ixcldice, ixnumliq, ixnumice, ndb_evp, cmfdb_evp, ii, kk, ks, lchnk, & 
                           rho, thv_mean_b, cmf_db, eps_dn, del_dn,                                              &
                           kevp_rain_dn, kevp_snow_dn, eta2, rbuoy_dn, rdrag, rjet, nonzero, wdmin,              &
                           evp_qvR, evp_qvS, evp_tr, S_b_qvR, S_b_qvS, w_db )
  ! ------------------------------------------------------------------------------------------------------- ! 
  ! Compute 'evp_qvR, evp_qvS >= 0' [kg/kg] and 'S_b_qvR, S_b_qvS >= 0' [kg/kg/Pa], 'w_db >= wdmin [m/s]'   !
  ! ------------------------------------------------------------------------------------------------------- !
    implicit none
    integer,  intent(in)     :: ncnst, ixcldliq, ixcldice, ixnumliq, ixnumice
    integer,  intent(in)     :: ndb_evp, ii, kk, ks, lchnk, ievp_prep
    real(kind_phys), intent(in)     :: cmfdb_evp
    real(kind_phys), intent(in)     :: z_b, z_t, p_b, p_t
    real(kind_phys), intent(in)     :: w_dt, bogtop
    real(kind_phys), intent(in)     :: th_in, qv_in, ql_in, qi_in, tr_in(ncnst), qmin(ncnst)
    real(kind_phys), intent(in)     :: S_t_qvR_in, S_t_qvS_in
    real(kind_phys), intent(in)     :: flxrain_bot_upeesm, flxsnow_bot_upeesm, flxtrrs_bot_upeesm(ncnst), a_p_msfc
    real(kind_phys), intent(in)     :: rho, thv_mean_b, cmf_db, eps_dn, del_dn 
    real(kind_phys), intent(in)     :: kevp_rain_dn, kevp_snow_dn, eta2 
    real(kind_phys), intent(in)     :: rbuoy_dn, rdrag, rjet
    real(kind_phys), intent(in)     :: nonzero, wdmin
    real(kind_phys), intent(out)    :: evp_qvR, evp_qvS, evp_tr(ncnst)
    real(kind_phys), intent(out)    :: S_b_qvR, S_b_qvS
    real(kind_phys), intent(out)    :: w_db
    real(kind_phys)                 :: es, qs
    integer   mt, iter, id_exit, niter
    real(kind_phys)  tmp1, tmp2
    real(kind_phys)  evp_max
    real(kind_phys)  S_t_qvR, S_t_qvS
    real(kind_phys)  S_qvR, S_qvS
    real(kind_phys)  evp_qvR_pre, evp_qvS_pre
    real(kind_phys)  S_b_qvR_pre, S_b_qvS_pre
    real(kind_phys)  t_in, thv_in, tw_in, qw_in
    real(kind_phys)  t, qv, th, thv, bogbot, wd2  
    real(kind_phys)  subsat_db
    real(kind_phys)  dp, dz
    real(kind_phys)  lambda
    real(kind_phys)  rndb_evp

    ! ----------------------- !
    ! Compute basic variables !
    ! ----------------------- !

    niter = 1
    if( ievp_prep .eq. -1 .or. ievp_prep .eq. -5 ) niter = 0   ! Forward Method
    lambda   = 0.5_kind_phys 
    rndb_evp = real(ndb_evp,kind_phys)
    dp       = p_b - p_t              ! [  Pa ] >= 0.
    dz       = z_t - z_b              ! [  z  ] >= 0.
    t_in     = th_in*exnf(p_b)
    call findsp_single( qv_in, t_in, p_b, tw_in, qw_in, ii, kk, lchnk )

    ! -------------------------------------------------------------------------------------------------------------- !
    ! Two sufficient-necessary downdraft state variables used for computing evaporation rate at the base interface : !
    !   (1) subsat_db                                                                                                !
    !   (2) w_db                                                                                                     !
    ! which should continuously updated within the iteration loop.                                                   !
    ! -------------------------------------------------------------------------------------------------------------- !

    S_qvR      = S_t_qvR_in
    S_qvS      = S_t_qvS_in
    evp_qvR    = S_qvR * dp
    evp_qvS    = S_qvS * dp
  ! evp_qvR = max( 0._kind_phys, min( evp_qvR, eta2 * flxrain_bot_upeesm / max( nonzero, cmf_db ) / max( 1._kind_phys, rndb_evp ) ) )
  ! evp_qvS = max( 0._kind_phys, min( evp_qvS, eta2 * flxsnow_bot_upeesm / max( nonzero, cmf_db ) / max( 1._kind_phys, rndb_evp ) ) )
    evp_qvR = max( 0._kind_phys, min( evp_qvR, eta2 * flxrain_bot_upeesm / max( nonzero, cmfdb_evp ) ) )
    evp_qvS = max( 0._kind_phys, min( evp_qvS, eta2 * flxsnow_bot_upeesm / max( nonzero, cmfdb_evp ) ) )
    evp_max = max( qw_in - qv_in, 0._kind_phys )
    if( ( evp_qvR + evp_qvS ) .gt. evp_max ) then
        tmp1 = evp_qvR * evp_max / ( evp_qvR + evp_qvS )
        tmp2 = evp_qvS * evp_max / ( evp_qvR + evp_qvS )
        evp_qvR = tmp1
        evp_qvS = tmp2
    endif
    S_t_qvR    = evp_qvR / max( dp, nonzero )
    S_t_qvS    = evp_qvS / max( dp, nonzero )
    t          = t_in - ( xlv / cp ) * evp_qvR - ( xls / cp ) * evp_qvS
    qv         = qv_in + evp_qvR + evp_qvS 
    if( ievp_prep .eq. 1 ) then
        t      = t_in
        qv     = qv_in
    endif

    call qsat(t,p_b,es,qs)
    subsat_db  = min( 1._kind_phys, max( 0._kind_phys, 1._kind_phys - qv / max( qs, nonzero ) ) )
    th         = t / exnf(p_b)
    thv        = th * ( 1._kind_phys + zvir * qv - ql_in - qi_in )
    bogbot     = rbuoy_dn * ( 1._kind_phys - thv / thv_mean_b )
    call progup_wu2( -( rdrag*eps_dn - rjet*del_dn ), rho, p_t, p_b, -bogtop, -bogbot, w_dt**2._kind_phys, 0._kind_phys, wd2 )
    w_db       = max( wdmin, sqrt( max( wd2, nonzero ) ) )

    S_b_qvR    = ( 1._kind_phys / ( rho * g ) / max( w_db, wdmin ) ) * kevp_rain_dn * subsat_db * & 
                 sqrt( max( 0._kind_phys, flxrain_bot_upeesm / max( a_p_msfc, nonzero ) ) )
    S_b_qvS    = ( 1._kind_phys / ( rho * g ) / max( w_db, wdmin ) ) * kevp_snow_dn * subsat_db * & 
                 sqrt( max( 0._kind_phys, flxsnow_bot_upeesm / max( a_p_msfc, nonzero ) ) )

    ! ---------------------------------------------------------------------------------------------------------------- !
    ! Perform implicit iteration                                                                                       !
    ! The requires output from the iteration loop :                                                                    !
    !   (1) 'evp_qvR, evp_qvS [kg/kg]    >= 0'                                                                         ! 
    !   (2) 'S_b_qvR, S_b_qvS [kg/kg/Pa] >= 0'                                                                         !
    ! where the above (1) and (2) are fully consistent ( evp_qvR = 0.5_kind_phys * ( S_t_qvR + S_b_qvR ) * dp, evp_qvS = .. ) !
    ! regardless of the iteration number. However, consistency between 't = t_in-(xlv/cp)*evp_qvR-(xls/cp)*evp_qvS,    !
    ! qv = qv_in + evp_qvR + evp_qvS' and 'S_b_qvR, S_b_qvS' at the base interface can be obtained as iteration is     !
    ! executed many times.                                                                                             !
    ! By setting 'do iter = 1, 1', I can use centered difference instead of forward difference which is a default.     !
    ! ---------------------------------------------------------------------------------------------------------------- !

    id_exit = 0
    do iter = 1, niter

       ! ----------------------------------------------------- !
       ! Compute 'raw' evaporation rate at the base interface. !
       ! ----------------------------------------------------- !

       S_b_qvR = ( 1._kind_phys / ( rho * g ) / max( w_db, wdmin ) ) * kevp_rain_dn * subsat_db * & 
                 sqrt( max( 0._kind_phys, flxrain_bot_upeesm / max( a_p_msfc, nonzero ) ) )
       S_b_qvS = ( 1._kind_phys / ( rho * g ) / max( w_db, wdmin ) ) * kevp_snow_dn * subsat_db * & 
                 sqrt( max( 0._kind_phys, flxsnow_bot_upeesm / max( a_p_msfc, nonzero ) ) )

       ! --------------------------------------------------------------------------------------- !
       ! Impose a limiter on the computed 'raw' evaporation rate at the base interface.          !
       ! Use 'evp_max = qw_in - qv_in' and 'eta2 * flxrain_bot_upeesm...' which do not include   ! 
       ! evaporation of precipitation yet.                                                       !
       ! --------------------------------------------------------------------------------------- !  

       S_qvR   = 0.5_kind_phys * ( S_t_qvR + S_b_qvR )
       S_qvS   = 0.5_kind_phys * ( S_t_qvS + S_b_qvS )
       if( ievp_prep .eq. 1 ) then
           S_qvR   = S_b_qvR
           S_qvS   = S_b_qvS
       endif
       evp_qvR = S_qvR * dp
       evp_qvS = S_qvS * dp
     ! evp_qvR = max( 0._kind_phys, min( evp_qvR, eta2 * flxrain_bot_upeesm / max( nonzero, cmf_db ) / max( 1._kind_phys, rndb_evp ) ) )
     ! evp_qvS = max( 0._kind_phys, min( evp_qvS, eta2 * flxsnow_bot_upeesm / max( nonzero, cmf_db ) / max( 1._kind_phys, rndb_evp ) ) )
       evp_qvR = max( 0._kind_phys, min( evp_qvR, eta2 * flxrain_bot_upeesm / max( nonzero, cmfdb_evp ) ) )
       evp_qvS = max( 0._kind_phys, min( evp_qvS, eta2 * flxsnow_bot_upeesm / max( nonzero, cmfdb_evp ) ) )
       evp_max = max( qw_in - qv_in, 0._kind_phys )
       if( ( evp_qvR + evp_qvS ) .gt. evp_max ) then
           tmp1 = evp_qvR * evp_max / ( evp_qvR + evp_qvS )
           tmp2 = evp_qvS * evp_max / ( evp_qvR + evp_qvS )
           evp_qvR = tmp1
           evp_qvS = tmp2
       endif
     ! Mar.15.2014. By commenting in (activating) below 'if' block, compute 'non-zero' 'S_b_qvR,S_b_qvS'
     !              for mixing downdraft generated in the current layer. This is important to obtain reasonably
     !              strong downdraft vertical velocity. 
!?     if( ks .ne. kk .and. dp .gt. 0._kind_phys ) then
       if( ks .ne. kk ) then
       S_qvR = evp_qvR / max( dp, nonzero )
       S_qvS = evp_qvS / max( dp, nonzero )
       S_b_qvR = 2._kind_phys * S_qvR - S_t_qvR    ! IMPORTANT : This must be allowed to be negative to prevent negative precipitation flux.   
       S_b_qvS = 2._kind_phys * S_qvS - S_t_qvS    ! IMPORTANT : This must be allowed to be negative to prevent negative precipitation flux.
       if( ievp_prep .eq. 1 ) then
           S_b_qvR = 1._kind_phys * S_qvR          ! IMPORTANT : This must be allowed to be negative to prevent negative precipitation flux.   
           S_b_qvS = 1._kind_phys * S_qvS          ! IMPORTANT : This must be allowed to be negative to prevent negative precipitation flux.
       endif
       endif

       ! ------------------------------------------------------------------------------------------------------------- !
       ! Compute 'implicit' precipitation rate at the top interface by averaging the 'current' precipitation rate with !
       ! the 'previous' precipitation rate computed at the previous iteration loop.                                    !
       ! Since both 'current' and 'previous' precipitation rates satisfies the limiters,                               !
       ! the average 'implicit' precipitation rate also satisfies the limiter automatically.                           !
       ! ------------------------------------------------------------------------------------------------------------- !

       if( iter .gt. 1 ) then
           S_b_qvR = lambda * S_b_qvR + ( 1._kind_phys - lambda ) * S_b_qvR_pre
           S_b_qvS = lambda * S_b_qvS + ( 1._kind_phys - lambda ) * S_b_qvS_pre
           S_qvR   = 0.5_kind_phys * ( S_t_qvR + S_b_qvR )
           S_qvS   = 0.5_kind_phys * ( S_t_qvS + S_b_qvS )
           if( ievp_prep .eq. 1 ) then
               S_qvR   = S_b_qvR
               S_qvS   = S_b_qvS
           endif
           evp_qvR = S_qvR * dp
           evp_qvS = S_qvS * dp
         ! if( kk .ge. 10 .and. kk .le. 12 ) then
         !     write(6,*)
         !     write(6,*) 'UNICON : Convergence test within the subroutine evap_prep_up'
         !     write(6,*) 'kk, iter, abs( evp_qvR + evp_qvS - evp_qvR_pre - evp_qvS_pre )    = ', kk, iter, abs( evp_qvR + evp_qvS - evp_qvR_pre - evp_qvS_pre ) 
         !     write(6,*) 'kk, iter, evp_qvR, evp_qvS, S_b_qvR, S_b_qvS, w_db, subsat_db, dp = ', kk, iter, evp_qvR, evp_qvS, S_b_qvR, S_b_qvS, w_db, subsat_db, dp  
         !     write(6,*)
         ! endif
           if( abs( evp_qvR + evp_qvS - evp_qvR_pre - evp_qvS_pre ) .lt. 1.e-6_kind_phys ) then
               id_exit = 1 
           endif 
       endif

       S_b_qvR_pre = S_b_qvR 
       S_b_qvS_pre = S_b_qvS
       evp_qvR_pre = evp_qvR
       evp_qvS_pre = evp_qvS

       ! ----------------------------------------------------------------------------------------------------------------- !
       ! Update state variable at the top interface.                                                                       ! 
       ! At this stage, 'exql = 0.5_kind_phys * ( S_b_ql + S_t_ql ) * dp' is exactly satisfied.                                   !
       ! However, our 'S_t_ql' becomes inconsistent with the below updated 'ql = ql_in - exql' since 'S_t_ql' was computed !
       ! using 'ql = ql_in - exql(old)' where 'exql(old)' differs from 'exql'.                                             ! 
       ! This inconsistency will be removed as iteration goes on.                                                          !
       ! ----------------------------------------------------------------------------------------------------------------- !
     
       t         = t_in - ( xlv / cp ) * evp_qvR - ( xls / cp ) * evp_qvS
       qv        = qv_in + evp_qvR + evp_qvS 
       call qsat(t,p_b,es,qs)
       subsat_db = min( 1._kind_phys, max( 0._kind_phys, 1._kind_phys - qv / max( qs, nonzero ) ) )
       th        = t / exnf(p_b)
       thv       = th * ( 1._kind_phys + zvir * qv - ql_in - qi_in )
       bogbot    = rbuoy_dn * ( 1._kind_phys - thv / thv_mean_b )
       call progup_wu2( -( rdrag*eps_dn - rjet*del_dn ), rho, p_t, p_b, -bogtop, -bogbot, w_dt**2._kind_phys, 0._kind_phys, wd2 )
       w_db      = max( wdmin, sqrt( max( wd2, nonzero ) ) )
       if( id_exit .eq. 1 ) goto 10

    enddo
 10 S_b_qvR = max( 0._kind_phys, S_b_qvR ) ! Reset to non-negative value before sending to output.
    S_b_qvS = max( 0._kind_phys, S_b_qvS ) ! Reset to non-negative value before sending to output.

    ! ------------------------------- !
    ! Treatment of In-Cumulus Tracers !
    ! ------------------------------- !

    do mt = 1, ncnst
       if( mt .eq. 1 ) then
           evp_tr(mt) = evp_qvR + evp_qvS
       elseif( mt .eq. ixcldliq .or. mt .eq. ixcldice .or. mt .eq. ixnumliq .or. mt .eq. ixnumice ) then
           evp_tr(mt) = 0._kind_phys
       else
           evp_tr(mt) = flxtrrs_bot_upeesm(mt) * ( evp_qvR + evp_qvS ) / & 
                        max( ( flxrain_bot_upeesm + flxsnow_bot_upeesm ) , nonzero )
       endif
       evp_tr(mt) = max( evp_tr(mt), qmin(mt) - tr_in(mt) )
    enddo

    return

  end subroutine evap_prep_dn

!--------------------------------------------------------------------------------------------------

subroutine progup_thlqt(eps_mix,eps_dia,qsrcg,pb,pt,qmid,gamq,qub,qut)

   ! -------------------------------------------------------------------------------------------------------------------------- !
   ! Compute cumulus updraft properties at the top interface for 'thl,qt'                                                       !
   ! for individual updraft plume. This subroutine is directly from the                                                         !
   ! progupeff_thlqt but edge computation is removed for multiple plume                                                         !
   ! approach.                                                                                                                  !
   !                                                                                                                            !
   !  (1) qsrcg = ['q'/Pa]                                                                                                      !
   !                                                                                                                            !
   !      Both in case that 'q' decreases as   updraft rises (i.e., precipitation production  process with q = ql, qi)          ! 
   !      and 'q' increases as downdraft sinks (i.e., precipitation evaporation process with q = qv),                           !
   !      the specified 'qsrcg' should be positive.                                                                             !
   !      That is, we should specify 'qsrcg' in pressure-coordinate system by considering dp = -dz.                             !
   !                                                                                                                            !
   !  (2) eps_mix, eps_dia = [1/Pa]                                                                                             !
   !                                                                                                                            !
   !      Here, 'eps_mix' is an inverse of damping length scale by 'entrainment mixing' while                                   !
   !      the 'eps_dia' is an inverse of damping length scale by 'diabatic forcing'                                             ! 
   !      ('eps_dia = +  cat * max( (qc_cu - qcrit) / qc_cu, 0._kind_phys ) [1/Pa]' which is auto-conversion efficiency for updraft or ! 
   !      ('eps_dia = - (kevp / (rho*g*wd)) * sqrt( Fp / ap ) / qsd' [1/Pa]' which is evaporation efficiency for downdraft ).   !
   !      Be careful that for convective updraft, both the specified 'eps_mix, eps_dia > 0' while                               !
   !      for downdraft process, both 'eps_mix, eps_dia < 0' as is shown in the calling routine in the main body.               !
   !                                                                                                                            !
   ! Done.                                                                                                                      !                                         
   ! -------------------------------------------------------------------------------------------------------------------------- ! 

   real(kind_phys), intent(in)                     :: eps_mix, eps_dia, qsrcg
   real(kind_phys), intent(in)                     :: pb, pt
   real(kind_phys), intent(in)                     :: qub
   real(kind_phys), intent(in)                     :: qmid, gamq
   real(kind_phys), intent(out)                    :: qut
   real(kind_phys)  eps, qb, a, b, c, fp, dp
   
   ! Environmental conservative scalar at the base interface 

   dp = pb - pt
   qb = qmid + 0.5_kind_phys * dp * gamq 

   ! Ensemble-mean updraft value at the top interface
   ! Aug.18.2010. I chekced that below 'qut' exactly reproduce the UW 
   ! in a raw form, but not the Taylor-extended form. I verified that
   ! my Taylor extended form is wrong but the UW has the correct 
   ! Taylor-expanded form. Thus, I am using the UW-Taylor extended form.
   ! As a result, my nelow formula produced the exactly same results as UW in
   ! all cases.
   
   eps =   eps_mix + eps_dia
   a   = - eps    
   b   = - eps * gamq
   c   = - eps * ( qb - gamq * pb )

   if(abs(a*dp).gt.1.e-3_kind_phys) then
      ! I checked that below exactly reproduced the UW.
      !prp   qut = ( qub - ( b * pb + c ) / a + ( b / a / a ) ) * exp( a * dp ) + ( b * pt + c ) / a - ( b / a / a ) 
      fp  = 1._kind_phys - eps_dia / eps
      qut = ( qub - ( ( b * pb + c ) / a - ( b / a / a ) ) * fp - qsrcg / a ) * exp( a * dp ) + &
         ( ( b * pt + c ) / a - ( b / a / a ) ) * fp + qsrcg / a
   else   
      ! Aug.18.2010. Below original form is wrong. Thus, I am using the
      ! UW Taylor extended form.
      ! Jul.12.2011. I may need to correct below Taylor-expansion formula using both 'q_b' and 'q_t'
      !              similar to the modification in progup_wu2. 
      ! Apr.18.2012. I realized that 'dp * a * ( qub - qb )' produces results non-trivially different from
      !              the 'dp * ( a * ( qub - qb ) )'. This was the reason why '057a' produces different
      !              results from '057' in bluefire. 
      !              To be consistent with the treatment in progup_uv, I should use the UW.Correct.1.
      !              But in order to reproduce the '056' simulation, let's use UW.Correct.2. for the
      !              time being.  This should be carefully chosen in future to be consistent with the
      !              progup_uv.
      !              I checked that below UW.Correct.2. exactly reproduces '056'.
      ! qut = qub - dp * ( 0.5_kind_phys * b * ( pt + pb ) + c )            ! Original, Wrong. 
      !prp   qut = qub + dp * ( a * ( qub - qb ) - qsrcg )                ! New = UW. Correct.1.
      qut = qub + dp * ( a * ( qub - qb ) - qsrcg - eps_dia * qb ) ! New = UW. Correct.1.
      ! qut = qub + dp *   a * ( qub - qb ) - dp * qsrcg             ! New = UW. Correct.2.
   endif

   ! Below is from UW shallow convection
   ! Produced similar result as above.
   ! On Jul.21.2010., change to the UW formula for consistency with
   ! the description paper. Also change 1.e-4 to 1.e-3 for Taylor expansion.   
   
   ! if( eps*dp .gt. 1.e-3_kind_phys ) then
   !     qut = ( qmid + gamq / eps - gamq * dp / 2._kind_phys ) - &
   !           ( qmid + gamq * dp / 2._kind_phys - qub + gamq / eps ) * exp( -eps * dp )
   ! else
   !     qut = qub + ( qmid + gamq * dp / 2._kind_phys - qub ) * eps * dp
   ! endif

end subroutine progup_thlqt

!--------------------------------------------------------------------------------------------------

subroutine progup_uv(eps,PGFuv,pb,pt,qmid,gamq,gamqPGF,qub,qut)

   ! ---------------------------------------------------------------------- !
   ! Compute cumulus updraft properties at the top interface for 'u,v'      !
   ! for individual updraft plume. This subroutine is directly from the     !
   ! progupeff_thlqt but edge computation is removed for multiple plume     !
   ! approach.                                                              !
   ! This is same as 'progupeff_thlqt_single' except that coef. c is        !
   ! re-defined by including PGFc effect.                                   !
   ! Done.                                                                  !
   ! ---------------------------------------------------------------------- ! 

   real(kind_phys), intent(in)                     :: eps, PGFuv
   real(kind_phys), intent(in)                     :: pb, pt
   real(kind_phys), intent(in)                     :: qub
   real(kind_phys), intent(in)                     :: qmid, gamq, gamqPGF 
   real(kind_phys), intent(out)                    :: qut
   real(kind_phys)  qb, a, b, c, dp

   ! Environmental conservative scalar at the base interface 
   ! Aug.18.2010. I checked that below 'qut' exactly reproduce the UW 
   ! in a raw form, but not the Taylor-extended form. I verified that
   ! my Taylor extended form is wrong but the UW has the correct 
   ! Taylor-expanded form. Thus, I am using the UW-Taylor extended form.
   ! As a result, my nelow formula produced the exactly same results as UW in
   ! all cases.

   ! Apr.5.2011. The 'gamq' multiplied to 'PGFuv' may need to be the one 
   ! using the inter-layer slope, not the within the layer slope. Thus,
   ! I am using gamqPGF for PGFc effect.

   dp = pb - pt
   qb = qmid + 0.5_kind_phys * dp * gamq 

   ! Ensemble-mean updraft value at the top interface

   a = - eps  
   b = - eps * gamq
   c = - eps * ( qb - gamq * pb ) + PGFuv * gamqPGF

   if( abs( a * dp ) .gt. 1.e-3_kind_phys ) then
      qut = ( qub - ( b * pb + c ) / a + ( b / a / a ) ) * exp( a * dp ) + ( b * pt + c ) / a - ( b / a / a ) 
   else   
      ! Aug.18.2010. Below original form is wrong. Thus, I am using the
      ! UW Taylor extended form.
      ! Apr.18.2012. For consistency with progup_thlqt, I use UW.Correct.2. below.
      !              The diffrence between UW.Correct.1. and UW.Correct.2. is from numerical truncation error.
      ! qut = qub - dp * ( 0.5_kind_phys * b * ( pt + pb ) + c )          ! Original, Wrong.
      qut = qub + dp * ( a * ( qub - qb ) - PGFuv * gamqPGF )    ! New = UW. Correct.1.
      ! qut = qub + dp *   a * ( qub - qb ) - dp * gamqPGF * PGFuv ! New = UW. Correct.2.
   endif

   ! Below is from UW shallow convection
   ! Produced similar result as above.
   ! On Jul.21.2010., change to the UW formula for consistency with
   ! the description paper. Also change 1.e-4 to 1.e-3 for Taylor expansion.   
    
   ! if( eps * dp .gt. 1.e-3_kind_phys ) then
   !     qut = ( qmid + ( 1._kind_phys - PGFuv ) * gamq / eps - gamq * dp / 2._kind_phys ) - &
   !           ( qmid + gamq * dp / 2._kind_phys - qub + ( 1._kind_phys - PGFuv ) * gamq / eps ) * exp( -eps * dp )
   ! else
   !     qut = qub + ( qmid + gamq * dp / 2._kind_phys - qub ) * eps * dp - PGFuv * gamq * dp
   ! endif

end subroutine progup_uv

!--------------------------------------------------------------------------------------------------

subroutine progup_wu2(eps,rho,pb,pt,bogbot,bogtop,wwub,wwe,wwut)

   ! ----------------------------------------------------------------- !
   ! Compute squared ensemble-mean cumulus updraft vertical at the top ! 
   ! interface. Note that 'bogbot,bogtop' ( which are non-dimensional  !
   ! variables ) already contains buoyancy coefficient, rbuoy in it.   ! 
   ! Eqn: dw2/dp + a*w2 = c                                            !
   ! Note that even the case of 'wwub=0.' is reasonably treated within !
   ! this subroutine.                                                  ! 
   ! Note that following 'compute_dp' should be chosen in a consistent !
   ! way as the one in this subroutine.                                !
   ! Mar.11.2013. Add 'wwe [m2/s2]' which is the square of vertical    !
   !              velocity of mixing environmental airs that is        !
   !              assumed to be a height-independent constant in the   !
   !              layer considered.                                    !   
   ! worg, bogbot, bogtop : Have no unit.                              !
   ! Done.                                                             !
   ! ----------------------------------------------------------------- !

   real(kind_phys), intent(in)                     :: eps, rho
   real(kind_phys), intent(in)                     :: pb, pt 
   real(kind_phys), intent(in)                     :: bogbot, bogtop
   real(kind_phys), intent(in)                     :: wwub
   real(kind_phys), intent(in)                     :: wwe
   real(kind_phys), intent(out)                    :: wwut
   !real(kind_phys)  a, c
   real(kind_phys)  dp, expfac, gammaB, worg
   !real(kind_phys)  delbog

   ! Coefficients

   dp = pb - pt
   !a  = - 2._kind_phys * eps  
   !c  = - ( bogbot + bogtop ) / rho
   worg = rho * eps * wwe

   ! Option 1. Ensemble-mean updraft value at the top interface without
   ! considering buoyancy slope in each layer.
   ! For consistency with compute_dp, good to use this. 
   ! Aug.18.2010. I chekced that below 'qut' exactly reproduce the UW 
   ! if 'bogbot = bogtop' in the UW form.
   
   ! if( abs( a * dp ) .gt. 1.e-3_kind_phys ) then
   !    wwut = ( wwub - c / a ) * exp( a * dp ) + c / a 
   ! else   
   !  ! Aug.18.2010. Below original form is wrong. Thus, I am using the
   !  ! correct UW Taylor extended form.
   !  ! wwut = wwub - c * dp                                ! Original Wrong
   !    wwut = wwub * ( 1._kind_phys - 2._kind_phys * eps * dp ) - c * dp ! New Correct.
   ! endif

   ! Option 2. Ensemble-mean updraft value at the top interface
   !           with consideration of buoyancy slope in each layer.
   !           Below is the most perfect formula I should use.
   !           My below formula produced the exactly same results as UW in
   !           all cases.
   ! Jul.9.2011. Below is not correct for downdraft momentum equation.
   !             I should correct this. Probably, I should correct ( add (-) sign )
   !             to the input argument of bogbot and bogtop.
   ! Jul.10.2011. I modified below Taylor expansion formula by using 'eps .gt. 1.e-5' and
   !              by changing 'bogbot' to '0.5_kind_phys * ( bogbot + bogtop )' for full consistency.
   ! Below is old before Jul.10.2011.
   ! if( eps * dp .gt. 1.e-3_kind_phys ) then
   !     expfac = exp( -2._kind_phys * eps * dp )
   !     gammaB = ( bogtop - bogbot ) / dp
   !     wwut = wwub * expfac + ( gammaB * dp + (1._kind_phys-expfac)*(bogbot + gammaB/(-2._kind_phys*eps)) )/(rho*eps)
   !  else
   !     wwut = wwub * ( 1._kind_phys - 2._kind_phys * eps * dp ) + 2._kind_phys * bogbot * dp / rho 
   ! endif
   ! Below is new after Jul.10.2011.
   ! Mar.11.2013. Add 'worg' only in the below block of formula.
   !              Thus, I should only use below block. 

   if( abs(dp) .lt. 1.e-10_kind_phys ) then
      wwut = wwub
   else
      if( abs(eps) .gt. 1.e-6_kind_phys ) then
         expfac = exp( -2._kind_phys * eps * dp )
         gammaB = ( bogtop - bogbot ) / dp
         wwut = wwub * expfac + ( gammaB * dp + (1._kind_phys-expfac)*(bogbot + gammaB/(-2._kind_phys*eps) + worg) )/(rho*eps)
      else
         wwut = wwub * ( 1._kind_phys - 2._kind_phys * eps * dp ) + ( bogbot + bogtop + 2._kind_phys * worg ) * dp / rho 
      endif
   endif

   ! Below is from UW shallow convection
   ! Produced similar result as above.
   ! On Jul.21.2010., change to the UW formula for consistency with the description paper.
   ! However, Taylor-expanded formula is correctly revised.   

   ! delbog = bogtop - bogbot
   ! expfac = exp( -2._kind_phys * eps * dp )
   ! if( eps * dp .gt. 1.e-3_kind_phys ) then
   !     wwut = wwub * expfac + ( delbog + (1._kind_phys-expfac)*(bogbot + delbog/(-2._kind_phys*eps*dp)))/(rho*eps)
   ! else
   !   ! Below is original Taylor formula.
   !   ! wwut = wwub + dp * ( bogbot + bogtop ) / rho
   !   ! Below is corrected Taylor formula on Jul.21.2010.
   !     wwut = wwub * ( 1._kind_phys - 2._kind_phys * eps * dp ) + 2._kind_phys * bogbot * dp / rho 
   ! endif

end subroutine progup_wu2

!--------------------------------------------------------------------------------------------------

real(kind_phys) function compute_dp(eps,rho,pb,pt,bogbot,bogtop,wwub,wwe)

   ! ---------------------------------------------------------------------- !
   ! Compute vertical distance that convective updraft with non-zero        ! 
   ! entrainment can move.                                                  !
   ! For simplicity, use the height-independent average buoyancy.           !
   ! Note 'compute_dp = [Pa] > 0'. eps = [1/Pa]. bogbot,bgtop = [ no unit ] !
   ! Note that above 'progup_wu2' should be chosen in a consistent way      !
   ! as the one in this subroutine.                                         !
   ! Mar.11.2013. Add 'wwe [m2/s2]' which is the square of vertical         !
   !              velocity of mixing environmental airs that is             !
   !              assumed to be a height-independent constant in the        !
   !              layer considered.                                         !   
   ! worg, bogbot, bogtop : Have no unit.                                   !
   ! Done.                                                                  !
   ! ---------------------------------------------------------------------- !

   real(kind_phys), intent(in)                     :: eps, rho
   real(kind_phys), intent(in)                     :: pb, pt 
   real(kind_phys), intent(in)                     :: bogbot, bogtop
   real(kind_phys), intent(in)                     :: wwub
   real(kind_phys), intent(in)                     :: wwe
   !real(kind_phys)  a, c
   real(kind_phys)  dp, dpen, gammaB, worg
   real(kind_phys)  x0, x1, f, fs, s00
   integer   iteration

   !a = 2._kind_phys * eps  
   !c = ( bogbot + bogtop ) / rho
   dpen = pb - pt
   gammaB = ( bogtop - bogbot ) / dpen
   worg = rho * eps * wwe

   ! Option 1. Without considering buoyancy slope in each layer.  

   ! if( c .ge. 0._kind_phys ) then
   !     dp = dpen 
   ! else
   !     if( a .lt. 1.e-12_kind_phys ) then
   !         dp = - wwub / c
   !     else
   !         dp = - ( 1._kind_phys / a ) * log( c / ( c - a * wwub ) )
   !     endif
   ! endif 
   
   ! Option 2. With considering buoyancy slope in each layer.  
   !           Below is from the UW code.
   !           Below is the most perfect choice.
   
   ! Jul.10.2011. Rigorously speaking, below Taylor expansion should also be modified following
   !              the same day's modification in progup_wu2 above, using 'eps' not the 
   !              product of 'eps * dpen', although this modification is likely to have
   !              a very minor effect. 
   ! Mar.11.2013. Add 'worg' in the below block. However, I have not added 'worg' in the
   !              Taylor-expended block, which should be done in future. 
   !              Probably, I should redefine s00 = ( bogbot + worg ) / rho - eps * wwub.
   !              Since this seems to be quite straightforwatd, I also modified this Taylor block. 

   ! s00 = bogbot / rho - eps * wwub
   s00 = ( bogbot + worg ) / rho - eps * wwub
   if( eps * dpen .le. 1.e-8_kind_phys ) then
      if( s00 .ge. 0._kind_phys ) then
         x0 = dpen
      else
         x0 = max( 0._kind_phys, min( dpen , -0.5_kind_phys * wwub / s00 ) )
      endif
   else
      if( s00 .ge. 0._kind_phys ) then
         x0 = dpen
      else
         x0 = 0._kind_phys
      endif
      do iteration = 1, 5
         f  = exp(-2._kind_phys*eps*x0)*(wwub-(bogbot-gammaB/(2._kind_phys*eps)+worg)/(eps*rho)) + &
            (gammaB*x0+bogbot-gammaB/(2._kind_phys*eps)+worg)/(eps*rho)
         fs = -2._kind_phys*eps*exp(-2._kind_phys*eps*x0)*(wwub-(bogbot-gammaB/(2._kind_phys*eps)+worg)/(eps*rho)) + &
            (gammaB)/(eps*rho)
         ! Sep.28.2010. Very rarely, fs = 0 happens. So, I added below fixer.
         if( fs .ge. 0._kind_phys ) fs = max(fs, nonzero)
         if( fs .lt. 0._kind_phys ) fs = min(fs,-nonzero)
         x1 = x0 - f/fs
         x0 = x1
         x0 = min( dpen, max( 0._kind_phys, x0 ) )
      end do
   endif
   dp = x0

   compute_dp = min( dpen, max( 0._kind_phys, dp ) )

end function compute_dp

!--------------------------------------------------------------------------------------------------

subroutine buosort_downdraft( cuL, cuU, enL, xdown_min, xdown_max )

   ! -------------------------------------------------------------- !
   ! Perform buoyancy sorting of downdrafts and find the ranges of  ! 
   ! mixing fraction, (xdown_min, xdown_max), where mixtures within !
   ! this range will move down into the base interface.             !
   ! It must be and is that x = 0 is corresponding to cuL, while    !
   !                        x = 1 is corresponding to cuU.          !
   ! Done.                                                          !
   ! -------------------------------------------------------------- !

   real(kind_phys), intent(in)   :: cuL, cuU, enL 
   real(kind_phys), intent(out)  :: xdown_min, xdown_max

   if( cuU.gt.cuL ) then
      xdown_min = 0._kind_phys
      xdown_max = min(1._kind_phys,max(0._kind_phys,(enL-cuL)/(cuU-cuL)))
      return
   elseif( cuU.lt.cuL ) then
      xdown_min = min(1._kind_phys,max(0._kind_phys,(enL-cuL)/(cuU-cuL)))
      xdown_max = 1._kind_phys
      return
   elseif( cuU.eq.cuL ) then
      if( cuU.lt.enL ) then
         xdown_min = 0._kind_phys
         xdown_max = 1._kind_phys
         return
      elseif( cuU.gt.enL ) then
         xdown_min = 0._kind_phys
         xdown_max = 0._kind_phys
         return
      else
         ! Below can be any values between 0 <= x <= 1
         ! Apr.21.2011. xdown_min is changed to 0 instead of 0.5 as described in the text.
         ! xdown_min = 0.5_kind_phys
         xdown_min = 0.0_kind_phys
         xdown_max = 0.5_kind_phys
         return
      endif
      return
   endif

end subroutine buosort_downdraft

!--------------------------------------------------------------------------------------------------

subroutine compute_PDF( PDFtype, zLL, zUU, zbar, zmass, zmass_L )

   ! ---------------------------------------------------------------------------------------- !
   ! Compute mass-flux weighted (or equally, probability weighted) mean z                     !
   ! (zbar) and corrresponding mass fraction (zmass) within a given single                    !
   ! interval between (zL,zU).                                                                !
   ! The PDFtype is as follows:                                                               ! 
   !    'PDFupG' : Updraft PDF with Gamma distribution   ( 0 < z = alpha < alpha_max = 3 )    !
   !    'PDFupW' : Updraft PDF with Weibull distri.      ( 0 < z = alpha < alpha_max = 3 )    !
   !    'PDFdnU' : Downdraft PDF with Uniform distri.    ( 0 < z = alpha < 2 )                !     
   !    'PDFdnT' : Downdraft PDF with Triangular distri. ( 0 < z = alpha < 2 )                ! 
   !    'PDFbsU' : Microscopic buoyancy sorting PDF with Uniform distri.    ( 0 < z = x < 1 ) !
   !    'PDFbsT' : Microscopic buoyancy sorting PDF with Triangular distri. ( 0 < z = x < 1 ) !
   !    'PDFbsQ' : Microscopic buoyancy sorting PDF with Quadratic distri.  ( 0 < z = x < 1 ) !
   !    'PDFbsN' : Microscopic buoyancy sorting PDF with p=0.5 in symmetric beta distri.      !
   !               ( 0 < z = x < 1 )                                                          !
   ! Note that                                                                                ! 
   !     x = 0 : cumulus core,                                                                !
   !     x = 1 : environmental mean                                                           !
   !     alpha = 0 : cumulus edge toward environment                                          !
   !     alpha = 2 or 3 : cumulus core-edge away from the environment                         !
   ! As the most general application, I should use the following                              !
   ! symmetric beta distribution for microscopic buoyancy sorting                             !
   !        P(x) = (x*(1-x))**(p-1) / B(p,p)   , p > 0, 0 <= x <= 1                           !          
   !        B(p,p) = Gamma(p)*Gamma(p)/Gamma(2*p)                                             !
   !        B(0.1,0.1) = 19.7146, B(0.5,0.5) = 3.1416                                         !
   !        B(1,1) = 1, B(2,2) = 1/6                                                          !
   ! This shows how much time the system has for mixing for a given                           !
   ! model time step dt. If 'p' increases, it means that the system has                       !
   ! more time for buoyancy sorting mixing. This 'p' can parameterized                        !  
   ! as a function of timeres/dt.                                                             !          
   ! Nov.04.2014.                                                                             !
   ! Currently, only 'PDFbsU' and 'PDFbsQ' computes 'zmass_L'. So, only choose these two      !
   ! options. Computation of 'zmass_L' for the other cases should be done later.              !            
   ! Done.                                                                                    !
   ! ---------------------------------------------------------------------------------------- !

   character(len=6), intent(in)   :: PDFtype 
   real(kind_phys),  intent(in)   :: zLL, zUU
   real(kind_phys),  intent(out)  :: zbar, zmass, zmass_L
   real(kind_phys)   zmax, Pn, zL, zU

   ! Take such that zU is always equal or larger than zL.

   zL = min(zLL,zUU) 
   zU = max(zLL,zUU) 
  
   ! Main computation.
   ! 'zmass' is fractional area (mass) surrounded by [zL,zU]
   ! 'zbar'  is PDF-weighted mean z in [zL,zU]
   ! 'zmass_L' is fractional mass of the convective updraft (on the zL side) 

   select case (PDFtype)
   case ('PDFupG')
      zmax = alpha_max
      Pn   = 1._kind_phys - exp(-2._kind_phys*zmax)*(1._kind_phys+2._kind_phys*zmax)
      ! Currently, below analytical formula is obtained only when gma = 2._kind_phys
      zmass = (exp(-2._kind_phys*zL)*(1._kind_phys+2._kind_phys*zL)-exp(-2._kind_phys*zU)*(1._kind_phys+2._kind_phys*zU))/Pn
      if( abs(zL-zU) .le. 1.e-10_kind_phys ) then
         zbar = zL
         zmass = 0._kind_phys
         return
      else
         zbar = exp(-2._kind_phys*zL)*(zL**2+zL+0.5_kind_phys)-exp(-2._kind_phys*zU)*(zU**2+zU+0.5_kind_phys)
         zbar = 2._kind_phys*zbar/Pn/zmass
         return
      endif
   case ('PDFdnU')
      zmax = 2._kind_phys
      zmass = 0.5_kind_phys*(zU-zL)
      zbar  = 0.5_kind_phys*(zU+zL)
      return
   case ('PDFbsU')
      zmax = 1._kind_phys
      zmass = 1.0_kind_phys*(zU-zL)
      zbar  = 0.5_kind_phys*(zU+zL)
      zmass_L = zU-0.5_kind_phys*zU**2-(zL-0.5_kind_phys*zL**2)
      return
   case ('PDFbsQ')
      zmax = 1._kind_phys
      zmass = (3._kind_phys*zU**2-2._kind_phys*zU**3)-(3._kind_phys*zL**2-2._kind_phys*zL**3)
      if( zmass .gt. 0._kind_phys ) then
         zbar = 0.5_kind_phys*((4._kind_phys*zU**3-3._kind_phys*zU**4)-(4._kind_phys*zL**3-3._kind_phys*zL**4))/zmass
      else
         zbar = 0.5_kind_phys*(zL+zU)
      endif
      zmass_L = 3._kind_phys*zU**2-4._kind_phys*zU**3+1.5_kind_phys*zU**4-(3._kind_phys*zL**2-4._kind_phys*zL**3+1.5_kind_phys*zL**4)
      return
   case ('PDFbsN')
      zmax = 1._kind_phys
      zmass = 0.3183_kind_phys*( asin(2._kind_phys*zU-1._kind_phys) - asin(2._kind_phys*zL-1._kind_phys) )
      if( zL .eq. 0._kind_phys .and. zU .eq. 1._kind_phys ) zmass = 1._kind_phys
      zmass = max(0._kind_phys,min(1._kind_phys,zmass)) 
      if( zmass .gt. 0._kind_phys ) then
         zbar = ( 0.5_kind_phys*asin(2._kind_phys*zU-1._kind_phys) - sqrt(zU*(1._kind_phys-zU)) - 0.5_kind_phys*asin(2._kind_phys*zL-1._kind_phys) + sqrt(zL*(1._kind_phys-zL)) ) / &
            ( asin(2._kind_phys*zU-1._kind_phys) - asin(2._kind_phys*zL-1._kind_phys) )
      else
         zbar = 0.5_kind_phys*(zL+zU)
      endif
      return

   end select

end subroutine compute_PDF

!--------------------------------------------------------------------------------------------------

subroutine compute_epsdelnod( PDFtype, xc, epsnod, delnod )

   ! ------------------------------------------------------------------ !
   ! Compute non-dimensional fraction entrainment and detrainment rate  !
   ! for a given 'xc' and normalized symmetric PDF(x) within 0 < x < 1. !
   ! This is for microscopic buoyancy sorting.                          ! 
   ! The PDFtype is as follows:                                         !
   !  'PDFbsU' : PDF with Uniform    distri. ( 0 < x < 1 )              !
   !  'PDFbsT' : PDF with Triangular distri. ( 0 < x < 1 )              !
   !  'PDFbsQ' : PDF with Quadratic  distri. ( 0 < x < 1 )              !
   !  'PDFbsN' : PDF with p=0.5 in symmetric beta distri. ( 0 < x < 1 ) !
   ! Note that                                                          ! 
   !     x = 0 : cumulus core,                                          !
   !     x = 1 : environmental mean                                     !
   ! In order to avoid infinite area at x = 1, either triangular or     !
   ! quadratic distribution is better than the uniform distribution.    !
   ! As the most general application, I should use the following        !
   ! symmetric beta distribution:                                       !
   !        P(x) = (x*(1-x))**(p-1) / B(p,p)   , p > 0, 0 <= x <= 1     !          
   !        B(p,p) = Gamma(p)*Gamma(p)/Gamma(2*p)                       !
   !        B(0.1,0.1) = 19.7146, B(0.5,0.5) = 3.1416                   !
   !        B(1,1) = 1, B(2,2) = 1/6                                    !
   ! This shows how much time the system has for mixing for a given     !
   ! model time step dt. If 'p' increases, it means that the system has !
   ! more time for buoyancy sorting mixing. This 'p' can parameterized  ! 
   ! as a function of timeres/dt.                                       !          
   ! Done.                                                              !
   ! ------------------------------------------------------------------ !

   character(len=6), intent(in)  :: PDFtype 
   real(kind_phys), intent(in)   :: xc
   real(kind_phys), intent(out)  :: epsnod, delnod
  
   ! Main computation.
   ! 'epsnod' is non-dimensional fractional entrainment rate
   ! 'delnod' is non-dimensional fractional detrainment rate

   select case (PDFtype)
   case ('PDFbsU')
      epsnod = xc**2
      delnod = (1._kind_phys-xc)**2 
      return
   case ('PDFbsT')
      if( xc .le. 0.5_kind_phys ) then
         epsnod = (8._kind_phys/3._kind_phys)*xc**3 
         delnod =  1._kind_phys-4._kind_phys*xc**2+(8._kind_phys/3._kind_phys)*xc**3
         return
      else
         epsnod = -1._kind_phys/3._kind_phys + 4._kind_phys*xc**2-(8._kind_phys/3._kind_phys)*xc**3 
         delnod = (8._kind_phys/3._kind_phys)*(1._kind_phys-xc)**3
         return
      endif
      return
   case ('PDFbsQ')
      epsnod = xc**3*(4._kind_phys-3._kind_phys*xc)
      delnod = (1._kind_phys-6._kind_phys*xc**2+8._kind_phys*xc**3-3._kind_phys*xc**4)
      return
   case ('PDFbsN')
      epsnod =         2._kind_phys*0.3183_kind_phys*(0.5_kind_phys*asin(2._kind_phys*xc-1._kind_phys)-sqrt(xc*(1._kind_phys-xc))-0.5_kind_phys*asin(-1._kind_phys)) 
      delnod = 1._kind_phys - 2._kind_phys*0.3183_kind_phys*(0.5_kind_phys*asin(2._kind_phys*xc-1._kind_phys)+sqrt(xc*(1._kind_phys-xc))-0.5_kind_phys*asin(-1._kind_phys))
      if( xc .eq. 0._kind_phys ) then
         epsnod = 0._kind_phys
         delnod = 1._kind_phys
      elseif( xc .eq. 1._kind_phys ) then                                            
         epsnod = 1._kind_phys
         delnod = 0._kind_phys
      endif
      epsnod = max(0._kind_phys,min(1._kind_phys,epsnod))
      delnod = max(0._kind_phys,min(1._kind_phys,delnod)) 
      return
   end select

end subroutine compute_epsdelnod

!--------------------------------------------------------------------------------------------------

subroutine buosorts_UW(rbuoy,p,w_cu,thl_cu,qt_cu,w_eg,thl_eg,qt_eg,thv_en,cridis,xc,xs,thv_cu,thv_eg,thvxs)

   ! ---------------------------------------------------------------------- !
   ! Buoyancy Sorting Algorithm from UWShCu but can treat the case of       !
   ! different 'thv_eg' and 'thv_en' for handling organized non-uniform     ! 
   ! environmental airs being entrained.                                    !
   ! This codes assumes that gamv_en = gamv_cu = 0                          !
   ! Except this, this code is quite general and correct for any 'cridis',  !
   ! and computationally efficient. So, I should use this code in my GCM.   !
   ! For microscopic buoyancy sorting,                                      !
   !     x = 0 : cumulsu ensemble-mean                                      !
   !     x = 1 : non-uniform environmental value                            !
   ! Done.                                                                  !
   ! May.03.2011. In order to save computation time, I changed the variable !
   !              names: thvxsat -> thvxs, thlxsat -> thlxs, qtxsat -> qtxs !
   !              xsat -> xs which does not change the answer.              !
   ! Mar.11.2013. I added 'w_eg' to treat the effect of organized flow      !
   !              within PBL.                                               !
   ! ---------------------------------------------------------------------- ! 

   real(kind_phys), intent(in)          :: rbuoy    
   real(kind_phys), intent(in)          :: w_cu, w_eg
   real(kind_phys), intent(in)          :: thl_cu, thl_eg
   real(kind_phys), intent(in)          :: qt_cu, qt_eg
   real(kind_phys), intent(in)          :: cridis
   real(kind_phys), intent(in)          :: p, thv_en
   real(kind_phys), intent(out)         :: xc, xs, thv_cu, thv_eg, thvxs
   integer   id_check, status, kk
   real(kind_phys)  thlxs, qtxs, x_cu, x_en, thv_x0, thv_x1
   real(kind_phys)  th, qv, ql, qi, qse
   real(kind_phys)  aquad, bquad, cquad, xc1, xc2, excess_cu, excess_eg, xs1, xs2
   real(kind_phys)  es
   real(kind_phys)  qs
   real(kind_phys)  qsat_arg
   real(kind_phys)  exn

   ! ---------------------------------------------------------------- !
   ! Calculate environmental and cumulus saturation.                  !
   ! Note that in order to calculate saturation excess, we should use !
   ! liquid water temperature instead of temperature  as the argument !
   ! of "qsat". But note normal argument of "qsat" is temperature.    !
   ! ---------------------------------------------------------------- !

   exn        = (p/p00)**rovcp
   call conden(p,thl_eg,qt_eg,th,qv,ql,qi,qse,id_check)
   thv_eg     = th * ( 1._kind_phys + zvir*qv - ql - qi )
   qsat_arg   = thl_eg*exn
   call qsat(qsat_arg, p, es, qs)
   excess_eg  = qt_eg - qs

   call conden(p,thl_cu,qt_cu,th,qv,ql,qi,qse,id_check)
   thv_cu     = th * ( 1._kind_phys + zvir * qv - ql - qi )
   qsat_arg   = thl_cu*exn
   call qsat(qsat_arg, p, es, qs)
   excess_cu  = qt_cu - qs

   if( (excess_cu*excess_eg).lt.0._kind_phys ) then
      ! May.03.2011. Regardless of the relative magnitude of 'excess_cu,excess_eg', it should be
      !              xs = excess_cu / ( excess_cu - excess_eg ). In the mother code, this 'xs'
      !              is also directly used to compute 'xe_min, xe_max'. So, below mistake can 
      !              influence actual model computation. I fixed this bug by commenting out
      !              4 lines in the below if block.  
      ! if( excess_cu .gt. excess_eg ) then 
      xs = excess_cu / ( excess_cu - excess_eg );
      ! else
      !     xs = excess_eg / ( excess_eg - excess_cu );
      ! endif
      ! May.03.2011. In order to be fully compatible with the mother routine that computes
      !              'xe_min, xe_max', I should set 'xs = 0' for all these cases. This also
      !              removes the distinuity when excess_cu=excess_eg=0. Thus, I removed
      !              below 5 lines and reset xs = 0.
      ! else
      !     xs = 0._kind_phys
      ! endif
      ! Below block is old code before May.03.2011.  
      ! May.16.2011. I restored to the 'old code, since 'xs = 1._kind_phys' is also fully compatible with the
      !              mother subroutine and correct even within this subroutine.
   elseif( excess_cu.le.0._kind_phys .and. excess_eg.le.0._kind_phys ) then
      xs = 0._kind_phys
   elseif( excess_cu.ge.0._kind_phys .and. excess_eg.ge.0._kind_phys ) then
      xs = 1._kind_phys
   endif

   ! ----------------------------------------------------------------- !
   ! Case 1 : When both cumulus and env. are unsaturated or saturated. !
   ! ----------------------------------------------------------------- !
  
   thvxs = thv_cu
   ! May.03.2011. In order to save computation time, I re-write below if line using the above newly computed 'xs',
   !              which should produce identical results. 
   if( xs .eq. 0._kind_phys .or. xs .eq. 1._kind_phys ) then
      ! if( ( excess_cu .le. 0._kind_phys .and. excess_eg .le. 0._kind_phys ) .or. ( excess_cu .ge. 0._kind_phys .and. excess_eg .ge. 0._kind_phys ) ) then
      ! Below is the original UW code assuming 'thv_eg=thv_en'
      ! xc = min(1._kind_phys,max(0._kind_phys,1._kind_phys-2._kind_phys*rbuoy*g*cridis/w_cu**2._kind_phys*(1._kind_phys-thv_cu/thv_en)))
      ! Below is the revised code considering the difference between 'thv_eg' and 'thv_en'
      thv_x0 = thv_cu;
      thv_x1 = thv_eg;
      ! aquad =  w_cu**2;
      ! bquad =  2._kind_phys*rbuoy*g*cridis*(thv_x1 - thv_x0)/thv_en - 2._kind_phys*w_cu**2;
      aquad =  (w_cu-w_eg)**2;
      bquad =  2._kind_phys*rbuoy*g*cridis*(thv_x1 - thv_x0)/thv_en - 2._kind_phys*w_cu*(w_cu-w_eg);
      cquad =  2._kind_phys*rbuoy*g*cridis*(thv_x0 - thv_en)/thv_en +       w_cu**2;
      if( ( bquad**2-4._kind_phys*aquad*cquad ) .ge. 0._kind_phys ) then
         call roots(aquad,bquad,cquad,xs1,xs2,status)
         xc = min(1._kind_phys,max(0._kind_phys,min(1._kind_phys,min(xs1,xs2))))
      else
         xc = 1._kind_phys;
      endif
   else
      ! -------------------------------------------------- !
      ! Case 2 : When either cumulus or env. is saturated. !
      ! -------------------------------------------------- !
      ! May.03.2011. I commented out below 'xs' since it is already computed above.
      !              This will save computation time.
      ! xs  = excess_cu / ( excess_cu - excess_eg );
      thlxs = thl_cu + xs * ( thl_eg - thl_cu );
      qtxs  = qt_cu  + xs * ( qt_eg - qt_cu );
      call conden(p,thlxs,qtxs,th,qv,ql,qi,qse,id_check)
      thvxs = th * ( 1._kind_phys + zvir * qv - ql - qi )
      ! -------------------------------------------------- !
      ! kk=1 : Cumulus Segment, kk=2 : Environment Segment !
      ! -------------------------------------------------- !
      do kk = 1, 2
         if( kk .eq. 1 ) then
            thv_x0 = thv_cu;
            thv_x1 = ( 1._kind_phys - 1._kind_phys / xs ) * thv_cu + ( 1._kind_phys / xs ) * thvxs;
         else
            thv_x1 = thv_eg;
            thv_x0 = ( xs / ( xs - 1._kind_phys ) ) * thv_eg + ( 1._kind_phys/( 1._kind_phys - xs ) ) * thvxs;
         endif
         ! aquad =  w_cu**2._kind_phys;
         ! bquad =  2._kind_phys*rbuoy*g*cridis*(thv_x1 - thv_x0)/thv_en - 2._kind_phys*w_cu**2._kind_phys;
         aquad =  (w_cu-w_eg)**2._kind_phys;
         bquad =  2._kind_phys*rbuoy*g*cridis*(thv_x1 - thv_x0)/thv_en - 2._kind_phys*w_cu*(w_cu-w_eg);
         ! Below is the original UW code assuming 'thv_eg=thv_en'
         ! cquad =  2._kind_phys*rbuoy*g*cridis*(thv_x0 - thv_eg)/thv_en +       w_cu**2._kind_phys;
         ! Below is the revised code considering the difference between 'thv_eg' and 'thv_en'
         cquad =  2._kind_phys*rbuoy*g*cridis*(thv_x0 - thv_en)/thv_en +       w_cu**2._kind_phys;
         if( kk .eq. 1 ) then
            if( ( bquad**2._kind_phys-4._kind_phys*aquad*cquad ) .ge. 0._kind_phys ) then
               call roots(aquad,bquad,cquad,xs1,xs2,status)
               x_cu = min(1._kind_phys,max(0._kind_phys,min(xs,min(xs1,xs2))))
            else
               x_cu = xs;
            endif
         else
            if( ( bquad**2._kind_phys-4._kind_phys*aquad*cquad) .ge. 0._kind_phys ) then
               call roots(aquad,bquad,cquad,xs1,xs2,status)
               x_en = min(1._kind_phys,max(0._kind_phys,max(xs,min(xs1,xs2))))
            else
               x_en = 1._kind_phys;
            endif
         endif
      enddo
      if( x_cu .eq. xs ) then
         xc = max(x_cu, x_en);
      else
         xc = x_cu;
      endif
   endif

end subroutine buosorts_UW

!--------------------------------------------------------------------------------------------------

subroutine positive_moisture( &
   cp, xlv, xls, pcols, ncol, mkx, dt, qvmin, qlmin, qimin, dp, qv, ql, qi, t, s, qvten, &
   qlten, qiten, sten )

   ! ------------------------------------------------------------------------------- !
   ! Author : Sungsu Park. AMP/CGD/NCAR.                                             !
   ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
   ! force them to be larger than minimum value by (1) condensating water vapor      !
   ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
   ! layer. '2._kind_phys' is multiplied to the minimum values for safety.                  !
   ! Update final state variables and tendencies associated with this correction.    !
   ! If any condensation happens, update (s,t) too.                                  !
   ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
   ! input tendencies.                                                               !
   ! Be careful the order of k : '1': near-surface layer, 'mkx' : top layer          ! 
   ! May.03.2011. Additional refinement is added in the lowest model layer for       !
   !              complete treatment.                                                !
   ! ------------------------------------------------------------------------------- !

   integer,  intent(in)     :: pcols, ncol, mkx
   real(kind_phys), intent(in)     :: cp, xlv, xls
   real(kind_phys), intent(in)     :: dt, qvmin, qlmin, qimin
   real(kind_phys), intent(in)     :: dp(pcols,mkx)
   real(kind_phys), intent(inout)  :: qv(pcols,mkx), ql(pcols,mkx), qi(pcols,mkx), t(pcols,mkx), s(pcols,mkx)
   real(kind_phys), intent(inout)  :: qvten(pcols,mkx), qlten(pcols,mkx), qiten(pcols,mkx), sten(pcols,mkx)
   integer   i, k
   real(kind_phys)  dql, dqi, dqv, sum, aa, dum 

   do i = 1, ncol
      do k = mkx, 1, -1    ! From the top to the 1st (lowest) layer from the surface
         dql        = max(0._kind_phys,1._kind_phys*qlmin-ql(i,k))
         dqi        = max(0._kind_phys,1._kind_phys*qimin-qi(i,k))
         qlten(i,k) = qlten(i,k) +  dql/dt
         qiten(i,k) = qiten(i,k) +  dqi/dt
         qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
         sten(i,k)  = sten(i,k)  + xlv * (dql/dt) + xls * (dqi/dt)
         ql(i,k)    = ql(i,k) +  dql
         qi(i,k)    = qi(i,k) +  dqi
         qv(i,k)    = qv(i,k) -  dql - dqi
         s(i,k)     = s(i,k)  +  xlv * dql + xls * dqi
         t(i,k)     = t(i,k)  + (xlv * dql + xls * dqi)/cp
         dqv        = max(0._kind_phys,1._kind_phys*qvmin-qv(i,k))
         qvten(i,k) = qvten(i,k) + dqv/dt
         qv(i,k)    = qv(i,k)    + dqv
         if( k .ne. 1 ) then 
            qv(i,k-1)    = qv(i,k-1)    - dqv*dp(i,k)/dp(i,k-1)
            qvten(i,k-1) = qvten(i,k-1) - dqv*dp(i,k)/dp(i,k-1)/dt
         endif
         qv(i,k) = max(qv(i,k),qvmin)
         ql(i,k) = max(ql(i,k),qlmin)
         qi(i,k) = max(qi(i,k),qimin)
      end do
      ! May.03.2011. Below block is additionally added for completeness.
      ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally
      ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
      ! preserves column moisture.
      if( dqv .gt. 0._kind_phys ) then
         sum = 0._kind_phys
         do k = 1, mkx
            if( qv(i,k) .gt. 2._kind_phys*qvmin ) sum = sum + qv(i,k)*dp(i,k)
         enddo
         aa = dqv*dp(i,1)/max(1.e-20_kind_phys,sum)
         if( aa .lt. 0.5_kind_phys ) then
            do k = 1, mkx
               if( qv(i,k) .gt. 2._kind_phys*qvmin ) then
                  dum        = aa*qv(i,k)
                  qv(i,k)    = qv(i,k) - dum
                  qvten(i,k) = qvten(i,k) - dum/dt
               endif
            enddo
         else
!            write(iulog,*) 'Full positive_moisture is impossible in UNICON'
         endif
      endif
   end do

end subroutine positive_moisture

!--------------------------------------------------------------------------------------------------

subroutine positive_tracer( pcols, ncol, mkx, dt, trmin, dp, tr, trten )

   ! ------------------------------------------------------------------------------- !
   ! If any 'tr < trmin' are developed in any layer, force them to be larger than    !
   ! minimum value by transporting water vapor from the very lower layer.            !
   ! Update final state variables and tendencies associated with this correction.    !
   ! Note that 'tr' is the final state variables after applying corresponding        !
   ! input tendencies.                                                               !
   ! Be careful the order of k : '1': near-surface layer, 'mkx' : top layer          !
   ! May.03.2011. Additional refinement is added in the lowest model layer for       !
   !              complete treatment.                                                ! 
   ! ------------------------------------------------------------------------------- !

   integer,  intent(in)     :: pcols, ncol, mkx
   real(kind_phys), intent(in)     :: dt, trmin
   real(kind_phys), intent(in)     :: dp(pcols,mkx)
   real(kind_phys), intent(inout)  :: tr(pcols,mkx)
   real(kind_phys), intent(inout)  :: trten(pcols,mkx)
   integer   i, k
   real(kind_phys)  dtr, sum, aa, dum 

   do i = 1, ncol
      do k = mkx, 1, -1    ! From the top to the 1st (lowest) layer from the surface
         dtr = max(0._kind_phys,1._kind_phys*trmin-tr(i,k))
         trten(i,k) = trten(i,k) + dtr/dt
         tr(i,k)    = tr(i,k)    + dtr
         if( k .ne. 1 ) then 
            tr(i,k-1)    = tr(i,k-1)    - dtr*dp(i,k)/dp(i,k-1)
            trten(i,k-1) = trten(i,k-1) - dtr*dp(i,k)/dp(i,k-1)/dt
         endif
         tr(i,k) = max(tr(i,k),trmin)
      end do
      ! May.03.2011. Below block is additionally added for completeness.
      ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally
      ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
      ! preserves column moisture.
      if( dtr .gt. 0._kind_phys ) then
         sum = 0._kind_phys
         do k = 1, mkx
            if( tr(i,k) .gt. 2._kind_phys*trmin ) sum = sum + tr(i,k)*dp(i,k)
         enddo
         aa = dtr*dp(i,1)/max(1.e-20_kind_phys,sum)
         if( aa .lt. 0.5_kind_phys ) then
            do k = 1, mkx
               if( tr(i,k) .gt. 2._kind_phys*trmin ) then
                  dum        = aa*tr(i,k)
                  tr(i,k)    = tr(i,k) - dum
                  trten(i,k) = trten(i,k) - dum/dt
               endif
            enddo
         else
            ! write(iulog,*) 'Full positive_tracer is impossible in UNICON'
         endif
      endif
   end do

end subroutine positive_tracer

!--------------------------------------------------------------------------------------------------

subroutine findsp_single( q, t, p, tsp, qsp, i_in, k_in, lchnk )

   ! Wrapper for the findsp subroutine in wv_saturation.

   integer , intent(in) :: lchnk     
   integer , intent(in) :: i_in, k_in      
   real(kind_phys), intent(in) :: q      ! Water vapor [kg/kg]
   real(kind_phys), intent(in) :: t      ! Temperature [K]
   real(kind_phys), intent(in) :: p      ! Pressure    [Pa]

   real(kind_phys), intent(out) :: tsp   ! Saturation temp [K]
   real(kind_phys), intent(out) :: qsp   ! Saturation mixing ratio [kg/kg]

   logical :: use_ice = .true.
   integer :: status
   !------------------------------------------------------------------------------

   call findsp(q, t, p, use_ice, tsp, qsp, status)

   ! Currently, only 2 and 8 seem to be treated as fatal errors.
   if (status == 2) then
!      write(iulog,*) ' findsp not converging at i,k,lchnk = ', i_in, k_in, lchnk
!      write(iulog,*) ' t, q, p ', t, q, p
!      write(iulog,*) ' tsp, qsp ', tsp, qsp
!      call endrun('findsp_single:: not converging')
       call abort
   else if (status == 8) then
!      write(iulog,*) ' the enthalpy is not conserved at i,k,lchnk = ', i_in, k_in, lchnk
!      write(iulog,*) ' t, q, p ', t, q, p
!      write(iulog,*) ' tsp, qsp ', tsp, qsp
!      call endrun('findsp_single:: enthalpy is not conserved')
       call abort
   endif
   
end subroutine findsp_single

!==================================================================================================
! Internal utilities
!==================================================================================================

subroutine roots(a,b,c,r1,r2,status)

   ! --------------------------------------------------------- !
   ! Subroutine to solve the second order polynomial equation. !
   ! I should check this subroutine later.                     !
   ! Done.                                                     !
   ! --------------------------------------------------------- !

   real(kind_phys), intent(in)  :: a
   real(kind_phys), intent(in)  :: b
   real(kind_phys), intent(in)  :: c
   real(kind_phys), intent(out) :: r1
   real(kind_phys), intent(out) :: r2
   integer , intent(out) :: status
   real(kind_phys)              :: q, rmin, rmax

   r1 = 0._kind_phys
   r2 = 0._kind_phys
   status = 0
   if(a .eq. 0) then                              ! Form b*x + c = 0
      if(b .eq. 0) then                           ! Failure: c = 0
         status = 1
      else                                        ! b*x + c = 0
         r1 = -c/b
         r2 = r1
      endif
   else
      if(b .eq. 0._kind_phys) then                       ! Form a*x**2 + c = 0
         if(a*c .gt. 0._kind_phys) then                  ! Failure: x**2 = -c/a < 0
            status = 2  
         else                                     ! x**2 = -c/a 
            r1 = sqrt(-c/a)
            r2 = -r1
         endif
      else                                        ! Form a*x**2 + b*x + c = 0
         if((b**2 - 4._kind_phys*a*c) .lt. 0._kind_phys) then   ! Failure, no real roots
            status = 3
         else
            q  = -0.5_kind_phys*(b + sign(1.0_kind_phys,b)*sqrt(b**2 - 4._kind_phys*a*c))
            r1 =  q/a
            r2 =  c/q
         endif
      endif
   endif
   
   rmin = min(r1,r2)
   rmax = max(r1,r2)
   r1 = rmin
   r2 = rmax

end subroutine roots

!--------------------------------------------------------------------------------------------------

end module unicon_utils
