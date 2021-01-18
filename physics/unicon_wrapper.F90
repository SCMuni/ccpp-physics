!>\file unicon_wrapper.F90
!! This file is unicon wrapper with CCPP interface coupling to FV3

module unicon_wrapper

use machine,  only: kind_phys
use unicon,   only: unicon_init, compute_unicon, positive_moisture, positive_tracer
use physcons  , g => con_g, cp => con_cp, xlv => con_hvap, xlf=>con_hfus,             &
                zvir=>con_fvirt, rair=>con_rd, mwh2o=>con_amw, mwdry=>con_amd

implicit none

public :: unicon_wrapper_init, unicon_wrapper_run, unicon_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine unicon_wrapper_init(mpirank, mpiroot, errmsg, errflg)

         implicit none

         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! initialize ccpp error handling variables
         errmsg = ''
         errflg = 0

         if (mpirank==mpiroot) then
            write(0,*) '-----------------------------------------------------------------'
            write(0,*) ' --- WARNING --- the CCPP unicon convection scheme is currently under development, use at your own risk --- WARNING ---'
            write(0,*) '-----------------------------------------------------------------'
         end if

      end subroutine unicon_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_unicon_wrapper_finalize Argument Table
!!
      subroutine unicon_wrapper_finalize()
      end subroutine unicon_wrapper_finalize

!> \defgroup unicon_group UNICON wrapper Module
!! This is the unicon wrapper
!>\defgroup unicon_wrapper UNICON wrapper Module  
!> \ingroup unicon_group
!! This is the UNICON wrapper Module
!! \section arg_table_unicon_wrapper_run Argument Table
!! \htmlinclude unicon_wrapper_run.html
!!
!>\section unicon_wrapper UNICON General Algorithm
!> @{
   subroutine unicon_wrapper_run(im,km,ntracer,dt,kdt,                                  &
                                 ps0, phii, p0, phil, dp0,     u0,  v0, qv0, ql0, qi0,  &
                                 tr0, t0,  kpblh, pblh, qflx, shflx, taux, tauy,        &
                                 area, xland,raincv,                        &
                                 cush, cushavg, cuorg, awk_PBL, delta_thl_PBL,          &
                                 delta_qt_PBL, delta_u_PBL, delta_v_PBL, delta_tr_PBL,  &
                                 cu_cmfum, cu_cmfr, cu_thlr, cu_qtr, cu_ur, cu_vr,      &
                                 cu_qlr, cu_qir, cu_trr, cu_cmfrd, cu_thlrd, cu_qtrd,   &
                                 cu_urd, cu_vrd, cu_qlrd, cu_qird,                      &
                                 errmsg,errflg)


   integer, parameter :: &
   mit  = 18,     &! Number of previous time steps saved for convective organization
   nseg = 1,      &! Number of updraft segments [ # ]
   nsub = 1,      &! Number of subgrid, even number [ # ], if nsub = 1 don't use stochastic profile
   ! For advecting organization-related variables
   n_org = 5                      ! Number of constituents
!==================================================================================================
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   integer, intent(in)   :: im                      !  Number of columns 
   integer, intent(in)   :: km                      !  k = 1 : Lowest layer, k = km : Top layer
   integer, intent(in)   :: ntracer                 !  Number of tracers
   integer, intent(in)   :: kdt                     !  Number of time step index
   real(kind_phys), intent(in) :: dt             ! Time step in seconds: 2*delta_t [s]
   real(kind_phys), intent(in) :: ps0(im,0:km)   ! Environmental pressure at the interface [Pa]
   real(kind_phys), intent(in) :: phii(im,km+1)  ! Environmental geop-height at the interface [m2/s2]
   real(kind_phys), intent(in) :: p0(im,km)      ! Environmental pressure at the mid-point [m]
   real(kind_phys), intent(in) :: phil(im,km)    ! Environmental geop-height   at the mid-point [m2/s2]
   real(kind_phys), intent(in) :: dp0(im,km)     ! Environmental layer pressure thickness  [Pa] > 0
   real(kind_phys), intent(inout) :: u0(im,km)      ! Environmental zonal wind [m/s]
   real(kind_phys), intent(inout) :: v0(im,km)      ! Environmental meridional wind [m/s]
   real(kind_phys), intent(inout) :: qv0(im,km)     ! Environmental water  vapor specific humidity [kg/kg]
   real(kind_phys), intent(inout) :: ql0(im,km)     ! Environmental liquid water specific humidity [kg/kg]
   real(kind_phys), intent(inout) :: qi0(im,km)     ! Environmental ice          specific humidity [kg/kg]
   real(kind_phys), intent(inout) :: tr0(im,km,ntracer) ! Environmental tracers [ #/kg, kg/kg ]
   real(kind_phys), intent(inout) :: t0(im,km)      ! Environmental temperature [K]
   integer, intent(in)   :: kpblh(im)              !  Layer index with PBL top in it or at the base interface
   real(kind_phys), intent(in) :: pblh(im)       ! PBL top height [ m ]

   real(kind_phys), intent(in) :: qflx(im)      ! Upward water vapor flux into atmosphere at surface [ kg/m2/s ]
   real(kind_phys), intent(in) :: shflx(im)     ! Upward sensible heat flux into atmosphere at surface [ J/m2/s ]
   real(kind_phys), intent(in) :: taux(im)      ! Upward zonal      wind stress into atmosphere at surface [ kg m/s /m2/s ] 
   real(kind_phys), intent(in) :: tauy(im)      ! Upward meridional wind stress into atmosphere at surface [ kg m/s /m2/s ] 
   real(kind_phys), intent(in) :: area(im)      ! Physical area of individual grid box [ meter*meter ]    
   real(kind_phys), intent(in) :: xland(im)     ! landmask: sea/land/ice=0/1/2

   real(kind_phys), intent(out) :: raincv(im)   ! convective rainfall amount on physics timestep [ m ]


   real(kind_phys), intent(inout) :: cush(im)  ! Cumulus top height [ m ]    --memory
   real(kind_phys), intent(inout) :: cushavg(im)    ! Mean cumulus top height weighted by updraft masss flux at surface [ m ]
   real(kind_phys), intent(inout) :: cuorg(im)  ! Covective organization parameter [ 0-1 ] --memory
   real(kind_phys), intent(inout) :: awk_PBL(im)   ! Wake area within PBL [ 0 - 1 ]
   real(kind_phys), intent(inout) :: delta_thl_PBL(im)   ! Difference of thl between off-wake region and grid-mean value averaged over the PBL [ K ]
   real(kind_phys), intent(inout) :: delta_qt_PBL(im)  ! Difference of qt  between off-wake region and grid-mean value averaged over the PBL [ kg/kg ]
   real(kind_phys), intent(inout) :: delta_u_PBL(im)   ! Difference of u   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   real(kind_phys), intent(inout) :: delta_v_PBL(im)   ! Difference of v   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   real(kind_phys), intent(inout) :: delta_tr_PBL(im,ntracer)        ! Difference of tr  between off-wake region and grid-mean value averaged over the PBL [ kg/kg, #/kg ]

   real(kind_phys), intent(inout) :: cu_cmfum(im,km)  !  The amount of mass involved in the updraft buoyancy sorting at the previous time step [ kg/s/m2 ]
   real(kind_phys), intent(inout) :: cu_cmfr(im,km)   !  The amount of detrained mass from convective updraft and downdraft at the previous time step [ kg/s/m2 ]
   real(kind_phys), intent(inout) :: cu_thlr(im,km)   ! Mass-flux weighted mean 'thl' of detrained mass from convective updraft and downdraft at the previous time step [ K ]
   real(kind_phys), intent(inout) :: cu_qtr(im,km)    ! Mass-flux weighted mean 'qt'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   real(kind_phys), intent(inout) :: cu_ur(im,km)     ! Mass-flux weighted mean 'u'   of detrained mass from convective updraft and downdraft at the previous time step [ m/s ]
   real(kind_phys), intent(inout) :: cu_vr(im,km)     ! Mass-flux weighted mean 'v'   of detrained mass from convective updraft and downdraft at the previous time step [ m/s ]
   real(kind_phys), intent(inout) :: cu_qlr(im,km)    ! Mass-flux weighted mean 'ql'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   real(kind_phys), intent(inout) :: cu_qir(im,km)    ! Mass-flux weighted mean 'qi'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   real(kind_phys), intent(inout) :: cu_trr(im,km,ntracer) ! Mass-flux weighted mean 'tr'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   real(kind_phys), intent(inout) :: cu_cmfrd(im,km)  ! The amount of detrained mass from convective downdraft at the previous time step [ kg/s/m2 ]
   real(kind_phys), intent(inout) :: cu_thlrd(im,km)  ! Mass-flux weighted mean 'thl' of detrained mass from convective downdraft at the previous time step [ K ]
   real(kind_phys), intent(inout) :: cu_qtrd(im,km)   ! Mass-flux weighted mean 'qt'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   real(kind_phys), intent(inout) :: cu_urd(im,km)    ! Mass-flux weighted mean 'u'   of detrained mass from convective downdraft at the previous time step [ m/s ]
   real(kind_phys), intent(inout) :: cu_vrd(im,km)    ! Mass-flux weighted mean 'v'   of detrained mass from convective downdraft at the previous time step [ m/s ]
   real(kind_phys), intent(inout) :: cu_qlrd(im,km)   ! Mass-flux weighted mean 'ql'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   real(kind_phys), intent(inout) :: cu_qird(im,km)   ! Mass-flux weighted mean 'qi'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   real(kind_phys) :: cu_trrd(im,km,ntracer)  ! Mass-flux weighted mean 'tr'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]


!??? hli  real(kind_phys) :: dpdry0(im,km)              ! Environmental layer dry pressure thickness  [Pa] > 0
   real(kind_phys) :: dpdry0(im,km)              ! Environmental layer dry pressure thickness  [Pa] > 0  -- (moist)
!   -- remove hydrometeor mass
   real(kind_phys) :: ast0(im,km)                ! Stratiform fractional area at the layer mid-point [ fraction ] --0
   real(kind_phys) :: tke0(im,0:km)              !  TKE at the interface [ m2/s2 ]  (qgrs(ntke)) 
   real(kind_phys) :: bprod0(im,0:km)            ! Buoyancy production at the interface [ m2/s3 ] --0 rom tke
!??? hli   integer,   :: ipbl(im)     !  If 1(0), PBL is CL(STL)  -- 1
!(convective)
   integer :: ipbl(im)     !  If 1(0), PBL is CL(STL)
   real(kind_phys) :: wstar(im)    ! Turbulent convective velocity scale within PBL [ m/s ]  -sqrt(tke(i,1)) 
   real(kind_phys) :: tkes(im)     ! Turbulent kinetic energy at surface directly from the UW PBL scheme [ m/s ]    (qgrs(ntke))
   real(kind_phys) :: went(im) ! Entrainment rate at the PBL top interface directly from the UW PBL scheme [ m/s]   --0
   real(kind_phys) :: aflx(im,ntracer)  ! Upward tracer fluxes into atmosphere at surface [ #/m2/s, kg/m2/s ] --0
   real(kind_phys) :: ocnfrac(im)  ! Ocean Fraction [ fraction ]
   real(kind_phys) :: landfrac(im) ! Land  Fraction [ fraction ]    
   real(kind_phys) :: icefrac(im)  ! Ice Fraction [ fraction ]
   real(kind_phys) :: sgh(im)      ! Standard deviation of subgrid topographic height [ meter ]  --0
   real(kind_phys) :: sgh30(im)    ! Standard deviation of subgrid topographic height at 30 s horizontal area [ meter ] --0

   ! Aug.08.2013. Evaporation of stratiform precipitation
   real(kind_phys) :: am_evp_st(im,km)  ! Evaporation area of stratiform precipitation [fraction]       --0
   real(kind_phys) :: evprain_st(im,km) ! Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0. --0
   real(kind_phys) :: evpsnow_st(im,km) ! Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0. --0
   ! Formal output variables

   real(kind_phys) :: am_u(im,km)  ! Updraft fractional area [ fraction ] 
   real(kind_phys) :: qlm_u(im,km) ! Area-weighted in-cloud LWC within updraft fractional area [ kg / kg ]
   real(kind_phys) :: qim_u(im,km) ! Area-weighted in-cloud IWC within updraft fractional area [ kg / kg ]
   real(kind_phys) :: am_d(im,km)  ! Downdraft fractional area [ fraction ] 
   real(kind_phys) :: qlm_d(im,km) ! Area-weighted in-cloud LWC within downdraft fractional area [ kg / kg ]
   real(kind_phys) :: qim_d(im,km) ! Area-weighted in-cloud IWC within downdraft fractional area [ kg / kg ]
   real(kind_phys) :: cmf_u(im,km) ! Upward convective mass flux at the interface [ kg / s / m2 ]
   real(kind_phys) :: slflx(im,0:km)  ! Net upward convective flux of liquid static energy [ J / s / m2 ]
   real(kind_phys) :: qtflx(im,0:km)  ! Net upward convective flux of total specific humidity [ kg / s / m2 ]

   real(kind_phys) :: qvten(im,km)    ! Tendency of water vapor specific humidity [ kg / kg / s ]
   real(kind_phys) :: qlten(im,km)    ! Tendency of liquid water mixing ratio [ kg / kg / s ]
   real(kind_phys) :: qiten(im,km)    ! Tendency of ice mixing ratio [ kg / kg / s ]
   real(kind_phys) :: sten(im,km)     ! Tendency of dry static energy [ J / kg / s ]
   real(kind_phys) :: uten(im,km)     ! Tendency of zonal wind [ m / s / s ]
   real(kind_phys) :: vten(im,km)     ! Tendency of meridional wind [ m / s / s ]
   real(kind_phys) :: trten(im,km,ntracer)   ! Tendency of tracers [ # / kg / s, kg / kg / s ]
   real(kind_phys) :: qrten(im,km)    ! Production rate of rain by lateral expels of cumulus condensate [kg/kg/s]
   real(kind_phys) :: qsten(im,km)    ! Production rate of snow by lateral expels of cumulus condensate [kg/kg/s]
   real(kind_phys) :: precip(im)      ! Precipitation flux at surface in flux unit [ m / s ]
   real(kind_phys) :: snow(im)        ! Snow flux at surface in flux unit [ m / s ]
   real(kind_phys) :: rice2(im)       ! Precipitation flux at surface in flux unit [ m / s ]
   real(kind_phys) :: evapc(im,km)    ! Evaporation rate of convective precipitation within environment [ kg/kg/s ]
   real(kind_phys) :: rqc(im,km)      ! Production rate of raw detrained LWC+IWC  [kg/kg/s] > 0
   real(kind_phys) :: rqc_l(im,km)    ! Production rate of raw detrained LWC      [kg/kg/s] > 0
   real(kind_phys) :: rqc_i(im,km)    ! Production rate of raw detrained IWC      [kg/kg/s] > 0
   real(kind_phys) :: rliq(im)        ! Vertical integral of 'rqc_out'   in flux unit [m/s]
   real(kind_phys) :: rice(im)        ! Vertical integral of 'rqc_l_out' in flux unit [m/s]
   real(kind_phys) :: rnc_l(im,km)    ! Production rate of raw detrained droplet number of cloud liquid droplets [#/kg/s] > 0
   real(kind_phys) :: rnc_i(im,km)    ! Production rate of raw detrained droplet number of cloud    ice droplets [#/kg/s] > 0
   real(kind_phys) :: cnt(im)         ! Cloud top  interface index ( ki = kpen )
   real(kind_phys) :: cnb(im)         ! Cloud base interface index ( ki = krel-1 )
   real(kind_phys) :: cmf_det(im,km)  ! Detrained mass flux only from convective updraft (not from environmental air) and downdraft [ kg / s / m2 ] 
   real(kind_phys) :: ql_det(im,km)   ! Detrained LWC without mixing with the environment ( flux-convergence & subsidence-detrainment consistent ) [ kg / kg ]
   real(kind_phys) :: qi_det(im,km)   ! Detrained LWC without mixing with the environment ( flux-convergence & subsidence-detrainment consistent ) [ kg / kg ]
   !jihoons20170412. optional input argument, random alpha value for stochastic convection.
!   real(kind_phys) ,optional    :: alpha_rand(im)


   ! -------------------------------------------------------- !
   ! Local variables !
   ! --------------- !
   ! fields in physics buffer
   integer   :: lchnk                               !  Numer of chucks
   real(kind_phys) :: tr0_c(im,km,ntracer)      !  Environmental tracers [ # / kg, kg / kg ]
   integer  :: i, k                      !  Column   index ; Vertical index for local fields
   integer :: seed, iend
   real(kind_phys) :: xls    ! Latent heat of sublimation
   real(kind_phys), parameter :: alpha_crit  = 2._kind_phys
   real(kind_phys) :: tmp
  !real(kind_phys) :: alpha(im,km,nsub/2+1)            !  alpha
   real(kind_phys) :: alpha(im)            !  alpha
   real(kind_phys) :: zs0(im,0:km)               ! Environmental height   at the interface [m]
   real(kind_phys) :: z0(im,km)                  ! Environmental height   at the mid-point [m]
   real(kind_phys) :: s0(im,km)                  ! Environmental dry static energy [J/kg]


   ! --------- !
   ! Main body !
   ! --------- !
   errmsg = ''
   errflg = 0


   xls   = xlv + xlf


   do k = 1, km
     do i = 1, im
       s0(i,k)=cp*t0(i,k)+g*z0(i,k)
     enddo
   enddo




  !seed = (state%pmid(1,km-2) - int(state%pmid(1,km-2)))      * 1000000000
   seed = (p0(1,km-2) - int(p0(1,km-2)))      * 1000000000
   do i = 1, im
      tmp=100._kind_phys
      do while (abs(tmp) .gt. alpha_crit)
        if(seed .eq. 0) then
          seed = 123456
        end if
        tmp=kind_phys_normal_01(seed)
      end do
      alpha(i)=abs(tmp)  !always positive
!      alpha(i,1,1)=abs(tmp)  !always positive
   end do
   
   call unicon_init(xlv, cp, xlf, zvir, rair, g, mwh2o/mwdry)


   do i=1,im
    zs0(i,km)=phii(i,km+1)/g
    do k=1,km
    zs0(i,k-1)=phii(i,k)/g
    z0 (i,k  )=phil(i,k)/g
    enddo
   enddo

   do i=1,im
     if(xland(i)==0)then
       ocnfrac(i) =100.
       landfrac(i)=0.
       icefrac(i) =0.
     elseif(xland(i)==1)then
       ocnfrac(i) =0.
       landfrac(i)=100.
       icefrac(i) =0.
     elseif(xland(i)==2)then
       ocnfrac(i) =0.
       landfrac(i)=0.
       icefrac(i) =100.
     endif
   enddo

   lchnk = 1
   iend  = 1

   s0(:,:)     =0.      ! Environmental dry static energy [J/kg]
   ast0(:,:)   =0.      ! Stratiform fractional area at the layer mid-point [ fraction ]
   tke0(:,:)   =0.      ! TKE at the interface [ m2/s2 ]
   bprod0(:,:) =0.      ! Buoyancy production at the interface [ m2/s3 ]
!??? hli   integer,   :: ipbl(im)     !  If 1(0), PBL is CL(STL)   use 0
   wstar(:)    =0.      ! Turbulent convective velocity scale within PBL [ m/s ]
   tkes(:)     =0.      ! Turbulent kinetic energy at surface directly from the UW PBL scheme [ m/s ]
   went(:)     =0.      ! Entrainment rate at the PBL top interface directly from the UW PBL scheme [ m/s]
   sgh(:)      =0.      ! Standard deviation of subgrid topographic height [ meter ] 
   sgh30(:)    =0.      ! Standard deviation of subgrid topographic height at 30 s horizontal area [ meter ]

   am_evp_st(:,: )=0.      ! Evaporation area of stratiform precipitation [fraction]
   evprain_st(:,:)=0.      ! Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0.
   evpsnow_st(:,:)=0.      ! Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0.
   cush(:)        =0.      ! Cumulus top height [ m ]
   cushavg(:)     =0.      ! Mean cumulus top height weighted by updraft masss flux at surface [ m ]
   cuorg(:)       =0.      ! Covective organization parameter [ 0-1 ]
   awk_PBL(:)     =0.      ! Wake area within PBL [ 0 - 1 ]
   delta_thl_PBL(:)=0.     ! Difference of thl between off-wake region and grid-mean value averaged over the PBL [ K ]
   delta_qt_PBL(:) =0.     ! Difference of qt  between off-wake region and grid-mean value averaged over the PBL [ kg/kg ]
   delta_u_PBL(:)  =0.     ! Difference of u   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   delta_v_PBL(:)  =0.     ! Difference of v   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
   delta_tr_PBL(:,:)=0.    ! Difference of tr between off-wake region and grid-mean value averaged over the PBL [ kg/kg, #/kg]
   cu_cmfum (:,:)  =0.     !  The amount of mass involved in the updraft buoyancy sorting at the previous time step [ kg/s/m2 ]
   cu_cmfr(:,:)    =0.     !  The amount of detrained mass from convective updraft and downdraft at the previous time step [ kg/s/m2 ]
   cu_thlr(:,:)    =0.     ! Mass-flux weighted mean 'thl' of detrained mass from convective updraft and downdraft at the previous time step [ K ]
   cu_qtr(:,:)     =0.     ! Mass-flux weighted mean 'qt'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   cu_ur(:,:)      =0.     ! Mass-flux weighted mean 'u'   of detrained mass from convective updraft and downdraft at the previous time step [ m/s ]
   cu_vr(:,:)      =0.     ! Mass-flux weighted mean 'v'   of detrained mass from convective updraft and downdraft at the previous time step [ m/s ]
   cu_qlr(im,km)   =0      ! Mass-flux weighted mean 'ql'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   cu_qir(:,:)     =0.     ! Mass-flux weighted mean 'qi'  of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   cu_trr(:,:,:)   =0.     ! Mass-flux weighted mean 'tr' of detrained mass from convective updraft and downdraft at the previous time step [ kg/kg ]
   cu_cmfrd(:,:)   =0.     ! The amount of detrained mass from convective downdraft at the previous time step [ kg/s/m2 ]
   cu_thlrd(:,:)   =0.     ! Mass-flux weighted mean 'thl' of detrained mass from convective downdraft at the previous time step [ K ]
   cu_qtrd(:,:)    =0.     ! Mass-flux weighted mean 'qt'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   cu_urd(:,:)     =0.     ! Mass-flux weighted mean 'u'   of detrained mass from convective downdraft at the previous time step [ m/s ]
   cu_vrd(:,:)     =0.     ! Mass-flux weighted mean 'v'   of detrained mass from convective downdraft at the previous time step [ m/s ]
   cu_qlrd(:,:)    =0.     ! Mass-flux weighted mean 'ql'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   cu_qird(:,:)    =0.     ! Mass-flux weighted mean 'qi'  of detrained mass from convective downdraft at the previous time step [ kg/kg ]
   cu_trrd(:,:,:)  =0.     ! Mass-flux weighted mean 'tr' of detrained mass from convective downdraft at the previous time step [ kg/kg ]

   call compute_unicon( im       , km        , iend    , ntracer    , mit     , dt,   kdt , &
                        ps0       , zs0        , p0        , z0       , dp0     , dpdry0  , &
                        t0        , s0         , qv0       , ql0      , qi0     , tr0     , & 
                        u0        , v0         , ast0      , tke0     , bprod0  ,           &
                        ipbl      , kpblh      , pblh      , wstar    , tkes    , went    , & 
                        qflx      , shflx      , taux      , tauy ,     aflx,               &
                        ocnfrac   , landfrac   , icefrac   , sgh      , sgh30   ,           &
                        area      , am_evp_st  , evprain_st, evpsnow_st,                    &
                        cush      , cushavg    , cuorg     ,                                &
                        awk_PBL                , delta_thl_PBL        , delta_qt_PBL      , & 
                        delta_u_PBL            , delta_v_PBL          , delta_tr_PBL      , &
                        cu_cmfum  , cu_cmfr    , cu_thlr   , cu_qtr   , cu_ur   , cu_vr   , &
                        cu_qlr    , cu_qir     , cu_trr    ,                                & 
                        cu_cmfrd  , cu_thlrd   , cu_qtrd   , cu_urd   , cu_vrd  ,           &
                        cu_qlrd   , cu_qird    , cu_trrd   ,                                & 
                        am_u      , qlm_u      , qim_u     ,                                &
                        am_d      , qlm_d      , qim_d     ,                                &
                        cmf_u     , slflx      , qtflx     ,                                & 
                        qvten     , qlten      , qiten     , trten    ,                     &
                        sten      , uten       , vten      ,                                &
                        qrten     , qsten      ,                                            &
                        rqc_l     , rqc_i      , rqc       , rnc_l    , rnc_i   ,           &
                        rliq      , rice2      , precip    , snow     , evapc   ,           &
                        cnt       , cnb        , cmf_det   , ql_det   , qi_det  ,           &
                        lchnk, alpha  )
!                       lchnk, alpha(:,1,1)  )

   ! -------------------------------------------------------- !  

   !qlten(:im,:km) = qlten(:im,:km) - rqc_l(:im,:km)
   !qiten(:im,:km) = qiten(:im,:km) - rqc_i(:im,:km)
   !trten(:im,:km,ixnumliq) = trten(:im,:km,ixnumliq) - rnc_l(:im,:km)
   !trten(:im,:km,ixnumice) = trten(:im,:km,ixnumice) - rnc_i(:im,:km)

   ! --------------------------------------------- !
   ! Prevent negative cloud condensate and tracers !
   ! --------------------------------------------- !

   !qv0_c(:im,:km)  = qv0(:im,:km) + qvten(:im,:km)*dt
   !ql0_c(:im,:km)  = ql0(:im,:km) + qlten(:im,:km)*dt
   !qi0_c(:im,:km)  = qi0(:im,:km) + qiten(:im,:km)*dt
   !t0_c(:im,:km)   =  t0(:im,:km) + (1._kind_phys/cp)*sten(:im,:km)*dt
   !s0_c(:im,:km)   =  s0(:im,:km) +  sten(:im,:km)*dt
   !do mt = 1, ntracer
   !   tr0_c(:im,:km,mt)  = tr0(:im,:km,mt) + trten(:im,:km,mt)*dt
   !enddo

   ! Note : Since 'positive_moisture' will only perform condensation not the evaporation, 
   !        we don't need to modify corresponding 'nl,ni'.
   !        Thus, current version is completely OK.

   !do mt = 1, ntracer
   !   trten_ori(:im,:km,mt) = trten(:im,:km,mt)
   !enddo

   !call positive_moisture( &
   !   cp, xlv, xls, im, im, &
   !   km, dt, qmin(1), qmin(2), qmin(3), &
   !   dp0, qv0_c, ql0_c, qi0_c, t0_c, &
   !   s0_c, qvten, qlten, qiten, sten )

   !do mt = 1, ntracer

   !   if( cnst_get_type_byind(mt) .eq. 'wet' ) then
   !      pdel0(:im,:km) = dp0(:im,:km)
   !   else
   !      pdel0(:im,:km) = dpdry0(:im,:km)
   !   endif

   !   if (cnst_is_mam_num(mt)) then
   !      trmin = 1.e-5_kind_phys
   !   else
   !      trmin = qmin(mt)
   !   end if

   !   call positive_tracer( im, im, km, dt, trmin, pdel0, tr0_c(:,:,mt), trten(:,:,mt) )

   !enddo

!  -- update u,v,t,q and tracers --

!   do k = 1, km
!      k_inv = km + 1 - k
!      ptend%q(:im,k_inv,1)        = qvten(:im,k)       ! Convective tendency of specific humidity
!      ptend%q(:im,k_inv,ixcldliq) = qlten(:im,k)       ! Convective tendency of liquid water mixing ratio
!      ptend%q(:im,k_inv,ixcldice) = qiten(:im,k)       ! Convective tendency of ice mixing ratio
!      ptend%s(:im,k_inv)      = sten(:im,k)        ! Convective tendency of static energy
!      ptend%u(:im,k_inv)      = uten(:im,k)        ! Convective tendency of zonal wind
!      ptend%v(:im,k_inv)      = vten(:im,k)        ! Convective tendency of meridional wind
!   end do
    do k=1,km
     do i=1,im
      t0(i,k)=t0(i,k)+dt*sten(i,k) ! dh=cp*dT+g*dz
     enddo
    enddo

    do k=1,km
     do i=1,im
      u0 (i,k)=u0 (i,k)+uten (i,k)*dt
      v0 (i,k)=v0 (i,k)+vten (i,k)*dt
      qv0(i,k)=qv0(i,k)+qvten(i,k)*dt
      ql0(i,k)=ql0(i,k)+qlten(i,k)*dt
      qi0(i,k)=qi0(i,k)+qiten(i,k)*dt
     enddo
    enddo

    do i=1,im
     raincv(i)=.001*precip(i)*dt
    enddo

  end subroutine unicon_wrapper_run


  real(kind_phys) function kind_phys_normal_01 (seed)
   implicit none

   real(kind_phys) :: r1
   real(kind_phys) :: r2
   real(kind_phys), parameter :: kind_phys_pi = 3.141592653589793_kind_phys
   integer (kind=4) :: seed
   real(kind_phys) :: x

   r1 = kind_phys_uniform_01(seed)
   r2 = kind_phys_uniform_01(seed)

   x=sqrt(-2._kind_phys *log(r1))*cos(2._kind_phys*kind_phys_pi*r2)

   kind_phys_normal_01 = x
   return
  end function

  real(kind_phys) function kind_phys_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      kind_phys_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
  implicit none

  integer  ::  k
  integer  :: seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
  kind_phys_uniform_01 = real ( seed, kind_phys ) *4.656612875e-10_kind_phys

  return
  end function

!> @}
  end module unicon_wrapper
