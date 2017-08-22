!---------------------------------------------------------
!  The code is run through a bash script
!---------------------------------------------------------
program main
!
      use environment
      use control
      use parametersDPG
      use GMP
      use data_structure3D
      use physics
      use uhm2
      use problem
!
      implicit none
      character(len=1024) :: argv
      real*8  :: err,rnorm,rvoid,factor
      integer :: mdle,i,kref,idec,nvoid,niter,ibeg,iend,nreflag,istep,nrdof_old,nrdof_new,nstop,idec_solve
      integer, dimension(2) :: flag
!
!----------------------------------------------------------------------
!
!  ...initialize environment
      call begin_environment
!
!  ...read in HP3D input files location (if option is not present, the default value is used)
!
!                             option label      // explanation                // default value     // parameter
      call get_option_string( '-file-control'    , 'Control file'              , './files/control'  , FILE_CONTROL)
      call get_option_string( '-file-geometry'   , 'Geometry file'             , './files/cube', FILE_GEOM   )
      call get_option_string( '-file-phys'       , 'Physics file'              , './files/physics'  , FILE_PHYS   )
      call get_option_string( '-file-refinement' , 'Refinement files location' , '../../files/ref'  , FILE_REFINE )
      call get_option_string( '-file-history'    , 'History file'              , './files/history'  , FILE_HISTORY)
      call get_option_string( '-file-err'        , 'Error file'                , './files/dump_err' , FILE_ERR    )
!
!  ...read in problem dependent parameters (defined in module parametersDPG,DPGH1)
!
!                              option label      // explanation                // default value     // parameter
      call get_option_int(    '-nord-add'           , 'NORD_ADD'                  , 2                  , NORD_ADD    )
      call get_option_int(    '-order-approx'       , 'ORDER_APPROX'              , 2                  , ORDER_APPROX)
      call get_option_int(    '-orderx'             , 'NPX'                       , 0                  , NPX         )
      call get_option_int(    '-ordery'             , 'NPY'                       , 0                  , NPY         )
      call get_option_int(    '-orderz'             , 'NPZ'                       , 0                  , NPZ         )
      call get_option_int(    '-isol'               , 'ISOL'                      , 1                  , ISOL        )
      call get_option_int(    '-comp'               , 'ICOMP_EXACT'               , 1                  , ICOMP_EXACT )
!
      call get_option_int(    '-inner-product'      , 'INNER_PRODUCT'             , 1                  , INNER_PRODUCT)
      call get_option_real(   '-mu'                 , 'MU'                        , 1.d0               , MU          )
      call get_option_real(   '-epsilon'            , 'EPSILON'                   , 1.d0               , EPSILON     )
      call get_option_real(   '-sigma'              , 'SIGMA'                     , 0.d0               , SIGMA       )
      call get_option_real(   '-omega'              , 'OMEGA'                     , 1.d0               , OMEGA       )
       call get_option_int(    '-ibc'               , 'IBCFlag'                   , 0                  , IBCFlag     )
      GAMMA_IMP = dsqrt(1.d0-(PI**2/OMEGA**2))

!
!  ...finalize
      call end_environment
!
!  ...print fancy header
      write(*,*)'                      '
      write(*,*)'// --  PRIMAL DPG FOR MAXWELL EQUATION -- //'
      write(*,*)'                      '
!
!  ...initialize problem
      call initialize

!     Kyungjoo's magic...
      UHM_VERBOSE            = .FALSE.
      UHM_SOLVER_REPORT_SHOW = .FALSE.
!
      call get_command(argv)
      call uhm_initialize(argv)
!
      call uhm_option_begin
      call uhm_option_end
!
      call uhm_time_begin(UHM_TIME_LEVEL_FRONT)
      call uhm_direct_lu_piv_initialize( &
!              UHM_DOUBLE, NR_RHS_PROB, 256, UHM_SOLVER_PTR)
              UHM_DOUBLE, MY_NR_RHS, 256, UHM_SOLVER_PTR)
!
!  ...display menu in infinite loop
 10   continue
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      write(*,*) 'Quit ...................................0'
      write(*,*) '                                         '
      write(*,*) 'Geometry graphics (X11) ................1'
      write(*,*) 'HP3D graphics (X11) ....................3'
      write(*,*) 'Print Data Structure arrays ............5'
      write(*,*) 'Dumpout Data Structure .................7'
      write(*,*) '                                         '
      write(*,*) ' --  Geometry & Refinements  --          '
      write(*,*) 'Geometry error ........................12'
      write(*,*) 'Interactive H-refinements .............31'
      write(*,*) 'Uniform     H-refinements .............32'
      write(*,*) 'Adaptive    H-refinements .............33'
      write(*,*) '                                         '
      write(*,*) 'Solve (frontal) .......................40'
      write(*,*) 'Solve (MUMPS) .........................50'
      write(*,*) 'Solve (UHM) ...........................60'
      write(*,*) '                                         '
      write(*,*) 'Compute Hcurl error ..................100'
      write(*,*) 'Compute residual .....................110'
      write(*,*) 'Adaptive DPG refinements .............120'
      write(*,*) 'Compute BC data interpolation error...130'
      write(*,*) '                                         '
      write(*,*) 'My tests..... ........................200'
      write(*,*) 'My own tests..... ....................210'
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
!
      read(*,*) idec
!
!----------------------------------------------------------------------
!
      select case(idec)
!  ...quit
      case( 0) ; call finalize ; stop
!
!  ...GMP graphics
      case( 1) ; call graphg
!
!  ...hp3d graphics
      case( 3) ; call graphb
!
!  ...print data structure
      case( 5) ; call result
!
!  ...dump out data structure
      case( 7)
        write(*,*)'dumping out GMP...'
        call dumpout_GMP
        write(*,*)'dumping out physics...'
        call dumpout_physics
        write(*,*)'dumping out HP3D...'
        call dumpout_hp3d('./files/dumpc3Dhp')
!
!  ...geometry error
      case(12) ; call geom_error(err,rnorm)
!
!  ...interactive refinements
      case(31)
        write(*,*)'Active elements:'
        mdle=0
        do i=1,NRELES
          call nelcon(mdle, mdle)
!
          select case(NODES(mdle)%type)
          case('mdlp') ; write(*,7000) mdle
          case('mdlb') ; write(*,7001) mdle
          case('mdln') ; write(*,7002) mdle
          case('mdld') ; write(*,7003) mdle
          endselect
 7000     format(' mdle = ',i6,' ; PRISM')
 7001     format(' mdle = ',i6,' ; BRICK')
 7002     format(' mdle = ',i6,' ; TET')
 7003     format(' mdle = ',i6,' ; PYRAMID')
!
        enddo
        call display_ref_kind
        write(*,7010)
 7010   format(' mdle,kref =')
        read(*,*) mdle,kref
!
!       refine element
        call refine(mdle,kref)
!       recover 1-irregular mesh, update geometry and Dirichlet dof's
        call close_mesh ; call update_gdof ; call update_ddof
!
!     uniform global H-refinements
      case(32)
!       refine elements
        call global_href
!       recover 1-irregular mesh, update geometry and Dirichlet dof's
        call close_mesh ; call update_gdof ; call update_Ddof
!
!     adaptive H-refinements
      case(33) ; call adaptivity_geom(0.3d0, nvoid,rvoid)
!
!     frontal solve
      case(40)
        call solve1(MY_NR_RHS)
!
!     MUMPS solve
      case(50)
        call mumps_solve_seq(MY_NR_RHS)
!
!     UHM solve
      case(60)
        call uhm_solve
        call uhm_solver_flush(UHM_SOLVER_PTR)
!
!     compute Hcurl error for the E-field only
      case(100)
        flag(1)=1; flag(2)=0
        call compute_error(flag,1)
!
!  ...compute residual
      case(110) ; !!!!!call compute_residual
!
!  ...adaptive DPG refinements
      case(120)
  333   write(*,7011)
 7011   format('main: SET INITIAL, FINAL STEP,', &
               ' REFINEMENT FLAG (1-href,2-pref,3-hpref), FACTOR, idec_solve')
        read(*,*) ibeg,iend,nreflag,factor,idec_solve
        if (nreflag.ne.1.and.nreflag.ne.2.and.nreflag.ne.3) go to 333
        istep=ibeg-1
        do while(istep.lt.iend)
          istep=istep+1
          nrdof_old = NRDOFSE
          !!!!!!!call adapt_DPG(idec_solve,istep,nreflag,factor,nstop)
          if (nstop.eq.1) exit
          nrdof_new = NRDOFSE
          if (nrdof_new.eq.nrdof_old) then
            istep=istep-1
            factor = factor*0.25d0
            write(*,7023) factor
 7023       format('main: NEW factor = ',e12.5)
            if (factor.lt.0.000001d0) exit
          endif
        enddo
!
!     compute BC data interpolation error
      case(130)
        !!!!call compute_BC_interp_error
!
      case(200) ; call my_tests
!
      case(210) ; call my_own_tests
      endselect
!
!  ...go back to menu
      goto 10
!
!
endprogram main

      subroutine my_tests
      end subroutine my_tests

!
