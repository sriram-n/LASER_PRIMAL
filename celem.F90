!===============================================================================================================
!  REMARK : unless you are an advanced user, there is no reason to modify this routine.
!===============================================================================================================
!> Purpose : routine computes modified stiffness matrix and load vector for an element
!
!! @param[in]  Mdle   - an element (middle node) number
!! @param[in]  Idec   - 1 = info on nodes only ; 2 = compute element matrices as well
!! @param[out] Nrdofs - number of element local dof for each physical attribute (needed to dimension work
!!                      arrays in calling routines)
!! @param[out] Nrdofm - number of modified element dof in the expanded mode (needed to dimension work arrays
!!                      in calling routines)
!! @param[out] Nrdofc - number of modified element dof for the coupled problem after compression
!! @param[out] Nodm   - actual (unconstrained) nodes in the order : middle, mid-face, mid-edge, vertex nodes
!! @param[out] NdofmH - the corresponding number of H1      dof
!! @param[out] NdofmE - the corresponding number of H(curl) dof
!! @param[out] NdofmV - the corresponding number of H(div)  dof
!! @param[out] NdofmQ - the corresponding number of L2      dof
!! @param[out] Nrnodm - number of the modified element nodes
!! @param[out] Bload  - 1D array containing the modified load vector
!! @param[out] Astif  - 1D array containing the modified stiffness matrix
!---------------------------------------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine celem(Mdle,Idec, Nrdofs,Nrdofm,Nrdofc,Nodm,NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
!
      use physics
      use data_structure3D
!
      implicit none
      integer,                      intent(in)  :: Mdle,Idec
      integer, dimension(NR_PHYSA), intent(out) :: Nrdofs
      integer,                      intent(out) :: Nrdofm,Nrdofc
      integer, dimension(MAXNODM) , intent(out) :: Nodm
      integer, dimension(MAXNODM) , intent(out) :: NdofmH,NdofmE,NdofmV,NdofmQ
      integer,                      intent(out) :: Nrnodm
      VTYPE  ,                      intent(out) :: Bload(*),Astif(*)

      integer, dimension(NR_PHYSA) :: nbcond
!
!---------------------------------------------------------------------------------------------------------------
!
!  ...this is a hack to eliminate middle node dof for traces in celem_system by using BC flag
!  ...  We assume the following physical attributes for the Laser problem:
!           (1) Hcurl for Maxwell field, (2) Hcurl for Maxwell trace
!

      nbcond(1)=0; nbcond(2)=1;
      call encod(nbcond,10,NR_PHYSA, NODES(Mdle)%bcond)
      call set_index(NODES(Mdle)%case,NODES(Mdle)%bcond, NODES(Mdle)%index)
!
!     redirect to the system routine
      call celem_system(Mdle,Idec, Nrdofs,Nrdofm,Nrdofc,Nodm,NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
!
!  ...reset the BC flag to zero
      NODES(Mdle)%bcond = 0
      call set_index(NODES(Mdle)%case,NODES(Mdle)%bcond, NODES(Mdle)%index)
!
!
!
end subroutine celem
