!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_output
!**********************************************************************
!                                                                     *
!     Module for data output                                          *
!                                                                     *
!**********************************************************************
use mpi
use mg_intstate,only:mg_intstate_type

public output_spec1_2d
public output_spec1_2dd
public output_spec2_2d
public output_vertical_2d

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine output_spec1_2d                      &
!***********************************************************************
!                                                                      !
!   Outpyt a 2D array                                                  !
!                                                                      !
!***********************************************************************
(obj_mgbf,V)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
implicit none
type(mg_intstate_type)::obj_mgbf
real(r_kind),dimension(1:obj_mgbf%nm,1:obj_mgbf%mm),intent(in):: V
real(i_kind):: m,n
!-----------------------------------------------------------------------

    do m=1,obj_mgbf%mm
      write(100+obj_mgbf%mype,'(f9.3)') (V(n,m),n=1,obj_mgbf%nm)
    enddo

!-----------------------------------------------------------------------
                        endsubroutine output_spec1_2d


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine output_vertical_2d                   &
!***********************************************************************
!                                                                      !
!   Output of verical a 2D array                                       !
!                                                                      !
!***********************************************************************
(obj_mgbf,V,my0)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
implicit none
type(mg_intstate_type)::obj_mgbf
real(r_kind),dimension(1:obj_mgbf%nm,1:obj_mgbf%LM),intent(in):: V
integer(i_kind):: my0
real(i_kind):: n,l
!-----------------------------------------------------------------------

if(my==my0) then
    do L=1,obj_mgbf%Lm
      write(500+nx,'(f9.3)') (V(n,L),n=1,obj_mgbf%nm)
    enddo
endif
call barrierMPI

!-----------------------------------------------------------------------
                        endsubroutine output_vertical_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine output_spec2_2d                      &
!***********************************************************************
!                                                                      !
!   Outlyt a 2D array                                                  !
!                                                                      !
!***********************************************************************
(obj_mgbf,V)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
implicit none
type(mg_intstate_type)::obj_mgbf
real(r_kind),dimension(-1:obj_mgbf%imL+1,-1:obj_mgbf%jmL+1),intent(in):: V
real(i_kind):: j,i
!-----------------------------------------------------------------------

    do j=obj_mgbf%jmL+1,-1,-1
      write(100+obj_mgbf%mype,'(f9.3)') (V(i,j),i=-1,obj_mgbf%imL+1)
    enddo

!-----------------------------------------------------------------------
                        endsubroutine output_spec2_2d

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine output_spec1_2dd                     &
!***********************************************************************
!                                                                      !
!   Outpyt a 2D array                                                  !
!                                                                      !
!***********************************************************************
(obj_mgbf,V,nmax,mmax)
!-----------------------------------------------------------------------
use kinds, only: r_kind,i_kind
implicit none
type(mg_intstate_type)::obj_mgbf
real(r_kind),dimension(1:nmax,1:mmax),intent(in):: V
integer(i_kind),intent(in):: nmax,mmax
real(i_kind):: m,n
!-----------------------------------------------------------------------

    do m=1,mmax
      write(100+obj_mgbf%mype,'(f9.3)') (V(n,m),n=1,nmax)
    enddo

!-----------------------------------------------------------------------
                        endsubroutine output_spec1_2dd


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_output
