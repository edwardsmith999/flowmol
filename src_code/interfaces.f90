!======================================================================
!   Fortran 90 interfaces for external subroutines ( not in modules)
!
!======================================================================

module interfaces
        implicit none
        
interface error_abort

       subroutine error_abort_s(msg)
               character(len=*), intent(in), optional :: msg
       end subroutine error_abort_s

       subroutine error_abort_si(msg,i)
               character(len=*), intent(in) :: msg
               integer, intent(in)          :: i
       end subroutine error_abort_si

end interface error_abort


interface SubcommSum

        subroutine SubcommSumVect(A, na, ixyz)

                integer, intent(in) :: na, ixyz !Direction of sub-comm
                real(kind(0.d0)) A(na)

        end subroutine SubcommSumVect

        subroutine SubcommSumInt(A, ixyz)

			integer, intent(in) :: ixyz !Direction of sub-comm
	        integer	A
        end
    
        subroutine SubcommSumIntVect(A, na, ixyz)

			integer, intent(in) :: na, ixyz !Direction of sub-comm
         	integer	A(na)
        end

end interface SubcommSum

end module interfaces
