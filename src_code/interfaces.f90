!======================================================================
!   Fortran 90 interfaces for external subroutines ( not in modules)
!
!======================================================================

module interfaces
        implicit none
        
interface error_abort

       subroutine error_abort(msg)
               character(len=*), intent(in), optional :: msg
       end subroutine error_abort

end interface error_abort


interface SubcommSum

        subroutine SubcommSumVect(A, na, ixyz)
	        use messenger

                integer, intent(in) :: na, ixyz !Direction of sub-comm
                double precision A(na)

        end subroutine SubcommSumVect

        subroutine SubcommSumInt(A, ixyz)
	        use messenger

                integer, intent(in) :: ixyz !Direction of sub-comm
	        integer	A
        end
    
        subroutine SubcommSumIntVect(A, na, ixyz)
	        use messenger

                integer, intent(in) :: na, ixyz !Direction of sub-comm
         	integer	A(na)
        end

end interface SubcommSum

end module interfaces
