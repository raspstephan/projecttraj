! Fortran subroutines for use in traj_tools python module
! Build using f2py -c futils.f90 -m futils

module futils

    contains

    subroutine minxspan(data, crit, span, istart, istop, startval, stopval)
        implicit none
        
        ! Declaration of I/O variables
        real, dimension(:), intent(in) :: data
        real, intent(in) :: crit
        integer, intent(inout) :: span, istart, istop
        real, intent(inout) :: startval, stopval
!f2py intent(in,out) :: span, istart, istop, startval, stopval
        
        ! Declaration of internal variables
        integer :: i, j, datalen, found
        
        
        ! Initial values
        datalen = size(data)
        
        !span = datalen + 1
        i = 1
        ! print *, maxval(data(i:)), data(i), maxval(data(i:)) - minval(data(i:)), crit
        do while ((i <= datalen) .and. ((maxval(data(i:)) - minval(data(i:))) > crit))
            ! print *, 'here'
            if (data(i) < data(i + 1)) then 
                j = 1
                found = 0
                do while ((i + j <= datalen) .and. (found == 0))
                    if (((data(i + j) - data(i)) >= crit) .and. (j < span)) then
                        istart = i
                        istop = i + j
                        span = j
                        startval = data(i)
                        stopval = data(i + j)
                        !print *, startval, stopval
                        found = 1
                    endif
                    j = j + 1
                end do
            endif
            i = i + 1
        end do
        
        ! Adjust for zero-based indexing in Python
        istart = istart - 1
        istop = istop - 1

    end subroutine minxspan

end module futils