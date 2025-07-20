module interp
    use error
    implicit none
    integer, parameter :: kind=8
    private
    public :: kind, lin_int
    contains
    
    
    real(kind=kind) function lin_int(x, xs, ys) 
        real(kind), intent(in) :: x !!Point to interpolate at. 
        real(kind), intent(in) :: xs(:) !!x coordinates. Must be sorted.
        real(kind), intent(in) :: ys(:) !!y coordinates
        integer :: sz,ii
        sz = size(xs)
        if(sz /= size(ys)) then
            error stop "size of input arrays must be the same."
        endif
        do ii =1,sz
            if(xs(ii) < x) then
                lin_int = lerp(x, xs(ii), ys(ii), xs(ii+1), ys(ii+1))
                return
            endif
        end do
        !!We dont allow for extrapolation. only interpolation
        error stop "input value not within bounds of input array"
    end function

    pure elemental real(kind=kind) function lerp(x, x0, y0, x1, y1) !!Linear interpolation between two points
        real(kind), intent(in) :: x !!coordinate to linear interp at 
        real(kind), intent(in) :: x0,y0,x1,y1 !!Known coordinates
        lerp = y0 + (x-x0) * (y1-y0)/(x1-x0)
    end function

end module