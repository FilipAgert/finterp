program main
    use interp
    implicit none

    real(kind) :: xs(4), ys(4)
    real(kind) :: y
    real(kind) :: x(10), splined(10), lined(10)
    type(cubic_spline) :: spline
    type(lin_interpolator) :: lin_interp
    integer ::ii
    x = linspace(2.0_kind,4.0_kind,10)
    xs = [2.0_kind, 3.0_kind, 4.0_kind,5.0_kind]
    ys = [3.0_kind,2.0_kind,4.0_kind,5.0_kind]
    call spline%setup(xs,ys)
    call lin_interp%setup(xs,ys)
    splined = spline%eval(x)
    lined = lin_interp%eval(x)
    write(*,'(A)') "x        spline     linear"
    do ii = 1,size(x)
         write(*,'(3F10.1)')x(ii), splined(ii),lined(ii)
    end do

end program