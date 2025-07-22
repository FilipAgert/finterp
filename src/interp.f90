module interp
    implicit none
    integer, parameter :: kind=8
    private
    public :: kind, linspace, cubic_spline, lin_interpolator

    type, abstract :: interpolator
        real(kind), private, allocatable :: xs(:) !!x coordinates. Must be sorted.
        real(kind), private, allocatable :: ys(:) !!y coordinates
        logical, private :: is_setup = .false.
        contains
        procedure :: eval => eval_interp !!Evaluate the interpolator at a given point
        procedure(eval_idx_interface),deferred, private:: eval_idx!!Evaluate the interpolator at a given point
        procedure(setup_interpolator), deferred :: setup !!Setup the interpolator with given xs and ys
    end type interpolator

    type, extends(interpolator) :: cubic_spline
        real(kind), private, allocatable :: coeffs(:)

        contains
        procedure, pass(self), private :: eval_idx => eval_cubic_spline
        procedure, pass(self) :: setup => setup_cubic_spline
    end type cubic_spline


    type, extends(interpolator) :: lin_interpolator
        contains
        procedure, pass(self), private :: eval_idx => eval_linear_interpolator
        procedure, pass(self) :: setup => setup_linear_interpolator
    end type lin_interpolator


    abstract interface
        subroutine setup_interpolator(self, xs, ys)
            import:: interpolator, kind
            class(interpolator), intent(inout) :: self
            real(kind), intent(in) :: xs(:) !!x coordinates. Must be sorted.
            real(kind), intent(in) :: ys(:) !!y coordinates
        end subroutine

        !!Idx is the index of the last element in xs smaller than or equal to x.
        pure elemental real(kind) function eval_idx_interface(self, x, idx)
            import :: kind, interpolator
            class(interpolator), intent(in) :: self
            real(kind), intent(in) :: x
            integer, intent(in) :: idx

        end function
    end interface 
    contains
    pure elemental real(kind) function eval_interp(self, x)
        class(interpolator), intent(in) :: self
        real(kind), intent(in) :: x
        integer :: idx
        if(.not. self%is_setup) error stop "interpolator not setup, call setup first"
        if(x > maxval(self%xs,1) .or. x < minval(self%xs,1)) error stop "input value not within bounds of input array"
        idx = maxloc(self%xs - x, 1, mask=self%xs <= x) !!Idx of the element smaller than or equal to x
        eval_interp = self%eval_idx(x, idx)
    end function
    subroutine setup_cubic_spline(self, xs, ys)
        class(cubic_spline), intent(inout) :: self
        real(kind), intent(in) :: xs(:) !!x coordinates. Must be sorted.
        real(kind), intent(in) :: ys(:) !!y coordinates
        integer :: sz
        sz = size(xs)
        if(sz /= size(ys)) error stop "size of input arrays must be the same."
        allocate(self%xs(sz))
        allocate(self%ys(sz))
        self%xs = xs
        self%ys = ys
        allocate(self%coeffs((sz-1)*4)) !!4 coefficients for each segment
        self%coeffs = q_spline_coeffs(xs,ys)
        self%is_setup = .true.
    end subroutine setup_cubic_spline

    pure real(kind) elemental function eval_cubic_spline(self, x, idx)
        class(cubic_spline), intent(in) :: self
        real(kind), intent(in) :: x
        integer, intent(in) :: idx
        eval_cubic_spline = self%coeffs((idx-1)*4 + 1) * x**3 + &
                            self%coeffs((idx-1)*4 + 2) * x**2 + &
                            self%coeffs((idx-1)*4 + 3) * x + &
                            self%coeffs((idx-1)*4 + 4)
    end function 
    
    function q_spline_coeffs(xs, ys)
        !! Interpolates the function value at a given point x using a cubic spline.
        !! Given two arrays `xs` and `ys` (of the same size), representing
        !! sample points of a function `y = f(x)`, this function estimates
        !! the value at a point `x` using a cubic spline.
        !!
        !! The input array `xs` must be sorted in ascending order.
        !! The value of `x` must lie within the bounds of `xs`.
        !!
        !! Author: Filip Agert, 2025
        real(kind), intent(in) :: xs(:) !!x coordinates. Must be sorted.
        real(kind), intent(in) :: ys(:) !!y coordinates
        integer :: sz
        real(kind) :: b((size(xs)-1)*4)!!RHS of the linear system
        real(kind) :: q_spline_coeffs((size(xs)-1)*4)
        real(kind) :: A((size(xs)-1)*4,(size(xs)-1)*4)!!Matrix of the linear system
        integer :: ii, startidx
        sz = size(xs)
        if(sz /= size(ys)) error stop "size of input arrays must be the same."
        A = 0
        b = 0
        !! We have piecewise functions in between each pair of points
        !! ax^3 + bx^2 + cx + d
        !! We require that the function is continous, and that the first and second derivative are continous at each point.
        !! 3a1x^2 + 2b1x + c1 - 3a2x^2 - 2b2x - c2 = 0
        !! 6a1x + 2b1 - 6a2x - 2b2 = 0
        do ii = 1,sz-1
            startidx = 4*ii - 3
            A(startidx, startidx:startidx+3) = [xs(ii)**3,xs(ii)**2,xs(ii),1.0_kind]!function value at point to the left.
            b(startidx) = ys(ii)

            A(startidx+1, startidx:startidx+3) = [xs(ii+1)**3,xs(ii+1)**2,xs(ii+1),1.0_kind] !function value at point to the right.
            b(startidx+1) = ys(ii+1)

            !First derivative at point to the right minus first derivative at point to the left. (Continous first derivative)
            if(ii /= sz-1) then
                A(startidx+2, startidx:startidx+7) = [3*xs(ii+1)**2, 2*xs(ii+1), 1.0_kind, 0.0_kind, -3*xs(ii+1)**2, -2*xs(ii+1), -1.0_kind, 0.0_kind] 
                A(startidx+3, startidx:startidx+7) = [6*xs(ii+1), 2.0_kind, 0.0_kind,0.0_kind, -6*xs(ii+1), -2.0_kind, 0.0_kind, 0.0_kind] 
            endif
        end do
        A(size(A,1)-1, 1:2) = [6*xs(1),2.0_kind] !second derivative at the first point
        A(size(A,1), size(A,2)-3:size(A,2)-2) = [6*xs(sz),2.0_kind] !second derivative at the last point

        !Solve linear system.
        call solve_linsys(A,b,size(A,1))

        q_spline_coeffs = b
    end function

    function linspace(start, stop, num)
        real(kind), intent(in) :: start, stop
        integer, intent(in) :: num
        real(kind) :: linspace(num)
        integer :: i
        if (num < 1) error stop "num must be at least 1"
        do i = 1, num
            linspace(i) = start + (stop - start) * (i - 1) / (num - 1)
        end do
    end function

    
    
    real(kind=kind) function lin_int(x, xs, ys) 
        !! Linearly interpolates the function value at a given point x.
        !!
        !! Given two arrays `xs` and `ys` (of the same size), representing
        !! sample points of a function `y = f(x)`, this function estimates
        !! the value at a point `x` using linear interpolation.
        !!
        !! The input array `xs` must be sorted in ascending order.
        !! The value of `x` must lie within the bounds of `xs`.
        !!
        !! Author: Filip Agert, 2025
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

    subroutine solve_linsys(A,b,dim) !Solve Ax=b for x
        real(kind), dimension(dim,dim), intent(in) :: A
        real(kind), dimension(dim), intent(inout) :: b !!In: b. Out: Solution
        integer,intent(in) :: dim
        real(kind), dimension(dim,1) :: Bmat
        integer, dimension(dim) :: ipiv
        integer ::info
        external :: dgesv
        Bmat(:,1) = b

        ! Solve the linear system A * x = b, with 
        !   XTX (Np, Np).   P: (Np)    Xt(Np,Nv)     BE(Nv)
        !Call LAPACK routine dgesv to solve the linear system
        call dgesv(dim,1,A,dim, ipiv,Bmat,dim,info)
        if(info /= 0) then
            print*, "Error when solving linsys, info = ", info
        endif
        b = Bmat(:,1)
    end subroutine

    subroutine setup_linear_interpolator(self, xs, ys)
        class(lin_interpolator), intent(inout) :: self
        real(kind), intent(in) :: xs(:) !!x coordinates. Must be sorted.
        real(kind), intent(in) :: ys(:) !!y coordinates
        integer :: sz
        sz = size(xs)
        if(sz /= size(ys)) error stop "size of input arrays must be the same."
        allocate(self%xs(sz))
        allocate(self%ys(sz))
        self%xs = xs
        self%ys = ys
        self%is_setup = .true.
    end subroutine setup_linear_interpolator

    pure real(kind) elemental function eval_linear_interpolator(self, x, idx)
        class(lin_interpolator), intent(in) :: self
        real(kind), intent(in) :: x
        integer, intent(in) :: idx
        eval_linear_interpolator = lerp(x, self%xs(idx), self%ys(idx), self%xs(idx+1), self%ys(idx+1))

    end function 
    
end module