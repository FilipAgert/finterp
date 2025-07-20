program main
    use interp
    implicit none

    real(kind) :: xs(1), ys(2)
    real(kind) :: x
    x = 1
    x = lin_int(x,xs,ys)

    write(*,*) "Hello world!"

end program