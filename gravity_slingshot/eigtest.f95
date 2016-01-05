program eigtest
    complex A(3,3)
    real eigs(3)
    A(1,1) = cmplx(1,0)
    A(1,2) = cmplx(0,2)
    A(1,3) = cmplx(3,0)
    A(2,1) = cmplx(0,-2)
    A(2,2) = cmplx(5,0)
    A(2,3) = cmplx(1,-1)
    A(3,1) = cmplx(3,0)
    A(3,2) = cmplx(1,1)
    A(3,3) = cmplx(7,0)
    call heevd(A, eigs)
    write(*,*) eigs
end
