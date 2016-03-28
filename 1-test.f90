real*8 function test(x)
    implicit none
    real*8 :: x
    test=x+2
    return
end function test

program main
    implicit none
  real*8,external::test
    real*8 a
    a=6
    a=test(a)
    write(*,*) "The Answer is : ",a
 end program main

real function p(args)
    implicit none
    real :: args
    
end function p
