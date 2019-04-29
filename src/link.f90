! BINOMIAL
subroutine binomial_variance(x,n,varmu)
implicit none
integer :: n
double precision :: x(n),varmu(n)
varmu = x * (1.d0 - x)
end subroutine binomial_variance

subroutine logitlink(x,n,mu)
implicit none
integer :: n
double precision :: x(n),mu(n)
mu = log(x / (1.d0 - x))
end subroutine logitlink

subroutine logitlinkinv(x,n,eta)
implicit none
integer :: n,i
double precision :: x(n),eta(n),tmp
double precision, parameter :: thresh = 2.220446049250313E-16
double precision, parameter :: invthresh = 4503599627370496.d0
do i = 1, n
    if(x(i).lt.(-30.d0)) then
        tmp = thresh
    else if(x(i).gt.30.d0) then
        tmp = invthresh
    else
        tmp = exp(x(i))
    end if
    eta(i) = tmp / (1.d0 + tmp)
end do
end subroutine logitlinkinv

subroutine logitmueta(x,n,eta)
implicit none
integer :: n,i
double precision :: x(n),eta(n)
double precision, parameter :: thresh = 2.220446049250313E-16
do i = 1, n
    if((x(i).lt.(-30.d0)).or.(x(i).gt.30.d0)) then
        eta(i) = thresh
    else
        eta(i) = exp(x(i)) / ((1.d0 + exp(x(i)))**2)
    end if
end do
end subroutine logitmueta

subroutine probitlink(x,n,mu)
implicit none
integer :: n,i
double precision :: x(n),qnm,mu(n)
do i = 1, n
    mu(i) = qnm(x(i))
end do
end subroutine probitlink

subroutine probitlinkinv(x,n,eta)
implicit none
integer :: n,i
double precision :: x(n),pnm,eta(n)
double precision, parameter :: thresh = 8.12589066470190d6
do i = 1, n
    eta(i) = x(i)
    if(eta(i).le.(-thresh)) eta(i) = -thresh
    if(eta(i).ge.thresh) eta(i) = thresh
    eta(i) = pnm(eta(i))
end do
end subroutine probitlinkinv

subroutine probitmueta(x,n,eta)
implicit none
integer :: n,i
double precision :: x(n),dnm,eta(n)
double precision, parameter :: thresh = 2.220446049250313E-16
do i = 1, n
    eta(i) = dnm(x(i))
    if(eta(i).le.thresh) eta(i) = thresh
end do
end subroutine probitmueta

!POISSON
subroutine poisson_variance(x,n,varmu)
implicit none
integer :: n
double precision :: x(n),varmu(n)
varmu = x
end subroutine poisson_variance

subroutine loglink(x,n,mu)
implicit none
integer :: n
double precision :: x(n),mu(n)
mu = log(x)
end subroutine loglink

subroutine loglinkinv(x,n,eta)
implicit none
integer :: n,i
double precision :: x(n),eta(n)
double precision, parameter :: thresh = 2.220446049250313E-16
eta = exp(x)
do i = 1, n
    if(eta(i).le.thresh) eta(i) = thresh
end do
end subroutine loglinkinv

subroutine logmueta(x,n,eta)
implicit none
integer :: n,i
double precision :: x(n),eta(n)
double precision, parameter :: thresh = 2.220446049250313E-16
eta = exp(x)
do i = 1, n
    if(x(i).le.thresh) x(i) = thresh
end do
end subroutine logmueta

! GAMMA
subroutine gamma_variance(x,n,varmu)
implicit none
integer :: n
double precision :: x(n),varmu(n)
varmu = x**2
end subroutine gamma_variance

subroutine inverselink(x,n,mu)
implicit none
integer :: n
double precision :: x(n),mu(n)
mu = 1.d0 / x
end subroutine inverselink

subroutine inverselinkinv(x,n,eta)
implicit none
integer :: n
double precision :: x(n),eta(n)
eta = 1.d0 / x
end subroutine inverselinkinv

subroutine inversemueta(x,n,eta)
implicit none
integer :: n
double precision :: x(n),eta(n)
eta = -1.d0 / (x**2)
end subroutine inversemueta

subroutine identitylink(x,n,mu)
implicit none
integer :: n
double precision :: x(n),mu(n)
mu = x
end subroutine identitylink

subroutine identitylinkinv(x,n,eta)
implicit none
integer :: n
double precision :: x(n),eta(n)
eta = x
end subroutine identitylinkinv

subroutine identitymueta(x,n,eta)
implicit none
integer :: n
double precision :: x(n),eta(n)
x = x
eta = 1.d0
end subroutine identitymueta
