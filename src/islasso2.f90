! without sandwich formula
subroutine islasso3(X, y, n, p, theta, se, cov, lambda, alpha, pi, estpi, h, &
& itmax, itmaxse, tol, sigma2, trace, adaptive, offset, conv, stand, intercept, &
& eta, mu, res, dev, weights, hi, edf, grad2)

implicit none

!external variables
integer :: n, p, estpi, itmax, itmaxse, trace, conv, adaptive, stand, intercept
double precision :: X(n,p), y(n), theta(p), se(p), cov(p,p), lambda(p), alpha
double precision :: pi(p), h, tol, sigma2, offset(n), eta(n), mu(n), res(n), dev
double precision :: weights(n), hi(p), edf, grad2(p)

! internal variables
integer :: i, j, k, info
double precision :: xm(p), xse(p), xtw(p,n), xtx(p,p), xtwy(p), lambdaadapt(p)
double precision :: theta0(p), se0(p), grad(p), hess(p,p)
double precision :: invH(p,p), tempMat(p,p), cov1(p,p), redf, s2, ind, ind2

ind = 1.d0
ind2 = 1.d0

conv = 0
info = 0

xm = 0.d0
xse = 1.d0

hi = 1.d0
lambdaadapt = lambda / hi
s2 = sigma2

h = 0.1d0

if(stand.eq.1) then
    call standardize(X, xm, xse, n, p, intercept)
    lambdaadapt = lambdaadapt / xse
end if

! cmputing starting values
do i = 1, p
    xtw(i,:) = X(:,i) * weights
end do
call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
call DGEMV('N', p, n, 1.d0, xtw, p, y - offset, 1, 0.d0, xtwy, 1)

theta0 = theta
se0 = se

! main body
do i = 1, itmax
    if(trace.eq.9) call islasso_trace2_3(i)
    call rchkusr()

    if(adaptive.eq.1) lambdaadapt = lambda / (hi + 0.0000d1)

    do j = 1, itmax
        ! computing mixture parameters c
        if((estpi.eq.1).and.(i.gt.1)) then
            call logitlinkinv(abs(theta / se), p, pi)
            pi = 2.d0 * pi - 1.d0
        end if
        
        ! computing IWLS
        call fn1(theta, se, lambdaadapt, xtx, pi, p, hess, alpha)
        grad = xtwy
        call solve(hess, grad, p, info)
        if(info.ne.0) then
            conv = 2
            exit
        end if
        theta = grad

        ind2 = MAXVAL(abs(theta - theta0))
        if(ind2.le.tol) then
            exit
        end if

        theta0 = theta
    end do
    if((conv.eq.2).or.(itmaxse.eq.0)) exit
    
    call fn1(theta, se, lambdaadapt, xtx, pi, p, hess, alpha)
    call inv(p, hess, cov1, info)
    cov = cov + h * (cov1 - cov)
    do k = 1, p
        se(k) = sqrt(cov(k,k))
    end do
    
    ! checking possible convergence criterion
    ind = MAXVAL(abs(se - se0))
    if(trace.eq.2) call islasso_trace1_2_2(tol, i, MAXVAL(lambda), dev, redf, s2, &
        & ind, ind2)
    if(trace.eq.1) call islasso_trace1_7_2(tol, i, MAXVAL(lambda), dev, redf, s2, &
        & ind, ind2)

    if(ind.le.(tol*10)) then
        if((trace.eq.1).or.(trace.eq.2)) call islasso_trace1_8(1)
        exit
    end if
    
    ! conv = 1 if i >= itmax
    if(i.ge.itmaxse) then
        conv = 1
        exit
    end if
    
    se0 = se
end do

! updating output components
itmax = i
lambda = lambdaadapt
tol = ind

! if standardized beta and se, then return to the original scale
if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    do j = 1, p
        X(:,j) = X(:,j) * xse(j) + xm(j)
        xtw(j,:) = xtw(j,:) * xse(j) + xm(j)
        se(j) = sqrt(cov(j,j))
    end do
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
    lambda = lambda * xse
end if

call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
eta = eta + offset
mu = eta
res = y - mu
dev = sum(weights * (res**2))

! updating components for variance covariance matrix
call fn2(theta, se, lambda, xtx, pi, p, hess, alpha)
call inv(p, hess, invH, info)
if(info.ne.0) then
    conv = 2
end if
call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, tempMat, p)
do k = 1, p
    hi(k) = tempMat(k,k)
end do
edf = sum(hi)
redf = n - edf
if(sigma2.lt.0) sigma2 = dev / redf
cov = cov / sigma2
call gradient(theta, se, lambda, xtw, res, pi, n, p, grad2, alpha)

end subroutine islasso3


subroutine fn1(theta,se,lambda,xtx,pi,p,hess,alpha)
implicit none
integer :: p,i
double precision :: theta(p),se(p),lambda(p),xtx(p,p),pi(p),hess(p,p)
double precision :: theta2(p),pnm,temp1,alpha
hess = xtx
theta2 = theta
do i = 1, p
    if(abs(theta2(i)).lt.0.000000000001d0) theta2(i) = 0.000001d0
    temp1 = theta2(i) / se(i)
    hess(i,i) = hess(i,i) + lambda(i) * alpha * ( pi(i) * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + (1.d0 - pi(i)) * &
        & (2.d0 * pnm(temp1, 0.d0, 0.00001d0) - 1.d0) ) / theta2(i) + lambda(i) * (1.d0 - alpha)
end do
end subroutine fn1


subroutine fn2(theta,se,lambda,xtx,pi,p,hess,alpha)
implicit none
integer :: p,i
double precision :: theta(p),se(p),lambda(p),xtx(p,p),pi(p),hess(p,p)
double precision :: dnm,temp1,alpha
hess = xtx
do i = 1, p
    temp1 = theta(i) / se(i)
    hess(i,i) = hess(i,i) + 2.d0 * lambda(i) * alpha * ( pi(i) * dnm(temp1, 0.d0, 1.d0) + &
        & (1.d0 - pi(i)) * dnm(temp1, 0.d0, 0.00001d0) ) / se(i) + (1.d0 - alpha) * lambda(i)
end do
end subroutine fn2

subroutine islasso_red(X, y, n, p, theta, se, lambda, alpha, pi, &
& itmax, tol, offset, conv, weights)

implicit none

!external variables
integer :: n,p,itmax,conv
double precision :: X(n,p),y(n),theta(p),se(p),lambda(p),alpha
double precision :: pi(p),tol,offset(n),weights(n)

! internal variables
integer :: i, j, info
double precision :: ind2,xtw(p,n),xtx(p,p),xtwy(p),theta0(p),hess(p,p),grad(p)

ind2 = 1.d0
conv = 0
info = 0

! cmputing starting values
do i = 1, p
    xtw(i,:) = X(:,i) * weights
end do
call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
call DGEMV('N', p, n, 1.d0, xtw, p, y - offset, 1, 0.d0, xtwy, 1)

theta0 = theta

! main body
do j = 1, itmax
    ! computing IWLS
    call fn1(theta, se, lambda, xtx, pi, p, hess, alpha)
    grad = xtwy
    call solve(hess, grad, p, info)
    if(info.ne.0) then
        conv = 2
        exit
    end if
    theta = grad

    ind2 = MAXVAL(abs(theta - theta0))
    if(ind2.le.tol) then
        exit
    end if

    theta0 = theta
end do

! updating output components
itmax = j
tol = ind2
end subroutine islasso_red

! with sandwich formula
subroutine islasso(X, y, n, p, theta, se, cov, lambda, alpha, pi, estpi, h, &
& itmax, itmaxse, tol, sigma2, trace, adaptive, offset, conv, stand, intercept, &
& eta, mu, res, dev, weights, hi, edf, grad2)

implicit none

!external variables
integer :: n, p, estpi, itmax, itmaxse, trace, conv, adaptive, stand, intercept
double precision :: X(n,p), y(n), theta(p), se(p), cov(p,p), lambda(p), alpha
double precision :: pi(p), h, tol, sigma2, offset(n), eta(n), mu(n), res(n), dev
double precision :: weights(n), hi(p), edf, grad2(p)

! internal variables
integer :: i, j, k, info
double precision :: xm(p), xse(p), xtw(p,n), xtx(p,p), xtwy(p), lambdaadapt(p)
double precision :: theta0(p), se0(p), grad(p), hess(p,p)
double precision :: invH(p,p), tempMat(p,p), cov1(p,p), redf, s2, ind, ind2

ind = 1.d0
ind2 = 1.d0

conv = 0
info = 0

xm = 0.d0
xse = 1.d0

hi = 1.d0
lambdaadapt = lambda / hi
s2 = sigma2

h = 0.1d0

! standardizing
if(stand.eq.1) then
    call standardize(X, xm, xse, n, p, intercept)
    lambdaadapt = lambdaadapt / xse
end if

! cmputing starting values
do i = 1, p
    xtw(i,:) = X(:,i) * weights
end do
call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
call DGEMV('N', p, n, 1.d0, xtw, p, y - offset, 1, 0.d0, xtwy, 1)

theta0 = theta
se0 = se

! main body
do i = 1, itmax
    if(trace.eq.9) call islasso_trace2_3(i)
    call rchkusr()

    if(adaptive.eq.1) lambdaadapt = lambda / (hi + 0.0000d1)

    do j = 1, itmax
        ! computing mixture parameters c
        if((estpi.eq.1).and.(i.gt.1)) then
            call logitlinkinv(abs(theta / se), p, pi)
            pi = 2.d0 * pi - 1.d0
        end if
        
        ! computing IWLS
        call fn1(theta, se, lambdaadapt, xtx, pi, p, hess, alpha)
        grad = xtwy
        call solve(hess, grad, p, info)
        if(info.ne.0) then
            conv = 2
            exit
        end if
        theta = grad

        ind2 = MAXVAL(abs(theta - theta0))
        if(ind2.le.tol) then
            exit
        end if

        theta0 = theta
    end do
    if((conv.eq.2).or.(itmaxse.eq.0)) exit
    
    call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
    eta = eta + offset
    mu = eta
    res = y - mu
    dev = sum(weights * (res**2))

    ! updating components for variance covariance matrix
    call fn2(theta, se, lambdaadapt, xtx, pi, p, hess, alpha)
    call inv(p, hess, invH, info)
    if(info.ne.0) then
        conv = 2
        exit
    end if
    call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, tempMat, p)
    call DGEMM('N', 'N', p, p, p, 1.d0, tempMat, p, invH, p, 0.d0, cov1, p)
    cov = cov + h * (cov1 - cov)
    
    do k = 1, p
      hi(k) = tempMat(k,k)
    end do
    edf = sum(hi)
    redf = n - edf
    if(sigma2.lt.0) s2 = dev / redf

    do k = 1, p
        se(k) = sqrt(s2 * cov(k,k))
    end do
    
    ! checking possible convergence criterion
    ind = MAXVAL(abs(se - se0))
    if(trace.eq.2) call islasso_trace1_2_2(tol, i, MAXVAL(lambda), dev, redf, s2, &
        & ind, ind2)
    if(trace.eq.1) call islasso_trace1_7_2(tol, i, MAXVAL(lambda), dev, redf, s2, &
        & ind, ind2)

    if(ind.le.(tol*10)) then
        if((trace.eq.1).or.(trace.eq.2)) call islasso_trace1_8(1)
        exit
    end if
    
    ! conv = 1 if i >= itmax
    if(i.ge.itmaxse) then
        conv = 1
        exit
    end if
    
    se0 = se
end do

! updating output components
itmax = i
lambda = lambdaadapt
tol = ind

! if standardized beta and se, then return to the original scale
if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    do j = 1, p
        X(:,j) = (X(:,j) + xm(j)) * xse(j)
        xtw(j,:) = (xtw(j,:) + xm(j)) * xse(j)
        se(j) = sqrt(cov(j,j))
    end do
    lambda = lambda * xse
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
    
    call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
    eta = eta + offset
    mu = eta
    res = y - mu
    dev = sum(weights * (res**2))

    ! updating components for variance covariance matrix
    call fn2(theta, se, lambda, xtx, pi, p, hess, alpha)
    call inv(p, hess, invH, info)
    if(info.ne.0) then
        conv = 2
    end if
    call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, tempMat, p)
end if

!do k = 1, p
!    hi(k) = tempMat(k,k)
!end do
!edf = sum(hi)
!redf = n - edf
!if(sigma2.lt.0) sigma2 = dev / redf
!cov = cov / sigma2
sigma2 = s2
call gradient(theta, se, lambda, xtw, res, pi, n, p, grad2, alpha)

end subroutine islasso
