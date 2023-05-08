! without sandwich formula
subroutine islasso_glm2(X, y, n, p, theta, se, cov, lambda, alpha, pi, estpi, &
& h, itmax, itmaxse, tol, sigma2, trace, adaptive, offset, conv, stand, intercept, &
& eta, mu, dev, weights, hi, edf, fam, link, grad2)

implicit none

! external variables
integer :: n, p, estpi, itmax, itmaxse, trace, adaptive, conv, stand, intercept, fam, link
double precision :: X(n,p), y(n), theta(p), se(p), cov(p,p),lambda(p), alpha
double precision :: pi(p), h, tol, sigma2, offset(n), eta(n), mu(n), dev
double precision :: weights(n), hi(p), edf, grad2(p)

! internal variables
integer :: i, j, k, info ! linkinv, mueta, variance, trace2, stand2,
double precision :: xm(p), xse(p), theta0(p), varmu(n) !, eta2(n), mu2(n)
double precision :: mu_eta_val(n), res(n), z(n), w(n), se0(p), redf
double precision :: lambda2(p), s2, ind, ind2, cov1(p,p)
double precision :: xtw(p,n), xtx(p,p), xtwy(p), hess(p,p), grad(p), invH(p,p), tempMat(p,p)

hi = 1.d0

ind = 1.d0
ind2 = 1.d0

!variance = 4
!mueta = 3
!linkinv = 2
!trace2 = 0
!stand2 = 0

xm = 0.d0
xse = 1.d0

call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
eta = eta + offset
call family(fam, link, 2, eta, n, mu)
call family(fam, link, 4, mu, n, varmu)
call family(fam, link, 3, eta, n, mu_eta_val)
res = (y - mu) / mu_eta_val
!z = (eta - offset) + res
!w = weights * (mu_eta_val**2) / varmu

s2 = sigma2
lambda2 = lambda / hi

theta0 = theta
se0 = se

h = 0.1d0

if(stand.eq.1) then
    call standardize(X, xm, xse, n, p, intercept)
    lambda2 = lambda2 / xse
end if

do i = 1, itmax
    call rchkusr()

    if(adaptive.eq.1) lambda2 = lambda / (hi + 0.0000d1)

    z = (eta - offset) + res
    w = weights * (mu_eta_val**2) / varmu

    ! cmputing starting values
    do k = 1, p
        xtw(k,:) = X(:,k) * w
    end do
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
    call DGEMV('N', p, n, 1.d0, xtw, p, z, 1, 0.d0, xtwy, 1)

    do j = 1, itmax
        ! computing mixture parameters c
        if(estpi.eq.1) then
            call logitlinkinv(abs(theta / se), p, pi)
            pi = 2.d0 * pi - 1.d0
        end if
        
        ! computing IWLS
        call fn1(theta, se, lambda2, xtx, pi, p, hess, alpha)
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

!    eta2 = eta
    call family(fam, link, 2, eta, n, mu)
!    mu2 = mu
    call family(fam, link, 4, mu, n, varmu)
!    eta2 = eta
    call family(fam, link, 3, eta, n, mu_eta_val)
    res = (y - mu) / mu_eta_val
    
    call fn1(theta, se, lambda2, xtx, pi, p, hess, alpha)
    call inv(p, hess, cov1, info)
    cov = cov + h * (cov1 - cov)
    do k = 1, p
        se(k) = sqrt(cov(k,k))
    end do
    
    ! checking possible convergence criterion
    ind = MAXVAL(abs(se - se0))
    if(trace.eq.2) call islasso_trace2_2_2(tol, i, MAXVAL(lambda), dev, redf, s2, &
        & ind, ind2)
    if(trace.eq.1) call islasso_trace2_7_2(tol, i, MAXVAL(lambda), dev, redf, s2, &
        & ind, ind2)

    if(ind.le.(tol*10)) then
        if(trace.eq.9) call islasso_trace2_6(i)
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
lambda = lambda2
tol = ind

! if standardized beta and se, then return to the original scale
if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    do j = 1, p
        X(:,j) = (X(:,j) + xm(j)) * xse(j)
        se(j) = sqrt(sigma2 * cov(j,j))
    end do
    lambda = lambda * xse
    
    call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
    eta = eta + offset

    call family(fam, link, 2, eta, n, mu)
    call family(fam, link, 4, mu, n, varmu)
    call family(fam, link, 3, eta, n, mu_eta_val)
    res = (y - mu) / mu_eta_val
    z = (eta - offset) + res
    w = weights * (mu_eta_val**2) / varmu
    do k = 1, p
        xtw(k,:) = X(:,k) * w
    end do
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
end if

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

end subroutine islasso_glm2


subroutine islasso_glm_red(X, y, n, p, theta, se, lambda, alpha, pi, &
& itmax, tol, offset, conv, weights, fam, link)

implicit none

!external variables
integer :: n,p,itmax,conv,fam,link
double precision :: X(n,p),y(n),theta(p),se(p),lambda(p),alpha
double precision :: pi(p),tol,offset(n),weights(n)

! internal variables
integer :: i,j,info
double precision :: ind2,eta(n),mu(n),varmu(n),mu_eta_val(n),res(n),theta0(p)
double precision :: z(n),w(n),xtw(p,n),xtx(p,p),xtwy(p),hess(p,p),grad(p)

ind2 = 1.d0
conv = 0
info = 0

!variance = 4
!mueta = 3
!linkinv = 2

call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
eta = eta + offset
call family(fam, link, 2, eta, n, mu)
call family(fam, link, 4, mu, n, varmu)
call family(fam, link, 3, eta, n, mu_eta_val)
res = (y - mu) / mu_eta_val

theta0 = theta

z = (eta - offset) + res
w = weights * (mu_eta_val**2) / varmu

! cmputing starting values
do i = 1, p
    xtw(i,:) = X(:,i) * w
end do
call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
call DGEMV('N', p, n, 1.d0, xtw, p, z, 1, 0.d0, xtwy, 1)

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
itmax = i
tol = ind2
end subroutine islasso_glm_red

! with sandwich formula
subroutine islasso_glm(X, y, n, p, theta, se, cov, lambda, alpha, pi, estpi, &
& h, itmax, itmaxse, tol, sigma2, trace, adaptive, offset, conv, stand, intercept, &
& eta, mu, dev, weights, hi, edf, fam, link, grad2)

implicit none

! external variables
integer :: n, p, estpi, itmax, itmaxse, trace, adaptive, conv, stand, intercept, fam, link
double precision :: X(n,p), y(n), theta(p), se(p), cov(p,p),lambda(p), alpha
double precision :: pi(p), h, tol, sigma2, offset(n), eta(n), mu(n), dev
double precision :: weights(n), hi(p), edf, grad2(p)

! internal variables
integer :: i, j, k, info ! linkinv, mueta, variance, trace2, stand2,
double precision :: xm(p), xse(p), theta0(p), varmu(n) !, eta2(n), mu2(n)
double precision :: mu_eta_val(n), res(n), z(n), w(n), se0(p), redf
double precision :: lambda2(p), s2, ind, ind2, cov1(p,p)
double precision :: xtw(p,n), xtx(p,p), xtwy(p), hess(p,p), grad(p), invH(p,p), tempMat(p,p)

hi = 1.d0

ind = 1.d0
ind2 = 1.d0

!variance = 4
!mueta = 3
!linkinv = 2
!trace2 = 0
!stand2 = 0

xm = 0.d0
xse = 1.d0

if(stand.eq.1) then
    call standardize(X, xm, xse, n, p, intercept)
    lambda2 = lambda2 / xse
end if

call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
eta = eta + offset
call family(fam, link, 2, eta, n, mu)
call family(fam, link, 4, mu, n, varmu)
call family(fam, link, 3, eta, n, mu_eta_val)
res = (y - mu) / mu_eta_val
!z = (eta - offset) + res
!w = weights * (mu_eta_val**2) / varmu

s2 = sigma2
lambda2 = lambda / hi

theta0 = theta
se0 = se

h = 0.1d0

do i = 1, itmax
    call rchkusr()

    if(adaptive.eq.1) lambda2 = lambda / (hi + 0.0000d1)

    z = (eta - offset) + res
    w = weights * (mu_eta_val**2) / varmu

    ! cmputing starting values
    do k = 1, p
        xtw(k,:) = X(:,k) * w
    end do
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
    call DGEMV('N', p, n, 1.d0, xtw, p, z, 1, 0.d0, xtwy, 1)

    do j = 1, itmax
        ! computing mixture parameters c
        if(estpi.eq.1) then
            call logitlinkinv(abs(theta / se), p, pi)
            pi = 2.d0 * pi - 1.d0
        end if
        
        ! computing IWLS
        call fn1(theta, se, lambda2, xtx, pi, p, hess, alpha)
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

!    eta2 = eta
    call family(fam, link, 2, eta, n, mu)
!    mu2 = mu
    call family(fam, link, 4, mu, n, varmu)
!    eta2 = eta
    call family(fam, link, 3, eta, n, mu_eta_val)
    res = (y - mu) / mu_eta_val
    dev = sum(weights * (res**2))
    
    ! updating components for variance covariance matrix
    call fn2(theta, se, lambda2, xtx, pi, p, hess, alpha)
    call inv(p, hess, invH, info)
    if(info.ne.0) then
        conv = 2
        exit
    end if
    call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, tempMat, p)
    call DGEMM('N', 'N', p, p, p, 1.d0, tempMat, p, invH, p, 0.d0, cov1, p)
    cov = cov + h * (cov1 - cov)
    do k = 1, p
        se(k) = sqrt(cov(k,k))
    end do
    
    ! checking possible convergence criterion
    ind = MAXVAL(abs(se - se0))
    if(trace.eq.2) call islasso_trace2_2_2(tol, i, MAXVAL(lambda), dev, redf, s2, &
        & ind, ind2)
    if(trace.eq.1) call islasso_trace2_7_2(tol, i, MAXVAL(lambda), dev, redf, s2, &
        & ind, ind2)

    if(ind.le.(tol*10)) then
        if(trace.eq.9) call islasso_trace2_6(i)
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
lambda = lambda2
tol = ind

! if standardized beta and se, then return to the original scale
if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    do j = 1, p
        X(:,j) = X(:,j) * xse(j) + xm(j)
        se(j) = sqrt(cov(j,j))
    end do
    lambda = lambda * xse
    
    call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
    eta = eta + offset
    call family(fam, link, 2, eta, n, mu)
    call family(fam, link, 4, mu, n, varmu)
    call family(fam, link, 3, eta, n, mu_eta_val)
    res = (y - mu) / mu_eta_val
    dev = sum(weights * (res**2))
    z = (eta - offset) + res
    w = weights * (mu_eta_val**2) / varmu
    do k = 1, p
        xtw(k,:) = X(:,k) * w
    end do
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
    call fn2(theta, se, lambda, xtx, pi, p, hess, alpha)
    call inv(p, hess, invH, info)
    if(info.ne.0) then
        conv = 2
    end if
    call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, tempMat, p)
end if

do k = 1, p
    hi(k) = tempMat(k,k)
end do
edf = sum(hi)
redf = n - edf
if(sigma2.lt.0) sigma2 = dev / redf
cov = cov / sigma2
call gradient(theta, se, lambda, xtw, res, pi, n, p, grad2, alpha)

end subroutine islasso_glm
