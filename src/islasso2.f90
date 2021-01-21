subroutine islasso2(X, y, n, p, theta, se, cov, lambda, alpha, pi, estpi, h, &
& itmax, tol, sigma2, trace, adaptive, offset, conv, stand, intercept, eta, &
& mu, res, dev, weights, hi, edf, grad2)

implicit none

!external variables
integer :: n, p, estpi, itmax, trace, conv, adaptive, stand, intercept
double precision :: X(n,p), y(n), theta(p), se(p), cov(p,p), lambda(p), alpha
double precision :: pi(p), h, tol, sigma2, offset(n), eta(n), mu(n), res(n), dev
double precision :: weights(n), hi(p), edf, grad2(p)

! internal variables
integer :: i, j, info
double precision :: xm(p), xse(p), f0, xtw(p,n), xtx(p,p), lambdaadapt(p)
double precision :: theta0(p), se0(p), grad(p), hess(p,p), f0_old, ind3, ind4
double precision :: invH(p,p), tempMat(p,p), cov1(p,p), redf, s2, ind, ind2, h2

conv = 0
info = 0

xm = 0.d0
xse = 1.d0

hi = 1.d0
lambdaadapt = lambda / hi
s2 = sigma2

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

! main body
do i = 1, itmax
    if(trace.eq.9) call islasso_trace2_3(i)
    call rchkusr()

    theta0 = theta
    se0 = se
    if(adaptive.eq.1) lambdaadapt = lambda / (hi + 0.0000d1)

    ! computing mixture parameters c
    if((estpi.eq.1).and.(i.gt.1)) then
        if(f0.lt.1.d0) then
            call logitlinkinv((theta / se)**2, p, pi)
            pi = 2.d0 * pi - 1.d0
        end if
    end if

    call gradient(theta, se, lambdaadapt, xtw, res, pi, n, p, grad, alpha)
    f0 = sqrt(sum(grad**2))
    f0_old = f0

    ! computing NR step with step halving search with armijo rule
    call hessian(theta, se0, lambdaadapt, xtx, pi, p, hess, alpha)
    call gradient(theta, se0, lambdaadapt, xtw, res, pi, n, p, grad, alpha)
    call solve(hess, grad, p, info)
    if(info.ne.0) then
        conv = 2
        exit
    end if
    call armijo(theta, se0, grad, f0, alpha, h, X, y, offset, n, p, &
            & lambdaadapt, eta, res, pi, xtw)
    mu = eta
    dev = sum(weights * (res**2))

    ! updating components for variance covariance matrix
    call hessian(theta, se0, lambdaadapt, xtx, pi, p, hess, alpha)
    call inv(p, hess, invH, info)
    if(info.ne.0) then
        theta = theta0
        se = se0
        f0 = f0_old
        conv = 2
        exit
    end if
    call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, tempMat, p)
    do j = 1, p
        hi(j) = tempMat(j,j)
    end do
    call DGEMM('N', 'N', p, p, p, 1.d0, tempMat, p, invH, p, 0.d0, cov1, p)
    edf = sum(hi)
    redf = n - edf
    if(sigma2.lt.0) s2 = dev / redf
    !if(i.eq.1) then
    !    h2 = 0.1d0
    !    cov = cov + h2 * (cov1 - cov)
    !    do j = 1, p
    !        se(j) = sqrt(1.d0 * cov(j,j))
    !    end do
    !    call gradient(theta, se, lambdaadapt, xtw, res, pi, n, p, grad, alpha)
    !    f0 = sqrt(sum(grad**2))
    !else
        call armijo2(theta, se, cov, cov1, s2, f0, alpha, h2, n, p, lambdaadapt, res, pi, xtw)
    !end if

    ! checking possible convergence criterion
    !ind = MAXVAL(abs(se - se0))
    !ind2 = MAXVAL(abs(theta - theta0))
    ind = sqrt(sum((se - se0)**2)) / p
    ind2 = sqrt(sum((theta - theta0)**2)) / p
    ind4 = 0.5d0 * (ind + ind2)
    ind3 = abs(f0 - f0_old)/(0.1d0 + abs(f0))
    if(trace.eq.2) call islasso_trace1_2(tol, i, MAXVAL(lambda), f0, dev, redf, s2, &
        & h, h2, ind, ind2, ind3, ind4)
    if(trace.eq.1) call islasso_trace1_7(tol, i, MAXVAL(lambda), f0, dev, redf, s2, &
        & h, h2, ind3)

    !ind = f0
    if(ind3.le.tol) then
        if(trace.eq.9) call islasso_trace2_6(i)
        if((trace.eq.1).or.(trace.eq.2)) call islasso_trace1_8()
        exit
    end if

    ! conv = 1 if i >= itmax
    if(i.ge.itmax) then
        conv = 1
        exit
    end if
end do

! updating output components
itmax = i
if(sigma2.lt.0) sigma2 = s2
lambda = lambdaadapt
call gradient(theta, se, lambda, xtw, res, pi, n, p, grad2, alpha)

! if standardized beta and se, then return to the original scale
if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    do j = 1, p
        se(j) = sqrt(sigma2 * cov(j,j))
    end do
    lambda = lambda * xse
end if
end subroutine islasso2
