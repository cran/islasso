subroutine islasso_glm(X, y, n, p, theta, se, cov, lambda, alpha, pi, estpi, &
& h, itmax, tol, sigma2, trace, adaptive, offset, conv, stand, intercept, eta, &
& mu, dev, weights, hi, edf, fam, link, grad2)

implicit none

! external variables
integer :: n, p, estpi, itmax, trace, adaptive, conv, stand, intercept, fam, link
double precision :: X(n,p), y(n), theta(p), se(p), cov(p,p),lambda(p), alpha
double precision :: pi(p), h, tol, sigma2, offset(n), eta(n), mu(n), dev
double precision :: weights(n), hi(p), edf, grad2(p)

! internal variables
integer :: i, j, itmax2, conv2, trace2 ! linkinv, mueta, variance, trace2, stand2,
double precision :: xm(p), xse(p), theta0(p), varmu(n), ind3, dev0, pi2(p), ind4 !, eta2(n), mu2(n)
double precision :: mu_eta_val(n), res(n), z(n), w(n), se0(p), cov0(p,p), s20
double precision :: lambda2(p), s2, hh, offset2(n), ind, ind2, f0, grad0(p), edf0

hi = 1.d0

!variance = 4
!mueta = 3
!linkinv = 2
!trace2 = 0
!stand2 = 0

xm = 0.d0
xse = 1.d0
dev0 = 1.d0

eta = eta + offset
call family(fam, link, 2, eta, n, mu)
!    mu2 = mu
call family(fam, link, 4, mu, n, varmu)
!    eta2 = eta
call family(fam, link, 3, eta, n, mu_eta_val)
res = (y - mu) / mu_eta_val
z = (eta - offset) + res
w = weights * (mu_eta_val**2) / varmu

offset2 = 0.d0
s2 = sigma2
lambda2 = lambda

if(stand.eq.1) then
    call standardize(X, xm, xse, n, p, intercept)
    lambda2 = lambda2 / xse
end if

trace2 = 0
if((trace.eq.1).or.(trace.eq.3)) then
    if(trace.eq.3) trace2 = 9
    call islasso_trace2_5(tol, MAXVAL(lambda))
end if
do i = 1, itmax
    call rchkusr()

    theta0 = theta
    se0 = se
    cov0 = cov
    grad0 = grad2
    edf0 = edf
    s20 = s2

    if((estpi.eq.1).and.(i.gt.1)) then
        call logitlinkinv((theta / se)**2, p, pi)
        pi = 2.d0 * pi - 1.d0
    end if
    pi2 = pi
    if(i.eq.1) pi2 = 1.d0

    itmax2 = itmax
    s2 = sigma2
    hh = h
    conv2 = conv

    call islasso2(X, z, n, p, theta, se, cov, lambda2, alpha, pi2, 0, hh, &
        & itmax2, tol, s2, trace2, adaptive, offset2, conv2, 0, intercept, &
        & eta, mu, res, dev, w, hi, edf, grad2)

    if((hh.lt.0.000001d0).or.((conv2).eq.2)) then
        theta = theta0
        se = se0
        cov = cov0
        grad2 = grad0
        edf = edf0
        s2 = s20
        conv = 3
        if(trace.eq.2) call islasso_trace2_2(tol, i, itmax2, MAXVAL(lambda), f0, dev, n - edf, s2, ind, ind2, ind3, ind4)
        if((trace.eq.1).or.(trace.eq.3)) call islasso_trace2_7(i, itmax2, dev, n - edf, s2, ind3, ind4)
        exit
    end if
    
!    eta = eta + offset
!    eta2 = eta
    call family(fam, link, 2, eta, n, mu)
!    mu2 = mu
    call family(fam, link, 4, mu, n, varmu)
!    eta2 = eta
    call family(fam, link, 3, eta, n, mu_eta_val)
    res = (y - mu) / mu_eta_val
    z = (eta - offset) + res
    w = weights * (mu_eta_val**2) / varmu

    f0 = sum(grad2**2)

    !ind = MAXVAL(abs(se - se0))
    !ind2 = MAXVAL(abs(theta - theta0))
    ind = sqrt(sum((se - se0)**2)) / p
    ind2 = sqrt(sum((theta - theta0)**2)) / p
    ind4 = 0.5d0 * (ind + ind2)
    ind3 = abs(ind4 - dev0) / (0.1d0 + abs(ind4)) !0.5d0 * (ind + ind2)
    if(trace.eq.2) call islasso_trace2_2(tol, i, itmax2, MAXVAL(lambda), f0, dev, n - edf, s2, ind, ind2, ind3, ind4)
    if((trace.eq.1).or.(trace.eq.3)) call islasso_trace2_7(i, itmax2, dev, n - edf, s2, ind3, ind4)
    
    ! conv = 0 if abs(se - se0) < tol
    if(ind3.lt.tol) then
        if(trace.ge.1) call islasso_trace1_8()
        exit
    end if
    dev0 = ind4
    ! conv = 1 if i >= itmax
    if(i.ge.itmax) then
        conv = 1
        exit
    end if
end do

itmax = i
!lambda = lambda2
sigma2 = s2
pi = pi2

if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    do j = 1, p
        se(j) = sqrt(sigma2 * cov(j,j))
    end do
end if
end subroutine islasso_glm
