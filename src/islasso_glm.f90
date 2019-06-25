subroutine islasso_glm(X,y,n,p,theta,se,cov,lambda,alpha,pi,estpi,h,itmax,tol,sigma2,trace,adaptive, &
& offset,conv,stand,intercept,eta,mu,dev,weights,hi,edf,fam,link,grad2)
implicit none
integer :: n,p,estpi,itmax,trace,adaptive,conv,stand,intercept,fam,link
double precision :: X(n,p),y(n),theta(p),se(p),cov(p,p),lambda(p),alpha,pi(p),h,tol,sigma2
double precision :: offset(n),eta(n),mu(n),dev,weights(n),hi(p),edf,grad2(p)
! internal variables
integer :: linkinv,mueta,variance,trace2,stand2,i,itmax2,conta
double precision :: xm(p),xse(p),thetaold(p),mu2(n),varmu(n),eta2(n),mu_eta_val(n),res(n),z(n),w(n)
double precision :: se0(p),pi2(p),lambda2(p),s2,hh,offset2(n),ind,ind2

variance = 4
mueta = 3
linkinv = 2
trace2 = 0
offset2 = 0.d0
stand2 = 0
ind = 0.d0
xm = 0.d0
xse = 1.d0
varmu = 0.d0
mu_eta_val = 0.d0
z = 0.d0
w = 0.d0
hi = 1.d0
conta = 0

if(stand.eq.1) call standardize(X,xm,xse,n,p,intercept)

do i = 1, itmax
    thetaold = theta
    mu2 = mu
    call family(fam,link,variance,mu2,n,varmu)
    eta2 = eta
    call family(fam,link,mueta,eta2,n,mu_eta_val)
    res = (y - mu) / mu_eta_val
    z = (eta - offset) + res
    w = (weights * (mu_eta_val**2) / varmu)

    se0 = se
    itmax2 = itmax / 2
    pi2 = pi
    lambda2 = lambda
    s2 = sigma2
    hh = h

    call islasso2(X,z,n,p,theta,se,cov,lambda2,alpha,pi2,estpi,hh,itmax2,tol,s2,trace2, &
& adaptive,offset2,conv,stand2,intercept,eta,mu,res,dev,w,hi,edf,grad2)

    eta = eta + offset
    eta2 = eta
    call family(fam,link,linkinv,eta2,n,mu)

    ind = MAXVAL(abs(se - se0))
    ind2 = MAXVAL(abs(theta - thetaold))
    if(trace.eq.1) call islasso_trace1_7(tol, i, lambda(2), dev, n - edf, s2, hh, ind, ind2)
    if(trace.ge.2) call islasso_trace1_2(tol, i, lambda(2), dev, n - edf, s2, hh, ind, ind2)
    if(ind.lt.tol) then
        if(trace.ge.1) call islasso_trace1_8()
        exit
    end if
    if(i.ge.itmax) then
        conv = 1
        exit
    end if
    if((hh.le.tol).and.(ind2.lt.tol)) then
        conta = conta + 1
        if(conta.ge.10) then
            if(trace.ge.1) call islasso_trace1_8()
            conv = 2
            exit
        end if
    else
        conta = 0
    end if
end do

itmax = i
pi = pi2
lambda = lambda2
cov = cov / s2

mu2 = mu
call family(fam,link,variance,mu2,n,varmu)
dev = sum(weights * ((y - mu)**2) / varmu)
if((sigma2).lt.0) sigma2 = dev / (n - edf)

!eta2 = eta
!call family(fam,link,mueta,eta2,n,mu_eta_val)
!w = (weights * (mu_eta_val**2)) / varmu
!call prod1(X,w,xtw,xtx,n,p)
!call hessian(theta,se,lambda,xtx,pi,p,hess)
!call inv2(p,hess,invH,conv)
!call prod2(xtx,tempMat,invH,cov,hi,p)
!edf = sum(hi)

if(stand.eq.1) call check_out(theta,cov,xm,xse,p,intercept)
cov = sigma2 * cov
do i = 1, p
    se(i) = sqrt(cov(i,i))
end do
end subroutine islasso_glm
