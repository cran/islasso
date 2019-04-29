subroutine islasso2(X,y,n,p,theta,se,cov,lambda,pi,estpi,h,itmax,tol,sigma2,trace,adaptive,offset, &
& conv,stand,intercept,eta,mu,res,dev,weights,hi,edf)
implicit none
integer :: n,p,estpi,itmax,trace,conv,adaptive,stand,intercept
double precision :: X(n,p),y(n),theta(p),se(p),cov(p,p),lambda(p),pi(p),h,tol,sigma2,offset(n)
double precision :: eta(n),mu(n),res(n),dev,weights(n),hi(p),edf
integer :: i,j
double precision :: lambdaadapt(p),se0(p),grad(p),hess(p,p),xtw(p,n),invH(p,p),tempMat(p,p)
double precision :: redf,s2,ind,xm(p),xse(p),xtx(p,p),cov1(p,p)!,diag(p,p)
lambdaadapt = lambda
grad = 0.d0
hess = 0.d0
xtw = 0.d0
invH = 0.d0
tempMat = 0.d0
redf = 0.d0
s2 = 0.d0
ind = 0.d0
xm = 0.d0
xse = 1.d0
xtx = 0.d0
cov1 = 0.d0
hi = 1.d0
conv = 0
!do i = 1, p
!    diag(i,i) = 1.d0
!end do
if(stand.eq.1) call standardize(X,xm,xse,n,p,intercept)
!xtw = TRANSPOSE(X)
!do j = 1, n
!    xtw(:,j) = xtw(:,j) * weights(j)
!end do
!xtx = MATMUL(xtw,X)
call prod1(X,weights,xtw,xtx,n,p)
if(adaptive.eq.1) lambdaadapt = lambda/hi
if(trace.eq.1) call islasso_trace1_1()
do i = 1, itmax
    if(adaptive.eq.1) lambdaadapt = lambda/(hi+0.0000d1)
    se0 = se
    if(estpi.eq.1) then
        pi = (theta/se)**2
        call logitlinkinv(pi, p, pi)
        pi = 2.d0 * pi - 1.d0
        !pi = 1.98d0 * pi - 0.98d0
    end if
    call hessian(theta,se0+0.0000d1,lambdaadapt,xtx,pi,p,hess)
    call gradient(theta,se0+0.0000d1,lambdaadapt,xtw,res,pi,n,p,grad)
    call solve3(hess,grad,p,conv)
    if(conv.ne.0) then
        call islasso_trace1_4(conv)
        exit
    end if
    theta = theta - h * grad
    call hessian(theta,se0,lambdaadapt,xtx,pi,p,hess)
    call inv3(p, hess, invH, conv)
!invH = diag
!call inv2(hess, invH, p, conv)
    if(conv.ne.0) then
        call islasso_trace1_3(conv)
        exit
    end if
    call prod2(xtx,tempMat,invH,cov1,hi,p)
    !tempMat = MATMUL(invH, xtx)
    !do j = 1, p
    !    hi(j) = tempMat(j,j)
    !end do
    !eta = MATMUL(X, theta) + offset
    call linear_predictor(X,theta,eta,offset,n,p)
    mu = eta
    res = y - mu
    dev = sum(res**2)
    if((sigma2).lt.0) then
        edf = sum(hi)
        redf = n - edf
        s2 = dev/redf
    else
        s2 = sigma2
    end if
    !cov = cov + h * (s2 * MATMUL(tempMat, invH) - cov)
    cov = cov + h * (s2 * cov1 - cov)
    do j = 1, p
        se(j) = sqrt(cov(j,j))
    end do
    ind = MAXVAL(abs(se - se0))
    if(trace.eq.1) call islasso_trace1_2(tol, i, ind)
    if(ind.lt.tol) exit
    if(i.ge.itmax) then
        conv = 1
        exit
    end if
end do
itmax = i
edf = sum(hi)
sigma2 = s2
lambda = lambdaadapt
if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    do i = 1, p
        se(i) = sqrt(cov(i,i))
    end do
end if
end subroutine islasso2
