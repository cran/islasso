subroutine islasso2(X,y,n,p,theta,se,cov,lambda,alpha,pi,estpi,h,itmax,tol,sigma2,trace,adaptive, &
& offset,conv,stand,intercept,eta,mu,res,dev,weights,hi,edf,grad2)
implicit none
integer :: n,p,estpi,itmax,trace,conv,adaptive,stand,intercept
double precision :: X(n,p),y(n),theta(p),se(p),cov(p,p),lambda(p),alpha,pi(p),h,tol,sigma2
double precision :: offset(n),eta(n),mu(n),res(n),dev,weights(n),hi(p),edf,grad2(p)
! internal variables
integer :: i,j,conta
double precision :: xm(p),xse(p),f0,xtw(p,n),xtx(p,p),lambdaadapt(p),thetaold(p),se0(p),grad(p)
double precision :: hess(p,p),dtg,hh,invH(p,p),tempMat(p,p),cov1(p,p),redf,s2,ind,ind2

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
conta = 0

! standardizing
if(stand.eq.1) call standardize(X,xm,xse,n,p,intercept)

! cmputing starting values
call linear_predictor(X,theta,eta,offset,n,p)
mu = eta
res = y - mu
f0 = sum(weights * (res**2)) + MAXVAL(lambda) * (alpha*sum(abs(theta)) + &
    & 0.5d0 * (1.d0-alpha)*sum(theta**2))
call prod1(X,weights,xtw,xtx,n,p)
if(adaptive.eq.1) lambdaadapt = lambda/hi
!if(trace.eq.1) call islasso_trace1_1()

! main body
do i = 1, itmax
    thetaold = theta
    if(adaptive.eq.1) lambdaadapt = lambda/(hi+0.0000d1)
    se0 = se

    ! computing mixture parameters c
    if((estpi.eq.1)) then
        pi = (theta/se)**2
        call logitlinkinv(pi, p, pi)
        pi = 2.d0 * pi - 1.d0
    end if

    ! computing NR step with step halving search with armijo rule
    call hessian(theta,se0,lambdaadapt,xtx,pi,p,hess,alpha)
    call gradient(theta,se0,lambdaadapt,xtw,res,pi,n,p,grad,alpha)
    grad2 = grad

    call solve3(-hess,grad,p,conv)
    if(conv.ne.0) then
        call islasso_trace1_4(conv)
        exit
    end if
    dtg = dot_product(grad, grad2)
    hh = h
    call armijo(theta,grad,dtg,f0,alpha,hh,X,y,weights,offset,n,p,MAXVAL(lambda),eta,res,tol)

    ! updating components for variance covariance matrix
    call hessian(theta,se0,lambdaadapt,xtx,pi,p,hess,alpha)
    call inv3(p, hess, invH, conv)
    if(conv.ne.0) then
        call islasso_trace1_3(conv)
        exit
    end if
    call prod2(xtx,tempMat,invH,cov1,hi,p)
    mu = eta
    dev = sum(weights * (res**2))
    if((sigma2).lt.0) then
        edf = sum(hi)
        redf = n - edf
        s2 = dev/redf
    else
        s2 = sigma2
    end if
    cov = cov + 0.1d0 * (s2 * cov1 - cov)
    do j = 1, p
        se(j) = sqrt(cov(j,j))
    end do

    ! checking possible convergence criterion
    ! conv = 0 if abs(se - se0) < tol
    ! conv = 1 if i >=itmax
    ! conv = 2 if error in step halving (h < tol) and abs(theta - thetaold) < tol
    ind = MAXVAL(abs(se - se0))
    ind2 = MAXVAL(abs(theta - thetaold))
    if(trace.eq.1) call islasso_trace1_7(tol, i, lambda(2), dev, redf, s2, hh, ind, ind2)
    if(trace.ge.2) call islasso_trace1_2(tol, i, lambda(2), dev, redf, s2, hh, ind, ind2)
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

! updating output components
itmax = i
edf = sum(hi)
sigma2 = s2
lambda = lambdaadapt
h = hh
call gradient(theta,se,lambdaadapt,xtw,res,pi,n,p,grad2,alpha)

! if standardized beta and se are returned to original scale
if(stand.eq.1) then
    call check_out(theta,cov,xm,xse,p,intercept)
    do i = 1, p
        se(i) = sqrt(cov(i,i))
    end do
end if
end subroutine islasso2
