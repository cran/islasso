subroutine inv(n, A, invA, info)
integer :: info, n
double precision :: A(n, n), invA(n, n)
integer :: i, j
invA = A
call DPOTRF('U', n, invA, n, info)
if (info.ne.0) call islasso_trace1_3_2(info)
call DPOTRI('U', n, invA, n, info)
if (info.ne.0) call islasso_trace1_3_2(info)
do j = 1, n
    do i = j + 1, n
        invA(i, j) = invA(j, i)
    end do
end do
end subroutine inv

subroutine inv2(n, A, invA, info)
integer :: info,n
double precision :: A(n, n),invA(n, n)
integer :: i,j,ipiv(n)
double precision :: work(n),nwork(n)
invA = A
!call DPOTRF('U', n, invA, n, info)
call DSYTRF('U', n, invA, n, ipiv , work , n, info)
if (info.ne.0) call islasso_trace1_3_2(info)
!call DPOTRI('U', n, invA, n, info)
call DSYTRI('U', n, invA, n, ipiv, nwork, info)
if (info.ne.0) call islasso_trace1_3_2(info)
do j = 1, n
    do i = j + 1, n
        invA(i, j) = invA(j, i)
    end do
end do
end subroutine inv2

subroutine inv3(n, A, invA, info)
integer :: info,n
double precision :: A(n, n),invA(n, n)
integer :: ipiv(n)
integer, parameter :: lwork = 100000
double precision :: work(lwork)
invA = A
call DGETRF(n, n, invA, n, ipiv, info)
if (info.ne.0) call islasso_trace1_3_2(info)
call DGETRI(n, invA, n, ipiv, work, lwork, info)
if (info.ne.0) call islasso_trace1_3_2(info)
!lwork = MIN(lwmax, INT(work(1)))
!call DGETRI(n, invA, n, ipiv, work, lwork, info)
!if (info.ne.0) call islasso_trace1_3_2(info)
!do j = 1, n
!do i = j + 1, n
!invA(i, j) = invA(j, i)
!end do
!end do
end subroutine inv3

subroutine solve(A, b, n, info)
integer :: info, n
double precision :: A(n, n), b(n, 1)
call DPOSV('U', n, 1, A, n, b, n, info)
if (info.ne.0) call islasso_trace1_4_2(info)
end subroutine solve

subroutine solve2(A, b, n, info)
integer :: info, n
double precision :: A(n, n), b(n, 1)
integer :: ipiv(n),lwork
integer, parameter :: lwmax = 100000
double precision :: work(lwmax)
!call DPOSV('U', n, 1, A, n, b, n, info)
lwork = -1
CALL DSYSV('U', n, 1, A, n, ipiv, b, n, work, lwork, info)
lwork = MIN(lwmax, INT(work(1)))
CALL DSYSV('U', n, 1, A, n, ipiv, b, n, work, lwork, info)
if (info.ne.0) call islasso_trace1_4_2(info)
end subroutine solve2

subroutine solve3(A, b, n, info)
integer :: info, n
double precision :: A(n, n), b(n, 1)
integer :: ipiv(n)
call DGESV(n, 1, A, n, ipiv, b, n, info)
if (info.ne.0) call islasso_trace1_4_2(info)
end subroutine solve3

!subroutine crossp(x,xtx,n,p)
!integer :: n,p
!double precision :: x(n,p),xtx(p,p)
!integer :: i,j
!!$omp parallel do shared(p, x, xtx) private(i, j)
!do i = 1, p
!    xtx(i,i) = dot_product(x(:,i),x(:,i))
!    do j = (i+1), p
!        xtx(i,j) = dot_product(x(:,i),x(:,j))
!        xtx(j,i) = xtx(i,j)
!    end do
!end do
!!$omp end parallel do
!end subroutine crossp

subroutine crossp(A, xtx, n, p)
integer :: n, p
double precision :: A(n, p), xtx(p, p)
call DGEMM('T', 'N', p, p, n, 1.d0, A, n, A, n, 0.d0, xtx, p)
end subroutine crossp

subroutine tcrossp(x, txx, n, p)
integer :: n,p
double precision :: x(n,p), txx(n,n)
call DGEMM('N', 'T', n, n, p, 1.d0, x, n, x, n, 0.d0, txx, n)
!do i = 1, n
!    txx(i,i) = dot_product(x(i,:),x(i,:))
!    do j = (i+1), n
!        txx(i,j) = dot_product(x(i,:),x(j,:))
!        txx(j,i) = txx(i,j)
!    end do
!end do
end subroutine tcrossp

!subroutine prod1(x,w,xtw,xtx,n,p)
!integer :: n,p
!double precision :: x(n,p),w(n),xtw(p,n),xtw2(n,p),xtx(p,p)
!integer :: i,j
!do i = 1, p
!    !xtw(i,:) = x(:,i) * w
!    xtw2(:,i) = x(:,i) * w
!    !xtx(i,i) = dot_product(xtw(i,:),x(:,i))
!    xtx(i,i) = dot_product(xtw2(:,i),x(:,i))
!    do j = (i+1), p
!        !xtx(i,j) = dot_product(xtw(i,:),x(:,j))
!        xtx(i,j) = dot_product(xtw2(:,i),x(:,j))
!        xtx(j,i) = xtx(i,j)
!    end do
!end do
!xtw = transpose(xtw2)
!end subroutine prod1

subroutine prod1(x,w,xtw,xtx,n,p)
integer :: n,p
double precision :: x(n,p),w(n),xtw(p,n),xtw2(n,p),xtx(p,p)
integer :: i
do i = 1, p
    xtw2(:,i) = x(:,i) * w
end do
xtw = transpose(xtw2)
call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, x, n, 0.d0, xtx, p)
end subroutine prod1

!subroutine prod2(xtx,tempMat,invH,cov1,hi,p)
!integer :: p
!double precision :: xtx(p,p),tempMat(p,p),invH(p,p),cov1(p,p),hi(p)
!double precision :: tempMat2(p,p),invH2(p,p)
!integer :: i,j,k
!invH2 = transpose(invH)
!do i = 1, p
!    do j = 1, p
!        !tempMat(i,j) = dot_product(invH(i,:),xtx(:,j))
!        tempMat2(j,i) = dot_product(invH2(:,i),xtx(:,j))
!        !tempMat(j,i) = tempMat(i,j)
!    end do
!    hi(i) = tempMat2(i,i)
!    cov1(i,i) = dot_product(tempMat2(:,i),invH(:,i))
!    do k = (i+1), p
!    cov1(i,k) = dot_product(tempMat2(:,i),invH(:,k))
!        cov1(k,i) = cov1(i,k)
!    end do
!end do
!tempMat = transpose(tempMat2)
!end subroutine prod2

subroutine armijo(beta,dir,dtg,f0,alpha,h,x,y,w,offset,n,p,lambda,eta,res,eps)
integer :: n,p
double precision :: beta(p),dir(p),dtg,f0,alpha,x(n,p),y(n),w(n),lambda,eps
double precision :: offset(n),eta(n),res(n),fn,betan(p),ind,c_1,fac,h
h = 1.d0
c_1 = 0.01d0
fac = 0.9d0
eta = 0.d0
betan = beta + h * dir
ind = 0.d0
call linear_predictor(x,betan,eta,offset,n,p)
res = y - eta
fn = sum(w * (res**2)) + lambda * (alpha*sum(abs(betan)) + (1.d0-alpha)*sum(betan**2))
ind = f0 + c_1 * h * dtg
do while (fn.gt.ind)
    h = fac * h
    betan = beta + h * dir
    call linear_predictor(x,betan,eta,offset,n,p)
    res = y - eta
    fn = sum(w * (res**2)) + lambda * (alpha*sum(abs(betan)) + 0.5d0*(1.d0-alpha)*sum(betan**2))
    ind = f0 + c_1 * h * dtg
    if(h.le.eps) exit
end do
beta = betan
f0 = fn
end subroutine armijo


subroutine prod2(xtx,tempMat,invH,cov1,hi,p)
integer :: p
double precision :: xtx(p,p),tempMat(p,p),invH(p,p),cov1(p,p),hi(p)
integer :: i
call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, tempMat, p)
call DGEMM('N', 'N', p, p, p, 1.d0, tempMat, p, invH, p, 0.d0, cov1, p)
do i = 1, p
    hi(i) = tempMat(i,i)
end do
end subroutine prod2


subroutine linear_predictor(x,beta,eta,offset,n,p)
integer :: n,p
double precision :: x(n,p),beta(p),eta(n),offset(n)
integer :: i,j
eta = 0.d0
do i = 1, n
    eta(i) = offset(i)
    do j = 1, p
        if(beta(j).ne.0.d0) eta(i) = eta(i) + x(i,j) * beta(j)
    end do
end do
end subroutine linear_predictor

!subroutine linear_predictor2(x,beta,offset,n,p)
!integer :: n,p
!double precision :: x(n,p),beta(p),offset(n)
!call DGEMV('N', n, p, 1.d0, x, n, beta, 1, 1.d0, offset, 1)
!end subroutine linear_predictor2

subroutine setdiff(p,iprofile,ind_noprofile)
integer :: p,iprofile,ind_noprofile(p-1)
integer :: i,k
k = 0
do i = 1, p
    if(i.ne.iprofile) then
        k = k + 1
        ind_noprofile(k) = i
    end if
end do
end subroutine setdiff
