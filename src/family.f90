subroutine family(fam,link,func,x,n,y)
implicit none
integer :: fam,link,func,n
double precision :: x(n),y(n)
select case(fam)
    case(1)
        select case(link)
            case(1)
                select case(func)
                    case(1)
                        call logitlink(x,n,y)
                    case(2)
                        call logitlinkinv(x,n,y)
                    case(3)
                        call logitmueta(x,n,y)
                    case(4)
                        call binomial_variance(x,n,y)
                end select
            case(2)
                select case(func)
                    case(1)
                        call probitlink(x,n,y)
                    case(2)
                        call probitlinkinv(x,n,y)
                    case(3)
                        call probitmueta(x,n,y)
                    case(4)
                        call binomial_variance(x,n,y)
                end select
        end select
    case(2)
        select case(link)
            case(1)
                select case(func)
                    case(1)
                        call loglink(x,n,y)
                    case(2)
                        call loglinkinv(x,n,y)
                    case(3)
                        call logmueta(x,n,y)
                    case(4)
                        call poisson_variance(x,n,y)
                end select
        end select
    case(3)
        select case(link)
            case(1)
                select case(func)
                    case(1)
                        call inverselink(x,n,y)
                    case(2)
                        call inverselinkinv(x,n,y)
                    case(3)
                        call inversemueta(x,n,y)
                    case(4)
                        call gamma_variance(x,n,y)
                end select
            case(2)
                select case(func)
                    case(1)
                        call loglink(x,n,y)
                    case(2)
                        call loglinkinv(x,n,y)
                    case(3)
                        call logmueta(x,n,y)
                    case(4)
                        call gamma_variance(x,n,y)
                end select
            case(3)
                select case(func)
                    case(1)
                        call identitylink(x,n,y)
                    case(2)
                        call identitylinkinv(x,n,y)
                    case(3)
                        call identitymueta(x,n,y)
                    case(4)
                        call gamma_variance(x,n,y)
                end select
        end select
end select
end subroutine family
