function val = Rq(p,q,ps,qs)
    val = (qs*(p^2 +1) - q*(ps*p + 1))/(sqrt((qs^2 + ps^2 +1) * (p^2 + q^2 +1)) * (p^2+q^2+1));
end