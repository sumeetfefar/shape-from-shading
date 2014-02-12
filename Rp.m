function val = Rp(p,q,ps,qs)
    val = (ps*(q*q +1) - p*(qs*q + 1))/(sqrt((qs*qs + ps*ps +1) * (p*p + q*q +1)) * (p*p+q*q+1));
end