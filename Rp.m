function val = Rp(p,q,ps,qs)
    val = (ps*(q*q +1) - p*(qs*q + 1))/(sqrt((qs^2 + ps^2 +1) * (p^2 + q^2 +1)) * (p^2+q^2+1));
end