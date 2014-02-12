function val = Rval(p,q,s)
    val = (s(1)*p + s(2)*q + 1)/sqrt( (s(2)*s(2) + s(1)*s(1) +1) * (p*p + q*q +1));
end
