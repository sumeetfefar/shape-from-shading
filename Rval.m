function val = Rval(p,q,s)
    val = (s(1)*p + s(2)*q + 1)/sqrt( (s(2)^2 + s(1)^2 +1) * (p^2 + q^2 +1));
end