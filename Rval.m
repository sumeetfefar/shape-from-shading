function val = Rval(p,q,s)
    val = (s(1)*p + s(2)*q + 1)/sqrt( (s(2)^2 + s(1)^2 +1) * (p^2 + q^2 +1));
end

function val = Rp(p,q,s)
    val = s(1)*(q^2 +1) - p*(s(2)*q + 1)/(sqrt((s(2)^2 + s(1)^2 +1) * (p^2 + q^2 +1)) * (p^2+q^2+1));
end

function val = Rq(p,q,s)
    val = s(2)*(p^2 +1) - q*(s(1)*p + 1)/(sqrt((s(2)^2 + s(1)^2 +1) * (p^2 + q^2 +1)) * (p^2+q^2+1));
end