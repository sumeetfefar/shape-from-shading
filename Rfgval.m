function val = Rfgval(f,g,s)
    val = ( (4 - f^2 - g^2)*(4 - s(1)^2 - s(2)^2) + 16*f*s(1) + 16*g*s(2) )/( (4 + f^2 + g^2) * (4 + s(1)^2 + s(2)^2) );
end
