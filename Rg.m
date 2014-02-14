function val = Rg(f,g,fs,gs)
    val = (-2*g*((4 - fs^2 - gs^2)*8 + 16*f*fs + 16*g*gs) + 16*gs*(4 + f^2 + g^2))/((4 + f^2 + g^2)^2 * (4 + fs^2 + gs^2));
end
