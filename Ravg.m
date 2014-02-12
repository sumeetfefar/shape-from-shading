function val = Ravg(A,i,j)
    val = ( A(i,j-1) + A(i-1,j) + A(i,j+1) + A(i+1,j) )/6 +(A(i-1,j-1)+ A(i-1,j+1)+A(i+1,j-1)+A(i+1,j+1))/12;
end