function u = potenUpper(A,u)
 n = size(A);
 B = inv (u*eye(n) -A);
 u = trace(B);
endfunction