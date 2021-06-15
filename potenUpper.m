function u = potenUpper(A,u)
 n = size(A)(1);
 B = inv (u*eye(n,n) -A);
 u = trace(B);
 
endfunction