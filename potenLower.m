function l = potenLower(A,l)
 n = size(A)(1);
 B = inv (A - l*eye(n,n));
 l = trace(B);
 
endfunction