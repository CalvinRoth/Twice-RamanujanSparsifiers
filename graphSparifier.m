% program to find the sparsification of a graph
% Takes as input G and w where
%% G := G is the adjacency matrix of input graph 
%% w := is the weight matrix for the edges
%% d := there will be at most d(n-1) edges in ouput graph
% returns H the laplicain of a sparsification of graph 
function H = graphSparifier(G, W, d)
    n = size(G)(1);
    
  deltL = 1;
  deltU = (sqrt(d)+1)/(sqrt(d)-1);
  
  epL = 1/sqrt(d);
  epU = (sqrt(d)-1)/(d + sqrt(d));
  
  l = -n/epL;
  u = n/epU; 
  
  [B, m] =  adjToIncidence(G);
  
  A  = zeros(n,n);
  
  barrierU = epU;
  barrierL = epL;
  
  vert = buildV(G, B, W, n, m);
  s = zeros(m,m);
  tic
  for i = 1:d*n
    upperT1 = inv(((u+deltU)*eye(n,n) - A));
    upperT2 = upperT1* upperT1;
    lowerT1 = inv( A - (l + deltU) * eye(n,n));
    lowerT2 = lowerT1 ** 2;
    for j = 1:size(vert)(2)
      v = vert(:, j);

      UA = v' * upperT2 * v;
      UA = UA / ( potenUpper(A, u) - potenUpper(A, u + deltU));
      UA += v' * upperT1 * v;
      UA = real(UA);
      LA = v' * lowerT2 * v;
      LA = LA / (potenLower(A, l + deltL) - potenLower(A, l));
      LA -= v' * lowerT1 * v;
      LA = real(LA);
      
      if UA < LA
        [i]
        midpoint = (UA + LA)/2;
        s(j,j) += 1/midpoint;
        A = A + (1/midpoint) * v * v';
        l += deltL;
        u += deltU;
        break;
      endif
    endfor

endfor
toc
kappa = (d+1+2*sqrt(d))/(d+1-2*sqrt(d))
eig(vert*s*vert' - vert*vert')
eigs(kappa*vert*vert'-vert*s*vert')
"Done with main loop"
  tic
  Lh = B' * sqrtm(W)* s * sqrtm(W) * B;
  
  H = Lh
  for i = 1:size(H)(1)
    for j = 1:size(H)(2)
      if(H(i,j) != 0)
        H(i,j) = 1;
      endif 
    endfor
  endfor
 
  toc
endfunction