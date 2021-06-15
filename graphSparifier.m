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
  
  epL - epU
  l = -n/epL;
  u = n/epU; 
  
  [B, m] =  adjToIncidence(G);
  
  A  = zeros(n,n);
  
  barrierU = epU;
  barrierL = epL;
  
  vert = buildV(G, B, W, n, m);
  s = zeros(m,m);
  
  
  for i = 1:d*n
    for j = 1:size(vert)(1)
      v = vert(:, j);
      augT1 = inv(((u+deltU)*eye(n,n) - A));
      augT2 = augT1 * augT1;
      UA = v' * augT2 * v;
      UA = UA / ( potenUpper(A, u) - potenUpper(A, u + deltU));
      UA += v' * augT1 * v;
      
      augT1 = inv( A - (l + deltU) * eye(n,n));
      augT2 = augT1 * augT1;
      LA = v' * augT2 * v;
      LA = LA / (potenLower(A, l + deltL) - potenLower(A, l));
      LA -= v' * augT1 * v;
      
      if UA < LA
        [i,j];
        midpoint = (UA + LA)/2;
        s(j,j) += 1/midpoint;
        A = A + (1/midpoint) * v * v';
        l += deltL;
        u += deltU;
        break;
      endif
    endfor
  endfor
  
  Lh = B' * sqrtm(W)* s * sqrtm(W) * B;
  Lh( ~any(Lh,2), : ) = [];  %rows
  Lh( :, ~any(Lh,1) ) = [];  %columns
  
  H = Lh
  for i = 1:size(H)(1)
    for j = 1:size(H)(2)
      if(H(i,j) != 0)
        H(i,j) = 1;
      endif 
    endfor
  endfor
 
  
endfunction