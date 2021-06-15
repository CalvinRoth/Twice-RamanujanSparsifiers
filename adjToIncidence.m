# Takes as input the adjacency graph for a matrix and return the incidence matrix. 
function [incid, nedges] = adjToIncidence(G)
  n = size(G)(1);
  edgePairs = zeros(0, 2);
  nedges = 0;
  for i = 1:n
    for j = i+1:n
      if (G(i,j) == 1)
        nedges += 1;
        edgePairs(end+1, :) = [i,j];
      endif
    endfor
  endfor
  
  incid = zeros(nedges, n);
  for i = 1:nedges
    incid(i, edgePairs(i, 1)) = 1;
    incid(i, edgePairs(i,2)) = -1;
  endfor
  
endfunction