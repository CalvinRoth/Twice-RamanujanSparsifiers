n = 25;
G = zeros(n,n);
W = zeros(0,1);
d = 5;
pr = 0.33; %probabality any two nodes are connected


for i =1:n
  for j=i+1:n
    if rand() < pr
      G(i,j) = 1;
      G(j,i) = 1;
      W(end+1,1) = 1;
    endif
  endfor
endfor
W = diag(W);

H = graphSparifier(G,W,d);

tic
% Plotting 
angle = 360/20;
radius = 3;

xy = zeros(n, 2);

for i = 1:n
  xy(i, 1) = radius * cos(angle * i);
  xy(i,2) = radius * sin(angle * i);
endfor
gplot(G, xy)
pause 
m = size(H)(1);
nedges = 0;
for i = 1:m
  for j = i+1:m
    if(H(i,j) != 0)
       nedges += 1;
    endif
  endfor
endfor
nedges
angle = 360/m;
radius = 3;

xy = zeros(m, 2);
for i = 1:m
  xy(i, 1) = radius * cos(angle * i);
  xy(i,2) = radius * sin(angle * i);
endfor

gplot(H, xy)
toc