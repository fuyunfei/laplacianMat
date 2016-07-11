function[L] = cotlapMatrix (vertex,face)  %%  output the cot laplacian Matrix 
n = size(vertex,1)
vertex=vertex';
 face=face';
W = sparse(n,n);
for i=1:3
   i1 = mod(i-1,3)+1;
   i2 = mod(i  ,3)+1;
   i3 = mod(i+1,3)+1;
   pp = vertex(:,face(i2,:)) - vertex(:,face(i1,:));
   qq = vertex(:,face(i3,:)) - vertex(:,face(i1,:));
   % normalize the vectors   
   pp = pp ./ repmat( sqrt(sum(pp.^2,1)), [3 1] );
   qq = qq ./ repmat( sqrt(sum(qq.^2,1)), [3 1] );
   % compute angles
   ang = acos(sum(pp.*qq,1));
   u = cot(ang);
   u = clamp(u, 0.01,100);
   W = W + sparse(face(i2,:),face(i3,:),u,n,n);
   W = W + sparse(face(i3,:),face(i2,:),u,n,n);
end
% Compute the symmetric Laplacian matrix.
d = full( sum(W,1) );
D = spdiags(d(:), 0, n,n);
L = D - W;
end  

