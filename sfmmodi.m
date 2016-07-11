close all ; clear all; 
load('pre.mat');  % combined with flow.mat art.mat face and vertece
 %% save obj <-- v ,f  
saveobjmesh('unproject.obj',unproject(:,:,1),unproject(:,:,2),unproject(:,:,3));
%% input the ply file 
fileID = fopen('1.txt','r');
formatSpec='%f %f %f %f %f %f %f\n';
sizeA = [7 Inf];
ply = fscanf(fileID,formatSpec,sizeA);
ply=ply';
color=ply(:,4:6);
ply=ply(:,1:3);
fclose(fileID);
ply_f=594.512858;  
ply_r= [0.999897 -0.010837 0.009367
        0.010407 0.998942 0.044801
        -0.009843 -0.044699 0.998952];
ply_t= [0.189848 0.067492 -0.992342];
plyrt= [ply_r ply_t'] * [ ply ones(size(ply,1),1)]';
plyrt=plyrt';
 
ux=flow(:,:,1);
uy=flow(:,:,2);

est_A=[-ply_f 0 480;0 -ply_f 480; 0 0 1];
ply_art=est_A*plyrt';
sfm2unpro=zeros(size(plyrt,1),2);

for i=1:size(ply_art,2)
    ply_art(1,i)=ply_art(1,i)/ply_art(3,i); 
    ply_art(2,i)=ply_art(2,i)/ply_art(3,i); 
    x=int16(ply_art(1,i))-ux(int16(ply_art(1,i)))/2;
    y=int16(ply_art(2,i))-uy(int16(ply_art(2,i)))/2;
    sfm2unpro(i,1)=y;
    sfm2unpro(i,2)=x;
end

sfm2unpro=[plyrt sfm2unpro]; %  
 
%%  L for unproj 
[vertex,face]=read_obj_file('unproject.obj');
L=cotlapMatrix (vertex,face);
unproj_rt=[[est_R est_T'] * [vertex ones(size(vertex,1),1)]']';

n=size(vertex,1);
L=sparse(L);
%we want to construct the matrix of the system for L-primes
L_prime = [   L     sparse(n,n)    sparse(n,n)   % the x-part
	           sparse(n,n)     L    sparse(n,n)   % the y-part
               sparse(n,n)     sparse(n,n)     L   ]; % the z-part
%============================for loop start 
 delta = L*vertex;
 Lring = triangulation2adjacency(face,vertex);
 for i= 1:n
  
  disp(i);
  ring = [i find(Lring(i,:))];   % the neighbors  
  V = vertex(ring,:)';
  V = [V
       ones(1,length(ring))];
  A=[];
  % ... Fill C in
  for r=1:length(ring)
    A(r,:) =                [V(1,r)   0           V(3,r)  (-1)*V(2,r) 1  0   0  ];
    A(length(ring)+r,:) =   [V(2,r)  (-1)*V(3,r)    0      V(1,r) 0  1   0  ];
    A(2*length(ring)+r,:) = [V(3,r)  V(2,r) (-1)*V(1,r)     0     0  0   1  ];
  end;
    
  Ainv = pinv(A);
  s =  Ainv(1,:);
  h1 =  Ainv(2,:);
  h2 = Ainv(3,:);
  h3 =  Ainv(4,:);
  
  delta_i = delta(i,:)';
  delta_ix = delta_i(1);
  delta_iy = delta_i(2);
  delta_iz = delta_i(3);
  % T*delta gives us an array of coefficients
  
%   temp1=delta_ix*s + delta_iy*(-1)*h3+delta_iz*h2;
%   temp2=delta_ix*h3 + delta_iy*s+ delta_iz*(-1)*h1;
%   temp3=delta_ix*(-1)*h2+ delta_iy*h1 +delta_iz*s; 
%    Tdelta =  [temp1 
%              temp2
%              temp3];
    Delta_forT=[delta_ix 0 delta_iz -delta_iy 
     delta_iy -delta_iz 0 delta_ix
     delta_iz delta_iy -delta_ix 0]; 
   Tdelta= Delta_forT*[s;h1; h2; h3]; 

  % updating the weights in Lx_prime, Ly_prime, Lw_prime
  L_prime(i,[ring (ring + n) (ring+2*n)]) = L_prime(i,[ring (ring + n) (ring+2*n)]) +(-1)*Tdelta(1,:);                                  
  L_prime(i+n,[ring (ring + n) (ring+2*n)]) = L_prime(i+n,[ring (ring + n) (ring+2*n)]) +(-1)*Tdelta(2,:)  ;      
  L_prime(i+2*n,[ring (ring + n) (ring+2*n)]) = L_prime(i+2*n,[ring (ring + n) (ring+2*n)])+(-1)*Tdelta(3,:);

end; 

w= 1;
static_anchors = [1,20,50,200]; % sine_circle
handle_anchors = [200];             % sine_circle
% building the least-squares system matrix
A_prime = L_prime;
rhs = zeros(3*n,1);
anch_pos = [];

anchors = [static_anchors handle_anchors];
for j=1:length(anchors)
    A_prime(anchors(j),:) = 0 ;     
    A_prime(anchors(j),anchors(j)) = w*1;
       A_prime(anchors(j)+n,:) = 0 ;  
       A_prime(anchors(j)+n,anchors(j)+n) =  w*1 ;
          A_prime(anchors(j)+2*n,:) = 0 ;  
          A_prime(anchors(j)+2*n, anchors(j)+2*n) =  w*1 ; 
    rhs(anchors(j))=w*vertex(anchors(j),1);      
    rhs(anchors(j)+n)=w*vertex(anchors(j),2);
    rhs(anchors(j)+2*n)=w*vertex(anchors(j),3);
%%    
%   A_prime = [A_prime
% 	     w*((1:(3*n))==anchors(j))
% 	     w*((1:(3*n))==(anchors(j)+n));
%          w*((1:(3*n))==(anchors(j)+2*n))];
%   rhs = [rhs
% 	 w*vertex(anchors(j),1)
% 	 w*vertex(anchors(j),2)
%      w*vertex(anchors(j),3)];
%      
%   anch_pos = [anch_pos
% 	      vertex(anchors(j),1:3)];
end;

rhs(handle_anchors)=rhs(handle_anchors)+1*w;
rhs(handle_anchors+n)=rhs(handle_anchors+n)+1*w;
rhs(handle_anchors+2*n)=rhs(handle_anchors+2*n)+1*w;

%rhs((lenr-2):lenr) = vertex(handle_anchors,:)'+[5 5 0]'*w; %[x_input y_input z_input]'

% solving for v-primes
xyz_col =  A_prime\rhs;
vertex1 = [xyz_col(1:n) xyz_col((n+1):(2*n)) xyz_col((2*n+1):(3*n))];
write_obj('xyz1.obj',vertex1,face,0);

%% show  anchor and mesh 
axis equal; 
trisurf(face,vertex1(:,1),vertex1(:,2),vertex1(:,3),'FaceColor',[0.26,0.33,1.0 ]);
light('Position',[-1.0,-1.0,100.0],'Style','infinite');
lighting phong; hold on; 
index = anchors ;
for index_i = 1: length(index)
 i=index(index_i);
scatter3(vertex(i,1), vertex(i,2), vertex(i,3) ,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 0 .75]);  hold on; 
scatter3(vertex1(i,1), vertex1(i,2), vertex1(i,3) ,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75]);  hold on; 
end 
scatter3(vertex1(i,1), vertex1(i,2), vertex1(i,3) ,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 0 0]);  hold on; 
    