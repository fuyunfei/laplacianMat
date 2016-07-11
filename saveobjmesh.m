function saveobjmesh(name,x,y,z,nx,ny,nz)
% SAVEOBJMESH Save a x,y,z mesh as a Wavefront/Alias Obj file
% SAVEOBJMESH(fname,x,y,z,nx,ny,nz)
%     Saves the mesh to the file named in the string fname
%     x,y,z are equally sized matrices with coordinates.
%     nx,ny,nz are normal directions (optional)

  
  normals=1;
  if (nargin<5) normals=0; end
  l=size(x,1); h=size(x,2);  

  n=zeros(l,h);
  
  fid=fopen(name,'w');
  nn=1;
  for i=1:l
    for j=1:h
      if(x(i,j)~=0)
      n(i,j)=nn; 
      fprintf(fid, 'v %f %f %f\n',x(i,j),y(i,j),z(i,j));
      nn=nn+1;
      end
    end
  end
  fprintf(fid,'g mesh\n');
  
  for i=1:(l-1)
    for j=1:(h-1)
      if (x(i,j)*x(i+1,j)*x(i+1,j+1)~=0) 	
          fprintf(fid,'f %d %d %d \n',n(i,j),n(i+1,j),n(i+1,j+1));
      end
          if (x(i,j)*x(i+1,j+1)*x(i,j+1)~=0)
          fprintf(fid,'f %d %d %d \n',n(i,j),n(i+1,j+1),n(i,j+1));
          end
    end
  end
  fprintf(fid,'g\n\n');
  fclose(fid);
  
