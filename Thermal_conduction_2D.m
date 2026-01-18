% Askisi 3.2
clear; clc;
nx = 81;
ny = 81;
n = nx*ny;
nex=nx-1;
ney=ny-1;
nel = nex*ney;

%Define matrices
T = zeros(n,1);
s = zeros(n,1);
a = zeros(n,n);
b = zeros(n,1);
z = zeros(nel,4);
iflagdir = ones(n,1);

ml = [4 2 2 1; 2 4 1 2; 2 1 4 2; 1 2 2 4];
llx = [2 1 -2 -1; 1 2 -1 -2; -2 -1 2 1; -1 -2 1 2];
lly = [2 -2 1 -1; -2 2 -1 1; 1 -1 2 -2;-1 1 -2 2];

%Parameters
Lx=1;
Ly=1;

dx = Lx/nex;
dy = Ly/ney;

%Define heat source
for i=1:nx
    for j=1:ny
        x=(i-1)*dx;
        y=(j-1)*dy;
        k=j+(i-1)*ny;
        s(k)=sin(pi*x)*sin(pi*y);
    end
end

%Form local element matrices
ll=dy/dx/6*llx+dx/dy/6*lly;
ml=ml*dx*dy/36;

%Form connectivity matrix
for ie=1:nel
    z(ie,1)=floor((ie-1)/ney)*ny+mod(ie-1,ney)+1;
    z(ie,2)=z(ie,1)+1;
    z(ie,3)=z(ie,1)+ny;
    z(ie,4)=z(ie,3)+1;
end

%Dirichlet boundary conditions
%left and right line
for j=1:ny
    iflagdir(j)=0;
    iflagdir(j+(nx-1)*ny)=0;
end

%bottom and top line
for i=1:nx
    iflagdir(1+(i-1)*ny)=0;
    iflagdir(ny+(i-1)*ny)=0;
end

%Form linear system
%Local element matrix assembly
for ie=1:nel
    for i=1:4
        i1=z(ie,i);
        for j=1:4
            j1=z(ie,j);
            a(i1,j1)=a(i1,j1)+ll(i,j);
            b(i1)=b(i1)+ml(i,j)*s(j1);
        end
    end
end

%Dirichlet  conditions
for i=1:n
    if iflagdir(i) == 0
       a(i,:)=0;
       a(i,i)=1;
       b(i)=0;
    end
end

% Solve linear system A.T = b
T = a\b;

%Store part of the solution
fid=fopen('out4.txt','w');
for i=1:nx
    x=dx*(i-1);
    k=(ny-1)/2+1+(i-1)*ny;
    fprintf(fid,'%10.6f %10.6f \n', x,T(k));
end
fclose(fid);
G=load("out4.txt");
%