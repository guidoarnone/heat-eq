function u=ADI_Murrone(u_cero,nx,nt)
tic
%Creamos las grillas. Aunque la grilla_x sea igual a la grilla_y, las
%creamos por separado
grilla_t = nt:nt:0.5-nt;
grilla_x = nx:nx:2-nx;
grilla_y = nx:nx:2-nx;

[X,Y] = meshgrid(grilla_x,grilla_y);

%Creamos las matrices
A = diag(4*ones(1,length(grilla_x)))+diag(-(2-5*nx)*ones(1,length(grilla_x)-1),-1)+diag(-(2+5*nx)*ones(1,length(grilla_x)-1),1);
A = sparse(A);
B = diag(-4*ones(1,length(grilla_x)))+diag((2-5*nx)*ones(1,length(grilla_x)-1),1)+diag((2+5*nx)*ones(1,length(grilla_x)-1),-1);
B = sparse(B);
u = zeros(length(grilla_x),length(grilla_y)); 

%valor inicial

for i=1:length(grilla_x)
    for j=1:length(grilla_y)
      u(i,j) = u_cero(grilla_x(i),grilla_y(j));
    end
end

I = eye(size(A));
I = sparse(I);
f = ones(length(grilla_x),length(grilla_x)); % la fuente
r=nt/(nx^2);


%empezamos la iteracion del metodo


for i=1:length(grilla_t)+1
   u_estrella = (I+0.25*r*A)\(u*((I+0.25*r*B))+(nt/2)*f);
   u = ((I+0.25*r*A')'\(((I+0.25*r*B')*u_estrella+(nt/2)*f))')' ;
   
   %comandos graficadores
   %surf(X,Y,u);
   %zlim([0,1])
   %shading interp;
   %view(0,90);
   %pause(0.1)
end

toc

end

