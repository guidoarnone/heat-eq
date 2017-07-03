function U = hoteq(alpha, h, dt, f_index)
  
  # Datos iniciales
  N = 1/h ;
  M = 1/dt;
  mu = alpha*N^2;
  eta = M + 4*mu;

  #Inicializamos la matriz del sistema y la grilla
  A = zeros(N^2,N^2);
  U = zeros(N,N);


  #Generamos las coordenadas de la cuadricula y el tiempo
  for i = 1:N
    x(i) = h*(i-1);
    y(i) = h*(i-1);
  end

  for i = 1:M
    t(i) = dt*(i-1);
  end

  for i = 1:N
    for j = 1:N
     u(coord(i,j, N)) = g(x(i),y(j));
    end 
  end

 #Resolvemos para cada t un sistema lineal, en base a la
 #matriz anterior, basandonos en cierta relacion de
 #recurrencia.
 for k = 1:M
   lu = (-1*M)*u;
   for i = 1:N
     for j = 1:N
       lu(coord(i,j, N)) -= f(x(i),y(j),t(k), f_index);
     end 
   end
   for i = 1:N
     for j = 1:N
       for l = 1:(N^2)
           A(coord(i,j,N), l) = 0; 
       end
       if i > 1 && i < N && j > 1 && j < N
         A(coord(i,j,N),coord(i,j,N)) = eta;
         A(coord(i,j,N),coord(i,j+1,N)) = mu;
         A(coord(i,j,N),coord(i,j-1,N)) = mu;
         A(coord(i,j,N),coord(i+1,j,N)) = mu;
         A(coord(i,j,N),coord(i-1,j,N)) = mu;
       end
     end
   end
   u = (A\lu')';
   U = to_matrix(u,N)
 end      
  
endfunction
