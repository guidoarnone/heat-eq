function A = hoteq(alpha, h, dt, f_index)
  
  # Datos iniciales
  N = 1/h;
  M = 1/dt;
  mu = alpha*N^2;
  eta = M + 4*mu;

  #Inicializamos A
  A = zeros(N^2,N^2);

  #Generamos las coordenadas de la cuadricula y el tiempo
  for i = 1:N
    x(i) = h*i;
    y(i) = h*i;
  end

  for i = 1:M
    t(i) = dt*i;
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
 end  	

endfunction
