function u = hoteq(alpha, h, dt, f_index)
  
  # Datos iniciales
  N = 1/h;
  M = 1/dt;
  mu = alpha*N^2;
  eta = M + 4*mu;
   
  #Inicializamos A
  A = zeros((N+1)^2,(N+1)^2);

  #Generamos las coordenadas de la cuadricula y el tiempo
  for i = 1:N+1 #x e y empiezan desde 0
    x(i) = h*(i-1);
    y(i) = h*(i-1);
  end

  for i = 1:M+1     #t empieza desde 0
    t(i) = dt*(i-1);
  end
  
  for i = 1:N+1
	  for j = 1:N+1
    
         lu(coord(i,j,N)) = g(x(i),y(j)); 
    end
  end
  

 #Resolvemos para cada t un sistema lineal, en base a la
 #matriz anterior, basandonos en cierta relacion de
 #recurrencia.
 
 for k = 1:M+1
   u = (-1*M)*lu;
   for i = 1:N+1
  	 for j = 1:N+1
       lu(coord(i,j,N)) -= f(x(i),y(j),t(k), f_index);
     end 
   end
 end
 
 
   for i = 1:N+1
     for j = 1:N+1
       for k = 1:((N+1)^2)
           A(coord(i,j,N), k) = 0;
       
       end
         
       if i > 1 && i < N+1 && j > 1 && j < N+1
                  
           A(coord(i,j,N),coord(i,j,N)) = eta;
           A(coord(i,j,N),coord(i,j+1,N)) = mu;
           A(coord(i,j,N),coord(i,j-1,N)) = mu;
           A(coord(i,j,N),coord(i+1,j,N)) = mu;
           A(coord(i,j,N),coord(i-1,j,N)) = mu;
         
           
       end
     end
   end
   
 u = (A\lu')';
 
 end  	 	

endfunction
