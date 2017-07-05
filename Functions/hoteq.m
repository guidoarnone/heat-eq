function U = hoteq(alpha, h, dt, f_index)
  
  # Datos iniciales
  N = 1/h + 1;
  M = 1/dt + 1;
  K = N-2; 
  mu = (alpha*dt)/(h^2);

  #  En base a la discretizacion del espacio, generamos un sistema de ecuaciones,
  #donde la temperatura en los puntos de la discretizacion en el tiempo actual son
  #las incognitas, la temperatura en el paso anterior es el dato conocido, y la matriz
  #depende de una relacion mu entre dt, h y alpha.
  #  La matriz resulta en un principio no inversible, pues los puntos de la 
  #discretizacion que se ecuentran en los bordes dan la ecuacion u^n_ij = 0, agregando
  #filas de ceros. Esto se resuelve ignorando los puntos del borde de la region, que no aportan informacion
  #pues su temperatura es siempre nula, y resolviendo el mismo problema para una submatriz
  #que tiene el mismo formato, y es inversible. Como tenemos que pensar la discretizacion del
  #cuadrado como un vector columna de incognitas, armamos las funciones auxiliares coord y to_matrix. La primera es
  #simplemente una biyeccion de {1,..k}x{1,..k} a {1,..,k^2}, y la segunda toma un vector de k^2 coordenadas y lo
  #pasa a una matriz de k x k.
  #  Finalmente, resolvemos el sistema iterativamente tantas veces como sea pedido, y luego pasamos de la submatriz conseguida
  #a la que se nos pide simplemente agregando ceros en los bordes.
  #Los parametros son homónimos a los de la ecuacion original. El último, f_index, es para indicar que funcion f de las propuestas
  #se debe utilizar. En otro archivo se encuentra la funcion g. Segun vimos, la relacion que deben cumplir los coeficientes para
  # que la resolucion sea estable es mu = alpha*dt/(h^2) < 1/4. Por ejemplo, se puede hacer:
  #
  #octave:1> hoteq(1,0.1,0.001,1)
  #
  #En este caso discretizamos el cuadrado en 100 puntos y el tiempo en 1000 instantes. Utilizamos la funcion 1. Por defecto, g esta definida 
  #como g(x,y) = 0, pero se puede utilizar por ejemplo g(x,y) = 10000*sin(pi*x)*y*(1-y) y descomentar mas abajo la impresion de la matriz
  #inicial, para comparar la evolucion de la temperatura.

  U = zeros(N,N);
  H = zeros(K,K); 
  A = zeros(K^2,K^2); 


  #Generamos las coordenadas de la cuadricula y el tiempo
  for i = 1:N
    x(i) = h*(i-1);
    y(i) = h*(i-1);
  end

  for i = 1:M
    t(i) = dt*(i-1);
  end

  #Inicializamos el vector de incognitas, que para el paso inicial es conocido.
  for i = 1:K
    for j = 1:K
     #Corremos en uno las posiciones en el lado derecho de la igualdad, pues las coordenadas
     #(i,j) de H, la cual ahora pensamos en forma de columna via u, son las coordenadas
     #(i+1,j+1) de la matriz original U
     u(coord(i,j,K)) = g(x(i+1),y(j+1));
    end 
  end

 #Llenamos la matriz del sistema en base a los coeficientes
 #anteriormente calculados
 for i = 1:K
    for j = 1:K
        A(coord(i,j,K),coord(i,j,K)) = (-1) -(4*mu);
      if i > 1
        A(coord(i,j,K),coord(i-1,j,K)) = mu;
      end
      if i < K
        A(coord(i,j,K),coord(i+1,j,K)) = mu;
      end
      if j > 1 
        A(coord(i,j,K),coord(i,j-1,K)) = mu;
      end
      if j < K
        A(coord(i,j,K),coord(i,j+1,K)) = mu;
      end
    end
 end

 #Descomentar para imprimir la matriz con la temperatura inicial
 #H = to_matrix(u,K);
 #for i = 2:(N-1)
 #  for j = 2:(N-1)
 #    U(i,j) = H(i-1,j-1);
 #  end
 #end
 #U
 
 #Resolvemos para cada t un sistema lineal, en base a la
 #solucion anterior del sistema basandonos en cierta relacion de
 #recurrencia.
 for s = 1:M
   lu = (-1)*u;
   for i = 1:K
     for j = 1:K
       #Corremos en uno las posiciones en el lado derecho de la igualdad, pues las coordenadas
       #(i,j) de H, la cual ahora pensamos en forma de columna via u, son las coordenadas
       #(i+1,j+1) de la matriz original U
       lu(coord(i,j,K)) -= dt*f(x(i+1),y(j+1),t(s), f_index);
     end 
   end
   u = (A\lu')';
   #Descomentar para imprimir la matriz con la temperatura en cada paso
   #H = to_matrix(u,K);
   #for i = 2:(N-1)
   #  for j = 2:(N-1)
   #    U(i,j) = H(i-1,j-1);
   #  end
   #end 
 end      

 H = to_matrix(u,K);
 for i = 2:(N-1)
   for j = 2:(N-1)
     U(i,j) = H(i-1,j-1);
   end
 end 

endfunction
