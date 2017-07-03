function U = to_matrix(u, N)
 #Pasamos de la columna de incognitas a la forma matricial
 for i = 1:N
   for j = 1:N
     U(i,j) = u(coord(i,j,N));
   end
  end 
endfunction