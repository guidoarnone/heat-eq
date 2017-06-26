function z=f6(x,y,t)

r=rem(t,2) #Resto de la division de t por 2.

#Observacion: Si t pertenece al intervalo [2n-2,2n-1] entonces el resto de dividir t por 2 debe ser menor o igual que 1.
#Si t pertenece al intervalo [2n-1,2n] para el mismo n que antes,entonces su resto debe ser mayor o igual que 1, pues si t fuese multiplo de 2 la condicion anterior lo contemplaria (2n-2=0 mod 2).

if x>= (1/4) && x<= (1/2) && y>= (1/4) && y<= (1/2) && r>=0 && r<=1;

  k=1;

else

  k=0;

endif

if x>= (1/2) && x<= (3/4) && y>= (1/2) && y<= (3/4) && r<=1;

 l=1;

else
 l=0;
endif

z=k+l;

endfunction


