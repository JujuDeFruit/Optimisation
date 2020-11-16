function [callback] = MoindresCarres(a, b, x, y)
  % Initialisation
  callback = 0;
  N = size(x, 1);
  % On effectue la méthode des moindres carrées pour chaque valeur de a et b sur
  % les éléments x et y. 
  for i=1:N
      callback = callback + (a*x(i) + b - y(i))^2;
  endfor
endfunction