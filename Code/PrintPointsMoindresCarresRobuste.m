function [] = PrintPointsMoindresCarresRobuste(robuste, x, y, space, dimension, sigma)
  
  a = linspace(space(1), space(2), dimension);
  b = a;
  
  moindresCarres = ones(dimension);

  for i=1:dimension
    for j=1:dimension
      moindresCarres(i,j) = MoindresCarres(a(i), b(j), x, y);
    endfor
  endfor
  
  mC = MinMatrix(moindresCarres, space, dimension);
  
  if sigma != 1
    rob = QuasiNewton(10^-4, 10, x, y, 2.5, 1.5, sigma);
    robuste = [rob(1, end); rob(2, end)];
  endif;
  
  figure;
  plot(x,y,'o');xlabel('a'); ylabel('b');
  hold on;
  plot(x,mC(2) + x*mC(1),'-r');
  plot(x, robuste(2) + x*robuste(1), '-g');
  hold off;
  legend({'Points', 'Moindres carrés', 'Robuste'});
  title(strcat('Espace de a et b : [', int2str(space(1)), ', ', int2str(space(2)), '] avec sigma = ', int2str(sigma)));
endfunction
