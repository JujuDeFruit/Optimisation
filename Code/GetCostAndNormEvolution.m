function [cost, norm] = GetCostAndNormEvolution(a, b, x, y)
  N = size(a, 2);
  cost = [];
  norm = [];
  for i=1:N
    % On r�cup�re toutes les p�nalisation de cauchy pour les diff�rents points 
    % (a b)' pour voir l'�volution.
    cost(end + 1) = PenalisationDeCauchy(1, a(i), b(i), x, y);
    
    % De m�me, on r�cup�re la norme de chaque gradient en tout point (a b)'.
    aGradient = AGradient(1, a(i), b(i), x, y);
    bGradient = BGradient(1, a(i), b(i), x, y);
    norm(end + 1) = sqrt(aGradient^2 + bGradient^2);
  endfor
endfunction