function [cost, norm] = GetCostAndNormEvolution(a, b, x, y)
  N = size(a, 2);
  cost = [];
  norm = [];
  for i=1:N
    % On récupère toutes les pénalisation de cauchy pour les différents points 
    % (a b)' pour voir l'évolution.
    cost(end + 1) = PenalisationDeCauchy(1, a(i), b(i), x, y);
    
    % De même, on récupère la norme de chaque gradient en tout point (a b)'.
    aGradient = AGradient(1, a(i), b(i), x, y);
    bGradient = BGradient(1, a(i), b(i), x, y);
    norm(end + 1) = sqrt(aGradient^2 + bGradient^2);
  endfor
endfunction