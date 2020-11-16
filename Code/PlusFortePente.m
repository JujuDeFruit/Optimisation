function [listPoints] = PlusFortePente(epsilon, NMax, x, y, InitialA, InitialB)
  % Initialisation
  % - Gradient initial (le gradient est stocké sous forme de vecteur).
  aGradientInitial = AGradient(1, InitialA, InitialB, x, y);
  bGradientInitial = BGradient(1, InitialA, InitialB, x, y);
  grad = [aGradientInitial; bGradientInitial];
  
  % - Norme initiale.
  norme = sqrt(aGradientInitial^2 + bGradientInitial^2);
  
  % - Point initial sous forme de vecteur.
  point = [InitialA; InitialB];
  
  % Liste des points à retourner.
  listPoints = [point];
  
  k=1;
  
  % Epsilon correspond à la précision demandée.
  while(k < NMax && norme > epsilon)
    % Calcul de la direction d. 
    d = -grad;
    
    % On appelle l'algorithme de Fletcher-Lemarechal pour trouver la valeur de 
    % Alpha.
    alpha = FletcherLeMarechal(point, d, 10, grad, x, y, 10^-3, 0, 10^9, 10^-3, 0.99, 20);
    
    % On passe au point suivant.
    point = point + alpha * d;
    
    % On pousse le nouveau point à la fin de la liste des points à retourner.
    listPoints(:, end + 1) = point;
    
    % Calcul du gradient du nouveau point.
    aGradient = AGradient(1, point(1), point(2), x, y);
    bGradient = BGradient(1, point(1), point(2), x, y);
    grad = [aGradient; bGradient];
    
    % Calcul de la nouvelle norme.
    norme = sqrt(grad(1)^2 + grad(2)^2);
    
    k++;
  endwhile
endfunction
