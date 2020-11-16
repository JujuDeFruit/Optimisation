function [approx] = QuasiNewton(epsilon, NMax, x, y, InitialA, InitialB, sigma = 1)
  % Initialisation
  % - Matrice identité 2x2.
  I = [1 0; 0 1];
  
  % - H initiale est la matrice identité de dimension 2, à l'initialisation.
  H = I;
  
  % - Point initial sous forme de vecteur. 
  point = [InitialA; InitialB];
  
  % - Calcul du gradient initial pour trouver la direction d initiale. 
  aGrad = AGradient(sigma, InitialA, InitialB, x, y);
  bGrad = BGradient(sigma, InitialA, InitialB, x, y);
  
  % - Initialisation du tableau retourné qui contient l'ensemble des points 
  % (a; b) sous forme de vecteur 2x1.
  approx = [point];
  
  % - Calcul de la norme initiale
  norm = sqrt(aGrad^2 + bGrad^2);
  
  k = 0;
  
  while (k < NMax && norm > epsilon)

    % Nouveau gradient (matrice 2x1).
    grad = [aGrad; bGrad];
    
    % Calcul de la nouvelle direction (matrice 2x1).
    d = -H*grad;
  
    % On appelle l'alogrithme de Fletcher-Lemarechal pour trouver la valeur de 
    % Alpha
    alpha = FletcherLeMarechal(point, d, 10, grad, x, y, 1, 0, 10^6, 10^-3, 0.99, 20);
    
    % Calcul du nouveau point.
    point = point + alpha * d;
    
    k++;
    
    % Calcul du nouveau gradient.
    aGrad = AGradient(sigma, point(1), point(2), x, y);
    bGrad = BGradient(sigma, point(1), point(2), x, y);
    
    % Calcul du nouveau H.
    % Y_k-1.
    yk = [aGrad; bGrad] - grad;
    % d_k-1 barre.
    d_ = alpha*d;
    
    % On découpe le calcul de H en H1, H2 et H3 pour que ce soit + lisible.
    H1 = I-((d_*yk')/(d_'*yk));
    H2 = I-((yk*d_')/(d_'*yk));
    H3 = (d_*d_')/(d_'*yk);
    
    H = (H1*H*H2)+H3;
    
    % Calcul de la nouvelle norme.
    norm = sqrt(aGrad^2 + bGrad^2);
    
    % On pousse le nouveau point dans la liste des points à retourner.
    approx = [approx, point];
      
  endwhile
endfunction