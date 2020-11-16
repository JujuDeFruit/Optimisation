function [distance] = GetDistance(approx)
  
  % Récuparation des valeurs de A et de B.
  aApprox = approx(1, :);
  bApprox = approx(2, :);
  
  % Taille maximale pour connaître le nombre d'itération à effectuer.
  N = size(aApprox, 2);
  
  distance = [];
  
  for i=1:N-1
    % La distance correspond à la distance séparant le point actuel du point 
    % suivant. 
    % C'est pour ça que l'itération max est à N - 1.
    distance = [distance, sqrt((aApprox(i+1) - aApprox(i))^2 + (bApprox(i+1) - bApprox(i))^2)];
    % On pousse la nouvelle distance à la fin du tableau.
  endfor
endfunction
