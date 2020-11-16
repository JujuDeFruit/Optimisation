function [distance] = GetDistance(approx)
  
  % R�cuparation des valeurs de A et de B.
  aApprox = approx(1, :);
  bApprox = approx(2, :);
  
  % Taille maximale pour conna�tre le nombre d'it�ration � effectuer.
  N = size(aApprox, 2);
  
  distance = [];
  
  for i=1:N-1
    % La distance correspond � la distance s�parant le point actuel du point 
    % suivant. 
    % C'est pour �a que l'it�ration max est � N - 1.
    distance = [distance, sqrt((aApprox(i+1) - aApprox(i))^2 + (bApprox(i+1) - bApprox(i))^2)];
    % On pousse la nouvelle distance � la fin du tableau.
  endfor
endfunction
