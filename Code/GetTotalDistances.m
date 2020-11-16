function [totalDistances] = GetTotalDistances(distancesMatrix)
  % Distance totale initiale de 0. On crée un tableau qui référence la distance 
  % totale à chaque itération. Ainsi, plus l'indice dans le tableau est élévé, 
  % plus la valeur est élevée.
  totalDistances = [0];
  N = size(distancesMatrix, 2);
  
  for i=1:N
    % On pousse à la fin du tableau la valeur sommée des distances passées.
    totalDistances(i + 1) = totalDistances(end) + distancesMatrix(i);
  endfor
endfunction