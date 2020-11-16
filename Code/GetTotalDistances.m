function [totalDistances] = GetTotalDistances(distancesMatrix)
  % Distance totale initiale de 0. On cr�e un tableau qui r�f�rence la distance 
  % totale � chaque it�ration. Ainsi, plus l'indice dans le tableau est �l�v�, 
  % plus la valeur est �lev�e.
  totalDistances = [0];
  N = size(distancesMatrix, 2);
  
  for i=1:N
    % On pousse � la fin du tableau la valeur somm�e des distances pass�es.
    totalDistances(i + 1) = totalDistances(end) + distancesMatrix(i);
  endfor
endfunction