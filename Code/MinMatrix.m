function [point] = MinMatrix(matrix, space, dimension)
  [M,I] = min(matrix);
  [mini, column] = min(M);
  row = I(column);
  echantillonA = (space(2) - space(1))/dimension;
  echantillonB = (space(2) - space(1))/dimension;
  aMin = space(1) + row * echantillonA;
  bMin = space(1) + column * echantillonB;
  point = [aMin; bMin];
endfunction
