function [] = Affichage(x, y, a, b, approx, robuste, iterations, nom)
  aApprox = approx(1, :);
  bApprox = approx(2, :);

  figure;
  contour(a, b, robuste, 20); title('Représentation 2D de la fonction robuste'); xlabel('b'); ylabel('a');
  hold on;
  plot(bApprox, aApprox, '-r');
  hold off;
  legend({'Pente', nom});

  figure;
  subplot(2,2,1);
  distance = GetDistance(approx);
  plot(distance); 
  xlabel('Nombre d''iterations'); ylabel('Distance'); 
  title('Distance en fonction du nombre d''itérations');
  axis([0 iterations]);

  subplot(2,2,2);
  [cost, norm] = GetCostAndNormEvolution(aApprox, bApprox, x, y);
  plot(cost);
  hold on;
  plot(norm);
  hold off;
  xlabel('Nombre d''iterations');
  axis([0 iterations]);
  title('Évolution de la fonction coût et de la norme du gradient');
  legend({'Fonction coût', 'Norme du gradient'});

  subplot(2,2,3);
  totalDistances = GetTotalDistances(distance);
  plot(totalDistances); 
  xlabel('Nombre d''iterations'); ylabel('Distance totale'); 
  title('Distance totale en fonction du nombre d''itérations');
  axis([0 iterations]);
endfunction
