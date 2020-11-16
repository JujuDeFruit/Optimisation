function [alpha] = FletcherLeMarechal(point, d, NMax, grad, x, y, alpha0, alphaL, alphaR, beta1, beta2, lambda)
  
  % Initialisation
  i = 1;
  stop = false;
  alphaArray = [alpha0];
  alpha = alpha0;
  infini = alphaR;
  
  while(i <= NMax && stop == false)
    % Calculons gamma (voir diapo 17 chap. 3).
    gamma = -beta1*grad'*d;
    
    % Cauchy au point initial (voir diapo 17 du chap. 3). f représente la 
    % fonction de Cauchy.
    CauchyX = PenalisationDeCauchy(1, point(1), point(2), x, y);
    
    % Point correspond au point initial auxquel on ajoute alpha * direction
    point = point + alpha * d;
    % Cauchy au point x+ad.
    CauchyXPlusAlphaD = PenalisationDeCauchy(1, point(1), point(2), x, y);
    
    %% Vérification de la condition CW1
    if(CauchyXPlusAlphaD <= CauchyX - alpha * gamma)
      CW1 = true;
    else 
      CW1 = false;
    endif;
    
    % Calculons le gradient de x+ad. pour vérifier la condition CW2 
    % page 18 chap.3.
    gradientOnXPlusAlphaD = [AGradient(1, point(1), point(2), x, y); BGradient(1, point(1), point(2), x, y)];
 
    % Calculons la dégérescence.
    deg = (gradientOnXPlusAlphaD'*d)/(grad'*d);
    
    % Vérification de la condition CW2.
    if(deg <= beta2)
      CW2 = true;
    else
      CW2 = false;
    endif;
    
    % Appliquons les changements si les deux conditions CW1 et CW2 ne sont pas 
    % respectées. On applique les nouvelles valeurs de alpha si une des deux 
    % conditions au moins n'est pas remplie. 
    if (CW1 == false || CW2 == false)
      if (CW1 == false)
        % Pas trop long
        alphaR = alpha;
        alpha = (alphaL + alphaR)/2;
      else 
        % Pas trop court
        alphaL = alpha;
        % Infini est la valeur initiale de AlphaR.
        if (alphaR < infini)
          alpha = (alphaL + alphaR)/2;
        else 
          alpha = alpha * lambda;
        endif
      endif
      % On pousse la nouvelle valeur de alpha dans le tableau des valeurs de 
      % Alpha.
      alphaArray(end + 1) = alpha;
    else
      stop = true;
    endif    
    i++;
  endwhile
  
  % On retourne le min des alpha.
  alpha = min(alphaArray);
endfunction