close all; 

file = open('data.mat');

x = file.x;
y = file.y_noisy;

%% Estimation au sens des moindres carrés.

dimension = 150;
space = [-25; 40];
a = linspace(space(1), space(2), dimension);
b = a;

%% Question 1:
moindresCarres = ones(dimension);

for i=1:dimension
  for j=1:dimension
    moindresCarres(i,j) = MoindresCarres(a(i), b(j), x, y);
  endfor
endfor

figure(1);
contour(a,b,moindresCarres); 
xlabel('b'); ylabel('a');
title('Moindres carrés');
figure(2);
mesh(a,b,moindresCarres);

%% Question 2:
pointMoindreCarrees = MinMatrix(moindresCarres, space, dimension);
aMin = pointMoindreCarrees(1)
bMin = pointMoindreCarrees(2)

figure;
contour(a, b, moindresCarres, 100); xlabel('b'); ylabel('a');
hold on;
plot(bMin, aMin, '*r');
hold off;
legend({'Moindres carrés', 'Minimum'});
title('Moindres carrés avec en  son centre le minimum calculé en a = 12.267 et b = 1');

figure(3);
plot(x,bMin + x*aMin,'-r');
hold on;
plot(x,y,'o');
hold off;
axis([-1 5 -1 16]);
xlabel('a'); ylabel('b');
title('Moindres carrées pour a = 12.267 et b = 1');
grid('on');
legend({'Moindre carrés', 'Points'});

%% Question 3:

% Création de la nouvelle matrice X nx2.
xExtend = ones(size(x,1), 2);
xExtend(:, 1) = x;

 Calcul du résultat quadratique.
quad = inv(xExtend'*xExtend)*xExtend'*y;
aMinBetterAccuracy = quad(1,1)
bMinBetterAccuracy = quad(2,1)

%% Question 4:
figure(4);
contour(a, b, moindresCarres, 100); xlabel('b'); ylabel('a');
hold on;
plot(bMinBetterAccuracy, aMinBetterAccuracy, '*r');
hold off;
legend({'Moindres carrés', 'Minimum précis'});
title('Moindres carrés avec en  son centre le minimum calculé en a = 12.071 et b = 0.96401');

 Question 5:
figure(5);
plot(x, bMinBetterAccuracy + x * aMinBetterAccuracy, '-r');
hold on;
plot(x, y, 'o');
hold off;
xlabel('a'); ylabel('b');
legend({'Moindres carrés plus précis', 'Points'});
title('Représentation de l''approximation des moindres carrés');

%% Estimation robuste.

dimension = 150;

space = [0; 5];

a = linspace(space(1), space(2), dimension);
b = a;

%% Question 6:
p = @(r) (1/2)*log(1+r.^2);
pSeconde = @(r) (1-(r.^2))/(1+(r.^2)^2);

figure(6);
fplot(p, [-10 10]);
hold on;
fplot(pSeconde, [-10 10]);
hold off;
xlabel('r'); ylabel('p');
title('Représentation de la pénalisation de Cauchy pour sigma = 1');
legend({'Fonction de pénalisation', 'Dérivée seconde'});

%% Question 7:

robuste = zeros(dimension);

for i=1:dimension
  for j=1:dimension
    robuste(i, j) = PenalisationDeCauchy(1, a(i), b(j), x, y);
  endfor
endfor

pointRobuste = MinMatrix(robuste, space, dimension);
aMinRobuste = pointRobuste(1)
bMinRobuste = pointRobuste(2)

figure(7);
contour(a, b, robuste, 40); xlabel('b'); ylabel('a');
legend({'Fonction coût'});
title('Modélisation de l''estimation robuste');

%% Question 8:
aGradient = zeros(dimension);
bGradient = aGradient;

for i=1:dimension
  for j=1:dimension
    aGradient(i, j) = AGradient(1, a(i), b(j), x, y);
    bGradient(i, j) = BGradient(1, a(i), b(j), x, y);
  endfor
endfor

 On crée deux vecteurs qui vont permettrent de servir de repère à la fonction 
 quiver. Ils vont balayer tous les points de l'espace en y référencant le 
 gradient correspondant.
redimensionnerMatrix = zeros(1, dimension^2);
redimensionnerMatrix2 = redimensionnerMatrix;

for i=1:dimension
  redimensionnerMatrix(dimension*(i-1)+1:dimension*i) = linspace(a(i), a(i), dimension);
  redimensionnerMatrix2(dimension*(i-1)+1:dimension*i) = a;
endfor

 On met les matrices gradient à plat sous forme de vecteur.
AGrad = reshape(aGradient, [1, dimension^2]);
BGrad = reshape(bGradient, [1, dimension^2]);

hold on;
quiver(redimensionnerMatrix, redimensionnerMatrix2, AGrad, BGrad);
axis equal;
title('Gradient (a b) de la fonction coût');
hold off;

%% Question 9:
iterations = 2000;
precision = 10^-4;

approxPlusFortePente1 = PlusFortePente(precision, iterations, x, y, 5, 5);
plusFortePente1 = [approxPlusFortePente1(1, end); approxPlusFortePente1(2, end)]
Affichage(x, y, a, b, approxPlusFortePente1, robuste, iterations, 'Plus forte pente 1');

%% Question 10:
approxPlusFortePente2 = PlusFortePente(precision, iterations, x, y, 3, 4);
plusFortePente2 = [approxPlusFortePente2(1, end); approxPlusFortePente2(2, end)]
Affichage(x, y, a, b, approxPlusFortePente2, robuste, iterations, 'Plus forte pente 2');

approxPlusFortePente3 = PlusFortePente(precision, iterations, x, y, 2.5, 1.5);
plusFortePente3 = [approxPlusFortePente3(1, end); approxPlusFortePente3(2, end)]
Affichage(x, y, a, b, approxPlusFortePente3, robuste, iterations, 'Plus forte pente 3');

%% Question 11:

iterations = 10;
approxNewton1 = QuasiNewton(precision, iterations, x, y, 10, 10);
quasiNewton1 = [approxNewton1(1, end); approxNewton1(2, end)]
Affichage(x, y, a, b, approxNewton1, robuste, iterations, 'Quasi-Newton 1');

approxNewton2 = QuasiNewton(precision, iterations, x, y, 2.5, 1.5);
quasiNewton2 = [approxNewton2(1, end); approxNewton2(2, end)];
Affichage(x, y, a, b, approxNewton2, robuste, iterations, 'Quasi-Newton 2');

approxNewton3 = QuasiNewton(precision, iterations, x, y, 5, 5);
quasiNewton3 = [approxNewton3(1, end); approxNewton3(2, end)]
Affichage(x, y, a, b, approxNewton3, robuste, iterations, 'Quasi-Newton 3');

%% Question 12:
PrintPointsMoindresCarresRobuste(quasiNewton2, x, y, space, dimension, 1);
PrintPointsMoindresCarresRobuste(quasiNewton2, x, y, [0, 1], dimension, 1);
PrintPointsMoindresCarresRobuste(quasiNewton2, x, y, [-20, 100], dimension, 1);

PrintPointsMoindresCarresRobuste(quasiNewton1, x, y, space, dimension, 1);

PrintPointsMoindresCarresRobuste(quasiNewton2, x, y, space, dimension, 1);
PrintPointsMoindresCarresRobuste(quasiNewton2, x, y, space, dimension, 100);
PrintPointsMoindresCarresRobuste(quasiNewton2, x, y, space, dimension, 100000000000);