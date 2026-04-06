function [S, H] = bss_pca(X);

% Separation aveugle de sources.
%
% Analyse en composantes principales (principal component analysis, PCA).
%
% SYNTAX: [S, H] = bss_pca(X);
%
%       S  : sources estimees (normalisees a une puissance unite)
%       H  : matrice de melange estimee
%       X  : signaux observes.

% Nombre de capteurs et d'echantillons
[p, T] = size(X);

% Decomposition en valeurs singulieres de la matrice de donnees
[V, S, U] = svd(X', 0);

% Matrice de melange estimee
H = U*S/sqrt(T);

% Sources estimees (normalisees a puissance unite)
S = sqrt(T)*V';
