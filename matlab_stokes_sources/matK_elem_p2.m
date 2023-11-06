function [Kel] = matK_elem_p2(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice de rigidité elementaire en P2 Lagrange
%
% SYNOPSIS [Kel] = matK_elem_p2(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de rigidité elementaire (matrice 6x6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% calcul de la matrice de rigidité
% -------------------------------
% Initialisation
Kel = zeros(6); % A COMPLETER

% Points et poids de quadrature
S_hat = [1/6, 1/6;
         2/3, 1/6;
         1/6, 2/3];
poids = 1/6;

% A COMPLETER
a = x2 - x1;
b = x3 - x1;
c = y2 - y1;
d = y3 - y1;
determinant = a*d - b*c;
grad_w = zeros(2,3,6);
grad_w(1,:,1) = -3+4*S_hat(:,1)+4*S_hat(:,2);
grad_w(2,:,1) = -3+4*S_hat(:,1)+4*S_hat(:,2);
grad_w(1,:,2) = 4*S_hat(:,1)-1;
grad_w(2,:,2) = 0;
grad_w(1,:,3) = 0;
grad_w(2,:,3) = 4*S_hat(:,2)-1;
grad_w(1,:,4) = 4-4*S_hat(:,2)-8*S_hat(:,1);
grad_w(2,:,4) = -4*S_hat(:,1);
grad_w(1,:,5) = 4*S_hat(:,2);
grad_w(2,:,5) = 4*S_hat(:,1);
grad_w(1,:,6) = -4*S_hat(:,2);
grad_w(2,:,6) = 4-4*S_hat(:,1)-8*S_hat(:,2);
inv_BT = 1/determinant*[d, -c;-b, a];
for I=1:6
    for J=1:6
        s = 0;
        for q=1:3
            s = s + dot(inv_BT*grad_w(:,q,I),inv_BT*grad_w(:,q,J));
        end
        Kel(I,J) = poids * s * abs(determinant);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
