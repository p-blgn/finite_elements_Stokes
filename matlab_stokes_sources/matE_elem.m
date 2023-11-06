function [Eel] = matE_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice elementaire du bloc rectangulaire (p, dv1/dx)
%
% SYNOPSIS [Eel] = matE_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Eel matrice elementaire rectangulaire (matrice 6x3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% calcul de la matrice elementaire du bloc rectangulaire (p, dv1/dx)
% ------------------------------------------------------------------
Eel = zeros(6,3);

% A COMPLETER

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
lambdas = zeros(3,3);
for i=1:3
    lambdas(1,i) = 1-S_hat(i,1)-S_hat(i,2);
    lambdas(2,i) = S_hat(i,1);
    lambdas(3,i) = S_hat(i,2);
end
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
    for J=1:3
        s = 0;
        for q=1:3
            s = s - lambdas(J,q)*dot(inv_BT(1,:),grad_w(:,q,I));
        end
        Eel(I,J) = poids * s * abs(determinant);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
