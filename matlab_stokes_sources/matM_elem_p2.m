function [Mel] = matM_elem_p2(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice de masse elementaire en P2 Lagrange
%
% SYNOPSIS [Mel] = matM_elem_p2(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Mel matrice de masse elementaire (matrice 6x6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% calcul de la matrice de masse
% -------------------------------
% Initialisation
Mel = zeros(6); % A COMPLETER

% Points et poids de quadrature
S_hat = [0.0915762135098, 0.0915762135098;
         0.8168475729805, 0.0915762135098;
         0.0915762135098, 0.8168475729805;
         0.1081030181681, 0.4459484909160;
         0.4459484909160, 0.1081030181681;
         0.4459484909160, 0.4459484909160];
poids = [0.05497587183, 0.05497587183, 0.05497587183, 0.1116907948, 0.1116907948, 0.1116907948];

% A COMPLETER
a = x2 - x1;
b = x3 - x1;
c = y2 - y1;
d = y3 - y1;
lambdas = zeros(3,6);
for i=1:6
    lambdas(1,i) = 1-S_hat(i,1)-S_hat(i,2);
    lambdas(2,i) = S_hat(i,1);
    lambdas(3,i) = S_hat(i,2);
end
w = zeros(6);
w(1,:) = lambdas(1,:).*(2*lambdas(1,:)-1);
w(2,:) = lambdas(2,:).*(2*lambdas(2,:)-1);
w(3,:) = lambdas(3,:).*(2*lambdas(3,:)-1);
w(4,:) = lambdas(1,:).*lambdas(2,:)*4;
w(5,:) = lambdas(2,:).*lambdas(3,:)*4;
w(6,:) = lambdas(3,:).*lambdas(1,:)*4;
determinant = a*d - b*c;
for I=1:6
    for J=1:6
        s = 0;
        for q=1:6
            s = s + poids(q)*w(I,q)*w(J,q);
        end
        Mel(I,J) = abs(determinant) * s;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
