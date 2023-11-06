function val = g1(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la contrainte g1 en (x,y)
%
% SYNOPSIS val = g1(x,y)
%
% INPUT * x,y : les 2 coordonnees du point
%
% OUTPUT - val : valeur de g1 au point considéré
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val = 4*(1-y)*y;
end