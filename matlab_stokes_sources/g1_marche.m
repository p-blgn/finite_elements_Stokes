function val = g1_marche(x,y,ref)
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
if ref==2
    val = -16*(y-1)*(y-0.5);
else
    val = 0;
end
end