% =====================================================
%
% une routine pour la mise en oeuvre des EF Taylor-Hood
% (P^1-P^1) pour l'equation de Stokes
%
% | - nu*\Delta u + \Grad p = f,   dans \Omega
% |                 - div u = 0,   dans \Omega
% |                       u = g,   sur \Gamma_D
% |           nu*du/dn - pn = 0,   sur \Gamma_N
%
% =====================================================

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomMarcheP1.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,~,~,~,~] = lecture_msh(nom_maillage);

nu = 1;      % viscosite

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------

% on ne connait pas encore le nombre de sommets, attention ce n'est pas la moitie du nombre des noeuds !

KK = sparse(Nbpt, Nbpt);    % matrice de rigidite : bloc (u_i, v_i)
EEtmp = sparse(Nbpt, Nbpt);  % bloc (p, v1) carre temporaire, ce bloc sera rectangulaire apres extraction des noeuds correspondant aux sommets
FFtmp = sparse(Nbpt, Nbpt);  % bloc (p, v2) carre temporaire, ce bloc sera rectangulaire apres extraction des noeuds correspondant aux sommets
isavertex = zeros(Nbpt, 1); % tableau booleen "is a vertex" reperant les noeuds qui sont des sommets (et pas des milieux)

% boucle sur les triangles
% ------------------------
for l=1:Nbtri

  % Coordonnees des sommets du triangle
  % A COMPLETER
  S1 = Coorneu(Numtri(l,1),:);
  S2 = Coorneu(Numtri(l,2),:);
  S3 = Coorneu(Numtri(l,3),:);

  % On met 1 dans les entrees de isavertex correspondants aux 3 sommets
  isavertex(Numtri(l,1)) = 1;
  isavertex(Numtri(l,2)) = 1;
  isavertex(Numtri(l,3)) = 1;

  % Calcul des matrices elementaires du triangle l
  Mel = matM_elem(S1, S2, S3);
  Kel = matK_elem(S1, S2, S3);
  Eel = matE_elem_P1(S1, S2, S3); 
  Fel = matF_elem_P1(S1, S2, S3);

  % On fait l'assemblage des blocs (u_i, v_i)
  % A COMPLETER
  for i=1:3
      I = Numtri(l,i);
      for j=1:3
          J = Numtri(l,j);
          KK(I,J) = KK(I,J) + Kel(i,j);
      end
  end

  % On fait l'assemblage des blocs (p, v1) et (p, v2)
  % A COMPLETER
  for i=1:3
      I = Numtri(l,i);
      for j=1:3
          J = Numtri(l,j);
          EEtmp(I,J) = EEtmp(I,J) + Eel(i,j);
          FFtmp(I,J) = FFtmp(I,J) + Fel(i,j);
      end
  end

end % for l

% find(isavertex~=0) donne la liste des numeros des sommets
numsommets = find(isavertex==1);
Ns = length(numsommets);   % nombre de sommets

EE = sparse(Nbpt, Ns);      % bloc rectangulaire (p, v1)
FF = sparse(Nbpt, Ns);      % bloc rectangulaire (p, v2)

% On extrait de EEtmp et FFtmp les colonnes correspondantes a des sommets
% (par ailleurs les autres colonnes sont a priori nulles)
% A COMPLETER

for i=1:Ns
    EE(:,i) = EEtmp(:,numsommets(i));
    FF(:,i) = FFtmp(:,numsommets(i));
end

% Matrice elements finis par blocs
% --------------------------------
AA = sparse(2*Nbpt + Ns, 2*Nbpt + Ns);

% A COMPLETER
AA(1:Nbpt, 1:Nbpt) = nu*KK;
AA((Nbpt+1):(2*Nbpt), (Nbpt+1):(2*Nbpt)) = nu*KK;
AA(1:Nbpt, (2*Nbpt+1):(2*Nbpt + Ns)) = EE;
AA((Nbpt+1):(2*Nbpt), (2*Nbpt+1):(2*Nbpt + Ns)) = FF;
AA((2*Nbpt+1):(2*Nbpt + Ns), 1:Nbpt) = EE';
AA((2*Nbpt+1):(2*Nbpt + Ns), (Nbpt+1):(2*Nbpt)) = FF';

% Pseudo-elimination
% ------------------
% on n'impose Dirichlet que sur les inconnues P1 de la vitesse
LL = zeros(2*Nbpt+Ns,1);  % vecteur second membre
[tilde_AA, tilde_LL] = elimine_stokes_marche(AA,LL,Refneu,Coorneu);

% Resolution du systeme lineaire
% ----------
UU = tilde_AA\tilde_LL;

% visualisation
% -------------
U1 = UU(1:Nbpt);
U2 = UU((Nbpt+1):(2*Nbpt));
P = UU((2*Nbpt+1):(2*Nbpt + Ns));

affiche_ordre1(U1, Numtri, Coorneu,sprintf('U1 - %s', nom_maillage));
affiche_ordre1(U2, Numtri, Coorneu,sprintf('U2 - %s', nom_maillage));

% Dans la i-eme entree de DEnumnoeudAnumsommet on trouve le numero en tant que sommet du i-eme noeud  
% (et on trouve 0 si l'i-eme noeud n'est pas un sommet) 
DEnumnoeudAnumsommet = zeros(Nbpt, 1);
DEnumnoeudAnumsommet(numsommets) = 1:Ns;   
% Pour l'affichage on decale les numeros dans Numtri
NumtriDecale = zeros(Nbtri,3); 
NumtriDecale = DEnumnoeudAnumsommet(Numtri(:,1:3)); 
affiche_ordre1(P, NumtriDecale, Coorneu(numsommets,:),sprintf('pression - %s', nom_maillage));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022