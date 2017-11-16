%function [errL2, errH1] = principal_periodique(eta)
eta = 0.1;
% =====================================================
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions periodiques
% sur le maillage nom_maillage.msh
%
% | -div(A grad u ) u + u= f,   dans \Omega
% |         u periodique,   sur le bord
% =====================================================

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre_per.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  Nl = [Numtri(l,1), Numtri(l,2), Numtri(l,3)];
  S1 = Coorneu(Nl(1),:);
  S2 = Coorneu(Nl(2),:);
  S3 = Coorneu(Nl(3),:);
  % calcul des matrices elementaires du triangle l 
  
   Kel=matKvar_elem(S1, S2, S3);
           
   Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemblage de la matrice globale et du second membre
  for i = 1:3
      for j = 1:3
          KK(Nl(i),Nl(j)) = KK(Nl(i),Nl(j)) + Kel(i,j);
          MM(Nl(i),Nl(j)) = MM(Nl(i),Nl(j)) + Mel(i,j);
      end
  end
end % for l
MM = eta * MM;

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = f(Coorneu(:,1),Coorneu(:,2));
LL = MM * FF;

% Projection sur l espace V_p
% matrice de projection 
PP = eye(Nbpt);
Nbords = sum(Refneu(:,:)==2);
PP(1:4,1:4) = ones(4)/4;
for k = 1:Nbords
    SO = 4              +  k;
    ES = 4 +     Nbords +  k;
    NO = 4 + 3 * Nbords - (k-1);
    WQ = 4 + 4 * Nbords - (k-1);
    [PP(SO,SO), PP(SO,NO), PP(NO,NO), PP(NO,SO)] = deal(1/2);
    [PP(ES,ES), PP(ES,WS), PP(WS,WS), PP(WS,ES)] = deal(1/2);
end
PPt = PP';
%
AA = MM+KK;
AAp = PP * AA * PPt;
LLP = PP * LL;

% inversion
% ----------
UUp = AAp\LLp;

% Expression de la solution dans toute la base
% ———————
UU = PPt * UUp;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Periodique - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
EE = UU_exact - UU;
% Calcul de l erreur L2
errL2 = EE.' * MM * EE;
UL2 = UU_exact.' * MM * UU_exact;
errL2 = errL2 / UL2;
% Calcul de l erreur H1
errH1 = EE.' * KK * EE;
UH1 = UU_exact.' * KK * UU_exact;
errH1 = errH1 / UH1;
% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end
