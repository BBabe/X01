function Mel = matM_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM_elem :
% calcul la matrices de masse elementaire en P1 lagrange
%
% SYNOPSIS Mel = matM_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Mel matrice de masse elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

x12 = x1-x2; x23 = x2-x3; x31 = x3-x1;
y12 = y1-y2; y23 = y2-y3; y31 = y3-y1;

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

% calcul de la matrice de masse
% -----------------------------
Mel = zeros(3,3);
for i=1:3
	for j=1:3
		% A COMPLETER
		Mel(i,j) = lambda(i,(x1+x2)/2,(y1+y2)/2)...
            * lambda(j,(x1+x2)/2,(y1+y2)/2);
		Mel(i,j) = Mel(i,j) + lambda(i,(x1+x3)/2,(y1+y3)/2)...
            * lambda(j,(x1+x3)/2,(y1+y3)/2);
		Mel(i,j) = Mel(i,j) + lambda(i,(x2+x3)/2,(y2+y3)/2)...
            * lambda(j,(x3+x2)/2,(y3+y2)/2);
	end; % j
end; % i
Mel = Mel / (6*abs(D));

function z = lambda(i,x,y)
    switch i
        case 1
            z = y23*(x-x3) - x23*(y-y3);
        case 2
            z = y31*(x-x1) - x31*(y-y1);
        case 3
            z = y12*(x-x2) - x12*(y-y2);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
