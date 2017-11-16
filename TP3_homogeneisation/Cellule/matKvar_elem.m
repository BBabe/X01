function  [Kel] = matKvar_elem(S1,S2,S3)

M = [1/3 1/3 ; 1/5 1/5 ; 1/5 3/5 ;3/5 1/5];
w = [-9/32 25/96 25/96 25/96];

syms x y
fref = x*S2.' + y*S3.' + (1-x-y)*S1.';
jac = jacobian(fref,[x y]);
jacdet = abs(det(jac));
jac = jac.';
jac = inv(jac);

frefA = subs(fref,{x,y},{M(:,1).',M(:,2).'});
jacdetA = subs(jacdet,{x,y},{M(:,1).',M(:,2).'});
jacA = subs(jac,{x,y},{M(:,1).',M(:,2).'});

gradw = [-1 1 0;-1 0 1];
Kel = zeros(3,3);
for i=1:3
    for j=1:3
        for p = 1:4
            Kel(i,j) = Kel(i,j) + w(p)...
                *(transpose(mat_A(frefA(:,p))*jacA(:,[p 4+p])*gradw(:,i))...
                  * (jacA(:,[p 4+p])*gradw(:,j)))*jacdetA(p);
        end
    end; % j
end; % i

end