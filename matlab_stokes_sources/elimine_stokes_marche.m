function [tilde_AA,tilde_LL] = elimine_stokes_marche(AA,LL,Refneu,Coorneu)
N = length(Refneu);
M = length(LL);
tilde_AA = AA;
tilde_LL = LL;
for i=1:N
    if Refneu(i)==1 || Refneu(i)==2 || Refneu(i)==3
        tilde_AA(i,:)=zeros(1,M);
        tilde_AA(N+i,:)=zeros(1,M);
        tilde_AA(i,i) = 1;
        tilde_AA(N+i,N+i) = 1;
        tilde_LL(i) = g1_marche(Coorneu(i,1),Coorneu(i,2),Refneu(i));
        tilde_LL(N+i) = g2_marche(Coorneu(i,1),Coorneu(i,2));
    end
end