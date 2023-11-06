function [tilde_AA tilde_LL]=elimine(AA, LL, Refneu)
N=size(AA,1);
for i=1:N
    if(Refneu(i)==1) %si le noeud est sur le bord
        AA(i,:)=zeros(1,N); %on annule la ligne
        AA(:,i)=zeros(N,1); %on annule la colonne
        AA(i,i)=1;
        LL(i)=0;
    end
end
tilde_AA=AA;
tilde_LL=LL;
end

