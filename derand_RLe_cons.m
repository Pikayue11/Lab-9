%  de_rand binomial
%   RL random local, r1 by tournament, Kaelo, Ali 2007
function y=derand_RLe_cons(P,hodf,hodviol,F,CR,expt)
N=length(P(:,1));
d=length(P(1,:));
% prd1=size(P)
% prd=expt(1)
y=P(expt(1),1:d);
vyb=nahvyb_expt(N,3,expt);	% three random points without expt

r123=P(vyb,:);
hodf123=hodf(vyb);
hodviol123=hodviol(vyb);

trivybrane=[r123 hodf123 hodviol123];
trivybrane=sortrows(trivybrane,d+1);
trivybrane=sortrows(trivybrane,d+2);

r1=trivybrane(1,1:d);
if rand  < 0.5
    r2=trivybrane(2,1:d);
    r3=trivybrane(3,1:d);
else
    r2=trivybrane(2,1:d);
    r3=trivybrane(3,1:d);
end

v=r1+F*(r2-r3);
change=find(rand(1,d)<CR);
if isempty(change) % at least one element is changed
    change=1+fix(d*rand(1));
end
y(change)=v(change);
%