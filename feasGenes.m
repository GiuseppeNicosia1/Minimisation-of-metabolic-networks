function t = feasGenes(child, essentialCouples, V)
%  ATTENZIONE, MODIFICATA PER NON INCLUDERE LE RIACCENSIONI ***
child = logical(child(1:V));
M = essentialCouples(child,:);
t = false(1,V);
for i = 1:size(M,1)
    t = or(t,M(i,:));
end
temp = 1:V;
t2 = child == 0; % ***
t = temp(~t & t2); % ***
% t = temp(~t);
end
