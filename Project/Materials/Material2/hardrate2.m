function da = hardrate2(peteff, dmat)
% m2hardrate computes the hardening rate for the given state
%
% INPUT
%       peteff:     the plastic effective strain
%       dmat:       matieral parameters,
%                       - alpha
%                       - n
%                       - sy0:  Initial yield stress
%
% OUTPUT
%       da:         the hardening rate
alpha = dmat(1);
n = dmat(2);
sy0 = dmat(3);

da = n*alpha*sy0*(peteff)^(n-1);
end