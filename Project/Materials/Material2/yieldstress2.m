function sy = yieldstressm2(peteff, dmat)
% m2yieldstress computes the yield stress for the given state
%
% INPUT
%       peteff:     the plastic effective strain
%       dmat:       matieral parameters,
%                       - alpha
%                       - n
%                       - sy0:  Initial yield stress 
%
% OUTPUT
%       da:         the yield stress
alpha = dmat(1);
n = dmat(2);
sy0 = dmat(3);

sy = sy0 + sy0*alpha*(peteff)^(n);
end