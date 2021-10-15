function da = hardrate1(peteff, dmat)
% m1hardrate computes the hardening rate for the given state
%
% INPUT
%       peteff:     the plastic effective strain
%       dmat:       matieral parameters,
%                       - kinf
%                       - h
%                       - sy0:  Initial yield stress (Not used, but for
%                       consistancy it can be good to specify)
%
% OUTPUT
%       da:         the hardening rate
kinf = dmat(1);
h = dmat(2);

da = h*exp(-h/kinf*peteff);
end