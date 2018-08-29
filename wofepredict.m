function P = wofepredict(mdl,Cx)

% Prediction
%
% Syntax
%
%     P = wofepredict(mdl,Cx)
%   
% Description
%
%     wofepredict takes the weights of evidence model mdl and a cell array
%     of predictor grids Cx. GRIDobjs in Cx must be in the same order as
%     they have been used to train the model.
%
% Input arguments
%
%     mdl    weights of evidence model
%     Cx     cell array of GRIDobjs. All GRIDobjs must be spatially
%            aligned.
%
% Output arguments
%
%     P      GRIDobj with posterior probabilities
%
%
% See also: wofe, wofetabulate
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 07. May, 2015 

if numel(Cx) ~= numel(mdl.pred)
    error(['number of predictor variables in the model differs from the\newline' ...
           'number of grids supplied as second argument']);
end

logitD  = log(mdl.P_D/(1-mdl.P_D));

P = GRIDobj(Cx{1});

logitD  = repmat(logitD,numel(P.Z),1);

for r=1:numel(mdl.pred);
    [~,IX] = histc(double(Cx{r}.Z(:)),double(mdl.pred(r).levels));
    I = IX~=0;
    Wpl = mdl.Wpl(r).W(IX(I));
    Wpl = Wpl(:);
    
    Wmi = mdl.Wmi(r).W;
    Wmi = repmat(Wmi,mdl.pred(r).nrlevels,1);
    Wmi(eye(mdl.pred(r).nrlevels,'uint8')>0) = 0;
    Wmi = sum(Wmi,2);
    Wmi = Wmi(IX(I));
    
    logitD(I) = logitD(I) + Wpl + Wmi;
    logitD(~I) = nan;
end

P.Z = reshape(exp(logitD)./(1+exp(logitD)),P.size);
    
    
