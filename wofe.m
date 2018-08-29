function mdl = wofe(Cx,Y,varargin)

% Weights of Evidence
%
% Syntax
%
%     mdl = wofe(Cx,Y)
%     mdl = wofe(Cx,Y,pn,pv,...)
%
% Description
%
%     Weights of evidence is a method for combining evidence in support of
%     a hypothesis. 
%
% Input
%
%     Cx      npred*2 cell vector with npred variables (GRIDobj) and 
%             vectors defining the edges for quantizing continuous metric 
%             variables, such that Cx = {GRIDobj1 edges1 GRIDobj2 edges2
%             ...}. If the GRIDobj contains ordinal or nominal data, supply
%             the corresponding edges vector as empty vector or string.
%     Y       logical GRIDobj containing the response variable.
%
% Parameter name/value pairs
%
%     mask    logical grid (GRIDobj) that is used as mask
%
% Output
%
%     mdl      structure array with output
%             .P_D    : prior probability. This is an empirical prior and
%                       is estimated from the group relative frequencies. 
%             .P_BD   : conditional probability for each variable. For each
%                       variable, there is a 4xn corresponding to the n 
%                       classes and the probabilities p(B|D), P(B|~D), 
%                       P(~B|D), and P(~B|~D).
%             .Wp     : positive weights
%             .Wm     : negative weights
%             .s2Wp   : variance of Wp
%             .s2Wm   : variance of Wm
%             .C      : contrast (Wp-Wm)
%             .P_DB   : posterior probability (D given B)
%             .P_DnB  : posterior probability (D given ~B)
%
% See also: wofepredict, wofetabulate
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 07. May, 2015


% Cy is a binary variable as a GRIDobj
% Cx is a cell array that contains GRIDobjs and additional information
% in the order
%     {GRIDobj1      categorical levels ... 
%      GRIDobj2      categorical levels ...
%      GRIDobj3      []}

% In addition following parameter value pairs
% 'mask' - GRIDobj
% 


%% Load predictor variables

% Input parsing
p = inputParser;
p.FunctionName = 'wofe';
addRequired(p,'Cx',@(x) iscell(x) && mod(numel(x),2)==0)
addRequired(p,'Y',@(x) isa(x,'GRIDobj') && islogical(x.Z));
addParameter(p,'mask',[]);
parse(p,Cx,Y,varargin{:});

% input argument
y = Y.Z(:);
npred = numel(Cx)/2;

% check if input arguments align spatially with each other.
for r = 1:npred
    tf = validatealignment(Y,Cx{r*2-1});
    if ~tf
        error('TopoToolbox:wofe',...
             ['Predictor variable #%d fails to align\newline'...
              'with the response variable.'],r);
    end
end

% preallocate array that holds the predictor variables 
X = zeros(prod(Cx{1}.size),npred);

% quantize grids if required and write information to mdl
for r=1:npred
    
    if ~isempty(Cx{r*2})
        % quantize
        levels = Cx{r*2};
        levels = levels(:);
        minx = min(Cx{r*2-1});
        maxx = max(Cx{r*2-1});
        
        % adjust levels to cover the full range of the data
        if levels(1)>minx
            levels(1) = minx;
        end
        
        if levels(end) < maxx
            levels(end+1) = maxx+1;
        elseif levels(end) == maxx
            levels(end) = maxx+1;
        end          
        
        X(:,r) = Cx{r*2-1}.Z(:);
        mdl.pred(r).levels = levels;
        levels(end) = maxx;
        mdl.pred(r).levelcenters = levels(1:end-1)+diff(levels)/2;
        mdl.pred(r).nrlevels = numel(levels)-1;
        
    else
        % no need to quantize since data was supplied as binary or
        % multiclass maps
        levels = unique(Cx{r*2-1}.Z(~isnan(Cx{r*2-1}.Z(:))));
        levels(isnan(levels)) = [];
        mdl.pred(r).levels = levels;
        mdl.pred(r).levelcenters = levels;
        mdl.pred(r).nrlevels = numel(levels);
        X(:,r) = Cx{r*2-1}.Z(:);
    end
    
    % add name and units
    mdl.pred(r).name = Cx{r*2-1}.name;
    mdl.pred(r).unit = Cx{r*2-1}.zunit;
    
end

% Find nans, infs and restrict data to mask if supplied
I = ~any(isnan(X) | isinf(X),2) | isnan(y) | isinf(y);
if ~isempty(p.Results.mask)
    I = I & p.Results.mask.Z(:);
end

X = X(I,:);
y = y(I);
y = y>0;
nrc = numel(y);

%% Calculate probabilities, weights and contrasts

% Prior probability
N_D     = nnz(y);
mdl.P_D = N_D./nrc;

logitD  = log(mdl.P_D/(1-mdl.P_D));

% Conditional probability
for r = 1:npred
    [~,IX] = histc(X(:,r),double(mdl.pred(r).levels));
    % create dummy binary variables
    x = sparse((1:nrc)',IX,1,nrc,mdl.pred(r).nrlevels);
    
    for r2 = 1:mdl.pred(r).nrlevels
        % Conditional probabilities
        mdl.P_BD(r).P(1,r2) = nnz(x(:,r2) & y)/N_D;    % P(B|D)
        mdl.P_BD(r).P(2,r2) = nnz(x(:,r2) & ~y)/(nrc-N_D);   % P(B|~D)
        mdl.P_BD(r).P(3,r2) = nnz(~x(:,r2) & y)/N_D;   % P(~B|D)
        mdl.P_BD(r).P(4,r2) = nnz(~x(:,r2) & ~y)/(nrc-N_D);  % P(~B|~D)
        
        % Positive weights (if no data in binary pattern, set weights to
        % zero)
        mdl.Wpl(r).W(1,r2)   = log(mdl.P_BD(r).P(1,r2)/mdl.P_BD(r).P(2,r2));
        if isinf(mdl.Wpl(r).W(1,r2))
            mdl.Wpl(r).W(1,r2) = 0;
        end
        mdl.Wpl(r).s2W(1,r2) = 1/(mdl.P_BD(r).P(1,r2)*N_D) + 1/(mdl.P_BD(r).P(2,r2)*(nrc-N_D));
        if isinf(mdl.Wpl(r).s2W(1,r2))
            mdl.Wpl(r).s2W(1,r2) = nan;
        end
        
        % Negative weights
        mdl.Wmi(r).W(1,r2)   = log(mdl.P_BD(r).P(3,r2)/mdl.P_BD(r).P(4,r2));
        mdl.Wmi(r).s2W(1,r2) = 1/(mdl.P_BD(r).P(3,r2)*N_D) + 1/(mdl.P_BD(r).P(4,r2)*(nrc-N_D));
        
        % Contrasts
        mdl.C(r).C(1,r2)   = mdl.Wpl(r).W(1,r2) - mdl.Wmi(r).W(1,r2);
        mdl.C(r).sC(1,r2)  = sqrt(mdl.Wpl(r).s2W(1,r2) + mdl.Wmi(r).s2W(1,r2));
        
        % Logit D|B
        LogitDBpl = logitD + mdl.Wpl(r).W(1,r2);
        LogitDBmi = logitD + mdl.Wmi(r).W(1,r2);
        
        % Conditional probabilities for each binary pattern
        mdl.P_DBpl(r).P(1,r2) = exp(LogitDBpl)/(1+exp(LogitDBpl));
        mdl.P_DBmi(r).P(1,r2) = exp(LogitDBmi)/(1+exp(LogitDBmi));
  
    end
end
    
%% and this is it











    



