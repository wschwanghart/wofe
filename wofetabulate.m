function T = wofetabulate(mdl)

% write WofE model to a table
%
% Syntax
%
%     T = wofetabulate(mdl)
%
% Description
%
%
% Input arguments
%
%     mdl     WofE model return by wofe
%
% Output arguments
%
%     T       cell array
%
% See also: wofe, wofepredict
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 07. May, 2015

T{1} = 'Causative factor';
T{2} = 'Unit';
T{3} = 'Binary predictors';
T{4} = 'W+';
T{5} = 's(W+)';
T{6} = 'W-';
T{7} = 's(W-)';
T{8} = 'C';
T{9} = 's(C)';

rowcount = 2;
for r=1:numel(mdl.pred);

    % Name of predictor variable
    name   = mdl.pred(r).name;
    % remove parts in brackets
    name((strfind(name,'(')-1):end) = '';    
    T{rowcount,1} = name; %#ok<*AGROW>
    T{rowcount,2} = mdl.pred(r).unit;
    
    for r2 = 1:mdl.pred(r).nrlevels;
        if mdl.pred(r).nrlevels < numel(mdl.pred(r).levels)
            % Variable was quantized
            T{rowcount,3} = [num2str(mdl.pred(r).levels(r2),3) '~' num2str(mdl.pred(r).levels(r2+1),3)];
        else
            T{rowcount,3} = num2str(mdl.pred(r).levels(r2));
        end
        T{rowcount,4} = mdl.Wpl(r).W(r2);
        T{rowcount,5} = sqrt(mdl.Wpl(r).s2W(r2));
        T{rowcount,6} = mdl.Wmi(r).W(r2);
        T{rowcount,7} = sqrt(mdl.Wmi(r).s2W(r2));
        T{rowcount,8} = mdl.C(r).C(r2);
        T{rowcount,9} = mdl.C(r).sC(r2);
        
        rowcount = rowcount + 1;

    end
end
        