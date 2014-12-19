function  [Opts]=OptsSet(nSig,x_org)

Opts.org       =   x_org;
Opts.nSig      =   nSig;                                 % Variance of the noise image
Opts.SearchWin =   30;                                   % Non-local patch searching window
Opts.delta     =   0.2;                                  % Optsameter between each iter
Opts.c         =   2*sqrt(2);                            % Constant num for the weight vector
Opts.Innerloop =   2;                                    % InnerLoop Num of between re-blockmatching
Opts.ReWeiIter =   3;
if nSig<=20
    Opts.PatchSize    =   6;                            % Patch size
    Opts.ArrayNo        =   70;                           % Initial Non-local Patch number
    Opts.Iter          =   8;                            % total iter numbers
    Opts.lamada        =   0.54;                         % Noise estimete Optsameter
elseif nSig <= 40
    Opts.PatchSize       =   6;
    Opts.ArrayNo        =   90;
    Opts.Iter          =   5;
    Opts.lamada        =   0.56; 
elseif nSig<=60
    Opts.PatchSize       =   8;
    Opts.ArrayNo        =   120;
    Opts.Iter          =   14;
    Opts.lamada        =   0.58; 
else
    Opts.PatchSize       =   9;
    Opts.ArrayNo        =   140;
    Opts.Iter          =   14;
    Opts.lamada        =   0.58; 
end

Opts.step      =   floor((Opts.PatchSize)/2-1);                   
