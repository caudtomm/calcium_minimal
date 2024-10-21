
function [labels, NCopt, NCopt2, param] = adaptAPclust(sw)
% programs for adaptive Affinity Propagation clustering; an improvement
% of Affinity Propagation clusteirng (see Frey & Dueck, Science, Feb. 2007)
% Note: Statistics Toolbox of Matlab needs to be installed
% WANG Kaijun: wangkjun@yahoo.com, July, Sept. 2007.

id = 6;        % selecting a data set, rows - data points, columns - dimensions
algorithm = 1;  % 1 --- adaptive AP, 0 --- original AP
nrun = 50000;   % max iteration times, default 50000
nrun2 = 2000;   % max iteration times for original AP
nconv = 50;     % convergence condition, default 50
pstep = 0.01;   % decreasing step of preferences: pstep*pmedian, default 0.01
lam = 0.5;      % damping factor, default 0.5
cut = 3;        % after clustering, drop an cluster with number of samples < cut
% splot = 'plot'; % observing a clustering process when it is on
splot = 'noplot';

param.algorithm = algorithm;
param.nrun = nrun;
param.nconv = nconv;
param.pstep = pstep;
param.damping = lam;
param.cut = cut;
param.metric = 'correlation';

% initialization
type = 2;       % 1: Euclidean distances % 2: Pearson correlation coefficients
simatrix = 0;   % 0: data as input; 1: similarity matrix as input
data_load      % loading a data file or similarity matrix

NCopt = 0;
NCopt2 = 0;


disp(' '); disp(['==> Clustering is running, please wait ...']);
if algorithm
   tic;
   if simatrix
      [labels,NCs,labelid,iend,Sp,Slam,NCfixs] = adapt_apcluster(M,type,...
        p,pstep,simatrix,'convits',nconv,'maxits',nrun,'dampfact',lam,splot);
   else
      [labels,NCs,labelid,iend,Sp,Slam,NCfixs] = adapt_apcluster(data,type,...
        p,pstep,simatrix,'convits',nconv,'maxits',nrun,'dampfact',lam,splot);
   end
   if simatrix; nrow = size(M,1); end
  [NC,Sil,Silmin] = solution_evaluation(data,M,labels,NCs,...
      NCfixs,simatrix,nrow,type,cut);
  trun = toc;
  fprintf('\n## Running time = %g seconds \n', trun);
  fprintf('## Running iterations = %g \n', iend);
  
  % finding an optimal clustering solution
  solution_findK
  
else
    tic;
    if ~simatrix
       M = simatrix_make(data,type,nrow);
    end
    if ~length(p)
        dn = find(M(:,3)>-realmax);
        p = median(M(dn,3));         % Set preference to similarity median
    end
    [labels,netsim,iend,unconverged] = apcluster(M,p,'convits',...
        nconv,'maxits',nrun2,'dampfact',lam,splot);
    trun = toc;
    fprintf('\n## Running time = %g seconds \n', trun);
    fprintf('## Running iterations = %g \n', iend);
    
    % finding an clustering solution
    solution_findK
end

