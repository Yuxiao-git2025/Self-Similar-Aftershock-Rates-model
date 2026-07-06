% =========================================================================
% XA.Yu., (2026/07/05)
% SSAR magnitude correlation analysis
clc;
close all;
addpath(pwd,'Data\');

CataFile='ExampleCata.txt';         % The Mc=1.5 and total 376311 events
% CataFile='SyntheticCata.txt';       % The Mc=1.5 and total 27251 events
Cata=SSAR_ReadCatalog(CataFile);      % Time unit is second (not day)
SSAR_PlotMT(Cata(:,1),Cata(:,2),600); % Plot M-T

%% ========================================================================
% >> Analysis settings
rng('shuffle');
McGrids=[1.6 1.8];
m0Grids=-3.0:0.1:3.0;
Nshuffles=300;
ExpdtGrids=[2 2.8];    % i.e., dt<10^Expdt
TimeSigma=1;           % multiples of σ when show results

% >> Conditioning type:
% 'subseq'       : all subsequent pairs
% 'within_dt'    : subsequent pairs with dt < y
% 'TMD'          : subsequent mother-daughter pairs (need more information)
% 'within_dt_MD' : subsequent mother-daughter pairs with dt < y
% conditionType='within_dt';
conditionType='within_dt_MD';

% >> Randomization type:
% 'subcatalog' : randomize within selected pairs, non-trivial correlation
% 'fullcatalog': randomize from full catalog, trivial + non-trivial
randomType='subcatalog';
% randomType='fullcatalog';

% >> Randomize format:
% 'L2' means randomize subsequent/daughter magnitude
% 'L1' means randomize preceding/mother magnitude
shuffleType='L2';


% =========================================================================
% >> Main loop
if length(ExpdtGrids)*length(McGrids)>8
    errStd('# Please adjust the range of parameters');
end
Result=zeros(length(m0Grids),5,length(ExpdtGrids)*length(McGrids));
IDcount=1;

figure;
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for ID1=1:length(ExpdtGrids)
    Expdt=ExpdtGrids(ID1);
    TimeLimit=10^Expdt;
    fprintf('# [%d] dt<10^%.1f =%.2f s\n',ID1,Expdt,TimeLimit);
    for ID2=1:length(McGrids)
        Mc=McGrids(ID2);
        % Filter catalog by magnitude cut-off threshold
        CataNew=Cata(Cata(:,2)>=Mc,:);
        CataNew(isnan(CataNew))=-2;
        fprintf('#     [%d-%d] Mc=%.1f\n',ID1,ID2,Mc);
        fprintf('#           Events ratio=%.2f\n',size(CataNew,1)/size(Cata,1)*100);

        % Build observed magnitude differences
        % pairArr: [m_as(i+1),M(i),dm(=m_as-M),dt(=t_as-t),Ismember]
        [dmag,pairArr,mask]=SSAR_BuildPairs(CataNew,TimeLimit,conditionType);
        fprintf('#           Select ratio=%.2f\n',sum(mask)/(size(CataNew,1)-1)*100);
        fprintf('#           Using=%d\n',sum(mask));

        % Compute dP
        [m0List, deltaP, errStd, probObs, probRandMean]=SSAR_CalDeltaP( ...
            dmag, pairArr, CataNew, m0Grids, Nshuffles, randomType, shuffleType);

        % Save result
        Result(:,:,IDcount)=[m0List(:), deltaP(:), errStd(:), probObs(:), probRandMean(:)];

        % Plot result
        SSAR_PlotResults(m0List, deltaP, errStd, Mc, Expdt, IDcount, TimeSigma);
        IDcount=IDcount+1;
    end
end
fprintf('# (The plotting will show %d*sigma) \n',TimeSigma);

