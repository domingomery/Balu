% Bex_fscombination
%
% Toolbox: Balu
%    Combination of feature selection algorithms.
%
%    This example shows how to combine different feature selection algorthm
%    in orde to obtain the highest performance.
%
%    The evaluation of the performance is using a simple LDA classifier
%
%    Example:
%    load datareal
%    Bex_fscombination
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl


function Bex_fscombination(f,d)

op_sfs.m      = 40;                    % 40 features will be selected
op_sfs.s      = 1;                     % 100% of sample will be used
op_sfs.show   = 1;                     % display results
op_sfs.b.name = 'fisher';              % SFS with Fisher
s             = Bfs_balu(f,d,op_sfs);  % index of selected features
XSFS40        = f(:,s);                % SFS first 40 features
XPCA40        = Bft_pca(f,40);         % first 40 principal components

disp('1. Feature selection with SFS...')
X             = XSFS40(:,1:6);
lda_performance(X,d)                   % function at the end of this code

disp('2. Feature selection with PCA only...')
X             = XPCA40(:,1:6);
lda_performance(X,d)

disp('3. Feature selection with PCA of SFS...')
X             = Bft_pca(XSFS40,6);     % first 6 principal components
lda_performance(X,d)                   % function at the end of this code

disp('4. Feature selection with SFS of PCA and SFS...')
X1            = [XPCA40 XSFS40];       
op_sfs.m      = 6;                     % 6 features will be selected
s             = Bfs_balu(X1,d,op_sfs); % index of selected features
X             = X1(:,s);
lda_performance(X,d)                   % function at the end of this code

disp('5. Feature selection with SFS of PCA of SFS...')
X1            = Bft_pca(XSFS40,20);    % X1 from last step
X2            = [X1 XSFS40];
s             = Bfs_balu(X2,d,op_sfs); % index of selected features
X  = X2(:,s);
lda_performance(X,d)                   % function at the end of this code


disp('6. Feature selection with SFS of PCA and all features...')
% (c) Esteban Cortazar, 2010
X1            = [XPCA40 f];
s             = Bfs_balu(X1,d,op_sfs); % index of selected features
X             = X1(:,s);
lda_performance(X,d)                   % function at the end of this code

disp('7. Feature selection with LSEF of SFS...')
op2.m         = 6;                     % 5 features will be selected
op2.show      = 0;                     % display results
s             = Bfs_lsef(XSFS40,op2);  % index of selected features
X             = X1(:,s);
lda_performance(X,d)                   % function at the end of this code

disp('8. Feature selection with FOSMOD of SFS...')
s             = Bfs_fosmod(XSFS40,op2);% data from last step
X             = X1(:,s);
lda_performance(X,d)                   % function at the end of this code

disp('9. Feature selection with PLSR of SFS...')
X             = Bft_plsr(XSFS40,d,6);
lda_performance(X,d)                   % function at the end of this code

end



function lda_performance(X,d)
op_lda.p      = [];
ds            = Bcl_lda(X,d,X,op_lda);
Bev_performance(d,ds)
enterpause
end


