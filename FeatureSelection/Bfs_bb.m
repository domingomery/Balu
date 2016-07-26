% selec = Bfs_bb(X,d,options)
%
% Toolbox: Balu
%    Feature selection using Branch & Bound for fatures X according to
%    ideal classification d. optins.m features will be selected.
%    options.b.method = 'fisher' uses Fisher objetctive function.
%    options.b.method = 'sp100' uses as criteria Sp @Sn=100%.
%    options.b can be any Balu classifier structure (see example).
%    options.show = 1 display results.
%    selec is the indices of the selected features.
%
% Example 1: Branch and Bound search from using Fisher dicriminant
%            (compare this example with exhaustive serach)
%    load datareal
%    s1 = [279 235 268 230 175 165 207 160 269 157]; %indices using Example
%                                                    %of Bfs_sfs.
%    X1 = f(:,s1);                    % preselected features
%    Xn1 = fn(s1,:);
%    op.m = 3;                        % 3 features will be selected
%    op.show = 1;                     % display results
%    op.b.name = 'fisher';            % SFS with Fisher
%    s = Bfs_bb(X1,d,op);             % index of selected features
%    X2 = X1(:,s);                    % selected features
%    Xn2 = Xn1(s,:)                   % list of feature names
%
% Example 2: Branch and Bound using a KNN classifier
%            (compare this example with exhaustive serach)
%    load datareal
%    s1 = [279 235 268 230 175 165 207 160 269 157]; %indices using Example
%                                                    %of Bfs_sfs.
%    X1 = f(:,s1);                    % preselected features
%    Xn1 = fn(s1,:);
%    op.m = 4;                        % 3 features will be selected
%    op.show = 1;                     % display results
%    op.b.name = 'knn';               % SFS with KNN
%    op.b.options.k = 5;              % 5 neighbors
%    s = Bfs_bb(X1,d,op);             % index of selected features
%    X2 = X1(:,s);                    % selected features
%    Xn2 = Xn1(s,:)                   % list of feature names
%
% (c) Irene Zuccar, PUC-DCC, Jul. 2011
% http://grima.ing.puc.cl


function selec = Bfs_bb(X,d,options)

f       = X;
p       = options.m;
NumFeat = size(f,2);
Bend    = false;

% Initial node in Bnode
Bnode   = nchoosek(1:NumFeat,NumFeat);
% Score for initial node
JBnode  = Bfa_score(f,d,options);

if isnan(JBnode)
    Bnode = [1,Bnode];
else
    Bnode = [JBnode,Bnode];
end
% Add initial node to stack Bopened.
Bopened = [];
Bclosed = [];
Bopened = [Bnode;Bopened];

bestJ   = 0;
selec   = [];
r       = 0;

while (~Bend)
    % Take Bnode from stack and store in Bclosed
    [Bopened,Bclosed,Bnode] = Bfs_bbpop(Bopened,Bclosed);
    JBnode=Bnode(1,1);
    % Prune criteria
    if JBnode > bestJ
        b = size(Bnode,2)-1;
        % If it has not reached the leaves then create the children of Bnode
        if b>p
            % It generates the children that have not already generated,
            % it computes the score and it stores and sorts in Bopened, 
            % it takes the best one.
            h = Bfs_bbchild(Bnode,f,d,Bopened,Bclosed,options);
            Bopened = Bfs_bbpush(Bopened,h);
        end
        % If Bnode is a leave > update prune criteria, and store the
        % selected features
        if b == p
            bestJ = JBnode;
            selec = Bnode;
        end
    end
    % If empty stack then end search in tree
    a = size(Bopened,1);
    if a == 0
        Bend = true;
    end
    % The best result in 10000 iterations
    r=r+1;
    if r==10000
        Bend = true;
        error('Branch and bound for feature selection interrupted: too many iterations.')
    end
end
selec = selec(2:end)';
end



function resp = Bfs_bbis(Bnode,Bopened,Bclosed)
resp = false;
m = size(Bopened,1);
for i = 1:m
    if Bopened(i,2:end) == Bnode(2:end)
        resp = true;
        return;
    end
end
m = size(Bclosed,1);
for i = 1:m
    if Bclosed(i,2:end) == Bnode(2:end)
        resp=true;
        return;
    end
end
end

function resp = Bfs_bbchild(padre,f,d,Bopened,Bclosed,options)

NumFeat = size(Bclosed,2)-1;
nPadre  = size(padre,2)-1;
Bchild    = nchoosek(padre(1,2:nPadre+1),nPadre-1);
Bchild    = [zeros(nPadre,1),Bchild];
Bchild    = padarray(Bchild,[0,NumFeat-nPadre+1],0,'post');
resp=[];
for i=1:nPadre
    if ~Bfs_bbis(Bchild(i,:),Bopened,Bclosed)
        fchild = f(:,Bchild(i,2:nPadre));
        % JBchild = Bfa_jfisher(fchild,d);
        JBchild = Bfa_score(fchild,d,options);
        Bchild(i,1) = JBchild;
        resp = [Bchild(i,:);resp];
    end
end
resp = -sortrows(-resp);
end

function [Bopened,Bclosed,Bnode] = Bfs_bbpop(Bopened,Bclosed)

Bnode = Bopened(1,:);
Bclosed = [Bnode;Bclosed];
aux=Bnode(2:end);
pos = find(aux==0);
q = size(pos,2);
if q~=0
    Bnode(pos+1:end) = [];
end
Bopened(1,:) = [];
end

function Bopened = Bfs_bbpush(Bopened,h)
Bopened = [h;Bopened];
end



