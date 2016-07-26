function bcl = Bcl_build(varargin)

v = varargin;
n = nargin;
if compare(class(v),'cell')==0
    v = v{1};
    n = length(v);
end



for k=1:n
    s = lower(char(v(k)));
    ns = length(s);
    ok = 0;
    if compare(s(1:3),'knn')==0
        bcl(k).name = 'knn';
        bcl(k).options.k = str2num(s(4:end));
        ok = 1;
    end
    
    if compare(s(1:3),'svm')==0
        bcl(k).name = 'svmplus';
        bcl(k).options.kernel = str2num(s(4:end));
        ok = 1;
    end 
    
    if ns>5
        if compare(s(1:5),'nnglm')==0
            bcl(k).name = 'nnglm';
            bcl(k).options.method = str2num(s(6:end));
            bcl(k).options.iter   = 12;
            
            ok = 1;
        end
    end
    
    if not(ok)
        switch lower(s)
            case 'lda'
                bcl(k).name = 'lda';
                bcl(k).options.p = [];     
                
            case 'qda'
                bcl(k).name = 'qda';
                bcl(k).options.p = [];     
            case 'maha'
                bcl(k).name = 'maha';
                bcl(k).options = [];       
                
            case 'dmin'
                bcl(k).name = 'dmin';
                bcl(k).options = [];       
            case 'pnn'
                bcl(k).name = 'pnn';
                bcl(k).options.p = [];     
                
            case 'adaboost'
                bcl(k).name = 'adaboost';
                bcl(k).options.iter = 10;  
                
            case 'boosting'
                bcl(k).name = 'boosting';
                bcl(k).options.s = 0.3;    
            otherwise
                error('Bcl_build does not recognize %s as a classifier method.',s)
        end
    end
end