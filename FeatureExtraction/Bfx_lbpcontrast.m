% Pietikainen, M. et al (2011): Computer Vision Using Local Binary
% Patterns, Springer.


function J = Bfx_lbpcontrast(I,options)

if ~exist('options','var')
    m = 3;
else
    m = options.m;
end
n  = (m-1)/2;

[N,M,P] = size(I);

if P==1 % 2D
    m0 = (m^2+1)/2;
    ii = [1:m0-1 m0+1:m^2];
    J = zeros(N,M);
    for i=n+1:N-n
        for j=n+1:M-n
            s = I(i-n:i+n,j-n:j+n);
            if std2(s)>0
                st = s(ii);
                it = st>=s(m0);
                Jplus  = mean(st(it));
                Jminus = mean(st(not(it)));
                if isnan(Jplus)
                    Jplus = 0;
                    Jminus = 0;
                end
                if isnan(Jminus)
                    Jminus = 0;
                    Jplus = 0;
                end
                J(i,j) = Jplus-Jminus;
            end
        end
    end
else
    m0 = (m^3+1)/2;
    ii = [1:m0-1 m0+1:m^3];
    J = zeros(N,M,P);
    for i=n+1:N-1
        for j=n+1:M-1
            for k=n+1:P-1
                s = I(i-n:i+n,j-n:j+n,k-n:k+n);
                if std2(s)>0
                    st = s(ii);
                    it = st>=s(m0);
                    Jplus  = mean(st(it));
                    Jminus = mean(st(not(it)));
                    if isnan(Jplus)
                        Jplus = 0;
                        Jminus = 0;
                    end
                    if isnan(Jminus)
                        Jminus = 0;
                        Jplus = 0;
                    end
                    J(i,j,k) = Jplus-Jminus;
                end
            end
        end
    end
end
