%function Bio_latextable(row_names,col_names,fmt,T)
%
% Toolbox: Balu
%    Code for a latex table.
%    
%    row_names is a cell with the names of the rows.
%    col_names is a cell with the names of the columns.
%    fmt is a cell with the format of each column.
%    T is the table.
%
% Example:
%
% col_names = {'cols','col 1','col 2','col 3','col4'};
% row_names = {'row 1','row 2','row 3'};
% fmt = {'%5.2f','%6.4f','%3.1f','%7.4f'};
% T = rand(3,4);
% Bio_latextable(row_names,col_names,fmt,T)
%
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%
function Bio_latextable(row_names,col_names,fmt,T)

disp('%')
disp('% Latex code starts here:')
disp('%')

[N,M] = size(T);

s = char(ones(M+1,1)*' c ')';
fprintf('\\begin{table}\n')
fprintf('\\caption{Please write caption here.}\n')
fprintf('\\begin{center}\n');
fprintf('\\begin{tabular}{ %s }\n',s(:));
fprintf('\\hline\n');
for j=1:M
    fprintf('  %s &',col_names{j});
end
fprintf('  %s \\\\\n',col_names{M+1});
fprintf('\\hline\n');
for i=1:N
    fprintf('  %s ',row_names{i});
    for j=1:M
        s = ['fprintf(' char(39) ' & ' fmt{j} char(39) ',T(i,j));'];
        eval(s);
    end
    fprintf('\\\\\n');
end
fprintf('\\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\\label{Tab:LABEL}\n');
fprintf('\\end{center}\n');
fprintf('\\end{table}\n');
