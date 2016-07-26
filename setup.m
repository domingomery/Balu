clear all
close all
clc

BaluPath = pwd;

io = 'InputOutput';

addpath(fullfile(BaluPath));
addpath(fullfile(BaluPath,io));

I = imresize(imread('LogoBalu.png'),0.5);
imshow(I)
verBalu = 'Balu 4.0.1';
title(['Installing ' verBalu ' ...']);



d = dir;

n = length(d);

ft = Bio_statusbar(['Installing ' verBalu ]);

for i=1:n
    ft = Bio_statusbar(i/n,ft);
    st = d(i).name;
    if and(exist(st,'dir'),length(st)>2)
        pause(0.2)
        fprintf('Adding directory %s...\n',st);
        if strcmp(io,st)~=1
            addpath(fullfile(BaluPath,st));
        end
    end
end
delete(ft)
savepath

title(['Installing ' verBalu 'installed']);

fprintf('%s installed succefully!\n',verBalu);

pause(2)
close all

