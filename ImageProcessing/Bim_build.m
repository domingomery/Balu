function bim = Bim_build(bdir,bext,bresize)

bfiles = dir([bdir '*.' bext]);

bim.path          = bdir;
bim.prefix        =  '*';   
bim.extension     =  bext;
bim.imgmin        = 1;
bim.imgmax        = length(bfiles);
bim.show          = 1;
if exist('bresize','var')
    bim.resize = bresize;
end
end