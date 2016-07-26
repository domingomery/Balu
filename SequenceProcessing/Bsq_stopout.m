function [Yfilt,ifilt,iw] = Bsq_stopout(Y,options)

voc         = options.voc;
i_go        = options.i_go;
kd_voc      = options.kd_voc;
kd_go       = options.kd_go;
stopth      = options.stopth;

iy          = (1:size(Y,1))';
ifilt       = false(size(Y,1),1);
[iw,dist]   = vl_kdtreequery(kd_voc,voc,Y','NumNeighbors',1);
id          = sort(dist);
stopth_val  = id(round(stopth*length(id)));
ii          = dist<stopth_val;
Y0          = Y(ii,:);
iy0         = iy(ii);
i0          = iw(ii);
[~,dist_go] = vl_kdtreequery(kd_go,i_go,double(i0),'NumNeighbors',1);
Yfilt       = Y0(dist_go==0,:);
iy2         = iy0(dist_go==0);
ifilt(iy2)  = true;







