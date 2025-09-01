
close all
bigx=im_mosaic(squeeze(imHASTE(:,101:300,1:11,kim)),3,4,100);
bigx=imresize(bigx,[600 800],'bicubic');
close all
Y=im_tight(bigx,1,[0 0.8*max(bigx(:))])
seqname=cell2mat(fnp(kim+9));
load(seqname);
seq.read(strcat(seqname(1:end-3),'seq'))
tstart=0;
fig=seq.plot('TimeRange',[0 1],'timeDisp','ms','stacked',1);
if(strcmp(acqP.T2prep,'on')), k0=round((acqP.TEeff-acqP.TEprep+acqP.TE)/acqP.TE);
    else
        k0=round(acqP.TEeff/acqP.TE);
    end
figure
plot(acqP.flip,'o-'),axis([xlim 0 180]);
set(gcf,'Position',[100 100 400 300])
title('flip angles')
figure
set(gcf,'Position',[100 100 400 300])
plot(acqP.PEorder,'.-')
hold on, plot([k0 k0],[min(acqP.PEorder(:)) max(acqP.PEorder(:))],'k-')
title('phase encoding order')