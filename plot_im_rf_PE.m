close all
plotflag='0011';
pnc=pwd;
for kim=1:4

    bigx=im_mosaic(squeeze(imTSEfast(:,:,:,kim)),3,4,150);
    %bigx=imresize(bigx,[600 800],'bicubic');
    if (plotflag(1)=='1')
        close
        Y=im_tight(bigx,0.4,[0 1*max(bigx(:))])
    end
    cd(pnp)
    seqname=cell2mat(fnp(kim));
    load(seqname);
    seq.read(strcat(seqname(1:end-3),'seq'))
    cd(pnc)
    tstart=0;
    if (plotflag(2)=='1')
        fig=seq.plot('TimeRange',[0.1 0.2],'timeDisp','ms','stacked',1);
        set(gcf,'Position',[100 100 350 240])
    end
    if(strcmp(acqP.T2prep,'on')), k0=round((acqP.TEeff-acqP.TEprep+acqP.TE)/acqP.TE);
    else
        k0=round(acqP.TEeff/acqP.TE);
    end
    if (plotflag(3)=='1')
        figure
        plot(acqP.flip,'o-'),axis([xlim 0 180]);
        set(gcf,'Position',[100 100 400 300])
        title('flip angles')
    end
    if (plotflag(4)=='1')
        figure
        set(gcf,'Position',[100 100 320 240])
        plot(acqP.PEorder,'.-')
        hold on, plot([k0 k0],[min(acqP.PEorder(:)) max(acqP.PEorder(:))],'k-')
        title('phase encoding order')
    end
end