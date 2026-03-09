% very basic and crude non-cartesian recon using griddata()
%%
disp('select rawdata')
[fnr, pnr] = uigetfile('*.dat', 'Select rawdata files','MultiSelect', 'on');
if(iscell(fnr)), nir=numel(fnr); else nir=1; fnr={fnr}; end
%%
disp('select protocol(s)')
[fnp, pnp] = uigetfile('*.mat', 'Select protocol files','MultiSelect', 'on');
if(iscell(fnp)), nip=numel(fnp); else nip=1; fnp={fnp}; end
%% Load datata
nuFFTflag=0;
FFTflag=1;
recnorm=0;r1=1;r2=5;ncut=0.5;
mp=0;
clear imsos imall ref_*
for kim=1
    clear imcoil acqP acq kt* raw* ref*
%% load data and parameters
    mp=mp+1
    cd(pnp)
    seqname=cell2mat([fnp(kim)]);
    load(seqname)
    acqP.nDummy=1; acqP.nrep=2; acqP.PEorder=[acqP.PEorder;acqP.PEorder];
    grappaflag=0; pocsflag=0;
    if(acqP.GR_fac>1), grappaflag=1; end
    if(acqP.HF_fac>0.5), pocsflag=1; end
    wfac=0.5;
    pocsflag=0;             % if you want to test reco without POCS
    %grappaflag=0;           % if you want to test reco without GRAPPA
    
    cd(pnr)
    rawname=cell2mat(fnr(:,kim));
    twix_obj = mapVBVD(rawname);
    %rawdata =double(twix_obj.image.unsorted());
    rawdata = double(twix_obj{2}.image.unsorted());

    %% resort rawdata
    raw0=permute(rawdata,[1 3 2]);
    
    siraw=size(raw0);
    nex=sum(acqP.nexc);
    echocount=acqP.necho;

    if(acqP.nDummy==1)
        raw=reshape(raw0,[siraw(1) echocount acqP.NSlices*acqP.nrep nex siraw(3)]);
    else
        raw=reshape(raw0,[siraw(1) echocount acqP.NSlices*acqP.nrep nex+1 siraw(3)]);
    end
    
    if(acqP.nDummy==0)
        rawp=raw(:,:,:,1,:);
        rawi=raw(:,:,:,2:end,:);
        %% process Dummy-scan
        [~,slorder]=sort(acqP.OSlices,'ascend');
        siD=size(rawp);
        ind=find(abs(rawp)==max(abs(rawp(:))));
        pos=mod(ind,siD(1));
        for ksl=1:acqP.NSlices
            tempr=sos(real(squeeze(rawp(:,:,ksl,1,:))));
            tempi=sos(imag(squeeze(rawp(:,:,ksl,1,:))));
            ref_raw(:,1:siD(2),ksl,mp)=tempr+1i*tempi;
            ref_angle(1:siD(2),ksl,mp)=angle(rawp(pos,:,ksl,mp));
        end
        ref_raw(:,:,:,mp)=ref_raw(:,:,slorder,mp);
        ref_angle(:,:,mp)=180/pi*ref_angle(:,slorder,mp);
    else
        rawi=raw;
    end
    if(strcmp(acq.cpmg_mod,'alt'))
        rawi(:,2:2:end,:,:,:)=-rawi(:,2:2:end,:,:,:);
    end
    %end    % used to extract Dummyscan without image reconstruction

    %% separate navigators
    if(sum(acqP.navmode)>0)
        rawnav=0*rawi(:,1,:,:,:);
    end
    
    if(acqP.navmode(2)>0),
        if(strcmp(acqP.T2prep,'on')), k0=round((acqP.TEeff-acqP.TEprep+acqP.TE)/acqP.TE);
        else
            k0=round(acqP.TEeff/acqP.TE);
        end
        rawnav(:,2,:,:,:)=rawi(:,k0,:,:,:); rawi=rawi(:,[1:k0-1 k0+1:end],:,:,:); acqP.PEorder=acqP.PEorder([1:k0-1 k0+1:end],:);acqP.necho=acqP.necho-1; end
    if(acqP.navmode(1)>0), rawnav(:,1,:,:,:)=rawi(:,1,:,:,:);rawi=rawi(:,2:end,:,:,:); end
    if(acqP.navmode(3)>0), rawnav(:,3,:,:,:)=rawi(:,end,:,:,:); rawi=rawi(:,1:end-1,:,:,:); end

    %% normalize data
    if(recnorm==1)
        ind=find(acqP.navmode>0);
        nnav=length(ind);
        [mph,amp]=pg_cpmg_f([1 0 0 0],1,acqP.flip,0,acqP.TE*1000,r1,r2,0);
        signorm=squeeze(abs(mph(2,1,1:acqP.necho-nnav)));
        ind=find(signorm<ncut);
        signorm(ind)=ncut;
        disp('normalize echoes')
        for k=1:acqP.necho-nnav
            rawi(:,k,:,:,:)=rawi(:,k,:,:,:)./signorm(k);
        end
    end
    if(recnorm==2)
        ind=find(acqP.navmode>0);
        nnav=length(ind);
        if(acqP.nDummy==0)
            signorm=zeros([acqP.necho-nnav 1])+ncut;
            temp=squeeze(rawp);
            ref=sos(sos(temp));
            peaks=ref(pos,:);
            peaks=peaks./max(peaks(:));
            signorm(1:length(peaks))=peaks;
            if(acqP.navmode(2)>0),
                signorm=signorm([1:acqP.navk-1 acqP.navk+1:end]);
            end
            if(acqP.navmode(1)>0),
                signorm=signorm(2:end);
            end
            if(acqP.navmode(3)>0),
                signorm=signorm([1:end-1]);
            end
            disp('normalize echoes')
            for k=1:acqP.necho-nnav
                rawi(:,k,:,:,:)=rawi(:,k,:,:,:)./signorm(k);
            end
        else
            disp('no reference scan available')
        end
    end
    %acqP.PEorder=temp;
    PEorder=acqP.PEorder(:);
    PEorder=PEorder-min(PEorder(:))+1;
    PEcount=acqP.PEcount(:);
    pemax=find(PEorder==0);
    %%

    PEall=floor(acqP.PEorder);
    if(acqP.navmode(2)>0), PEall=PEall([1:acqP.PEindex-1 acqP.PEindex+1:end],:); end
    if(acqP.navmode(1)>0), PEall=PEall(2:end,:); end
    if(acqP.navmode(3)>0), PEall=PEall(1:end-1,:); end
    siPEall=size(PEall);

    PEind=PEall+acqP.Ny/2;
    if(max(PEind(:))>acqP.Ny)
        siPE=max(PEind(:));
    else
        siPE=acqP.Ny;
    end
    
    nwav=length(acqP.wav);
    rawsl_tot=squeeze(zeros([siraw(1) siPE acqP.NSlices*acqP.nrep siraw(3) nwav]));
    rawis=permute(rawi,[1 2 4 3 5]);
    siraws=size(rawis);
    rawis=reshape(rawis,[siraws(1) siraws(2)*siraws(3) siraws(4) siraws(5)]);



    %% FT-reco
    PEind=PEind(1:siPEall(1)/acqP.nrep,:);
    ks=1;
    if(length(acqP.wav)==1)
        for ksl=1:acqP.NSlices
            for krep=1:acqP.nrep
                rawsl_tot(:,PEind(:),ks,:)=rawis(:,:,ks,:); ks=ks+1;
            end
        end
        rawslices=squeeze(rawsl_tot(:,end-acqP.Ny+1:end,:,:));
    else
        PEOS_reco;
    end
    %%
    if(grappaflag==1)
        for ksl=1:acqP.NSlices*acqP.nrep
            rawsl=squeeze(rawsl_tot(:,:,ksl,:));
            try
                GRAPPA_reco
                rawslices(:,:,ksl,:)=rawsl_rec;
            catch
                disp(strcat('GRAPPA didnt work for :',fnp(kim),', slice',num2str(ksl)))
            end
        end
    end

    %%
    if(pocsflag==1)
        for ksl=1:acqP.NSlices*acqP.nrep
            rawsl=squeeze(rawslices(:,:,ksl,:));
            rawsl=permute(rawsl,[3 1 2]);
            rawsl(:,:,1)=0;
            try
                [im, rawpocs] = pocs(rawsl);
                rawpocs=permute(rawpocs,[2 3 1]);
                rawslices(:,:,ksl,:)=rawpocs;
            catch
                disp(strcat('pocs didnt work for ',fnp(kim),', slice',num2str(ksl)))
            end
        end
    end
    %%
    zf=0;
    rawtemp=rawslices;
    sirs=size(rawslices);
    if(zf==1)

        rawzf=zeros([sirs(1) sirs(1) sirs(3) sirs(4)]);
        rawzf(:,floor((sirs(1)-sirs(2))/2)+1:floor((sirs(1)-sirs(2))/2)+sirs(2),:,:)=rawtemp;
        rawslices=rawzf;
    end
    if(length(sirs)==3), clear temp; temp(:,:,1,:)=rawslices; rawslices=temp; sirs=size(rawslices); end;
    %% FFT-reco
    PEindex=sort(PEorder(:));
    if(PEindex(1)==0), PEindex=PEindex+1; end
    
    if(FFTflag==1)
        for ksl=1:acqP.NSlices*acqP.nrep
            for kcoil=1:sirs(4)
                %imcoil(:,:,kcoil)=squeeze(fftshift(ifft2(fftshift(filt.*squeeze(rawslices(:,:,ksl,kcoil))))));
                imcoil(:,:,kcoil)=squeeze(fftshift(ifft2(fftshift(squeeze(rawslices(:,:,ksl,kcoil))))));
            end
            temp=sos(squeeze(imcoil(:,:,:)));
            imsos(:,:,ksl)=temp(end:-1:1,:)';
        end
    end
    %%
    [~,slorder]=sort(acqP.OSlices,'ascend');
    imsos1=0*imsos;
    for k=1:acqP.nrep
        imsos1(:,:,k:acqP.nrep:end)=imsos(:,:,acqP.nrep*slorder-acqP.nrep+k);
    end
    imall(:,:,:,mp)=imsos1(:,:,:);

    if(nuFFTflag==1)
        myTSE_reco_nufft
        imall_n(:,:,:,mp)=imsos_nuFFT;
    end
end
%% extract navigators
if(sum(acqP.navmode)>0),
    sipro=size(rawnav);
    projnavc=0*rawnav;projnav=zeros([sipro(1:4)]);
    for knav=1:sipro(2)
        for ksl=1:acqP.NSlices*acqP.nrep
            for kproj=1:sipro(4)
                for kcoil=1:sipro(5)
                    temp=squeeze(rawnav(:,knav,ksl,kproj,kcoil));
                    projnavc(:,knav,ksl,kproj,kcoil)=fftshift(fft(temp));
                end
                projr=sos(real(projnavc(:,knav,ksl,kproj,:)));
                proji=sos(imag(projnavc(:,knav,ksl,kproj,:)));
                projnav(:,knav,ksl,kproj)=projr+i*proji;
            end
        end
        projnav(:,knav,slorder,:)=projnav(:,knav,:,:);
    end
end

%% sort images from interleaved acquisition
temp1=imsos(:,:,1:acqP.NSlices);
temp2=imsos(:,:,acqP.NSlices+1:2*acqP.NSlices);
imall=temp1(:,:,slorder);
imall(:,:,acqP.NSlices+1:2*acqP.NSlices)=temp2(:,:,slorder);
%%

script nufftreco
x=1
