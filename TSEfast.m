%% Create a TSE sequence and export for execution
%
% For descroption of parameters see manual

%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m for slew. Unspecificied hardware
% parameters will be assigned default values.
tic
runseq=0;
clear count* acqP acq
plotflag='100001';
acq.dG=100e-6; acq.dGs=100e-6;
pulseflag_ex=0;
pulseflag_ref=0;
%%
system = mr.opts('MaxGrad', 45, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);
sysphas = mr.opts('MaxGrad', 70, 'GradUnit', 'mT/m', ...
    'MaxSlew', 125, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);
B0=3;
system.maxB1=2000;
%%
% A new sequence object is created by calling the class constructor.
nseq=22;                                                                     %  comment for interleaved acquisition
for kseq=12:21                                                               %  comment for interleaved acquisition
    seq=mr.Sequence(system);                                                 %  comment for interleaved acquisition

    if(runseq==1), load('TSE_.mat');
    else

        acqP.T2prep='off';
        acqP.fatsat='off';
        acqP.flipback='off';
        acqP.sat_ppm=-3.45;
        acqP.sat_freq=acqP.sat_ppm*1e-6*B0*system.gamma;
        %% Sequence events
        % Some sequence parameters are defined using standard MATLAB variables
        acqP.fov=200e-3;
        acqP.Nx=200;        acqP.Ny=200;        acqP.PEfac=acqP.Nx/acqP.Ny;     %=1:same fov; =acqP.Nx/acqP.Ny:same res
        acqP.necho=200;
        acqP.samplingTime= 2.8e-3;                                              %should be a multiple of Nx in microseconds
        acqP.NSlices=21;     acqP.sliceGAP=3;  acqP.sliceThickness=4e-3;
        acqP.flipref=40;   acqP.flipflag=4;
        acqP.TE=4e-3;      acqP.TEeff=80e-3;   acqP.TEprep=30e-3; acqP.dTE=0.0e-3;
        acqP.TR=16;       acqP.TI=0e-3;       acqP.TImod=1;
        acqP.HF_fac=0.9;      acqP.PEref=20;      acqP.GR_fac=2;
        acqP.PEtype='centric';   acqP.PEover=0;
        acqP.nDummy=0;  % =0: dummy scan will be performed

        % acqP-parameters related to special acquisition modes.
        acqP.GRramp=0;
        acqP.nrep=1;
        acqP.navmode=[0 0 0];  acqP.navk=5;                                     % navigator echoes at beginning and end of echotrain
        acqP.wav=[1];    acqP.wavmode='cont';    acqP.PEinc=0;
        acqP.GDs=0;
        acqP.GDx=0;
        acqP.GDy=0;
        acqP.tD=1e-3;
        acq.dGD=400e-6;
        %%  acq-Parameters
        acq.accfac=1;
        acq.sigpy='off';
        acq.cpmg_mod='const';
        acq.tSpS=0.3e-3;
        acq.tEx=0.7e-3;
        acq.tBwPex=2;
        acq.ExApo=0.25;
        tExwd=acq.tEx+system.rfRingdownTime+system.rfDeadTime;
        acq.nPinit=1;
        acq.tRef=0.5e-3;
        acq.tBwPref=0.5;
        acq.RefApo=0.4;
        tRefwd=acq.tRef+system.rfRingdownTime+system.rfDeadTime;
        acq.tRefSE=4e-3;
        acq.tBwPrefSE=8;
        tRefwdSE=acq.tRefSE+system.rfRingdownTime+system.rfDeadTime;
        acq.rfTIdur=5e-3;
        acq.rfTItBwP=8;
        acq.GSexFac=1;
        acq.GSrefFac=1;
        acq.GSSEfac=1;
        acq.GSIRFac=1;
        acq.fspR=2;
        acq.fspS=3;
        acq.spoilermode='amp';
        %slice order factor
        if(nseq>1), myTSE_para;  end
    end
    %% Calculate dependant parameters
    if(strcmp((acqP.PEtype),'paired')&&(acqP.Ny/acqP.necho>1.5))
        disp('paired PE-mode only valid for single-shot acquisition, sequence construction terminated')
        return
    end

    if(strcmp(acqP.T2prep,'on')), k0=round((acqP.TEeff-acqP.TEprep+acqP.TE)/acqP.TE);
    else
        k0=round(acqP.TEeff/acqP.TE);
    end
    if(acqP.navmode(1)>0), k0=k0-1; end
    if(acqP.navmode(2)>0), if(acqP.navk<k0),k0=k0-1; end; end
    if(k0<1)
        disp('TEeff was adpated to include navigator echoes')
        k0=1;
    end
    if(strcmp(acqP.T2prep,'on')&&(acqP.navmode(1)==2))
        disp('PE-navigator in first echo is not compatible with acqP.T2prep=on, navigator will be in read direction')
        acqP.navmode(1)=0;
    end
    tExwd=acq.tEx+system.rfRingdownTime+system.rfDeadTime;
    tRefwd=acq.tRef+system.rfRingdownTime+system.rfDeadTime;
    tRefwdSE=acq.tRefSE+system.rfRingdownTime+system.rfDeadTime;
    tSpex=0.5*(acqP.TE-tExwd-tRefwd)+acqP.dTE;
    readoutTime=acqP.samplingTime+2*system.adcDeadTime;
    acq.tSpR=0.5*(acqP.TE-readoutTime-tRefwd);
    acqP.samplingTime = readoutTime - 2*system.adcDeadTime;

    slofac=0;                                                               % determines direction of slice order
    acqP.OSlices=-[1:2:acqP.NSlices 2:2:acqP.NSlices]+ceil(acqP.NSlices/2);
    if(slofac==1)
        acqP.OSlices=acqP.OSlices(end:-1:1);
    end
    if(acqP.tD==0),     acqP.GDs=0;    acqP.GDx=0;    acqP.GDy=0; acqP.tD=10*acq.dG; end
    if strcmp(acq.cpmg_mod,'alt')
        rfex_phase=0;rfref_phase=pi; rfref_phase0=pi;
    else
        rfex_phase=pi/2; rfref_phase=0;
    end
    %%
    %%% Base gradients
    %%% Slice selection
    %
    % First, the slice selective RF pulses (and corresponding slice gradient)
    % are generated using the |makeSincPulse| function.
    % Gradients are recalculated such that their flattime covers the pulse plus
    % the rfdead- and rfringdown- times.
    % %
    flipex=90*pi/180;

    [rfex, gzex] = mr.makeSincPulse(flipex,system,'Duration',acq.tEx,...
        'sliceThickness',acqP.sliceThickness,'apodization',acq.ExApo,'timeBwProduct',acq.tBwPex,'PhaseOffset',rfex_phase,'use','excitation');
    rfex.delay=system.rfRingdownTime;
    GSex = mr.makeTrapezoid('z',system,'amplitude',gzex.amplitude,'FlatTime',tExwd,'riseTime',acq.dG);
    GSex.amplitude=acq.GSexFac*GSex.amplitude;
    GSex.area=acq.GSexFac*GSex.area;
    %if(plotflag(1)=='1'), plotPulse(rfex,GSex); end
    GSexfirst=GSex.amplitude;
    acq.GSex=GSex.amplitude;

    [rfref0, gzref] = mr.makeSincPulse(pi,system,'Duration',acq.tRef,...
        'sliceThickness',acqP.sliceThickness,'apodization',0,'timeBwProduct',acq.tBwPref,'PhaseOffset',rfref_phase,'use','refocusing');
    rfref.delay=system.rfRingdownTime;
    GSref = mr.makeTrapezoid('z',system,'amplitude',gzref.amplitude,'FlatTime',tRefwd,'riseTime',acq.dG);
    GSref.amplitude=acq.GSrefFac*GSref.amplitude;
    GSref.area=acq.GSrefFac*GSref.area;
    acq.GSref=GSref.amplitude;
    refenvelope0=rfref0.signal;

    [rfref, gzref] = mr.makeSincPulse(pi,system,'Duration',acq.tRef,...
        'sliceThickness',acqP.sliceThickness,'apodization',acq.RefApo,'timeBwProduct',acq.tBwPref,'PhaseOffset',rfref_phase,'use','refocusing');
    rfref.delay=system.rfRingdownTime;
    GSref = mr.makeTrapezoid('z',system,'amplitude',gzref.amplitude,'FlatTime',tRefwd,'riseTime',acq.dG);
    GSref.amplitude=acq.GSrefFac*GSref.amplitude;
    GSref.area=acq.GSrefFac*GSref.area;
    acq.GSref=GSref.amplitude;
    refenvelope=rfref.signal;
    %%
    if(strcmp(acq.sigpy,'on'))
        [rfrefSE, gzrefSE] = mr.makeSLRpulse(pi,'Duration',acq.tRefSE,'PhaseOffset',rfref_phase,...
            'SliceThickness',acqP.sliceThickness,'timeBwProduct',0.5*acq.tBwPrefSE,'passbandRipple',1,'stopbandRipple',1e-2,'filterType','ms','system',system,'use','refocusing');

    else
        [rfrefSE, gzrefSE] = mr.makeSincPulse(pi,system,'Duration',acq.tRefSE,'PhaseOffset',rfref_phase,...
            'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',acq.tBwPrefSE,'use','refocusing');
    end
    GSrefSE = mr.makeTrapezoid('z',system,'amplitude',gzref.amplitude,'FlatTime',tRefwdSE,'riseTime',acq.dG);
    GSrefSE.amplitude=acq.GSSEfac*GSrefSE.amplitude;
    GSrefSE.area=acq.GSSEfac*GSrefSE.area;
    rfrefSE.delay=system.rfDeadTime;
    %if(plotflag(1)=='1'), plotPulse(rfref,GSref); end
    GSrefirst=GSref.amplitude;
    GSSErefirst=GSrefSE.amplitude;
    AGSex=GSex.area/2;
    AGSref=GSref.area/2;
    if(strcmp(acq.spoilermode,'amp'))
        GSspr = mr.makeTrapezoid('z',system,'amplitude',GSref.amplitude*acq.fspS,'duration',acq.tSpS,'riseTime',acq.dG);
    else
        GSspr = mr.makeTrapezoid('z',system,'area',GSref.area*acq.fspS,'duration',acq.tSpS,'riseTime',acq.dG);
    end
    GSspex = mr.makeTrapezoid('z',system,'area',GSspr.area-AGSex,'duration',tSpex,'riseTime',acq.dG);
    GSspr_end = mr.makeTrapezoid('z',system,'amplitude',1e6,'duration',4e-3,'riseTime',5*acq.dG);
    GSspr_end.delay=5*acq.dG;
    %% Fatsat
    if(strcmp(acqP.fatsat,'on'))
        rf_fst=4e-3;
        if (B0<2), rf_fst=1e-5*floor(1e5*10e-3*1.5/B0); end
        rf_fs = mr.makeGaussPulse(110*pi/180,'system',system,'Duration',rf_fst,...
            'bandwidth',abs(acqP.sat_freq),'freqOffset',acqP.sat_freq,'use','other');
        gz_fs = mr.makeTrapezoid('z',system,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm
        fs_dur=gz_fs.delay+gz_fs.riseTime+gz_fs.flatTime+gz_fs.fallTime;
    else
        fs_dur=0;
    end

    %% Inversion recovery
    if(acqP.TI>0)

        if(strcmp(acq.sigpy,'on'))
            [rfIR,gzIR]=mr.makeAdiabaticPulse('hypsec','sliceThickness',acqP.sliceThickness,'duration',acq.rfTIdur,...
                'bandwidth',acq.tBwPex/(1000*acq.tEx),'dwell',1e-6,'beta',1600, 'mu',2.5,'use','inversion','system',system);
        else
            [rfIR, gzIR] = mr.makeSincPulse(pi,system,'Duration',acq.rfTIdur,...
                'sliceThickness',acqP.sliceThickness,'apodization',0.5,'timeBwProduct',acq.rfTItBwP,'PhaseOffset',0,'delay',system.rfDeadTime,'use','other');
        end
        rfIR.delay=rfIR.delay+acq.dG;
        GSIR = mr.makeTrapezoid('z',system,'amplitude',gzIR.amplitude,'FlatTime',rfIR.shape_dur+system.rfRingdownTime+system.rfDeadTime,'riseTime',acq.dG);
        GSIR.amplitude=acq.GSIRFac*GSIR.amplitude;
        GSIRsp=mr.makeTrapezoid('z',system,'amplitude',0.5*system.maxGrad,'duration',4e-3,'riseTime',acq.dG,'delay',mr.calcDuration(GSIR)-acq.dG);
        GSIRtot=mr.addGradients({GSIR,GSIRsp},'system',system);
        GSIRtot_times=[0; GSIRtot.tt+acq.dG; GSIRtot.shape_dur+2*acq.dG];
        GSIRtot_gradients=[0; GSIRtot.waveform; 0];
        GSIRtot=mr.makeExtendedTrapezoid('z','times',GSIRtot_times,'amplitudes',GSIRtot_gradients);
        TI_dur=GSIRtot.shape_dur; TI_t2P=rfIR.delay+rfIR.shape_dur/2;
        acqP.nDummy=0;
    else
        TI_dur=0;

    end



    %%
    %%% Readout gradient
    % To define the remaining encoding gradients we need to calculate the
    % $k$-space sampling. The Fourier relationship
    %
    % $$\Delta k = \frac{1}{acqP.fov}$$
    %
    % Therefore the area of the readout gradient is $n\Delta k$.
    deltak=1/acqP.fov;
    kWidth = acqP.Nx*deltak;

    GRacq = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime,'riseTime',acq.dG);
    facs=factor(round(acqP.samplingTime*10^6));facs=facs(end:-1:1);
    ind=find(cumprod(facs)>=acqP.Nx*acq.accfac*system.adcSamplesDivisor);
    npoints=prod(facs(1:ind(1)));
    adc = mr.makeAdc(npoints,'Duration',acqP.samplingTime, 'Delay', system.adcDeadTime);%,'Delay',GRacq.riseTime);
    adcSE=adc;
    if(strcmp(acq.spoilermode,'amp'))
        GRspr = mr.makeTrapezoid('x',system,'amplitude',GRacq.amplitude*acq.fspR*(1+acqP.GRramp),'duration',acq.tSpR,'riseTime',acq.dG);
    else
        GRspr = mr.makeTrapezoid('x',system,'area',GRacq.area*acq.fspR*(1+acqP.GRramp),'duration',acq.tSpR,'riseTime',acq.dG);
    end
    GRsprSE = mr.makeTrapezoid('x',system,'amplitude',GRacq.amplitude*acq.fspR*(1+acqP.GRramp),'duration',acq.tSpS,'riseTime',acq.dG);
    AGRspr=GRspr.area;%GRacq.area/2*acq.fspR;




    %%
    %%% Phase encoding
    % To move the $k$-space trajectory away from 0 prior to the readout a
    % prephasing gradient must be used. Furthermore rephasing of the slice
    % select gradient is required.

    nex=floor(acqP.Ny/acqP.necho);

    Ny0=acqP.Ny;
    acqP.PEorder=[];acqP.PEcount=[];acqP.PEindex=[];
    clear nexc;

    nav=length(acqP.wav);
    for kav=1:nav
        %  for kav=1
        [PEorder,PEindex,NS] = PE_order(acqP.PEtype,floor(acqP.Ny*acqP.wav(kav)),acqP.necho,k0,[],acqP.HF_fac*acqP.wav(kav),acqP.PEref,acqP.GR_fac,acqP.PEover);
        if(NS==1), acqP.necho=length(PEindex); end
        PEcount=0*PEorder+kav;
        try
            acqP.PEorder=[acqP.PEorder PEorder]+acqP.PEinc*kav;
        catch
            disp(strcat('warning ! ETL of ',num2str(kav),'th average too short'));
            siP=size(acqP.PEorder); siPE=size(PEorder); PEnew=zeros([siP(1) siPE(2)]);
            PEnew(1:siPE(1),1:siPE(2))=PEorder; PEorder=PEnew;
            PEcount=0*PEorder+kav;
            acqP.PEorder=[acqP.PEorder PEorder]+acqP.PEinc*kav;
        end
        acqP.PEcount=[acqP.PEcount PEcount];
        acqP.PEindex=[acqP.PEindex acqP.PEindex];
        nexc(kav)=NS;
    end

    nex=sum(nexc);
    if(nav>1)
        if(strcmp(acqP.wavmode,'cont'))
            [PEtemp,ind]=sort(acqP.PEorder(:));
            count=acqP.PEcount(:);
            count=count(ind);
            [PEorder,PEindex,NS] = PE_order(acqP.PEtype,floor(acqP.Ny*acqP.wav(kav)),acqP.necho,k0,PEtemp',acqP.HF_fac*acqP.wav(kav),acqP.PEref,acqP.GR_fac,acqP.PEover);
            acqP.PEorder=PEorder;
            size(acqP.PEorder);
            acqP.PEindex=PEindex;
            [~,PEind]=sort(acqP.PEorder(:));
            [~,PEindb]=sort(PEind);
            acqP.PEcount=count(PEindb);
        end
    end
    acqP.PEorder0=acqP.PEorder;
    if(acqP.navmode(1)>0)
        siPE=size(acqP.PEorder);
        temp=zeros([siPE(1)+1 siPE(2)]);
        temp(2:siPE(1)+1,:)=acqP.PEorder;
        acqP.PEorder=temp;
        acqP.necho=acqP.necho+1;
        k0=k0+1;
    end
    if(acqP.navmode(2)>0)
        siPE=size(acqP.PEorder);
        temp=zeros([siPE(1)+1 siPE(2)]);
        temp(1:acqP.navk-1,:)=acqP.PEorder(1:acqP.navk-1,:);
        temp(acqP.navk+1:end,:)=acqP.PEorder(acqP.navk:end,:);
        acqP.PEorder=temp;
        acqP.necho=acqP.necho+1;
        k0=k0+1;
    end
    if(acqP.navmode(3)>0)
        siPE=size(acqP.PEorder);
        acqP.PEorder(siPE(1)+1,siPE(2))=0;
        acqP.necho=acqP.necho+1;
    end
    if(plotflag(1)=='1')
        figure
        hold on
        plot(acqP.PEorder,'.-')
        temp=ylim;
        plot([k0 k0],[temp(1) 0.8*temp(2)],'k-')
        text(k0,0.9*temp(2),'k0')
        title(strcat('PE-order (k vs. Necho) ,',acqP.PEtype))
    end
    acqP.nexc=nexc;
    temp=[];
    for k=1:acqP.nrep
        temp=[temp;acqP.PEorder];
    end
    acqP.PEorder=temp;
    phaseAreas= acqP.PEorder*deltak*acqP.PEfac;
    acq.phaseAreas=phaseAreas;

    %% split gradients and recombine into blocks
    % lets start with slice selection....
    GS1times=[0 GSex.riseTime];
    GS1amp=[0 GSexfirst];
    GS1 = mr.makeExtendedTrapezoid('z','times',GS1times,'amplitudes',GS1amp);
    GS1rew_amp=[GSexfirst 0];
    GS1rew = mr.makeExtendedTrapezoid('z','times',GS1times,'amplitudes',GS1rew_amp);

    GS2times=[0 GSex.flatTime];
    GS2amp=[GSexfirst GSexfirst];
    GS2 = mr.makeExtendedTrapezoid('z','times',GS2times,'amplitudes',GS2amp);

    GS3times=[0 GSspex.riseTime GSspex.riseTime+GSspex.flatTime GSspex.riseTime+GSspex.flatTime+GSspex.fallTime];
    GS3amp=[GSexfirst GSspex.amplitude GSspex.amplitude GSrefirst];
    GS3 = mr.makeExtendedTrapezoid('z','times',GS3times,'amplitudes',GS3amp);

    GS3rewamp=[GSrefirst GSspex.amplitude GSspex.amplitude GSexfirst];
    GS3rew = mr.makeExtendedTrapezoid('z','times',GS3times,'amplitudes',GS3rewamp);

    GS4times=[0 GSref.flatTime];
    GS4amp=[GSrefirst GSrefirst];
    GS4 = mr.makeExtendedTrapezoid('z','times',GS4times,'amplitudes',GS4amp);

    GS4SEtimes=[0 GSrefSE.flatTime];
    GS4SEamp=[GSSErefirst GSSErefirst];
    GS4SE = mr.makeExtendedTrapezoid('z','times',GS4SEtimes,'amplitudes',GS4SEamp);

    GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5amp=[GSrefirst GSspr.amplitude GSspr.amplitude 0];
    GS5SEtimes=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5SEamp=[GSSErefirst GSspr.amplitude GSspr.amplitude 0];

    if(GSspr.flatTime==0)
        GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.fallTime];
        GS5amp=[GSrefirst GSspr.amplitude  0];
        GS5SEtimes=[0 GSspr.riseTime GSspr.riseTime+GSspr.fallTime];
        GS5SEamp=[GSSErefirst GSspr.amplitude 0];
    end
    GS5 = mr.makeExtendedTrapezoid('z','times',GS5times,'amplitudes',GS5amp);
    GS5SE = mr.makeExtendedTrapezoid('z','times',GS5SEtimes,'amplitudes',GS5SEamp);
    GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS7amp=[0 GSspr.amplitude GSspr.amplitude GSrefirst];
    GS7SEamp=[0 GSspr.amplitude GSspr.amplitude GSSErefirst];
    GS7SEtimes=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];

    if(GSspr.flatTime==0)
        GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
        GS7amp=[0 GSspr.amplitude GSrefirst];
        GS7SEtimes=[0 GSspr.riseTime GSspr.riseTime+GSspr.fallTime];
        GS7SEamp=[0 GSspr.amplitude GSSErefirst];
    end


    GS7 = mr.makeExtendedTrapezoid('z','times',GS7times,'amplitudes',GS7amp);
    GS7on=GS7;
    GS7SE = mr.makeExtendedTrapezoid('z','times',GS7SEtimes,'amplitudes',GS7SEamp);
    GSarea=calcArea(GS4)+calcArea(GS5)+calcArea(GS7);
    if((max(acqP.navmode)==3))
        GSpredur=system.gradRasterTime*floor(readoutTime/4/system.gradRasterTime);
        GSread = mr.makeTrapezoid('z',system,'Amplitude',5e5,'Duration',readoutTime-2*GSpredur,'riseTime',acq.dG);
        GSprew=mr.makeTrapezoid('z',system,'Area',-GSread.area/2,'Duration',GSpredur,'riseTime',acq.dG);
        GSacqtimes=[0 acq.dG GSpredur-acq.dG GSpredur GSpredur+acq.dG GSpredur+acq.dG+GSread.flatTime readoutTime-GSpredur readoutTime-GSpredur+acq.dG readoutTime-acq.dG readoutTime]
        GSacqamp=[0 GSprew.amplitude GSprew.amplitude 0 GSread.amplitude GSread.amplitude 0 GSprew.amplitude GSprew.amplitude 0 ]
        %plot(GSacqtimes,GSacqamp)
        GSacq = mr.makeExtendedTrapezoid('z','times',GSacqtimes,'amplitudes',GSacqamp);
    end

    if(strcmp(acqP.T2prep,'on'))
        GS3p1=mr.makeExtendedTrapezoid('z','times',[0 GSex.fallTime],'amplitudes',[GSex.amplitude 0]);
        GS3p2=mr.makeTrapezoid('z',system,'area',-AGSex+GSarea/4,'duration',tSpex,'riseTime',acq.dG);
        % and now the readout gradient....
        GS5p1 = mr.makeTrapezoid('z',system,'Area',GSarea/4,'Duration',acq.tSpS);
    end

    GR5times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
    GR5amp=[0 GRspr.amplitude GRspr.amplitude GRacq.amplitude*(1+acqP.GRramp)];
    GR7times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
    GR7amp=[GRacq.amplitude*(1+acqP.GRramp) GRspr.amplitude GRspr.amplitude 0];
    if(GRspr.flatTime==0)
        GR5times=[0 GRspr.riseTime GRspr.riseTime+GRspr.fallTime];
        GR5amp=[0 GRspr.amplitude GRacq.amplitude*(1+acqP.GRramp)];
        GR7times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
        GR7amp=[GRacq.amplitude*(1+acqP.GRramp) GRspr.amplitude 0];
    end


    GR5 = mr.makeExtendedTrapezoid('x','times',GR5times,'amplitudes',GR5amp);

    GR6times=[0 readoutTime/2 readoutTime];
    GR6amp=[GRacq.amplitude*(1+acqP.GRramp) GRacq.amplitude*(1-acqP.GRramp) GRacq.amplitude*(1+acqP.GRramp)];

    GR6 = mr.makeExtendedTrapezoid('x','times',GR6times,'amplitudes',GR6amp);
    GR6on=GR6;

    GR7 = mr.makeExtendedTrapezoid('x','times',GR7times,'amplitudes',GR7amp);
    GR7on=GR7;

    GRarea=calcArea(GR5)+calcArea(GR6)+calcArea(GR7);
    GRpreph = mr.makeTrapezoid('x',system,'Area',GRarea/2,'duration',tSpex,'riseTime',acq.dGs);
    GR3=GRpreph;
    if(any(acqP.navmode(:)>1))
        try
            GRrew = mr.makeTrapezoid('x',system,'area',GRarea/2,'duration',GRspr.flatTime+2*acq.dG,'riseTime',acq.dG);
        catch
            ('only navigators in read out direction allowed')
            acqP.navmode=ceil(acqP.navmode/3);
        end
    end

    GRpre_s=mr.makeTrapezoid('x',system,'Area',GRarea/2,'duration',tSpex,'riseTime',2*acq.dG);
    GRpre_s.delay=5*acq.dG;
    % and now the PE-navigator
    GPacq = mr.makeTrapezoid('y',system,'FlatArea',kWidth,'FlatTime',readoutTime,'riseTime',acq.dG);
    if(any(acqP.navmode(:)==2))
        GPRpre = mr.makeTrapezoid('y',system,'Area',-GPacq.area/2,'Duration',GRspr.flatTime+2*acq.dG,'riseTime',acq.dG);
    end
    nPEecho=acqP.necho;

    %% Calculate filltimes
    acq.Extime=mr.calcDuration(GS1)+mr.calcDuration(GS2)+mr.calcDuration(GS3);
    acq.Reftime=mr.calcDuration(GS4)+mr.calcDuration(GS5)+mr.calcDuration(GS7)+readoutTime;
    tEx2P=rfex.delay+rfex.shape_dur/2;
    if(acqP.flipback==0)
        tend=mr.calcDuration(GS4)+mr.calcDuration(GS5)+mr.calcDuration(GSspr_end);
    else
        tend=mr.calcDuration(GS2)+mr.calcDuration(GS3)+mr.calcDuration(GS4)+mr.calcDuration(GS1rew)+mr.calcDuration(GSspr_end);
    end
    tETrain=acq.Extime+acqP.necho*acq.Reftime+tend+fs_dur+TI_dur;
    if(strcmp(acqP.T2prep,'on'))
        tETrain=tETrain+acqP.TEeff;
    end
    if(acqP.NSlices*tETrain>(acqP.TR)),
        acqP.TR=acqP.NSlices*tETrain;
        disp(strcat('Warning!!! acqP.TR too short, adapted to include all slices to : ',num2str(acqP.TR*1000),' ms'));
    end
    TR_tot=(acqP.TR-acqP.NSlices*tETrain)/acqP.NSlices;
    TIfill=0;
    IRSlices=acqP.OSlices;
    if(acqP.TI>0)
        nTI=floor(acqP.TI/(tETrain+TR_tot));
        IRSlices=circshift(acqP.OSlices,-2);
        if(nTI==0), acqP.TImod=1; end
        if(acqP.TImod==0)
            acqP.TI=nTI*tETrain;
            disp(strcat('acqP.TI adapted t0 : ',num2str(acqP.TI*1000),' ms'));
        else
            TIfill=acqP.TI-nTI*tETrain;
            if(TIfill>TR_tot), TR_tot=(acqP.TI-nTI*tETrain)/(nTI+1);
                TIfill=TR_tot-tEx2P-TI_dur+TI_t2P;                  % TIfill is corrected for pulse timing
                acqP.TR=acqP.NSlices*(tETrain+TR_tot);
                disp(strcat('acqP.TR adapted t0 : ',num2str(acqP.TR*1000),' ms'));
            end
        end
    end
    TRfill=TR_tot-TIfill;

    % round to gradient raster
    TIfill=system.gradRasterTime * round(TIfill / system.gradRasterTime);
    TRfill=system.gradRasterTime * round(TRfill / system.gradRasterTime);
    delayTR = mr.makeDelay(TRfill);
    delayTI = mr.makeDelay(TIfill);
    %% Fill times for SE-TSE
    if(strcmp(acqP.T2prep,'on'))
        TEfill1=acqP.TEprep/2-GS2.shape_dur/2-GS3p1.shape_dur-mr.calcDuration(GS3p2)-mr.calcDuration(GRsprSE)-GS4SE.shape_dur/2;
        TEfill2=acqP.TEprep/2-GS5SE.shape_dur-GS4SE.shape_dur/2-GS4.shape_dur/2-acqP.TE/2+tRefwd;

        if(TEfill1<=0),
            fprintf(2,'\n\nTEprep too short, sequence creation terminated\n')
            return
        end
        if(TEfill1<acqP.tD), acqP.tD=TEfill1-acq.dG;
            disp(strcat('Warninng, acqP.tD has been reduced to:',num2str(acqP.tD*1000),' ms'))
        end

        delayTE1 = mr.makeDelay(TEfill1);
        delayTE2 = mr.makeDelay(TEfill2);
    end
    %% Diffusion gradients
    if(strcmp(acqP.T2prep,'on'))
        GDs=mr.makeTrapezoid('z',system,'amplitude',acqP.GDs/1000*system.gamma,'duration',acqP.tD,'riseTime',acq.dGD,'delay',acq.dG);
        GDx=mr.makeTrapezoid('x',system,'amplitude',acqP.GDx/1000*system.gamma,'duration',acqP.tD,'riseTime',acq.dGD,'delay',acq.dG);
        GDy=mr.makeTrapezoid('y',system,'amplitude',acqP.GDy/1000*system.gamma,'duration',acqP.tD,'riseTime',acq.dGD,'delay',acq.dG);
        GDxr=GDx; GDxr.delay=TEfill2-acqP.tD-acq.dG;
        GDyr=GDy; GDyr.delay=TEfill2-acqP.tD-acq.dG;
        GDsr=GDs; GDsr.delay=TEfill2-acqP.tD-acq.dG;
    end

    %% calculate flip angles
    if(acqP.flipflag<=1)
        rf=acqP.flipref+zeros([1 acqP.nrep.*acqP.necho]); end

    if(acqP.flipflag==1),  rf(1)=90+acqP.flipref/2;end

    if(acqP.flipflag==2),
        kz=k0;
        if(kz<5),kz=5; end;
        flipend=acqP.flipref;
        [rf,~] = TRAPS_flip(acqP.flipref,(acqP.nrep).*acqP.necho,kz,'opt',2,acqP.flipref,flipend,[kz kz kz acqP.necho]);
    end

    if(acqP.flipflag==3), flipend=60;
        kz=k0;
        if(kz<5),kz=5; end
        [rf,~] = TRAPS_flip(acqP.flipref,(acqP.nrep).*acqP.necho,k0,'90+a/2',2,acqP.flipref,flipend,[kz kz kz acqP.necho]);
        rf(1)=acqP.flipref;
    end

    if(acqP.flipflag==4), flipend=60;
        kz=k0;
        if(kz<5),kz=5; end
        [rf,~] = TRAPS_flip(acqP.flipref,(acqP.nrep)*acqP.necho,k0,'opt',2,acqP.flipref,flipend,[6 6 kz acqP.necho]);
    end

    if(acqP.flipflag==5),
        ka=k0;
        if(ka==1); ka=2; end
        flipend=60;
        [rf,~] = TRAPS_flip(acqP.flipref,(acqP.nrep)*acqP.necho,k0,'opt',2,180,60,[ceil(ka/2) ka ka acqP.necho]);

    end

    rf=rf(1:acqP.nrep*acqP.necho);
    if(strcmp(acqP.T2prep,'on'))
        rf=[rf];
    end
    acqP.flip=rf(1:acqP.nrep*acqP.necho);
    if(plotflag(2)=='1')
        figure
        plot(acqP.flip,'.-'), axis([xlim 0 180]);
        title('flip angle along echotrain')
    end
    %%
    % the following lines can be usedto calculate the total RF-power of the
    % sequence in a.u. For comparison of different sequences the calculated
    % values can be referenced to a reference sequence (e.g. a sequence using
    % 180

    [RFP,RFPtr,RFPrel] = calcRFP(rfex,rfref,acqP);
    acqP.pow=RFPrel;
    disp(strcat('rel. power compared to 180° refocussing :',num2str(acqP.pow)));
    rfpower=acqP.pow;
    %% Combine gradients
    GS5.delay=tRefwd;
    GS7.delay=acqP.TE-acq.tSpS;
    GSesp=mr.addGradients({GS4,GS5,GS7});

    if(strcmp(acqP.T2prep,'on'))
        GS7.delay=acqP.TE-acq.tSpS-tRefwd;
        GSespSE=mr.addGradients({GS5p1,GS7});
    end
    %%
    GR5.delay=tRefwd;
    GR6.delay=tRefwd+acq.tSpR;
    GR7.delay=tRefwd+acq.tSpR+readoutTime;
    GResp=mr.addGradients({GR5,GR6,GR7});
    %%
    GR5.delay=0;
    GR6.delay=acq.tSpR;
    GR7.delay=acq.tSpR+readoutTime;
    GRespSE=mr.addGradients({GR5,GR6,GR7});
    %%
    GPRpre.delay=0;
    GPacq.delay=mr.calcDuration(GPRpre)-acq.dG;
    GPResp=mr.addGradients({GPRpre,GPacq});
    GPRpre.delay=mr.calcDuration(GPResp)-acq.dG;
    GPResp=mr.addGradients({GPResp,GPRpre});
    GPResp.delay=tRefwd;
    %%
    GSResp=GPResp;
    GSResp.channel='z';
    GSResp=mr.addGradients({GSesp,GSResp});
    %%
    if(any(acqP.navmode(:)>1))
        GRrew.delay=0;
        GRresp=GRrew;
        GRrew.delay=acq.tSpR+readoutTime;
        GRresp=mr.addGradients({GRresp,GRrew});
        GRresp.delay=tRefwd;
    end
    %%
    adc.delay=GS4.shape_dur+acq.tSpR;
    adcSE.delay=GS4.shape_dur+acq.tSpR-tRefwd;
    kech0=1;

    %% Define sequence blocks
    % Next, the blocks are put together to form the sequence
    for kex=acqP.nDummy:nex                                                 %  comment for interleaved acquisition
       for s=1:acqP.NSlices
            rfex.freqOffset=GSex.amplitude*acqP.sliceGAP*acqP.sliceThickness*acqP.OSlices(s);
            rfref.freqOffset=GSref.amplitude*acqP.sliceGAP*acqP.sliceThickness*acqP.OSlices(s);
            rfrefSE.freqOffset=GSrefSE.amplitude*acqP.sliceGAP*acqP.sliceThickness*acqP.OSlices(s);
            rfex.phaseOffset=rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex); % align the phase for off-center slices
            if(strcmp(acq.cpmg_mod,'alt')), rfref_phase=rfref_phase0; end
            rfref.phaseOffset=rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref); % dito
            rfrefSE.phaseOffset=rfref_phase-2*pi*rfrefSE.freqOffset*mr.calcRfCenter(rfrefSE); % dito
            seq.addBlock(delayTR);
            if(acqP.TI>0),
                rfIR.freqOffset=GSIR.amplitude*acqP.sliceGAP*acqP.sliceThickness*IRSlices(s);
                rfIR.phaseOffset=-2*pi*rfref.freqOffset*mr.calcRfCenter(rfIR);
                seq.addBlock(rfIR,GSIRtot);
                seq.addBlock(delayTI);
            end
            if(strcmp(acqP.fatsat,'on')), seq.addBlock(rf_fs,gz_fs); end
            if(strcmp(acqP.T2prep,'on'))
                kech0=2; kech=1;
                if (kex>0)
                    phaseArea=acq.phaseAreas(kech,kex);
                end
                GPpreSE = mr.makeTrapezoid('y',sysphas,'Area',phaseArea);
                GPrewSE = mr.makeTrapezoid('y',sysphas,'Area',-phaseArea);
                GPrewon=GPrew;
                GPpreSE.delay=0;
                GPrewSE.delay=acqP.TE-mr.calcDuration(GPrewSE)-tRefwd;
                GRPespSE=mr.addGradients({GPpreSE,GPrewSE});
                seq.addBlock(GS1);
                seq.addBlock(GS2,rfex);
                seq.addBlock(GS3p1);
                seq.addBlock(GS3p2,GRpre_s);
                seq.addBlock(delayTE1,GDx,GDy,GDs);
                seq.addBlock(GS7SE,GRsprSE);
                seq.addBlock(GS4SE,rfrefSE);
                if(strcmp(acq.cpmg_mod,'alt')), rfref_phase=mod(rfref_phase+pi,2*pi); end
                seq.addBlock(GS5SE,GRsprSE);
                seq.addBlock(delayTE2,GDxr,GDyr,GDsr);
                seq.addBlock(GSespSE,GRespSE,GRPespSE,adcSE);

            else
                seq.addBlock(GS1);
                seq.addBlock(GS2,rfex);
                seq.addBlock(GS3,GR3);

            end
            phaseArea=0;
            for kech=kech0:acqP.necho*acqP.nrep
                if (kex>0)
                    phaseArea=acq.phaseAreas(kech,kex);
                end
                GPpre = mr.makeTrapezoid('y',sysphas,'Area',phaseArea);
                GPrew = mr.makeTrapezoid('y',sysphas,'Area',-phaseArea);
                GPrewon=GPrew;
                GPpre.delay=GS4.shape_dur;
                GPrew.delay=acqP.TE-mr.calcDuration(GPrew);
                GPesp=mr.addGradients({GPpre,GPrew});
                if(kech<=acq.nPinit)
                    rfref.signal=refenvelope0*acqP.flip(kech)/180;
                else
                    rfref.signal=refenvelope*acqP.flip(kech)/180;
                end

                if(strcmp(acq.cpmg_mod,'alt')), rfref.phaseOffset=rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref); end
                if((kech==1)&&acqP.navmode(1)==2)
                    seq.addBlock(GRresp,GPResp,GSesp,adc,rfref);

                elseif((kech==1)&&acqP.navmode(1)==3)
                    seq.addBlock(GRresp,GSResp,adc,rfref);

                elseif((kech==acqP.navk)&&acqP.navmode(1)==2)
                    seq.addBlock(GRresp,GPResp,GSesp,adc,rfref);
                elseif((kech==acqP.navk)&&acqP.navmode(1)==3)
                    seq.addBlock(GRresp,GSResp,adc,rfref);

                elseif((kech==acqP.necho)&&acqP.navmode(3)==2)
                    seq.addBlock(GRresp,GPResp,GSesp,adc,rfref);
                elseif((kech==acqP.necho)&&acqP.navmode(3)==3)
                    seq.addBlock(GRresp,GSResp,adc,rfref);


                else
                    seq.addBlock(GResp,GSesp,GPesp,adc,rfref);
                end
                if(strcmp(acq.cpmg_mod,'alt')), rfref_phase=mod(rfref_phase+pi,2*pi); end


            end
            GS4.delay=0;
            GS5.delay=0;
            if(strcmp(acqP.flipback,'on'))
                rfex.phaseOffset=pi+rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex);
                seq.addBlock(GS3rew);
                seq.addBlock(GS2,rfex);
                seq.addBlock(GS1rew);
            else
                seq.addBlock(GS4);
                seq.addBlock(GS5);
            end
            seq.addBlock(GSspr_end);



       end                                                                  %  comment for interleaved acquisition 

    end                                                                     %  comment for interleaved acquisition                                                                    
    acqP.necho=nPEecho; %reset number of echoes
    acqP.Nx=npoints;    %set Nx to number of samples
    %return                                                                 % Uncomment for interleaved acquisition 
    %% check whether the timing of the sequence is correct
    [ok, error_report]=seq.checkTiming;

    if (ok)
        fprintf('Timing check passed successfully\n');
    else
        fprintf('Timing check failed! Error listing follows:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end
    toc
    tic
    %% k-space trajectory calculation
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
    %%
    % plot k-spaces
    if(plotflag(3)=='1')
        figure; plot(t_ktraj,ktraj'); title('k-space components as functions of time'); end% plot the entire k-space trajectory
    if(plotflag(4)=='1')
        figure; plot(ktraj(1,:),ktraj(2,:),'b',...
            ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
        axis('equal'); % enforce aspect ratio for the correct trajectory display
        title('2D k-space (x-y)');
    end
    if(plotflag(5)=='1')
        figure; plot(ktraj(1,:),ktraj(3,:),'b',...
            ktraj_adc(1,:),ktraj_adc(3,:),'r.'); % a 2D plot
        axis('equal'); % enforce aspect ratio for the correct trajectory display
        title('2D k-space(x-z)');
    end

    %% Write to file

    % The sequence is written to file in compressed form according to the file
    % format specification using the |write| method.
    if(~exist('seqno')), seqno=1; end
    seqname=strcat('TSE_',num2str(seqno));
    seq.write(strcat(seqname,'.seq'));
    seqno=seqno+1;
    save(seqname,'system','pulseflag_ex','pulseflag_ref','acqP','acq');
    %%

    if(plotflag(6)=='1')
        fig=seq.plot('TimeRange',[TRfill TRfill+0.02],'timeDisp','ms','stacked',1);
        %fig=seq.plot('TimeRange',[0 0.1],'timeDisp','ms','stacked',1);
    end
end
% [~, pns_n, pns_c, tpns]=seq.calcPNS('MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
% [pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('MP_GradSys_P034_X60.asc'); % Cima.X PNS
% [~, pns_n, pns_c, tpns]=seq.calcPNS('MP_GPA_K2309_2250V_951A_AS82_prisma.asc'); % prisma
toc
