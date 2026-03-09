clear count* acqP acq
seq=mr.Sequence(system);
nseq=2;nex=16;
for kex=1:nex
    for kseq=1:nseq
        TSEplus_ileave;
    end
end

%% check whether the timing of the sequence is correct
disp('final calculations')
tic
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end


%% k-space trajectory calculation

[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%%
% plot k-spaces
if(plotflag(3)=='1')
    figure; plot(t_ktraj,ktraj'); axis([TRfill TRfill+acqP.TR/acqP.NSlices ylim]); title('k-space components as functions of time'); end% plot the entire k-space trajectory
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
end

toc
