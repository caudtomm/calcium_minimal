function RF_NT_NoiseTest1(traceResult)
% traceResult : rows are cells


%tracemat=traceResult.timeTraceMat;
tracemat = traceResult;

framerate=7.8125;

dtracemat=diff(tracemat(:,45:end)');

% mean(abs(dtracemat))
% mean(mean(abs(dtracemat)))


noiselevels=median(abs(dtracemat))/sqrt(framerate)
mean(noiselevels)

figure(789);clf;
subplot(221);
imagesc(tracemat);
%setfig3;
colorbar;
title(['Mean noise: ',num2str(mean(noiselevels))]);

subplot(222);
imagesc(dtracemat);
%setfig3;
colorbar;

subplot(223);
plot(noiselevels);
%setfig3;
xlabel('ROI #');
ylabel('Noise level')

subplot(224);
histogram(noiselevels)
xlabel('Noise level')


end