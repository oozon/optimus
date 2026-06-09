load sound_transmission.mat

p_ref = 1; % reference pressure (Pa)
barycentre_pressures = [1, barycentre_pressures]; % adding a pressure with mag 1 
frequencies = [0, frequencies]; % adding frequency 0 

barycentre_SPLs = 20*log10(abs(barycentre_pressures)./p_ref);
barycentre_unwrapped_phase = unwrap(angle(barycentre_pressures));
barycentre_magnitude = abs(barycentre_pressures);

N_points = 16385; %based on the Nature paper 
Fs = 44100; %sampling frequency 
max_freq = Fs/2;
 

fq = linspace(0, max_freq, N_points);

interp_phase = interp1(frequencies, barycentre_unwrapped_phase, fq, 'spline');
interp_mag = interp1(frequencies, barycentre_magnitude, fq, 'spline');

delay = 0;

%linearise phase using constant delay 
linear_phase = -2 * pi * delay * fq/Fs;
phase_final = linear_phase + interp_phase;

%generate complex number with linearised, interpolated phase and
%interpolated magnitude 
complex_resp = interp_mag .* exp(1i * phase_final);

filter_order = 1000;
% ensure frequencies cover the band and H has same length
H = complex_resp; % complex target
w = 2*pi * fq / Fs;   % rad/sample in [0, pi]

na = 0;

% invfreqz fit (returns numerator b and denominator a)
[b_fit, a_fit] = invfreqz(H, w, filter_order, na);

% b_fit is the FIR coefficients; a_fit should be 1 for na=0

[h,f] = freqz(b_fit, a_fit, fq, Fs);

% figure;
% subplot(2,1,1);
% plot(f, 20*log10(abs(h)), 'b', LineWidth=2);
% hold on;
% plot(frequencies, barycentre_SPLs, 'r--', LineWidth=2); % raw samples
% ylabel('Magnitude (dB)');
% legend('FIR SPLs', 'Raw SPLs');
% grid on;
% title('Overlayed FIR and measured magnitude/frequency plot')
% hold off;
% 
% subplot(2,1,2);
% plot(f, unwrap(angle(h)), 'b', LineWidth=2);
% hold on;
% plot (frequencies, barycentre_unwrapped_phase, 'r--', LineWidth=2);
% xlabel('Frequency (Hz)');
% ylabel('Phase (rad)');
% legend('FIR phase', 'Raw unwrapped phase');
% grid on;
% title('Overlayed FIR and measured unwrapped phase/frequency plot')


%convolution of MRI audio

cd '/Users/ody-ozono/Library/CloudStorage/OneDrive-Personal/Desktop/School + Work/ibsc/Research Project/MRI SOUNDS FOR PROJECT'


suffix{1} = '1.5T ABDO T1 FAST GRE';
suffix{2} = '1.5T ABDO T2 SSFSE';
suffix{3} = '3T PELVIC T1 LAVA-FLEX';
suffix{4} = '3T PELVIC T2 DWI-EPI';

p_ref = 2e-5;

for i=3:3

    [y, Fs] = audioread([suffix{i},'.wav']);   

    if i==1
        l_peak = 106;
    elseif i==2
        l_peak = 116;
    elseif i==3
        l_peak = 126;
    elseif i==4
        l_peak = 131;
    end


    v_max = max(abs(y(:,1)));

    a = (p_ref*10^(l_peak/20))/v_max;
    y_scaled = a * y(:,1);
    SPL_before = 20*log10(abs(y_scaled)/p_ref);
% 
%     % audiowrite([suffix{i},' scaled.wav'], y_scaled, Fs, 'BitsPerSample', 32) %make scaled files
% 
    y_filtered = filter(b_fit, a_fit, y_scaled); %put the raw scaled y data in the filter 
    SPL_after = 20*log10(abs(y_filtered)/p_ref);
% 
%     % audiowrite([suffix{i},' scaled+filtered.wav'], y_filtered, Fs, 'BitsPerSample', 32); %write the filtered
% 
%     % figure;
%     % 
%     % spectrogram(y_scaled(:,1), hamming(256), 128, 256, Fs, 'yaxis');
%     % title(['Spectrogram: ', suffix{i}, ': ex-utero']);
%     % 
%     % figure;
%     % 
%     % spectrogram(y_filtered(:,1), hamming(256), 128, 256, Fs, 'yaxis');
%     % title(['Spectrogram: ', suffix{i}, ': in-utero']);
%     % 
%     % 
%     t = (0:length(y)-1) / Fs;
%     figure; 
%     plot(t, SPL_before)
%     title(['Time history: ', suffix{i}, ': ex-utero'])
%     xlabel('Time (s)');
%     ylabel('Amplitude (dB)');
%     fontsize(gcf, 23, 'points')
% 
%     figure; 
%     plot(t, SPL_after)
%     title(['Time history: ', suffix{i}, ': in-utero'])
%     xlabel('Time (s)');
%     ylabel('Amplitude (dB)');
%     fontsize(gcf, 20, 'points')
%     % 
%     % 
%     % % Plot Power Spectral Density (PSD) for before and after signals
%     % nfft = 2^14; % FFT points for good resolution
%     % window = hamming(1024);
%     % noverlap = 512;
%     % 
%     % % PSD of scaled (before)
%     % [pxx_before, f_psd] = pwelch(y_scaled, window, noverlap, nfft, Fs);
%     % pxx_before_dB = 10*log10(pxx_before);
%     % 
%     % % PSD of filtered (after)
%     % [pxx_after, ~] = pwelch(y_filtered, window, noverlap, nfft, Fs);
%     % pxx_after_dB = 10*log10(pxx_after);
%     % 
%     % figure;
%     % plot(f_psd, pxx_before_dB, 'b', 'LineWidth', 1.5);
%     % hold on;
%     % plot(f_psd, pxx_after_dB, 'r', 'LineWidth', 1.5);
%     % hold off;
%     % xlim([0, Fs/2]);
%     % xlabel('Frequency (Hz)');
%     % ylabel('Power/Frequency (dB/Hz)');
%     % title(['Power Spectral Density: ', suffix{i}]);
%     % legend('Ex-utero (scaled)', 'In-utero (scaled+filtered)', 'Location', 'best');
%     % grid on;
%     % 
%     peak_before = max(SPL_before)
%     peak_after = max(SPL_after)
%     % 
%     % 
%     % suffix{i}
%     % % 
%     % % y_scaled_rms = rms(y_scaled);
%     % % y_filtered_rms = rms(y_filtered);
%     % % spl_rms_before = 20 * log10 (y_scaled_rms/p_ref) %RMS BEFORE 
%     % % spl_rms_after = 20 * log10 (y_filtered_rms/p_ref) %RMS AFTER
%     % % 
    [p, cf] = poctave(y_scaled, Fs,'BandsPerOctave', 3);
    [p1, cf1] = poctave(y_filtered, Fs,'BandsPerOctave', 3);
% 
    % add labels and title to the octave-band plots and plot both on same axes
    figure;
    semilogx(cf, p, '-b', 'LineWidth', 1.5);
    hold on;
    semilogx(cf, p1, '-r', 'LineWidth', 1.5);
    hold off;
    xlabel('Centre frequency (Hz)');
    ylabel('Level (dB)');
    title(['Octave spectrum: ', suffix{i}, ' (Ex-utero vs In-utero)']);
    grid on;
    legend({'Ex-utero (scaled)','In-utero (scaled+filtered)'}, 'Location', 'best');
    fontsize(gcf, 23, 'points')
% 

%     % % Define the 20-second window
%     % start_time = 0; % in seconds
%     % end_time = 20;  % in seconds
%     % indices = (start_time*Fs + 1) : (end_time*Fs);
%     % y_20s_scaled = y_scaled(indices, :); % Extracts 20 seconds
%     % y_20s_filtered = y_filtered(indices, :);
%     % 
%     % % Calculate Leq
%     % % Assuming y is calibrated to Pascal
%     % rms20s_scaled = rms(y_20s_scaled); 
%     % rms20s_filtered = rms(y_20s_filtered);
%     % leq_20s_scaled = 20 * log10(rms20s_scaled / p_ref);
%     % leq_20s_filtered = 20 * log10(rms20s_filtered / p_ref);
%     % 
%     % disp(['Leq for 20s (scaled): ', num2str(leq_20s_scaled), ' dB']);
%     % disp(['Leq for 20s (filtered): ', num2str(leq_20s_filtered), ' dB']);
% 
% 
%     % SPL_calc = splMeter('FrequencyWeighting','Z-weighting','SampleRate',48000);
%     % [~,Leq,~,~] = SPL_calc(y_scaled);
%     % 
%     % figure;
%     % plot(Leq);
%     % xlabel('Time (samples)');
%     % ylabel('Leq (dB)');
%     % title(['Leq of ', suffix{i}]);
%     % grid on;
% 
% 
end
