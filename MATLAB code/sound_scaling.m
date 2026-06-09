cd '/Users/ody-ozono/Library/CloudStorage/OneDrive-Personal/Desktop/School + Work/ibsc/Research Project/MRI SOUNDS FOR PROJECT'


suffix{1} = '1.5T ABDO T1 FAST GRE';
suffix{2} = '1.5T ABDO T2 SSFSE';
suffix{3} = '3T PELVIC T1 LAVA-FLEX';
suffix{4} = '3T PELVIC T2 DWI-EPI';

p_ref = 2e-5;

for i=1:length(suffix)

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
    SPL = 20*log10(abs(y_scaled)/p_ref);

    audiowrite([suffix{i},' scaled.wav'], y_scaled, Fs, 'BitsPerSample', 32)
    
        
end