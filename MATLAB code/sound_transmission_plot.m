load '/Users/ody-ozono/Library/CloudStorage/OneDrive-Personal/Desktop/School + Work/ibsc/Research Project/Sound Stuff matlab/sound_transmission.mat'

% Calculate the sound pressure level (SPL) in dB
p_ref = 1; % reference pressure (Pa)
barycentre_SPLs = 20*log10(abs(barycentre_pressures)./p_ref);
barycentre_unwrapped_phase = unwrap(angle(barycentre_pressures));
barycentre_magnitude = abs(barycentre_pressures);

p_rms = rms(barycentre_magnitude);
L_rms_dB = 20*log10(p_rms ./ p_ref);


figure
semilogx (frequencies, barycentre_SPLs, '-b') % use if you want to plot sound pressure levels (dB) instead of acoustic pressure (Pa)

hold on

semilogx(frequencies, repmat(L_rms_dB, size(frequencies)), '--r')

hold off

xlabel('Frequency (Hz)');
ylabel('In-utero sound pressure level (dB re 1 Pa)');
title('SPL Uterus Barycentre');
grid on;
legend('SPLs', 'L_{rms}');
