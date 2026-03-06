% for Re = [100, 400, 1000]
%     fprintf('Starting Re = %d...\n', Re);
%     LidDrivenCavity; 
% end

for lax_factor = 1e-3:1e-3:5e-3
    fprintf('Starting blending factor = %d...\n', lax_factor);
    LidDrivenCavity; 
end