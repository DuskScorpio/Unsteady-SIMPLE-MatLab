% for Re = [100, 400, 1000]
%     fprintf('Starting Re = %d...\n', Re);
%     LidDrivenCavity; 
% end
for Re = [100, 400, 1000]
    for lax_factor = 5e-2:5e-2:0.5
        fprintf('Starting blending factor = %d...\n', lax_factor);
        LidDrivenCavity; 
    end
end