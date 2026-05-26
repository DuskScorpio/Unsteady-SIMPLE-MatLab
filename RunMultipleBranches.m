function RunMultipleBranches(branches, varargin)
    % RunMultipleBranches - Run LidDrivenCavity script on multiple branches sequentially
    % Usage: RunMultipleBranches({'Lax', 'Leap-Frog'})
    %        RunMultipleBranches({'Lax', 'Leap-Frog', 'master'})
    
    if nargin == 0
        branches = {'master', 'Leap-Frog'};
    end
    
    if ischar(branches)
        branches = {branches};
    end
    
    originalBranch = strtrim(string(system('git rev-parse --abbrev-ref HEAD')));
    
    fprintf('\n=== Running simulation on %d branches ===\n', length(branches));
    
    for i = 1:length(branches)
        branch = branches{i};
        fprintf('\n>>> Switching to branch: %s\n', branch);
        system(sprintf('git checkout %s', branch));
        
        fprintf('>>> Running LidDrivenCavity on %s\n', branch);
        try
            LidDrivenCavity();
            fprintf('✓ %s completed successfully\n', branch);
        catch ME
            fprintf('✗ Error on %s: %s\n', branch, ME.message);
        end
    end
    
    % Return to original branch
    fprintf('\n>>> Returning to original branch: %s\n', originalBranch);
    system(sprintf('git checkout %s', originalBranch));
    fprintf('\n=== All simulations completed ===\n');
end
