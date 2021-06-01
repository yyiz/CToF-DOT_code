% This is a basic function to parse varargin for code. 
%
% Calling Convention:
% optargs = parseVarargin(var_args, name_list, opt_args)
%
% Contact: ar89@rice.edu
%
function opt_args = parseVarargin(var_args, name_list, opt_args)
    num_args = length(var_args);

    % varargin should always be even since all inputs are name-value pairs
    if mod(num_args, 2) == 1
        error('All inputs should be name value pairs.');
    end
    
    % iterate through name-value pairs
    for ind = 1:2:num_args
        
        % checks to see if valid names
        if ~ischar(var_args{ind})                                    % if name isn't a string, produce error
            error('Name designation must be type char.')
        elseif ~contains(name_list,var_args{ind},'IgnoreCase',true)  % if we can't find it, then produce error
            error('%s is not a known input.',var_args{ind});
        else                                                         % else, put in arguments if type is same
            log_ind = strcmpi(name_list, var_args{ind});
            if isa(var_args{ind+1}, class(opt_args{log_ind}))
                opt_args{log_ind} = var_args{ind+1};
            else
                error('%s must be type %s.', var_args{ind}, class(opt_args{log_ind}));
            end
        end
    end
end