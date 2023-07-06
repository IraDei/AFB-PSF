function [varginMap, varginList, vararginKeySet] = parse_varargin(segfix, varargin)
%PARSE_VARARGIN 
%   This automatic parsing function discriminates new param set via prefix
%   'vprm_' in string member varargin params.
    % segmentation fix init
    sf_len = size(segfix,2);
    abs_filter = 0;
    
    varginSet = {}; varginList = {};
    keySet = {};    % key set maps param title to param items 
    varginSet_size = [];
    n_var_collect = 0;
    n_vgin = nargin - 1;
    for i = 1:n_vgin
        accpet = 1;
        prm = varargin{i};
        if(ischar(prm))
            if abs_filter
                break;
            end
            
            lb_len = size(prm,2);
            if(lb_len>=sf_len && (strcmp(prm(1,1:sf_len),segfix) || strcmp(segfix, '*_')))
                % Specify prefix '*_' to parse all varargin sections
                if(n_var_collect>0)
                    varginList = {varginList, [segfix keySet{n_var_collect}], varginSet{n_var_collect}{:}};
                end
                n_var_collect = n_var_collect + 1;
                varginSet{n_var_collect} = {};
                
                varginSet_size(n_var_collect) = 0;
                % set title of current param block and init prm-block size
                if(sf_len<lb_len)
                    keySet{n_var_collect} = prm(sf_len+1:lb_len);
                else
                    keySet{n_var_collect} = prm(1:lb_len);
                    abs_filter = 1;
                end
                accpet = 0;
            end
        end
        
        if accpet && n_var_collect>0
            % numeric values are accepted as member of current param block
            grp_prm_ind = varginSet_size(n_var_collect) + 1;
            varginSet{n_var_collect}{grp_prm_ind} = prm;
            varginSet_size(n_var_collect) = grp_prm_ind;
        end
    end
    
    if(n_var_collect>0)
        varginList = [varginList, [segfix keySet{n_var_collect}], varginSet{n_var_collect}{:}];
    end
    
    % construct mapping from keySet to valSet
    if(~isempty(keySet) && ~isempty(varginSet))
        varginMap = containers.Map(keySet, varginSet, 'UniformValues', false);
    else
        varginMap = containers.Map; % return a empty map
    end
    vararginKeySet = keySet;
end

