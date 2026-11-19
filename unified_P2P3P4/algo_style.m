function s = algo_style(id)
persistent TABLE
if isempty(TABLE)
    ids = {
        'SRDO'
        'SRDO_NOLR'
        'SRDO_HP'
        'LSHADE'
        'IPOP_CMAES'
        'PSO'
        'GA'
        'SAR'
        'DBO'
        'MSO'
        'SFOA'
        'GLO'
        'HBO'
        'BNCO'
    };
    disps = {
        'SRDO'
        'SRDO (w/o LR)'
        'SRDO-HP'
        'L-SHADE'           
        'IPOP-CMA-ES'      
        'PSO'
        'GA'
        'SAR'
        'DBO'
        'MSO'
        'SFOA'
        'GLO'
        'HBO'
        'BNCO'
    };
    colors = [
        0.84 0.15 0.16
        0.00 0.00 0.00
        0.50 0.50 0.50
        0.12 0.47 0.71   
        1.00 0.50 0.05  
        0.30 0.75 0.93  
        1.00 0.85 0.10   
        0.17 0.73 0.60   
        0.90 0.56 0.76  
        0.68 0.86 0.54  
        0.50 0.20 0.60   
        0.93 0.69 0.25  
        0.40 0.40 0.80   
        0.75 0.30 0.60   
    ];
    lss = {
        '--' '-' ':' '-.' '-' '--' ':' '-.' '-' '--' ':' '-.' '-' '--'
    };
    mks = {
        'o' 's' '^' 'd' 'v' '<' '>' 'p' 'h' '*' '+' 'x' 'o' 's'
    };

    TABLE = struct();
    for i = 1:numel(ids)
        TABLE.(ids{i}) = struct('disp', disps{i}, 'color', colors(i,:), ...
                                'ls', lss{i}, 'mk', mks{i});
    end
end

if nargin < 1 || isempty(id)
    s = TABLE;
    return;
end
if iscell(id)
    s = arrayfun(@(x) TABLE.(x), id, 'UniformOutput', false);
else
    if ~isfield(TABLE, id)
        error('algo_style: Unknown algorithm identifier ''%s''', id);
    end
    s = TABLE.(id);
end
end
