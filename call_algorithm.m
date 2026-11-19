function [best_fit, best_pos, curve] = call_algorithm(name, N, Max_iter, lb, ub, dim, fobj)
    switch name
        case 'PSO'
            [best_fit, best_pos, curve] = PSO(N, Max_iter, lb, ub, dim, fobj);
        case 'GA'
            [best_fit, best_pos, curve] = GA(N, Max_iter, lb, ub, dim, fobj);
        case 'SAR'
            [best_fit, best_pos, curve] = SAR(N, Max_iter, lb, ub, dim, fobj);
        case 'DBO'
            [best_fit, best_pos, curve] = DBO(N, Max_iter, lb(1), ub(1), dim, fobj);
        case 'MSO'
            [best_fit, best_pos, curve] = MSO(N, Max_iter, lb, ub, dim, fobj);
        case 'SFOA'
            [best_pos, best_fit, curve] = SFOA(N, Max_iter, lb, ub, dim, fobj);
        case 'GLO'
            [best_fit, best_pos, curve] = GLO(N, Max_iter, lb, ub, dim, fobj);
        case 'HBO'
            [best_fit, best_pos, curve] = HBO(N, Max_iter, lb, ub, dim, fobj);
        case 'BNCO'
            [best_fit, best_pos, curve] = BNCO(N, Max_iter, lb, ub, dim, fobj);
        case 'SRDO'
            [best_fit, best_pos, curve] = SRDO(N, Max_iter, lb, ub, dim, fobj);
        otherwise
            error('Unknown algorithm: %s', name);
    end
    curve = curve(:)'; 
end