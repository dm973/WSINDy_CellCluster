function [f_learned,h_learned,d_learned] = gen_force_fcn_xv(W,fx_fcn_cell,fv_fcn_cell,hx_fcn_cell,hv_fcn_cell,dx_fcn_cell,dv_fcn_cell)
    
    J_fx = length(fx_fcn_cell);
    J_fv = length(fv_fcn_cell);
    J_hx = length(hx_fcn_cell);
    J_hv = length(hv_fcn_cell);
    J_dx = length(dx_fcn_cell);
    J_dv = length(dv_fcn_cell);
    
    Wf = W(1:J_fx*J_fv);
    Wh = W(J_fx*J_fv+1:J_fx*J_fv+J_hx*J_hv);
    Wd = W(J_fx*J_fv+J_hx*J_hv+1:J_fx*J_fv+J_hx*J_hv+J_dx*J_dv);
    
    if any(Wf)
        f_learned = build_F(Wf,fx_fcn_cell,fv_fcn_cell);
    else
        f_learned = [];
    end
    if any(Wh)
        h_learned = build_F(Wh,hx_fcn_cell,hv_fcn_cell);
    else
        h_learned = [];
    end
    if any(Wd)
        d_learned = build_F(Wd,dx_fcn_cell,dv_fcn_cell);
    else
        d_learned = [];
    end

end

function F = build_F(w,fx_fcn_cell,fv_fcn_cell)
    if ~all(w==0)
        if ~exist('fv_fcn_cell','var')
            fv_fcn_cell = {@(v) v*0+1};
        end        
        if isempty(fv_fcn_cell)
            fv_fcn_cell = {@(v) v*0+1};
        end
        F= @(x,v) x*0+v*0;
        for j=1:length(fx_fcn_cell)
            for i=1:length(fv_fcn_cell)
                if w((j-1)*length(fv_fcn_cell)+i)~=0
                    F = @(x,v) F(x,v) + w((j-1)*length(fv_fcn_cell)+i)*fx_fcn_cell{j}(x).*fv_fcn_cell{i}(v);
                end
            end
        end
    else
        F = [];
    end
end

