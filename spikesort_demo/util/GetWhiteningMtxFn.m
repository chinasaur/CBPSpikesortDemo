function WhiteningMtxFn = GetWhiteningMtxFn(pars)
    if (~isfield(pars, 'whitening_mtx_fn'))
        WhiteningMtxFn = @(len) eye(len);
    else
        WhiteningMtxFn = pars.whitening_mtx_fn;
    end
end
