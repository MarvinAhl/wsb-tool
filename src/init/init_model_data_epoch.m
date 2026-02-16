function model_data = init_model_data_epoch(et, is_4bp)
    model_data = init_model_data(is_4bp);

    [phi0, th0] = epoch_to_BCRFBP_angles(et);

    model_data.phi0 = phi0;
    model_data.th0 = th0;
end