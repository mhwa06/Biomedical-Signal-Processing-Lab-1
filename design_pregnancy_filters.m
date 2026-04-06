function Hd = design_pregnancy_filters(fs)
%DESIGN_PREGNANCY_FILTERS Create the filter chain for the noisy abdominal lead.

    cfg = assignment_config();

    [bHP, aHP] = butter(4, cfg.preg.hpCutoff/(fs/2), 'high');
    [bBS, aBS] = butter(2, cfg.preg.notch/(fs/2), 'stop');
    [bLP, aLP] = butter(4, cfg.preg.lpCutoff/(fs/2), 'low');

    Hd = struct();
    Hd.bHP = bHP; Hd.aHP = aHP;
    Hd.bBS = bBS; Hd.aBS = aBS;
    Hd.bLP = bLP; Hd.aLP = aLP;
end
