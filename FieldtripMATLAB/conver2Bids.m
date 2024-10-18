sub = {'OA', 'OB'};
age = [nan nan];
sex = {'m', 'm'};

for subindx = 1:numel(sub)
    cfg = [];
    cfg.method = 'copy';
    cfg.datatype = 'eeg';
    cfg.bidsroot = 'bids';
    cfg.sub = sub{subindx};
    cfg.participants.age = age(subindx);
    cfg.participants.sex = sex{si}
end