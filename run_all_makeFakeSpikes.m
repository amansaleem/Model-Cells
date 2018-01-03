function run_all_makeFakeSpikes(is_common)

idx = 0;
is_common = 0;

if is_common
    flag_indep = 0;
else
    flag_indep = 1;
end

idx = idx + 1;
expt_list(idx).animal     = 'M130920_BALL';
expt_list(idx).iseries    = 1025;
expt_list(idx).expt_list  = 102:103;

idx = idx + 1;
expt_list(idx).animal     = 'M130918_BALL';
expt_list(idx).iseries    = 1030;
expt_list(idx).expt_list  = 103:105;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 530;
expt_list(idx).expt_list  = 104:106;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 531;
expt_list(idx).expt_list  = 103:106;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 601;
expt_list(idx).expt_list  = 103:106;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 602;
expt_list(idx).expt_list  = 102:106;

idx = idx + 1;
expt_list(idx).animal     = 'M140502_BALL';
expt_list(idx).iseries    = 603;
expt_list(idx).expt_list  = 107:110;

idx = idx + 1;
expt_list(idx).animal     = 'M140502_BALL';
expt_list(idx).iseries    = 604;
expt_list(idx).expt_list  = 107:110;

for idx = 1:length(expt_list)
    display(['Processing ' num2str(idx)]);
    es = VRLoadMultipleExpts( ...
        expt_list(idx).animal, expt_list(idx).iseries, ...
        expt_list(idx).expt_list);
    
    [es pop] = genFakePopSpikes(es, 35, 'P', is_common, flag_indep);
    
    es.animal = expt_list(idx).animal;
    es.iseries = expt_list(idx).iseries;
    es.expt_list = expt_list(idx).expt_list;
    es.pop = pop;
    
    if is_common
        filename = ['./Data/es_model_common_exp_' num2str(idx)];
    else
        filename = ['./Data/es_model_indep_exp_' num2str(idx)];
    end
    save(filename, 'es')
end