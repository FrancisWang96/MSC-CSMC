function GOsim = importGO(DataSetName)

switch (DataSetName)
    case 'a'
        Xdata=csvread('./Datasets/a/simdata/a_sim_AIC.csv');
    case'cell'
        Xdata=csvread('./Datasets/cell/simdata/cell_sim_AIC.csv');
    case 'sporulation'
        Xdata=csvread('./Datasets/sporulation/simdata/sporulation_sim_AIC.csv');

    case'serum'
        Xdata=csvread('./Datasets/serum/simdata/serum_sim_AIC.csv');

    case'gal'
        Xdata=csvread('./Datasets/gal/simdata/Gal_sim_AIC.csv');
end
GOsim = Xdata;
end

















