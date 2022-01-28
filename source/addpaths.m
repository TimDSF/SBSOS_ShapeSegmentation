local = true;

if local
    setenv('MOSEKLM_LICENSE_FILE', '../licenses/mosek.lic');

    addpath(genpath('../libraries/HausdorffDist'));
    addpath(genpath('../libraries/mosek/9.3/toolbox'));
    addpath(genpath('../libraries/sedumi'));
    addpath(genpath('../libraries/SieveSDP-master'));
    addpath(genpath('../libraries/SOSTOOLS-3.300'));
    addpath(genpath('../libraries/Sparse_BSOS-master'));

else
    setenv('MOSEKLM_LICENSE_FILE', '../licences/mosek.lic');

    addpath(genpath('../libraries/HausdorffDist'));
    addpath(genpath('../libraries/mosek/9.3/toolbox'));
    addpath(genpath('../libraries/sedumi'));
    addpath(genpath('../libraries/SieveSDP-master'));
    addpath(genpath('../libraries/SOSTOOLS-3.300'));
    addpath(genpath('../libraries/Sparse_BSOS-master'));
end

clear local;