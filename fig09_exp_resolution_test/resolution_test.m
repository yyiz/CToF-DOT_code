
defaultPath = "../dat/all_res_tests";

if (testSep == 2)
    m_file = "res_test_6p5mm_skull_meas_06-Dec-2020_22-02-17";
    J_file = "res_test_6p5mm_skull_Jacobian_06-Dec-2020_22-01-46";
    dat_dir = "res_2mm_06-Dec-2020_res_test";

    mpath_full = sprintf("%s/%s/%s", defaultPath, dat_dir, m_file);
    Jpath_full = sprintf("%s/%s/%s", defaultPath, dat_dir, J_file);

    fistaOpts.lam1 = 2e7;
    fistaOpts.lam2 = 0;
    fistaOpts.maxItr = 100;
    fistaOpts.tol = 0;
    fistaOpts.nonneg = true; 
    fistaOpts.showFigs = false;
    
    load(Jpath_full);
    load(mpath_full);
    
    m = sum(m_diff, 2);
    J = sum(Jacobian, 3);

elseif (testSep == 1)
    m_file = "res_test_6p5mm_skull_meas_06-Dec-2020_23-03-11";
    J_file = "res_test_6p5mm_skull_Jacobian_06-Dec-2020_22-56-08";
    dat_dir = "res_1mm_06-Dec-2020_res_test";

    mpath_full = sprintf("%s/%s/%s", defaultPath, dat_dir, m_file);
    Jpath_full = sprintf("%s/%s/%s", defaultPath, dat_dir, J_file);

    fistaOpts.lam1 = 0;
    fistaOpts.lam2 = 0;
    fistaOpts.maxItr = 90;
    fistaOpts.tol = 0;
    fistaOpts.nonneg = true; 
    fistaOpts.showFigs = false;
    
    load(Jpath_full);
    load(mpath_full);
    
    m = sum(m_diff, 2);
    J = sum(Jacobian, 3);
    
elseif (testSep == 0.5)
    m_file = "res_test_6p5mm_skull_meas_10-Dec-2020_00-07-03";
    J_file = "res_test_6p5mm_skull_Jacobian_09-Dec-2020_23-59-58";
    dat_dir = "res_0p5mm_12_9_20_res_test";
    
    mpath_full = sprintf("%s/%s/%s", defaultPath, dat_dir, m_file);
    Jpath_full = sprintf("%s/%s/%s", defaultPath, dat_dir, J_file);
    
    load(Jpath_full);
    load(mpath_full);
    
    fistaOpts.lam1 = 5e9;
    fistaOpts.lam2 = 0;
    fistaOpts.maxItr = 100;
    fistaOpts.tol = 0;
    fistaOpts.nonneg = true; 
    fistaOpts.showFigs = false;
    
    adjustExpData;
else
    error("Invalid testSep: %f", testSep);
end

fwdType = "mvp";
sizeX = npts;

[f_fista, fT_fista, fistaStep] = mat2Handle(J, fwdType);

fprintf("Running FISTA...");
tic; [fistaRecon, fistaErr] = fista(m,f_fista,fT_fista,fistaStep,sizeX,fistaOpts); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

%% Display results

f1 = figure();
plot(fistaRecon, 'r', 'LineWidth', 2);

if (testSep ~= 0.5)
    meas_fp = load(mpath_full);
    J_fp = load(Jpath_full);
    line_start = meas_fp.xstart;
    line_start_ind = (line_start - J_fp.xstart) / J_fp.lineW + 1;
    line_w = meas_fp.lineW;
    if ~isfield(meas_fp, 'lineSep')
        line_sep = 1;
    else
        line_sep = meas_fp.lineSep;
    end
    rectWidth = 1.0;
    line1_start = line_start_ind;
    line2_start = line1_start + line_w + line_sep;

    xtick_vec = get(gca, 'XTick');
    xtick_vec = (xtick_vec - 1) * line_w + J_fp.xstart;
    xtick_str = sprintfc('%d',xtick_vec);
    xticklabels(xtick_str);
else
    J_lineW = 0.5;
    J_xstart = -(npts * J_lineW)/2;
    meas_fp = load(mpath_full);
    J_fp = load(Jpath_full);
    line_start = meas_fp.xstart;
    line_start_ind = (line_start - J_xstart) / J_lineW + 1;
    line_w = meas_fp.lineW;
    if ~isfield(meas_fp, 'lineSep')
        line_sep = 1;
    else
        line_sep = meas_fp.lineSep;
    end
    
    rectWidth = 1.0;
    line1_start = line_start_ind;
    line2_start = line1_start + line_w + line_sep;

    xtick_vec = get(gca, 'XTick');
    xtick_vec = (xtick_vec - 1) * line_w + J_xstart;
    xtick_str = sprintfc('%0.1f',xtick_vec);
    xticklabels(xtick_str);
end


ylabel("$\mu_a$ ($mm^{-1}$)", 'interpreter', 'latex');
xlabel("Lateral position (mm)", 'FontName', 'Times New Roman');

set(gca, 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman');

testSepStr = strrep(sprintf("%0.1d", testSep), '.', 'p');
saveFname = sprintf("%s/resTest_recon_%smm", savePath, testSepStr); 

%% Save results

if ~(exist(savePath, 'dir'))
    mkdir(savePath);
end

matname = sprintf("%s.mat", saveFname);
if (exist(matname, 'file'))
    fprintf("File exists. Overwrite?\n");
    keyboard;
end

save(matname, "-v7.3", "m_diff", "J", "fwdType", "sizeX",...
    "fistaOpts", "fistaStep", "fistaRecon", "fistaErr",...
    "fistaRuntime");
export_fig(f1, sprintf("%s.png", saveFname), '-m3', '-transparent', '-png');
savefig(f1, sprintf("%s.fig", saveFname));
fprintf("Saved file to %s\n", saveFname);


