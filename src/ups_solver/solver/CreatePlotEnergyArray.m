%% Create the arrys to record the energies at each iteration.
% input: options.maxit       the maximum iteration.
% output: some arrays        used to record energies.
function plot_energy = CreatePlotEnergyArray(options)
plot_energy.tab_primal      = zeros(1,options.maxit);
plot_energy.tab_dual        = zeros(1,options.maxit);
plot_energy.tab_primal_abs  = zeros(1,options.maxit);
plot_energy.tab_dual_abs    = zeros(1,options.maxit);
plot_energy.tab_theta_min   = zeros(1,options.maxit);
plot_energy.tab_res_theta   = zeros(1,options.maxit);
plot_energy.tab_res_z       = zeros(1,options.maxit*options.LinS.maxit);
plot_energy.tab_energy      = zeros(1,4*options.maxit+1);
plot_energy.tab_objective   = zeros(1,4*options.maxit+1);
plot_energy.tab_no_smooth   = zeros(1,4*options.maxit+1);
plot_energy.tab_res_rho_gap = zeros(1,options.maxit);
plot_energy.tab_res_s       = zeros(1,options.maxit);
plot_energy.tab_z           = zeros(1,options.maxit);
plot_energy.tab_theta       = zeros(1,options.maxit);
plot_energy.tab_u           = zeros(1,options.maxit);
plot_energy.tab_rho         = zeros(1,options.maxit);
plot_energy.tab_s           = zeros(1,options.maxit);
end