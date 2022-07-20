using DifferentialEquations
using Pertussis
using MAT
using Plots
demopath = "DataSource/demo_para_ontario.mat"
str = "pop_groups"
pop_groups = matread(demopath)[str]

casepath = "DataSource/ontario_incidence.mat"
casestr = "Ontario_monthly_incidence"
Ontario_monthly_incidence = matread(casepath)[casestr]
function solve_pertussis(theta)
    Tspan = (2011, 2020)
    n_agegroup = 20
    ED_0 = [Ontario_monthly_incidence[1, 1] * pop_groups[1] / sum(pop_groups[1:3]); Ontario_monthly_incidence[1, 1] * pop_groups[2] / sum(pop_groups[1:3]); Ontario_monthly_incidence[1, 1] * pop_groups[3] / sum(pop_groups[1:3]); Ontario_monthly_incidence[1, 2]; Ontario_monthly_incidence[1, 3]; Ontario_monthly_incidence[1, 4]; Ontario_monthly_incidence[1, 5] * pop_groups[7] / sum(pop_groups[7:12]); Ontario_monthly_incidence[1, 5] * pop_groups[8] / sum(pop_groups[7:12]); Ontario_monthly_incidence[1, 5] * pop_groups[9] / sum(pop_groups[7:12]); Ontario_monthly_incidence[1, 5] * pop_groups[10] / sum(pop_groups[7:12]); Ontario_monthly_incidence[1, 5] * pop_groups[11] / sum(pop_groups[7:12]); Ontario_monthly_incidence[1, 5] * pop_groups[12] / sum(pop_groups[7:12]); Ontario_monthly_incidence[1, 6] * pop_groups[13] / sum(pop_groups[13:15]); Ontario_monthly_incidence[1, 6] * pop_groups[14] / sum(pop_groups[13:15]); Ontario_monthly_incidence[1, 6] * pop_groups[15] / sum(pop_groups[13:15]); Ontario_monthly_incidence[1, 7] * pop_groups[16] / sum(pop_groups[16:20]); Ontario_monthly_incidence[1, 7] * pop_groups[17] / sum(pop_groups[16:20]); Ontario_monthly_incidence[1, 7] * pop_groups[18] / sum(pop_groups[16:20]); Ontario_monthly_incidence[1, 7] * pop_groups[19] / sum(pop_groups[16:20]); Ontario_monthly_incidence[1, 7] * pop_groups[20] / sum(pop_groups[16:20])] ./ eta_EDtoD
    D_0 = 0 * ones(n_agegroup)
    H_0 = 0 * ones(n_agegroup)
    Sr_0 = 0 * ones(n_agegroup)
    EU_0 = 1 * ones(n_agegroup)
    U_0 = 0 * ones(n_agegroup)
    RC_0 = ED_0 .* eta_EDtoD .* 2
    RI_0 = 0 * ones(n_agegroup)
    RD_0 = 0 * ones(n_agegroup)
    VC_0 = [0; 2.36 / 12 * 0.857 * pop_groups[8, 2]; 1.49 / 12 * 0.928 * pop_groups[8, 3]; 0; 0; 0.677 * 0.925 * pop_groups[8, 5]; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
    VD_0 = 0 * ones(n_agegroup)
    S_0 = pop_groups[8, :] - ED_0 - D_0 - H_0 - Sr_0 - EU_0 - U_0 - RC_0 - RI_0 - RD_0 - VC_0 - VD_0
    StatVar_0 = [S_0; ED_0; D_0; H_0; Sr_0; EU_0; U_0; RC_0; RI_0; RD_0; VC_0; VD_0]
    prob = ODEProblem(pertussisage, StatVar_0, Tspan, theta)
    solve(prob)
end

theta = [0.9; 0.5; 0.5; 0.5; 0.1; 0.01; 0.5; 0.15; 0.03; 0.03; 0.003; 0.003; 1;0];

sol=solve_pertussis(theta)

plot(sol.t,sol[41,:])

