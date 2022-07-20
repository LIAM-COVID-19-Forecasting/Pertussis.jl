using Pertussis
using DifferentialEquations

u0 = 0.1
tspan = (0.0, 10.0)
p = 1.0
pertussisage1(u, p, t) = sin(2π * (t)) + p * cos(2π * (t))
prob = ODEProblem(pertussisage1, u0, tspan, p)
sol = solve(prob)
plot(sol)
function pertussisage2(du, u, p, t)
    S = u[1]
    I = u[2]
    N = sum(u)
    β, γ = p
    du[1] = -β * S * I / N + γ * I
    du[2] = β * S * I / N - γ * I
end
u01 = [0.9, 0.1]
p = [0.3, 0.1]
tspan = (0.0, 20.0)
prob2 = ODEProblem(pertussisage2, u01, tspan, p)
sol = solve(prob2)

using Plots
plot(sol)


contactmatpath = "DataSource/C20_Y_2004_to_2040.mat"

contactmatstr = "C20_Y_2004_to_2040"

contactfunc = contactmatrix(contactmatpath, contactmatstr)

size(contactfunc)

x = 2021:0.1:2050

plot(x, contactfunc[1, 4].(x))

AA = [contactfunc[i, j](1) for i in 1:20, j in 1:20]

heatmap(AA)

using MAT

heatmap(matread(contactmatpath)[contactmatstr][2])

AA == matread(contactmatpath)[contactmatstr][2]

ii = 3
matread(contactmatpath)[contactmatstr][ii+1] == [contactfunc[i, j](ii) for i in 1:20, j in 1:20]


# demo
# Years_demo=Years_deaths_ontario;
# r=Aging_Rates;
# B=births_table_ontario;
# mu=Death_Rates_calc;
# nu=PC_Migration_Rates;

demopath = "DataSource/demo_para_ontario.mat"
demostr = "births_table_ontario"
birthyear = demorate(demopath, demostr)

demopath1 = "DataSource/demo_para_ontario_monthly.mat"
demostr1 = "births_table_ontario_even"
birthyear1 = demorate(demopath1, demostr1)
x = 0:10
plot(x, birthyear1.(x))



agegroup_str = ["0-1m" "2-11m" "1-4y" "5-9y" "10-14y" "15-19y" "20-24y"
    "25-29y" "30-34y" "35-39y" "40-44y" "45-49y" "50-54y" "55-59y" "60-64y"
    "65-69y" "70-74y" "75-79y" "80-84y" "85+y"];
using LaTeXStrings
paramstr = [L"\sigma_1" L"\sigma_{2p}" L"\sigma_{3p}" L"\sigma_{4p}" L"\sigma^r_scale"
    L"\rho" L"p_1" L"p_{2p}" L"p_{3p}" L"p_{4p}" L"t_2" L"t_3" L"\epsilon_p"]


rate_DTaP_2 = [5.53; 4.21; 4.21; 3.86; 3.59; 2.89; 4.21; 2.36; 1.73; 1.73;2.36;4.21; 4.69; 2.89; 2.89; 4.69; 4.69]
rate_DTaP_3 = [1.59; 1.59; 1.59; 1.59; 1.59; 1.43; 1.49; 1.43; 1.29; 1.88; 1.41;1.84; 1.87; 1.95; 1.95; 2.66; 2.66]
psi = psifun(rate_DTaP_2, rate_DTaP_3)