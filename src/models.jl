@doc raw"""
    pertussisage(rhs, StatVar, theta, t)

main model of pertussis.
agegroup_str = ["0-1m" "2-11m" "1-4y" "5-9y" "10-14y" "15-19y" "20-24y"
    25-29y" "30-34y" "35-39y" "40-44y" "45-49y" "50-54y" "55-59y" "60-64y"
    "65-69y" "70-74y" "75-79y" "80-84y" "85+y"]

paramstr: $[\sigma_1; \sigma_{2p}; \sigma_{3p};\sigma_{4p}; \sigma^r_scale;\rho; p_1; p_{2p}; p_{3p}; p_{4p}; t_{2}; t_3; \epsilon_p;\epsilon_t]$
"""
function pertussisage(rhs, StatVar, theta, t)
    n_agegroup = 20#number of age groups
    Years = collect(2011:1:2020)#years considered in this study

    #disease related parameter values
    eta_EDtoD = 365 / 9 * ones(n_agegroup)
    eta_EUtoU = 365 / 8 * ones(n_agegroup)
    gamma_DtoH = 365 / 15.6 * ones(n_agegroup)
    gamma_HtoRC = 365 / 5 * ones(n_agegroup)
    gamma_DtoRC = 365 / 21 * ones(n_agegroup)
    gamma_UtoRC = 365 / 15 * ones(n_agegroup)

    #waning of natural immunity
    tau_RCtoRI = 1 / 10 * ones(n_agegroup)
    tau_RItoRD = 1 / 15 * ones(n_agegroup)
    tau_RDtoSr = 1 / 5 * ones(n_agegroup)

    #vaccine-related parameter values
    #waning of immunity
    #Years_vacc=collect(2004:1:2020)
    tau_VCtoVD = ones(n_agegroup)
    tau_VCtoVD[1] = 1 / 1.48
    tau_VCtoVD[3:5] .= 1 / 4.955
    tau_VCtoVD[6:end] .= 1 / 4.02
    tau_VDtoSr = copy(tau_VCtoVD)

    #vaccine efficacy
    alpham = [0.914; 0; 0.086]# Materanl booster efficacy

    alpha = Array{Vector{Float64}}(undef, 4, 4)# each element=row (year)*column (age group) 
    for j = 1:4
        for i = 1:4
            alpha[i, j] = zeros(n_agegroup)
        end
    end
    alpha[1, 1][6] = 0.925
    alpha[1, 3][6] = 0.075
    alpha[1, 1][8] = 0.889
    alpha[1, 3][8] = 0.111
    alpha[2, 1] = copy(alpha[1, 1])
    alpha[3, 1] = copy(alpha[1, 1])
    alpha[4, 1] = copy(alpha[1, 1])
    alpha[2, 2] = copy(alpha[1, 3])
    alpha[3, 3] = copy(alpha[1, 3])
    alpha[4, 4] = copy(alpha[1, 3])

    beta = Array{Vector{Float64}}(undef, 4, 4)# each element=matrix row (year)*matrix column (age group) 
    for j = 1:4
        for i = 1:4
            beta[i, j] = zeros(n_agegroup)
        end
    end
    beta[1, 1][2] = 0.857
    beta[1, 3][2] = 0.143
    beta[1, 1][3] = 0.928
    beta[1, 3][3] = 0.072
    beta[2, 1] = copy(beta[1, 1])
    beta[3, 1] = copy(beta[1, 1])
    beta[4, 1] = copy(beta[1, 1])
    beta[2, 2] = copy(beta[1, 3])
    beta[3, 3] = copy(beta[1, 3])
    beta[4, 4] = copy(beta[1, 3])

    #vaccine coverage
    rate_DTaP_2 = [5.53; 4.21; 4.21; 3.86; 3.59; 2.89; 4.21; 2.36; 1.73; 1.73; 2.36; 4.21; 4.69; 2.89; 2.89; 4.69; 4.69]
    rate_DTaP_3 = [1.59; 1.59; 1.59; 1.59; 1.59; 1.43; 1.49; 1.43; 1.29; 1.88; 1.41; 1.84; 1.87; 1.95; 1.95; 2.66; 2.66]
    psi = psifun(rate_DTaP_2, rate_DTaP_3)


    coverage_adolescent = [0.627; 0.627; 0.627; 0.627; 0.627; 0.689; 0.665; 0.677; 0.699; 0.604; 0.574; 0.65; 0.632; 0.721; 0.741; 0.762; 0.762]
    coverage_adult = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0.093; 0.0965; 0.1; 0.2155; 0.331; 0.331; 0.34]
    phi = phifun(coverage_adolescent, coverage_adult)

    #n_statvar=12*n_agegroup; #number of state variables


    sigma = [theta[1] * ones(3); theta[2] * theta[1] * ones(2); theta[3] * theta[1] * ones(10); theta[4] * theta[1] * ones(5)]
    sigmar = theta[5] * sigma
    rho = theta[6] * ones(20)
    p = [theta[7] * ones(3); theta[7] * theta[8] * ones(2); theta[7] * theta[9] * ones(10); theta[7] * theta[10] * ones(5)]

    #infection due to travel
    rate = zeros(n_agegroup)
    rate[4:15] .= theta[11]
    rate[16:end] .= theta[12]
    epsilon_p = theta[13]
    epsilon_t = theta[13]

    #
    # Variables
    S = StatVar[1:n_agegroup]
    ED = StatVar[n_agegroup+1:2*n_agegroup]
    D = StatVar[2*n_agegroup+1:3*n_agegroup]
    H = StatVar[3*n_agegroup+1:4*n_agegroup]
    SR = StatVar[4*n_agegroup+1:5*n_agegroup]
    EU = StatVar[5*n_agegroup+1:6*n_agegroup]
    U = StatVar[6*n_agegroup+1:7*n_agegroup]
    RC = StatVar[7*n_agegroup+1:8*n_agegroup]
    RI = StatVar[8*n_agegroup+1:9*n_agegroup]
    RD = StatVar[9*n_agegroup+1:10*n_agegroup]
    VC = StatVar[10*n_agegroup+1:11*n_agegroup]
    VD = StatVar[11*n_agegroup+1:12*n_agegroup]
    N = S + ED + D + H + SR + EU + U + RC + RI + RD + VC + VD


    # Contact Matrix function
    contactmatpath = "DataSource/C20_Y_2004_to_2040.mat"
    contactmatstr = "C20_Y_2004_to_2040"
    contactfunc = contactmatrix(contactmatpath, contactmatstr)
    # Demo rate
    demopath = "DataSource/demo_para_ontario.mat"
    #demopathmonthly = "DataSource/demo_para_ontario_monthly.mat"
    Years_demo = demorate(demopath, "Years_deaths_ontario")
    r = [6.0; 1.2; 0.25; ones(16); 0]
    B = demorate(demopath, "births_table_ontario")
    # nu_monthly=demorate(demopathmonthly,"PC_Migration_Rates")

    birthfitness = [1; zeros(n_agegroup - 1)]

    for i in 1:n_agegroup
        λi = p[i] * (1 + epsilon_p * cos(2 * pi * (t - Years[1]) + pi)) * sum([365 * contactfunc[i, j](t) * (D[j] + rho[j] * U[j]) / N[j] for j in 1:n_agegroup])
        traveli = rate[i] * (1 + epsilon_t * cos(2 * pi * (t - Years[1]) + pi))
        nui = demorate(demopath, "PC_Migration_Rates", startyear=2004, age=i)
        mui = demorate(demopath, "Death_Rates_calc", startyear=2004, age=i)

        # S
        inputSi = i > 1 ? (1 - phi[i](t)) * r[i-1] * S[i-1] : 0
        rhs[i] = birthfitness[i] * B(t) * (1 - phim(t)) + inputSi + (nui(t) - mui(t) - λi - r[i] - psi[i](t) - traveli) * S[i]

        # Define E^D equations
        inputEDit = i > 1 ? r[i-1] * ED[i-1] : 0
        rhs[n_agegroup+i] = inputEDit + sigma[i] * λi * S[i] + sigmar[i] * λi * SR[i] - eta_EDtoD[i] * ED[i] + (nui(t) - mui(t) - r[i]) * ED[i]

        # Define D equations
        inputDit = i > 1 ? r[i-1] * D[i-1] : 0
        rhs[2*n_agegroup+i] = inputDit + sigma[i] * traveli * S[i] + eta_EDtoD[i] * ED[i] - (gamma_DtoH[i] + gamma_DtoRC[i]) * D[i] + (nui(t) - mui(t) - r[i]) * D[i]

        # Define H equations
        inputHit = i > 1 ? r[i-1] * H[i-1] : 0
        rhs[3*n_agegroup+i] = inputHit + gamma_DtoH[i] * D[i] - gamma_HtoRC[i] * H[i] + (nui(t) - mui(t) - r[i]) * H[i]

        # Define S^R equations
        inputSRit = i > 1 ? (1 - phi[i](t)) * r[i-1] * SR[i-1] + phi[i](t) * r[i-1] * (alpha[1, 3][i-1] * S[i-1] + alpha[3, 3][i-1] * SR[i-1]) : 0
        rhs[4*n_agegroup+i] = birthfitness[i] * B(t) * phim(t) * alpham[3] + inputSRit + psi[i](t) * (beta[1, 3][i] * S[i] + beta[3, 3][i] * SR[i]) + (nui(t) - mui(t) - r[i] - λi - psi[i](t)) * SR[i] + tau_RDtoSr[i] * RD[i] + tau_VDtoSr[i] * VD[i]

        # Define E^U equations
        inputEUit = i > 1 ? r[i-1] * EU[i-1] : 0
        rhs[5*n_agegroup+i] = inputEUit + (1 - sigma[i]) * λi * S[i] + (1 - sigmar[i]) * λi * SR[i] + λi * RD[i] + λi * VD[i] - eta_EUtoU[i] * EU[i] + (nui(t) - mui(t) - r[i]) * EU[i]

        # Define U equations
        inputUit = i > 1 ? r[i-1] * U[i-1] : 0
        rhs[6*n_agegroup+i] = inputUit + (1 - sigma[i]) * traveli * S[i] + eta_EUtoU[i] * EU[i] - gamma_UtoRC[i] * U[i] + (nui(t) - mui(t) - r[i]) * U[i]

        # Define R^C equations
        inputRCit = i > 1 ? r[i-1] * RC[i-1] : 0
        rhs[7*n_agegroup+i] = inputRCit + gamma_DtoRC[i] * D[i] + gamma_UtoRC[i] * U[i] + gamma_HtoRC[i] * H[i]
        +λi * RI[i] - tau_RCtoRI[i] * RC[i] + (nui(t) - mui(t) - r[i]) * RC[i]

        # Define R^I equations
        inputRIit = i > 1 ? r[i-1] * RI[i-1] : 0
        rhs[8*n_agegroup+i] = inputRIit + tau_RCtoRI[i] * RC[i] - (tau_RItoRD[i] + λi) * RI[i] + (nui(t) - mui(t) - r[i]) * RI[i]

        # Define R^D equations
        inputRDit = i > 1 ? ((1 - phi[i](t)) + phi[i](t) * alpha[4, 4][i-1]) * r[i-1] * RD[i-1] : 0
        rhs[9*n_agegroup+i] = inputRDit + beta[4, 4][i] * psi[i](t) * RD[i] + tau_RItoRD[i] * RI[i] - tau_RDtoSr[i] * RD[i] - λi * RD[i] - psi[i](t) * RD[i] + (nui(t) - mui(t) - r[i]) * RD[i]


        # Define V^C equations
        inputVCit = i > 1 ? r[i-1] * VC[i-1] + r[i-1] * phi[i](t) * (alpha[1, 1][i-1] * S[i-1] + alpha[2, 1][i-1] * VD[i-1] + alpha[3, 1][i-1] * SR[i-1] + alpha[4, 1][i-1] * RD[i-1]) : 0
        rhs[10*n_agegroup+i] = birthfitness[i] * B(t) * phim(t) * alpham[1] + inputVCit + psi[i](t) * (beta[1, 1][i] * S[i] + beta[2, 1][i] * VD[i] + beta[3, 1][i] * SR[i] + beta[4, 1][i] * RD[i]) - tau_VCtoVD[i] * VC[i] + (nui(t) - mui(t) - r[i]) * VC[i]


        # Define V^D equations
        inputVDit = i > 1 ? (1 - phi[i](t)) * r[i-1] * VD[i-1] + r[i-1] * phi[i](t) * (alpha[1, 2][i-1] * S[i-1] + alpha[2, 2][i-1] * VD[i-1] + alpha[3, 2][i-1] * SR[i-1] + alpha[4, 2][i-1] * RD[i-1]) : 0
        rhs[11*n_agegroup+i] = birthfitness[i] * B(t) * phim(t) * alpham[2] + inputVDit + psi[i](t) * (beta[1, 2][i] * S[i] + beta[2, 2][i] * VD[i] + beta[3, 2][i] * SR[i] + beta[4, 2][i] * RD[i]) + tau_VCtoVD[i] * VC[i] - tau_VDtoSr[i] * VD[i] + (nui(t) - mui(t) - r[i] - λi - psi[i](t)) * VD[i]
    end
    return rhs
end