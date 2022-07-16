@doc raw"""
    contactmatrix(matpath,str)

Read contact matrix mat file by using `MAT.jl`, and then using `DataInterpolations.jl` to obtain the
contact matrix function.

    mathpath is the path of the mat file. str is the key of the dictionary.
"""
function contactmatrix(matpath, str)
    datasource = matread(matpath)
    datasourcematrix = datasource[str]
    yearlen = size(datasourcematrix)[2]
    matrixsize = size(datasourcematrix[1])
    contactmatrix = Array[Function](undef, matrixsize)
    for j in 1:matrixsize[2]
        for i in 1:matrixsize[1]
            contactarray = [datasourcematrix[k][i, j] for k in 1:yearlen]
            timespan = 0:(yearlen-1)
            contactmatrix[i, j] = x -> QuadraticInterpolation(contactarray, timespan)(x)
        end
    end
    return contactmatrix
end

@doc raw"""
    demorate(demopath,str)

Read birth death age migration rate mat file by using `MAT.jl`, and then using `DataInterpolations.jl` to obtain the
birth rate function.

    demopath is the path of the mat file. str is the key of the dictionary.
"""
function demorate(demopath, str)
    datasource = matread(demopath)
    datasourcevector = datasource[str][:, 1]
    timespan = 0:(length(datasourcevector)-1)
    demorate = x -> QuadraticInterpolation(datasourcevector, timespan)(x)
    return demorate
end

@doc raw"""
    pertussisage(du, u, p, t)

main model of pertussis.
"""
function pertussisage(rhs, StatVar, theta, t)
 beta
 phim phi psi

# Contact Matrix function
contactmatpath = "DataSource/C20_Y_2004_to_2040.mat"
contactmatstr = "C20_Y_2004_to_2040"
contactfunc = contactmatrix(contactmatpath, contactmatstr)
# Demo rate
demopath = "DataSource/demo_para_ontario.mat"
demopathmonthly = "DataSource/demo_para_ontario_monthly.mat"
Years_demo=demorate(demopath,"Years_deaths_ontario") 
r=demorate(demopath,"Aging_Rates")
B=demorate(demopath,"births_table_ontario")
mu=demorate(demopath,"Death_Rates_calc")
nu=demorate(demopath,"PC_Migration_Rates")
nu_monthly=demorate(demopathmonthly,"PC_Migration_Rates")


agegroup_str=["0-1m";"2-11m";"1-4y";"5-9y";"10-14y";"15-19y";"20-24y";
    "25-29y";"30-34y";"35-39y";"40-44y";"45-49y";"50-54y";"55-59y";"60-64y";
    "65-69y";"70-74y";"75-79y";"80-84y";"85+y"]#age group division   
n_agegroup=length(agegroup_str)#number of age groups
Years=collect(2011:1:2020)#years considered in this study
Years_c=collect(2004:1:2021)
paramstr = [L"\sigma_1"; L"\sigma_[2p]"; L"\sigma_[3p]";L"\sigma_[4p]"; L"\sigma^r_scale";
    L"\rho"; L"p_1"; L"p_[2p]"; L"p_[3p]"; L"p_[4p]"; L"t_2"; L"t_3"; L"\epsilon_p"]
nparam=length(paramstr)

#disease related parameter values
eta_EDtoD=365/9*ones(n_agegroup,1);
eta_EUtoU=365/8*ones(n_agegroup,1);
gamma_DtoH=365/15.6*ones(n_agegroup,1);
gamma_HtoRC=365/5*ones(n_agegroup,1);
gamma_DtoRC=365/21*ones(n_agegroup,1);
gamma_UtoRC=365/15*ones(n_agegroup,1);

#waning of natural immunity
tau_RCtoRI=1/10*ones(n_agegroup,1);
tau_RItoRD=1/15*ones(n_agegroup,1);
tau_RDtoSr=1/5*ones(n_agegroup,1);

#vaccine-related parameter values
#waning of immunity
Years_vacc=collect(2004:1:2020)
tau_VCtoVD=ones(n_agegroup)
tau_VCtoVD[1]=1/1.48
tau_VCtoVD[3:5].=1/4.955
tau_VCtoVD[6:end].=1/4.02
tau_VDtoSr=copy(tau_VCtoVD)

#vaccine efficacy
alpham=[0.914;0;0.086];# Materanl booster efficacy

alpha=Array[any](undef,4,4);# each element=row (year)*column (age group) 
for j=1:4
    for i=1:4
        alpha[i,j]=zeros(length(Years_vacc),n_agegroup);
    end
end
alpha[1,1][:,6].=0.925;
alpha[1,3][:,6].=0.075;
alpha[1,1][:,8].=0.889;
alpha[1,3][:,8].=0.111;
alpha[2,1]=copy(alpha[1,1]);
alpha[3,1]=copy(alpha[1,1]);
alpha[4,1]=copy(alpha[1,1]);
alpha[2,2]=copy(alpha[1,3]);
alpha[3,3]=copy(alpha[1,3]);
alpha[4,4]=copy(alpha[1,3]);

beta=Array[any](undef,4,4);# each element=matrix row (year)*matrix column (age group) 
for j=1:4
    for i=1:4
        beta[i,j]=zeros(length(Years_vacc),n_agegroup);
    end
end
beta[1,1][:,2].=0.857;
beta[1,3][:,2].=0.143;
beta[1,1][:,3].=0.928;
beta[1,3][:,3].=0.072;
beta[2,1]=copy(beta[1,1]);
beta[3,1]=copy(beta[1,1]);
beta[4,1]=copy(beta[1,1]);
beta[2,2]=copy(beta[1,3]);
beta[3,3]=copy(beta[1,3]);
beta[4,4]=copy(beta[1,3]);

#vaccine coverage
phim=zeros(length(Years_vacc),1);#maternal booster
phim(find(Years_vacc==2018))=0.4;
phim(find(Years_vacc>=2019&Years_vacc<=2020))=0.395;

psi=zeros(length(Years_vacc),n_agegroup);
rate_DTaP_2=[5.53;4.21;4.21;3.86;3.59;2.89;4.21;2.36;1.73;1.73;2.36;...
    4.21;4.69;2.89;2.89;4.69;4.69];
rate_DTaP_3=[1.59;1.59;1.59;1.59;1.59;1.43;1.49;1.43;1.29;1.88;1.41;...
    1.84;1.87;1.95;1.95;2.66;2.66];
psi(find(Years_vacc>=2004&Years_vacc<=2020),2)=repmat(rate_DTaP_2,1);
psi(find(Years_vacc>=2004&Years_vacc<=2020),3)=repmat(rate_DTaP_3,1);

phi=zeros(length(Years_vacc),n_agegroup);
coverage_adolescent=[0.627;0.627;0.627;0.627;0.627;0.689;0.665;0.677;...
    0.699;0.604;0.574;0.65;0.632;0.721;0.741;0.762;0.762];
coverage_adult=[0;0;0;0;0;0;0;0;0;0;0.093;0.0965;0.1;0.2155;0.331;0.331;0.34];
phi(find(Years_vacc>=2004&Years_vacc<=2020),6)=repmat(coverage_adolescent,1);
phi(find(Years_vacc>=2004&Years_vacc<=2020),8)=repmat(coverage_adult,1);



   
n_statvar=12*n_agegroup; #number of state variables
rhs=zeros(n_statvar,length(time));
lambda=cell(length(time),1);# Define lambda

sigma=[theta(1)*ones(3,1);theta(2)*theta(1)*ones(2,1);theta(3)*theta(1)*ones(10,1);theta(4)*theta(1)*ones(5,1)];
sigmar=theta(5)*sigma;
rho=theta(6)*ones(20,1);
p=[theta(7)*ones(3,1);theta(7)*theta(8)*ones(2,1);theta(7)*theta(9)*ones(10,1);theta(7)*theta(10)*ones(5,1)];

#infection due to travel
rate=zeros(1,n_agegroup);
rate(1,4:15)=theta(11);
rate(1,16:end)=theta(12);
epsilon_p=theta(13);

#
# Variables
S=StatVar(1:n_agegroup);
ED=StatVar(n_agegroup+1:2*n_agegroup);
D=StatVar(2*n_agegroup+1:3*n_agegroup);
H=StatVar(3*n_agegroup+1:4*n_agegroup);
Sr=StatVar(4*n_agegroup+1:5*n_agegroup);
EU=StatVar(5*n_agegroup+1:6*n_agegroup);
U=StatVar(6*n_agegroup+1:7*n_agegroup);
RC=StatVar(7*n_agegroup+1:8*n_agegroup);
RI=StatVar(8*n_agegroup+1:9*n_agegroup);
RD=StatVar(9*n_agegroup+1:10*n_agegroup);
VC=StatVar(10*n_agegroup+1:11*n_agegroup);
VD=StatVar(11*n_agegroup+1:12*n_agegroup);
N=S+ED+D+H+Sr+EU+U+RC+RI+RD+VC+VD;

# Define equations
for i=1:length(time) 
    t = time(i);
    
    # Define lambda
    lambda[i]=zeros(n_agegroup,1);
    C=cell2mat(c(find(Years_c<=time(i),1,"last")));
    for j=1:n_agegroup
        lambda[i](j)=p(j)*(1+epsilon_p*cos(2*pi*(t-Years(1))+pi))*sum((365*C(j,:))".*(D+rho.*U)./N);
    end
    
    # Define S equations
    rhs(1:n_agegroup,i)=[1;zeros(n_agegroup-1,1)]*time_rates(t,Years_demo,B)...
        *(1-time_rates(t,Years_vacc,phim))...
        +diag(-r"-time_rates(t,Years_demo,mu)+time_rates(t,Years_demo,nu)...
        -lambda[i]"-time_rates(t,Years_vacc,psi))*S...
        +diag(r(1:end-1)".*(ones(1,n_agegroup-1)-time_rates(t,Years_vacc,phi(:,2:end))),-1)*S...
        -diag(rate)*S;
    
    # Define E^D equations
    rhs(n_agegroup+1:2*n_agegroup,i)=diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-eta_EDtoD")*ED...
        +diag(r(1:end-1),-1)*ED...
        +diag(sigma".*lambda[i]")*S...
        +diag(sigmar".*lambda[i]")*Sr;
    
    # Define D equations
    rhs(2*n_agegroup+1:3*n_agegroup,i)=diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-gamma_DtoH"-gamma_DtoRC")*D...
        +diag(r(1:end-1),-1)*D...
        +diag(eta_EDtoD)*ED...
        +diag(sigma".*rate)*S;
    
    # Define H equations
    rhs(3*n_agegroup+1:4*n_agegroup,i)=diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-gamma_HtoRC")*H...
        +diag(r(1:end-1),-1)*H...
        +diag(gamma_DtoH)*D;
 
    # Define S^r equations
    rhs(4*n_agegroup+1:5*n_agegroup,i)=[1;zeros(n_agegroup-1,1)]*time_rates(t,Years_demo,B)...
        *time_rates(t,Years_vacc,phim)*alpham(3)...
        +diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-time_rates(t,Years_vacc,psi)-lambda[i]"...
        +time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[3,3]))*Sr...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[1,3]))*S...
        +diag(tau_RDtoSr)*RD...
        +diag(time_rates(t,Years_vacc,tau_VDtoSr))*VD...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[1,3](:,2:end)),-1)*S...
        +diag(r(1:end-1)".*(ones(1,n_agegroup-1)-time_rates(t,Years_vacc,phi(:,2:end)))...
        +r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[3,3](:,2:end)),-1)*Sr;
    
    # Define E^U equations
    rhs(5*n_agegroup+1:6*n_agegroup,i)=diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-eta_EUtoU")*EU...
        +diag(r(1:end-1),-1)*EU...
        +diag((ones(n_agegroup,1)-sigma)".*lambda[i]")*S...
        +diag((ones(n_agegroup,1)-sigmar)".*lambda[i]")*Sr...
        +diag(lambda[i]")*(RD+VD);
    
    # Define U equations
    rhs(6*n_agegroup+1:7*n_agegroup,i)=diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-gamma_UtoRC")*U...
        +diag(r(1:end-1),-1)*U...
        +diag(eta_EUtoU)*EU...
        +diag((ones(n_agegroup,1)-sigma)".*rate)*S;
    
    # Define R^C equations
    rhs(7*n_agegroup+1:8*n_agegroup,i)=diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-tau_RCtoRI")*RC...
        +diag(r(1:end-1),-1)*RC...
        +diag(gamma_DtoRC)*D+diag(gamma_UtoRC)*U+diag(gamma_HtoRC)*H...
        +diag(lambda[i])*RI;
    
    # Define R^I equations
    rhs(8*n_agegroup+1:9*n_agegroup,i)=diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-tau_RItoRD"-lambda[i]")*RI...
        +diag(r(1:end-1),-1)*RI...
        +diag(tau_RCtoRI)*RC;
    
    # Define R^D equations
    rhs(9*n_agegroup+1:10*n_agegroup,i)=diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-tau_RDtoSr"-lambda[i]"...
        -time_rates(t,Years_vacc,psi)+time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[4,4]))*RD...
        +diag(tau_RItoRD)*RI... 
        +diag(r(1:end-1)".*(ones(1,n_agegroup-1)-time_rates(t,Years_vacc,phi(:,2:end)))...
        +r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[4,4](:,2:end)),-1)*RD;

    
    # Define V^C equations
    rhs(10*n_agegroup+1:11*n_agegroup,i)=[1;zeros(n_agegroup-1,1)]*time_rates(t,Years_demo,B)...
        *time_rates(t,Years_vacc,phim)*alpham(1)...
        +diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-time_rates(t,Years_vacc,tau_VCtoVD))*VC...
        +diag(r(1:end-1),-1)*VC...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[1,1]))*S...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[2,1]))*VD...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[3,1]))*Sr...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[4,1]))*RD...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[1,1](:,2:end)),-1)*S...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[2,1](:,2:end)),-1)*VD...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[3,1](:,2:end)),-1)*Sr...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[4,1](:,2:end)),-1)*RD;
    
    # Define V^D equations
    rhs(11*n_agegroup+1:12*n_agegroup,i)=[1;zeros(n_agegroup-1,1)]*time_rates(t,Years_demo,B)...
        *time_rates(t,Years_vacc,phim)*alpham(2)...
        +diag(-r"-time_rates(t,Years_demo,mu)...
        +time_rates(t,Years_demo,nu)-time_rates(t,Years_vacc,tau_VDtoSr)...
        -lambda[i]"-time_rates(t,Years_vacc,psi))*VD...
        +diag(r(1:end-1)".*(ones(1,n_agegroup-1)-time_rates(t,Years_vacc,phi(:,2:end))),-1)*VD...
        +diag(time_rates(t,Years_vacc,tau_VCtoVD))*VC...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[1,2]))*S...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[2,2]))*VD...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[3,2]))*Sr...
        +diag(time_rates(t,Years_vacc,psi).*time_rates(t,Years_vacc,beta[4,2]))*RD...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[1,2](:,2:end)),-1)*S...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[2,2](:,2:end)),-1)*VD...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[3,2](:,2:end)),-1)*Sr...
        +diag(r(1:end-1)".*time_rates(t,Years_vacc,phi(:,2:end)).*time_rates(t,Years_vacc,alpha[4,2](:,2:end)),-1)*RD;
end
end