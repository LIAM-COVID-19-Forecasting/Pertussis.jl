using Pertussis
using DifferentialEquations
using Plots

# Contactmatrix
contactmatpath = "DataSource/C20_Y_2004_to_2040.mat"

contactmatstr = "C20_Y_2004_to_2040"

contactfunc = contactmatrix(contactmatpath, contactmatstr)

size(contactfunc)

x = 2021:0.1:2050

plot(x, contactfunc[1, 4].(x))

AA = [contactfunc[i, j](2004) for i in 1:20, j in 1:20]

heatmap(AA)

using MAT

heatmap(matread(contactmatpath)[contactmatstr][2])

AA == matread(contactmatpath)[contactmatstr][1]

ii = 3
matread(contactmatpath)[contactmatstr][ii] == [contactfunc[i, j](ii + 2003) for i in 1:20, j in 1:20]


# psifun, phifun

rate_DTaP_2 = [5.53; 4.21; 4.21; 3.86; 3.59; 2.89; 4.21; 2.36; 1.73; 1.73; 2.36; 4.21; 4.69; 2.89; 2.89; 4.69; 4.69]
rate_DTaP_3 = [1.59; 1.59; 1.59; 1.59; 1.59; 1.43; 1.49; 1.43; 1.29; 1.88; 1.41; 1.84; 1.87; 1.95; 1.95; 2.66; 2.66]
psi = psifun(rate_DTaP_2, rate_DTaP_3)
x = 2004:0.1:2020
plot(x, psi[2].(x))

coverage_adolescent = [0.627; 0.627; 0.627; 0.627; 0.627; 0.689; 0.665; 0.677; 0.699; 0.604; 0.574; 0.65; 0.632; 0.721; 0.741; 0.762; 0.762]
coverage_adult = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0.093; 0.0965; 0.1; 0.2155; 0.331; 0.331; 0.34]
phi = phifun(coverage_adolescent, coverage_adult)
x = 2004:0.1:2020
plot(x, phi[6].(x))
