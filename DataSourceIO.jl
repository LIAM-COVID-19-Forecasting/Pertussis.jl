using MAT
using DataFrames
using Plots
using CSV
file = matread("./DataSource/C20_Y_2004_to_2021.mat")
file2 = matread("DataSource/C20_Y_2004_to_2040.mat")
file3 = matread("DataSource/ontario_incidence.mat")
file4 = matread("DataSource/demo_data_ontario.mat")
file5 = matread("DataSource/demo_para_ontario.mat")
file3["Ontario_yearly_incidence"]
file3["Ontario_monthly_incidence"]
file4["deaths_table_ontario"]

plot(file3["Ontario_monthly_incidence"])

scatter(file3["Ontario_monthly_incidence"][:, 1])

cols = ["0_4", "5_9", "10_14", "15_19", "20_49", "50_64", "65+"]

# incidence yearly
incidence_year = DataFrame(file3["Ontario_yearly_incidence"], cols)
CSV.write("DataSource/incidence_year.csv", incidence_year)

# incidence monthly 
incidence_month = DataFrame(file3["Ontario_monthly_incidence"], cols)
CSV.write("DataSource/incidence_month.csv", incidence_month)

