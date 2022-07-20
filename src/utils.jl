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
    contactmatrix = Array{Function,2}(undef, matrixsize)
    for j in 1:matrixsize[2]
        for i in 1:matrixsize[1]
            contactarray = [datasourcematrix[k][i, j] for k in 1:yearlen]
            timespan = 2004:(2003+yearlen)
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
function demorate(demopath, str; startyear=2011, age=1)
    datasource = matread(demopath)
    datasourcevector = datasource[str][:, age]
    timespan = startyear:(length(datasourcevector)+startyear-1)
    demorate = x -> QuadraticInterpolation(datasourcevector, timespan)(x)
    return demorate
end

@doc raw"""
    vaccrate(demopath,str)

Read birth death age migration rate mat file by using `MAT.jl`, and then using `DataInterpolations.jl` to obtain the
birth rate function.

    demopath is the path of the mat file. str is the key of the dictionary.
"""
function vaccrate(demopath, str) end


@doc raw"""
    phim(t)

maternal booster
"""
function phim(t)
    if t < 2018
        y = 0
    elseif t >= 2018 && t < 2019
        y = 0.4
    else
        y = 0.395
    end
    return y
end

@doc raw"""
    phifun()

coverage_adolescent=[0.627;0.627;0.627;0.627;0.627;0.689;0.665;0.677;...
    0.699;0.604;0.574;0.65;0.632;0.721;0.741;0.762;0.762];
coverage_adult=[0;0;0;0;0;0;0;0;0;0;0.093;0.0965;0.1;0.2155;0.331;0.331;0.34];
"""
function phifun(coverage_adolescent, coverage_adult)
    phimatrix = Array{Function,1}(undef, 20)
    for i in 1:20
        phimatrix[i] = x -> 0
    end
    yearlen = length(coverage_adolescent)
    timespan = 2004:(2003+yearlen)
    phimatrix[6] = x -> QuadraticInterpolation(coverage_adolescent, timespan)(x)
    phimatrix[8] = x -> QuadraticInterpolation(coverage_adult, timespan)(x)
    return phimatrix
end

@doc raw"""
    psifun()

rate_DTaP_2=[5.53;4.21;4.21;3.86;3.59;2.89;4.21;2.36;1.73;1.73;2.36;...
    4.21;4.69;2.89;2.89;4.69;4.69];
rate_DTaP_3=[1.59;1.59;1.59;1.59;1.59;1.43;1.49;1.43;1.29;1.88;1.41;...
    1.84;1.87;1.95;1.95;2.66;2.66];
"""
function psifun(rate_DTaP_2, rate_DTaP_3)
    psimatrix = Array{Function,1}(undef, 20)
    for i in 1:20
        psimatrix[i] = x -> 0
    end
    yearlen = length(rate_DTaP_2)
    timespan = 2004:(2003+yearlen)
    psimatrix[2] = x -> QuadraticInterpolation(rate_DTaP_2, timespan)(x)
    psimatrix[3] = x -> QuadraticInterpolation(rate_DTaP_3, timespan)(x)
    return psimatrix
end

