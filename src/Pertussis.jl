module Pertussis

# Write your package code here.
using MAT
using DataInterpolations
include("utils.jl")
include("models.jl")

export pertussisage
export contactmatrix
export demorate
export phifun
export psifun

end
