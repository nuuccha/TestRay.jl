
struct Ray  # Defined as in Mikhelson "optics for astronomical telescopes" p93
    x::Float64  #Position related to the origin
    y::Float64
    z::Float64
    cx::Float64 # direction cosines
    cy::Float64
    cz::Float64
    n::Float64 # refraction index of current
end

struct Surf  # Defined as in Mikhelson "optics for astronomical telescopes" p94
    r::Float64  # ROC
    d::Float64  # thickness BEFORE surface
    n::Float64   # refractive index
#    e2::Float64  #  e square (excentricity)
end

struct Paraxial  # paraxial lens
    f::Float64  # focal length   
end

struct Pupil
    distance::Float64
    radius::Float64
    y::Float64
    z::Float64
end

struct Glass
    nr::Float64
    nu::Float64
    typ::String
end

BK7=Glass(1.51680, 64.17, "D")
F2=Glass(1.62004, 36.37, "D")
F5=Glass(1.60342, 38.03, "D")
K10=Glass(1.50137, 60.41, "D")
LAFN7=Glass(1.7495, 34.95, "D")
NF2=Glass(1.62005, 36.43, "D")
NSF6=Glass(1.80518, 25.36, "D")
NPK52A=Glass(1.497, 81.61, "D")

