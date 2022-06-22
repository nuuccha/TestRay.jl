""" Defined as in Mikhelson "optics for astronomical telescopes" p93 
    Ray.x  Ray.y Ray.z is the position related to the origin 
    Ray.cx Ray.cy Ray.cz vector cosines, Ray.n is the refraction index
    if Ray.empty == 0 then the ray is processed, otherwise ignored"""
struct Ray   
    x::Float64  
    y::Float64
    z::Float64
    cx::Float64 
    cy::Float64
    cz::Float64
    n::Float64 
    empty::Int32  
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

