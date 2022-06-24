abstract type Glass end
struct AbbeGlass <: Glass
    nr::Float64
    nu::Float64
    typ::String
end

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
    empty::Number  # set to 0 to process the ray, set to 1 to ignore the ray
end

""" Surface is an abstract type with members  
such as Surf, Paraxial, etc. (to be extended) With this class
the optical system can be defined as an array of abstract surfaces, 
with the specifics defined on the fly """ 
abstract type Surface end

""" normal refracting surface """
mutable struct Surf <: Surface  # Defined as in Mikhelson "optics for astronomical telescopes" p94
    r::Float64  # ROC
    d::Float64  # thickness BEFORE surface
    n::Float64   # refractive index
    e2::Float64  #  e square (excentricity)
    stop::Number  # set to 1 for a stop
end
""" Refracting surface with AbbeGlass """
struct Abbe_Surf <: Surface  # Defined as in Mikhelson "optics for astronomical telescopes" p94
    r::Float64  # ROC
    d::Float64  # thickness BEFORE surface
    glass::AbbeGlass   # refractive index
    lambda::Number
    e2::Float64  #  e square (excentricity)
    stop::Number  # set to 1 for a stop
end
""" An attempt of paraxial surface, tricky, "sine condition???" """
struct Paraxial  <: Surface # paraxial lens
    f::Float64  # focal length   
end
""" Separate structure for pupil, has to be defined by ray propagation though the 
front and back part of the system """
struct Pupil <: Surface
    distance::Float64
    radius::Float64
    y::Float64
    z::Float64
end


""" Glasses with Abber Numbers """
BK7=AbbeGlass(1.51680, 64.17, "D")
F2=AbbeGlass(1.62004, 36.37, "D")
F5=AbbeGlass(1.60342, 38.03, "D")
K10=AbbeGlass(1.5014, 56.409, "D")
LAFN7=AbbeGlass(1.7495, 34.95, "D")
NF2=AbbeGlass(1.62005, 36.43, "D")
NSF6=AbbeGlass(1.80518, 25.36, "D")
NPK52A=AbbeGlass(1.497, 81.61, "D")
AIR=AbbeGlass(1.0, 10000., "D")
SF6=AbbeGlass(1.8052, 25.432, "D")
