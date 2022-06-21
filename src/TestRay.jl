module TestRay
using Plots
using Statistics

include("structures.jl")

include("routines.jl")

s=Surf(33.97, 20000.0, 1.62166)
#s=Surf(33.97, 2000.0, 1.0)
s1=Surf(0.0, 6.1878, 1.0)
s2=Surf(-62.67, 6.9862, 1.579)
s3=Surf(29.58, 2.495, 1.0)
s4=Surf(1625.0, 7.4852, 1.5791)
s5=Surf(27.34, 2.1957, 1.6577)
s6=Surf(-46.366, 8.189, 1.0)
s7=Surf(0.0, 83.732, 1.0 )
#=
s=Surf(0.0, 20000.0, 1.0)
s1=Surf(0, 5, 1.0)
s2=Surf(0, 97, 1.0)
=#
OptSys = (s,s1,s2,s3,s4,s5,s6,s7)

l=Surf(0, 20000.0, 1.5)
l1=Surf(-50, 20.0, 1.0)
l2=Surf(-20, 100., 1.0)
OptSys=(l, l1, l2)
#s=Surf(0.0, 20000.0, 1.0)
#OptSys = (s,)
entr_pupil = Pupil(20000.0-70, 10.0 ,0.0 ,0.0)
pencil = make__hex_pencil(1000.0, 0.0, entr_pupil, 41)
#=
#example of a composite pencil made of several pencils
pencil2 = make__hex_pencil(0.0, -1000.0, entr_pupil, 3)
pencil1 = make__hex_pencil(0.0, 0.0, entr_pupil, 3)
pencil3 = make__hex_pencil(1000.0, 0.0, entr_pupil, 3)
pencil = vcat(pencil2,pencil1)
pencil=vcat(pencil,pencil3)
=#

# Defining the focus position by diff propaating two axial rays 

    println("Focus position = ",focus_position(pencil,OptSys))
    
    #Propagating the pencil through the system:
    @time Propagate(pencil, OptSys)  # does work !!!
    
    # extracting the ray coordinates
    xx = Array{Float64}(undef,size(pencil))
    yy = Array{Float64}(undef,size(pencil))
        for i in 1:size(pencil)[1]
            #println(pencil[i,j])
            rr = pencil[i]
             xx[i] = rr.y
             yy[i] = rr.z
        end
    println(Pencil_rms(pencil))
    (display(scatter(xx, yy)))
    

end # module
