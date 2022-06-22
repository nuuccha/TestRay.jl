module TestRay
using Plots
using Statistics
using Optim

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
entr_pupil = Pupil(20000.0, 15.0 ,0.0 ,0.0)
xx = Array{Float64}(undef, 6)
xx[1] = 100
xx[2] = -100
xx[3] = 100
xx[4] = 100
xx[5] = 100
xx[6] = -100
xx1 = copy(xx)

function FF(xx::Array{Float64})
    pencil = make__hex_pencil(3000.0, 0.0, entr_pupil, 21)
    s=Surf(xx[1], 20000.0, 1.62166)
    s1=Surf(0.0, 6.1878, 1.0)
    s2=Surf(xx[2], 6.9862, 1.579)
    s3=Surf(xx[3], 2.495, 1.0)
    s4=Surf(xx[4], 7.4852, 1.5791)
    s5=Surf(xx[5], 2.1957, 1.6577)
    s6=Surf(xx[6], 8.189, 1.0)
    s7=Surf(0.0, 85, 1.0 )
    OptSys1 = (s,s1,s2,s3,s4,s5,s6,s7)
    Propagate(pencil, OptSys1)  # does work !!!
    rr1 = (Pencil_rms(pencil))[1]
    pencil = make__hex_pencil(0.0, 0.0, entr_pupil, 21)
    s=Surf(xx[1], 20000.0, 1.62166)
    s1=Surf(0.0, 6.1878, 1.0)
    s2=Surf(xx[2], 6.9862, 1.579)
    s3=Surf(xx[3], 2.495, 1.0)
    s4=Surf(xx[4], 7.4852, 1.5791)
    s5=Surf(xx[5], 2.1957, 1.6577)
    s6=Surf(xx[6], 8.189, 1.0)
    s7=Surf(0.0, 85, 1.0 )
    OptSys1 = (s,s1,s2,s3,s4,s5,s6,s7)
    Propagate(pencil, OptSys1)  # does work !!!
    rr2 = (Pencil_rms(pencil))[1]
    return rr1 + rr2
end
#result = optimize(FF,xx,NelderMead())
    
        result = optimize(FF,xx)
        
        #xx = copy(xx1)
    
    for ii in 1:10
        xx1=Optim.minimizer(result)
        global result = Optim.optimize(FF,xx1)
    end
    #=
        xx1=Optim.minimizer(result)
    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)


    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)
    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)

    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)
    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)

    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)
    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)

    xx1=Optim.minimizer(result)
    result = Optim.optimize(FF,xx1)
    xx1=Optim.minimizer(result)
    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)


    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)
    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)

    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)
    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)

    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)
    xx1=Optim.minimizer(result)
    result = optimize(FF,xx1)
    xx1=Optim.minimizer(result)

    =#
    xx1=Optim.minimizer(result)
    println(xx1, Optim.minimum(result))
pencil = make__hex_pencil(3000.0, 0.0, entr_pupil, 11)
    s=Surf(xx1[1], 20000.0, 1.62166)
    s1=Surf(0.0, 6.1878, 1.0)
    s2=Surf(xx1[2], 6.9862, 1.579)
    s3=Surf(xx1[3], 2.495, 1.0)
    s4=Surf(xx1[4], 7.4852, 1.5791)
    s5=Surf(xx1[5], 2.1957, 1.6577)
    s6=Surf(xx1[6], 8.189, 1.0)
    s7=Surf(0.0, 85, 1.0 )
    OptSys1 = (s,s1,s2,s3,s4,s5,s6,s7)
    Propagate(pencil, OptSys1)  # does work !!!
    # extracting the ray coordinates
    xx = Array{Float64}(undef,size(pencil))
    yy = Array{Float64}(undef,size(pencil))
        for i in 1:size(pencil)[1]
            #println(pencil[i,j])
            rr = pencil[i]
             xx[i] = rr.y
             yy[i] = rr.z
        end
    #println(Pencil_rms(pencil))
    (display(scatter(xx, yy)))

#println(FF(xx1))
#=
pencil = make__hex_pencil(0.0, 0.0, entr_pupil, 41)

ff = Array{Float64}(undef, 100)
dd = Array{Float64}(undef, 100)
for ii in 1:100
    dist = 83+0.01*ii
    pencil = make__hex_pencil(4000.0, 0.0, entr_pupil, 11)
    s7=Surf(0.0, dist, 1.0 )
    OptSys1 = (s,s1,s2,s3,s4,s5,s6,s7)
    Propagate(pencil, OptSys1)  # does work !!!
    ff[ii] = (Pencil_rms(pencil))[1]
    dd[ii] = dist
    ii = ii + 1
end
(display(plot(dd,ff)))
=#
#=
l=Surf(0, 20000.0, 1.5)
l1=Surf(-50, 20.0, 1.0)
l2=Surf(-20, 100., 1.0)
OptSys=(l, l1, l2)
#s=Surf(0.0, 20000.0, 1.0)
#OptSys = (s,)

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
    @time Propagate(pencil, OptSys1)  # does work !!!
    
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
    
=#
end # module
