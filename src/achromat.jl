#module TestRay
using Plots
using Statistics
using Optim

include("structures.jl")

include("routines.jl")
lambda = [450, 550, 630]
field = [0.0 0.0; 1000.0 0.0; 0.0 1000.0]
# array for all surfaces and all lambdas
s = Array{Surface,2}(undef, 5, size(lambda)[1])
ss = Array{Surface,2}(undef, 5, size(lambda)[1])

# definition of all surfaces with glasses
for ilam in 1:3
    ss[1,ilam]=Abbe_Surf(1000.0, 20000.0, K10, lambda[ilam], 0., 0)
    ss[2,ilam]=Abbe_Surf(1000.0, 5, SF6, lambda[ilam], 0., 0)
    ss[3,ilam]=Abbe_Surf(1000, 5. , F2, lambda[ilam], 0. ,0)
    ss[4,ilam]=Abbe_Surf(1000, 5. , AIR, lambda[ilam], 0. ,0)
    ss[5,ilam]=Abbe_Surf(1000, 100. , AIR, lambda[ilam], 0. ,0)
end
# conversion of all surfaces into Surf
for ilam in 1:size(lambda)[1]
    for i = 1:size(ss)[1]
        s[i,ilam] =ConvFromAbbe(ss[i,ilam])
    end 
end

entr_pupil = Pupil(20000.0, 5.0 ,0.0 ,0.0)
# array of variables, here only 6, but can be any size
xx = Array{Float64}(undef, 4)
xx[1] = 1. / 2000
xx[2] = 1. / 2000
xx[3] = 1. / 1000
xx[4] = 1. / 3000
# bounding arrays in terms of curvature
u = copy(xx)
l = copy(xx)
for i in 1:4
    u[i] = 0.0001
    l[i] = -0.0001
end

function FF(xx::Array{Float64})
    rr=0.
    for ifield in 1:3
        for ilam in 1:3
            #println(field[ifield,1])
            pencil = make__hex_pencil(field[ifield,1], field[ifield,2],entr_pupil, 11)
            s[1,ilam].r = 1 / xx[1]
            s[2,ilam].r = 1 / xx[2]
            s[3,ilam].r = 1 / xx[3]
            s[4,ilam].r = 1 / xx[4]
            Propagate(pencil, s[ : , ilam])  # does work !!!
            rr = rr + (Pencil_rms(pencil))[1]
        end
    end
    return rr
end
    
        result = optimize(FF,xx, l, u, NelderMead(),Optim.Options(g_tol = 1e-12))
        
        #xx = copy(xx1)
    
    for ii in 1:1
        xx2=Optim.minimizer(result)
        #xx2 = xx2 .+ rand(Float64,6)
        global result = Optim.optimize(FF,xx2, l, u, NelderMead(),Optim.Options(g_tol = 1e-12))
        
    end
    
    xx1=Optim.minimizer(result)
    println(1 ./ xx1, "  ", Optim.minimum(result))
    gr()
    for ifield  in 1:3
        for ilam in  1:3
            #println(field[ifield,1])
            pencil = make__hex_pencil(field[ifield,1], field[ifield,2],entr_pupil, 11)
            s[1,ilam].r = 1 ./ xx1[1]
            s[2,ilam].r = 1 ./ xx1[2]
            s[3,ilam].r = 1 ./ xx1[3]
            s[4,ilam].r = 1 ./ xx1[4]

            Propagate(pencil, s[ : , ilam])  # does work !!!
            show_pencil(pencil)
            println(Pencil_rms(pencil)[1])
            sleep(1.3)
        end
    end
 