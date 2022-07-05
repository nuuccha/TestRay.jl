#module TestRay
using Plots
using Statistics
using Optim

include("structures.jl")

include("routines.jl")
lambda = [450, 550, 630]
field = [0.0 0.0; 5000.0 0.0; 0.0 5000.0]
# array for all surfaces and all lambdas
s = Array{Surface,2}(undef, 8, size(lambda)[1])
ss = Array{Surface,2}(undef, 8, size(lambda)[1])

# definition of all surfaces with glasses
for ilam in 1:size(lambda)[1]
    ss[1,ilam]=Abbe_Surf(33.97, 20000.0, F2, lambda[ilam], 10.0, 0., 0)
    ss[2,ilam]=Abbe_Surf(0.0, 6.1878, AIR, lambda[ilam], 10.0, 0., 0)
    ss[3,ilam]=Abbe_Surf(-62.67, 6.9862, BK7, lambda[ilam], 10.0, 0. ,0)
    ss[4,ilam]=Abbe_Surf(29.58, 2.495, AIR, lambda[ilam], 10.0, 0., 0)
    ss[5,ilam]=Abbe_Surf(1625.0, 7.4852, F5, lambda[ilam], 10.0, 0., 0)
    ss[6,ilam]=Abbe_Surf(27.34, 2.1957, LAFN7, lambda[ilam], 10.0, 0., 0)
    ss[7,ilam]=Abbe_Surf(-46.366, 8.189, AIR, lambda[ilam], 10.0, 0., 0)
    ss[8,ilam]=Abbe_Surf(0.0, 83.732, AIR , lambda[ilam], 10.0, 0., 0)
end
# conversion of all surfaces into Surf
for ilam in 1:size(lambda)[1]
    for i = 1:size(ss)[1]
        s[i,ilam] =ConvFromAbbe(ss[i,ilam])
    end 
end

entr_pupil = Pupil(20000.0, 10.0 ,0.0 ,0.0)
# array of variables, here only 6, but can be any size
xx = Array{Float64}(undef, 6)
xx[1] = 1. / ss[1,1].r
xx[2] =1. / ss[3,1].r
xx[3] =1. / ss[4,1].r
xx[4] = 1. / ss[5,1].r
xx[5] = 1. / ss[6,1].r
xx[6] = 1. / ss[7,1].r
# bounding arrays in terms of curvature
u = copy(xx)
l = copy(xx)
for i in 1:size(xx)[1]
    u[i] = 0.03
    l[i] = -0.03
end


function FF(xx::Array{Float64})
    rr=0.
    for ifield in 1:3
        for ilam in 1:3
            #println(field[ifield,1])
            pencil = make__hex_pencil(field[ifield,1], field[ifield,2],entr_pupil, 5)
            
                #println(field[ifield,1])
            s[1,ilam].curv = xx[1]
            s[3,ilam].curv = xx[2]
            s[4,ilam].curv = xx[3]
            s[5,ilam].curv = xx[4]
            s[6,ilam].curv = xx[5]
            s[7,ilam].curv = xx[6]
            
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
#gr()

for ifield  in 1:3
for ilam in  1:3
    #println(field[ifield,1])
    pencil = make__hex_pencil(field[ifield,1], field[ifield,2],entr_pupil, 31)
    
        #println(field[ifield,1])
        #s[ii,ilam].curv = xx1[ii]
        s[1,ilam].curv = xx1[1]
        s[3,ilam].curv = xx1[2]
        s[4,ilam].curv = xx1[3]
        s[5,ilam].curv = xx1[4]
        s[6,ilam].curv = xx1[5]
        s[7,ilam].curv = xx1[6]

    Propagate(pencil, s[ : , ilam])  # does work !!!
    #show_pencil(pencil)
    println(Pencil_rms(pencil)[1])
    #sleep(1.3)
end
end
Plot_Sys_2D(s[1:size(s)[1],1])

#=
        result = optimize(FF,xx, l, u, NelderMead())
        
        #xx = copy(xx1)
    
    for ii in 1:10
        xx2=Optim.minimizer(result)
        #xx2 = xx2 .+ rand(Float64,6)
        global result = Optim.optimize(FF,xx2, l, u, NelderMead())
        
    end
    
    xx1=Optim.minimizer(result)
    println(1.0 ./ xx1, "  ", Optim.minimum(result))
    gr()
    for ifield  in 1:3
        for ilam in  1:3
            #println(field[ifield,1])
            pencil = make__hex_pencil(field[ifield,1], field[ifield,2],entr_pupil, 11)
            for ii in 1:size(xx1)[1]
                #println(field[ifield,1])
                s[ii,ilam].curv = xx1[ii]
            end
            Propagate(pencil, s[ : , ilam])  # does work !!!
            show_pencil(pencil)
            println(Pencil_rms(pencil)[1])
            sleep(1)
        end
    end
 
    =#

#=
    pencil = make__hex_pencil(0., 0.0, entr_pupil, 11)
    #s=Surf(xx1[1], 20000.0, 1.62166, 0., 0)
    sos = Abbe_Surf(xx1[1], 20000.0, NF2, 533, 0., 0)
    s = ConvFromAbbe(sos) 
    s1=Surf(0.0, 6.1878, 1.0, 0., 0)
    s2=Surf(xx1[2], 6.9862, 1.579, 0., 0)
    s3=Surf(xx1[3], 2.495, 1.0, 0., 0)
    s4=Surf(xx1[4], 7.4852, 1.5791, 0., 0)
    s5=Surf(xx1[5], 2.1957, 1.6577, 0., 0)
    s6=Surf(xx1[6], 8.189, 1.0, 0., 0)
    s7=Surf(0.0, 85, 1.0 , 0., 0)
    OptSys1 = (s,s1,s2,s3,s4,s5,s6,s7)
    Propagate(pencil, OptSys1)  # does work !!!
    println(Pencil_rms(pencil)[1])
    show_pencil(pencil)
    =#
    #=
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
=#

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
#end # module
