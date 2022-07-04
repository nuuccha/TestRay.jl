
""" Primitive system plot, requires Luxor.jl """
function Plot_Sys_2D(s::Array{Surface})
    (length, left, right) = SysThickness(s)
    Drawing(640, 400, "Lens.png")
    scal = 640.0 / (right -left) /1.3
    #scal = 1.
    origin(100,200)
    #bb=BoundingBox()
    t=0.0
    background("white")
    sethue("tomato")
    for i in 1:size(s)[1]
        #println(s[i].d," ",s[i].r," ",s[i].aper)
        t=t + s[i].d 
        tt = scal*t
        if isdefined(s[i],:r)
            r = scal * s[i].r
        else
            if s[i].curv == 0
                r=0
                else
                r = scal / s[i].curv
            end 
        end   
        ap=scal * s[i].aper
        
        if r == 0
            line(Point(tt, -ap), Point(tt,ap), :stroke)
        end
        if r > 0
            arc2r(Point(tt+r, 0), Point(tt + ap^2/2.0/r,ap), Point(tt+ap^2/2.0/r, -ap), :stroke)
        end
        if r < 0
            arc2r(Point(tt+r, 0), Point(tt + ap^2/2.0/r,-ap), Point(tt+ap^2/2.0/r, ap), :stroke)
        end

        
    end
    pencil = make_flat_pencil(1, entr_pupil,5)
    t=0
    for i in 1:4
        t1 = t
        t = t + s[i].d 
        pencil_ = copy(pencil)
        
        #Propagate(pencil,(s[i,1],))
        for j in 1:4
            pencil[j] = Propagate(pencil[j],s[j,1])
            println(pencil[j])
            line(Point(scal *t1, scal *pencil_[j].y), Point(scal *t,scal *pencil[j].y), :stroke)
        end
        

    end
    finish()
    preview()
end
