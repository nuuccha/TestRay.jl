function Propagate(R::Ray, S::Paraxial)
    # propagation through paraxial element
    #=
    hy = R.cy / sqrt(1.0 - R.cy*R.cy)
    hz = R.cz / sqrt(1.0 - R.cz*R.cz)
    hy = hy * S.f
    hz = hz * S.f

    ccy= -(R.y - hy) / S.f
    ccz= -(R.z - hz) / S.f

    cy=ccy / sqrt(1.0 + ccy*ccy) 
    cz=ccz / sqrt(1.0 + ccz*ccz)
    =#
    #cx=sqrt(1-R.cy^2 - R.cz^2)  
    hy = ((R.cy)) * S.f
    hz = ((R.cz)) * S.f 
    cy= ((-(R.y - hy) / S.f))
    cz= ((-(R.z - hz) / S.f))
    cx=sqrt(1.0-cy^2 - cz^2)
    R1 = Ray(R.x+1*(R.y^2+R.z^2)/2.0/S.f , R.y, R.z, cx, cy, cz, R.n, 0)
    if (R.y^2+R.z^2) > S.aper^2 R1.empty = 1 end
    #println(acos(cx))
    #println(ccz)
    return R1
end

""" original Donald P Feder JOSA 41 p630 Sep 1951
    Optical calculations with automatic computer machinery
    The function propagates ray distance d to the surface
    and then calculates the coordinates of refracted ray in the
    intersection point
    Rcx=sqrt(1.0 - R.cz^2 - R.cy^2)# calculated assumnig unit ray vector """
function Propagate(R::Ray, S::Surf)
    
    if R.empty != 0
        return Ray(R.x,R.y,R.z,R.cx,R.cy,R.cz,S.n, 1) 
    end
    ss = sqrt(R.cx^2 + R.cy^2 + R.cz^2)
    Rcx=R.cx/ss
    Rcy=R.cy/ss
    Rcz=R.cz/ss

    mu=R.n/S.n  # ratio of the refraction indices
    c=S.curv    # curvagure of the surface
        #println(c*(R.x^2+R.y^2+R.z^2)-2.0*R.x)  # test on distances
        #println(R.cx, " ", R.cy, " ", R.cz)  # test on angles
        
    e= S.d*Rcx - R.x*Rcx -R.y*Rcy - R.z*Rcz
    M1x = R.x + e*Rcx - S.d
    M12 = R.x^2 + R.y^2 + R.z^2 - e^2 + S.d^2 - 2.0 * S.d * R.x
        #println(R.x,"  ",R.y,"  ",R.z,"  ",Rcx,"  ",Rcy,"  ",Rcz,"  ")
        
    if (Rcx^2 - c*(c*M12 -2.0 * M1x)) < 0.0
        R1 = R1 = Ray(R.x,R.y,R.z,R.cx,R.cy,R.cz,S.n, 1)
        return R1
    end
        #println("!!!!!!!!!!!!!!!!!!!!!!!!!")
    ksi1 = sqrt(Rcx^2 - c*(c*M12 -2.0 * M1x))
        #println(ksi1)
    L = e + (c*M12 - 2*M1x)/(Rcx + ksi1)
    x = R.x + L*Rcx - S.d
    y = R.y + L*Rcy
    z = R.z + L*Rcz
    if  (1.0 - mu^2 * (1.0 - ksi1^2)) < 0.0
        R1 = Ray(R.x,R.y,R.z,R.cx,R.cy,R.cz,S.n, 1)
        return R1
    end
    ksi11 = sqrt(1.0 - mu^2 * (1.0 - ksi1^2))
    g1 = ksi11 - mu*ksi1
    cx = mu * Rcx - g1 * c * x + g1
    cy = mu * Rcy - g1 * c * y
    cz = mu * Rcz - g1 * c * z
    ll = sqrt(cx*cx + cy*cy + cz*cz)
        #if ll > 1 println(ll) end
    if (y^2+z^2) > S.aper^2 # cheking the vignetting 
        R1 = Ray(x,y,z,cx/ll,cy/ll,cz/ll,S.n, 1)
        else
        R1 = Ray(x,y,z,cx/ll,cy/ll,cz/ll,S.n, 0)
    end 
    return R1
   #println(c*(R1.x^2+R1.y^2+R1.z^2)-2*R1.x)  # test on distances
   #println(R1.cx^2+R1.cy^2+R1.cz^2)  # test on angles
end


""" 
Forms hex grid with origin at 0. 0. and horizontal step of 1.0
The Hex Order is defined by number of elements contained in one of the hexagon sides
For HexOrder = 5 the total number of hexagons in the grid will be 61   
"""
function HexagonalGrid(HexOrder=5, StepX = 1.0, xc=0, yc=0)
    numpoints = Int(round(HexOrder^2*3.1))
    GridCoord = Array{Float64, 2}(undef, numpoints ,2)
    StepY=StepX*sqrt(3.0)/2.0
    RowMax = HexOrder*2-1
    ij=0  # counter of generated points
    for i in 1:HexOrder # Double loop for the top part of the hexagon
        for j in append!(collect(0:(HexOrder-i)),collect((1-HexOrder):-1))
            ij += 1
            GridCoord[ij,1] = (j)*StepX +(i-1)*StepX/2. +xc
            GridCoord[ij,2] = (i-1)*StepY +yc
            #println(GridCoord[ij,1])
            #println(GridCoord[ij,2])
        end
    end
    
    for i in 2:HexOrder # Double loop for the bottom part of the hexagon
        for j in append!(collect(0:(HexOrder-i)),collect((1-HexOrder):-1))
            ij += 1
            GridCoord[ij,1] = (j)*StepX + (i-1)*StepX/2. +xc
            GridCoord[ij,2] = (i-1)*(-StepY) +yc
            #println(GridCoord[ij,1])
            #println(GridCoord[ij,2])
        end
    end
    #GridCoord = GridCoord ./ sqrt(GridCoord[ij,1]^2 + GridCoord[ij,2]^2)
    #println(GridCoord)
    return GridCoord[1:ij,:]
end

""" Hexagonal ray pencil pointed from x0 y0 to the entrance pupil 
    Hex_order is the number of rays in the side of hexagon
    for Hex_Order = 5 a hexagon pencil with 61 ray will be produced"""
function make__hex_pencil(y0, z0, pupil::Pupil, Hex_Order=5)
    # directions to the pupil center:
    dd = sqrt((pupil.y-y0)^2+(pupil.z-z0)^2+pupil.distance^2)
    dir_y=(pupil.y - y0) / dd
    dir_z=(pupil.z - z0) / dd
    # the pencil angle
    maxang = (pupil.radius) / dd
    angles = HexagonalGrid(Hex_Order)
    radii = HexagonalGrid(Hex_Order)
    np=size(radii)[1] 
    #println(np)
    pencil = Vector{Ray}(undef,np)
    angles = angles * maxang / (2.0 * (Hex_Order-1.))
    radii = radii * pupil.radius / (2.0 * (Hex_Order-1.))
    for i in 1:np
        radii[i,1] = radii[i,1]+pupil.y
        radii[i,2] = radii[i,2]+pupil.z
        dd = sqrt((radii[i,1]-y0)^2+(radii[i,2]-z0)^2+pupil.distance^2)
        angles[i,1]=(radii[i,1] - y0) / dd
        angles[i,2]=(radii[i,2] - z0) / dd
        cx=sqrt(1.0 - (angles[i,1])^2 - (angles[i,2])^2)
        pencil[i] = Ray(0.0, y0, z0, cx, angles[i,1] , angles[i,2] , 1.0, 0)
    end
    return(pencil)
end

""" Flat ray pencil pointed from y0 to the entrance pupil 
    Hex_order is the number of rays in the side of hexagon
    for Hex_Order = 5 a hexagon pencil with 61 ray will be produced"""
function make_flat_pencil(y0, pupil::Pupil, N_Rays = 5)
    pencil = Vector{Ray}(undef,N_Rays)
    angles = Array{Float64, 2}(undef, N_Rays ,2)
    dp = 2.0 * pupil.radius /(N_Rays -1.) # step inside pupil
    
    for i in 1:N_Rays
        rr = pupil.y + (i-1)*dp - pupil.radius # coord inside pupil
        dd = sqrt((rr-y0)^2+pupil.distance^2) # dist to the point
        angles[i,1] = (rr-y0)/dd # angle 
        angles[i,2] = 0.0 # z angle = 0 
        cx=sqrt(1.0 - (angles[i,1])^2 - (angles[i,2])^2) # x projection
        pencil[i] = Ray(0.0, y0, 0.0, cx, angles[i,1] , angles[i,2] , 1.0, 0)
    end
    #println(angles)
    return(pencil)
end




function Propagate(Single_Ray::Ray, OptSys::Tuple)
    for OptSurf in OptSys
        Single_Ray = Propagate(Single_Ray, OptSurf)
    end
    return(Single_Ray)
end

function focus_position(pencil::Vector{Ray},OSys::Tuple)
    r = pencil[1]
    r1 = pencil[2]
    r=Propagate(r, OSys)
    r1=Propagate(r1, OSys)
    p1=(r1.y - r.y + r1.cy*r1.x - r.cy*r.x) / (r.cy-r1.cy)
    #=
    ii=0
    for cc in pencil
        ii +=1
    end
    r1= pencil[ii]
    r1=Propagate(r1, OSys)
    p2=(r1.y - r.y + r1.cy*r1.x - r.cy*r.x) / (r.cy-r1.cy)
    println(p1,"  ",p2)
    =#
    return p1
    end


""" propagation of the ray pencil through the system """
function Propagate(Ray_pencil::Vector{Ray}, OptSys::Tuple)
    #println(size(Ray_pencil)[1], " ", size(Ray_pencil)[2])
        for OptSurf in OptSys
            for i in 1:size(Ray_pencil)[1]
                Ray_pencil[i] = Propagate(Ray_pencil[i],OptSurf)
            end
        end
        return(Ray_pencil)
    end

    """ Pencil RMS """
    function Pencil_rms(Ray_pencil::Vector{Ray})
    numpoints =size(Ray_pencil)[1]
    GC = Array{Float64,2}(undef, numpoints,3)
    ij=0    
        for i in 1:numpoints
            rrr = Ray_pencil[i]
            if  rrr.empty == 0 # checking for empty ray
                ij=ij+1
                GC[ij,1] = Ray_pencil[i].y
                GC[ij,2] = Ray_pencil[i].z
            end
        end
        GC[1:ij,1] = GC[1:ij,1] .- mean(GC[1:ij,1])
        GC[1:ij,2] = GC[1:ij,2] .- mean(GC[1:ij,2])
        GC[1:ij,3] =sqrt.(GC[1:ij,1].^2 + GC[1:ij,2].^2)
        if ij < 3  # checking for ray total to be more than 3
            #println("Warning: Less than 3 rays, check vignetting")
            return (1000.,1000.,1000)
        else
        #return("RMS y=",std(GC[:,1]),"RMS z=",std(GC[:,2]),"RMS radial=",std(GC[:,3]))
            return(std(GC[1:ij,3]),std(GC[1:ij,1]),std(GC[1:ij,2]))
        end  
    end  
    
    """ propagation of the ray pencil through the system """
    function Propagate(Ray_pencil::Vector{Ray}, OptSys::Vector{Surface})
        #println(size(Ray_pencil)[1], " ", size(Ray_pencil)[2])
            for OptSurf in OptSys
                for i in 1:size(Ray_pencil)[1]
                    Ray_pencil[i] = Propagate(Ray_pencil[i],OptSurf)
                end
            end
            return(Ray_pencil)
        end

""" Refractive index from the glass data n and nu, using Abbe approximation """
    function n_Abbe(lambda, glass::AbbeGlass)
        if glass.typ != "E"
            lambda1=486.1
            lambda2=589.2
            lambda3=656.3
        else
            lambda1=480
            lambda2=546.073
            lambda3=643.8
        end
    
    # cc1 is the result of  approx in the form n = a/lambda^pp +b
        pp=2.4
        lambda1 = lambda1^pp
        lambda2 = lambda2^pp
        lambda3 = lambda3^pp
        lambda = lambda^pp
       
        nr=glass.nr
        nu=glass.nu
        a = (nr-1.0)*lambda1*lambda3/(nu*(lambda3-lambda1))
        b = nr*lambda2*nu*(lambda3-lambda1)-(nr-1.0)*lambda1*lambda3
        b=b/(lambda2*nu*(lambda3-lambda1))
        refr_ind = ((a / lambda) + b )
        return refr_ind
    
    end
    
    function show_pencil(Ray_pencil::Vector{Ray})
        xx = Array{Float64}(undef,size(Ray_pencil)[1])
        yy = Array{Float64}(undef,size(Ray_pencil)[1])
        for i in 1:size(Ray_pencil)[1]
            #println(pencil[i,j])
            
             xx[i] = Ray_pencil[i].y
             yy[i] = Ray_pencil[i].z
        end
    #println(Pencil_rms(pencil))
    (display(scatter(xx, yy)))
    #(display(scatter!(xx, yy))) # to overlay graphs
    end

function ConvFromAbbe(sur::Abbe_Surf) 
    n = n_Abbe(sur.lambda, sur.glass)
    if sur.r == 0 
        curv = 0.
    else
        curv = 1.0 / sur.r
    end
    return Surf(curv, sur.d, n, sur.e2, sur.aper, sur.stop)
end

function TupleSize(a::Tuple)
    ii=0
    for asa in a
        ii += 1
    end
    return ii
end    

""" primitive thickness of the system """
function SysThickness(s::Array{Surface})
    t= 0.0
    tmax = 0.0
    tmin= 0.0
    for surf in s
        t=t+surf.d
        if t > tmax tmax = t end
        if t < tmin tmin = t end
    end
    return(t, tmin, tmax)
end
