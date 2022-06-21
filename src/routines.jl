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
    R1 = Ray(R.x+1*(R.y^2+R.z^2)/2.0/S.f , R.y, R.z, cx, cy, cz, R.n)
    #println(acos(cx))
    #println(ccz)
    return R1
end

function Propagate(R::Ray, S::Surf)
    
    # original Donald P Feder JOSA 41 p630 Sep 1951
    # Optical calculations with automatic computer machinery
    #The function propagates ray distance d to the surface
    # and then calculates the coordinates of refracted ray in the
    # intersection point
    #Rcx=sqrt(1.0 - R.cz^2 - R.cy^2)# calculated assumnig unit ray vector
    ss = sqrt(R.cx^2 + R.cy^2 + R.cz^2)
    Rcx=R.cx/ss
    Rcy=R.cy/ss
    Rcz=R.cz/ss

    mu=R.n/S.n  # ratio of the refraction indices
    if S.r == 0.0
        c = 0.0   # curvature of the surface
    else
        c = 1.0 / S.r
    end

    #println(c*(R.x^2+R.y^2+R.z^2)-2.0*R.x)  # test on distances
    #println(R.cx, " ", R.cy, " ", R.cz)  # test on angles
    
    e= S.d*Rcx - R.x*Rcx -R.y*Rcy - R.z*Rcz
    M1x = R.x + e*Rcx - S.d
    M12 = R.x^2 + R.y^2 + R.z^2 - e^2 + S.d^2 - 2.0 * S.d * R.x
    #println(R.x,"  ",R.y,"  ",R.z,"  ",Rcx,"  ",Rcy,"  ",Rcz,"  ")
    #=
    try
        ksi1 = sqrt(Rcx^2 - c*(c*M12 -2.0 * M1x))
    catch e
        println(R.x,"  ",R.y,"  ",R.z,"  ",Rcx,"  ",Rcy,"  ",Rcz,"  ")
    end
    =#
    ksi1 = sqrt(Rcx^2 - c*(c*M12 -2.0 * M1x))
    #println(ksi1)
    L = e + (c*M12 - 2*M1x)/(Rcx + ksi1)
    x = R.x + L*Rcx - S.d
    y = R.y + L*Rcy
    z = R.z + L*Rcz

    ksi11 = sqrt(1.0 - mu^2 * (1.0 - ksi1^2))
    g1 = ksi11 - mu*ksi1
    cx = mu * Rcx - g1 * c * x + g1
    cy = mu * Rcy - g1 * c * y
    cz = mu * Rcz - g1 * c * z
    ll = sqrt(cx*cx + cy*cy + cz*cz)
    #if ll > 1 println(ll) end
    R1 = Ray(x,y,z,cx/ll,cy/ll,cz/ll,S.n)
    #println(g1)
    #println(c)
    #println(y)
    return R1
    

   #println(c*(R1.x^2+R1.y^2+R1.z^2)-2*R1.x)  # test on distances
   #println(R1.cx^2+R1.cy^2+R1.cz^2)  # test on angles
end


""" 
Forms hex grid with origin at 0.0 and horizontal step of 1.0
The Hex Order is defined by number of elements contained in one of hexagon sides   
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

function make__hex_pencil(y0, z0, pupil::Pupil, Hex_Order=5)
    dd = sqrt((pupil.y-y0)^2+(pupil.z-z0)^2+pupil.distance^2)
    dir_y=(pupil.y - y0) / dd
    dir_z=(pupil.z - z0) / dd
    maxang = (pupil.radius) / dd
    angles = HexagonalGrid(Hex_Order)
    np=size(angles)[1] 
    #println(np)
    pencil = Vector{Ray}(undef,np)
    angles = angles * maxang / (2.0 * (Hex_Order-1.))
    for i in 1:np
        angles[i,1] = angles[i,1]+dir_y
        angles[i,2] = angles[i,2]+dir_z
            cx=sqrt(1.0 - (angles[i,1])^2 - (angles[i,2])^2)
            pencil[i] = Ray(0.0, y0, z0, cx, angles[i,1] , angles[i,2] , 1.0)
    end
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


""" propagation of the ray pencil through тхе system """
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
        for i in 1:numpoints
            GC[i,1] = Ray_pencil[i].y
            GC[i,2] = Ray_pencil[i].z
        end
        GC[:,1] = GC[:,1] .- mean(GC[:,1])
        GC[:,2] = GC[:,2] .- mean(GC[:,2])
        GC[:,3] =sqrt.(GC[:,1].^2 + GC[:,2].^2)
        #return("RMS y=",std(GC[:,1]),"RMS z=",std(GC[:,2]),"RMS radial=",std(GC[:,3]))
        return(std(GC[:,3]),std(GC[:,1]),std(GC[:,2]))  
    end  
    
    
""" Refractive index from the glass nema, using Abbe approximation """
    function n_Abbe(lambda, glass::Glass)
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
        cc1 = ((a/lambda) +b )
        return cc1 
    
    end
    