

mutable struct s1
    aaa::Any
    bbb::Any
    ccc::Any
end

s=s1(2,3,4);
fd=s1(5,6,8)

"""
This macro returns the string containing the actual name of the structure
"""
macro return_my_name(s)
    String(s)
end
""" this macro modifies the num-th field of a mutable structure to value val""" 
macro mutate_struct_field(s,num,val) 
   
    quote
        val = $(esc(val))
        num = $(esc(num))
        cc = String(@return_my_name($s)*"."* String((fieldnames(typeof($s))[num])))
        dd = eval(Meta.parse(cc*"="*string(val)))
    end
end

for ii in 1:3
    #cc = ii
    #cc1 = ii^2
    println(ii)
    @time @mutate_struct_field(fd,ii,sin(ii))
    #@mutate_struct_field(s,ii,ii)
    println(ii,"  ",fd)
    #println(ii,"  ",s)
end

#=

macro fieldN(s,num, val)
    quote
    a =  (fieldnames(typeof(s))[$num])
    b = Symbol("s.",a)
    eval($b = $val)
    end
end
=#
#@fieldN(s,2)
#println(b)
#s
#=
a =  fieldnames(typeof(s))[num]
println(a)
b = Symbol("s.",a)
println(b)
b=10
println(getfield(s,2))

=#
#=
macro structd(s,num,val)
    quote
    a =  (fieldnames(typeof(s))[$num])
    b = Symbol("s.",a)
    eval($b = $val)
    end
end

println(@structd(s,3,10))
=#

#println(ddd)
#=

macro sayhello(name)
    return :( println("Hello, ", $name) )
end
#println(String(@structd(s,2)))
#cc =eval(@structd(s,2))




struct s2
    a::Any
    b::Any
    c::Any
    d::Any
end
function len(a)
    l=0
    for c in (fieldnames(typeof(a)))
        l += 1
    end
    return l
end

function gcopy(a::s2)
    l=0
    for c in (fieldnames(typeof(a)))
        l += 1
        bb[l]=getfield(a,l)
    end
    return bb
end

a = s1("a1",pi,1)
b = s2("b1",1,28,pi^2)
#ccc = gcopy(b)
#println(ccc)
println(len(a))
println(len(b))

println(a)
println(b)
for c in fieldnames(typeof(b))
    println(c)
end

ff = ([1,2.0],[3,4])
ff[1][2] = 5
ff
ff1 = collect(ff)
ff1
using Printf
macro timeit(ex)
    quote
        local t0 = time()
        local val = $(esc(ex))
        local t1 = time()
        println("elapsed time in seconds: ")
        @printf "%.3f" t1-t0
        #@printf (val)
        val
    end
end
@timeit (for i in 1:10000 println(sqrt(11.)) end)
macroexpand(@timeit, 10.)
=#