# Generating structures
# Bus_Data 13 items
mutable struct Data_Bus
    busID
    bustype
    Pd
    Qd
    Gs
    Bs
    area
    Vm
    Va
    baseKV
    zone
    Vmax
    Vmin
end
# Gen_data 21 items
mutable struct Data_Gen
    bus
    Pg
    Qg
    Qmax
    Qmin
    Vg
    mBase
    status
    Pmax
    Pmin
    Pc1
    Pc2
    Qc1min
    Qc1max
    Qc2min
    Qc2max
    ramp_agc
    ramp_10
    ramp_30
    ramp_q
    apf
end
# Gen_Cost 7 items
mutable struct Data_GenCost
    type
    startup
    shutdown
    n
    ck2
    ck1
    ck0
end
# Branch_Data 13 items
mutable struct Data_Branch
    fbus
    tbus
    r
    x
    b_s
    rateA
    rateB
    rateC
    ratio
    angle
    status
    angmin
    angmax
    branchnum
    a
    g
    b
end
# Area_Data 2 items
mutable struct Data_Area
    area
    refbus
end
# new most used values
mutable struct Data_System
    nbus
    nbranch
    Branchfrom
    Branchto
    Busnumber
    Sbase
    ngen
    Gen_Bus
    neig_bus
    in_lines
    out_lines
    GenBusnumber
    in_lines_bus
    out_lines_bus
    nareas
    Res_per_area
end

system = Data_System(size(Cases[system_name]["bus"])[1],
                size(Cases[system_name]["branch"])[1],
                zeros(Float64,1,size(Cases[system_name]["branch"])[1]),
                zeros(Float64,1,size(Cases[system_name]["branch"])[1]),
                collect(1:size(Cases[system_name]["bus"])[1]),
                Cases[system_name]["baseMVA"],
                size(Cases[system_name]["gen"])[1],
                [Int[] for i in 1:size(Cases[system_name]["bus"])[1]],
                [Int[] for i in 1:size(Cases[system_name]["bus"])[1]],
                [Int[] for i in 1:size(Cases[system_name]["bus"])[1]],
                [Int[] for i in 1:size(Cases[system_name]["bus"])[1]],
                zeros(1,size(Cases[system_name]["gen"])[1]),
                [Int[] for i in 1:size(Cases[system_name]["bus"])[1]],
                [Int[] for i in 1:size(Cases[system_name]["bus"])[1]],
                0.00,
                [[] for i in 1:maximum(unique(Cases[system_name]["bus"][:,7]))])

Gen = Data_Gen(Int.(Cases[system_name]["gen"][:,1]),
                Cases[system_name]["gen"][:,2],
                Cases[system_name]["gen"][:,3],
                Cases[system_name]["gen"][:,4],
                Cases[system_name]["gen"][:,5],
                Cases[system_name]["gen"][:,6],
                Cases[system_name]["gen"][:,7],
                Cases[system_name]["gen"][:,8],
                Cases[system_name]["gen"][:,9],
                Cases[system_name]["gen"][:,10],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen],
                [[] for i in 1:system.ngen])
if Cases[system_name]["version"] == '2'
        push!(Gen.Pc1,Cases[system_name]["gen"][:,11])
        push!(Gen.Pc2,Cases[system_name]["gen"][:,12])
        push!(Gen.Qc1min,Cases[system_name]["gen"][:,13])
        push!(Gen.Qc1max,Cases[system_name]["gen"][:,14])
        push!(Gen.Qc2min,Cases[system_name]["gen"][:,15])
        push!(Gen.Qc2max,Cases[system_name]["gen"][:,16])
        push!(Gen.ramp_agc,Cases[system_name]["gen"][:,17])
        push!(Gen.ramp_10,Cases[system_name]["gen"][:,18])
        push!(Gen.ramp_30,Cases[system_name]["gen"][:,19])
        push!(Gen.ramp_q,Cases[system_name]["gen"][:,20])
        push!(Gen.apf,Cases[system_name]["gen"][:,21])
    # end
end

Bus = Data_Bus(Cases[system_name]["bus"][:,1],
             Cases[system_name]["bus"][:,2],
             Cases[system_name]["bus"][:,3],
             Cases[system_name]["bus"][:,4],
             Cases[system_name]["bus"][:,5],
             Cases[system_name]["bus"][:,6],
             Cases[system_name]["bus"][:,7],
             Cases[system_name]["bus"][:,8],
             Cases[system_name]["bus"][:,9],
             Cases[system_name]["bus"][:,10],
             Cases[system_name]["bus"][:,11],
             Cases[system_name]["bus"][:,12],
             Cases[system_name]["bus"][:,13])

Branch = Data_Branch(Cases[system_name]["branch"][:,1],
                Cases[system_name]["branch"][:,2],
                Cases[system_name]["branch"][:,3],
                Cases[system_name]["branch"][:,4],
                Cases[system_name]["branch"][:,5],
                Cases[system_name]["branch"][:,6],
                Cases[system_name]["branch"][:,7],
                Cases[system_name]["branch"][:,8],
                Cases[system_name]["branch"][:,9],
                Cases[system_name]["branch"][:,10],
                Cases[system_name]["branch"][:,11],
                Cases[system_name]["branch"][:,12],
                Cases[system_name]["branch"][:,13],
                collect(1:system.nbranch),
                ones(1,system.nbranch),
                zeros(1,system.nbranch),
                zeros(1,system.nbranch))

Gen_cost = Data_GenCost(Cases[system_name]["gencost"][:,1],
                  Cases[system_name]["gencost"][:,2],
                  Cases[system_name]["gencost"][:,3],
                  Cases[system_name]["gencost"][:,4],
                  Cases[system_name]["gencost"][:,5],
                  Cases[system_name]["gencost"][:,6],
                  Cases[system_name]["gencost"][:,7])
# Data Change
for i = 1:system.nbranch
  a = findall(isequal(Branch.fbus[i]),Bus.busID)
  system.Branchfrom[i] = a[1]
  b = findall(isequal(Branch.tbus[i]),Bus.busID)
  system.Branchto[i] = b[1]
  if Branch.ratio[i]!= 0
      Branch.a[i] = 1/Branch.ratio[i]
  end
end

Bus.Pd = Bus.Pd / system.Sbase
Bus.Qd = Bus.Qd / system.Sbase
Bus.Bs = Bus.Bs / system.Sbase
Bus.Gs = Bus.Gs / system.Sbase

Branch.b_s     = Branch.b_s / 2
Branch.angle = Branch.angle*pi/180
Branch.rateA = Branch.rateA / system.Sbase
Branch.rateB = Branch.rateB / system.Sbase
Branch.rateC = Branch.rateC / system.Sbase
Branch.angmin= Branch.angmin*pi/180
Branch.angmax= Branch.angmax*pi/180
Branch.g     = real(1 ./(Branch.r+Branch.x*im))
Branch.b     = imag(1 ./(Branch.r+Branch.x*im))

Gen.Qmax    =Gen.Qmax / system.Sbase
Gen.Qmin    =Gen.Qmin / system.Sbase
Gen.Pmax    =Gen.Pmax / system.Sbase
Gen.Pmin    =Gen.Pmin / system.Sbase
Gen.Pg      =Gen.Pg / system.Sbase
Gen.Qg      =Gen.Qg / system.Sbase

Gen_cost.ck1 = Gen_cost.ck1*system.Sbase
Gen_cost.ck2 = Gen_cost.ck2*system.Sbase^2

for i in 1:system.nbranch
    a = convert(Int, system.Branchfrom[i])
    b = convert(Int, system.Branchto[i])
    push!(system.neig_bus[a],b)
    push!(system.neig_bus[b],a)
end
for i in 1:system.ngen
    a = findall(x->x==Gen.bus[i],Bus.busID)[1]
    system.GenBusnumber[i]=a[1]
end
for i in 1:system.ngen
    a = convert(Int, system.GenBusnumber[i])
    push!(system.Gen_Bus[a],i)
end

#Identificacion de lineas que salen y que llegan a una determinada barra
for i in 1:system.nbranch
	a = convert(Int, system.Branchfrom[i])
	b = convert(Int, system.Branchto[i])
	    push!(system.out_lines[a],i)
	    push!(system.in_lines[b],i)
end
system.nareas = length(unique(Bus.area))
