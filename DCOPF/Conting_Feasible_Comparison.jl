using BilevelJuMP
using JuMP
using Gurobi
using MAT
using MATLAB
using Ipopt
using Mosek
using MosekTools
using Printf
using PowerModels
using PGLib
using DataFrames
using CSV
using Plots

path_data = "C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos"

path_org = pwd()

cd(path_data)

system_name = "pglib_opf_case30_ieee"
#system_name = "pglib_opf_case57_ieee"
#system_name = "pglib_opf_case118_ieee"


Cases = Dict()
Cases[system_name] = mat"run(append($system_name,'.m'))"

cd(path_org)
# Generating structures
include("C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\DCOPF\\CreateStructures_OPF.jl")

#Branch.rateA = Branch.rateA*100 Los duales igual porq estan congestionado

# # ===================DEMANDS & DISPATCH IMPORT=====================#
path_data1 = "C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\pDC_Model_$system_name.mat"
path_data2 = "C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\DC_Model_$system_name.mat"
path_data = "C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Matpower_CPF\\Generation_Matpower_$system_name.mat"

pDC_Model = matread(path_data1)
DC_Model = matread(path_data2)
PgDC = DC_Model["Pg"]
PgPDC = pDC_Model["Pg"]

Pg_DC_feasible = zeros(1000, system.ngen)
Pg_PDC_feasible = zeros(1000, system.ngen)

Generation_Matpower = matread(path_data)
PgAC = Generation_Matpower["DataMATPOWER"]["DispatchAC_Probe"]
DEMAND = Generation_Matpower["DataMATPOWER"]["Demand_Probe"]
#PgAC = Generation_Matpower["Data1"]["DataMATPOWER"]["DispatchAC_Probe"]
#DEMAND = Generation_Matpower["Data1"]["DataMATPOWER"]["Demand_Probe"]

RANGEPD= Int(size(DEMAND,2)./2)

# # ========================IMPORT DATA============================================#
Pd = DEMAND[1:1000,1:RANGEPD]
Qd = DEMAND[1:1000,RANGEPD+1:Int(size(DEMAND,2))]
## ===============================================================================#

# Reference bus
ref = findall(x->x==3, Bus.bustype)

#Phase Shifting
cosA=cos.(Branch.angle)
sinA=sin.(Branch.angle)

# Generator index per bus
DispGen_bus = [Int[] for i in 1:system.nbus]
for i in 1:system.ngen
    a = convert(Int, system.GenBusnumber[i])
    push!(DispGen_bus[a],i)
end

# Import from PowerModels
# network_data0 = PowerModels.parse_file("pglib_opf_case30_ieee.m")

# Results: Dictionary used to initialize variables
# flag_model: 0 => convex, 1 => nonconvex, 2 => second-order conic
# flag_solver: 0 => Ipopt, 1 => Gurobi, 2 => Mosek
# noise between 0 and 1

function ACOPF(system,Branch,Bus,Gen,e_0,f_0,Results,lambda0,flag_model,flag_solver,mipgap,l,Pd,Qd,Pg_Probe,Line)

    if flag_solver == 0
        m_ac = Model(Ipopt.Optimizer)
        set_optimizer_attribute(m_ac, MOI.Silent(),true)
        set_optimizer_attribute(m_ac, "max_iter",6000)   
        set_optimizer_attribute(m_ac, "tol", 1e-8)  
    elseif flag_solver == 1
        m_ac = Model(Gurobi.Optimizer)
        set_optimizer_attribute(m_ac, "BarHomogeneous", 1) 
        # set_optimizer_attribute(m_ac, MOI.Silent()rue)
        # set_optimizer_attribute(m_ac, "OutputFlag", 0) 
        if flag_model == 1
            set_optimizer_attribute(m_ac, "NonConvex", 2) 
            set_optimizer_attribute(m_ac, "MIPGap", mipgap)
            set_optimizer_attribute(m_ac, "MIPFocus",3)
        end
    elseif flag_solver == 2
        m_ac = Model(Mosek.Optimizer)
        set_optimizer_attribute(m_ac, MOI.Silent(), true)
        set_optimizer_attribute(m_ac, "MSK_DPAR_INTPNT_CO_TOL_INFEAS", 1e-6)
    end 

    c = Dict()
    s = Dict()

    for i = 1:system.nbus
        c[i,i] = @variable(m_ac)
    end
    for k = 1:system.nbranch
        i = Int(system.Branchfrom[k])
        j = Int(system.Branchto[k])
        c[i,j] = @variable(m_ac)
        s[i,j] = @variable(m_ac)
    end

    @variable(m_ac, e[1:system.nbus])
    @variable(m_ac, f[1:system.nbus])

    # Active/real power produced by unit j at time period t
    @variable(m_ac, Pg[1:system.ngen] >= 0)
    # Reactive power produced by generating unit j at time period t
    @variable(m_ac, Qg[1:system.ngen])
    # Active Power flow from i to j nodes
    @variable(m_ac, Pfr[Branch.branchnum])
    # Reactive Power flow from i to j nodes
    @variable(m_ac, Qfr[Branch.branchnum])
    #Active Power flow from j to i nodes
    @variable(m_ac, Pto[Branch.branchnum])
    # Reactive Power flow from j to i nodes
    @variable(m_ac, Qto[Branch.branchnum])

    @variable(m_ac, sl >= 0)

    @variable(m_ac, BaseCost >= 0)  

    @constraint(m_ac, f[ref[1]] == 0)

    # lambda0 = 0.35625*4.9025e5  #39_epri

    #Phase Shifting
    cosA=cos.(Branch.angle)
    sinA=sin.(Branch.angle)

    # Initialization
    for k = 1:system.nbus
        if "e" in keys(Results)
            set_start_value(e[k], Results["e"][k])
        else
            set_start_value(e[k], 1)
        end
        if "f" in keys(Results)
            set_start_value(f[k], Results["f"][k])
        else
            set_start_value(f[k], 0)
        end
        if "c" in keys(Results)
            set_start_value(c[k,k], Results["c"][k,k])
        else
            set_start_value(c[k,k], 1)
        end
    end
    for k = 1:system.nbranch
        i = Int(system.Branchfrom[k])
        j = Int(system.Branchto[k])
        if "c" in keys(Results)
            set_start_value(c[i,j], Results["c"][i,j])
        else
            set_start_value(c[i,j], 1)
        end
        if "s" in keys(Results)
            set_start_value(s[i,j], Results["s"][i,j])
            set_start_value(s[j,i], Results["s"][j,i])
        else
            set_start_value(s[i,j], 0)
            set_start_value(s[i,j], 0)
        end
        if "Pfr" in keys(Results)
            set_start_value(Pfr[k], Results["Pfr"][k])
        end
        if "Qfr" in keys(Results)
            set_start_value(Qfr[k], Results["Qfr"][k])
        end
        if "Pto" in keys(Results)
            set_start_value(Pto[k], Results["Pto"][k])
        end
        if "Qto" in keys(Results)
            set_start_value(Qto[k], Results["Qto"][k])
        end
    end
    for k = 1:system.ngen
        if "Pg" in keys(Results)
            set_start_value(Pg[k], Results["Pg"][k])
        else
            set_start_value(Pg[k], 0.5*(Gen.Pmin[k] + Gen.Pmax[k]))
        end
        if "Qg" in keys(Results)
            set_start_value(Qg[k], Results["Qg"][k])   
        else
            set_start_value(Qg[k], 0.5*(Gen.Qmin[k] + Gen.Qmax[k]))
        end
    end

    # Total cost

    @objective(m_ac, Min,sum((Pg_Probe[l,i]-Pg[i])^2 for i in 1:system.ngen)) 


    # If model nonconvex, then...
    if flag_model == 1
        if flag_solver == 0 
            @NLconstraint(m_ac, dual1[i = 1:system.nbus], c[i,i] == e[i]^2 + f[i]^2)
            @NLconstraint(m_ac, dual2[k = 1:system.nbranch], c[Int(system.Branchfrom[k]),Int(system.Branchto[k])] == e[Int(system.Branchfrom[k])]*e[Int(system.Branchto[k])] + f[Int(system.Branchfrom[k])]*f[Int(system.Branchto[k])])
            @NLconstraint(m_ac, dual3[k = 1:system.nbranch], s[Int(system.Branchfrom[k]),Int(system.Branchto[k])] == e[Int(system.Branchfrom[k])]*f[Int(system.Branchto[k])] - e[Int(system.Branchto[k])]*f[Int(system.Branchfrom[k])])
        else
            @constraint(m_ac, dual1[i = 1:system.nbus], c[i,i] == e[i]^2 + f[i]^2)
            @constraint(m_ac, dual2[k = 1:system.nbranch], c[Int(system.Branchfrom[k]),Int(system.Branchto[k])] == e[Int(system.Branchfrom[k])]*e[Int(system.Branchto[k])] + f[Int(system.Branchfrom[k])]*f[Int(system.Branchto[k])])
            @constraint(m_ac, dual3[k = 1:system.nbranch], s[Int(system.Branchfrom[k]),Int(system.Branchto[k])] == e[Int(system.Branchfrom[k])]*f[Int(system.Branchto[k])] - e[Int(system.Branchto[k])]*f[Int(system.Branchfrom[k])])
        end
    end

    # If model convex approx, then...
    if flag_model == 0
        f_0[ref[1]] == 0

        @constraint(m_ac, dual1[i = 1:system.nbus], c[i,i] >= e[i]^2 + f[i]^2)
        @constraint(m_ac, [i = 1:system.nbus], c[i,i] <= 2*(e_0[i]*e[i] + f_0[i]*f[i]) - e_0[i]^2 - f_0[i]^2 + sl)  

        for k = 1:system.nbranch        
            i = convert(Int, system.Branchfrom[k])
            j = convert(Int, system.Branchto[k]) 
            if k == Line
                println("Line out of operation $k")
            else
                if flag_solver >= 1 
                    @constraint(m_ac, [(2*(e[i] - e[j])*(e_0[i] - e_0[j]) +
                                    2*(f[i] - f[j])*(f_0[i] - f_0[j]) -
                                    ((e_0[i] - e_0[j])^2 + (f_0[i] - f_0[j])^2) + sl + 4*c[i,j] + 1)/2,
                                    (2*(e[i] - e[j])*(e_0[i] - e_0[j]) +
                                    2*(f[i] - f[j])*(f_0[i] - f_0[j]) -
                                    ((e_0[i] - e_0[j])^2 + (f_0[i] - f_0[j])^2) + sl + 4*c[i,j] - 1)/2,
                                    e[i] + e[j], f[i] + f[j]] in SecondOrderCone())
                    @constraint(m_ac, [(2*(e[i] + e[j])*(e_0[i] + e_0[j]) +
                                    2*(f[i] + f[j])*(f_0[i] + f_0[j]) -
                                    ((e_0[i] + e_0[j])^2 + (f_0[i] + f_0[j])^2) + sl - 4*c[i,j] + 1)/2,
                                    (2*(e[i] + e[j])*(e_0[i] + e_0[j]) +
                                    2*(f[i] + f[j])*(f_0[i] + f_0[j]) -
                                    ((e_0[i] + e_0[j])^2 + (f_0[i] + f_0[j])^2) + sl - 4*c[i,j] - 1)/2,
                                    e[i] - e[j],f[i] - f[j]] in SecondOrderCone()) 

                    @constraint(m_ac, [(2*(e[i] + f[j])*(e_0[i] + f_0[j]) +
                                        2*(e[j] - f[i])*(e_0[j] - f_0[i]) -
                                    ((e_0[i] + f_0[j])^2 + (e_0[j] - f_0[i])^2) + sl - 4*s[i,j] + 1)/2, 
                                    (2*(e[i] + f[j])*(e_0[i] + f_0[j]) +
                                        2*(e[j] - f[i])*(e_0[j] - f_0[i]) -
                                    ((e_0[i] + f_0[j])^2 + (e_0[j] - f_0[i])^2) + sl - 4*s[i,j] - 1)/2, 
                                    e[i] - f[j], e[j] + f[i]] in SecondOrderCone())
                    @constraint(m_ac, [(2*(e[i] - f[j])*(e_0[i] - f_0[j]) +
                                        2*(e[j] + f[i])*(e_0[j] + f_0[i]) -
                                    ((e_0[i] - f_0[j])^2 + (e_0[j] + f_0[i])^2) + sl + 4*s[i,j] + 1)/2, 
                                    (2*(e[i] - f[j])*(e_0[i] - f_0[j]) +
                                        2*(e[j] + f[i])*(e_0[j] + f_0[i]) -
                                    ((e_0[i] - f_0[j])^2 + (e_0[j] + f_0[i])^2) + sl + 4*s[i,j] - 1)/2, 
                                    e[i] + f[j], e[j] - f[i]] in SecondOrderCone())                      
                else
                    @constraint(m_ac, (e[i] + e[j])^2 + (f[i] + f[j])^2 - 4*c[i,j] <=
                                    2*(e[i] - e[j])*(e_0[i] - e_0[j]) +
                                    2*(f[i] - f[j])*(f_0[i] - f_0[j]) -
                                    ((e_0[i] - e_0[j])^2 + (f_0[i] - f_0[j])^2) + sl)
                    @constraint(m_ac, (e[i] - e[j])^2 + (f[i] - f[j])^2 + 4*c[i,j] <=
                                    2*(e[i] + e[j])*(e_0[i] + e_0[j]) +
                                    2*(f[i] + f[j])*(f_0[i] + f_0[j]) -
                                    ((e_0[i] + e_0[j])^2 + (f_0[i] + f_0[j])^2) + sl)

                    @constraint(m_ac, (e[i] - f[j])^2 + (e[j] + f[i])^2 + 4*s[i,j] <=
                                    2*(e[i] + f[j])*(e_0[i] + f_0[j]) +
                                    2*(e[j] - f[i])*(e_0[j] - f_0[i]) +
                                    ((e_0[i] + f_0[j])^2 + (e_0[j] - f_0[i])^2) + sl)    
                    @constraint(m_ac, (e[i] + f[j])^2 + (e[j] - f[i])^2 - 4*s[i,j] <=
                                    2*(e[i] - f[j])*(e_0[i] - f_0[j]) +
                                    2*(e[j] + f[i])*(e_0[j] + f_0[i]) -
                                    ((e_0[i] - f_0[j])^2 + (e_0[j] + f_0[i])^2) + sl)
                end
            end     
        end      
    end

    

    #P_balance
    @constraint(m_ac, [k in 1:system.nbus], sum(Pg[i] for i in DispGen_bus[k])
        - Pd[l,k] - c[k,k] * Bus.Gs[k]
        - sum(Pfr[i] for i in system.out_lines[k])
        - sum(Pto[i] for i in system.in_lines[k]) == 0)

    #Q_balance
    @constraint(m_ac, [k in 1:system.nbus], sum(Qg[i] for i in DispGen_bus[k])
        - Qd[l,k] + c[k,k] * Bus.Bs[k]
        - sum(Qfr[i] for i in system.out_lines[k])
        - sum(Qto[i] for i in system.in_lines[k]) == 0)

    for k = 1:system.nbranch
        i = convert(Int, system.Branchfrom[k])
        j = convert(Int, system.Branchto[k])
        if k == Line
            println("Line out of operation $k")
        else
        #Pfr
        @constraint(m_ac, Pfr[k] == Branch.g[k] * Branch.a[k]^2 * c[i,i]
                                    - Branch.a[k] * Branch.g[k] *(c[i,j]*cosA[k]+s[i,j]*sinA[k])
                                    + Branch.a[k] * Branch.b[k] * (s[i,j]*cosA[k]-c[i,j]*sinA[k]))
                                   
        #Qfr
        @constraint(m_ac, Qfr[k] == -(Branch.b[k]+Branch.b_s[k]) * Branch.a[k]^2 * c[i,i]
                                    + Branch.a[k] * Branch.g[k] * (s[i,j]*cosA[k]-c[i,j]*sinA[k])
                                    + Branch.a[k] * Branch.b[k] * (c[i,j]*cosA[k]+s[i,j]*sinA[k]))

        #Pto
        @constraint(m_ac, Pto[k] == Branch.g[k] * c[j,j]
                                    - Branch.a[k] * Branch.g[k] * (c[i,j]*cosA[k]+s[i,j]*sinA[k])
                                    - Branch.a[k] * Branch.b[k] * (s[i,j]*cosA[k]-c[i,j]*sinA[k]))


        #Qto
        @constraint(m_ac, Qto[k] == -(Branch.b[k] + Branch.b_s[k]) * c[j,j]
                                    + Branch.a[k] * Branch.b[k] * (c[i,j]*cosA[k]+s[i,j]*sinA[k])
                                    - Branch.a[k] * Branch.g[k] * (s[i,j]*cosA[k]-c[i,j]*sinA[k]))
                                    
        end                   
    end


    # #  ============= Variables Fixed ================================
    # Voltage magnitude bounds per node and period:
    @constraint(m_ac, [i in 1:system.nbus], c[i,i] <= Bus.Vmax[i]^2)
    @constraint(m_ac, [i in 1:system.nbus], c[i,i] >= Bus.Vmin[i]^2)

    #Power Branch Limits
    for k in 1:system.nbranch

        if k == Line
            println("Line out of operation $k")
        else
            if flag_solver == 0 
                @NLconstraint(m_ac, Pfr[k]^2 + Qfr[k]^2 <= Branch.rateA[k]^2)
                @NLconstraint(m_ac, Pto[k]^2 + Qto[k]^2 <= Branch.rateA[k]^2)
            else
                @constraint(m_ac, [Branch.rateA[k],Pfr[k],Qfr[k]] in SecondOrderCone())
                @constraint(m_ac, [Branch.rateA[k],Pto[k],Qto[k]] in SecondOrderCone())
            end
        end

    end

    # if flag_model == 2 
    for k = 1:system.nbranch
        i = Int(system.Branchfrom[k])
        j = Int(system.Branchto[k])
        if k == Line
            println("Line out of operation $k")
        else
            if flag_solver != 0 
                @constraint(m_ac, [c[i,i], c[j,j], sqrt(2)*c[i,j], sqrt(2)*s[i,j]] in RotatedSecondOrderCone())
            else
                @NLconstraint(m_ac, c[i,j]^2 + s[i,j]^2 >= c[i,i]*c[j,j])
            end
        end
    end
    # end

    # Active power bounds per unit:
    @constraint(m_ac, [i in 1:system.ngen], Pg[i] <= Gen.Pmax[i])
    @constraint(m_ac, [i in 1:system.ngen], Pg[i] >= Gen.Pmin[i])
    # Reactive power bounds per unit:
    @constraint(m_ac, [i in 1:system.ngen], Qg[i] <= Gen.Qmax[i])
    @constraint(m_ac, [i in 1:system.ngen], Qg[i] >= Gen.Qmin[i])


    Output = Dict()
    Output["neqcon"] = num_constraints(m_ac, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64}) +
                       num_constraints(m_ac, VariableRef, MOI.EqualTo{Float64})
    Output["nincon"] = num_constraints(m_ac, GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64}) +
                       num_constraints(m_ac, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64}) +
                       num_constraints(m_ac, VariableRef, MOI.GreaterThan{Float64}) +
                       num_constraints(m_ac, VariableRef, MOI.LessThan{Float64})
    Output["ncvar"] = num_variables(m_ac)

    start = time()
    status = optimize!(m_ac)
    elapsed = time() - start;

    Output["flag"] = has_values(m_ac)
    print(termination_status(m_ac))

    Output["time"] = elapsed
    Output["Pg"] = zeros(system.ngen,1)
    Output["Qg"] = zeros(system.ngen,1)
    Output["Pfr"] = zeros(system.nbranch,1)
    Output["Qfr"] = zeros(system.nbranch,1)
    Output["Pto"] = zeros(system.nbranch,1)
    Output["Qto"] = zeros(system.nbranch,1)
    Output["e"] = zeros(system.nbus,1)
    Output["f"] = zeros(system.nbus,1)
    Output["V"] = zeros(system.nbus,1)
    Output["theta"] = zeros(system.nbus,1)
    Output["c"] = zeros(system.nbus,system.nbus,1)
    Output["s"] = zeros(system.nbus,system.nbus,1)
    Output["dual1"] = zeros(system.nbus,1)
    Output["dual2"] = zeros(system.nbranch,1)
    Output["dual3"] = zeros(system.nbranch,1)
    Output["LMP_P"] = zeros(system.nbus,1)
    Output["slack"] = value(sl)
    Output["Objective"] = objective_value(m_ac)
    if flag_model == 3
        Output["lambda"] = value(lambda)
    end

    if has_values(m_ac)
        for i in 1:system.nbus
            Output["e"][i] = value(e[i])
            Output["f"][i] = value(f[i])
            Output["c"][i,i] = value(c[i,i])
            Output["V"][i] = sqrt(value(c[i,i]))
            Output["theta"][i] = atand(value(f[i])/value(e[i]))
        end
        for i in 1:system.nbranch     
            Output["Pfr"][i] = value(Pfr[i])   
            Output["Qfr"][i] = value(Qfr[i])   
            Output["Pto"][i] = value(Pto[i])   
            Output["Qto"][i] = value(Qto[i])   
        end
        for i in 1:system.ngen
            Output["Pg"][i] = value(Pg[i])
            Output["Qg"][i] = value(Qg[i])
        end
        for k in 1:system.nbranch
            i = Int(system.Branchfrom[k])
        	j = Int(system.Branchto[k])
            Output["c"][i,j] = value(c[i,j])
            Output["s"][i,j] = value(s[i,j])
        end
    end
    return(Output)
end

# ============MAIN EXECUTION FLOW=============

# Initialize variables for objective values and initial guesses

e_0 = ones(system.nbus, 1)
f_0 = zeros(system.nbus, 1)
ObjectiveAC = zeros(Int(size(PgAC,1)), 1)
ObjectivePDC = zeros(Int(size(PgAC,1)), 1)
ObjectiveDC = zeros(Int(size(PgAC,1)), 1)
MODEL = 1  # Set model to nonconvex by default
SOLVER = 1  # Use Gurobi solver by default


    for sample in 1:1000
        println("Sample PDCOPF",l)
        Solution_NCVX = ACOPF(system, Branch, Bus, Gen, e_0, f_0, [], 0, MODEL, SOLVER, 1e-4,sample,Pd,Qd,PgPDC,[])
        Results = Dict()
        ObjectivePDC[l,:] .= Solution_NCVX["Objective"]
        Pg_PDC_feasible[l,:] = Solution_NCVX["Pg"]
    end

    pDC_Model["Feasibility_Cost"] = copy(ObjectivePDC)
    pDC_Model["Pg_Feasibility_Cost"] = copy(Pg_PDC_feasible)


    for sample in 1:1000
        println("Sample DCOPF",sample)
        Solution_NCVX = ACOPF(system, Branch, Bus, Gen, e_0, f_0, [], 0, MODEL, SOLVER, 1e-4,sample,Pd,Qd,PgDC,[])
        Results = Dict()
        ObjectiveDC[l,:] .= Solution_NCVX["Objective"]
        Pg_DC_feasible[l,:] = Solution_NCVX["Pg"]
    end
    DC_Model["Feasibility"] = copy(ObjectiveDC)
    DC_Model["Pg_Feasibility"] = copy(Pg_DC_feasible)


#MAT.matwrite("DC_Model_$system_name.mat", DC_Model)
#MAT.matwrite("pDC_Model_$system_name.mat", pDC_Model)


