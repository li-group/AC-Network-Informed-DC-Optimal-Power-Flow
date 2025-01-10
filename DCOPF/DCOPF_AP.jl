
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

# Define data path
path_data = "C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos"
# Store original path
path_org = pwd()

# Change directory to the data path
cd(path_data)

# Define the power system case
system_name = "pglib_opf_case57_ieee"
# system_name = "pglib_opf_case30_ieee"
# system_name = "pglib_opf_case118_ieee"

# Dictionary to store case information
Cases = Dict()
Cases[system_name] = mat"run(append($system_name,'.m'))"

# Generate structures for the system
include("CreateStructures_OPF.jl")

# # ===================DEMANDS & DISPATCH IMPORT=====================#

# Read Matpower data
path_data = "C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Matpower_CPF\\Generation_Matpower_$system_name.mat"
DataMATPOWER = matread(path_data)

# Read beta values from CSV
betas = Matrix(CSV.read("C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\NeuralNetwork\\Beta_nn3_$system_name.csv", DataFrame))

# Define flag for operation mode: 0 => Normal, 1 => Contingencies
flag_operation = 0

# Extract required data from the MATPOWER structure
PgAC = DataMATPOWER["DataMATPOWER"]["DispatchAC_Probe"]
LMPS_AC = DataMATPOWER["DataMATPOWER"]["DualAC_Probe"] * 100
load_power = DataMATPOWER["DataMATPOWER"]["Demand_Probe"]

BBUS0 = DataMATPOWER["DataMATPOWER"]["BBUS"]
BF0 = DataMATPOWER["DataMATPOWER"]["BF"]
PBUSINJ = DataMATPOWER["DataMATPOWER"]["PBUSINJ"]
PFINJ = DataMATPOWER["DataMATPOWER"]["PFINJ"]

# Restore original directory path
cd(path_org)

# Reference bus index
ref = findall(x -> x == 3, Bus.bustype)

# Generate generator index per bus
DispGen_bus = [Int[] for i in 1:system.nbus]
for i in 1:system.ngen
    a = convert(Int, system.GenBusnumber[i])
    push!(DispGen_bus[a], i)
end

# Remove duplicate neighboring buses
system.neig_bus = [unique(inner_vector) for inner_vector in system.neig_bus]

# Generator bus indices
GenBus = permutedims(Int.(system.GenBusnumber))
GenBus = GenBus[:, 1]

# DCOPF PARAMETRIC FUNCTION
function DCOPF(system, Branch, flag_model, flag_solver, mipgap, beta, l, Pd, line, BBUS, BF)
    # Define solver based on flag
    if flag_solver == 0
        m_dc = Model(Ipopt.Optimizer)
        set_optimizer_attribute(m_dc, "max_iter", 10000)
    elseif flag_solver == 1
        m_dc = Model(Gurobi.Optimizer)
        set_optimizer_attribute(m_dc, MOI.Silent(), true)
        if flag_model == 1
            set_optimizer_attribute(m_dc, "MIPGap", mipgap)
            set_optimizer_attribute(m_dc, "TimeLimit", 60)
        end
    elseif flag_solver == 2
        m_dc = Model(Mosek.Optimizer)
        set_optimizer_attribute(m_dc, MOI.Silent(), true)
        set_optimizer_attribute(m_dc, "MSK_DPAR_INTPNT_CO_TOL_INFEAS", 1e-6)
    end 

    # Define variables for the optimization
    @variable(m_dc, theta[1:system.nbus])
    @variable(m_dc, Pg[1:system.ngen] >= 0)
    @variable(m_dc, BaseCost >= 0)
    @variable(m_dc, Pd[l, i] >= slack1[i = 1:system.nbus] >= 0)
    @variable(m_dc, slack2[i = 1:system.nbus] == 0)


    # Define the objective function
    @objective(m_dc, Min, BaseCost)

    @constraint(m_dc, BaseCost >= sum(Gen_cost.ck1[i]*Pg[i] for i in 1:system.ngen)
    +1e4*sum(slack1[i]+slack2[i] for i in 1:system.nbus)) 

    # Balance constraint adjustment if contingencies are modeled
    BBUS1 = copy(BBUS)
    BF1 = copy(BF)
    if flag_operation == 1
        for k in 1:system.nbranch    
            if k == line
                i = convert(Int, system.Branchfrom[k])
                j = convert(Int, system.Branchto[k])
                BBUS1[i, i] -= 1 / Branch.x[k]
                BBUS1[j, j] -= 1 / Branch.x[k]
                BBUS1[i, j] += 1 / Branch.x[k]
                BBUS1[j, i] += 1 / Branch.x[k]
            end
        end
    end

    # Bus balance constraint
    @constraint(m_dc, Balance[k = 1:system.nbus], sum(Pg[m] for m in DispGen_bus[k]) ==
                      Pd[l, k] * beta[l, k] + Bus.Gs[k] + PBUSINJ[k] + sum(BBUS1[k, j] * theta[j] for j in 1:system.nbus) - slack1[k])

    @constraint(m_dc, theta[ref] == 0)

    # Line limits
    for k in 1:system.nbranch  
        if k == line
            BF1[line, :] .= 0
        end
        @constraint(m_dc, -Branch.rateA[k] <= sum(BF1[k, i] * theta[i] for i = 1:system.nbus) + PFINJ[k] <= Branch.rateA[k])
    end

    # Generator power bounds
    @constraint(m_dc, [i in 1:system.ngen], Pg[i] <= Gen.Pmax[i])
    @constraint(m_dc, [i in 1:system.ngen], Pg[i] >= Gen.Pmin[i])

    # Data storage for results
    Output = Dict()
    start = time()
    optimize!(m_dc)
    elapsed = time() - start
    
    if has_values(m_dc)
        Output["cost"] = objective_value(m_dc)
        Output["time"] = elapsed
        Output["BaseCost"] = value(BaseCost)
        Output["theta"] = JuMP.value.(theta)
        Output["slack1"] = JuMP.value.(slack1)
        Output["LMPs"] = JuMP.dual.(Balance)
    end    
    Output["Pg"] = zeros(system.ngen, 1)
    if has_values(m_dc)
        for i in 1:system.ngen
            Output["Pg"][i] = value(Pg[i])
        end
    end
    return Output
end


# ERROR FUNCTION FOR COMPARING MODELS
function ERROR(AC, DC)
    relative_error = (abs.(AC - DC) ./ AC) * 100
    unsigned_error = (AC - DC) ./ AC * 100
    total_error = sum(relative_error)
    Output = Dict()
    Output["relative_error"] = relative_error
    Output["total_error"] = total_error
    Output["unsigned_error"] = unsigned_error
    return Output
end

# Revenue Adequacy calculation function
function RevenueAdequency(Dual, Pg, Slack)
    sum_LMPs_Pg = sum(Dual[:, GenBus] .* Pg, dims = 2)
    sum_LMPs_Pd = sum(Dual .* (Pd - Slack), dims = 2)
    sum_cw = sum(1e4 * Slack, dims = 2)
    RevenueAdequencyPDC = []
    for i in 1:SAMPLES
        if sum_LMPs_Pd[i] - sum_LMPs_Pg[i] - sum_cw[i] >= 0
            push!(RevenueAdequencyPDC, i)
        end
    end
    Output = Dict()
    Output["RA_OK"] = RevenueAdequencyPDC
    Output["RA_Result"] = sum_LMPs_Pd - sum_LMPs_Pg - sum_cw
    return Output
end

# Cost Recovery function
function CostRecovery(Dual, Pg, Gen_cost)
    LMPs_Pg = Dual[:, GenBus] .* Pg
    Cg_Pg = Pg .* Gen_cost'
    CostRecoveryPDC = []
    for i in 1:SAMPLES
        for j in 1:system.ngen
            if LMPs_Pg[i, j] - Cg_Pg[i, j] >= 0
                push!(CostRecoveryPDC, j)
            end
        end
    end
    Output = Dict()
    Output["CR_OK"] = CostRecoveryPDC
    Output["CR_Result"] = LMPs_Pg - Cg_Pg
    return Output
end



## ========================FUNCTION PARAMETERS=================================#
RANGE_FROM = 1
RANGE_TO = 1000
MODEL = 1
SOLVER = 1
SAMPLES = RANGE_TO-RANGE_FROM+1

# Extract relevant data from previously loaded files
betas = betas[RANGE_FROM:RANGE_TO, :]
Pd = load_power[RANGE_FROM:RANGE_TO, 1:system.nbus]
PgAC = PgAC[RANGE_FROM:RANGE_TO, 1:system.ngen]

# Initialize matrices to store data
PgDC = zeros(SAMPLES, system.ngen)
timeDC = zeros(SAMPLES, 1)
DualDC = zeros(SAMPLES, system.nbus)
Slack1DC = zeros(SAMPLES, system.nbus)
PdRenew = zeros(SAMPLES, system.nbus)

PgPDC = zeros(SAMPLES, system.ngen)
timePDC = zeros(SAMPLES, 1)
DualPDC = zeros(SAMPLES, system.nbus)
Slack1PDC = zeros(SAMPLES, system.nbus)
Slack2PDC = zeros(SAMPLES, system.nbus)

global ErrorCounter = 0
no_solution_sampleDC = []

# Define line index for contingencies if flag_operation == 1
if flag_operation == 1
    line = 1:system.nbranch  # contingency cases
else
    line = []
end

# Run the DCOPF model with or without contingencies
if flag_operation == 1
    # Contingency DC
    Contingencias_DC = Dict()
    for lines in line
        println("Contingency ", lines)
        for l in 1:SAMPLES
            println("Sample ", l)
            try
                Solution_DC = DCOPF(system, Branch, MODEL, SOLVER, 1e-4, ones(SAMPLES, system.nbus), l, Pd, lines, BBUS0, BF0)
                PgDC[l, :] = Solution_DC["Pg"]
                timeDC[l, :] .= Solution_DC["time"]
                DualDC[l, :] = Solution_DC["LMPs"]
                Slack1DC[l, :] .= Solution_DC["slack1"]
            catch
                push!(no_solution_sampleDC, l)                
            end
        end
        Contingencias_DC["Pg_$lines"] = copy(PgDC)
        Contingencias_DC["time_$lines"] = copy(timeDC)
        Contingencias_DC["Dual_$lines"] = copy(DualDC)
        Contingencias_DC["Slack1_$lines"] = copy(Slack1DC)
    end
else
    # Run the DCOPF without contingencies
    for l in 1:SAMPLES
        Solution_DC = DCOPF(system, Branch, MODEL, SOLVER, 1e-4, ones(SAMPLES, system.nbus), l, Pd, line, BBUS0, BF0)
        PgDC[l, :] = Solution_DC["Pg"]
        timeDC[l, :] .= Solution_DC["time"]
        DualDC[l, :] = Solution_DC["LMPs"]
        Slack1DC[l, :] .= Solution_DC["slack1"]
    end
    DC_Model["Pg"] = copy(PgDC)
    DC_Model["time"] = copy(timeDC)
    DC_Model["Dual"] = copy(DualDC)
    DC_Model["Slack1"] = copy(Slack1DC)
end

no_solution_sample = []

# DCOPF PARAMETRIC FUNCTION with Contingencies
if flag_operation == 1
    # Parametric DC with contingencies
    Contingencias_pDC = Dict()
    for lines in line
        println("Contingency ", lines)
        for l in 1:SAMPLES
            println("Sample ", l)
            try
                Solution_PDC = DCOPF(system, Branch, MODEL, SOLVER, 1e-4, betas, l, Pd, lines, BBUS0, BF0)
                PgPDC[l, :] = Solution_PDC["Pg"]
                timePDC[l, :] .= Solution_PDC["time"]
                DualPDC[l, :] = Solution_PDC["LMPs"]
                Slack1PDC[l, :] .= Solution_PDC["slack1"]

                # Adjust demand and rerun if slack is present
                if sum(Slack1PDC[l, :]) != 0
                    PdRenew = Pd[l, :] - Slack1PDC[l, :]
                    Solution_PDC = DCOPF(system, Branch, MODEL, SOLVER, 1e-4, betas, l, PdRenew, lines, BBUS0, BF0)
                    PgPDC[l, :] = Solution_PDC["Pg"]
                    timePDC[l, :] .= Solution_PDC["time"]
                    DualPDC[l, :] = Solution_PDC["LMPs"]
                    Slack1PDC[l, :] .= Solution_PDC["slack1"]
                end
            catch
                push!(no_solution_sample, l)
            end
        end

        Contingencias_pDC["Pg_$lines"] = copy(PgPDC)
        Contingencias_pDC["time_$lines"] = copy(timePDC)
        Contingencias_pDC["Dual_$lines"] = copy(DualPDC)
        Contingencias_pDC["Slack1_$lines"] = copy(Slack1PDC)
    end
else
    # Run the parametric DCOPF without contingencies
    for l in 1:SAMPLES
        try
            Solution_PDC = DCOPF(system, Branch, MODEL, SOLVER, 1e-4, betas, l, Pd, line, BBUS0, BF0)
            PgPDC[l, :] = Solution_PDC["Pg"]
            timePDC[l, :] .= Solution_PDC["time"]
            DualPDC[l, :] = Solution_PDC["LMPs"]
            Slack1PDC[l, :] .= Solution_PDC["slack1"]

            # Adjust demand and rerun if slack is present
            if sum(Slack1PDC[l, :]) != 0
                PdRenew[l, :] = Pd[l, :] - Slack1PDC[l, :]
                Solution_PDC = DCOPF(system, Branch, MODEL, SOLVER, 1e-4, betas, l, PdRenew, line, BBUS0, BF0)
                PgPDC[l, :] = Solution_PDC["Pg"]
                timePDC[l, :] .= Solution_PDC["time"]
                DualPDC[l, :] .= Solution_PDC["LMPs"]
                Slack1PDC[l, :] .= Solution_PDC["slack1"]
            end
        catch
            push!(no_solution_sample, l)
        end
        pDC_Model["Pg"] = copy(PgPDC)
        pDC_Model["time"] = copy(timePDC)
        pDC_Model["Dual"] = copy(DualPDC)
        pDC_Model["Slack1"] = copy(Slack1PDC)
    end
end            
  
    

# Cost comparison between AC, pDC, and DC models
if flag_operation == 1
    Contingencias_AC["Pg"] = DataMATPOWER["DataMATPOWER"]["ContingenciaAC"]
    for i in line
        Contingencias_AC["Pg_$i"] = DataMATPOWER["DataMATPOWER"]["ContingenciaAC"]["Pg_$i"][2:1001, :]
        Contingencias_pDC["BaseCost_$i"] = reshape(Contingencias_pDC["Pg_$i"] * Gen_cost.ck1, (length(Contingencias_pDC["Pg_$i"] * Gen_cost.ck1), 1)) 
        Contingencias_DC["BaseCost_$i"] = reshape(Contingencias_DC["Pg_$i"] * Gen_cost.ck1, (length(Contingencias_DC["Pg_$i"] * Gen_cost.ck1), 1)) 
        Contingencias_AC["BaseCost_$i"] = reshape(DataMATPOWER["DataMATPOWER"]["ContingenciaAC"]["Pg_$i"][2:1001, :] * Gen_cost.ck1, (length(DataMATPOWER["DataMATPOWER"]["ContingenciaAC"]["Pg_$i"][2:1001, :] * Gen_cost.ck1), 1)) 
    end

    # Calculate cost errors between the methods
    for i in line
        Contingencias_pDC["ErrorCost_$i"] = ERROR(Contingencias_AC["BaseCost_$i"], Contingencias_pDC["BaseCost_$i"])
        Contingencias_DC["ErrorCost_$i"] = ERROR(Contingencias_AC["BaseCost_$i"], Contingencias_DC["BaseCost_$i"])
    end

else
    # Cost calculation without contingencies
    BaseCostAC = reshape(PgAC * Gen_cost.ck1, (length(PgAC * Gen_cost.ck1), 1))
    pDC_Model["BaseCost"] = reshape(PgPDC * Gen_cost.ck1, (length(PgPDC * Gen_cost.ck1), 1))
    DC_Model["BaseCost"] = reshape(PgDC * Gen_cost.ck1, (length(PgDC * Gen_cost.ck1), 1))

    # Calculate cost errors between the methods
    ErrorDC = ERROR(BaseCostAC, DC_Model["BaseCost"])
    ErrorPDC = ERROR(BaseCostAC, pDC_Model["BaseCost"])
    pDC_Model["ErrorCost"] = ErrorPDC
    DC_Model["ErrorCost"] = ErrorDC

    # Calculate dispatch errors between AC, pDC, and DC models
    ErrorDCPg = ERROR(PgAC, PgDC)
    ErrorPDCPg = ERROR(PgAC, PgPDC)
    pDC_Model["ErrorPg"] = ErrorPDCPg
    DC_Model["ErrorPg"] = ErrorDCPg

    # Revenue adequacy calculation for pDC
    RevenueAdequency_PDC = RevenueAdequency(pDC_Model["Dual"], pDC_Model["Pg"], pDC_Model["Slack1"])
    pDC_Model["RA_OK"] = RevenueAdequency_PDC["RA_OK"]
    pDC_Model["RA_Result"] = RevenueAdequency_PDC["RA_Result"]

    # Revenue adequacy calculation for DC
    RevenueAdequency_DC = RevenueAdequency(DC_Model["Dual"], DC_Model["Pg"], DC_Model["Slack1"])
    DC_Model["RA_OK"] = RevenueAdequency_DC["RA_OK"]
    DC_Model["RA_Result"] = RevenueAdequency_DC["RA_Result"]

    # Cost recovery calculation for pDC
    CostRecovery_PDC = CostRecovery(pDC_Model["Dual"], pDC_Model["Pg"], Gen_cost.ck1)
    pDC_Model["CR_OK"] = CostRecovery_PDC["CR_OK"]
    pDC_Model["CR_Result"] = CostRecovery_PDC["CR_Result"]

    # Cost recovery calculation for DC
    CostRecovery_DC = CostRecovery(DC_Model["Dual"], DC_Model["Pg"], Gen_cost.ck1)
    DC_Model["CR_OK"] = CostRecovery_DC["CR_OK"]
    DC_Model["CR_Result"] = CostRecovery_DC["CR_Result"]
end

# Export data to MATLAB format for further analysis
if flag_operation == 1
    MAT.matwrite("Contingencias_pDC_$system_name.mat", Contingencias_pDC)
    MAT.matwrite("Contingencias_DC_$system_name.mat", Contingencias_DC)
    MAT.matwrite("Contingencias_AC_$system_name.mat", Contingencias_AC)
else
    MAT.matwrite("pDC_Model_$system_name.mat", pDC_Model)
    MAT.matwrite("DC_Model_$system_name.mat", DC_Model)
end

#MAT.matwrite("pDC_Model_unconstrain_$system_name.mat", pDC_Model)
#MAT.matwrite("DC_Model_unconstrain_$system_name.mat", DC_Model)
#MAT.matwrite("pDC_Model_new_$system_name.mat", pDC_Model)
#MAT.matwrite("DC_Model_new_$system_name.mat", DC_Model)

