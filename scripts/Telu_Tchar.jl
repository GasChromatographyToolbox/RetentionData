# script to simulate different conventional GC runs
# to correlate Telu and Tchar in relation to different parameters
using Pkg
Pkg.activate(".")
using CSV, DataFrames, DrWatson, GasChromatographySimulator

#variation of parameters:
# stationary phase
sp = ["SLB5ms","SPB50","Wax","Rxi17SilMS","DB5ms","Rxi5ms"]
# initial temperature:
Tinit = collect(-40.0:40.0:320.0)
# heating rate:
RT = [0.001, 0.01, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
# inlet pressure:
pin = collect(50.0:50.0:300.0)
# column length:
L = 30.0#collect(10.0:10.0:60.0)
# diameter and film thickness (phase ratio constant):
ddf = 0.25#[0.1, 0.25, 0.32, 0.53]
# gas
gas = "He"#["H2", "He", "N2"]

parameters = Dict(
					"L" => L, #m
	                "d" => ddf.*1e-3, #m
	                "df" => ddf.*1e-6, #m
	                "sp" => sp,
	                "gas" => gas,
	                "rate" => RT, #°C/min
	                "Tinit" => Tinit, #°C
					"Tend" => 1000.0, #°C
	                "constMode" => "Pressure", #"Pressure" or "Flow"
					##:constMode => "Flow",
	                "Fpin" => pin, # in kPa(g) or mL/min
					##:Fpin => 1.0, # in mL/min
	                "pout" => 0.0, # in kPa(a)
	                "solute_db_path" => "/Users/janleppert/Documents/GitHub/ThermodynamicData/Databases/",
	                "solute_db" => "newformat_Kcentric_2022-05-25T12:49:45.405.csv",
	                "abstol" => 1e-6,
	                "reltol" => 1e-3
	)
dicts = dict_list(parameters);
Ndict = dict_list_count(parameters)

function sortout_CASmissing(solutes, db)
	new_solutes = String[]
	for i=1:length(solutes)
		try
			cas = GasChromatographySimulator.CAS_identification([solutes[i]])[!, 2][1]
			push!(new_solutes, solutes[i])
		catch
			cas = Missing
		end
	end
	return new_solutes
end

function makesim(dict::Dict) #(adapted from 'script-sim-heating_rate.jl' in VGGC-project)
    @unpack L, d, df, sp, gas, rate, Tinit, Tend, constMode, Fpin, pout, solute_db_path, solute_db, abstol, reltol = dict

    Tsteps = [Tinit, Tend]
    Δtsteps = [0.0, (Tend-Tinit)./rate.*60.0]
	if constMode == "Pressure"
    	Fpinsteps = [Fpin, Fpin].*1000.0 .+ 101300.0 # conversion from kPa(g) to Pa(a)
		# calculate holdup-time at reference temperature 150°C
		tMref = GasChromatographySimulator.holdup_time(150.0+273.15, Fpin*1000.0+101300.0, pout, L, d, gas, control="Pressure")
	else
		Fpinsteps = [Fpin, Fpin]./60e6 # conversion from mL/min to m³/s
		# calculate holdup-time at reference temperature 150°C
		tMref = GasChromatographySimulator.holdup_time(150.0+273.15, Fpin/60e6, pout, L, d, gas, control="Flow")
	end
	poutsteps = [pout, pout].*1000.0
    sys = GasChromatographySimulator.Column(L,d,df,sp,gas)
    prog = GasChromatographySimulator.Program(Δtsteps,Tsteps,Fpinsteps,poutsteps,L)
	db = DataFrame(CSV.File(string(solute_db_path,"/",solute_db), header=1, silencewarnings=true))
	solutes = sortout_CASmissing(GasChromatographySimulator.all_solutes(sp, db), db)
	# the usage of a custom database for ChemicalIdentifiers seems not to work
    sub = GasChromatographySimulator.load_solute_database(db, sp, gas, solutes, zeros(length(solutes)), zeros(length(solutes)))
    opt = GasChromatographySimulator.Options(abstol=abstol, reltol=reltol, control=constMode)
    
    par = GasChromatographySimulator.Parameters(sys, prog, sub, opt)
    # run the simulation 
    pl, sol = GasChromatographySimulator.simulate(par)
	
    fulldict = copy(dict)
    fulldict["peaklist"] = pl
    fulldict["par"] = par
	fulldict["tMref"] = tMref
    return fulldict
end

Nsim = Int[]
Nexist = Int[]
Threads.@threads for i=1:length(dicts)
    file = datadir("sims", "Telu_Tchar", savename("Telu_Tchar", dicts[i], "jld2"; ignores=["solute_db_path", "solute_db"]))
    if isfile(file) == false # check if the file already exist, make the simulation only if the file does not exist
        f = makesim(dicts[i])
        @tagsave(file, f)
        push!(Nsim, 1)
    else
        push!(Nexist, 1)
    end
end
println("$(sum(Nsim))/$(Ndict) simulated, $(sum(Nexist))/$(Ndict) already exist.")