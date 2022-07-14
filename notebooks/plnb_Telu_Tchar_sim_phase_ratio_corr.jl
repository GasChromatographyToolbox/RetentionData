### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 31dc694a-9cfb-495b-baa4-c30e85634744
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	using DrWatson, DataFrames, CSV, Interpolations, Plots, LsqFit, Measurements
	using PlutoUI, Statistics, LambertW
	using GasChromatographySimulator
	plotly()
	TableOfContents()
end

# ╔═╡ 42fd0ee9-cc2f-4889-b5d6-4431829aa5a7
md"""
# Simulation

## Correlation of Elution Temperature $T_{elu}$, characteristic Temperature $T_{char}$ and dimensionless Heating Rate $r_T$

Simulations were made with the script 'script-sim-heating_rate-phaseratio.jl'. A simple temperature program heating from 30°C to 360°C (no initial holding time). Heating rate, inlet pressure, **stationary phases**, **solutes** and **phase ratio $\beta$** were varied between different simulations. 
"""

# ╔═╡ 7e517e2a-7ab7-4e2f-8c35-3f6b37c70e33
begin
	# settings (adapted from 'script-sim-heating_rate-phaseratio.jl' in VGGC-project)
	# program: single ramp, no holding times
	rate = [1.0, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0, 30.0, 50.0] # °C/min
	Tinit = 30.0 # °C
	Tend = 360.0 # °C
	
	pin = [100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0]# kPa(g)
	pout = [101.3] # kPa(a)

	parameters = Dict(
	                :L => 30.0,
	                :d => 2.5e-4,
	                :df => [0.05e-6, 0.1e-6, 0.25e-6, 0.5e-6, 1.0e-6],
	                :sp => ["SLB5ms","SPB50","Wax","FS5ms","Rxi17SilMS","DB5ms","Rxi5MS","ZB-PAH-CT"],
	                :gas => "He",
	                :rate => rate,
	                :Tinit => Tinit,
					:Tend => Tend,
	                #const pressure is needed for constant difference between elution temp. and Tchar
	                :pin => pin, # in kPa(g)
	                :pout => pout, # in kPa(a)
	                :solute_db_path => projectdir("Databases"),
	                :solute_db => "Database_append.csv",
	                :abstol => 1e-6,
	                :reltol => 1e-3
	)
end

# ╔═╡ 2f202810-dc23-4760-8bb4-f00e4fd96f23
Ndict = dict_list_count(parameters)

# ╔═╡ b7f9ec8d-c637-401d-a47a-6d3c30568abd
dicts = dict_list(parameters);

# ╔═╡ 2bb5eb7c-7507-43ab-8c5e-0086f0e21a1c
projectdir("Databases")

# ╔═╡ 42e7ca1a-b7d8-4cd8-9998-02f27cedbc59
function makesim(dict::Dict) #(adapted from 'script-sim-heating_rate.jl' in VGGC-project)
    @unpack L, d, df, sp, gas, rate, Tinit, Tend, pin, pout, solute_db_path, solute_db, abstol, reltol = dict

    Tsteps = [Tinit, Tend]
    Δtsteps = [0.0, (Tend-Tinit)./rate.*60.0]
    pinsteps = [pin, pin].*1000.0 .+ 101300.0
	poutsteps = [pout, pout].*1000.0
    sys = GasChromatographySimulator.constructor_System(L,d,df,sp,gas)
    prog = GasChromatographySimulator.constructor_Program(Δtsteps,Tsteps,pinsteps,poutsteps,L)
	db = DataFrame(CSV.File(string(solute_db_path,"/",solute_db), header=1, silencewarnings=true))
	solutes = GasChromatographySimulator.all_solutes(sp, db)
    sub = GasChromatographySimulator.load_solute_database(db, sp, gas, solutes, zeros(length(solutes)), zeros(length(solutes)))
    opt = GasChromatographySimulator.Options(OwrenZen5(), abstol, reltol, "inlet", true)
    
    par = GasChromatographySimulator.Parameters(sys, prog, sub, opt)
    # run the simulation 
    sol = GasChromatographySimulator.solve_system_multithreads(par; ng=true)
    pl = GasChromatographySimulator.peaklist(sol, par)
    fulldict = copy(dict)
    fulldict[:peaklist] = pl
    fulldict[:par] = par
    return fulldict
end

# ╔═╡ dc2c53dd-89cc-4e93-bd58-e20b2b5750e6
# 5670 simulations in 655s

# ╔═╡ 7c1b33e7-4e58-401a-b914-3c2004ce6ca8
begin
	# run the simulation: 379s for 2835 simulations
	f = Array{Any}(undef, length(dicts))
	for (i, d) in enumerate(dicts)
	    f[i] = makesim(d)
	end
end

# ╔═╡ e8d0e6dd-d0a7-4692-baeb-71712022f23a
md"""
## Correction of the parameter $T_{char}$ for the phase ratio $\beta$

... some text (see notebook `plnb_lnk_beta.jl`)
"""

# ╔═╡ bf31b8ee-6222-4596-b177-3db9a91ac677
function dict_array_to_dataframe(f)
	df = DataFrame(f[1])
	for i=2:length(f)
		push!(df, f[i])
	end
	return df
end

# ╔═╡ a57cfb19-7718-4336-bead-fff25b19287c
data = dict_array_to_dataframe(f)

# ╔═╡ 388a7f39-8faa-4cfa-9a38-360a9e124cde
function extract_sub_sort(data, i;  solute_db_path=datadir("Databases"))
		# extract the thermodynamic parameters Tchar, θchar and ΔCp and sort them according to retention time
		# load database and extract Tchar (and name)
	db = DataFrame(CSV.File(string(solute_db_path,"/",data.solute_db[i]), header=1, silencewarnings=true))
	solutes = GasChromatographySimulator.all_solutes(data.sp[i], db)
	sub = GasChromatographySimulator.load_solute_database(db, data.sp[i], data.gas[i], solutes, zeros(length(solutes)), zeros(length(solutes)))
	# sort Tchar according to TR
	sub_sort = Array{GasChromatographySimulator.Substance}(undef, length(sub))
	for j=1:length(sub)
		sub_sort[findfirst(sub[j].name.==data.peaklist[i].Name)] = sub[j]
	end
	return sub_sort
end

# ╔═╡ e30cb724-1f9c-448c-9b59-977e4e1da4fc
function Tchar_β(β, Tchar0, θchar0, ΔCp0, β0)
	# correction of Tchar (estimated for φ₀) for the actual phase ratio β of the simulated system
	R = 8.3145
	a = θchar0*ΔCp0/R
	b = Tchar0/a+1
	Δlnβ = log(β0)-log(β)
	y = -b*exp(Δlnβ*R/ΔCp0 - b)
	if y>=(-1/exp(1)) && y<0
		Tchar = -Tchar0*b/lambertw(y,-1)
		branch = -1
	else
		Tchar = -Tchar0*b/lambertw(y,0)
		branch = 0
	end
	# -> criteria for the goodness of the thermodynamic parameters:
	# -1/e<=y<0
	# otherwise there are problems
	return Tchar, branch
end

# ╔═╡ 2ccd805a-fcf3-44c6-a79e-6db1ff6e08d0
function linear_fit_β(data, i;  solute_db_path=projectdir("Databases"), α=0.05)
	# linear fit of TR over Tchar with slope = 1
	# -> intercept = mean difference TR-Tchar
	sub_sort = extract_sub_sort(data, i; solute_db_path=solute_db_path)
	@. lin_model(x,p) = x + p[1]
	p0 = [0.0]
	# filter for tR<tProg, for tR>tProg elution temperature is constant
	tProg = (data.Tend[i]-data.Tinit[i])/data.rate[i]*60 # duration of heating ramp in s
	b = findfirst(tProg.<=data.peaklist[i].tR) # last element with tR<tProg
	if isa(b, Nothing)
		b = length(sub_sort)
	else
		b = b-1
	end
	# filter for unretained solutes, use μR to identify them as outlyer
	# iterativly calculate relativ standard deviation (std/mean) and reduce the first index until relative error is <α
	μR = 1.0 ./ (1.0 .+ data.peaklist[i].kR)
	a = 1
	if std(μR[1:b])/mean(μR[1:b]) < α
		a = 1
	else
		for ii=2:b
			if std(μR[ii:b])/mean(μR[ii:b]) < α
				a = ii
				break
			end
		end
	# change in the following line Tchar_sort to recalculated Tchar for the corresponding phase ratio β
	end
	mean_μR = mean(μR[a:b])
	std_μR = std(μR[a:b])
	Tchar_β_sort = Array{Float64}(undef, length(sub_sort))
	Tchar_sort = Array{Float64}(undef, length(sub_sort))
	branch = Array{Int64}(undef, length(sub_sort))
	for j=1:length(sub_sort)
		Tchar_β_sort[j], branch[j] = Tchar_β(data.β[i],sub_sort[j].Tchar,sub_sort[j].θchar,sub_sort[j].ΔCp,1/(4*sub_sort[j].φ₀))
		Tchar_sort[j] = sub_sort[j].Tchar
	end
	lin_fit = curve_fit(lin_model, Tchar_β_sort[a:b].-273.15, data.peaklist[i].TR[a:b], p0)
	return lin_fit, a, b, mean_μR, std_μR, Tchar_β_sort, Tchar_sort, branch
end

# ╔═╡ 5f27ef08-9237-4379-b15a-65cfa311881b
begin
	#data=collect_results!(datadir("sims","results_heating_rate-phase_ratio.jld2"),datadir("sims","heating_rate-phase_ratio"))
	sort!(data, [:sp, :pin, :rate, :df])
	data.β = zeros(length(data.rate))
	data.a = zeros(Int, length(data.rate))
	data.b = zeros(Int, length(data.rate))
	data.ΔTelu = measurement(0.0, 0.0).*ones(length(data.rate))
	data.μelu = measurement(0.0, 0.0).*ones(length(data.rate))
	for i=1:length(data.rate)
		data.β[i] = (data.d[i] - 2*data.df[i])^2/(data.d[i]^2 - (data.d[i] - 2*data.df[i])^2)
		lin_fit, data.a[i], data.b[i], mean_μR, std_μR, Tchar_β_sort, Tchar_sort, branch = linear_fit_β(data, i; α=0.05)
		data.ΔTelu[i] = measurement(lin_fit.param[1], standard_errors(lin_fit)[1])
		data.μelu[i] = measurement(mean_μR, std_μR)
		data.peaklist[i].Tchar = Tchar_sort
		data.peaklist[i].Tchar_β = Tchar_β_sort
		data.peaklist[i].branch = branch
	end
end

# ╔═╡ b743e2e3-3dea-4897-baab-b438d6137daa
function holdup_time(T, pin, pout, L, d, df, gas)
	# holdup time for conventional GC (no gradient)
	Tst = 273.15
	# viscosity
	if gas=="He"
		ηst = 18.63e-6
		ξ₀ = 0.6958
		ξ₁ = -0.0071
	elseif gas=="H2"
		ηst = 8.382e-6
		ξ₀ = 0.6892
		ξ₁ = 0.005
	elseif gas=="N2"
		ηst = 16.62e-6
		ξ₀ = 0.7665
		ξ₁ = -0.0378
	elseif gas=="Ar"
		ηst = 21.04e-6
		ξ₀ = 0.8131
		ξ₁ = -0.0426
	else
		error("Unknown selection of gas. Choose one of these: He, H2, N2 or Ar.")
	end
	η = ηst*(T/Tst)^(ξ₀ + ξ₁*(T-Tst)/Tst)
	id = d-2*df
	tM = 128/3*L^2/id^2*η*(pin^3-pout^3)/(pin^2-pout^2)^2
	return tM
end

# ╔═╡ ef26356e-8f53-4eae-afa1-c6938500421b
function plot_chromatogram(peaklist, T_itp, Δtsteps, L; labelpeaks=true, tempprog=false, neg=false)
    if peaklist.tR[end]<sum(Δtsteps)
        tend = sum(Δtsteps)
    else
        tend = 1.1*peaklist.tR[end]
    end
    if neg==true
        sgn = -1.0
    else
        sgn = 1.0
    end
    gauss(t, tR, τR) = 1/sqrt(2*π*τR^2)*exp(-(t-tR)^2/(2*τR^2))
    gauss_sum(t) = sum(gauss.(t, peaklist.tR, peaklist.τR))
    time = 0:tend/10000:tend
    chrom = gauss_sum.(time)
    pc = plot(time, sgn.*chrom,
                title="Chromatogram",
                legend=false,
                grid=false,
                xlims=(0,tend))
    # peak annotations
    if labelpeaks==true
        peaklabel = Array{Plots.PlotText}(undef, size(peaklist)[1])
        for i=1:size(peaklist)[1]
            peaklabel[i] = Plots.text(string(i), 8)
        end
        scatter!(peaklist.tR, sgn.*gauss_sum.(peaklist.tR), 
                    series_annotations=peaklabel,
                    markersize=0)
    end
    if tempprog==true
        plot!(twinx(),time, [T_itp(0,time).-273.15,T_itp(L,time).-273.15],
                    label=["z=0" "z=L"],
                    legend=:topright,
                    ylims=(minimum([minimum(T_itp(0,time)),minimum(T_itp(L,time))])-273.15-5,1.2*(maximum([maximum(T_itp(0,time)),maximum(T_itp(L,time))])-273.15)),
                    xlims=(0,tend),
                    ylabel="temperature in °C",
                    color=[:red :blue],
                    grid=true)
    end
    xlabel!("time in s")

    return pc, time, chrom
end

# ╔═╡ 9f6eefa6-b4e0-426e-a81f-42342ab05f3c
begin
	md"""
	### Survey of single simulation
	
	#### Select Simulation
	
	$(@bind ii NumberField(1:size(data)[1])) of $(size(data)[1])
	
	"""
end

# ╔═╡ f7c4089c-bc10-4acc-b92c-63b74107f029
data.peaklist[ii]

# ╔═╡ fe226f03-840d-4325-bee0-c437b6a2da5f
data.sp[ii]

# ╔═╡ 30fc1f13-891c-499b-a285-f195f6849cc2
data.rate[ii]

# ╔═╡ 3a9a2dca-c20f-4d72-b1f8-b9d2123b7efd
data.β[ii]

# ╔═╡ 35d21a42-f8a3-4b55-accd-6596253a617b
data.rate[ii]*data.tM_ref[ii]/30/60

# ╔═╡ a814e1bc-3795-4c26-9db1-e953c4cf45d1
data.ΔTelu[ii]

# ╔═╡ 5a6b9a80-4be0-4a4f-84f4-c5e60e1e662a
begin
	Δtsteps = [0.0, (data.Tend[ii]-data.Tinit[ii])./data.rate[ii].*60.0]
	pc = plot_chromatogram(data.peaklist[ii], data.par[ii].prog.T_itp[ii], Δtsteps, data.L[ii])[1]
	
	p1 = scatter(data.peaklist[ii].Tchar.-273.15,data.peaklist[ii].TR,
					legend=:bottomright,
					legendfontsize=6,
					xlabel="Tchar in °C",
					ylabel="elution temperature Telu in °C",
					title="heating rate r=$(data.rate[ii])°C/min, pin=$(data.pin[ii])kPa(g), $(data.sp[ii]), β=$(round(data.β[ii],digits=2))",
					label="data org.",
					markersize=4
					)
	scatter!(data.peaklist[ii].Tchar_β.-273.15,data.peaklist[ii].TR,
					label="data β corr."
					)
	scatter!(data.peaklist[ii].Tchar_β[data.a[ii]:data.b[ii]].-273.15,data.peaklist[ii].TR[data.a[ii]:data.b[ii]], label="used data", marker=:xcross, markersize=4)
	plot!(p1, data.peaklist[ii].Tchar_β.-273.15, data.peaklist[ii].Tchar_β.-273.15, label="TR=Tchar")
	plot!(p1, data.peaklist[ii].Tchar_β[data.a[ii]:data.b[ii]].-273.15, Measurements.value.(data.peaklist[ii].Tchar_β[data.a[ii]:data.b[ii]].-273.15 .+ data.ΔTelu[ii]),label="ΔTelu=$(data.ΔTelu[ii])°C")
	
	p2 = scatter(data.peaklist[ii].Tchar_β.-273.15,data.peaklist[ii].kR,
					legend=false,
					xlabel="Tchar_β in °C",
					ylabel="retention factor kR",
					title="heating rate r=$(data.rate[ii])°C/min, pin=$(data.pin[ii])kPa(g), $(data.sp[ii]), β=$(round(data.β[ii],digits=2))",
					label="data"
					)
	
	p3 = scatter(data.peaklist[ii].Tchar_β.-273.15,1.0 ./ (1.0 .+ data.peaklist[ii].kR),
					legend=false,
					xlabel="Tchar_β in °C",
					ylabel="mobility μR",
					title="heating rate r=$(data.rate[ii])°C/min, pin=$(data.pin[ii])kPa(g), $(data.sp[ii]), β=$(round(data.β[ii],digits=2))",
					label="data"
					)
	plot!(p3, [data.peaklist[ii].Tchar_β[data.a[ii]], data.peaklist[ii].Tchar_β[data.b[ii]]].-273.15, [data.μelu[ii],data.μelu[ii]])
	
	pp = plot(pc, p1, p2, p3, lay=(2,2), title="")
	md"""
	
	heating rate r = $(data.rate[ii]) °C/min, pin = $(data.pin[ii]) kPa(g), pout = $(data.pout[ii]) kPa(a), $(data.sp[ii]) , β = $(round(data.β[ii],digits=2))
	
	$(embed_display(pp))
	
	$(embed_display(data.peaklist[ii]))
	"""
end

# ╔═╡ b1435002-009b-4508-92de-df4458783396
begin
	T_ref = 150+273.15
	data.tM_ref = zeros(Float64, length(data.rate))
	for i=1:length(data.rate)
		data.tM_ref[i] = holdup_time(T_ref, data.pin[i]*1000+101300, data.pout[i]*1000, data.L[i], data.d[i], data.df[i], data.gas[i])
	end
	#α = 0.05
	#model(x,p) = p[1] .+ p[2].*sqrt.(x)
	#p0 = [-100.0, 100.0]
	#fit = Array{Any}(undef, length(unique(data.β)))
	#x1 = minimum(data.rate.*data.tM_ref./(60*30))
	#x2 = maximum(data.rate.*data.tM_ref./(60*30))
	#x = x1:(x2-x1)/1000:x2
	p_ΔTelu_TM_ref_3 = scatter(legend=:bottomright,
				xlabel="dimless heating rate",
				ylabel="ΔTec in °C"
				)
	#p1 = Array{Measurement{Float64}}(undef, length(unique(data.β)))
	#p2 = Array{Measurement{Float64}}(undef, length(unique(data.β)))
	#R2 = Array{Float64}(undef, length(unique(data.β)))
	for k=1:length(unique(data.β))
		data_filter_k = filter([:β] => (β) -> β==unique(data.β)[k],
										data)
		scatter!(p_ΔTelu_TM_ref_3, data_filter_k.rate.*data_filter_k.tM_ref./(60*30), data_filter_k.ΔTelu,
					label="β=$(round(unique(data.β)[k], digits=1))"
					)
		#fit[k] = curve_fit(model, data_filter_k.rate.*data_filter_k.tM_ref./(60*30), Measurements.value.(data_filter_k.ΔTelu), p0)
		#p1[k] = measurement(fit[k].param[1], standard_errors(fit[k])[1])
		#p2[k] = measurement(fit[k].param[2], standard_errors(fit[k])[2])
		#p1bottom = confidence_interval(fit[k], α)[1][1]
		#p1top = confidence_interval(fit[k], α)[1][2]
		#p2bottom = confidence_interval(fit[k], α)[2][1]
		#p2top = confidence_interval(fit[k], α)[2][2]
		#R2[k] = Rsquare(fit[k], model(x, fit[k].param))
		#plot!(p_ΔTelu_TM_ref_3, 
		#			x, model(x, fit[k].param), 
		#			label="($(round(p1[k], digits=2)))°C + ($(round(p2[k], digits=2))) √x, R²=$(round(R2[k], digits=4))", 
		#			c=:red,
		#			ribbon=(model(x,[p1top,p2top]).-model(x,[p1bottom,p2bottom]))./2,
		#			fillalpha=0.25,
		#			legend=:outertop,
		#			size=(700,700)
		#			)
	end
	
	#plot_resid = plot(xlabel="dimless heating rate", ylabel="residuen")
	#for k=1:length(unique(data.β))
	#	data_filter_k = filter([:β] => (β) -> β==unique(data.β)[k],
	#									data)
	#	scatter!(plot_resid, data_filter_k.rate.*data_filter_k.tM_ref./(60*30), fit[k].resid)
	#end
	
	md"""
	
	$(embed_display(p_ΔTelu_TM_ref_3))
	
	
	"""
	#$(embed_display(plot_resid))
end

# ╔═╡ 480d9ab5-808c-43f8-83de-76ea40d804bb
function Rsquare(fit, y)
	sstot = sum((y.-mean(y)).^2)
	ssres = sum(fit.resid.^2)
	R2 = 1-ssres/sstot
	return R2
end

# ╔═╡ 1a8f46de-ac7f-46fd-9d7b-c381f7480a43
begin
	# fit of models to the data, differentiated for β and stationary phase
	model1(x,p) = p[1] .+ p[2].*sqrt.(x)
	p01 = [-100.0, 100.0]
	model2(x,p) = p[1] .+ p[2].*x.^p[3]
	p02 = [-100.0, 100.0, 0.5]
	model3(x,p) = p[1] .+ p[2].*sqrt.(x) .+ p[3]./x
	p03 = [-100.0, 100.0, 0.5]
	
	model = [model1, model2, model3]
	p0 = [p01, p02, p03]
	
	phases = Array{String}(undef, length(unique(data.sp)))
	for i=1:length(unique(data.sp))
		phases[i] = unique(data.sp)[i]
	end
	βs = unique(data.β)
	
	fit = Array{Any}(undef, length(βs), length(phases), 3)
	R² = Array{Float64}(undef, length(βs), length(phases), 3)
	
	for j=1:length(phases)
		for i=1:length(βs)
			data_filter = filter([:β, :sp] => (β, sp) -> β==βs[i] && sp==phases[j],
										data)
			for k=1:3
				fit[i,j,k] = curve_fit(model[k], data_filter.rate.*data_filter.tM_ref./(60*30), Measurements.value.(data_filter.ΔTelu), p0[k])
				R²[i,j,k] = Rsquare(fit[i,j,k], Measurements.value.(data_filter.ΔTelu))
			end
		end
	end
end

# ╔═╡ b6f64134-b483-453b-95fd-eb751a9d1e78
R²

# ╔═╡ bec53d21-4b00-48c2-b083-92c4e1e559f4
fit[1,1,1]

# ╔═╡ b34adec6-190d-48e5-ba0f-06d96f588eed
md"""

#### Select model for plots

model $(@bind select_model Select(["1", "2", "3"]))
"""

# ╔═╡ c597ad34-1da7-4654-a7f5-1b319df43678
md"""

#### Select stationary phase

$(@bind select_phase Select(phases))

"""

# ╔═╡ a752d071-5131-447f-bf02-a58bf50c57ca
begin
	p_ΔTec_phases = scatter(legend=:bottomright,
				xlabel="dimless heating rate",
				ylabel="ΔTec in °C"
				)
	for i=1:length(unique(data.β))
		data_filter = filter([:β, :sp] => (β, sp) -> β==unique(data.β)[i] && sp==select_phase, data)
		scatter!(p_ΔTec_phases, data_filter.rate.*data_filter.tM_ref./(60*30), data_filter.ΔTelu,
					label="β=$(round(unique(data.β)[i], digits=1))",
					title="$(select_phase)"
					)
		x1 = minimum(data_filter.rate.*data_filter.tM_ref./(60*30))
		x2 = maximum(data_filter.rate.*data_filter.tM_ref./(60*30))
		x = x1:(x2-x1)/100:x2
		j = findfirst(phases.==select_phase)
		k = parse(Int, select_model)
		plot!(p_ΔTec_phases, x, model[k](x, fit[i, j, k].param), label="R²=$(round(R²[i,j,k],digits=4))")
	end
	md"""
	
	$(embed_display(p_ΔTec_phases))
	
	
	"""
end

# ╔═╡ e35e94c9-125d-4860-8e67-1b7fad8c264a
md"""
Jumps in the graph, e.g. for Rxi17SilMS from dimless heating rates of 0.35 to 0.39, are a result of not ideal filtering of the solutes for the estimation of ΔTelu. For r<0.35 only some late eluted solutes (PAHs) are used, less retained solutes have a higher deviation of the mobility (which is used to filter the solutes) from these late eluted solutes. Therfore the estimation of ΔTelu is to low. -> **Is there another way to select the right solutes for the ΔTelu estimation?**

Other jumps, espeacialy for low β (61.8) could have other reasons (e.g. problem with calculation of ΔTchar_β, due to usage of the 0-branch of the Lambert-W-function) -> this has to be tested.
"""

# ╔═╡ 52361dbe-5713-4d39-9917-9d6028ac2c5c
begin
	phase_ratios = Array{String}(undef, length(unique(data.β)))
	for i=1:length(unique(data.β))
		phase_ratios[i] = string(unique(data.β)[i])
	end
	md"""

	#### Select phase ratio
	
	$(@bind select_β Select(phase_ratios))
	
	"""
end

# ╔═╡ 8bab719c-0770-4318-ad35-4e4aa32f9f98
begin
	p_ΔTec_β = scatter(legend=:bottomright,
				xlabel="dimless heating rate",
				ylabel="ΔTec in °C"
				)
	for j=1:length(phases)
		data_filter = filter([:sp, :β] => (sp, β) -> sp==phases[j] && β==parse(Float64,select_β), data)
		scatter!(p_ΔTec_β, data_filter.rate.*data_filter.tM_ref./(60*30), data_filter.ΔTelu,
					label=phases[j],
					title="β = $(select_β)"
					)
		x1 = minimum(data_filter.rate.*data_filter.tM_ref./(60*30))
		x2 = maximum(data_filter.rate.*data_filter.tM_ref./(60*30))
		x = x1:(x2-x1)/100:x2
		i = findfirst(βs.==parse(Float64,select_β))
		k = parse(Int, select_model)
		plot!(p_ΔTec_β, x, model[k](x, fit[i, j, k].param), label="R²=$(round(R²[i,j,k],digits=4))")
	end
	md"""
	
	$(embed_display(p_ΔTec_β))
	
	
	"""
end

# ╔═╡ 592614b9-c4e9-4977-9028-ce7ee4dd48ba
begin
	p_ΔTec_minβ = scatter(legend=:bottomright,
				xlabel="dimless heating rate",
				ylabel="ΔTec in °C"
				)
	for k=1:length(phases)
		data_filter_k = filter([:sp, :β] => (sp, β) -> sp==phases[k] && β==minimum(parse.(Float64,phase_ratios)), data)
		scatter!(p_ΔTec_minβ, data_filter_k.rate.*data_filter_k.tM_ref./(60*30), data_filter_k.ΔTelu,
					label=phases[k],
					title="β = $(minimum(parse.(Float64,phase_ratios)))"
					)
	end
	
	p_ΔTec_medianβ = scatter(legend=:bottomright,
				xlabel="dimless heating rate",
				ylabel="ΔTec in °C"
				)
	for k=1:length(phases)
		data_filter_k = filter([:sp, :β] => (sp, β) -> sp==phases[k] && β==median(parse.(Float64,phase_ratios)), data)
		scatter!(p_ΔTec_medianβ, data_filter_k.rate.*data_filter_k.tM_ref./(60*30), data_filter_k.ΔTelu,
					label=phases[k],
					title="β = $(median(parse.(Float64,phase_ratios)))"
					)
	end
	
	p_ΔTec_maximumβ = scatter(legend=:bottomright,
				xlabel="dimless heating rate",
				ylabel="ΔTec in °C"
				)
	for k=1:length(phases)
		data_filter_k = filter([:sp, :β] => (sp, β) -> sp==phases[k] && β==maximum(parse.(Float64,phase_ratios)), data)
		scatter!(p_ΔTec_maximumβ, data_filter_k.rate.*data_filter_k.tM_ref./(60*30), data_filter_k.ΔTelu,
					label=phases[k],
					title="β = $(maximum(parse.(Float64,phase_ratios)))"
					)
	end
	md"""
	#### Differences for stationary phases
	
	For small values of the phase ratio β (β=61.75, thick film df=1µm for d=0.25mm) the data for different stationary phases diverge.
	
	$(embed_display(p_ΔTec_minβ))
	
	For medium values of the phase ratio β (β=249.25, standard film thicknes df=0.25µm for d=0.25mm) the data for all stationary phases is similar.
	
	$(embed_display(p_ΔTec_medianβ))
	
	For high values of the phase ratio β (β=249.25, thin film thicknes df=0.05µm for d=0.25mm) the data for all stationary phases is similar.
	
	$(embed_display(p_ΔTec_maximumβ))
	
	Overall the data for FS5ms has the highest difference to the data of the other stationary phases. The quality of the data for the FS5ms stationary phase could be the reason.
	"""
end

# ╔═╡ 9349968b-0e9c-4ba6-a49f-546121d1b049
begin
	# model 1:
	p_m1_p1 = plot(xlabel="ln(β)", ylabel="parameter p1", legend=:bottomright, legendfontsize=6)
	p_m1_p2 = plot(xlabel="ln(β)", ylabel="parameter p2", legend=:topright, legendfontsize=6)
	p_m2_p1 = plot(xlabel="ln(β)", ylabel="parameter p1", legend=:bottomright, legendfontsize=6)
	p_m2_p2 = plot(xlabel="ln(β)", ylabel="parameter p2", legend=:topright, legendfontsize=6)
	p_m2_p3 = plot(xlabel="ln(β)", ylabel="parameter p3", legend=:topleft, legendfontsize=6)
	p_m3_p1 = plot(xlabel="ln(β)", ylabel="parameter p1", legend=:bottomright, legendfontsize=6)
	p_m3_p2 = plot(xlabel="ln(β)", ylabel="parameter p2", legend=:topright, legendfontsize=6)
	p_m3_p3 = plot(xlabel="ln(β)", ylabel="parameter p3", legend=:topleft, legendfontsize=6)
	par = Array{Measurements.Measurement{Float64}}(undef, length(βs), length(phases), 3, 3)
	for j=1:length(phases)
		for i=1:length(βs)
			for k=1:3
				for l=1:3
					if l==3 && k==1
						v = NaN
						u = NaN
					else
						v = fit[i,j,k].param[l]
						u = standard_errors(fit[i,j,k])[l]
					end
				par[i,j,k,l] = measurement(v, u) 
				end
			end
		end
	end
	
	for j=1:length(phases)
		scatter!(p_m1_p1, log.(βs), par[:,j,1,1], label=phases[j])
		scatter!(p_m1_p2, log.(βs), par[:,j,1,2], label=phases[j])
		scatter!(p_m2_p1, log.(βs), par[:,j,2,1], label=phases[j])
		scatter!(p_m2_p2, log.(βs), par[:,j,2,2], label=phases[j])
		scatter!(p_m2_p3, log.(βs), par[:,j,2,3], label=phases[j])
		scatter!(p_m3_p1, log.(βs), par[:,j,3,1], label=phases[j])
		scatter!(p_m3_p2, log.(βs), par[:,j,3,2], label=phases[j])
		scatter!(p_m3_p3, log.(βs), par[:,j,3,3], label=phases[j])
	end
	p_m1=plot(p_m1_p1, p_m1_p2) 
	p_m2=plot(p_m2_p1, p_m2_p2, p_m2_p3) 
	p_m3=plot(p_m3_p1, p_m3_p2, p_m3_p3) 
	md"""
	#### Plot of the model parameters over phase ratio β
	
	##### Model 1: $y = p_1 + p_2\sqrt{x}$
	
	$(embed_display(p_m1))
	"""
end

# ╔═╡ a0cc03bc-56dc-490a-8832-58a1075adf2c
md"""
##### Model 2: $y = p_1 + p_2 x^{p_3}$
	
$(embed_display(p_m2))
"""

# ╔═╡ 4945cb84-5fc4-46d4-820c-d6adf54c6a6c
md"""
##### Model 3: $y = p_1 + p_2\sqrt{x} + p_3/x$
	
$(embed_display(p_m3))
"""

# ╔═╡ Cell order:
# ╠═31dc694a-9cfb-495b-baa4-c30e85634744
# ╠═42fd0ee9-cc2f-4889-b5d6-4431829aa5a7
# ╠═7e517e2a-7ab7-4e2f-8c35-3f6b37c70e33
# ╠═2f202810-dc23-4760-8bb4-f00e4fd96f23
# ╠═b7f9ec8d-c637-401d-a47a-6d3c30568abd
# ╠═2bb5eb7c-7507-43ab-8c5e-0086f0e21a1c
# ╠═42e7ca1a-b7d8-4cd8-9998-02f27cedbc59
# ╠═dc2c53dd-89cc-4e93-bd58-e20b2b5750e6
# ╠═7c1b33e7-4e58-401a-b914-3c2004ce6ca8
# ╠═e8d0e6dd-d0a7-4692-baeb-71712022f23a
# ╠═bf31b8ee-6222-4596-b177-3db9a91ac677
# ╠═a57cfb19-7718-4336-bead-fff25b19287c
# ╠═5f27ef08-9237-4379-b15a-65cfa311881b
# ╠═f7c4089c-bc10-4acc-b92c-63b74107f029
# ╠═fe226f03-840d-4325-bee0-c437b6a2da5f
# ╠═30fc1f13-891c-499b-a285-f195f6849cc2
# ╠═3a9a2dca-c20f-4d72-b1f8-b9d2123b7efd
# ╠═35d21a42-f8a3-4b55-accd-6596253a617b
# ╠═a814e1bc-3795-4c26-9db1-e953c4cf45d1
# ╠═2ccd805a-fcf3-44c6-a79e-6db1ff6e08d0
# ╠═388a7f39-8faa-4cfa-9a38-360a9e124cde
# ╠═e30cb724-1f9c-448c-9b59-977e4e1da4fc
# ╠═b743e2e3-3dea-4897-baab-b438d6137daa
# ╠═ef26356e-8f53-4eae-afa1-c6938500421b
# ╠═9f6eefa6-b4e0-426e-a81f-42342ab05f3c
# ╠═5a6b9a80-4be0-4a4f-84f4-c5e60e1e662a
# ╠═b1435002-009b-4508-92de-df4458783396
# ╠═480d9ab5-808c-43f8-83de-76ea40d804bb
# ╠═1a8f46de-ac7f-46fd-9d7b-c381f7480a43
# ╠═b6f64134-b483-453b-95fd-eb751a9d1e78
# ╠═bec53d21-4b00-48c2-b083-92c4e1e559f4
# ╟─b34adec6-190d-48e5-ba0f-06d96f588eed
# ╟─c597ad34-1da7-4654-a7f5-1b319df43678
# ╟─a752d071-5131-447f-bf02-a58bf50c57ca
# ╟─e35e94c9-125d-4860-8e67-1b7fad8c264a
# ╟─52361dbe-5713-4d39-9917-9d6028ac2c5c
# ╟─8bab719c-0770-4318-ad35-4e4aa32f9f98
# ╟─592614b9-c4e9-4977-9028-ce7ee4dd48ba
# ╠═9349968b-0e9c-4ba6-a49f-546121d1b049
# ╠═a0cc03bc-56dc-490a-8832-58a1075adf2c
# ╠═4945cb84-5fc4-46d4-820c-d6adf54c6a6c
