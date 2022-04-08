### A Pluto.jl notebook ###
# v0.17.1

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

# ╔═╡ 63a79620-7509-11eb-2afd-9de5f7f069d9
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	using DrWatson, DataFrames, CSV, Interpolations, Plots, LsqFit, Measurements, PlutoUI
	using GasChromatographySimulator
	plotly()
	TableOfContents()
end

# ╔═╡ 76c45ca2-7509-11eb-2ec4-89575f222aad
begin
	d(x) = GasChromatographySimulator.gradient(x, 2.44e-4)
	df(x) = GasChromatographySimulator.gradient(x, 2.44e-7)
	gf(x) = GasChromatographySimulator.gradient(x, [[0.0, 0.0]])
	parameters = Dict(
					:L => 29.43,
					:d => d,
					:a_d => [[2.44e-4]],
					:df => df,
					:a_df => [[2.44e-7]],
					:sp => ["SLB5ms","SPB50","Wax","FS5ms","Rxi17SilMS","DB5ms","Rxi5MS"],
					:gas => "He",
					:rate => [1.0, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0, 30.0, 50.0], # °C/min,
					:Tsteps => [[30.0, 300.0]],
					#const pressure is needed for constant difference between elution temp. and Tchar
					:pin => [100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0], # in kPa(g)
					:poutsteps => [[1.0, 1.0 ] .* 101300],
					:gf => gf, # no gradient
					:a_gf => [[0.0, 0.0]],
					:solute_db_path => datadir("Databases"),
					:solute_db => "Database_append.csv",
					:abstol => 1e-6,
					:reltol => 1e-3
					)
	Ndict = dict_list_count(parameters)
	dicts = dict_list(parameters)
	md"""
	Define parameters of simulations.
	"""
end

# ╔═╡ d243af7d-c7c7-49e3-810b-4ddbd7e4c8d2
function all_solutes(sp, db)
	db_filter = filter([:Phase] => x -> x==sp, db)
	solutes = db_filter.Name
	return solutes
end

# ╔═╡ 44bcda5c-751f-11eb-32f7-53d72a578fe4
begin
	function makesim(dict::Dict)
		@unpack L, d, a_d, df, a_df, sp, gas, rate, Tsteps, pin, poutsteps, gf, a_gf, solute_db_path, solute_db, abstol, reltol = dict
		Tinit = Tsteps[1]
		Tend = Tsteps[2]
		Δtsteps = [0.0, (Tend-Tinit)./rate.*60.0]
		pinsteps = [pin, pin].*1000.0 .+ 101300.0
		sys = GasChromatographySimulator.constructor_System(L,a_d[1],a_df[1],sp,gas)
		prog = GasChromatographySimulator.constructor_Program(Δtsteps,
																Tsteps,
																pinsteps,
																poutsteps,
																sys.L)
		# here a determination of all available solutes for sp is needed!
		# or a version of load_solute_database(), which loads the data for all available solutes of a sp
		db = DataFrame(CSV.File(string(solute_db_path,"/",solute_db), header=1, silencewarnings=true))
		solutes = all_solutes(sp, db)
		sub = GasChromatographySimulator.load_solute_database(db, sp, gas, solutes, zeros(length(solutes)), zeros(length(solutes)))
		opt = GasChromatographySimulator.Options(OwrenZen5(), abstol, reltol, "inlet", true)
		par = GasChromatographySimulator.Parameters(sys,prog,sub,opt)

		# run the simulation 
		sol = GasChromatographySimulator.solve_system_multithreads(par; ng=true) 
		pl = GasChromatographySimulator.peaklist(sol, par)
		fulldict = copy(dict)
		fulldict[:peaklist] = pl
		return fulldict
	end
	md""" 
	Definition of the simulation function.
	"""
end

# ╔═╡ 6dffdff4-751f-11eb-269d-2981558563dc
begin
	f = Array{Any}(undef, length(dicts))
	for (i, d) in enumerate(dicts)
		f[i] = makesim(d)
	end
	data = DataFrame(f)
	md"""
	Run the simulations.
	"""
end

# ╔═╡ dc5213b0-75d2-11eb-28b7-0b9af9b307d0
begin
	Nsim = 0
	for i=1:Ndict
		Nsim = Nsim + length(data.peaklist[i].Name)
	end
end

# ╔═╡ 519525a4-75d3-11eb-35be-990ffbd63057
md"""
	Number of total simulations (number of solutes simulated with $(Ndict) different settings): $(Nsim)
	"""

# ╔═╡ 1ad8b6dc-7521-11eb-2660-e17f3c9a8e83
begin	
	function extract_Tchar_sort(data, i)
		# extract the Tchar-values and sort them according to retention time
		# load database and extract Tchar (and name)
		db = DataFrame(CSV.File(string(data.solute_db_path[i],"/",data.solute_db[i]), header=1, silencewarnings=true))
		solutes = all_solutes(data.sp[i], db)
		sub = GasChromatographySimulator.load_solute_database(db, data.sp[i], data.gas[i], solutes, zeros(length(solutes)), zeros(length(solutes)))
		# sort Tchar according to TR
		Tchar_sort = Array{Float64}(undef, length(sub))
		for j=1:length(sub)
			Tchar_sort[findfirst(sub[j].name.==data.peaklist[i].Name)] = sub[j].Tchar
		end
		return Tchar_sort
	end
	
	function plot_TR_over_Tchar(data, i)
		# plot TR over Tchar
		# load database and extract Tchar (and name)
		Tchar_sort = extract_Tchar_sort(data, i)
		# plot TR over Tchar
		p1 = scatter(Tchar_sort.-273.15,data.peaklist[i].TR,
					legend=:bottomright,
					xlabel="Tchar in °C",
					ylabel="elution temperature Telu in °C",
					title="heating rate r=$(data.rate[i])°C/min, pin=$(data.pin[i])kPa(g), $(data.sp[i])",
					label="data"
					)
		plot!(Tchar_sort.-273.15, Tchar_sort.-273.15, label="Telu=Tchar")
		# plot difference TR-Tchar over Tchar
		#p2 = scatter(Tchar_sort.-273.15, f[:peaklist].TR .- (Tchar_sort.-273.15))
		return p1
	end
	
	function plot_μR_over_Tchar(data, i)
		# plot TR over Tchar
		# load database and extract Tchar (and name)
		Tchar_sort = extract_Tchar_sort(data, i)
		# plot TR over Tchar
		p1 = scatter(Tchar_sort.-273.15,1.0 ./ (1.0 .+ data.peaklist[i].kR),
					legend=:bottomright,
					xlabel="Tchar in °C",
					ylabel="mobility μR",
					title="heating rate r=$(data.rate[i])°C/min, pin=$(data.pin[i])kPa(g), $(data.sp[i])",
					label="data"
					)
		return p1
	end
	
	function linear_fit(data, i)
		# linear fit of TR over Tchar with slope = 1
		# -> intercept = mean difference TR-Tchar
		Tchar_sort = extract_Tchar_sort(data, i)
		@. lin_model(x,p) = x + p[1]
		p0 = [0.0]

		# filter for tR<tProg, for tR>tProg elution temperature is constant
		tProg = (data.Tsteps[i][2]-data.Tsteps[i][1])/data.rate[i]*60 # duration of heating ramp in s
		b = findfirst(tProg.<=data.peaklist[i].tR) # last element with tR<tProg
		if isa(b, Nothing)
			b = length(Tchar_sort)
		else
			b = b-1
		end

		# filter for unretained solutes, use kR or μR to identify them as outlyer
		# iterativly fit a zero-slope, reduce the first index until relative error of the intersept is <1%
		@. zero_model(x,p) = 0*x + p[1]
		pp0 = [0.25]
		zero_fit = curve_fit(zero_model, Tchar_sort.-273.15, 1.0 ./ (1.0 .+ data.peaklist[i].kR), pp0)
		a = 1
		if standard_errors(zero_fit)[1]/zero_fit.param[1] < 0.01
			a = 1
		else
			for ii=2:b
				zero_fit = curve_fit(zero_model, Tchar_sort[ii:b].-273.15, 1.0 ./ (1.0 .+ data.peaklist[i].kR[ii:b]), pp0)
				if standard_errors(zero_fit)[1]/zero_fit.param[1] < 0.01
					a = ii
					break
				end
			end
			#a = a + 1
			zero_fit = curve_fit(zero_model, Tchar_sort[a:b].-273.15, 1.0 ./ (1.0 .+ data.peaklist[i].kR[a:b]), pp0)
		end

		lin_fit = curve_fit(lin_model, Tchar_sort[a:b].-273.15, data.peaklist[i].TR[a:b], p0)
		return lin_fit, a, b, zero_fit
	end
	
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
	md""" 
	Definition of some additional functions.
	"""
end

# ╔═╡ bf66134e-7520-11eb-08cc-b32ca71e220b
begin
	T_ref = 150+273.15
	data.tM_ref = zeros(Float64, length(data.rate))
	data.a = zeros(Int64, length(data.rate)) # first index of allowed solutes
	data.b = zeros(Int64, length(data.rate)) # last index of allowed solutes
	data.ΔTelu = measurement(0.0, 0.0).*ones(length(data.rate))
	data.μelu = measurement(0.0, 0.0).*ones(length(data.rate))
	for i=1:length(data.rate)
		data.tM_ref[i] = holdup_time(T_ref, data.pin[i]*1000+101300, data.poutsteps[i][1], data.L[i], data.a_d[i][1], data.a_df[i][1], data.gas[i])
		lin_fit, data.a[i], data.b[i], zero_fit = linear_fit(data,i)
		data.ΔTelu[i] = measurement(lin_fit.param[1], standard_errors(lin_fit)[1])
		data.μelu[i] = measurement(zero_fit.param[1], standard_errors(zero_fit)[1])
	end
	md"""
	Evaluation of data.
	"""
end

# ╔═╡ 055c9c1c-7522-11eb-1335-79d5591172fa
begin
	md"""
	## Plot ``T_{elu}`` over ``T_{char}`` 
	
	for simulation number: $(@bind ii Slider(1:size(data)[1], show_value=true))
	"""
end

# ╔═╡ c9bcaddc-7523-11eb-179b-d11cfc574b0e
begin
	pTR = plot_TR_over_Tchar(data, ii)
	Tchar_sort = extract_Tchar_sort(data, ii)
	plot!(pTR, Tchar_sort[data.a[ii]:data.b[ii]].-273.15, Measurements.value.((Tchar_sort[data.a[ii]:data.b[ii]].-273.15) .+ data.ΔTelu[ii]),
        label="ΔTelu=$(data.ΔTelu[ii])°C")
end

# ╔═╡ ec126bf0-7524-11eb-3a08-d70fd0a3fb33
begin
	# plot with Interpolations
	p_ΔTelu_TM_ref_3 = plot(legend=:bottomright,
				xlabel="heating rate in °C/tM_ref",
				ylabel="ΔTelu in °C"
				)
	for j=1:length(unique(data.sp))
		data_filter_j = filter(:sp => sp -> sp==unique(data.sp)[j], data)
		data_filter_j.rate_dimless = data_filter_j.rate.*data_filter_j.tM_ref./60
		sort!(data_filter_j, :rate_dimless)
		itp = LinearInterpolation((data_filter_j.rate_dimless, ), data_filter_j.ΔTelu, extrapolation_bc=Flat())
		x = data_filter_j.rate_dimless[1]:(data_filter_j.rate_dimless[end]-data_filter_j.rate_dimless[1])/100:data_filter_j.rate_dimless[end]
		# plot
		plot!(p_ΔTelu_TM_ref_3, x, Measurements.value.(itp(x)), label="$(unique(data.sp)[j])")
		# estimate rate_dimless, where ΔTelu = 0°C
		itp_inv = LinearInterpolation((Measurements.value.(data_filter_j.ΔTelu), ), data_filter_j.rate_dimless, extrapolation_bc=Flat())
		scatter!(p_ΔTelu_TM_ref_3, [itp_inv(0), itp_inv(0)], [0,0], label="r=$(round(itp_inv(0),digits=2))°C/tM_ref", markershape=:+)
		# estimate ΔTelu at rate_dimless=10
		scatter!(p_ΔTelu_TM_ref_3, [10, 10], [itp(10), itp(10)], label="ΔTelu=$(itp(10))°C", markershape=:x)
	end
	md"""
	## Plot ``\Delta T_{elu}`` over the normalized heating rate for different stationary phases.

	$(embed_display(p_ΔTelu_TM_ref_3))
	"""
end

# ╔═╡ 9ea9e414-75d4-11eb-11b4-4bcec0435589
begin
	# plot with Interpolations
	p_ΔTelu_TM_ref_all = plot(legend=:bottomright,
				xlabel="heating rate in °C/tM_ref",
				ylabel="ΔTelu in °C"
				)
	datacopy = copy(data)
	datacopy.rate_dimless = datacopy.rate.*datacopy.tM_ref./60
	sort!(datacopy, :rate_dimless)
	itp = LinearInterpolation((datacopy.rate_dimless, ), datacopy.ΔTelu, extrapolation_bc=Flat())
	x = datacopy.rate_dimless[1]:(datacopy.rate_dimless[end]-datacopy.rate_dimless[1])/100:datacopy.rate_dimless[end]
	# plot
	plot!(p_ΔTelu_TM_ref_all, x, Measurements.value.(itp(x)), label="all")
	# estimate rate_dimless, where ΔTelu = 0°C
	sort!(datacopy, :ΔTelu)
	itp_inv = LinearInterpolation((Measurements.value.(datacopy.ΔTelu), ), datacopy.rate_dimless, extrapolation_bc=Flat())
	scatter!(p_ΔTelu_TM_ref_all, [itp_inv(0), itp_inv(0)], [0,0], label="r=$(round(itp_inv(0),digits=2))°C/tM_ref", markershape=:+)
	# estimate ΔTelu at rate_dimless=10
	scatter!(p_ΔTelu_TM_ref_all, [10, 10], [itp(10), itp(10)], label="ΔTelu=$(itp(10))°C", markershape=:x)
	md"""
	## Plot ``\Delta T_{elu}`` over the normalized heating rate for all stationary phases.

	$(embed_display(p_ΔTelu_TM_ref_all))
	"""
end

# ╔═╡ Cell order:
# ╠═63a79620-7509-11eb-2afd-9de5f7f069d9
# ╠═76c45ca2-7509-11eb-2ec4-89575f222aad
# ╠═d243af7d-c7c7-49e3-810b-4ddbd7e4c8d2
# ╠═44bcda5c-751f-11eb-32f7-53d72a578fe4
# ╠═6dffdff4-751f-11eb-269d-2981558563dc
# ╟─dc5213b0-75d2-11eb-28b7-0b9af9b307d0
# ╟─519525a4-75d3-11eb-35be-990ffbd63057
# ╠═1ad8b6dc-7521-11eb-2660-e17f3c9a8e83
# ╟─bf66134e-7520-11eb-08cc-b32ca71e220b
# ╟─055c9c1c-7522-11eb-1335-79d5591172fa
# ╟─c9bcaddc-7523-11eb-179b-d11cfc574b0e
# ╟─ec126bf0-7524-11eb-3a08-d70fd0a3fb33
# ╟─9ea9e414-75d4-11eb-11b4-4bcec0435589
