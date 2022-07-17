### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ 5090d60e-01e0-11ed-30aa-a1eb336e73bf
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir());
    Pkg.add([
		Pkg.PackageSpec(name="CSV"),
		Pkg.PackageSpec(name="DataFrames"),
		Pkg.PackageSpec(name="DrWatson", version="2.9"),
		Pkg.PackageSpec(name="GasChromatographySimulator", version="0.3.6"),
		#Pkg.PackageSpec(name="Interpolations"),
		#Pkg.PackageSpec(name="LambertW"),
		Pkg.PackageSpec(name="LsqFit", version="0.12"),
		Pkg.PackageSpec(name="Measurements", version="2.7"),
        Pkg.PackageSpec(name="Plots", version="1"),
        Pkg.PackageSpec(name="PlutoUI", version="0.7.39"),
		Pkg.PackageSpec(name="Statistics"),
	]);
    using CSV, DataFrames, DrWatson, GasChromatographySimulator, LsqFit, Measurements, Plots, PlutoUI, Statistics;
end;

# ╔═╡ 68f02f89-3af7-453c-bc62-59630bf2b44f
begin
	plotly()
	#gr()
	TableOfContents()
end

# ╔═╡ ffab5c90-0ad2-46dc-953a-f995c04588fa
md"""
# Simulation of Correlation of Elution Temperature ``T_{elu}``, characteristic Temperature ``T_{char}`` depending on Heating Rate ``R_T``

In this notebook simulations of different programs (heating rates) are made and the correlation between ``T_{char}`` and ``T_{elu}`` are investigated.
"""

# ╔═╡ e09c4c7f-f7e6-4fcb-8e02-832e9c6146cf
md"""
## Settings for Simulation
"""

# ╔═╡ d39c4cde-27cc-417b-8c5e-df53cad4695c
begin
	# program: single ramp, no holding times
	# variation of parameters:
	# - over stationary phases:
	sp = ["SLB5ms","SPB50","Wax","Rxi17SilMS","DB5ms","Rxi5ms"]
	# - over selected key
	var_key = :rate
	variation = [0.001, 0.01, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
	var_label = "heating rate in °C/min"
	
	parameters = Dict(
					:L => 30.0, #m
	                :d => 2.5e-4, #m
	                :df => 0.25e-6, #m
	                :sp => sp,
	                :gas => "He",
	                :rate => variation, #°C/min
	                :Tinit => 40.0, #°C
					:Tend => 1000.0, #°C
	                :constMode => "Pressure", #"Pressure" or "Flow"
					##:constMode => "Flow",
	                :Fpin => 150.0, # in kPa(g) or mL/min
					##:Fpin => 1.0, # in mL/min
	                :pout => 0.0, # in kPa(a)
	                :solute_db_path => "/Users/janleppert/Documents/GitHub/ThermodynamicData/Databases/",
	                :solute_db => "newformat_Kcentric_2022-05-25T12:49:45.405.csv",
	                :abstol => 1e-6,
	                :reltol => 1e-3
	)
	dicts = dict_list(parameters);
end

# ╔═╡ dc2331f0-01f7-451e-83fb-6893b2c836b7
Ndict = dict_list_count(parameters)

# ╔═╡ edce9bf3-fd9f-4d3e-82f6-82cc9334422f
md"""
## Run simulation
"""

# ╔═╡ 25d3dae2-18d8-4131-a17c-a974eb8c00a2
md"""
## Model for the fit

### Four parameters

``T_{elu}(T_{char}) = \frac{1}{b} \ln{\left(\exp{\left(b\left(T_0 + m T_{char}\right)\right)}+\exp{\left(b T_1\right)}\right)}``

Parameters: 
- ``b`` [°C⁻¹] ... curvature/bend of the transition area between the two lines
- ``T_0`` [°C] ... intercept of the increasing line
- ``m`` ... slope of the increasing line (should be 1.0)
- ``T_1`` [°C] ... value/level of the constant line
"""

# ╔═╡ c2bf501e-8aaa-48b2-a4b4-a54288d09250
@. model(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*x)) + exp(p[1]*p[4]))

# ╔═╡ fc54a58c-01a3-4c25-86bb-d0484f3c59d8
md"""
### Constant slope ``m=1``

``T_{elu}(T_{char}) = \frac{1}{b} \ln{\left(\exp{\left(b\left(T_0 + T_{char}\right)\right)}+\exp{\left(b T_1\right)}\right)}``

Parameters: 
- ``b`` [°C⁻¹] ... curvature/bend of the transition area between the two lines
- ``T_0`` [°C] ... intercept of the increasing line
- ``T_1`` [°C] ... value/level of the constant line
"""

# ╔═╡ 0267473e-1c39-4a65-bac0-5a66588f3b0f
@. model_m1(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + 1.0*x)) + exp(p[1]*p[3]))

# ╔═╡ 28809c3a-852f-4700-8d7a-4c6cdf7af169
md"""
### Constant slope ``m=1`` and bend ``b=0.05°C⁻¹``

``T_{elu}(T_{char}) = \frac{1}{0.05} \ln{\left(\exp{\left(0.05\left(T_0 + T_{char}\right)\right)}+\exp{\left(0.05 T_1\right)}\right)}``

Parameters: 
- ``T_0`` [°C] ... intercept of the increasing line
- ``T_1`` [°C] ... value/level of the constant line
"""

# ╔═╡ 75e43083-c9ec-4371-8b17-6e0be57db885
@. model_m1_b005(x,p) = 1/0.05 * log(exp(0.05*(p[1] + 1.0*x)) + exp(0.05*p[2]))

# ╔═╡ 1a1bd6d9-7ad2-4cf5-9dcd-90accf339a7e
md"""
# Plot ``T_{elu}`` over ``T_{char}``
"""

# ╔═╡ 3a3ebd7b-8bad-47da-b5c2-9abebca8b23f
md"""
### Fit four parameters
"""

# ╔═╡ a205bcfc-ef60-4630-b48b-07d530976189
md"""
### Fit ``m=1``
"""

# ╔═╡ 5d314d84-ba24-4e9f-ba1e-e01ddf265ae7
md"""
### Fit ``m=1`` and ``b=0.05``°C⁻¹
"""

# ╔═╡ a2f9c71b-938f-4852-81d9-d33744d4db41
md"""
### ``R^2``
"""

# ╔═╡ 1add6cbe-9357-44a7-90d2-24604996ce71
md"""
## Combine all stationary phases

Plot of the data ``T_{elu}(T_{char})`` for all solute -- stationary phase combinations suggest an independent behaviour of the correlation between ``T_{elu}`` and ``T_{char}``.

#### Select variied value
	
$(@bind select_var Select(variation)) $(var_label)
"""

# ╔═╡ 7bb24245-2fe0-40e8-8f8b-27215c6a164d
md"""
### Fit four parameters
"""

# ╔═╡ 834896db-9162-4bca-b03b-bb2b74558395
md"""
### Fit ``m=1``
"""

# ╔═╡ ce00e529-943f-4b67-b423-986772469c52
md"""
### Fit ``m=1`` and ``b=0.05``°C⁻¹
"""

# ╔═╡ cab9ed7e-fd71-4ba8-ad8c-8ce16bb478bb
@. model_pot(x,p) = p[1] + p[2]*x^p[3]

# ╔═╡ f7e23edf-28b6-4a25-85df-e02585ca771a
@. model_log(x,p) = p[1] + p[2]*log(x)

# ╔═╡ 0865c20f-9df8-4330-80e2-8c75db19695e
@. model_log_lin(x,p) = p[1] + p[2]*log(x) + p[3]*x

# ╔═╡ c4e1f6c7-f17d-4a73-a81d-69b6f679cc06
@. model_log_lin_fix(x,p) = -273.15 + p[1]*log(x) + p[2]*x

# ╔═╡ d707240f-da1b-40cd-9f94-4be00ed625ee
@. model_pot_fix(x,p) = 40.0 + p[1]*sqrt(x)

# ╔═╡ d572f262-bb09-4485-97b9-0b4e0f680134
md"""
### Observation
Dependency of the four parameters from the heating rate suggest a simplified model with two constant parameters (``b, m``), one squarroot depending parameter $T_0$ and one linear depending parameter ``T_1``.
"""

# ╔═╡ 091828ef-c187-4705-ae17-5c6b32941654
md"""
### ``R^2``
"""

# ╔═╡ 95f36cbb-4677-4721-aa1d-44abf164bbe7
md"""
### MR-Fit for seven parameters
"""

# ╔═╡ 7fdd98f0-2633-452d-b98c-5423864d4a60
md"""
## Multivariate regression

Multivariate regression of ``T_{elu}(T_{char}, R_t) with the model:

``T_{elu}(T_{char}, R_T) = \frac{1}{b} \ln{\left(\exp{\left(b\left(T_0(R_T) + m T_{char}\right)\right)}+\exp{\left(b T_1(R_T)\right)}\right)}``

with ``b``and ``m`` constant and ``T_0`` and ``T_1`` depending on ``R_T``:

``T_0(R_T) = A_2 + B_2 R_T^{C_2}``

``T_1(R_T) = A_3 + B_3 R_T``
"""

# ╔═╡ d5dc0059-85fd-4fea-8a90-ab26d484fbec
model_T₀(x,A2,B2,C2) = A2 + B2*x^C2

# ╔═╡ 065c6b7a-1292-4005-b33a-502967040d51
model_T₁(x,A3,B3) = A3 + B3*x

# ╔═╡ 69b14db0-fc69-462d-8bd3-1c5780839a38
@. model_mr7(x,p) = 1/p[1] * log(exp(p[1]*(model_T₀.(x[:,2],p[3], p[4], p[5]) + p[2]*x[:,1])) + exp(p[1]*model_T₁.(x[:,2], p[6], p[7])))

# ╔═╡ 488478cd-54b5-44ec-ab08-94a72b523701
md"""
Fix values:
- ``m=1.0``
- ``C_2=0.5`` (sqrt-function)
- ``A_2=-273.15``
"""

# ╔═╡ d35a323c-6b95-4439-9144-66243ceec5c1
md"""
### MR-Fit for four parameters
"""

# ╔═╡ 9802ddab-9dce-4413-b3e6-696ccf99516a
@. model_mr4(x,p) = 1/p[1] * log(exp(p[1]*(model_T₀.(x[:,2],-273.15, p[2], 1/3) + 1.0*x[:,1])) + exp(p[1]*model_T₁.(x[:,2], p[3], p[4])))

# ╔═╡ 7214dd27-114c-4430-8e1c-456b8cb226b2
md"""
# End
"""

# ╔═╡ 76c3e597-c430-4b75-9e5c-f9878f89d9d4
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

# ╔═╡ b6f902b3-cf37-44ba-9be5-88d48f70e75f
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
    fulldict[:peaklist] = pl
    fulldict[:par] = par
	fulldict[:tMref] = tMref
    return fulldict
end

# ╔═╡ 7936ff29-e873-40b2-8b8a-bb1a1335cb6a
function dict_array_to_dataframe(f)
	df = DataFrame(f[1])
	for i=2:length(f)
		push!(df, f[i])
	end
	return df
end

# ╔═╡ ce7d30d7-0bcb-4219-a8c1-cf19d1403b6a
begin
	# run the simulation
	f = Array{Any}(undef, length(dicts))
	for (i, d) in enumerate(dicts)
	    f[i] = makesim(d)
	end
	data = dict_array_to_dataframe(f)
end

# ╔═╡ 7f9aca9e-dd4d-4dd0-ae3e-2eed9342bd88
begin
	md"""
	## Survey of single simulation
	
	#### Select Simulation
	
	$(@bind ii NumberField(1:size(data)[1])) of $(size(data)[1])
	
	"""
end

# ╔═╡ ad811db6-4678-4d1c-bd44-ebde12f693de
begin
	sort!(data, [:sp, :rate])
	TeluTchar = Array{DataFrame}(undef, length(data.rate))
	for i=1:length(data.rate)
		Tchars = Array{Float64}(undef, length(data.par[i].sub))
		Telus = Array{Float64}(undef, length(data.par[i].sub))
		Names = Array{String}(undef, length(data.par[i].sub))
		kelus = Array{Float64}(undef, length(data.par[i].sub))
		for j=1:length(data.par[i].sub)
			jj = findfirst(data.par[i].sub[j].name.==data.peaklist[i].Name)
			Tchars[j] = data.par[i].sub[j].Tchar - 273.15
			Telus[j] = data.peaklist[i].TR[jj]
			Names[j] = data.par[i].sub[j].name
			kelus[j] = data.peaklist[i].kR[jj]
		end
		TeluTchar[i] = filter!([:Telu] => x -> !isnan(x), sort!(DataFrame(Name=Names, Tchar=Tchars, Telu=Telus, kelu=kelus), [:Tchar]))  # sort Tchar and filter out NaN values
	end
end

# ╔═╡ 8c9b564b-9ebe-4914-a378-68cd4ca029cf
begin
	TeluTchar_allsp = Array{DataFrame}(undef, length(variation))
	for i=1:length(variation)
		data_f = filter([var_key] => x -> x == variation[i], data)
		Tchars = Float64[]
		Telus = Float64[]
		Names = String[]
		kelus = Float64[]
		sps = String[]
		for j=1:length(data_f.par)
			for k=1:length(data_f.par[j].sub)
				kk = findfirst(data_f.par[j].sub[k].name.==data_f.peaklist[j].Name)
				push!(Tchars, data_f.par[j].sub[k].Tchar - 273.15)
				push!(Telus, data_f.peaklist[j].TR[kk])
				push!(Names, data_f.par[j].sub[k].name)
				push!(kelus, data_f.peaklist[j].kR[kk])
				push!(sps, data_f.sp[j])
			end
		end
		TeluTchar_allsp[i] = filter!([:Telu] => x -> !isnan(x), sort!(DataFrame(Name=Names, sp=sps, Tchar=Tchars, Telu=Telus, kelu=kelus), [:Tchar])) # sort Tchar and filter out NaN values
	end
end

# ╔═╡ 31edd1fc-aa56-438e-99e8-eca3d1daf32e
function fitting(TeluTchar, model, p0)
	N = size(TeluTchar)[1]
	fits = Array{LsqFit.LsqFitResult}(undef, N)
	for i=1:N
		xdata = TeluTchar[i].Tchar
		ydata = TeluTchar[i].Telu
		fits[i] = curve_fit(model, xdata, ydata, p0)
	end
	return fits
end

# ╔═╡ 3e51232b-74fa-420b-a7ed-11c24792d277
function Rsquare(fit, y)
	sstot = sum((y.-mean(y)).^2)
	ssres = sum(fit.resid.^2)
	R2 = 1-ssres/sstot
	return R2
end

# ╔═╡ fcaa1bf3-d4cf-4488-b422-55dbcaf53894
function collect_param(fits, TeluTchar)
	N = size(fits)[1]
	M = length(fits[1].param)
	param = Array{Measurement{Float64}}(undef, N, M)
	R² = Array{Float64}(undef, N)
	for i=1:N
		sigmas = NaN.*ones(M)
		try
			sigmas = stderror(fits[i])
		catch
			sigmas = NaN.*ones(M)
		end
		for j=1:M
			param[i,j] = fits[i].param[j] ± sigmas[j]
		end
		R²[i] = Rsquare(fits[i], TeluTchar[i].Telu)
	end
	return param, R²
end

# ╔═╡ 39f72f16-adf8-4b63-a721-c9bd43c66aea
begin #fit for the single simulations
	fits = fitting(TeluTchar, model, [0.01, 0.0, 1.0, 40.0])
	param, R² = collect_param(fits, TeluTchar)
end

# ╔═╡ 20243f4c-4735-47e5-b8de-836944c0ff5a
begin
	p_b = scatter(1:size(param)[1], param[:,1], xlabel="data number", ylabel="b in °C⁻¹", legend=false)
	p_T₀ = scatter(1:size(param)[1], param[:,2], xlabel="data number", ylabel="T₀ in °C", legend=false)
	p_m = scatter(1:size(param)[1], param[:,3], xlabel="data number", ylabel="m", legend=false)
	p_T₁ = scatter(1:size(param)[1], param[:,4], xlabel="data number", ylabel="T₁ in °C", legend=false)
	p_param = plot(p_b, p_T₀, p_m, p_T₁, lay=(2,2), title="")
end

# ╔═╡ 4131601f-2b69-4c89-be2f-aad158b79cf8
begin #fit for the single simulations
	fits_m1 = fitting(TeluTchar, model_m1, [0.01, 0.0, 40.0])
	param_m1, R²_m1 = collect_param(fits_m1, TeluTchar)
end

# ╔═╡ 66d1577c-f7f7-451f-b3e2-17ddea89b38e
begin
	p_b_m1 = scatter(1:size(param_m1)[1], param_m1[:,1], xlabel="data number", ylabel="b in °C⁻¹", legend=false)
	p_T₀_m1 = scatter(1:size(param_m1)[1], param_m1[:,2], xlabel="data number", ylabel="T₀ in °C", legend=false)
	p_m_m1 = scatter(1:size(param_m1)[1], ones(size(param_m1)[1]), xlabel="data number", ylabel="m", legend=false)
	p_T₁_m1 = scatter(1:size(param_m1)[1], param_m1[:,3], xlabel="data number", ylabel="T₁ in °C", legend=false)
	p_param_m1 = plot(p_b_m1, p_T₀_m1, p_m_m1, p_T₁_m1, lay=(2,2), title="")
end

# ╔═╡ 10a9bd4d-412a-4d62-9a50-292baad6bd99
begin #fit for the single simulations
	fits_m1_b005 = fitting(TeluTchar, model_m1_b005, [0.0, 40.0])
	param_m1_b005, R²_m1_b005 = collect_param(fits_m1_b005, TeluTchar)
end

# ╔═╡ 2d9d2faa-78f6-4d1e-b633-5abc82fb66e0
begin
	tlim = (0.0, (data.Tend[ii]-data.Tinit[ii])./data.rate[ii].*60.0)
	pc = GasChromatographySimulator.plot_chromatogram(data.peaklist[ii], tlim)[1]
	
	p1 = scatter(TeluTchar[ii].Tchar,TeluTchar[ii].Telu,
					legend=:bottomright,
					legendfontsize=6,
					xlabel="Tchar in °C",
					ylabel="elution temperature Telu in °C",
					label="data",
					markersize=4
					)
	
	plot!(p1, TeluTchar[ii].Tchar, TeluTchar[ii].Tchar, label="Telu=Tchar", c=:black, s=:dash)
	plot!(p1, TeluTchar[ii].Tchar, model(TeluTchar[ii].Tchar, fits[ii].param), c=:red, w=3, label="fit")
	plot!(p1, TeluTchar[ii].Tchar, model_m1(TeluTchar[ii].Tchar, fits_m1[ii].param), c=:orange, w=2, label="fit m=1")
	plot!(p1, TeluTchar[ii].Tchar, model_m1_b005(TeluTchar[ii].Tchar, fits_m1_b005[ii].param), c=:darkgreen, w=2, label="fit m=1 b=0.05")
	
	p2 = scatter(TeluTchar[ii].Tchar,TeluTchar[ii].kelu,
					legend=false,
					xlabel="Tchar in °C",
					ylabel="retention factor at elution kelu",
					label="data"
					)
	
	p3 = scatter(TeluTchar[ii].Tchar,1.0 ./ (1.0 .+ TeluTchar[ii].kelu),
					legend=false,
					xlabel="Tchar in °C",
					ylabel="mobility at elution μR",
					label="data"
					)
	
	pp = plot(pc, p1, p2, p3, lay=(2,2), title="")
	md"""
	
	variied quantity: $(var_label) = $(data[ii, var_key])
	
	stationary phase : $(data.sp[ii])
	
	$(embed_display(pp))
	
	$(embed_display(TeluTchar[ii]))
	"""
end

# ╔═╡ 972b44fb-d35d-41cf-9219-64d13baf5416
begin
	p_b_m1_b005 = scatter(1:size(param_m1_b005)[1], 0.05.*ones(size(param_m1_b005)[1]), xlabel="data number", ylabel="b in °C⁻¹", legend=false)
	p_T₀_m1_b005 = scatter(1:size(param_m1_b005)[1], param_m1_b005[:,1], xlabel="data number", ylabel="T₀ in °C", legend=false)
	p_m_m1_b005 = scatter(1:size(param_m1_b005)[1], ones(size(param_m1_b005)[1]), xlabel="data number", ylabel="m", legend=false)
	p_T₁_m1_b005 = scatter(1:size(param_m1_b005)[1], param_m1_b005[:,2], xlabel="data number", ylabel="T₁ in °C", legend=false)
	p_param_m1_b005 = plot(p_b_m1_b005, p_T₀_m1_b005, p_m_m1_b005, p_T₁_m1_b005, lay=(2,2), title="")
end

# ╔═╡ f332180b-4b9a-4a1b-b13c-24a3c50e3461
begin
	p_R² = scatter(1:length(R²), R², xlabel="sim number", ylabel="R²", label="4 param")
	scatter!(p_R², 1:length(R²_m1), R²_m1, label="m=1")
	scatter!(p_R², 1:length(R²_m1_b005), R²_m1_b005, label="m=1, b=0.05")
	p_R²
end

# ╔═╡ a0c73457-bf0b-4553-ba1e-9b1fdff4eac2
begin #fit for all simulations of same heating rate
	fits_allsp = fitting(TeluTchar_allsp, model, [0.01, 0.0, 1.0, 40.0])
	param_allsp, R²_allsp = collect_param(fits_allsp, TeluTchar_allsp)
end

# ╔═╡ 340212ad-b8e5-428c-8757-8a79c61a5174
begin
	p_b_allsp = scatter(variation, param_allsp[:,1], xlabel=var_label, ylabel="b in °C⁻¹", legend=false)
	p_T₀_allsp = scatter(variation, param_allsp[:,2], xlabel=var_label, ylabel="T₀ in °C", legend=false)
	p_m_allsp = scatter(variation, param_allsp[:,3], xlabel=var_label, ylabel="m", legend=false)
	p_T₁_allsp = scatter(variation, param_allsp[:,4], xlabel=var_label, ylabel="T₁ in °C", legend=false)
	p_param_allsp = plot(p_b_allsp, p_T₀_allsp, p_m_allsp, p_T₁_allsp, lay=(2,2), title="")
end

# ╔═╡ 2f4a5e5c-76c7-4d4f-a38d-a509535c1cbc
begin #fit for all simulations of same heating rate
	fits_allsp_m1 = fitting(TeluTchar_allsp, model_m1, [0.01, 0.0, 40.0])
	param_allsp_m1, R²_allsp_m1 = collect_param(fits_allsp_m1, TeluTchar_allsp)
end

# ╔═╡ 26ddca0e-063e-4b78-a34e-d0d23b07a965
begin
	p_b_allsp_m1 = scatter(variation, param_allsp_m1[:,1], xlabel=var_label, ylabel="b in °C⁻¹", legend=false)
	p_T₀_allsp_m1 = scatter(variation, param_allsp_m1[:,2], xlabel=var_label, ylabel="T₀ in °C", legend=false)
	p_m_allsp_m1 = scatter(variation, ones(length(variation)), xlabel=var_label, ylabel="m", legend=false)
	p_T₁_allsp_m1 = scatter(variation, param_allsp_m1[:,3], xlabel=var_label, ylabel="T₁ in °C", legend=false)
	p_param_allsp_m1 = plot(p_b_allsp_m1, p_T₀_allsp_m1, p_m_allsp_m1, p_T₁_allsp_m1, lay=(2,2), title="")
end

# ╔═╡ 2efeae52-54f6-490b-80d6-fc7af9a1fbd8
begin #fit for all simulations of same heating rate
	fits_allsp_m1_b005 = fitting(TeluTchar_allsp, model_m1_b005, [0.0, 40.0])
	param_allsp_m1_b005, R²_allsp_m1_b005 = collect_param(fits_allsp_m1_b005, TeluTchar_allsp)
end

# ╔═╡ 952cc093-306a-460e-a421-2f9d3a3b6f32
begin
	ii_ = findfirst(select_var.==variation)
	#tlim_ = (0.0, (data.Tend[ii_rate]-data.Tinit[ii_rate])./data.rate[ii_rate].*60.0)
	#pc_ = GasChromatographySimulator.plot_chromatogram(data.peaklist[ii_rate], tlim_)[1]
	
	p1_ = scatter(TeluTchar_allsp[ii_].Tchar,TeluTchar_allsp[ii_].Telu,
					legend=:bottomright,
					legendfontsize=6,
					xlabel="Tchar in °C",
					ylabel="elution temperature Telu in °C",
					label="data",
					markersize=4
					)
	
	plot!(p1_, TeluTchar_allsp[ii_].Tchar, TeluTchar_allsp[ii_].Tchar, label="Telu=Tchar", c=:black, s=:dash)
	plot!(p1_, TeluTchar_allsp[ii_].Tchar, model(TeluTchar_allsp[ii_].Tchar, fits_allsp[ii_].param), c=:red, w=3, label="fit")
	plot!(p1_, TeluTchar_allsp[ii_].Tchar, model_m1(TeluTchar_allsp[ii_].Tchar, fits_allsp_m1[ii_].param), c=:orange, w=2, label="fit m=1")
	plot!(p1_, TeluTchar_allsp[ii_].Tchar, model_m1_b005(TeluTchar_allsp[ii_].Tchar, fits_allsp_m1_b005[ii_].param), c=:darkgreen, w=2, label="fit m=1 b=0.05")
	
	p2_ = scatter(TeluTchar_allsp[ii_].Tchar,TeluTchar_allsp[ii_].kelu,
					legend=false,
					xlabel="Tchar in °C",
					ylabel="retention factor at elution kelu",
					label="data"
					)
	
	p3_ = scatter(TeluTchar_allsp[ii_].Tchar,1.0 ./ (1.0 .+ TeluTchar_allsp[ii_].kelu),
					legend=false,
					xlabel="Tchar in °C",
					ylabel="mobility at elution μR",
					label="data"
					)
	
	pp_ = plot(p1_, p2_, p3_, lay=(2,2), title="")
	md"""

	
	$(embed_display(pp_))

	$(embed_display(TeluTchar_allsp[ii_]))
	"""
end

# ╔═╡ cd7c4b87-7774-4ef9-a090-a26940b13b5f
begin
	plotly()
	p_b_allsp_m1_b005 = scatter(variation, 0.05.*ones(length(variation)), xlabel=var_label, ylabel="b in °C⁻¹", legend=false)
	p_T₀_allsp_m1_b005 = scatter(variation, param_allsp_m1_b005[:,1], xlabel=var_label, ylabel="T₀ in °C", legend=false)
	p_m_allsp_m1_b005 = scatter(variation, ones(length(variation)), xlabel=var_label, ylabel="m", legend=false)
	p_T₁_allsp_m1_b005 = scatter(variation, param_allsp_m1_b005[:,2], xlabel=var_label, ylabel="T₁ in °C", legend=false)
	p_param_allsp_m1_b005 = plot(p_b_allsp_m1_b005, p_T₀_allsp_m1_b005, p_m_allsp_m1_b005, p_T₁_allsp_m1_b005, lay=(2,2), title="")
end

# ╔═╡ ca8b1960-d376-4010-94d7-55d6742b508b
Measurements.value.(param_allsp_m1_b005[:,1])

# ╔═╡ 00cb1a96-54fb-487e-8682-82777afee0ec
begin
	gr()
	fit_pot_T0 = curve_fit(model_pot, variation, Measurements.value.(param_allsp_m1_b005[:,1]), [1.0, 1.0, 1.0])
	R²_pot_T0 = Rsquare(fit_pot_T0, Measurements.value.(param_allsp_m1_b005[:,1]))

	fit_log_T0 = curve_fit(model_log, variation, Measurements.value.(param_allsp_m1_b005[:,1]), [1.0, 1.0])
	R²_log_T0 = Rsquare(fit_log_T0, Measurements.value.(param_allsp_m1_b005[:,1]))

	fit_log_lin_T0 = curve_fit(model_log_lin, variation, Measurements.value.(param_allsp_m1_b005[:,1]), [1.0, 1.0, 1.0])
	R²_log_lin_T0 = Rsquare(fit_log_lin_T0, Measurements.value.(param_allsp_m1_b005[:,1]))
	
	fit_log_lin_fix_T0 = curve_fit(model_log_lin_fix, variation, Measurements.value.(param_allsp_m1_b005[:,1]), [1.0, 1.0])
	R²_log_lin_fix_T0 = Rsquare(fit_log_lin_fix_T0, Measurements.value.(param_allsp_m1_b005[:,1]))
	
	p_fit_T0 = scatter(variation, param_allsp_m1_b005[:,1], xlabel=var_label, ylabel="T₀ in °C", legend=:bottomright)
	plot!(p_fit_T0, variation, model_pot(variation, fit_pot_T0.param), label="pot, R²=$(round(R²_pot_T0; digits=4))")
	plot!(p_fit_T0, variation, model_log(variation, fit_log_T0.param), label="log, R²=$(round(R²_log_T0; digits=4))")
	plot!(p_fit_T0, variation, model_log_lin(variation, fit_log_lin_T0.param), label="log_lin, R²=$(round(R²_log_lin_T0; digits=4))")
	plot!(p_fit_T0, variation, model_log_lin_fix(variation, fit_log_lin_fix_T0.param), label="log_lin_fix, R²=$(round(R²_log_lin_fix_T0; digits=4))")
end	

# ╔═╡ d6390804-daa5-41ae-a2ac-0b7aa67d50ad
fit_pot_T0.param

# ╔═╡ 89042ae0-9425-4675-884a-8e631d7b1994
fit_log_T0.param

# ╔═╡ 700fd250-cf92-4c55-a947-a94eea8074ff
fit_log_lin_T0.param

# ╔═╡ 6fb9f24c-bc74-4a5e-81e7-4963c6ba9257
fit_log_lin_fix_T0.param

# ╔═╡ d29cf2e2-3106-4c3f-9cfa-207e754d1a3c
begin
	plotly()
	fit_pot_T1 = curve_fit(model_pot, variation, Measurements.value.(param_allsp_m1_b005[:,2]), [1.0, 1.0, 1.0])
	R²_pot_T1 = Rsquare(fit_pot_T1, Measurements.value.(param_allsp_m1_b005[:,2]))

	fit_log_T1 = curve_fit(model_log, variation, Measurements.value.(param_allsp_m1_b005[:,2]), [1.0, 1.0])
	R²_log_T1 = Rsquare(fit_log_T1, Measurements.value.(param_allsp_m1_b005[:,2]))

	fit_log_lin_T1 = curve_fit(model_log_lin, variation, Measurements.value.(param_allsp_m1_b005[:,2]), [1.0, 1.0, 1.0])
	R²_log_lin_T1 = Rsquare(fit_log_lin_T0, Measurements.value.(param_allsp_m1_b005[:,2]))

	fit_pot_fix_T1 = curve_fit(model_pot_fix, variation, Measurements.value.(param_allsp_m1_b005[:,2]), [1.0])
	R²_pot_fix_T1 = Rsquare(fit_pot_fix_T1, Measurements.value.(param_allsp_m1_b005[:,2]))
	
	p_fit_T1 = scatter(variation, param_allsp_m1_b005[:,2], xlabel=var_label, ylabel="T₀ in °C", legend=:bottomright)
	plot!(p_fit_T1, variation, model_pot(variation, fit_pot_T1.param), label="pot, R²=$(round(R²_pot_T1; digits=4))")
	plot!(p_fit_T1, variation, model_log(variation, fit_log_T1.param), label="log, R²=$(round(R²_log_T1; digits=4))")
	plot!(p_fit_T1, variation, model_log_lin(variation, fit_log_lin_T1.param), label="log_lin, R²=$(round(R²_log_lin_T1; digits=4))")
	plot!(p_fit_T1, variation, model_pot_fix(variation, fit_pot_fix_T1.param), label="pot_fix, R²=$(round(R²_pot_fix_T1; digits=4))")
end	

# ╔═╡ 9f9c3e04-c8e0-4203-b6ab-af824c3b6b7b
fit_pot_T1.param

# ╔═╡ f080f8d0-7b80-4544-bf33-94639f8cd2ca
fit_log_T1.param

# ╔═╡ d5359d49-84c6-48af-9327-c3b2a903fd83
fit_log_lin_T1.param

# ╔═╡ 3fd3c176-ebd4-4698-a89e-e9ede8c74786
fit_pot_fix_T1.param

# ╔═╡ 79862725-39a0-479f-a146-c1c2d7ac559f
scatter(log.(variation), param_allsp_m1_b005[:,1], xlabel=string("log ", var_label), ylabel="T₀ in °C", legend=false)

# ╔═╡ 90ecdb78-ed4e-4619-8bd2-be7bf78cb079
scatter(sqrt.(variation), param_allsp_m1_b005[:,1], xlabel=string("sqrt(", var_label,")"), ylabel="T₀ in °C", legend=false)

# ╔═╡ 56ff6f01-5bb2-4a54-833a-686381c743b1
scatter((variation).^(1/3), param_allsp_m1_b005[:,1], xlabel=string("(", var_label,")^1/3"), ylabel="T₀ in °C", legend=false)

# ╔═╡ 6c6e508a-2e9f-4510-8aa4-5b60adf38006
begin
	p_R²_allsp = scatter(variation, R²_allsp, xlabel=var_label, ylabel="R²", label="4 param")
	scatter!(p_R²_allsp, variation, R²_allsp_m1, label="m=1")
	scatter!(p_R²_allsp, variation, R²_allsp_m1_b005, label="m=1, b=0.05")
	p_R²_allsp
end

# ╔═╡ 0502f3d8-b6b5-43d3-970c-1baf50340daf
function TeluTcharVar(variation, var_key, data)
	TeluTchar_allsp_allvar = Array{DataFrame}(undef, length(variation))
	var = Float64[]
	Tchars = Float64[]
	Telus = Float64[]
	Names = String[]
	kelus = Float64[]
	sps = String[]
	for i=1:length(variation)
		data_f = filter([var_key] => x -> x == variation[i], data)
		for j=1:length(data_f.par)
			for k=1:length(data_f.par[j].sub)
				kk = findfirst(data_f.par[j].sub[k].name.==data_f.peaklist[j].Name)
				push!(var, variation[i])
				push!(Tchars, data_f.par[j].sub[k].Tchar - 273.15)
				push!(Telus, data_f.peaklist[j].TR[kk])
				push!(Names, data_f.par[j].sub[k].name)
				push!(kelus, data_f.peaklist[j].kR[kk])
				push!(sps, data_f.sp[j])
			end
		end
	end	
	TeluTchar_allsp_allvar = filter!([:Telu] => x -> !isnan(x), sort!(DataFrame(Name=Names, sp=sps, Tchar=Tchars, Telu=Telus, kelu=kelus, var=var), [:Tchar])) # sort Tchar and filter out NaN values
	return TeluTchar_allsp_allvar
end

# ╔═╡ 58d0d92a-e97f-430b-a96b-0b7829e370f3
Telu_Tchar_Var = TeluTcharVar(variation, var_key, data)

# ╔═╡ 33c8a24c-25d8-4b05-94b1-c6eb5a26c310
scatter(Telu_Tchar_Var.Tchar, Telu_Tchar_Var.var, Telu_Tchar_Var.Telu, xlabel="Tchar in °C", ylabel=var_label, zlabel="elution temperature Telu in °C", legend=false)

# ╔═╡ 3a4f771e-24ed-47a2-b8d5-7efc2d0e4a3f
function fitting_mr(Telu_Tchar_Var, model, p0)
	xdata = [Telu_Tchar_Var.Tchar Telu_Tchar_Var.var]
	ydata = Telu_Tchar_Var.Telu
	fits = curve_fit(model, xdata, ydata,p0)

	M = length(fits.param)
	param = Array{Measurement{Float64}}(undef, M)
	sigmas = NaN.*ones(M)
	try
		sigmas = stderror(fits)
	catch
		sigmas = NaN.*ones(M)
	end
	for j=1:M
		param[j] = fits.param[j] ± sigmas[j]
	end
	R² = Rsquare(fits, Telu_Tchar_Var.Telu)
	return fits, param, R²
end

# ╔═╡ 2b241902-d2e5-4d66-b466-0eb1ccd0e485
fit_mr7, param_mr7, R²_mr7 = fitting_mr(Telu_Tchar_Var, model_mr7, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

# ╔═╡ 3637562b-4524-4b23-a248-6c4d9adc544e
md"""
#### Result
``b = `` $(param_mr7[1])

``m = `` $(param_mr7[2])

``A_2 = `` $(param_mr7[3])

``B_2 = `` $(param_mr7[4])

``C_2 = `` $(param_mr7[5])

``A_3 = `` $(param_mr7[6])

``B_3 = `` $(param_mr7[7])

``R^2 = `` $(R²_mr7)
"""

# ╔═╡ 47f187ae-56d9-4b52-a66f-23840926ad10
begin
	p_mr7 = scatter(legend=:topleft, xlabel="Tchar in °C", ylabel="elution temperature Telu in °C")
	for i=1:length(variation)
		scatter!(p_mr7, TeluTchar_allsp[i].Tchar, TeluTchar_allsp[i].Telu, label=variation[i])
		plot!(p_mr7, TeluTchar_allsp[i].Tchar, model_mr7([TeluTchar_allsp[i].Tchar variation[i].*ones(length(TeluTchar_allsp[i].Tchar))], fit_mr7.param), c=i, w=2, label="")

		#plot!(p_mr5, TeluTchar_allsp[i].Tchar, fit_mr5.param[2].+fit_mr5.param[3].*TeluTchar_allsp[i].Tchar, c=i, style=:dash, label="")
		#plot!(p_mr5, TeluTchar_allsp[i].Tchar, (fit_mr5.param[4]+fit_mr5.param[5]*variation[i]).*ones(length(TeluTchar_allsp[i].Tchar)), c=i, style=:dash, label="")
	end
	p_mr7
end

# ╔═╡ e2134ac3-eb9f-4f11-a37b-209ef5d1de82
fit_mr4, param_mr4, R²_mr4 = fitting_mr(Telu_Tchar_Var, model_mr4, [1.0, 1.0, 1.0, 1.0])

# ╔═╡ 044924b3-89fa-4580-8062-cc22b8aa7815
begin
	p_mr4 = scatter(legend=:topleft, xlabel="Tchar in °C", ylabel="elution temperature Telu in °C")
	for i=1:length(variation)
		scatter!(p_mr4, TeluTchar_allsp[i].Tchar, TeluTchar_allsp[i].Telu, label=variation[i])
		plot!(p_mr4, TeluTchar_allsp[i].Tchar, model_mr4([TeluTchar_allsp[i].Tchar variation[i].*ones(length(TeluTchar_allsp[i].Tchar))], fit_mr4.param), c=i, w=2, label="")
		#plot!(p_mr3, TeluTchar_allsp[i].Tchar, fit_mr3.param[2].+1.0.*TeluTchar_allsp[i].Tchar, c=i, style=:dash, label="")
		#plot!(p_mr3, TeluTchar_allsp[i].Tchar, (fit_mr3.param[3]+1.0*variation[i]).*ones(length(TeluTchar_allsp[i].Tchar)), c=i, style=:dash, label="")
	end
	p_mr4
end

# ╔═╡ 4c9b6aca-7739-4f55-b7d9-569194f5dea8
md"""
#### Result
``b = `` $(param_mr4[1])

``m = `` 1.0

``A_2 = `` -273.15

``B_2 = `` $(param_mr4[2])

``C_2 = `` 0.5

``A_3 = `` $(param_mr4[3])

``B_3 = `` $(param_mr4[4])

``R^2 = `` $(R²_mr4)
"""

# ╔═╡ Cell order:
# ╠═68f02f89-3af7-453c-bc62-59630bf2b44f
# ╠═ffab5c90-0ad2-46dc-953a-f995c04588fa
# ╟─e09c4c7f-f7e6-4fcb-8e02-832e9c6146cf
# ╠═d39c4cde-27cc-417b-8c5e-df53cad4695c
# ╠═dc2331f0-01f7-451e-83fb-6893b2c836b7
# ╟─edce9bf3-fd9f-4d3e-82f6-82cc9334422f
# ╠═ce7d30d7-0bcb-4219-a8c1-cf19d1403b6a
# ╠═25d3dae2-18d8-4131-a17c-a974eb8c00a2
# ╠═c2bf501e-8aaa-48b2-a4b4-a54288d09250
# ╟─fc54a58c-01a3-4c25-86bb-d0484f3c59d8
# ╠═0267473e-1c39-4a65-bac0-5a66588f3b0f
# ╟─28809c3a-852f-4700-8d7a-4c6cdf7af169
# ╠═75e43083-c9ec-4371-8b17-6e0be57db885
# ╟─1a1bd6d9-7ad2-4cf5-9dcd-90accf339a7e
# ╟─7f9aca9e-dd4d-4dd0-ae3e-2eed9342bd88
# ╟─ad811db6-4678-4d1c-bd44-ebde12f693de
# ╟─2d9d2faa-78f6-4d1e-b633-5abc82fb66e0
# ╟─3a3ebd7b-8bad-47da-b5c2-9abebca8b23f
# ╟─39f72f16-adf8-4b63-a721-c9bd43c66aea
# ╟─20243f4c-4735-47e5-b8de-836944c0ff5a
# ╟─a205bcfc-ef60-4630-b48b-07d530976189
# ╟─4131601f-2b69-4c89-be2f-aad158b79cf8
# ╟─66d1577c-f7f7-451f-b3e2-17ddea89b38e
# ╟─5d314d84-ba24-4e9f-ba1e-e01ddf265ae7
# ╟─10a9bd4d-412a-4d62-9a50-292baad6bd99
# ╟─972b44fb-d35d-41cf-9219-64d13baf5416
# ╟─a2f9c71b-938f-4852-81d9-d33744d4db41
# ╟─f332180b-4b9a-4a1b-b13c-24a3c50e3461
# ╠═1add6cbe-9357-44a7-90d2-24604996ce71
# ╟─8c9b564b-9ebe-4914-a378-68cd4ca029cf
# ╠═952cc093-306a-460e-a421-2f9d3a3b6f32
# ╟─7bb24245-2fe0-40e8-8f8b-27215c6a164d
# ╠═a0c73457-bf0b-4553-ba1e-9b1fdff4eac2
# ╟─340212ad-b8e5-428c-8757-8a79c61a5174
# ╟─834896db-9162-4bca-b03b-bb2b74558395
# ╟─2f4a5e5c-76c7-4d4f-a38d-a509535c1cbc
# ╟─26ddca0e-063e-4b78-a34e-d0d23b07a965
# ╟─ce00e529-943f-4b67-b423-986772469c52
# ╟─2efeae52-54f6-490b-80d6-fc7af9a1fbd8
# ╠═cd7c4b87-7774-4ef9-a090-a26940b13b5f
# ╠═cab9ed7e-fd71-4ba8-ad8c-8ce16bb478bb
# ╠═f7e23edf-28b6-4a25-85df-e02585ca771a
# ╠═0865c20f-9df8-4330-80e2-8c75db19695e
# ╠═c4e1f6c7-f17d-4a73-a81d-69b6f679cc06
# ╠═d707240f-da1b-40cd-9f94-4be00ed625ee
# ╠═ca8b1960-d376-4010-94d7-55d6742b508b
# ╠═00cb1a96-54fb-487e-8682-82777afee0ec
# ╠═d6390804-daa5-41ae-a2ac-0b7aa67d50ad
# ╠═89042ae0-9425-4675-884a-8e631d7b1994
# ╠═700fd250-cf92-4c55-a947-a94eea8074ff
# ╠═6fb9f24c-bc74-4a5e-81e7-4963c6ba9257
# ╠═d29cf2e2-3106-4c3f-9cfa-207e754d1a3c
# ╠═9f9c3e04-c8e0-4203-b6ab-af824c3b6b7b
# ╠═f080f8d0-7b80-4544-bf33-94639f8cd2ca
# ╠═d5359d49-84c6-48af-9327-c3b2a903fd83
# ╠═3fd3c176-ebd4-4698-a89e-e9ede8c74786
# ╠═79862725-39a0-479f-a146-c1c2d7ac559f
# ╠═90ecdb78-ed4e-4619-8bd2-be7bf78cb079
# ╠═56ff6f01-5bb2-4a54-833a-686381c743b1
# ╠═d572f262-bb09-4485-97b9-0b4e0f680134
# ╟─091828ef-c187-4705-ae17-5c6b32941654
# ╠═6c6e508a-2e9f-4510-8aa4-5b60adf38006
# ╠═58d0d92a-e97f-430b-a96b-0b7829e370f3
# ╠═33c8a24c-25d8-4b05-94b1-c6eb5a26c310
# ╠═95f36cbb-4677-4721-aa1d-44abf164bbe7
# ╟─7fdd98f0-2633-452d-b98c-5423864d4a60
# ╠═d5dc0059-85fd-4fea-8a90-ab26d484fbec
# ╠═065c6b7a-1292-4005-b33a-502967040d51
# ╠═69b14db0-fc69-462d-8bd3-1c5780839a38
# ╠═2b241902-d2e5-4d66-b466-0eb1ccd0e485
# ╠═3637562b-4524-4b23-a248-6c4d9adc544e
# ╠═47f187ae-56d9-4b52-a66f-23840926ad10
# ╠═488478cd-54b5-44ec-ab08-94a72b523701
# ╠═d35a323c-6b95-4439-9144-66243ceec5c1
# ╠═9802ddab-9dce-4413-b3e6-696ccf99516a
# ╠═e2134ac3-eb9f-4f11-a37b-209ef5d1de82
# ╠═044924b3-89fa-4580-8062-cc22b8aa7815
# ╠═4c9b6aca-7739-4f55-b7d9-569194f5dea8
# ╟─7214dd27-114c-4430-8e1c-456b8cb226b2
# ╠═b6f902b3-cf37-44ba-9be5-88d48f70e75f
# ╠═76c3e597-c430-4b75-9e5c-f9878f89d9d4
# ╠═7936ff29-e873-40b2-8b8a-bb1a1335cb6a
# ╠═31edd1fc-aa56-438e-99e8-eca3d1daf32e
# ╠═fcaa1bf3-d4cf-4488-b422-55dbcaf53894
# ╠═3e51232b-74fa-420b-a7ed-11c24792d277
# ╠═0502f3d8-b6b5-43d3-970c-1baf50340daf
# ╠═3a4f771e-24ed-47a2-b8d5-7efc2d0e4a3f
# ╠═5090d60e-01e0-11ed-30aa-a1eb336e73bf
