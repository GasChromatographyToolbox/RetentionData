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

# ╔═╡ 7d9ae1fa-047c-11ed-0af1-c1af45ffc663
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

# ╔═╡ f8ad60a7-d51c-40b0-8f82-1fc78dcd4b54
TableOfContents()

# ╔═╡ 562d78f0-6e8d-4b44-9cb3-70acd226f439
md"""
# Evaluation ``T_{elu}(T_{char})``
"""

# ╔═╡ 2a8dd6b6-dc9c-453b-9652-2226fc1f39ff
md"""
## Load data
"""

# ╔═╡ e5f30f41-3205-4746-969f-3d7c106e11bf
begin
	df = collect_results("/Users/janleppert/Documents/GitHub/ThermodynamicData/data/sims/Telu_Tchar/")
	df.dimless_rate = df.rate.*df.tMref./30.0
	df
end

# ╔═╡ 650a0f9b-77c8-4918-adf5-28fae3c7dd35
scatter(df.rate, df.dimless_rate)

# ╔═╡ e069b363-fc28-4665-9d27-f95371433dd2
names(df)

# ╔═╡ 3a85a6d5-d20c-4245-b8f7-5779d57e6f04
sp = sort!(unique(df.sp))

# ╔═╡ 56822d93-0c5f-44a1-a9fe-f034b9561d28
Tinit = sort!(unique(df.Tinit))

# ╔═╡ 17b52208-736e-4b21-a584-c49a9849e1ec
RT = sort!(unique(df.rate))

# ╔═╡ b7d6a89b-6ce2-426e-90d2-4ea289a407b8
pin = sort!(unique(df.Fpin))

# ╔═╡ 825ca4e1-02b4-4415-9a92-44dd92e94487
rT = sort!(unique(df.dimless_rate))

# ╔═╡ 39655593-9b8c-444e-8f7d-48cd92d7b9d8


# ╔═╡ 686250d4-235d-4699-88cc-5aaf1273b319
md"""
## Plot ``T_{elu}(T_{char}, T_{init}, r_T)``

$(@bind select_i_rT Slider(1:length(rT)))
"""

# ╔═╡ 8f5bed6b-726e-4b19-bd20-ec41151ce83d
md"""
## Plot ``T_{elu}(T_{char}, T_{init}, R_T, p_{in})``

$(@bind select_i_RT Slider(1:length(RT)))

$(@bind select_i_pin Slider(1:length(pin)))
"""

# ╔═╡ 00350329-7af5-4723-b858-f412a414c67f
# [x] filter out the triglycerides -> Tchar>700K
# fit the model with four parameters to every set of single combinations of Tinit, RT, pin resp. Tinit and rT for all stationary phases together

# ╔═╡ f67161c6-dc5f-44a9-b4ae-937f122b76df
md"""
## Model for the fit

### Four parameters

``T_{elu}(T_{char}) = \frac{1}{b} \ln{\left(\exp{\left(b\left(T_0 + m T_{char}\right)\right)}+\exp{\left(b T_1\right)}\right)}``

Parameters: 
- ``T_1`` [°C] ... value/level of the constant line
- ``T_0`` [°C] ... intercept of the increasing line
- ``b`` [°C⁻¹] ... curvature/bend of the transition area between the two lines
- ``m`` ... slope of the increasing line (should be 1.0)

"""

# ╔═╡ a261cc22-30c0-4d07-b8de-758618ec608d
@. model4(x,p) = 1/p[3] * log(exp(p[3]*(p[1] + p[4]*x)) + exp(p[3]*p[2]))

# ╔═╡ 689207e3-04f1-4521-8f87-ee8ebe5b09e7
md"""
### Three parameters (``b=``const.)

``T_{elu}(T_{char}) = \frac{1}{b} \ln{\left(\exp{\left(b\left(T_0 + m T_{char}\right)\right)}+\exp{\left(b T_1\right)}\right)}``

Parameters: 
- ``T_0`` [°C] ... intercept of the increasing line
- ``T_1`` [°C] ... value/level of the constant line
- ``m`` ... slope of the increasing line (should be 1.0)
"""

# ╔═╡ ce4b4f7f-69f0-4622-81fe-7b3f5c8089de
bconst = 0.1

# ╔═╡ 7e51c8ab-2e1a-4b0f-8ea3-9ae9399eb82b
@. model_b(x,p) = 1/bconst * log(exp(bconst*(p[1] + p[3]*x)) + exp(bconst*p[2]))

# ╔═╡ 7a9f8b03-3cb3-40a2-9dd8-7f4a5b95a018
md"""
### Three parameters (``m=1``)

``T_{elu}(T_{char}) = \frac{1}{b} \ln{\left(\exp{\left(b\left(T_0 + T_{char}\right)\right)}+\exp{\left(b T_1\right)}\right)}``

Parameters: 
- ``T_0`` [°C] ... intercept of the increasing line
- ``T_1`` [°C] ... value/level of the constant line
- ``b`` [°C⁻¹] ... curvature/bend of the transition area between the two lines
"""

# ╔═╡ 2214b214-ff3b-448a-8f7e-a4a754f8a5a3
mconst = 1.0

# ╔═╡ d969f868-81e5-494b-a528-ca86e07f8b2b
@. model_m1(x,p) = 1/p[3] * log(exp(p[3]*(p[1] + mconst*x)) + exp(p[3]*p[2]))

# ╔═╡ a8d54a2b-f9ab-4f95-aa0e-12819cc35ea7
md"""
### Two parameters (``b=``const., ``m=1``)

``T_{elu}(T_{char}) = \frac{1}{b} \ln{\left(\exp{\left(b\left(T_0 + T_{char}\right)\right)}+\exp{\left(b T_1\right)}\right)}``

Parameters: 
- ``T_0`` [°C] ... intercept of the increasing line
- ``T_1`` [°C] ... value/level of the constant line
"""

# ╔═╡ 13f2b66b-05e3-4ee2-a01a-fbd47d57233b
@. model_b_m1(x,p) = 1/bconst * log(exp(bconst*(p[1] + mconst*x)) + exp(bconst*p[2]))

# ╔═╡ cd2185f3-9f37-4371-9b94-833df7dd1cbd
md"""
## Fits
"""

# ╔═╡ d8550c75-51e3-4d25-ad73-c4d72eac4267
md"""
### For var. ``R_T`` and ``p_{in}``

$(@bind select_model Select([model4, model_b, model_m1, model_b_m1]))
"""

# ╔═╡ ebf5ee8e-430e-4634-8543-fbeec58fb0f8
md"""
Select pin: $(@bind select_pin Select(pin; default=150.0))
"""

# ╔═╡ 8897fb75-8a19-48d6-8ec9-4e16fb36f9c2
md"""
### For ``r_T``

$(@bind select_model_ Select([model4, model_b, model_m1, model_b_m1]))
"""

# ╔═╡ 599d46f2-c2e5-415f-9a05-60696ae57a52
function plot_parameters_i(param)	
	return plot(plot(1:length(param.R²), param.R², ylabel="R²"),
		plot(1:length(param.T₀), Measurements.value.(param.T₀), ylabel="T₀"),
		plot(1:length(param.T₁), Measurements.value.(param.T₁), ylabel="T₁"),
		plot(1:length(param.b), Measurements.value.(param.b), ylabel="b"),
		plot(1:length(param.m), Measurements.value.(param.m), ylabel="m"),
		legend=false
	)
end

# ╔═╡ c13c9d9c-7b81-4a16-90dd-f56c69e15f23
function plot_parameters_Tinit(param)
	p_T₀ = plot(xlabel="Tinit in °C", ylabel="T₀")
	p_T₁ = plot(xlabel="Tinit in K", ylabel="T₁")
	p_b = plot(xlabel="Tinit in K", ylabel="b")
	p_m = plot(xlabel="Tinit in K", ylabel="m")
	if "RT" in names(param)
		rate = unique(param.RT)
		rate_key = :RT
	elseif "rT" in names(param)
		rate = unique(param.rT)
		rate_key = :rT
	end
	for i=1:length(RT)
		param_f = filter([rate_key] => (x) -> x == rate[i], param)
		plot!(p_b, param_f.Tinit.+273.15, Measurements.value.(param_f.b), label=round(rate[i],sigdigits=3), m=:circle)

		plot!(p_T₀, param_f.Tinit, Measurements.value.(param_f.T₀), label=round(rate[i],sigdigits=3), m=:circle)
	
		plot!(p_m, param_f.Tinit.+273.15, Measurements.value.(param_f.m), label=round(rate[i],sigdigits=3), m=:circle)
	
		plot!(p_T₁, param_f.Tinit.+273.15, Measurements.value.(param_f.T₁), label=round(rate[i],sigdigits=3), m=:circle)
	end
	plot(p_T₀, p_T₁, p_b, p_m)
end

# ╔═╡ 4eac8fc5-3d30-44b6-828e-e40f17bd1d7e
md"""
## Show selected Fits
"""

# ╔═╡ 668bfd06-7308-4cc8-83e2-34d51f9730fd
md"""
### ``R_T`` and ``p_{in}``

``T_{init}``: $(@bind show_i_Tinit Slider(1:length(Tinit)))

``R_T``: $(@bind show_i_RT Slider(1:length(RT)))

``p_{in}``: $(@bind show_i_pin Slider(1:length(pin)))
"""

# ╔═╡ 789f8085-76c9-4fba-95d2-b6b6b603f868
function plot_fit(TeluTchar, param, fit, model, i_Tinit, i_RT, i_pin) 
	Tinit = unique(TeluTchar.Tinit)
	RT = unique(TeluTchar.RT)
	pin = unique(TeluTchar.pin)
	TeluTchar_f = filter([:Tinit, :RT, :pin] => (x,y,z) -> x == Tinit[i_Tinit] && y == RT[i_RT] && z == pin[i_pin], TeluTchar)
	param_f = filter([:Tinit, :RT, :pin] => (x,y,z) -> x == Tinit[i_Tinit] && y == RT[i_RT] && z == pin[i_pin], param)
	p_fit = scatter(TeluTchar_f.Tchar, TeluTchar_f.Telu, label="sim", title=string("fit Tinit=", Tinit[i_Tinit], "°C, RT=", RT[i_RT], "°C/min, pin=", pin[i_pin], "kPa(g)"))
	plot!(p_fit, TeluTchar_f.Tchar, model(TeluTchar_f.Tchar, fit[i_Tinit, i_RT, i_pin].param), label="R²=$(round(param_f.R²[1]; digits=4))")
	T₀ = param_f.T₀[1]
	T₁ = param_f.T₁[1]
	b = param_f.b[1]
	m = param_f.m[1]
	return p_fit, T₀, T₁, b, m  
end

# ╔═╡ 030af251-fe5b-46c8-b6d1-1b2efb7cd6c3
function plot_fit(TeluTchar, param, fit, model, i_Tinit, i_rT) 
	Tinit = unique(TeluTchar.Tinit)
	rT = unique(TeluTchar.rT)
	TeluTchar_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[i_Tinit] && y == rT[i_rT], TeluTchar)
	param_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[i_Tinit] && y == rT[i_rT], param)
	p_fit = scatter(TeluTchar_f.Tchar, TeluTchar_f.Telu, label="sim", title=string("fit Tinit=", Tinit[i_Tinit], "°C, rT=", round(rT[i_rT];sigdigits=3)))
	plot!(p_fit, TeluTchar_f.Tchar, model(TeluTchar_f.Tchar, fit[i_Tinit, i_rT].param), label="R²=$(round(param_f.R²[1]; digits=4))")
	T₀ = param_f.T₀[1]
	T₁ = param_f.T₁[1]
	b = param_f.b[1]
	m = param_f.m[1]
	return p_fit, T₀, T₁, b, m  
end

# ╔═╡ 65c91137-034b-43c5-a7d1-c25cc694b1d4
md"""
### ``r_T``

``T_{init}``: $(@bind show_i_Tinit_ Slider(1:length(Tinit)))

``r_T``: $(@bind show_i_rT_ Slider(1:length(rT)))
"""

# ╔═╡ c6143a16-c8d9-48ec-8588-7f430851576b
md"""
Data from Duong2022 (PAHs on Rxi17SilMS), Stultz2020 (Dioxins, Dibenzofurans on Rxi17SilMS) and some solutes from Brehmer2022 (PCBs on Rxi17SilMS) diverge from other data on the linear increasing side of ``T_{elu}(T_{char})`` with a smaller slope. Other data, e.g. Brehmer2022 (FAMES on Rxi17SilMS) and n-Alkanes on different stationary phases (Boswell2012 on DB5ms, Gaida2021 on Rxi5ms and Marquart on Rxi17SilMS) follows another slope of nearly ``m=1``.
"""

# ╔═╡ c3e32361-b417-4305-acf3-5862813d2d0c
md"""
The slope is depending on ``r_T`` and ``θ_{char}``.
"""

# ╔═╡ 19f116ba-bfd3-4579-9e45-b1afb61a9cbd
md"""
# MARK
"""

# ╔═╡ 5bc9a02e-62b1-48de-8287-33be309402fe
md"""
## Plot ``T_{elu}(T_{char}, θ_{char}, r_T, T_{init})``

rT: $(@bind select_i_rT__ Slider(1:length(rT)))

Tinit: $(@bind select_i_Tinit__ Slider(1:length(Tinit)))
"""

# ╔═╡ 6bef1367-bc12-40cd-8c6d-2d5934cdb9d6
md"""
## Fit ``r_T`` with filtered data
**For now: exclude Brehmer_PCB, Duong_FAME, Stultz_Dioxins from the data**
"""

# ╔═╡ 0e504ae7-4165-4261-a02d-85141e2f5b5f
md"""
``T_{init}``: $(@bind show_i_Tinit_1 Slider(1:length(Tinit)))

``r_T``: $(@bind show_i_rT_1 Slider(1:length(rT)))
"""

# ╔═╡ 16991c8c-b7c3-495b-93ce-70b949017b3c
md"""
### Fit ``m=1``
"""

# ╔═╡ 2e224d8b-6fb1-4b0b-aacb-28294f8e7b8d
md"""
# End
"""

# ╔═╡ 365ae9ff-9b94-457a-a687-d1c8581fc5ab
model(1.0,1.0)

# ╔═╡ 78e4ed32-9943-458a-aff8-3a9d0397e400
mm = model

# ╔═╡ 78941ab1-25bf-4b10-935b-4deff67fe5dc
mm==model_m1

# ╔═╡ 4fe853a2-100b-4c25-8559-1ea008c0cad2
function Rsquare(fit, y)
	sstot = sum((y.-mean(y)).^2)
	ssres = sum(fit.resid.^2)
	R2 = 1-ssres/sstot
	return R2
end

# ╔═╡ fd778a85-b7f5-4d87-b671-4ac3efee4b79
function fit_data(model, xdata, ydata, p0; lower=-Inf.*ones(length(p0)), upper=Inf.*ones(length(p0)))
	fit = curve_fit(model, xdata, ydata, p0, lower=lower, upper=upper)
	sigmas = NaN.*ones(length(p0))
	try
		sigmas = stderror(fit)
	catch
		sigmas = NaN.*ones(length(p0))
	end
	par = Array{Measurement{Float64}}(undef, length(p0))
	for i=1:length(p0)
		par[i] = fit.param[i] ± sigmas[i]
	end
	R2 = Rsquare(fit, ydata)
	return par, R2, fit
end

# ╔═╡ 7bc0108e-6842-4850-92a5-0cae5d7d4706
function fitting_RT_pin(TeluTchar, model)
	Tinit = sort!(unique(TeluTchar.Tinit))
	RT = sort!(unique(TeluTchar.RT))
	pin = sort!(unique(TeluTchar.pin))
	fits = Array{LsqFit.LsqFitResult}(undef, length(Tinit), length(RT), length(pin))
	Tinits = Float64[]
	RTs = Float64[]
	pins = Float64[]
	bs = Measurement{Float64}[]
	T₀s = Measurement{Float64}[]
	ms = Measurement{Float64}[]
	T₁s = Measurement{Float64}[]
	R²s = Float64[]
	converged = Bool[] 
	for i=1:length(Tinit)
		for j=1:length(RT)
			for k=1:length(pin)
				TeluTchar_f = filter([:Tinit, :RT, :pin] => (x,y,z) -> x == Tinit[i] && y == RT[j] && z == pin[k], TeluTchar)
				xdata = TeluTchar_f.Tchar
				ydata = TeluTchar_f.Telu
				if model == model4
					p0 = [Tinit[i], Tinit[i]+273.15, 0.05, 1.0]
					lb = [-Inf, -Inf, 0.0, 0.0]
					ub = [+Inf, +Inf, 0.5, 2.0]
					par, R2, fits[i,j,k] = fit_data(model, xdata, ydata, p0; lower=lb, upper=ub)
					push!(T₀s, par[1])
					push!(T₁s, par[2])
					push!(bs, par[3])
					push!(ms, par[4])
				elseif model == model_m1
					p0 = [Tinit[i], Tinit[i]+273.15, 0.05]
					lb = [-Inf, -Inf, 0.0]
					ub = [+Inf, +Inf, 0.5]
					par, R2, fits[i,j,k] = fit_data(model, xdata, ydata, p0; lower=lb, upper=ub)
					push!(T₀s, par[1])
					push!(T₁s, par[2])
					push!(bs, par[3])
					push!(ms, mconst)
				elseif model == model_b
					p0 = [Tinit[i], Tinit[i]+273.15, 1.0]
					lb = [-Inf, -Inf, 0.0]
					ub = [+Inf, +Inf, 2.0]
					par, R2, fits[i,j,k] = fit_data(model, xdata, ydata, p0; lower=lb, upper=ub)
					push!(T₀s, par[1])
					push!(T₁s, par[2])
					push!(bs, bconst)
					push!(ms, par[3])
				elseif model == model_b_m1
					p0 = [Tinit[i], Tinit[i]+273.15]
					lb = [-Inf, -Inf]
					ub = [+Inf, +Inf]
					par, R2, fits[i,j,k] = fit_data(model, xdata, ydata, p0; lower=lb, upper=ub)
					push!(T₀s, par[1])
					push!(T₁s, par[2])
					push!(bs, bconst)
					push!(ms, mconst)
				end
				push!(Tinits, Tinit[i])
				push!(RTs, RT[j])
				push!(pins, pin[k])
				push!(R²s, R2)
				push!(converged, fits[i,j,k].converged)
			end
		end
	end
	params = DataFrame(Tinit=Tinits, RT=RTs, pin=pins, T₀=T₀s, T₁=T₁s, b=bs, m=ms, R²=R²s, converged=converged)
	return params, fits
end

# ╔═╡ 2f750d6e-461f-4e6a-bc92-2939dbb55421
function fitting_rT(TeluTchar, model)
	Tinit = sort!(unique(TeluTchar.Tinit))
	rT = sort!(unique(TeluTchar.rT))
	fits = Array{LsqFit.LsqFitResult}(undef, length(Tinit), length(rT))
	Tinits = Float64[]
	rTs = Float64[]
	pins = Float64[]
	bs = Measurement{Float64}[]
	T₀s = Measurement{Float64}[]
	ms = Measurement{Float64}[]
	T₁s = Measurement{Float64}[]
	R²s = Float64[]
	converged = Bool[] 
	for i=1:length(Tinit)
		for j=1:length(RT)
			TeluTchar_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[i] && y == rT[j], TeluTchar)
			xdata = TeluTchar_f.Tchar
			ydata = TeluTchar_f.Telu
			if model == model4
				p0 = [Tinit[i], Tinit[i]+273.15, 0.05, 1.0]
				lb = [-Inf, -Inf, 0.0, 0.0]
				ub = [+Inf, +Inf, 0.5, 2.0]
				par, R2, fits[i,j] = fit_data(model, xdata, ydata, p0; lower=lb, upper=ub)
				push!(T₀s, par[1])
				push!(T₁s, par[2])
				push!(bs, par[3])
				push!(ms, par[4])
			elseif model == model_m1
				p0 = [Tinit[i], Tinit[i]+273.15, 0.05]
				lb = [-Inf, -Inf, 0.0]
				ub = [+Inf, +Inf, 0.5]
				par, R2, fits[i,j] = fit_data(model, xdata, ydata, p0; lower=lb, upper=ub)
				push!(T₀s, par[1])
				push!(T₁s, par[2])
				push!(bs, par[3])
				push!(ms, mconst)
			elseif model == model_b
				p0 = [Tinit[i], Tinit[i]+273.15, 1.0]
				lb = [-Inf, -Inf, 0.0]
				ub = [+Inf, +Inf, 2.0]
				par, R2, fits[i,j] = fit_data(model, xdata, ydata, p0; lower=lb, upper=ub)
				push!(T₀s, par[1])
				push!(T₁s, par[2])
				push!(bs, bconst)
				push!(ms, par[3])
			elseif model == model_b_m1
				p0 = [Tinit[i], Tinit[i]+273.15]
				lb = [-Inf, -Inf]
				ub = [+Inf, +Inf]
				par, R2, fits[i,j] = fit_data(model, xdata, ydata, p0; lower=lb, upper=ub)
				push!(T₀s, par[1])
				push!(T₁s, par[2])
				push!(bs, bconst)
				push!(ms, mconst)
			end
			push!(Tinits, Tinit[i])
			push!(rTs, rT[j])
			push!(R²s, R2)
			push!(converged, fits[i,j].converged)
		end
	end
	params = DataFrame(Tinit=Tinits, rT=rTs, T₀=T₀s, T₁=T₁s, b=bs, m=ms, R²=R²s, converged=converged)
	return params, fits
end

# ╔═╡ 5e26e413-612c-4728-9c08-15f45aa10e1f
function fitting_rT(TeluTchar)
	Tinit = sort!(unique(TeluTchar.Tinit))
	rT = sort!(unique(TeluTchar.rT))
	fits = Array{LsqFit.LsqFitResult}(undef, length(Tinit), length(rT))
	Tinits = Float64[]
	rTs = Float64[]
	bs = Measurement{Float64}[]
	T₀s = Measurement{Float64}[]
	ms = Measurement{Float64}[]
	T₁s = Measurement{Float64}[]
	R²s = Float64[]
	for i=1:length(Tinit)
		for j=1:length(rT)
			TeluTchar_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[i] && y == rT[j], TeluTchar)
			xdata = TeluTchar_f.Tchar
			ydata = TeluTchar_f.Telu
			p0 = [0.05, Tinit[i], 1.0, Tinit[i]+273.15]
			lb = [0.0, -Inf, 0.0, -Inf]
				ub = [0.5, +Inf, 2.0, +Inf]
			fits[i,j] = curve_fit(model4, xdata, ydata, p0, lower=lb, upper=ub)
			if fits[i,j].converged == true
				sigmas = NaN.*ones(4)
				params = NaN.*ones(4)
				try
					sigmas = stderror(fits[i,j])
					params = fits[i,j].param
				catch
					sigmas = NaN.*ones(4)
					params = NaN.*ones(4)
				end
				push!(Tinits, Tinit[i])
				push!(rTs, rT[j])
				push!(bs, params[1] ± sigmas[1])
				push!(T₀s, params[2] ± sigmas[2])
				push!(ms, params[3] ± sigmas[3])
				push!(T₁s, params[4] ± sigmas[4])
				push!(R²s, Rsquare(fits[i,j], ydata))
			end
		end
	end
	params = DataFrame(Tinit=Tinits, rT=rTs, b=bs, T₀=T₀s, m=ms, T₁=T₁s, R²=R²s)
	return params, fits
end

# ╔═╡ 1739ae21-a60d-413e-832d-cb12fad0ca34
function fitting(TeluTchar, model)
	if "RT" in names(TeluTchar)
		params, fits = fitting_RT_pin(TeluTchar, model)
	elseif "rT" in names(TeluTchar)
		params, fits = fitting_rT(TeluTchar, model)
	end
	return params, fits
end	

# ╔═╡ 804fd07b-b4f0-4119-a708-60d4973e84eb
function fitting_rT_m1(TeluTchar)
	Tinit = sort!(unique(TeluTchar.Tinit))
	rT = sort!(unique(TeluTchar.rT))
	fits = Array{LsqFit.LsqFitResult}(undef, length(Tinit), length(rT))
	Tinits = Float64[]
	rTs = Float64[]
	bs = Measurement{Float64}[]
	T₀s = Measurement{Float64}[]
	T₁s = Measurement{Float64}[]
	R²s = Float64[]
	for i=1:length(Tinit)
		for j=1:length(rT)
			TeluTchar_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[i] && y == rT[j], TeluTchar)
			xdata = TeluTchar_f.Tchar
			ydata = TeluTchar_f.Telu
			p0 = [0.05, Tinit[i], Tinit[i]+273.15]
			lb = [0.0, -Inf, -Inf]
			ub = [0.5, +Inf, +Inf]
			fits[i,j] = curve_fit(model_m1, xdata, ydata, p0, lower=lb, upper=ub)
			if fits[i,j].converged == true
				sigmas = NaN.*ones(3)
				params = NaN.*ones(3)
				try
					sigmas = stderror(fits[i,j])
					params = fits[i,j].param
				catch
					sigmas = NaN.*ones(3)
					params = NaN.*ones(3)
				end
				push!(Tinits, Tinit[i])
				push!(rTs, rT[j])
				push!(bs, params[1] ± sigmas[1])
				push!(T₀s, params[2] ± sigmas[2])
				push!(T₁s, params[3] ± sigmas[3])
				push!(R²s, Rsquare(fits[i,j], ydata))
			end
		end
	end
	params = DataFrame(Tinit=Tinits, rT=rTs, b=bs, T₀=T₀s, T₁=T₁s, R²=R²s)
	return params, fits
end

# ╔═╡ 59e2f85f-acdd-4ca5-a91c-78b9bb42364e
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

# ╔═╡ 3f7dad44-11bc-4535-b7cb-db5938f5d2d4
df.par[5].sub[1]

# ╔═╡ 98df485c-de5b-4ea9-ba12-7d0677a4e7dc
function TeluTchar(data)
	sp = sort!(unique(data.sp))
	Tinit = sort!(unique(data.Tinit))
	RT = sort!(unique(data.rate))
	pin = sort!(unique(data.Fpin))
	#TeluTchar_ = Array{DataFrame}(undef, length(Tinit), length(RT), length(pin), length(sp))
	Tchars = Float64[]
	Telus = Float64[]
	Names = String[]
	kelus = Float64[]
	sps = String[]
	Tinits = Float64[]
	RTs = Float64[]
	pins = Float64[]
	Anns = String[]
	θchars = Float64[]
	ΔCps = Float64[]
	for i1=1:length(Tinit)
		for i2=1:length(RT)
			for i3=1:length(pin)
				for i4=1:length(sp)
					data_f = filter([:Tinit, :rate, :Fpin, :sp] => (x1, x2, x3, x4) -> x1 == Tinit[i1] && x2 == RT[i2] && x3 == pin[i3] && x4 == sp[i4], data)
					for j=1:length(data_f.par)
						for k=1:length(data_f.par[j].sub)
							kk = findfirst(data_f.par[j].sub[k].name.==data_f.peaklist[j].Name)
							push!(Tinits, Tinit[i1])
							push!(RTs, RT[i2])
							push!(pins, pin[i3])
							push!(sps, sp[i4])
							push!(Tchars, data_f.par[j].sub[k].Tchar)
							push!(Telus, data_f.peaklist[j].TR[kk] + 273.15)
							push!(Names, data_f.par[j].sub[k].name)
							push!(kelus, data_f.peaklist[j].kR[kk])
							push!(Anns, split(data_f.par[j].sub[k].ann, ',')[1])
							push!(θchars, data_f.par[j].sub[k].θchar)
							push!(ΔCps, data_f.par[j].sub[k].ΔCp)
						end
					end
				end
			end
		end
	end	
	TeluTchar_ = filter!([:Telu, :Tchar] => (x,y) -> !isnan(x) && y < 700.0, sort!(DataFrame(Name=Names, Tinit=Tinits, RT=RTs, pin=pins, sp=sps, Tchar=Tchars, Telu=Telus, kelu=kelus, ann=Anns, θchar=θchars, ΔCp=ΔCps), [:Tchar])) # sort Tchar and filter out NaN values
	return TeluTchar_
end

# ╔═╡ fb743e94-7f59-45cf-b329-b980d3403af2
Telu_Tchar = TeluTchar(df)

# ╔═╡ dac5f3f5-1327-4b65-be25-46764829e5b4
# plot Telu(Tchar, Tinit) for fixed rT (resp. fixed RT and pin) and variation of sp
begin
	plotly()
	p___ = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu", title=string("RT = ", RT[select_i_RT], "°C/min, pin = ", pin[select_i_pin], "kPa(g)"))
	for i=1:length(sp)
		Telu_Tchar_f = filter([:sp, :RT, :pin] => (x,y,z) -> x == sp[i] && y == RT[select_i_RT] && z == pin[select_i_pin], Telu_Tchar)
		scatter!(p___, Telu_Tchar_f.Tchar, Telu_Tchar_f.Tinit, Telu_Tchar_f.Telu, label=sp[i])
	end
	p___
	md"""

	$(embed_display(p___))

	For high heating rates, ``R_T \geq 50`` and all inlet pressures (but extremely for low values), certain solutes on certain stationary phases (Rxi17SilMS, Rxi5ms, SPB50) show a to high ``T_{elu}`` value compared to other solutes, especially in the region, where the solutes should elute with nearly ``T_{init}``.

	-> exclude the simulations with ``R_T \geq 50``.
	"""
end

# ╔═╡ 7d0b670d-03de-40e3-9eef-6f38289d3230
begin
	f = filter([:RT, :pin, :Tinit, :Telu, :Tchar] => (x1,x2,x3,x5,x6) -> x1 == 100.0 && x2 == 50.0 && x3 == 320.0 && x5 >= 890.0 && x6 <= 400.0, Telu_Tchar)
	md"""

	Identification of the biggest outliers:

	``R_T = 100``°C/min, ``p_{in}=50``kPa(g), ``T_{init}=320``°C, ``T_{elu}>890``K, ``T_{char}<400``K

	$(f)
	"""
end

# ╔═╡ 3b21afa0-c6b2-4531-974b-f4f82d8c7b40
param_RT_pin, fit_RT_pin = fitting(Telu_Tchar, select_model);

# ╔═╡ 4dd137ec-45d7-4be3-934f-454adefb8203
plot_parameters_i(filter([:converged] => x -> x == true, param_RT_pin))

# ╔═╡ adea9d70-dcda-4792-b365-d15274855fd8
plot_parameters_Tinit(filter([:pin] => x -> x == select_pin,param_RT_pin))

# ╔═╡ f54b0681-3a3e-45cf-b274-4962989821dc
begin
	p_fit_RT_pin = plot_fit(Telu_Tchar, param_RT_pin, fit_RT_pin, model4, show_i_Tinit, show_i_RT, show_i_pin) 
	md"""
	#### Fit result
	
	``T_0 =`` $(p_fit_RT_pin[2]) K

	``T_1 =`` $(p_fit_RT_pin[3]) K

	``b =`` $(p_fit_RT_pin[4])
	
	``m =`` $(p_fit_RT_pin[5]) K⁻¹

	$(embed_display(p_fit_RT_pin[1]))
	"""
end

# ╔═╡ 7290ff57-ead8-4518-8506-0f845cf14d0f
"RT" in names(Telu_Tchar)

# ╔═╡ 2a0e8394-076a-4d0f-8458-b2f838c9bd77
function TeluTchar_dimless_rate(data)
	sp = sort!(unique(data.sp))
	Tinit = sort!(unique(data.Tinit))
	rT = sort!(unique(data.dimless_rate))
	#TeluTchar_ = Array{DataFrame}(undef, length(Tinit), length(RT), length(pin), length(sp))
	Tchars = Float64[]
	Telus = Float64[]
	Names = String[]
	kelus = Float64[]
	sps = String[]
	Tinits = Float64[]
	rTs = Float64[]
	Anns = String[]
	θchars = Float64[]
	ΔCps = Float64[]
	for i1=1:length(Tinit)
		for i2=1:length(rT)
			for i3=1:length(sp)
					data_f = filter([:Tinit, :dimless_rate, :sp] => (x1, x2, x3) -> x1 == Tinit[i1] && x2 == rT[i2] && x3 == sp[i3], data)
					for j=1:length(data_f.par)
						for k=1:length(data_f.par[j].sub)
							kk = findfirst(data_f.par[j].sub[k].name.==data_f.peaklist[j].Name)
							push!(Tinits, Tinit[i1])
							push!(rTs, rT[i2])
							push!(sps, sp[i3])
							push!(Tchars, data_f.par[j].sub[k].Tchar)
							push!(Telus, data_f.peaklist[j].TR[kk] + 273.15)
							push!(Names, data_f.par[j].sub[k].name)
							push!(kelus, data_f.peaklist[j].kR[kk])
							push!(Anns, split(data_f.par[j].sub[k].ann, ',')[1])
							push!(θchars, data_f.par[j].sub[k].θchar)
							push!(ΔCps, data_f.par[j].sub[k].ΔCp)
						end
					end
			end
		end
	end	
	TeluTchar_ = filter!([:Telu, :Tchar] => (x,y) -> !isnan(x) && y < 700.0, sort!(DataFrame(Name=Names, Tinit=Tinits, rT=rTs, sp=sps, Tchar=Tchars, Telu=Telus, kelu=kelus, ann=Anns, θchar=θchars, ΔCp=ΔCps), [:Tchar])) # sort Tchar and filter out NaN values
	return TeluTchar_
end

# ╔═╡ f9e7b6a4-3e68-4adc-abeb-b9a6abbaf068
Telu_Tchar_rT = TeluTchar_dimless_rate(df)

# ╔═╡ 922865aa-9d53-4231-b581-7a160a395420
begin
	plotly()
	p = plot(xlabel="Tchar", ylabel="rT", zlabel="Telu")
	for i=1:length(Tinit)
		Telu_Tchar_rT_f = filter([:Tinit] => x -> x == Tinit[i], Telu_Tchar_rT)
		scatter!(p, Telu_Tchar_rT_f.Tchar, Telu_Tchar_rT_f.rT, Telu_Tchar_rT_f.Telu, label=Tinit[i])
	end
	p
end

# ╔═╡ 4d669db9-1661-4aa6-87f4-2977ae11e989
begin
	plotly()
	p_ = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	for i=1:length(rT)
		Telu_Tchar_rT_f = filter([:rT] => x -> x == rT[i], Telu_Tchar_rT)
		scatter!(p_, Telu_Tchar_rT_f.Tchar, Telu_Tchar_rT_f.Tinit, Telu_Tchar_rT_f.Telu, label=round(rT[i];sigdigits=3))
	end
	p_
end

# ╔═╡ 4e46e5eb-74ca-400b-bca2-7329af1879ac
# plot Telu(Tchar, Tinit) for fixed rT (resp. fixed RT and pin) and variation of sp
begin
	#select_i_rT = findfirst(select_rT.==rT)
	plotly()
	p__ = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu", title=string("rT = ", round(rT[select_i_rT]; sigdigits=3)))
	for i=1:length(sp)
		Telu_Tchar_rT_f = filter([:sp, :rT] => (x,y) -> x == sp[i] && y == rT[select_i_rT], Telu_Tchar_rT)
		scatter!(p__, Telu_Tchar_rT_f.Tchar, Telu_Tchar_rT_f.Tinit, Telu_Tchar_rT_f.Telu, label=sp[i])
	end
	md"""

	$(embed_display(p__))

	For extreme high dimensionless heating rates, ``r_T>100``, certain solutes on certain stationary phases (Rxi17SilMS, Rxi5ms, SPB50) show a to high ``T_{elu}`` value compared to other solutes, especially in the region, where the solutes should elute with nearly ``T_{init}``.

	-> exclude the simulations with ``r_T>100``.
	"""
end

# ╔═╡ 2f374d9f-54a4-48ce-8aa0-22ae407fd5b3
param_rT, fit_rT = fitting(Telu_Tchar_rT, select_model_);

# ╔═╡ 1c37e89a-6a90-4926-bb37-c4f7f445c9f8
plot_parameters_i(filter([:converged] => x -> x == true, param_rT))

# ╔═╡ b373c92c-d465-4a0f-a925-ae50360c277e
plot_parameters_Tinit(filter([:converged] => x -> x == true, param_rT))

# ╔═╡ 7c07e06f-5c31-4453-8d0a-d89842ebe7d2
begin
	p_fit_rT = plot_fit(Telu_Tchar_rT, param_rT, fit_rT, model4, show_i_Tinit_, show_i_rT_) 
	md"""
	#### Fit result
	
	``T_0 =`` $(p_fit_rT[2]) K

	``T_1 =`` $(p_fit_rT[3]) K

	``b =`` $(p_fit_rT[4])
	
	``m =`` $(p_fit_rT[5]) K⁻¹

	$(embed_display(p_fit_rT[1]))
	"""
end

# ╔═╡ e97005f1-f165-4dfb-8be0-4d078eb27930
begin
	TeluTchar_rT_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[show_i_Tinit_] && y == rT[show_i_rT_], Telu_Tchar_rT)
	param_rT_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[show_i_Tinit_] && y == rT[show_i_rT_], param_rT)
	p_fit_sp = scatter(title=string("fit Tinit=", Tinit[show_i_Tinit_], "°C, rT=", round(rT[show_i_rT_];sigdigits=3)))
	for i=1:length(sp)
		TeluTchar_rT_f_sp = filter([:sp] => (x) -> x == sp[i], TeluTchar_rT_f)
		scatter!(p_fit_sp, TeluTchar_rT_f_sp.Tchar, TeluTchar_rT_f_sp.Telu, label=sp[i])
	end
	plot!(p_fit_sp, TeluTchar_rT_f.Tchar, select_model(TeluTchar_rT_f.Tchar, Measurements.value.([param_rT_f.T₀[1], param_rT_f.T₁[1], param_rT_f.b[1], param_rT_f.m[1]])), label="R²=$(round(param_rT_f.R²[1]; digits=4))")	
	md"""
	#### Plot Fit result for different stationary phases

	$(embed_display(p_fit_sp))
	"""
end

# ╔═╡ bd87d080-07fc-4c99-9e6c-1df8229cb46c
begin
	brehmer_PCB = filter([:ann, :Tchar, :Name] => (x,y,z) -> x == "Brehmer2022" && y>=530.0 && occursin("PCB", z), TeluTchar_rT_f)
	brehmer_FAME = filter([:ann, :Tchar, :Name] => (x,y,z) -> x == "Brehmer2022" && y>=530.0 && !occursin("PCB", z), TeluTchar_rT_f)
	duong_PAH = filter([:ann, :Tchar] => (x,y) -> x == "Duong2022"  && y>=530.0, TeluTchar_rT_f)
	stultz_Dioxins = filter([:ann, :Tchar] => (x,y) -> x == "Stultz2020"  && y>=530.0, TeluTchar_rT_f)
	boswell_alkanes = filter([:ann, :Tchar, :Name] => (x,y,z) -> x == "Boswell2012"  && y>=530.0 && occursin("C",z), TeluTchar_rT_f)
	boswell_nonalkanes = filter([:ann, :Tchar, :Name] => (x,y,z) -> x == "Boswell2012"  && y>=530.0 && !occursin("C",z), TeluTchar_rT_f)
	marquart_alkanes = filter([:ann, :Tchar,:Name] => (x,y,z) -> x == "Marquart2020"  && y>=530.0 && occursin("C",z), TeluTchar_rT_f)
	marquart_nonalkanes = filter([:ann, :Tchar,:Name] => (x,y,z) -> x == "Marquart2020"  && y>=530.0 && !occursin("C",z), TeluTchar_rT_f)
	gaida_alkanes = filter([:ann, :Tchar] => (x,y) -> x == "Gaida2021"  && y>=530.0, TeluTchar_rT_f)
	p1 = scatter(brehmer_PCB.Tchar, brehmer_PCB.θchar, xlabel="Tchar in K", ylabel="θchar in °C", label="Brehmer_PCB", palette=:tab10)
	scatter!(p1, brehmer_FAME.Tchar, brehmer_FAME.θchar, label="Brehmer_FAME")
	scatter!(p1, duong_PAH.Tchar, duong_PAH.θchar, label="Duong_PAH")
	scatter!(p1, stultz_Dioxins.Tchar, stultz_Dioxins.θchar, label="stultz_Dioxins")
	scatter!(p1, boswell_alkanes.Tchar, boswell_alkanes.θchar, label="boswell_alkanes")
	scatter!(p1, boswell_nonalkanes.Tchar, boswell_nonalkanes.θchar, label="boswell_nonalkanes")
	scatter!(p1, marquart_alkanes.Tchar, marquart_alkanes.θchar, label="marquart_alkanes")
	scatter!(p1, marquart_nonalkanes.Tchar, marquart_nonalkanes.θchar, label="marquart_nonalkanes")
	scatter!(p1, gaida_alkanes.Tchar, gaida_alkanes.θchar, label="gaida_alkanes")
	p2 = scatter(brehmer_PCB.Tchar, brehmer_PCB.ΔCp, xlabel="Tchar in K", ylabel="ΔCp in J mol⁻¹ K⁻¹", label="Brehmer_PCB", palette=:tab10)
	scatter!(p2, brehmer_FAME.Tchar, brehmer_FAME.ΔCp, label="Brehmer_FAME")
	scatter!(p2, duong_PAH.Tchar, duong_PAH.ΔCp, label="Duong_PAH")
	scatter!(p2, stultz_Dioxins.Tchar, stultz_Dioxins.ΔCp, label="stultz_Dioxins")
	scatter!(p2, boswell_alkanes.Tchar, boswell_alkanes.ΔCp, label="boswell_alkanes")
	scatter!(p2, boswell_nonalkanes.Tchar, boswell_nonalkanes.ΔCp, label="boswell_nonalkanes")
	scatter!(p2, marquart_alkanes.Tchar, marquart_alkanes.ΔCp, label="marquart_alkanes")
	scatter!(p2, marquart_nonalkanes.Tchar, marquart_nonalkanes.ΔCp, label="marquart_nonalkanes")
	scatter!(p2, gaida_alkanes.Tchar, gaida_alkanes.ΔCp, label="gaida_alkanes")

	p3 = scatter(brehmer_PCB.Tchar, brehmer_PCB.Telu, xlabel="Tchar in K", ylabel="Telu K", label="Brehmer_PCB", palette=:tab10)
	scatter!(p3, brehmer_FAME.Tchar, brehmer_FAME.Telu, label="Brehmer_FAME")
	scatter!(p3, duong_PAH.Tchar, duong_PAH.Telu, label="Duong_PAH")
	scatter!(p3, stultz_Dioxins.Tchar, stultz_Dioxins.Telu, label="stultz_Dioxins")
	scatter!(p3, boswell_alkanes.Tchar, boswell_alkanes.Telu, label="boswell_alkanes")
	scatter!(p3, boswell_nonalkanes.Tchar, boswell_nonalkanes.Telu, label="boswell_nonalkanes")
	scatter!(p3, marquart_alkanes.Tchar, marquart_alkanes.Telu, label="marquart_alkanes")
	scatter!(p3, marquart_nonalkanes.Tchar, marquart_nonalkanes.Telu, label="marquart_nonalkanes")
	scatter!(p3, gaida_alkanes.Tchar, gaida_alkanes.Telu, label="gaida_alkanes")
	plot!(p3, 530.0:1.0:680.0, select_model_(530.0:1.0:680.0, fit_rT[show_i_Tinit_, show_i_rT_].param), label="R²=$(round(param_rT_f.R²[1]; digits=4))", c=:black, s=:dash, w=2)
	
	plot(p1, p2, p3)
end

# ╔═╡ e950d80d-8b1a-427c-a4f0-361d62f9a450
begin
	plotly()
	p___a = plot(xlabel="Tchar", ylabel="θchar", zlabel="Telu", title=string("rT = ", rT[select_i_rT__], ", Tinit = ", Tinit[select_i_Tinit__], "°C"))
	for i=1:length(sp)
		Telu_Tchar_f = filter([:sp, :rT, :Tinit] => (x,y,z) -> x == sp[i] && y == rT[select_i_rT__] && z == Tinit[select_i_Tinit__], Telu_Tchar_rT)
		scatter!(p___a, Telu_Tchar_f.Tchar, Telu_Tchar_f.θchar, Telu_Tchar_f.Telu, label=sp[i])
	end
	p___a
end

# ╔═╡ 600773b7-8b3b-4e30-91ca-3ce0407dae68
begin
	plotly()
	p___b = plot(ylabel="Tchar", xlabel="θchar", zlabel="Telu", title=string("rT = ", rT[1], ", Tinit = ", Tinit[4], "°C"))
	for i=1:length(sp)
		Telu_Tchar_f = filter([:sp, :rT, :Tinit] => (x,y,z) -> x == sp[i] && y == rT[1] && z == Tinit[4], Telu_Tchar_rT)
		scatter!(p___b, Telu_Tchar_f.θchar, Telu_Tchar_f.Tchar, Telu_Tchar_f.Telu, label=sp[i])
	end
	p___c = plot(ylabel="Tchar", xlabel="θchar", zlabel="Telu", title=string("rT = ", rT[end], ", Tinit = ", Tinit[4], "°C"))
	for i=1:length(sp)
		Telu_Tchar_f = filter([:sp, :rT, :Tinit] => (x,y,z) -> x == sp[i] && y == rT[end] && z == Tinit[4], Telu_Tchar_rT)
		scatter!(p___c, Telu_Tchar_f.θchar, Telu_Tchar_f.Tchar, Telu_Tchar_f.Telu, label=sp[i])
	end

	
	md"""
	### Note
	
	For low dimensionless heating rates ``r_T`` the slope of ``T_{elu}(T_{char})`` depends strongly on ``θ_{char}``. With increasing ``θ_{char}`` the slope decreases.

	$(embed_display(p___b))
	
	In the 3D-plot of ``T_{elu}(T_{char}, θ_{char})`` for a low ``r_T`` value, the plane, in which the ``T_{elu}`` data is located is twisted and it is not parallel to the ``T_{char}`` axis. 

	$(embed_display(p___c))

	For increasing ``r_T`` the plane rotates to a parallel position to the ``T_{char}`` axis.
	"""
end

# ╔═╡ aa34de2c-dbf1-4551-aac7-94f361fa9b6a
begin
	Telu_Tchar_rT_data_excluded = filter([:ann, :Name] => (x,y) -> (x == "Brehmer2022" &&  !occursin("PCB", y)) && x != "Duong2022" && x != "Stultz2020", Telu_Tchar_rT)
end

# ╔═╡ 43d45ad1-0289-4992-a089-22c6b299931d
param_rT_data_excluded, fit_rT_data_excluded = fitting_rT(Telu_Tchar_rT_data_excluded);

# ╔═╡ f6677948-c480-4557-ac40-e38a76470daf
param_rT_data_excluded

# ╔═╡ fc86ff7a-bbab-4278-a1dd-5dfe3da87adf
fit_rT_data_excluded[1,6]

# ╔═╡ 7d5fdf32-1619-4835-976c-8d0513efb24c
begin
	p_b_rT_data_excluded = plot(xlabel="Tinit in K", ylabel="b")
	for i=1:length(rT)
		param_f = filter([:rT] => x -> x == rT[i], param_rT_data_excluded)
		plot!(p_b_rT_data_excluded, param_f.Tinit.+273.15, Measurements.value.(param_f.b), label=round(rT[i];sigdigits=3), m=:circle)
	end
	p_T₀_rT_data_excluded = plot(xlabel="Tinit in °C", ylabel="T₀")
	for i=1:length(rT)
		param_f = filter([:rT] => (x) -> x == rT[i], param_rT_data_excluded)
		plot!(p_T₀_rT_data_excluded, param_f.Tinit, Measurements.value.(param_f.T₀), label=round(rT[i];sigdigits=3), m=:circle)
	end
	p_m_rT_data_excluded = plot(xlabel="Tinit in K", ylabel="m")
	for i=1:length(rT)
		param_f = filter([:rT] => (x) -> x == rT[i], param_rT_data_excluded)
		plot!(p_m_rT_data_excluded, param_f.Tinit.+273.15, Measurements.value.(param_f.m), label=round(rT[i];sigdigits=3), m=:circle)
	end
	p_T₁_rT_data_excluded = plot(xlabel="Tinit in K", ylabel="T₁")
	for i=1:length(rT)
		param_f = filter([:rT] => (x) -> x == rT[i], param_rT_data_excluded)
		plot!(p_T₁_rT_data_excluded, param_f.Tinit.+273.15, Measurements.value.(param_f.T₁), label=round(rT[i];sigdigits=3), m=:circle)
	end
	plot(p_b_rT_data_excluded, p_T₀_rT_data_excluded, p_m_rT_data_excluded, p_T₁_rT_data_excluded)
end

# ╔═╡ 10266ff0-7e64-41c9-b420-d6c15257ceef
begin
	TeluTchar_rT_f1 = filter([:Tinit, :rT] => (x,y) -> x == Tinit[show_i_Tinit_1] && y == rT[show_i_rT_1], Telu_Tchar_rT_data_excluded)
	param_rT_f1 = filter([:Tinit, :rT] => (x,y) -> x == Tinit[show_i_Tinit_1] && y == rT[show_i_rT_1], param_rT_data_excluded)
	p_fit_rT_1 = scatter(title=string("fit Tinit=", Tinit[show_i_Tinit_1], "°C, rT=", round(rT[show_i_rT_1];sigdigits=3)), xlabel="Tchar in K", ylabel="Telu in K")
	#Ann = unique(Telu_Tchar.ann)
	for i=1:length(Ann)
		TeluTchar_rT_ff = filter([:ann] => (z) -> z == Ann[i], TeluTchar_rT_f1)
		scatter!(p_fit_rT_1, TeluTchar_rT_ff.Tchar, TeluTchar_rT_ff.Telu, label=Ann[i])
	end
	plot!(p_fit_rT_1, TeluTchar_rT_f1.Tchar, model(TeluTchar_rT_f1.Tchar, fit_rT_data_excluded[show_i_Tinit_1, show_i_rT_1].param), label="R²=$(round(param_rT_f1.R²[1]; digits=4))")
	md"""
	#### Fit result
	``b =`` $(param_rT_f1.b[1])
	
	``T_0 =`` $(param_rT_f1.T₀[1]) K
	
	``m =`` $(param_rT_f1.m[1]) K⁻¹
	
	``T_1 =`` $(param_rT_f1.T₁[1]) K

	``R^2 =`` $(round(param_rT_f1.R²[1];sigdigits=4)) 

	$(embed_display(p_fit_rT_1))
	"""
end

# ╔═╡ f917c1c9-8bb8-4e35-bb92-24df9e2945a1
param_rT_data_excluded_m1, fit_rT_data_excluded_m1 = fitting_rT_m1(Telu_Tchar_rT_data_excluded);

# ╔═╡ 9392b89e-6a64-4a7c-80dd-c6124820bdb3
begin
	p_b_rT_data_excluded_m1 = plot(xlabel="Tinit in K", ylabel="b")
	for i=1:length(rT)
		param_f = filter([:rT] => x -> x == rT[i], param_rT_data_excluded_m1)
		plot!(p_b_rT_data_excluded_m1, param_f.Tinit.+273.15, Measurements.value.(param_f.b), label=round(rT[i];sigdigits=3), m=:circle)
	end
	p_T₀_rT_data_excluded_m1 = plot(xlabel="Tinit in °C", ylabel="T₀")
	for i=1:length(rT)
		param_f = filter([:rT] => (x) -> x == rT[i], param_rT_data_excluded_m1)
		plot!(p_T₀_rT_data_excluded_m1, param_f.Tinit, Measurements.value.(param_f.T₀), label=round(rT[i];sigdigits=3), m=:circle)
	end
	p_m_rT_data_excluded_m1 = plot(xlabel="Tinit in K", ylabel="m")
	for i=1:length(rT)
		param_f = filter([:rT] => (x) -> x == rT[i], param_rT_data_excluded_m1)
		plot!(p_m_rT_data_excluded_m1, param_f.Tinit.+273.15, ones(length(param_f.Tinit)), label=round(rT[i];sigdigits=3), m=:circle)
	end
	p_T₁_rT_data_excluded_m1 = plot(xlabel="Tinit in K", ylabel="T₁")
	for i=1:length(rT)
		param_f = filter([:rT] => (x) -> x == rT[i], param_rT_data_excluded_m1)
		plot!(p_T₁_rT_data_excluded_m1, param_f.Tinit.+273.15, Measurements.value.(param_f.T₁), label=round(rT[i];sigdigits=3), m=:circle)
	end
	plot(p_b_rT_data_excluded_m1, p_T₀_rT_data_excluded_m1, p_m_rT_data_excluded_m1, p_T₁_rT_data_excluded_m1)
end

# ╔═╡ f518d526-34ca-4c15-8255-4ee79c8581dd
begin
	xx = param_rT_data_excluded_m1.Tinit[isnan.(param_rT_data_excluded_m1.b).==false]
	yy = param_rT_data_excluded_m1.rT[isnan.(param_rT_data_excluded_m1.b).==false]
	bb = param_rT_data_excluded_m1.b[isnan.(param_rT_data_excluded_m1.b).==false]
	p3D_b = scatter(xx, yy, Measurements.value.(bb), xlabel="Tinit in °C", ylabel="rT", zlabel="b")
	T₀T₀ = param_rT_data_excluded_m1.T₀[isnan.(param_rT_data_excluded_m1.T₀).==false]
	p3D_T₀ = scatter(xx, yy, Measurements.value.(T₀T₀), xlabel="Tinit in °C", ylabel="rT", zlabel="T₀")
	T₁T₁ = param_rT_data_excluded_m1.T₁[isnan.(param_rT_data_excluded_m1.T₁).==false]
	p3D_T₁ = scatter(xx, yy, Measurements.value.(T₁T₁), xlabel="Tinit in °C", ylabel="rT", zlabel="T₁")
	plot(p3D_b, p3D_T₀, p3D_T₁)
end

# ╔═╡ 2869dce2-0151-4ead-9ba4-39ac12876c95
isnan.(param_rT_data_excluded_m1.b).==false

# ╔═╡ 512a39c1-1fd0-4d6b-9d45-e45457dfabfb
fit_rT_data_excluded_m1[5,6]

# ╔═╡ Cell order:
# ╠═f8ad60a7-d51c-40b0-8f82-1fc78dcd4b54
# ╠═562d78f0-6e8d-4b44-9cb3-70acd226f439
# ╠═2a8dd6b6-dc9c-453b-9652-2226fc1f39ff
# ╠═e5f30f41-3205-4746-969f-3d7c106e11bf
# ╠═650a0f9b-77c8-4918-adf5-28fae3c7dd35
# ╠═e069b363-fc28-4665-9d27-f95371433dd2
# ╠═3a85a6d5-d20c-4245-b8f7-5779d57e6f04
# ╠═56822d93-0c5f-44a1-a9fe-f034b9561d28
# ╠═17b52208-736e-4b21-a584-c49a9849e1ec
# ╠═b7d6a89b-6ce2-426e-90d2-4ea289a407b8
# ╠═825ca4e1-02b4-4415-9a92-44dd92e94487
# ╠═fb743e94-7f59-45cf-b329-b980d3403af2
# ╠═f9e7b6a4-3e68-4adc-abeb-b9a6abbaf068
# ╠═922865aa-9d53-4231-b581-7a160a395420
# ╠═4d669db9-1661-4aa6-87f4-2977ae11e989
# ╠═39655593-9b8c-444e-8f7d-48cd92d7b9d8
# ╠═686250d4-235d-4699-88cc-5aaf1273b319
# ╟─4e46e5eb-74ca-400b-bca2-7329af1879ac
# ╠═8f5bed6b-726e-4b19-bd20-ec41151ce83d
# ╠═dac5f3f5-1327-4b65-be25-46764829e5b4
# ╠═7d0b670d-03de-40e3-9eef-6f38289d3230
# ╠═00350329-7af5-4723-b858-f412a414c67f
# ╠═f67161c6-dc5f-44a9-b4ae-937f122b76df
# ╠═a261cc22-30c0-4d07-b8de-758618ec608d
# ╟─689207e3-04f1-4521-8f87-ee8ebe5b09e7
# ╠═7e51c8ab-2e1a-4b0f-8ea3-9ae9399eb82b
# ╠═ce4b4f7f-69f0-4622-81fe-7b3f5c8089de
# ╟─7a9f8b03-3cb3-40a2-9dd8-7f4a5b95a018
# ╠═d969f868-81e5-494b-a528-ca86e07f8b2b
# ╠═2214b214-ff3b-448a-8f7e-a4a754f8a5a3
# ╟─a8d54a2b-f9ab-4f95-aa0e-12819cc35ea7
# ╠═13f2b66b-05e3-4ee2-a01a-fbd47d57233b
# ╟─cd2185f3-9f37-4371-9b94-833df7dd1cbd
# ╟─d8550c75-51e3-4d25-ad73-c4d72eac4267
# ╠═3b21afa0-c6b2-4531-974b-f4f82d8c7b40
# ╠═4dd137ec-45d7-4be3-934f-454adefb8203
# ╟─ebf5ee8e-430e-4634-8543-fbeec58fb0f8
# ╟─adea9d70-dcda-4792-b365-d15274855fd8
# ╠═8897fb75-8a19-48d6-8ec9-4e16fb36f9c2
# ╠═2f374d9f-54a4-48ce-8aa0-22ae407fd5b3
# ╠═1c37e89a-6a90-4926-bb37-c4f7f445c9f8
# ╠═599d46f2-c2e5-415f-9a05-60696ae57a52
# ╠═c13c9d9c-7b81-4a16-90dd-f56c69e15f23
# ╠═b373c92c-d465-4a0f-a925-ae50360c277e
# ╠═4eac8fc5-3d30-44b6-828e-e40f17bd1d7e
# ╠═668bfd06-7308-4cc8-83e2-34d51f9730fd
# ╠═f54b0681-3a3e-45cf-b274-4962989821dc
# ╠═789f8085-76c9-4fba-95d2-b6b6b603f868
# ╠═030af251-fe5b-46c8-b6d1-1b2efb7cd6c3
# ╠═65c91137-034b-43c5-a7d1-c25cc694b1d4
# ╟─7c07e06f-5c31-4453-8d0a-d89842ebe7d2
# ╟─e97005f1-f165-4dfb-8be0-4d078eb27930
# ╟─c6143a16-c8d9-48ec-8588-7f430851576b
# ╟─bd87d080-07fc-4c99-9e6c-1df8229cb46c
# ╟─c3e32361-b417-4305-acf3-5862813d2d0c
# ╟─19f116ba-bfd3-4579-9e45-b1afb61a9cbd
# ╟─5bc9a02e-62b1-48de-8287-33be309402fe
# ╟─e950d80d-8b1a-427c-a4f0-361d62f9a450
# ╟─600773b7-8b3b-4e30-91ca-3ce0407dae68
# ╠═6bef1367-bc12-40cd-8c6d-2d5934cdb9d6
# ╠═aa34de2c-dbf1-4551-aac7-94f361fa9b6a
# ╠═43d45ad1-0289-4992-a089-22c6b299931d
# ╠═f6677948-c480-4557-ac40-e38a76470daf
# ╠═fc86ff7a-bbab-4278-a1dd-5dfe3da87adf
# ╠═7d5fdf32-1619-4835-976c-8d0513efb24c
# ╠═0e504ae7-4165-4261-a02d-85141e2f5b5f
# ╠═10266ff0-7e64-41c9-b420-d6c15257ceef
# ╠═16991c8c-b7c3-495b-93ce-70b949017b3c
# ╠═f917c1c9-8bb8-4e35-bb92-24df9e2945a1
# ╠═9392b89e-6a64-4a7c-80dd-c6124820bdb3
# ╠═f518d526-34ca-4c15-8255-4ee79c8581dd
# ╠═2869dce2-0151-4ead-9ba4-39ac12876c95
# ╠═2e224d8b-6fb1-4b0b-aacb-28294f8e7b8d
# ╠═7290ff57-ead8-4518-8506-0f845cf14d0f
# ╠═365ae9ff-9b94-457a-a687-d1c8581fc5ab
# ╠═78e4ed32-9943-458a-aff8-3a9d0397e400
# ╠═78941ab1-25bf-4b10-935b-4deff67fe5dc
# ╠═512a39c1-1fd0-4d6b-9d45-e45457dfabfb
# ╠═fd778a85-b7f5-4d87-b671-4ac3efee4b79
# ╠═1739ae21-a60d-413e-832d-cb12fad0ca34
# ╠═7bc0108e-6842-4850-92a5-0cae5d7d4706
# ╠═2f750d6e-461f-4e6a-bc92-2939dbb55421
# ╠═5e26e413-612c-4728-9c08-15f45aa10e1f
# ╠═804fd07b-b4f0-4119-a708-60d4973e84eb
# ╠═59e2f85f-acdd-4ca5-a91c-78b9bb42364e
# ╠═4fe853a2-100b-4c25-8559-1ea008c0cad2
# ╠═3f7dad44-11bc-4535-b7cb-db5938f5d2d4
# ╠═98df485c-de5b-4ea9-ba12-7d0677a4e7dc
# ╠═2a0e8394-076a-4d0f-8458-b2f838c9bd77
# ╠═7d9ae1fa-047c-11ed-0af1-c1af45ffc663
