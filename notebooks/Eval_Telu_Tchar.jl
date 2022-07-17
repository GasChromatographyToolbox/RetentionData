### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

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

# ╔═╡ 00350329-7af5-4723-b858-f412a414c67f
# [x] filter out the triglycerides -> Tchar>700K
# fit the model with four parameters to every set of single combinations of Tinit, RT, pin resp. Tinit and rT for all stationary phases together

# ╔═╡ f67161c6-dc5f-44a9-b4ae-937f122b76df
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

# ╔═╡ a261cc22-30c0-4d07-b8de-758618ec608d
@. model(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*x)) + exp(p[1]*p[4]))

# ╔═╡ cd2185f3-9f37-4371-9b94-833df7dd1cbd
md"""
## Fit
"""

# ╔═╡ ba3efec4-ae7b-4389-bd3a-aa21fc3cfe7a
md"""
## Assumption

``T_1 \prop T_{init} + a₀ RT``

"""

# ╔═╡ 4fe853a2-100b-4c25-8559-1ea008c0cad2
function Rsquare(fit, y)
	sstot = sum((y.-mean(y)).^2)
	ssres = sum(fit.resid.^2)
	R2 = 1-ssres/sstot
	return R2
end

# ╔═╡ 7bc0108e-6842-4850-92a5-0cae5d7d4706
function fitting(TeluTchar)
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
	for i=1:length(Tinit)
		for j=1:length(RT)
			for k=1:length(pin)
				TeluTchar_f = filter([:Tinit, :RT, :pin] => (x,y,z) -> x == Tinit[i] && y == RT[j] && z == pin[k], TeluTchar)
				xdata = TeluTchar_f.Tchar
				ydata = TeluTchar_f.Telu
				p0 = [0.05, Tinit[i], 1.0, Tinit[i]+273.15]
				fits[i,j,k] = curve_fit(model, xdata, ydata, p0)
				sigmas = NaN.*ones(4)
				try
					sigmas = stderror(fits[i,j,k])
				catch
					sigmas = NaN.*ones(4)
				end
				push!(Tinits, Tinit[i])
				push!(RTs, RT[j])
				push!(pins, pin[k])
				push!(bs, fits[i,j,k].param[1] ± sigmas[1])
				push!(T₀s, fits[i,j,k].param[2] ± sigmas[2])
				push!(ms, fits[i,j,k].param[3] ± sigmas[3])
				push!(T₁s, fits[i,j,k].param[4] ± sigmas[4])
				push!(R²s, Rsquare(fits[i,j,k], ydata))
			end
		end
	end
	params = DataFrame(Tinit=Tinits, RT=RTs, pin=pins, b=bs, T₀=T₀s, m=ms, T₁=T₁s, R²=R²s)
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
						end
					end
				end
			end
		end
	end	
	TeluTchar_ = filter!([:Telu, :Tchar] => (x,y) -> !isnan(x) && y < 700.0, sort!(DataFrame(Name=Names, Tinit=Tinits, RT=RTs, pin=pins, sp=sps, Tchar=Tchars, Telu=Telus, kelu=kelus), [:Tchar])) # sort Tchar and filter out NaN values
	return TeluTchar_
end

# ╔═╡ fb743e94-7f59-45cf-b329-b980d3403af2
Telu_Tchar = TeluTchar(df)

# ╔═╡ dac5f3f5-1327-4b65-be25-46764829e5b4
# plot Telu(Tchar, Tinit) for fixed rT (resp. fixed RT and pin) and variation of sp
begin
	plotly()
	select_i_RT = 3
	select_i_pin = 4
	p___ = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	for i=1:length(sp)
		Telu_Tchar_f = filter([:sp, :RT, :pin] => (x,y,z) -> x == sp[i] && y == RT[select_i_RT] && z == pin[select_i_pin], Telu_Tchar)
		scatter!(p___, Telu_Tchar_f.Tchar, Telu_Tchar_f.Tinit, Telu_Tchar_f.Telu, label=sp[i])
	end
	p___
end

# ╔═╡ 4dd137ec-45d7-4be3-934f-454adefb8203
param, fit = fitting(Telu_Tchar);

# ╔═╡ d05208cc-443e-4764-a3af-00807cfacf4f
param

# ╔═╡ 88bbd9d5-d18a-4f61-af91-29b35c4b4e91
plot(1:length(param.R²), param.R²)

# ╔═╡ 1a8eac36-1ac5-4e36-84f0-fa554824bf28
plot(1:length(param.b), param.b)

# ╔═╡ e3cee02e-0b3d-41e3-9f5e-89dc0d5bb6b8
plot(1:length(param.T₀), Measurements.value.(param.T₀))

# ╔═╡ adea9d70-dcda-4792-b365-d15274855fd8
begin
	select_pin = 150.0
	
	p_b = plot(xlabel="Tinit in K", ylabel="b")
	for i=1:length(RT)
		param_f = filter([:RT, :pin] => (x,y) -> x == RT[i] && y == select_pin, param)
		plot!(p_b, param_f.Tinit.+273.15, Measurements.value.(param_f.b), label=RT[i], m=:circle)
	end
	p_T₀ = plot(xlabel="Tinit in °C", ylabel="T₀")
	for i=1:length(RT)
		param_f = filter([:RT, :pin] => (x,y) -> x == RT[i] && y == select_pin, param)
		plot!(p_T₀, param_f.Tinit, Measurements.value.(param_f.T₀), label=RT[i], m=:circle)
	end
	p_m = plot(xlabel="Tinit in K", ylabel="m")
	for i=1:length(RT)
		param_f = filter([:RT, :pin] => (x,y) -> x == RT[i] && y == select_pin, param)
		plot!(p_m, param_f.Tinit.+273.15, Measurements.value.(param_f.m), label=RT[i], m=:circle)
	end
	p_T₁ = plot(xlabel="Tinit in K", ylabel="T₁")
	for i=1:length(RT)
		param_f = filter([:RT, :pin] => (x,y) -> x == RT[i] && y == select_pin, param)
		plot!(p_T₁, param_f.Tinit.+273.15, Measurements.value.(param_f.T₁), label=RT[i], m=:circle)
	end
	plot(p_b, p_T₀, p_m, p_T₁)
end

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
						end
					end
			end
		end
	end	
	TeluTchar_ = filter!([:Telu, :Tchar] => (x,y) -> !isnan(x) && y < 700.0, sort!(DataFrame(Name=Names, Tinit=Tinits, rT=rTs, sp=sps, Tchar=Tchars, Telu=Telus, kelu=kelus), [:Tchar])) # sort Tchar and filter out NaN values
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
	select_i_rT = 10
	plotly()
	p__ = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	for i=1:length(sp)
		Telu_Tchar_rT_f = filter([:sp, :rT] => (x,y) -> x == sp[i] && y == rT[select_i_rT], Telu_Tchar_rT)
		scatter!(p__, Telu_Tchar_rT_f.Tchar, Telu_Tchar_rT_f.Tinit, Telu_Tchar_rT_f.Telu, label=sp[i])
	end
	p__
end

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
# ╠═4e46e5eb-74ca-400b-bca2-7329af1879ac
# ╠═dac5f3f5-1327-4b65-be25-46764829e5b4
# ╠═00350329-7af5-4723-b858-f412a414c67f
# ╠═f67161c6-dc5f-44a9-b4ae-937f122b76df
# ╠═a261cc22-30c0-4d07-b8de-758618ec608d
# ╠═cd2185f3-9f37-4371-9b94-833df7dd1cbd
# ╠═4dd137ec-45d7-4be3-934f-454adefb8203
# ╠═d05208cc-443e-4764-a3af-00807cfacf4f
# ╠═88bbd9d5-d18a-4f61-af91-29b35c4b4e91
# ╠═1a8eac36-1ac5-4e36-84f0-fa554824bf28
# ╠═e3cee02e-0b3d-41e3-9f5e-89dc0d5bb6b8
# ╟─adea9d70-dcda-4792-b365-d15274855fd8
# ╠═ba3efec4-ae7b-4389-bd3a-aa21fc3cfe7a
# ╠═7bc0108e-6842-4850-92a5-0cae5d7d4706
# ╠═59e2f85f-acdd-4ca5-a91c-78b9bb42364e
# ╠═4fe853a2-100b-4c25-8559-1ea008c0cad2
# ╠═98df485c-de5b-4ea9-ba12-7d0677a4e7dc
# ╠═2a0e8394-076a-4d0f-8458-b2f838c9bd77
# ╠═7d9ae1fa-047c-11ed-0af1-c1af45ffc663
