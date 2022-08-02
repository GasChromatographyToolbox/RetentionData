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

Load simulation results from the folder `/data/sims/Telu_Tchar/` (results from the script `/scripts/Telu_Tchar.jl`) and evaluate these for the correlation of ``T_{elu}`` with ``T_{char}`` and temperature program parameters ``T_{init}`` and ``r_T`` (resp. ``R_T`` and ``p_{in}``).
"""

# ╔═╡ 2a8dd6b6-dc9c-453b-9652-2226fc1f39ff
md"""
## Load data
"""

# ╔═╡ e5f30f41-3205-4746-969f-3d7c106e11bf
begin
	df = collect_results("/Users/janleppert/Documents/GitHub/ThermodynamicData/data/sims/Telu_Tchar/")
	df.dimless_rate = df.rate.*df.tMref./30.0./60.0
	df
end

# ╔═╡ 650a0f9b-77c8-4918-adf5-28fae3c7dd35
scatter(df.rate, df.dimless_rate, xlabel="heating rate in °C/min", ylabel="dimensionless heating rate")

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
bconst = 0.06

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
	for i=1:length(rate)
		param_f = filter([rate_key] => (x) -> x == rate[i], param)
		plot!(p_b, param_f.Tinit.+273.15, Measurements.value.(param_f.b), label=round(rate[i],sigdigits=3), m=:circle)

		plot!(p_T₀, param_f.Tinit, Measurements.value.(param_f.T₀), label=round(rate[i],sigdigits=3), m=:circle)
	
		plot!(p_m, param_f.Tinit.+273.15, Measurements.value.(param_f.m), label=round(rate[i],sigdigits=3), m=:circle)
	
		plot!(p_T₁, param_f.Tinit.+273.15, Measurements.value.(param_f.T₁), label=round(rate[i],sigdigits=3), m=:circle)
	end
	plot(p_T₀, p_T₁, p_b, p_m)
end

# ╔═╡ eb37f6b4-23b8-4163-8ec4-67324aba214a
function plot_parameters_rT(param; logrT=true)
	if logrT == true
		xlbl = "log(rT)"
	else
		xlbl = "rT"
	end
	p_T₀ = plot(xlabel=xlbl, ylabel="T₀")
	p_T₁ = plot(xlabel=xlbl, ylabel="T₁")
	p_b = plot(xlabel=xlbl, ylabel="b")
	p_m = plot(xlabel=xlbl, ylabel="m")
	rT = unique(param.rT)
	Tinit = unique(param.Tinit)
	for i=1:length(Tinit)
		param_f = filter([:Tinit] => (x) -> x == Tinit[i], param)
		if logrT == true
			x = log.(param_f.rT)
		else
			x = param_f.rT
		end
		plot!(p_b, x, Measurements.value.(param_f.b), label=round(Tinit[i],sigdigits=3), m=:circle)

		plot!(p_T₀, x, Measurements.value.(param_f.T₀), label=round(Tinit[i],sigdigits=3), m=:circle)
	
		plot!(p_m, x, Measurements.value.(param_f.m), label=round(Tinit[i],sigdigits=3), m=:circle)
	
		plot!(p_T₁, x, Measurements.value.(param_f.T₁), label=round(Tinit[i],sigdigits=3), m=:circle)
	end
	plot(p_T₀, p_T₁, p_b, p_m)
end

# ╔═╡ b5d6fb34-00ba-42c1-9d57-518bead48f4b
md"""
## Note

The parameter ``T_1`` is a linear function of ``T_{init}``:

``
T_1 = a_1 + b_1 T_{init}
``
"""

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

# ╔═╡ 5bc9a02e-62b1-48de-8287-33be309402fe
md"""
## Plot ``T_{elu}(T_{char}, θ_{char}, r_T, T_{init})``

rT: $(@bind select_i_rT__ Slider(1:length(rT)))

Tinit: $(@bind select_i_Tinit__ Slider(1:length(Tinit)))
"""

# ╔═╡ 6bef1367-bc12-40cd-8c6d-2d5934cdb9d6
md"""
## Fit ``r_T`` with filtered data
**For now: exclude Brehmer\_PCB, Duong\_FAME, Stultz\_Dioxins from the data**

**Also filter out following substances:**

- Pentanol, Rxi17SilMS	
- 2,3-Butanediol, Rxi5ms
- 2-Hexanone, Rxi17SilMS	
- Heptanol, Rxi17SilMS	

**limit the range of heating rates: 0.1<rT<1.0 (more typical values used in GC)**

**limit the ``T_{init}>200``°C**
"""

# ╔═╡ 8a730b2f-a51f-40ac-a3a9-054f0c903031
md"""
### Fit model4
"""

# ╔═╡ 16991c8c-b7c3-495b-93ce-70b949017b3c
md"""
### Fit ``m=1``
"""

# ╔═╡ 3844df54-1622-4313-92ee-36c6599f80d0
md"""
### Fit ``b=0.06``
"""

# ╔═╡ 69113eb5-cca4-40e7-8a5e-96e196a99fc0
md"""
### Fit ``m=1`` and ``b=0.06``
"""

# ╔═╡ aec7dae9-c73a-4cce-a440-2724799085e1
log(0.6)

# ╔═╡ e8e0b85f-286e-4665-a992-3f26f8eb016e
md"""
## Assumptions for ``T_0`` and ``T_1``

The constant linear branch ``T_1`` obviously depends directly on ``T_{init}``. For low and moderate heating rates, it linearly depend on ``rT``.

Assumption:

a) ``T_1 = a_1 + m_1 T_{init}``

b) ``T_1 = a_1 + m_1 T_{init} + n_1 r_T``

The intercept of the linear increasing branch, ``T_0`` depends on heating rate ``r_T``. A dependency on ``T_{init}`` is possible but not clear. For low to moderate heating rates ``T_0`` is a linear function of ``\log{r_T}``.

Assumption:

1) ``T_0 = a_0 + m_0 log(r_T)``

2) ``T_0 = a_0 + m_0 log(r_T) + n_0 r_T``

3) ``T_0 = a_0 + n_0 r_T``

"""

# ╔═╡ 34b5b40f-3be2-4b9f-ad4c-f4d99dac1299
md"""
## Multivariate fit ``T_{elu}(T_{char}, r_T, T_{init})``
"""

# ╔═╡ 78c48193-7021-4803-bf4b-bf4e03406c47
md"""
### 1a) model\_mv\_1a
"""

# ╔═╡ 494c2b8c-04b3-46ff-9c28-cc9b79f9a7d4
@. model_mv_1a(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*log(x[:,2]) + p[4]*x[:,1])) + exp(p[1]*(p[5] + p[6]*x[:,3])))

# ╔═╡ ac33cbec-938b-472c-863f-97b8c3d6093d
p0_1a = [0.05, -40.0, 16.0, 1.0, 0.0, 1.5]

# ╔═╡ ec7bb02e-b4cd-47aa-a694-07c14c0b5e86
md"""
### 1b) model\_mv\_1b 👍
"""

# ╔═╡ bd88682c-a94d-431b-a5dc-dbc7638a3e48
@. model_mv_1b(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*log(x[:,2]) + p[4]*x[:,1])) + exp(p[1]*(p[5] + p[6]*x[:,3] + p[7]*x[:,2])))

# ╔═╡ 47bf2c03-bfa0-4ecc-ace3-a64959c558d1
p0_1b = [0.05, -40.0, 16.0, 1.0, 0.0, 1.5, 1.5]

# ╔═╡ e0a8ba27-65c0-4425-a1f2-1f30fe768c30
md"""
### 1b0) model\_mv\_1b0, ``m=1``, ``m_1=1``, ``a_1=273.15``, ``a_0=0``
"""

# ╔═╡ 241a52ee-c094-4198-b964-67847ec313ec
@. model_mv_1b0(x,p) = 1/p[1] * log(exp(p[1]*(p[2]*log(x[:,2]) + x[:,1])) + exp(p[1]*(273.15 + x[:,3] + p[3]*x[:,2])))

# ╔═╡ 13e14de4-09d0-4013-8176-ac96c678a0a6
p0_1b0 = [0.05, 16.0, 1.5]

# ╔═╡ 59ae1c1b-a047-4bcb-b1a1-9f6e59d36695
md"""
### 1b1) model\_mv\_1b 👍, ``m=1``, ``m_1=1``, ``a_1=273.15``
"""

# ╔═╡ e0196c9b-ce1f-4417-b8e7-f1c9a942efb2
@. model_mv_1b1(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*log(x[:,2]) + x[:,1])) + exp(p[1]*(273.15 + x[:,3] + p[4]*x[:,2])))

# ╔═╡ 33efc8a9-b2ad-4cb9-8d90-245e3463c393
p0_1b1 = [0.05, -40.0, 16.0, 1.5]

# ╔═╡ 40805b6a-c69b-48b8-aefd-49a731c06de7
md"""
### 2a) model\_mv\_2a

Use of all heating rates.
"""

# ╔═╡ 3676b513-c1b6-46e9-b069-b14be86a78eb
@. model_mv_2a(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*log(x[:,2]) + p[4]*x[:,2] + p[5]*x[:,1])) + exp(p[1]*(p[6] + p[7]*x[:,3])))

# ╔═╡ 3793b682-2b9f-410a-86af-53791b42979d
p0_2a = [0.05, -40.0, 16.0, 20.0, 1.0, 0.0, 1.5]

# ╔═╡ 78935736-f4b1-41c5-bd43-328e41eb53e8
md"""
### 2b) model\_mv\_2b 👍
"""

# ╔═╡ fb0c53bd-7e57-417f-a2eb-584896eaf6b3
@. model_mv_2b(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*log(x[:,2]) + p[4]*x[:,2] + p[5]*x[:,1])) + exp(p[1]*(p[6] + p[7]*x[:,3] + p[8]*x[:,2])))

# ╔═╡ 01cd2d11-2e63-4bba-b7e5-84404cb55f0f
p0_2b = [0.05, -40.0, 16.0, 20.0, 1.0, 0.0, 1.5, 1.0]

# ╔═╡ d111ec3e-62d3-4737-bf56-c62eb3d034ac
md"""
### 2bmod) model\_mv\_2bmod
"""

# ╔═╡ d9668965-1e6f-4088-a686-2657922aa76e
@. model_mv_2bmod(x,p) = 1/p[1] * log(x[:,2]^p[3]*exp(p[1]*(p[2] + p[4]*x[:,2] + p[5]*x[:,1])) + exp(p[1]*(p[6] + p[7]*x[:,3] + p[8]*x[:,2])))

# ╔═╡ deba9143-50fa-4bf3-b7c7-a7a93c404f4c
p0_2bmod = [0.05, -40.0, 1.0, 20.0, 1.0, 0.0, 1.5, 1.0]

# ╔═╡ 40d8b557-65cf-4397-b429-2d020ba2c8c8
md"""
### 2b0) model\_mv\_2b0, ``m=1``, ``m_1=1``, ``a_1=273.15``, ``a_0=0.0``
"""

# ╔═╡ 62663709-8eed-4777-a121-61f5ca5c9973
@. model_mv_2b0(x,p) = 1/p[1] * log(exp(p[1]*(p[2]*log(x[:,2]) + p[3]*x[:,2] + 1.0*x[:,1])) + exp(p[1]*(273.15 + 1.0*x[:,3] + p[4]*x[:,2])))

# ╔═╡ ed9d99ee-01f3-4f59-9d26-5fe0841e0a40
p0_2b0 = [0.05, 16.0, 20.0, 1.0]

# ╔═╡ 983dfe18-83fb-4eff-9da6-e6e4ca79d506
md"""
### 2b1) model\_mv\_2b1 👍, ``m=1``, ``m_1=1``, ``a_1=273.15``
"""

# ╔═╡ 445c8480-19dd-4a5a-a022-97a758341a48
@. model_mv_2b1(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*log(x[:,2]) + p[4]*x[:,2] + 1.0*x[:,1])) + exp(p[1]*(273.15 + 1.0*x[:,3] + p[5]*x[:,2])))

# ╔═╡ 9eb0565c-b76f-4cdc-952d-c6139e4e6f63
p0_2b1 = [0.05, -40.0, 16.0, 20.0, 1.0]

# ╔═╡ 3286b6ee-2deb-4c2b-b1aa-3284b9697139
md"""
### 3a) model\_mv\_3a
"""

# ╔═╡ e36fbd07-42d2-44a3-8bcd-a80e19df18b1
@. model_mv_3a(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*x[:,2] + p[4]*x[:,1])) + exp(p[1]*(p[5] + p[6]*x[:,3])))

# ╔═╡ 840ec2d9-5389-4ba2-9a7b-9346763f93dc
p0_3a = [0.05, -40.0, 20.0, 1.0, 273.15, 1.0]

# ╔═╡ dd1670fe-9968-42a5-8f7a-3fb8c72c197e
md"""
### 3b) model\_mv\_3b
"""

# ╔═╡ f60a676a-3458-4681-9634-dabaede63383
@. model_mv_3b(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*x[:,2] + p[4]*x[:,1])) + exp(p[1]*(p[5] + p[6]*x[:,3] + p[7]*x[:,2])))

# ╔═╡ 2f6bfa25-acd7-41cf-8429-29c9df694a8d
p0_3b = [0.05, -40.0, 20.0, 1.0, 273.15, 1.0, 1.0]

# ╔═╡ 084f8e7f-ca98-4727-9346-2315b5f51fc2
md"""
### 3b1) model\_mv\_3b1, ``m=1``, ``m_1=1``, ``a_1=273.15``
"""

# ╔═╡ f00ff9f1-bb8a-4e66-a842-479a9b68e4e2
@. model_mv_3b1(x,p) = 1/p[1] * log(exp(p[1]*(p[2] + p[3]*x[:,2] + 1.0*x[:,1])) + exp(p[1]*(273.15 + 1.0*x[:,3] + p[4]*x[:,2])))

# ╔═╡ e66578f9-c40e-448f-aba7-12ff5e7738d8
p0_3b1 = [0.05, -40.0, 20.0, 1.0]

# ╔═╡ 7871ce05-73de-4b88-9b00-9a8a735ef3da
md"""
## Comparison
"""

# ╔═╡ 8e1c640d-75bd-46d1-a5cf-ec63f493cf10
md"""
Best models:

Model 1b:

``
\frac{1}{b}\left(\exp{\left(b\left(a_0 + m_0 log(r_T) + m T_{char}\right)\right)} + \exp{\left(b\left(a_1 + m_1 T_{init} + n_1 r_T\right)\right)}\right)
``

Model 2b:

``
\frac{1}{b}\left(\exp{\left(b\left(a_0 + m_0 log(r_T) + n_0 r_T + m T_{char}\right)\right)} + \exp{\left(b\left(a_1 + m_1 T_{init} + n_1 r_T\right)\right)}\right)
``

Some Parameters could be fixed:

``m=1``

``a_1=273.15`` (or ``a_1=0`` and ``T_{init}`` in K)

``m_1=1``

-> model 1b1 and 2b1 (3b1)

It would make sense to set 

``a_0 = 0``

-> models 1b0 and 2b0

Model 1b0 not good

**Model 2b0 looks good**.
"""

# ╔═╡ c31eee97-b94d-463a-bb91-41f18cd83742
ff(x) = 1 + 2*x

# ╔═╡ 6bbe03d7-cb42-4ff5-92a8-6441515ba8cc
gg(x) = 10 -0.6x

# ╔═╡ 4a14a4ba-162d-46bf-98c0-46107fff0b45
cc(x,k) = 1/k*log(exp(k*ff(x))+exp(k*gg(x))) 

# ╔═╡ 98eb4f63-c4fc-465f-928c-f7f805d49839
dd(x,k) = 1/k*log(exp(k*ff(x))-exp(k*gg(x))) 

# ╔═╡ b86e3b4d-b216-49b9-9e5d-73d57f24c72e
x_ = 3.5:0.1:9

# ╔═╡ ffb0843d-dcaf-4207-bd76-e1e21ebdc2cd
begin
	plot(x_, ff.(x_))
	plot!(x_, gg.(x_))
	plot!(x_, cc.(x_, -0.5))
	plot!(x_, cc.(x_, 0.5))
	plot!(x_, dd.(x_, +0.5))
end

# ╔═╡ 2e224d8b-6fb1-4b0b-aacb-28294f8e7b8d
md"""
# End
"""

# ╔═╡ 2c28f8ad-ab76-4131-8856-55a43b31fde6
function plot_fit_mv(TeluTchar, fit, model, i_Tinit, i_rT) 
	Tinit = unique(TeluTchar.Tinit)
	rT = unique(TeluTchar.rT)
	TeluTchar_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[i_Tinit] && y == rT[i_rT], TeluTchar)
	#param_f = filter([:Tinit, :rT] => (x,y) -> x == Tinit[i_Tinit] && y == rT[i_rT], param)
	p_fit = scatter(TeluTchar_f.Tchar, TeluTchar_f.Telu, label="sim", title=string("fit Tinit=", Tinit[i_Tinit], "°C, rT=", round(rT[i_rT];sigdigits=3)))
	plot!(p_fit, TeluTchar_f.Tchar, model([TeluTchar_f.Tchar TeluTchar_f.rT TeluTchar_f.Tinit], fit.param), label="mv_fit")
	#T₀ = param.T₀[1]
	#T₁ = param.T₁[1]
	#b = param.b[1]
	#m = param.m[1]
	return p_fit#, T₀, T₁, b, m  
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
		for j=1:length(rT)
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
	p_fit_RT_pin = plot_fit(Telu_Tchar, param_RT_pin, fit_RT_pin, select_model, show_i_Tinit, show_i_RT, show_i_pin) 
	md"""
	#### Fit result
	
	``T_0 =`` $(p_fit_RT_pin[2]) K

	``T_1 =`` $(p_fit_RT_pin[3]) K

	``b =`` $(p_fit_RT_pin[4])
	
	``m =`` $(p_fit_RT_pin[5]) K⁻¹

	$(embed_display(p_fit_RT_pin[1]))
	"""
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

	For high dimensionless heating rates, ``r_T>2``, certain solutes on certain stationary phases (Rxi17SilMS, Rxi5ms, SPB50) show a to high ``T_{elu}`` value compared to other solutes, especially in the region, where the solutes should elute with nearly ``T_{init}``..
	"""
end

# ╔═╡ 2f374d9f-54a4-48ce-8aa0-22ae407fd5b3
param_rT, fit_rT = fitting(Telu_Tchar_rT, select_model_);

# ╔═╡ 1c37e89a-6a90-4926-bb37-c4f7f445c9f8
plot_parameters_i(filter([:converged] => x -> x == true, param_rT))

# ╔═╡ b373c92c-d465-4a0f-a925-ae50360c277e
plot_parameters_Tinit(filter([:converged] => x -> x == true, param_rT))

# ╔═╡ 0fd2139b-3e50-48ad-9d76-6166b91d4a5b
plot_parameters_rT(filter([:converged] => x -> x == true, param_rT))

# ╔═╡ 7c07e06f-5c31-4453-8d0a-d89842ebe7d2
begin
	p_fit_rT = plot_fit(Telu_Tchar_rT, param_rT, fit_rT, select_model_, show_i_Tinit_, show_i_rT_) 
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
	filter!([:Name] => (x) -> !occursin("Pentanol",x) && !occursin("2,3-Butanediol",x) && !occursin("2-Hexanone",x) && !occursin("Heptanol",x), Telu_Tchar_rT_data_excluded)
	filter!([:rT] => x -> x > 0.1 && x < 1.0, Telu_Tchar_rT_data_excluded)
	filter!([:Tinit] => x -> x < 201.0, Telu_Tchar_rT_data_excluded)
	rT_f = sort!(unique(Telu_Tchar_rT_data_excluded.rT))
	Tinit_f = sort!(unique(Telu_Tchar_rT_data_excluded.Tinit))
end

# ╔═╡ 43d45ad1-0289-4992-a089-22c6b299931d
param_rT_data_excluded, fit_rT_data_excluded = fitting_rT(Telu_Tchar_rT_data_excluded, model4);

# ╔═╡ f6677948-c480-4557-ac40-e38a76470daf
param_rT_data_excluded

# ╔═╡ 58e58468-ee89-4f00-bed9-bbfd2f494be1
plot_parameters_Tinit(filter([:converged] => x -> x == true, param_rT_data_excluded))

# ╔═╡ b552ac0c-7f39-4fd6-bd37-13aae5a81d99
plot_parameters_rT(filter([:converged] => x -> x == true, param_rT_data_excluded))

# ╔═╡ 0e504ae7-4165-4261-a02d-85141e2f5b5f
md"""
``T_{init}``: $(@bind show_i_Tinit_1 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_1 Slider(1:length(rT_f)))
"""

# ╔═╡ 10266ff0-7e64-41c9-b420-d6c15257ceef
begin
	p_fit_rT_1 = plot_fit(Telu_Tchar_rT_data_excluded, param_rT_data_excluded, fit_rT_data_excluded, model4, show_i_Tinit_1, show_i_rT_1) 
	md"""
	#### Fit result
	``b =`` $(p_fit_rT_1[2])
	
	``T_0 =`` $(p_fit_rT_1[3]) K
	
	``m =`` $(p_fit_rT_1[4]) K⁻¹
	
	``T_1 =`` $(p_fit_rT_1[5]) K

	$(embed_display(p_fit_rT_1[1]))
	"""
end

# ╔═╡ f917c1c9-8bb8-4e35-bb92-24df9e2945a1
param_rT_data_excluded_m1, fit_rT_data_excluded_m1 = fitting_rT(Telu_Tchar_rT_data_excluded, model_m1);

# ╔═╡ 03121d84-5b1b-41b7-befa-c629b08e6f52
plot_parameters_Tinit(filter([:converged] => x -> x == true, param_rT_data_excluded_m1))

# ╔═╡ 101ffb81-ada0-48fa-9df8-272800af2565
plot_parameters_rT(filter([:converged] => x -> x == true, param_rT_data_excluded_m1))

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

# ╔═╡ be2237d4-be73-4368-9f38-88280732f225
md"""
``T_{init}``: $(@bind show_i_Tinit_2 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_2 Slider(1:length(rT_f)))
"""

# ╔═╡ c5d95d97-8903-4dbf-ad12-359d2a235e6c
begin
	p_fit_rT_2 = plot_fit(Telu_Tchar_rT_data_excluded, param_rT_data_excluded_m1, fit_rT_data_excluded_m1, model_m1, show_i_Tinit_2, show_i_rT_2) 
	md"""
	#### Fit result
	``b =`` $(p_fit_rT_2[2])
	
	``T_0 =`` $(p_fit_rT_2[3]) K
	
	``m =`` $(p_fit_rT_2[4]) K⁻¹
	
	``T_1 =`` $(p_fit_rT_2[5]) K

	$(embed_display(p_fit_rT_2[1]))
	"""
end

# ╔═╡ 4071fa66-095c-49e7-adfe-44c1fce9c59a
param_rT_data_excluded_b, fit_rT_data_excluded_b = fitting_rT(Telu_Tchar_rT_data_excluded, model_b);

# ╔═╡ a6b9948f-7364-4294-a744-fd13550a4796
param_rT_data_excluded_b

# ╔═╡ 69c0a9bc-4606-4e90-aed1-6ea8dbc6f039
plot_parameters_Tinit(filter([:converged] => x -> x == true, param_rT_data_excluded_b))

# ╔═╡ ca7c60d2-04ad-4980-b26f-c02c70ef84a8
plot_parameters_rT(filter([:converged] => x -> x == true, param_rT_data_excluded_b))

# ╔═╡ fbde7fa1-c584-4732-ba60-05837ee62d8f
md"""
``T_{init}``: $(@bind show_i_Tinit_3 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_3 Slider(1:length(rT_f)))
"""

# ╔═╡ 23ccbdb4-2c7c-499f-af74-a093304d60cd
begin
	p_fit_rT_3 = plot_fit(Telu_Tchar_rT_data_excluded, param_rT_data_excluded_b, fit_rT_data_excluded_b, model_b, show_i_Tinit_3, show_i_rT_3) 
	md"""
	#### Fit result
	``b =`` $(p_fit_rT_3[2])
	
	``T_0 =`` $(p_fit_rT_3[3]) K
	
	``m =`` $(p_fit_rT_3[4]) K⁻¹
	
	``T_1 =`` $(p_fit_rT_3[5]) K

	$(embed_display(p_fit_rT_3[1]))
	"""
end

# ╔═╡ 388a86d4-766f-4937-9a7a-4bd583b1c1b5
param_rT_data_excluded_b_m1, fit_rT_data_excluded_b_m1 = fitting_rT(Telu_Tchar_rT_data_excluded, model_b_m1);

# ╔═╡ 85696436-f843-4a2a-a582-a3d27c29d951
param_rT_data_excluded_b_m1

# ╔═╡ 59e7835e-e9f2-4888-b862-28ba1a019ab6
plot_parameters_Tinit(filter([:converged] => x -> x == true, param_rT_data_excluded_b_m1))

# ╔═╡ ed71d956-bc13-407c-8aab-02d700b711b7
plot_parameters_rT(filter([:converged] => x -> x == true, param_rT_data_excluded_b_m1))

# ╔═╡ 1a8865f8-db21-423f-8591-033f41197d2e
plot_parameters_rT(filter([:converged] => x -> x == true, param_rT_data_excluded_b_m1), logrT=false)

# ╔═╡ 32782159-a1ee-4570-bdbc-31a4506ad902
md"""
``T_{init}``: $(@bind show_i_Tinit_4 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_4 Slider(1:length(rT_f)))
"""

# ╔═╡ 997409e6-1292-4bc0-94f7-b1439d9beeaa
begin
	p_fit_rT_4 = plot_fit(Telu_Tchar_rT_data_excluded, param_rT_data_excluded_b_m1, fit_rT_data_excluded_b_m1, model_b_m1, show_i_Tinit_4, show_i_rT_4) 
	md"""
	#### Fit result
	``b =`` $(p_fit_rT_4[2])
	
	``T_0 =`` $(p_fit_rT_4[3]) K
	
	``m =`` $(p_fit_rT_4[4]) K⁻¹
	
	``T_1 =`` $(p_fit_rT_4[5]) K

	$(embed_display(p_fit_rT_4[1]))
	"""
end

# ╔═╡ 42411a67-69f4-4088-8b01-2815f223e4f9
xdata = [Telu_Tchar_rT_data_excluded.Tchar Telu_Tchar_rT_data_excluded.rT Telu_Tchar_rT_data_excluded.Tinit]

# ╔═╡ 76873708-0cd3-4048-a33b-29297fb58813
ydata = Telu_Tchar_rT_data_excluded.Telu

# ╔═╡ 8c9a6c3a-78ac-4a60-af68-7306dbb79b7b
par_mv_1a, R2_mv_1a, fit_mv_1a = fit_data(model_mv_1a, xdata, ydata, p0_1a)

# ╔═╡ 137174bf-fde5-436a-ab76-92168bca291f
par_mv_1a

# ╔═╡ 5fe03d15-d5cc-4c3b-af2a-4ed35b40bfad
R2_mv_1a

# ╔═╡ 0c741a3f-5f76-4c86-be39-e5a02c2fe7e1
par_mv_1b, R2_mv_1b, fit_mv_1b = fit_data(model_mv_1b, xdata, ydata, p0_1b)

# ╔═╡ 1d3474df-ede5-427b-8bf8-c50ecdfaedb5
par_mv_1b

# ╔═╡ 073912f2-5dcd-4c70-810d-9196709ee2d5
R2_mv_1b

# ╔═╡ e809c0f8-b939-413d-b6b8-1c76524c7d89
par_mv_1b0, R2_mv_1b0, fit_mv_1b0 = fit_data(model_mv_1b0, xdata, ydata, p0_1b0)

# ╔═╡ 7c3d5c74-3710-4671-891a-f2f75d23ee7e
par_mv_1b0

# ╔═╡ e0edb329-ba6b-4338-9386-6a0d34ff7c61
R2_mv_1b0

# ╔═╡ 5eab4c51-789d-4f78-825b-9514f948cde5
par_mv_1b1, R2_mv_1b1, fit_mv_1b1 = fit_data(model_mv_1b1, xdata, ydata, p0_1b1)

# ╔═╡ accf6e87-0b2c-4251-bff0-941d84053733
par_mv_1b1

# ╔═╡ 86e15050-d13f-4c6f-99d0-5ad2fc011937
R2_mv_1b1

# ╔═╡ d91ec9d0-d9db-47b4-973f-c9565825dc71
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_1a Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_1a Slider(1:length(rT_f)))
"""

# ╔═╡ ef742f31-4a9c-4d49-ad12-e3a81fa881cb
p_fit_mv_1a = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_1a, model_mv_1a, show_i_Tinit_mv_1a, show_i_rT_mv_1a)

# ╔═╡ bee6bac9-c137-4979-83e2-23462660d9bf
begin
	plotly()
	p_mv_1a_Tinit = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	Telu_Tchar_f = filter([:rT] => x -> x == rT[show_i_rT_mv_1a], Telu_Tchar_rT_data_excluded)
	scatter!(p_mv_1a_Tinit, Telu_Tchar_f.Tchar, Telu_Tchar_f.Tinit, Telu_Tchar_f.Telu, label=round(rT[show_i_rT_mv_1a];sigdigits=3))
	xTchar = sort!(unique(Telu_Tchar_f.Tchar))
	yTinit = sort!(unique(Telu_Tchar_f.Tinit))
	zdata = Array{Float64}(undef, length(xTchar), length(yTinit))
	for i=1:length(xTchar)
		for j=1:length(yTinit)
			zdata[i,j] = model_mv_1a([xTchar[i] rT[show_i_rT_mv_1a] yTinit[j]], fit_mv_1a.param)[1]
		end
	end
	plot!(p_mv_1a_Tinit, xTchar, yTinit, zdata', st=:surface)
end

# ╔═╡ 5c8678b9-d0f3-4db3-8f0d-42ff43cd511c
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_1b Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_1b Slider(1:length(rT_f)))
"""

# ╔═╡ bf37a7dc-2870-489c-9ff9-3b2423f4870c
p_fit_mv_1b = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_1b, model_mv_1b, show_i_Tinit_mv_1b, show_i_rT_mv_1b)

# ╔═╡ d5d09e6d-667c-4d5d-aaf3-0f4841c49ca1
begin
	plotly()
	p_mv_1b_Tinit = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	Telu_Tchar_ff = filter([:rT] => x -> x == rT[show_i_rT_mv_1b], Telu_Tchar_rT_data_excluded)
	scatter!(p_mv_1b_Tinit, Telu_Tchar_ff.Tchar, Telu_Tchar_ff.Tinit, Telu_Tchar_ff.Telu, label=round(rT[show_i_rT_mv_1b];sigdigits=3))
	xTchar_1b = sort!(unique(Telu_Tchar_ff.Tchar))
	yTinit_1b = sort!(unique(Telu_Tchar_ff.Tinit))
	zdata_1b = Array{Float64}(undef, length(xTchar_1b), length(yTinit_1b))
	for i=1:length(xTchar_1b)
		for j=1:length(yTinit_1b)
			zdata_1b[i,j] = model_mv_1b([xTchar_1b[i] rT[show_i_rT_mv_1b] yTinit_1b[j]], fit_mv_1b.param)[1]
		end
	end
	plot!(p_mv_1b_Tinit, xTchar_1b, yTinit_1b, zdata_1b', st=:surface)
end

# ╔═╡ 4b6b6c34-d103-4483-a483-00162a5f2583
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_1b0 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_1b0 Slider(1:length(rT_f)))
"""

# ╔═╡ f7f053fa-f0d9-4a40-bf31-3edc12f1f9df
p_fit_mv_1b0 = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_1b0, model_mv_1b0, show_i_Tinit_mv_1b0, show_i_rT_mv_1b0)

# ╔═╡ 351892bd-40ff-4162-b6be-222188f15eab
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_1b1 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_1b1 Slider(1:length(rT_f)))
"""

# ╔═╡ 4ebe8f12-c38e-48cc-9efb-562763684dc0
p_fit_mv_1b1 = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_1b1, model_mv_1b1, show_i_Tinit_mv_1b1, show_i_rT_mv_1b1)

# ╔═╡ 5a501229-fe26-42c5-85ad-366baddfbaee
xdata_2a = [Telu_Tchar_rT_data_excluded.Tchar Telu_Tchar_rT_data_excluded.rT Telu_Tchar_rT_data_excluded.Tinit]

# ╔═╡ 26706952-0716-42c7-804a-738801be929e
ydata_2a = Telu_Tchar_rT_data_excluded.Telu

# ╔═╡ 5630427e-22b4-4437-963c-a8c49ee76d57
par_mv_2a, R2_mv_2a, fit_mv_2a = fit_data(model_mv_2a, xdata_2a, ydata_2a, p0_2a)

# ╔═╡ 3c61bc26-c4cd-4cad-b405-272b0a962e49
par_mv_2a

# ╔═╡ 69a69280-f3b9-49ea-98fd-b439b495428a
R2_mv_2a

# ╔═╡ 96222088-76e8-4c36-9a04-fce4ee4a6c94
par_mv_2b, R2_mv_2b, fit_mv_2b = fit_data(model_mv_2b, xdata_2a, ydata_2a, p0_2b)

# ╔═╡ ddcfb67c-3bcf-4320-8603-5db5e42cccf1
par_mv_2b

# ╔═╡ d1a18fb3-2114-401c-91df-d2ce348aa1e1
R2_mv_2b

# ╔═╡ 19d9940d-7b58-4f22-9e28-55a51cc3dc56
par_mv_2b

# ╔═╡ 8b89812b-de28-4a4a-8d66-4979bfb41e27
R2_mv_2b

# ╔═╡ 5c239244-fdae-4b17-b6af-09054c3c26d1
par_mv_2bmod, R2_mv_2bmod, fit_mv_2bmod = fit_data(model_mv_2bmod, xdata_2a, ydata_2a, p0_2bmod)

# ╔═╡ b0bc5c1f-c8d7-4b72-b302-226c83b1001f
par_mv_2bmod

# ╔═╡ bbf1b7de-8d8f-49ba-8c88-3cd85d403560
R2_mv_2bmod

# ╔═╡ 525a9093-ce35-4f43-95fb-c2729cdb02b4
par_mv_2b0, R2_mv_2b0, fit_mv_2b0 = fit_data(model_mv_2b0, xdata_2a, ydata_2a, p0_2b0)

# ╔═╡ cbf3837e-dbda-4c20-83d5-74b936877a14
par_mv_2b0

# ╔═╡ 3a4c953f-a881-4d0d-b12a-015c39faf0bb
R2_mv_2b0

# ╔═╡ 146a2948-5a80-4e49-a3b0-615e7bd52769
par_mv_2b1, R2_mv_2b1, fit_mv_2b1 = fit_data(model_mv_2b1, xdata_2a, ydata_2a, p0_2b1)

# ╔═╡ c123731b-f488-485c-bf5f-d46ba56a9f9e
par_mv_2b1

# ╔═╡ 1c23b58a-4f4f-426c-ba4e-21e04858fdfb
R2_mv_2b1

# ╔═╡ 05a1c3cb-c778-4d9d-be80-88c1136320d9
DataFrame(model = ["mv_1b", "mv_1b0", "mv_1b1", "mv_2b", "mv_2b0", "mv_2b1"], 
	R2 = [R2_mv_1b, R2_mv_1b0, R2_mv_1b1, R2_mv_2b, R2_mv_2b0, R2_mv_2b1], 
	b = [par_mv_1b[1], par_mv_1b0[1], par_mv_1b1[1], par_mv_2b[1], par_mv_2b0[1], par_mv_2b1[1]], 
	a₀ = [par_mv_1b[2], 0.0, par_mv_1b1[2], par_mv_2b[2], 0.0, par_mv_2b1[2]], 
	m₀ = [par_mv_1b[3], par_mv_1b0[2], par_mv_1b1[3], par_mv_2b[3], par_mv_2b0[2], par_mv_2b1[3]], 
	n₀ = [0.0, 0.0, 0.0, par_mv_2b[4], par_mv_2b0[3], par_mv_2b1[4]],
	m = [par_mv_1b[4], 1.0, 1.0, par_mv_2b[5], 1.0, 1.0],
	a₁ = [par_mv_1b[5], 273.15, 273.15, par_mv_2b[6], 273.15, 273.15],
	m₁ = [par_mv_1b[6], 1.0, 1.0, par_mv_2b[7], 1.0, 1.0],
	n₁ = [par_mv_1b[7], par_mv_1b0[3], par_mv_1b1[4], par_mv_2b[8], par_mv_2b0[4], par_mv_2b1[5]])

# ╔═╡ 5a16fa66-5ddf-4f18-b398-5b47c881e3b3
par_mv_3a, R2_mv_3a, fit_mv_3a = fit_data(model_mv_3a, xdata_2a, ydata_2a, p0_3a)

# ╔═╡ 5f89244c-45f8-4115-bee7-2cc9ee8ef903
par_mv_3a

# ╔═╡ caa4e339-87c8-4ffa-94ae-1617f562ce0b
R2_mv_3a

# ╔═╡ 13f1725c-9355-410c-aa25-5f4de21e0bf4
par_mv_3b, R2_mv_3b, fit_mv_3b = fit_data(model_mv_3b, xdata_2a, ydata_2a, p0_3b)

# ╔═╡ 0013f6f5-f547-49a0-917b-89ee473f9439
par_mv_3b

# ╔═╡ f6b370b5-4944-44e3-a9a7-31e10798399e
R2_mv_3b

# ╔═╡ 1381f354-84ab-4496-9200-d99b10d98ef3
DataFrame(model = ["mv_1a", "mv_1b", "mv_2a", "mv_2b", "mv_3a", "mv_3b"], 
	R2 = [R2_mv_1a, R2_mv_1b, R2_mv_2a, R2_mv_2b, R2_mv_3a, R2_mv_3b], 
	b = [par_mv_1a[1], par_mv_1b[1], par_mv_2a[1], par_mv_2b[1], par_mv_3a[1], par_mv_3b[1]], 
	a₀ = [par_mv_1a[2], par_mv_1b[2], par_mv_2a[2], par_mv_2b[2], par_mv_3a[2], par_mv_3b[2]], 
	m₀ = [par_mv_1a[3], par_mv_1b[3], par_mv_2a[3], par_mv_2b[3], NaN, NaN], 
	n₀ = [NaN, NaN, par_mv_2a[4], par_mv_2b[4], par_mv_3a[3], par_mv_3b[3]],
	m = [par_mv_1a[4], par_mv_1b[4], par_mv_2a[5], par_mv_2b[5], par_mv_3a[4], par_mv_3b[4]],
	a₁ = [par_mv_1a[5], par_mv_1b[5], par_mv_2a[6], par_mv_2b[6], par_mv_3a[5], par_mv_3b[5]],
	m₁ = [par_mv_1a[6], par_mv_1b[6], par_mv_2a[7], par_mv_2b[7], par_mv_3a[6], par_mv_3b[6]],
	n₁ = [NaN, par_mv_1b[7], NaN, par_mv_2b[8], NaN, par_mv_3b[7]])

# ╔═╡ 7239b6aa-86fc-42c2-ad54-2fcee51f3d81
par_mv_3b1, R2_mv_3b1, fit_mv_3b1 = fit_data(model_mv_3b1, xdata_2a, ydata_2a, p0_3b1)

# ╔═╡ c1c961e9-41ea-43aa-bde5-0d94ba277f90
par_mv_3b1

# ╔═╡ 7b8e57db-5893-45da-83d2-bdd2c58b4650
R2_mv_3b1

# ╔═╡ add56d31-d978-4e14-9981-eff6a1fabb59
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_2a Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_2a Slider(1:length(rT_f)))
"""

# ╔═╡ 25743005-bb19-43de-9077-4f12b60d6136
p_fit_mv_2a = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_2a, model_mv_2a, show_i_Tinit_mv_2a, show_i_rT_mv_2a)

# ╔═╡ cea827e6-b1ae-4e50-85ae-d3e8566100ca
begin
	plotly()
	p_mv_2a_Tinit = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	Telu_Tchar_2a = filter([:rT] => x -> x == rT[show_i_rT_mv_2a], Telu_Tchar_rT_data_excluded)
	scatter!(p_mv_2a_Tinit, Telu_Tchar_2a.Tchar, Telu_Tchar_2a.Tinit, Telu_Tchar_2a.Telu, label=round(rT[show_i_rT_mv_2a];sigdigits=3))
	xTchar_2a = sort!(unique(Telu_Tchar_2a.Tchar))
	yTinit_2a = sort!(unique(Telu_Tchar_2a.Tinit))
	zdata_2a = Array{Float64}(undef, length(xTchar_2a), length(yTinit_2a))
	for i=1:length(xTchar_2a)
		for j=1:length(yTinit_2a)
			zdata_2a[i,j] = model_mv_2a([xTchar_2a[i] rT[show_i_rT_mv_2a] yTinit_2a[j]], fit_mv_2a.param)[1]
		end
	end
	plot!(p_mv_2a_Tinit, xTchar_2a, yTinit_2a, zdata_2a', st=:surface)
end

# ╔═╡ e9f02fe0-aa5d-4ed5-a913-19bccfe0ae7f
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_2b Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_2b Slider(1:length(rT_f)))
"""

# ╔═╡ 574ffed2-05ab-4ca9-8922-0dce17bd547d
p_fit_mv_2b = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_2b, model_mv_2b, show_i_Tinit_mv_2b, show_i_rT_mv_2b)

# ╔═╡ 3820b600-cfd2-4a7f-a2d3-fc7accc7d7f5
begin
	plotly()
	p_mv_2b_Tinit = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	Telu_Tchar_2b = filter([:rT] => x -> x == rT[show_i_rT_mv_2b], Telu_Tchar_rT_data_excluded)
	scatter!(p_mv_2b_Tinit, Telu_Tchar_2b.Tchar, Telu_Tchar_2b.Tinit, Telu_Tchar_2b.Telu, label=round(rT[show_i_rT_mv_2b];sigdigits=3))
	xTchar_2b = sort!(unique(Telu_Tchar_2b.Tchar))
	yTinit_2b = sort!(unique(Telu_Tchar_2b.Tinit))
	zdata_2b = Array{Float64}(undef, length(xTchar_2b), length(yTinit_2b))
	for i=1:length(xTchar_2b)
		for j=1:length(yTinit_2b)
			zdata_2b[i,j] = model_mv_2b([xTchar_2b[i] rT[show_i_rT_mv_2b] yTinit_2b[j]], fit_mv_2b.param)[1]
		end
	end
	plot!(p_mv_2b_Tinit, xTchar_2b, yTinit_2b, zdata_2b', st=:surface)
end

# ╔═╡ ddee026d-df9c-4cf0-b50d-9c8993973a45
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_2b0 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_2b0 Slider(1:length(rT_f)))
"""

# ╔═╡ 381f6f06-cbcc-475b-982e-f67b479b486d
p_fit_mv_2b0 = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_2b0, model_mv_2b0, show_i_Tinit_mv_2b0, show_i_rT_mv_2b0)

# ╔═╡ f945bc65-1b64-4157-bb0c-d96127f1d49a
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_2b1 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_2b1 Slider(1:length(rT_f)))
"""

# ╔═╡ 7f7d17ff-84ec-490a-a90f-f53bea8d8172
p_fit_mv_2b1 = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_2b1, model_mv_2b1, show_i_Tinit_mv_2b1, show_i_rT_mv_2b1)

# ╔═╡ ac9f8221-7258-4017-926c-fd8ecd4b9553
begin
	plotly()
	p_mv_2b1_Tinit = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	Telu_Tchar_2b1 = filter([:rT] => x -> x == rT[show_i_rT_mv_2b1], Telu_Tchar_rT_data_excluded)
	scatter!(p_mv_2b1_Tinit, Telu_Tchar_2b1.Tchar, Telu_Tchar_2b1.Tinit, Telu_Tchar_2b1.Telu, label=round(rT[show_i_rT_mv_2b1];sigdigits=3))
	xTchar_2b1 = sort!(unique(Telu_Tchar_2b1.Tchar))
	yTinit_2b1 = sort!(unique(Telu_Tchar_2b1.Tinit))
	zdata_2b1 = Array{Float64}(undef, length(xTchar_2b1), length(yTinit_2b1))
	residuen_2b1 = Array{Float64}(undef, length(xTchar_2b1), length(yTinit_2b1))
	for i=1:length(xTchar_2b1)
		for j=1:length(yTinit_2b1)
			zdata_2b1[i,j] = model_mv_2b1([xTchar_2b1[i] rT[show_i_rT_mv_2b1] yTinit_2b1[j]], fit_mv_2b1.param)[1]
			Telu = filter([:Tchar, :Tinit] => (x,y) -> x == xTchar_2b1[i] && y == yTinit_2b1[j], Telu_Tchar_2b1)
			residuen_2b1[i,j] = Telu.Telu[1] - zdata_2b1[i,j]
		end
	end
	plot!(p_mv_2b1_Tinit, xTchar_2b1, yTinit_2b1, zdata_2b1', st=:surface)
	p_mv_2b1_residuen = plot(xTchar_2b1, yTinit_2b1, residuen_2b1', st=:surface, xlabel="Tchar", ylabel="Tinit", zlabel="residuen(Telu)")
	p_mv_2b1_Tinit, p_mv_2b1_residuen
end

# ╔═╡ 637a3ec5-dab7-4fc0-b78f-2da6edec5547
Telu_ = filter([:Tchar, :Tinit] => (x,y) -> x == xTchar_2b1[2] && y == yTinit_2b1[5], Telu_Tchar_2b1)

# ╔═╡ d11cb6a8-6786-467b-a5f5-4ea2e163e199
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_3a Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_3a Slider(1:length(rT_f)))
"""

# ╔═╡ 4e75368a-2a81-4c43-aa48-21bc208ee66b
p_fit_mv_3a = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_3a, model_mv_3a, show_i_Tinit_mv_3a, show_i_rT_mv_3a)

# ╔═╡ 64ae3896-8b83-472c-9113-7ed879e751c0
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_3b Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_3b Slider(1:length(rT_f)))
"""

# ╔═╡ c9cdc785-5415-44de-9223-d817f6948553
p_fit_mv_3b = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_3b, model_mv_3b, show_i_Tinit_mv_3b, show_i_rT_mv_3b)

# ╔═╡ f0e1e4dd-140f-4bc9-80d2-5dcfe87f385e
md"""
``T_{init}``: $(@bind show_i_Tinit_mv_3b1 Slider(1:length(Tinit_f)))

``r_T``: $(@bind show_i_rT_mv_3b1 Slider(1:length(rT_f)))
"""

# ╔═╡ 81af5cb4-a41b-407f-8578-317aef73f41d
p_fit_mv_3b1 = plot_fit_mv(Telu_Tchar_rT_data_excluded, fit_mv_3b1, model_mv_3b1, show_i_Tinit_mv_3b1, show_i_rT_mv_3b1)

# ╔═╡ 90193caa-5e71-469c-8cd0-51ea21b034e6
begin
	plotly()
	p_mv_3b1_Tinit = plot(xlabel="Tchar", ylabel="Tinit", zlabel="Telu")
	Telu_Tchar_3b1 = filter([:rT] => x -> x == rT[show_i_rT_mv_3b1], Telu_Tchar_rT_data_excluded)
	scatter!(p_mv_3b1_Tinit, Telu_Tchar_3b1.Tchar, Telu_Tchar_3b1.Tinit, Telu_Tchar_3b1.Telu, label=round(rT[show_i_rT_mv_3b1];sigdigits=3))
	xTchar_3b1 = sort!(unique(Telu_Tchar_3b1.Tchar))
	yTinit_3b1 = sort!(unique(Telu_Tchar_3b1.Tinit))
	zdata_3b1 = Array{Float64}(undef, length(xTchar_3b1), length(yTinit_3b1))
	residuen_3b1 = Array{Float64}(undef, length(xTchar_3b1), length(yTinit_3b1))
	for i=1:length(xTchar_3b1)
		for j=1:length(yTinit_3b1)
			zdata_3b1[i,j] = model_mv_3b1([xTchar_3b1[i] rT[show_i_rT_mv_3b1] yTinit_3b1[j]], fit_mv_3b1.param)[1]
			Telu = filter([:Tchar, :Tinit] => (x,y) -> x == xTchar_3b1[i] && y == yTinit_3b1[j], Telu_Tchar_3b1)
			residuen_3b1[i,j] = Telu.Telu[1] - zdata_3b1[i,j]
		end
	end
	plot!(p_mv_3b1_Tinit, xTchar_3b1, yTinit_3b1, zdata_3b1', st=:surface)
	p_mv_3b1_residuen = plot(xTchar_3b1, yTinit_3b1, residuen_3b1', st=:surface, xlabel="Tchar", ylabel="Tinit", zlabel="residuen(Telu)")
	p_mv_3b1_Tinit, p_mv_3b1_residuen
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
# ╠═39655593-9b8c-444e-8f7d-48cd92d7b9d8
# ╠═686250d4-235d-4699-88cc-5aaf1273b319
# ╠═4e46e5eb-74ca-400b-bca2-7329af1879ac
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
# ╠═eb37f6b4-23b8-4163-8ec4-67324aba214a
# ╠═b373c92c-d465-4a0f-a925-ae50360c277e
# ╠═0fd2139b-3e50-48ad-9d76-6166b91d4a5b
# ╠═b5d6fb34-00ba-42c1-9d57-518bead48f4b
# ╠═4eac8fc5-3d30-44b6-828e-e40f17bd1d7e
# ╟─668bfd06-7308-4cc8-83e2-34d51f9730fd
# ╠═f54b0681-3a3e-45cf-b274-4962989821dc
# ╠═789f8085-76c9-4fba-95d2-b6b6b603f868
# ╟─65c91137-034b-43c5-a7d1-c25cc694b1d4
# ╟─7c07e06f-5c31-4453-8d0a-d89842ebe7d2
# ╟─e97005f1-f165-4dfb-8be0-4d078eb27930
# ╟─c6143a16-c8d9-48ec-8588-7f430851576b
# ╟─bd87d080-07fc-4c99-9e6c-1df8229cb46c
# ╟─c3e32361-b417-4305-acf3-5862813d2d0c
# ╟─5bc9a02e-62b1-48de-8287-33be309402fe
# ╟─e950d80d-8b1a-427c-a4f0-361d62f9a450
# ╟─600773b7-8b3b-4e30-91ca-3ce0407dae68
# ╠═6bef1367-bc12-40cd-8c6d-2d5934cdb9d6
# ╠═aa34de2c-dbf1-4551-aac7-94f361fa9b6a
# ╠═8a730b2f-a51f-40ac-a3a9-054f0c903031
# ╠═43d45ad1-0289-4992-a089-22c6b299931d
# ╠═f6677948-c480-4557-ac40-e38a76470daf
# ╠═58e58468-ee89-4f00-bed9-bbfd2f494be1
# ╠═b552ac0c-7f39-4fd6-bd37-13aae5a81d99
# ╠═0e504ae7-4165-4261-a02d-85141e2f5b5f
# ╟─10266ff0-7e64-41c9-b420-d6c15257ceef
# ╠═16991c8c-b7c3-495b-93ce-70b949017b3c
# ╠═f917c1c9-8bb8-4e35-bb92-24df9e2945a1
# ╠═03121d84-5b1b-41b7-befa-c629b08e6f52
# ╠═101ffb81-ada0-48fa-9df8-272800af2565
# ╠═be2237d4-be73-4368-9f38-88280732f225
# ╠═c5d95d97-8903-4dbf-ad12-359d2a235e6c
# ╠═f518d526-34ca-4c15-8255-4ee79c8581dd
# ╠═3844df54-1622-4313-92ee-36c6599f80d0
# ╠═4071fa66-095c-49e7-adfe-44c1fce9c59a
# ╠═a6b9948f-7364-4294-a744-fd13550a4796
# ╠═69c0a9bc-4606-4e90-aed1-6ea8dbc6f039
# ╠═ca7c60d2-04ad-4980-b26f-c02c70ef84a8
# ╠═fbde7fa1-c584-4732-ba60-05837ee62d8f
# ╠═23ccbdb4-2c7c-499f-af74-a093304d60cd
# ╠═69113eb5-cca4-40e7-8a5e-96e196a99fc0
# ╠═388a86d4-766f-4937-9a7a-4bd583b1c1b5
# ╠═85696436-f843-4a2a-a582-a3d27c29d951
# ╠═59e7835e-e9f2-4888-b862-28ba1a019ab6
# ╠═aec7dae9-c73a-4cce-a440-2724799085e1
# ╠═ed71d956-bc13-407c-8aab-02d700b711b7
# ╠═1a8865f8-db21-423f-8591-033f41197d2e
# ╠═32782159-a1ee-4570-bdbc-31a4506ad902
# ╠═997409e6-1292-4bc0-94f7-b1439d9beeaa
# ╠═e8e0b85f-286e-4665-a992-3f26f8eb016e
# ╠═34b5b40f-3be2-4b9f-ad4c-f4d99dac1299
# ╠═78c48193-7021-4803-bf4b-bf4e03406c47
# ╠═494c2b8c-04b3-46ff-9c28-cc9b79f9a7d4
# ╠═ac33cbec-938b-472c-863f-97b8c3d6093d
# ╠═42411a67-69f4-4088-8b01-2815f223e4f9
# ╠═76873708-0cd3-4048-a33b-29297fb58813
# ╠═8c9a6c3a-78ac-4a60-af68-7306dbb79b7b
# ╠═137174bf-fde5-436a-ab76-92168bca291f
# ╠═5fe03d15-d5cc-4c3b-af2a-4ed35b40bfad
# ╠═d91ec9d0-d9db-47b4-973f-c9565825dc71
# ╠═ef742f31-4a9c-4d49-ad12-e3a81fa881cb
# ╠═bee6bac9-c137-4979-83e2-23462660d9bf
# ╠═ec7bb02e-b4cd-47aa-a694-07c14c0b5e86
# ╠═bd88682c-a94d-431b-a5dc-dbc7638a3e48
# ╠═47bf2c03-bfa0-4ecc-ace3-a64959c558d1
# ╠═0c741a3f-5f76-4c86-be39-e5a02c2fe7e1
# ╠═1d3474df-ede5-427b-8bf8-c50ecdfaedb5
# ╠═073912f2-5dcd-4c70-810d-9196709ee2d5
# ╠═5c8678b9-d0f3-4db3-8f0d-42ff43cd511c
# ╠═bf37a7dc-2870-489c-9ff9-3b2423f4870c
# ╠═d5d09e6d-667c-4d5d-aaf3-0f4841c49ca1
# ╠═e0a8ba27-65c0-4425-a1f2-1f30fe768c30
# ╠═241a52ee-c094-4198-b964-67847ec313ec
# ╠═13e14de4-09d0-4013-8176-ac96c678a0a6
# ╠═e809c0f8-b939-413d-b6b8-1c76524c7d89
# ╠═7c3d5c74-3710-4671-891a-f2f75d23ee7e
# ╠═e0edb329-ba6b-4338-9386-6a0d34ff7c61
# ╠═4b6b6c34-d103-4483-a483-00162a5f2583
# ╠═f7f053fa-f0d9-4a40-bf31-3edc12f1f9df
# ╠═59ae1c1b-a047-4bcb-b1a1-9f6e59d36695
# ╠═e0196c9b-ce1f-4417-b8e7-f1c9a942efb2
# ╠═33efc8a9-b2ad-4cb9-8d90-245e3463c393
# ╠═5eab4c51-789d-4f78-825b-9514f948cde5
# ╠═accf6e87-0b2c-4251-bff0-941d84053733
# ╠═86e15050-d13f-4c6f-99d0-5ad2fc011937
# ╠═351892bd-40ff-4162-b6be-222188f15eab
# ╠═4ebe8f12-c38e-48cc-9efb-562763684dc0
# ╠═40805b6a-c69b-48b8-aefd-49a731c06de7
# ╠═3676b513-c1b6-46e9-b069-b14be86a78eb
# ╠═3793b682-2b9f-410a-86af-53791b42979d
# ╠═5a501229-fe26-42c5-85ad-366baddfbaee
# ╠═26706952-0716-42c7-804a-738801be929e
# ╠═5630427e-22b4-4437-963c-a8c49ee76d57
# ╠═3c61bc26-c4cd-4cad-b405-272b0a962e49
# ╠═69a69280-f3b9-49ea-98fd-b439b495428a
# ╠═add56d31-d978-4e14-9981-eff6a1fabb59
# ╠═25743005-bb19-43de-9077-4f12b60d6136
# ╠═cea827e6-b1ae-4e50-85ae-d3e8566100ca
# ╠═78935736-f4b1-41c5-bd43-328e41eb53e8
# ╠═fb0c53bd-7e57-417f-a2eb-584896eaf6b3
# ╠═01cd2d11-2e63-4bba-b7e5-84404cb55f0f
# ╠═96222088-76e8-4c36-9a04-fce4ee4a6c94
# ╠═ddcfb67c-3bcf-4320-8603-5db5e42cccf1
# ╠═d1a18fb3-2114-401c-91df-d2ce348aa1e1
# ╠═e9f02fe0-aa5d-4ed5-a913-19bccfe0ae7f
# ╠═574ffed2-05ab-4ca9-8922-0dce17bd547d
# ╠═3820b600-cfd2-4a7f-a2d3-fc7accc7d7f5
# ╠═d111ec3e-62d3-4737-bf56-c62eb3d034ac
# ╠═d9668965-1e6f-4088-a686-2657922aa76e
# ╠═deba9143-50fa-4bf3-b7c7-a7a93c404f4c
# ╠═5c239244-fdae-4b17-b6af-09054c3c26d1
# ╠═19d9940d-7b58-4f22-9e28-55a51cc3dc56
# ╠═b0bc5c1f-c8d7-4b72-b302-226c83b1001f
# ╠═8b89812b-de28-4a4a-8d66-4979bfb41e27
# ╠═bbf1b7de-8d8f-49ba-8c88-3cd85d403560
# ╠═40d8b557-65cf-4397-b429-2d020ba2c8c8
# ╠═62663709-8eed-4777-a121-61f5ca5c9973
# ╠═ed9d99ee-01f3-4f59-9d26-5fe0841e0a40
# ╠═525a9093-ce35-4f43-95fb-c2729cdb02b4
# ╠═cbf3837e-dbda-4c20-83d5-74b936877a14
# ╠═3a4c953f-a881-4d0d-b12a-015c39faf0bb
# ╠═ddee026d-df9c-4cf0-b50d-9c8993973a45
# ╠═381f6f06-cbcc-475b-982e-f67b479b486d
# ╠═983dfe18-83fb-4eff-9da6-e6e4ca79d506
# ╠═445c8480-19dd-4a5a-a022-97a758341a48
# ╠═9eb0565c-b76f-4cdc-952d-c6139e4e6f63
# ╠═146a2948-5a80-4e49-a3b0-615e7bd52769
# ╠═c123731b-f488-485c-bf5f-d46ba56a9f9e
# ╠═1c23b58a-4f4f-426c-ba4e-21e04858fdfb
# ╠═f945bc65-1b64-4157-bb0c-d96127f1d49a
# ╠═7f7d17ff-84ec-490a-a90f-f53bea8d8172
# ╠═ac9f8221-7258-4017-926c-fd8ecd4b9553
# ╠═637a3ec5-dab7-4fc0-b78f-2da6edec5547
# ╠═3286b6ee-2deb-4c2b-b1aa-3284b9697139
# ╠═e36fbd07-42d2-44a3-8bcd-a80e19df18b1
# ╠═840ec2d9-5389-4ba2-9a7b-9346763f93dc
# ╠═5a16fa66-5ddf-4f18-b398-5b47c881e3b3
# ╠═5f89244c-45f8-4115-bee7-2cc9ee8ef903
# ╠═caa4e339-87c8-4ffa-94ae-1617f562ce0b
# ╠═d11cb6a8-6786-467b-a5f5-4ea2e163e199
# ╠═4e75368a-2a81-4c43-aa48-21bc208ee66b
# ╠═dd1670fe-9968-42a5-8f7a-3fb8c72c197e
# ╠═f60a676a-3458-4681-9634-dabaede63383
# ╠═2f6bfa25-acd7-41cf-8429-29c9df694a8d
# ╠═13f1725c-9355-410c-aa25-5f4de21e0bf4
# ╠═0013f6f5-f547-49a0-917b-89ee473f9439
# ╠═f6b370b5-4944-44e3-a9a7-31e10798399e
# ╠═64ae3896-8b83-472c-9113-7ed879e751c0
# ╠═c9cdc785-5415-44de-9223-d817f6948553
# ╠═084f8e7f-ca98-4727-9346-2315b5f51fc2
# ╠═f00ff9f1-bb8a-4e66-a842-479a9b68e4e2
# ╠═e66578f9-c40e-448f-aba7-12ff5e7738d8
# ╠═7239b6aa-86fc-42c2-ad54-2fcee51f3d81
# ╠═c1c961e9-41ea-43aa-bde5-0d94ba277f90
# ╠═7b8e57db-5893-45da-83d2-bdd2c58b4650
# ╠═f0e1e4dd-140f-4bc9-80d2-5dcfe87f385e
# ╠═81af5cb4-a41b-407f-8578-317aef73f41d
# ╠═90193caa-5e71-469c-8cd0-51ea21b034e6
# ╠═7871ce05-73de-4b88-9b00-9a8a735ef3da
# ╠═1381f354-84ab-4496-9200-d99b10d98ef3
# ╠═8e1c640d-75bd-46d1-a5cf-ec63f493cf10
# ╠═05a1c3cb-c778-4d9d-be80-88c1136320d9
# ╠═c31eee97-b94d-463a-bb91-41f18cd83742
# ╠═6bbe03d7-cb42-4ff5-92a8-6441515ba8cc
# ╠═4a14a4ba-162d-46bf-98c0-46107fff0b45
# ╠═98eb4f63-c4fc-465f-928c-f7f805d49839
# ╠═b86e3b4d-b216-49b9-9e5d-73d57f24c72e
# ╠═ffb0843d-dcaf-4207-bd76-e1e21ebdc2cd
# ╠═2e224d8b-6fb1-4b0b-aacb-28294f8e7b8d
# ╠═fd778a85-b7f5-4d87-b671-4ac3efee4b79
# ╠═1739ae21-a60d-413e-832d-cb12fad0ca34
# ╠═7bc0108e-6842-4850-92a5-0cae5d7d4706
# ╟─2f750d6e-461f-4e6a-bc92-2939dbb55421
# ╠═5e26e413-612c-4728-9c08-15f45aa10e1f
# ╠═2c28f8ad-ab76-4131-8856-55a43b31fde6
# ╠═804fd07b-b4f0-4119-a708-60d4973e84eb
# ╠═59e2f85f-acdd-4ca5-a91c-78b9bb42364e
# ╠═030af251-fe5b-46c8-b6d1-1b2efb7cd6c3
# ╠═4fe853a2-100b-4c25-8559-1ea008c0cad2
# ╠═3f7dad44-11bc-4535-b7cb-db5938f5d2d4
# ╠═98df485c-de5b-4ea9-ba12-7d0677a4e7dc
# ╠═2a0e8394-076a-4d0f-8458-b2f838c9bd77
# ╠═7d9ae1fa-047c-11ed-0af1-c1af45ffc663
