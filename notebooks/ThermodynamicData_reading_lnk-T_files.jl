### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ ae71c5ce-9978-11ec-3f73-1197d152796a
using CSV, DataFrames, LambertW, ChemicalIdentifiers, LsqFit, Statistics, Measurements

# ╔═╡ 2028b40d-13a5-4c22-8335-8691bbe5a0ee
using Plots

# ╔═╡ d83d45da-f852-424e-b3aa-2a5881de28fe
using PlutoUI

# ╔═╡ aaa3793e-b9c8-4742-8005-33c336787bf7
include("/Users/janleppert/Documents/GitHub/ThermodynamicData/src/ThermodynamicData.jl")

# ╔═╡ 0e80aed9-f0ac-4077-b1f4-b298349d31c8
db_path = "/Users/janleppert/Documents/GitHub/ThermodynamicData/Databases"

# ╔═╡ 413fad2a-e73f-4933-9eed-b1aae8d376ce
"""
	load_lnkT_data(db_path)

Load the data files (.csv format) with `lnk-T` data from the folder `db_path` including all subfolders.
"""
function load_lnkT_data(db_path)
	all_csv = ThermodynamicData.collect_csv_paths(db_path)
	keyword1 = "log10k-T"
	keyword2 = "lnk-T"
	log10kT_csv = filter(contains(keyword1), all_csv)
	lnkT_csv = filter(contains(keyword2), all_csv)
	meta_data = ThermodynamicData.extract_meta_data([log10kT_csv; lnkT_csv])
	ThermodynamicData.load_csv_data!(meta_data)
	return meta_data
end

# ╔═╡ e2566bfc-1912-45bb-90c4-b05ee59c16af
all_csv = ThermodynamicData.collect_csv_paths(db_path)

# ╔═╡ ffa74790-bb57-4dc9-9bfc-8309d6d58f49
begin # filter for certain key words in the filename
	keyword1 = "lnk-T"
	keyword2 = "log10k-T"
	lnkT_csv = filter(contains(keyword1), all_csv)
	log10kT_csv = filter(contains(keyword2), all_csv)
end

# ╔═╡ 34c6c2bd-63a3-4af1-b6ac-1dccdf3ddb4c
[log10kT_csv; lnkT_csv]

# ╔═╡ 15047914-7463-4d22-934a-00f7c81ee44a
meta_data = ThermodynamicData.extract_meta_data([log10kT_csv; lnkT_csv])

# ╔═╡ baf00ccd-0a30-4cb4-9704-d5d9b23cfe69
ThermodynamicData.load_csv_data!(meta_data)

# ╔═╡ 422b832b-0b52-4c35-995b-a904069d42c5
contains(meta_data.filename[2], "log10k-T")

# ╔═╡ 29232c11-9e21-42aa-aca3-500563bbc9a9
log(2)/log10(2)

# ╔═╡ 9ea6d953-8556-4ddc-806b-02a74a3b94d8
Ts = names(meta_data.data[11])[2:end]

# ╔═╡ fbd168b6-8d91-4e60-b2ac-342e787ccb14
TT = Array{Float64}(undef, length(Ts))

# ╔═╡ a356e7c7-c1da-4907-b8b4-e7cf7c6f5c57
occursin(Ts[12], "_")

# ╔═╡ 035a63b2-1d9b-417b-ae3a-f27604a50898
occursin("_", Ts[12])

# ╔═╡ b946b6a7-4453-44a6-9419-1ecd235d2956
parse(Float64,split(Ts[13], "_")[1])

# ╔═╡ 5655ccd8-09db-40b9-9b6d-f88b5d328d43
"""
	T_column_names_to_Float(data::DataFrame)

Translates the column names of the dataframe `data` containing the values of isothermal temperature to Float64. For the case of identical temperatures, a number is appended, separated by an underscore `_`. If this is the case the string is splited at the `_` and only the first part is used.
"""
function T_column_names_to_Float(data::DataFrame)
	Ts = names(data)[2:end]
	T = Array{Float64}(undef, length(Ts))
	for i=1:length(Ts)
		if occursin("_", Ts[i])
			T[i] = parse(Float64,split(Ts[i], "_")[1])
		else
			T[i] = parse(Float64, Ts[i])
		end
	end
	return T
end

# ╔═╡ 77b226e1-4c36-47e7-b181-511667de9968
TT

# ╔═╡ f754b0d2-2f63-40fa-92c8-a692f52975c3
"""
	ABC(x, p)

The ABC-model for the relationship lnk(T). 

# Arguments:
- `x`: temperature `T` in K
- `p[1] = A - lnβ`
- `p[2] = B` in K
- `p[3] = C`
"""
@. ABC(x, p) = p[1] + p[2]/x + p[3]*log(x)

# ╔═╡ f99640cf-58e9-4bb9-997b-e76578edebfa
"""
	Kcentric(x, p)

The K-centric model for the relationship lnk(T).

# Arguments
- `x`: temperature `T` in K
- `p[1] = Tchar + Tst` in K
- `p[2] = θchar` in °C
- `p[3] = ΔCp/R`
"""
@. Kcentric(x, p) = (p[3] + p[1]/p[2])*(p[1]/x - 1) + p[3]*log(x/p[1])

# ╔═╡ 1d466a06-3865-477e-bc10-d811eb6a2732


# ╔═╡ 4a2d1ecb-02b6-44a5-953e-eedc6cef512f
T = 0.0:1.0:400.0

# ╔═╡ 996a3ca9-ede4-4ecf-aa24-38f3fb104cb2


# ╔═╡ 3624acaf-c9f8-4548-88da-14a22e3e054f
flagged_index = Int[]

# ╔═╡ 69e663e4-3814-419b-9225-1574eedb9b54
ABC0 = [-100.0, 10000.0, 10.0]

# ╔═╡ 14d85174-45f2-459f-82a3-498c9aa5b1b8
threshold = 0.001

# ╔═╡ 578bade5-7af7-4dae-815e-41ab1c9311d5
a = [3,4,6,7,8,9]

# ╔═╡ 0ac891c1-4ef5-4267-bfa3-13e79c316534
i_a = 8

# ╔═╡ be15dac4-309f-459a-b86b-c80fb12421bc
a[findall(a.!=i_a)]

# ╔═╡ e803eef6-5b7e-4f8b-a89c-0ea99117cf35
# evaluation of Fit

# ╔═╡ ba2dfc06-37df-44df-b8f4-50a506af4151
"""
	coeff_of_determination(fit, y)

Calculate the coefficient of determination R² for LsqFit.jl result `fit` and measured data `y`.
"""
function coeff_of_determination(fit, y)
	SSres = sum(fit.resid.^2)
	SStot = sum((y .- mean(y)).^2)
	R² = 1 - SSres/SStot
	return R²
end

# ╔═╡ bd6acff3-93b1-401f-9ec3-2d1d84485a4d
"""
	chi_square(fit)

Calculate chi-square χ² and the adjusted chi-square χ̄² for LsqFit.jl result `fit`.
"""
function chi_square(fit)
	χ² = sum(fit.resid.^2)
	χ̄² = sum(fit.resid.^2)/dof(fit)
	return χ², χ̄²
end

# ╔═╡ c25f377d-cd5c-4626-a09b-f0043e478b40
log(1000)

# ╔═╡ 09ed6ad1-f161-4288-9152-1be342267b5e
meta_data.beta0

# ╔═╡ 6743fce9-f000-4757-aca5-f347e4c38ef7
names(meta_data)

# ╔═╡ c61f2185-5c78-4cd9-81b5-2c9e472d5c37
R = 8.31446261815324

# ╔═╡ e270cf70-f13f-4fc5-91d4-b5181621f9da
T0 = 90.0

# ╔═╡ 9ad6c6da-27e4-4f22-b693-1689824aad0c
md""" 
# Next:
- [x] add fit-function version with residual-threshold, exclusion of data and re-fit
- [x] calculate R² (in fit dataframe), Χ² (only function defined)
- [x] extract the parameters A, B, C, Tchar, θchar and ΔCp (including errors, Measurements.jl)
- [x] calculate ΔHref and ΔSref
- [] export `AllParam` files

Separat new notebook
- [] read all `AllParam` files -> make one big dataframe
								-> flag low quality data
 								-> identify duplicates
 								-> ChemicalIdentifiers.jl
 								-> export one big data base resp. databases for every source (in old and new format)
"""

# ╔═╡ 1f252318-8645-4be0-a900-a9585f1c73ff
Tst = 273.15

# ╔═╡ e011d3c1-7070-4e4c-b804-d63424aa6f20
"""
	fit_models(data::Array{DataFrame}, β0::Array{Float64})

Fit the ABC-model and the K-centric-model at the lnk(T) data of the `data`-array of dataframes, using LsqFit.jl 
"""
function fit_models(data::Array{DataFrame}, β0::Array{Float64})

	fit = Array{DataFrame}(undef, length(data))
	for i=1:length(data)
		fitABC = Array{Any}(undef, length(data[i][!,1]))
		fitKcentric = Array{Any}(undef, length(data[i][!,1]))
		T = Array{Array{Float64}}(undef, length(data[i][!,1]))
		lnk = Array{Any}(undef, length(data[i][!,1]))
		name = Array{String}(undef, length(data[i][!,1]))
		ok_i = Array{Array{Int,1}}(undef, length(data[i][!,1]))
		R²_ABC = Array{Float64}(undef, length(data[i][!,1]))
		R²_Kcentric = Array{Float64}(undef, length(data[i][!,1]))
		for j=1:length(data[i][!,1])
			T[j] = T_column_names_to_Float(data[i])
			lnk[j] = collect(Union{Missing,Float64}, data[i][j,2:end])
			ii = findall(isa.(lnk[j], Float64))
			name[j] = data[i][!, 1][j]
			ABC0 = [-100.0, 10000.0, 10.0] # for p[1] = -100.0 - log(beta0) ?
			lb_ABC = [-Inf, 0.0, 0.0]
			ub_ABC = [-log(β0[i]), Inf, Inf]
			Kcentric0 = [200.0+Tst, 30.0, 10.0]
			lb_Kcentric = [0.0, 0.0, 0.0]
			ub_Kcentric = [Inf, Inf, 500.0]
			fitABC[j] = curve_fit(ABC, T[j][ii].+Tst, lnk[j][ii], ABC0, lower=lb_ABC, upper=ub_ABC)
			fitKcentric[j] = curve_fit(Kcentric, T[j][ii].+Tst, lnk[j][ii], Kcentric0, lower=lb_Kcentric, upper=ub_Kcentric)
			ok_i[j] = ii
			R²_ABC[j] = coeff_of_determination(fitABC[j], lnk[j][ii])
			R²_Kcentric[j] = coeff_of_determination(fitKcentric[j], lnk[j][ii])
		end
		fit[i] = DataFrame(Name=name, T=T, lnk=lnk, fitABC=fitABC, fitKcentric=fitKcentric, ok_i=ok_i, R²_ABC=R²_ABC, R²_Kcentric=R²_Kcentric)
	end
	return fit
end

# ╔═╡ 9025d205-965a-4843-a0b9-eeaa8d98d997
fit = fit_models(meta_data.data, meta_data.beta0)

# ╔═╡ 049e3387-eead-4806-acdc-77ff0e292d53
md"""
Select data set: $(@bind select_dataset Slider(1:length(fit); show_value=true))
"""

# ╔═╡ a10bec56-173e-460c-8eb5-040f3da1244c
md"""
Select substance: $(@bind select_substance Slider(1:length(fit[select_dataset].Name); show_value=true))
"""

# ╔═╡ 13949a69-a1fc-4b81-9eb6-66269e1bf853
fit[4]

# ╔═╡ 4c2eed13-f7fc-4e0e-ad64-e8bf2894e98e
fit[1].fitABC[14]

# ╔═╡ d7f8407c-5ad2-4df4-a222-53f2c478fd94
scatter(fit[1].fitABC[14].resid)

# ╔═╡ 9c66eca0-9ef7-4e2e-8c37-0d2c522df5b0
fit[1].T[14]

# ╔═╡ e7113c1b-e465-4983-a01c-52f00ed6b16f
maximum(abs.(fit[3].fitABC[2].resid))

# ╔═╡ e7dcf5dd-ab2d-4e45-ba9c-50d55896beb2
T_ = fit[1].T[10]

# ╔═╡ f404511b-41d1-4cac-a260-9409184becdd
lnk_ = fit[1].lnk[10]

# ╔═╡ 8c2ec059-872e-4e58-851c-e3628b477b33
ii_ = findall(isa.(fit[1].lnk[10], Float64))

# ╔═╡ caddeaba-0ad6-45a0-8af0-2fa61d09663c
ok_index = findall(ismissing.(lnk_[ii_]).==false)

# ╔═╡ 30f10519-ff73-4d46-8837-070d9c46cdbf
length(ok_index) <= 3

# ╔═╡ 5335f894-8b7a-41e5-a0e5-292bdde60881
fit[2].fitABC[5]

# ╔═╡ 0a5b2bc6-9dd0-4483-af64-50e889a7a1c4
sigma = stderror(fit[2].fitABC[5])

# ╔═╡ f7e45701-3ce6-49fa-942a-90102c82c784
margin_of_error = margin_error(fit[2].fitABC[5], 0.05)

# ╔═╡ bfbaf273-8993-4588-8822-220efc242bb9
margin_of_error./sigma

# ╔═╡ 9b1271c5-7686-4f23-9f31-4a7028afe1b1
confidence_inter = confidence_interval(fit[2].fitABC[5], 0.05)

# ╔═╡ c3e1408f-9cf7-4418-8b66-fe3e502a8278
covar = estimate_covar(fit[2].fitABC[5])

# ╔═╡ 1d273a65-577d-42ad-9534-7acf6ca1e176
SSres = sum(fit[2].fitABC[5].resid.^2) # Χ²

# ╔═╡ 70b59e3b-373f-455a-b58d-3083cacf3c3b
dof(fit[2].fitABC[5]) # degree of freedom

# ╔═╡ 17367fcc-8737-42a5-9546-4d7447420248
estimate_covar(fit[2].fitABC[5]).^2

# ╔═╡ cece68a2-86bf-49f1-bfb0-5dcc8a41eb5b
chi_square(fit[2].fitABC[5])

# ╔═╡ d9494de6-248e-4731-93ba-5cce35fcf0de
1.0./fit[2].fitABC[5].resid, 1.0./fit[2].fitABC[5].resid.^2, fit[2].fitABC[5].resid

# ╔═╡ 94bbe840-f321-4699-9a57-21a4291cf56e
fit[2].fitABC[5]

# ╔═╡ cf44ed56-ef17-4bb6-a131-82bfd405a3fe
fit[2].fitKcentric[5]

# ╔═╡ 14bd0e24-c16d-45bc-aaa3-d0775220374b
w1 = 1.0./fit[2].fitKcentric[5].resid.^2

# ╔═╡ 5ac44f86-f9c3-4a59-9764-11c3aa0a9e66
sum(w1)

# ╔═╡ 8ff63c84-7996-47e0-87ef-d00ea1eafe93
w1./sum(w1)

# ╔═╡ 252713b0-bebc-4a56-a154-44902fcee8bf
fit[3].fitABC[5].param

# ╔═╡ 88376f08-7ee9-43e2-a0d8-0e2dc4e63904
stderror(fit[3].fitABC[5])

# ╔═╡ e43a548b-a863-4ff1-ba76-7c4392175f5b
length(fit)

# ╔═╡ 97e9dd89-bf47-4c79-ab41-84e54bb0ec2c
aa = measurement(fit[3].fitABC[5].param[1], stderror(fit[3].fitABC[5])[1])

# ╔═╡ 646b3e92-6725-4363-9c76-e94c0d84af6c
aa + log(meta_data.beta0[3])

# ╔═╡ 4a1b443e-1965-432f-abb7-886edce3a7f0
typeof(aa)

# ╔═╡ 32f2b9b2-0982-41b9-b297-845f1829d301
fit[3].fitABC[5].param[1] ± stderror(fit[3].fitABC[5])[1]

# ╔═╡ 9efc51b7-10ab-4b15-bdab-be7bb3f55872
b = fit[1].fitKcentric[1].param

# ╔═╡ ea5a5cd1-e53a-47f8-ad8d-b46aa3d61a81
b[3] * R

# ╔═╡ 1c38490d-c6b6-443f-b856-2b33e8417277
"""
	fit_models!(meta_data::DataFrame)

Fit the ABC-model and the K-centric-model at the lnk(T) data of the `data`-array of dataframes, using LsqFit.jl and add the result in a new column (`fit`) of `meta_data 
"""
function fit_models!(meta_data::DataFrame)
	fit = fit_models(meta_data.data, meta_data.beta0)
	meta_data[!, "fit"] = fit
	return meta_data
end

# ╔═╡ 126dceaf-7e01-47a4-81ce-f83f799ac422
"""
	plot_lnk_fit(fit, i, j)

Plot the `lnk`-values over `T` of the selected dataset `i` of the selected substance `j`. Also the fit solution for the ABC-model and the K-centric-model are shown over the extended temperature range from 0°C to 400°C.
"""
function plot_lnk_fit(fit, i, j)
	T = 0.0:1.0:400.0
	plnk = scatter(fit[i].T[j], fit[i].lnk[j], label="data", title=fit[i].Name[j], xlabel="T in °C", ylabel="lnk")
	if "i_ABC" in names(fit[i])
		ex_i_ABC = fit[i].ex_i_ABC[j]
		scatter!(plnk, fit[i].T[j][ex_i_ABC], fit[i].lnk[j][ex_i_ABC], label="excluded ABC", m=:diamond, mcolor=:grey, msize=2)
		ex_i_Kcentric = fit[i].ex_i_Kcentric[j]
		scatter!(plnk, fit[i].T[j][ex_i_Kcentric], fit[i].lnk[j][ex_i_Kcentric], label="excluded Kcentric", m=:cross, mcolor=:orange, msize=2)
	end
	plot!(plnk, T, ABC(T.+Tst, fit[i].fitABC[j].param), label=string("ABC, R²=", round(fit[i].R²_ABC[j]; digits=4)))
	plot!(plnk, T, Kcentric(T.+Tst, fit[i].fitKcentric[j].param), label=string("Kcentric, R²=", round(fit[i].R²_Kcentric[j]; digits=4)))
	return plnk
end

# ╔═╡ edda5af7-5ca1-4456-b5eb-2a2ff3ad89d4
plot_lnk_fit(fit, select_dataset, select_substance)

# ╔═╡ 52fe47c1-d030-4936-b627-8b7886b6211b
"""
	plot_res_lnk_fit(fit, i, j)

Plot the residuum of `lnk`-values over `T` of the selected dataset `i` of the selected substance `j` to the fitted models.
"""
function plot_res_lnk_fit(fit, i, j)
	if "i_ABC" in names(fit[i])
		i_ABC = fit[i].i_ABC[j]
		i_Kcentric = fit[i].i_Kcentric[j]
	else
		i_ABC = findall(isa.(fit[i].lnk[j], Float64))
		i_Kcentric = findall(isa.(fit[i].lnk[j], Float64))
	end
	preslnk = scatter(fit[i].T[j][i_ABC], fit[i].fitABC[j].resid, 
							label=string("ABC, R²=", round(fit[i].R²_ABC[j]; digits=4)), title=fit[i].Name[j], xlabel="T in °C", ylabel="Resid(lnk)")
	scatter!(preslnk, fit[i].T[j][i_Kcentric], fit[i].fitKcentric[j].resid, 
							label=string("Kcentric, R²=", round(fit[i].R²_Kcentric[j]; digits=4)))
	if "ex_i_ABC" in names(fit[i])
		ex_i_ABC = fit[i].ex_i_ABC[j]
		ex_i_Kcentric = fit[i].ex_i_Kcentric[j]
		scatter!(preslnk, fit[i].T[j][ex_i_ABC], ABC(fit[i].T[j][ex_i_ABC].+Tst, fit[i].fitABC[j].param), 
							label="excluded ABC", m=:diamond, mcolor=:grey, msize=2)
		scatter!(preslnk, fit[i].T[j][ex_i_Kcentric], Kcentric(fit[i].T[j][ex_i_Kcentric].+Tst, fit[i].fitKcentric[j].param), 
							label="excluded Kcentric", m=:cross, mcolor=:orange, msize=2)
	end
	return preslnk
end

# ╔═╡ 125495fc-1609-4b80-a1f3-ade4104f7da2
plot_res_lnk_fit(fit, select_dataset, select_substance)

# ╔═╡ 522d5f26-9ad5-42ce-83f2-fdabfd5d8313
"""
	fit_ABC(T, lnk, β0, res_threshold)

Fit the ABC-model for `lnk` over `T` data. Data points with a residuum above `res_threshold` are excluded and the fit is recalculated, until all residua are below the threshold or only three data points are left.
"""
function fit_ABC(T, lnk, β0, res_threshold)

	flagged_index = Int[]
	ok_index = findall(ismissing.(lnk).==false)
	
	ABC0 = [-100.0, 10000.0, 10.0] # for p[1] = -100.0 - log(beta0) ?
	lb_ABC = [-Inf, 0.0, 0.0]
	ub_ABC = [-log(β0), Inf, Inf]
	fit = curve_fit(ABC, T[ok_index].+Tst, lnk[ok_index], ABC0, lower=lb_ABC, upper=ub_ABC)
	
	# check the residuen
	while maximum(abs.(fit.resid)) > res_threshold
		if length(ok_index) <= 3
			break
		end
		imax = findfirst(abs.(fit.resid).==maximum(abs.(fit.resid))) # exclude this data point
		push!(flagged_index, ok_index[imax])
		ok_index = ok_index[findall(ok_index.!=ok_index[imax])]
		fit = curve_fit(ABC, T[ok_index].+Tst, lnk[ok_index], ABC0, lower=lb_ABC, upper=ub_ABC)
	end
	return fit, flagged_index, ok_index
end

# ╔═╡ 04a0772d-d6e7-4661-88be-45f8e06d2aff
fit_ABC(T_, lnk_, 250.0, 0.001)

# ╔═╡ 761f3c9f-7561-499b-916d-767bc3f2aafb
"""
	fit_Kcentric(T, lnk, res_threshold)

Fit the Kcentric-model for `lnk` over `T` data. Data points with a residuum above `res_threshold` are excluded and the fit is recalculated, until all residua are below the threshold or only three data points are left.
"""
function fit_Kcentric(T, lnk, res_threshold)

	flagged_index = Int[]
	ok_index = findall(ismissing.(lnk).==false)
	
	Kcentric0 = [200.0+Tst, 30.0, 10.0]
	lb_Kcentric = [0.0, 0.0, 0.0]
	ub_Kcentric = [Inf, Inf, 500.0]
	fit = curve_fit(Kcentric, T[ok_index].+Tst, lnk[ok_index], Kcentric0, lower=lb_Kcentric, upper=ub_Kcentric)
	
	# check the residuen
	while maximum(abs.(fit.resid)) > res_threshold
		if length(ok_index) <= 3
			break
		end
		imax = findfirst(abs.(fit.resid).==maximum(abs.(fit.resid))) # exclude this data point
		push!(flagged_index, ok_index[imax])
		ok_index = ok_index[findall(ok_index.!=ok_index[imax])]
		fit = curve_fit(Kcentric, T[ok_index].+Tst, lnk[ok_index], Kcentric0, lower=lb_Kcentric, upper=ub_Kcentric)
	end
	return fit, flagged_index, ok_index
end

# ╔═╡ 45ff666f-3421-4367-bb50-4766df16efa6
"""
	fit_models_th(data::Array{DataFrame}, res_threshold)

Fit the ABC-model and the K-centric-model at the lnk(T) data of the `data`-array of dataframes, using LsqFit.jl 
"""
function fit_models_th(data::Array{DataFrame}, β0::Array{Float64}, res_threshold)
	fit = Array{DataFrame}(undef, length(data))
	for i=1:length(data)
		fitABC = Array{Any}(undef, length(data[i][!,1]))
		fitKcentric = Array{Any}(undef, length(data[i][!,1]))
		T = Array{Array{Float64}}(undef, length(data[i][!,1]))
		lnk = Array{Any}(undef, length(data[i][!,1]))
		name = Array{String}(undef, length(data[i][!,1]))
		excludedABC_i = Array{Array{Int,1}}(undef, length(data[i][!,1]))
		okABC_i = Array{Array{Int,1}}(undef, length(data[i][!,1]))
		excludedKcentric_i = Array{Array{Int,1}}(undef, length(data[i][!,1]))
		okKcentric_i = Array{Array{Int,1}}(undef, length(data[i][!,1]))
		R²_ABC = Array{Float64}(undef, length(data[i][!,1]))
		R²_Kcentric = Array{Float64}(undef, length(data[i][!,1]))
		for j=1:length(data[i][!,1])
			T[j] = T_column_names_to_Float(data[i])
			lnk[j] = collect(Union{Missing,Float64}, data[i][j,2:end])
			
			name[j] = data[i][!, 1][j]

			fitABC[j], excludedABC_i[j], okABC_i[j] = fit_ABC(T[j], lnk[j], β0[i], res_threshold)
			fitKcentric[j], excludedKcentric_i[j], okKcentric_i[j] = fit_Kcentric(T[j], lnk[j], res_threshold)

			R²_ABC[j] = coeff_of_determination(fitABC[j], lnk[j][okABC_i[j]])
			R²_Kcentric[j] = coeff_of_determination(fitKcentric[j], lnk[j][okKcentric_i[j]])
		end
		fit[i] = DataFrame(Name=name, T=T, lnk=lnk, fitABC=fitABC, fitKcentric=fitKcentric, i_ABC=okABC_i, i_Kcentric=okKcentric_i, ex_i_ABC=excludedABC_i, ex_i_Kcentric=excludedKcentric_i, R²_ABC=R²_ABC, R²_Kcentric=R²_Kcentric)
	end
	return fit
end

# ╔═╡ c7d96000-5a3c-4c46-85e0-ce0c101c4546
"""
	fit_models_th!(meta_data::DataFrame, β0::Array{Float64}, res_threshold)

Fit the ABC-model and the K-centric-model at the lnk(T) data of the `data`-array of dataframes, using LsqFit.jl and add the result in a new column (`fit`) of `meta_data 
"""
function fit_models_th!(meta_data::DataFrame, res_threshold)
	fit = fit_models_th(meta_data.data, meta_data.beta0, res_threshold)
	meta_data[!, "fit"] = fit
	return meta_data
end

# ╔═╡ 487afef9-9e37-48e5-bbd7-4031e6cc3a5e
fitt = fit_models_th(meta_data.data, meta_data.beta0, 0.01)

# ╔═╡ 21195058-acdd-4a49-accc-fcb6d1e11c0d
fitt[2]

# ╔═╡ b9d90721-44f1-40cd-b663-84905787063f
fitt[2].fitABC

# ╔═╡ 600feb4b-bf68-48f7-a46a-fa3cc205b808
fitt[2]

# ╔═╡ 600fd676-0677-49c8-920d-2d675ff2c64b
mean(fitt[2].lnk[5][fitt[2].i_ABC[5]])

# ╔═╡ 73511827-56ad-48a3-b514-8fa59c87366c
coeff_of_determination(fit[2].fitABC[5], fitt[2].lnk[5][fitt[2].i_ABC[5]])

# ╔═╡ 5204cc96-f477-4ae4-bb27-09c94b635a07
fitt[3].fitABC[5].param

# ╔═╡ dd7ef4db-5843-4bbc-a2c3-bb53a23917d8
stderror(fitt[3].fitABC[5])

# ╔═╡ a1b16e67-4253-4e4a-98df-66671410458b
fit_ = curve_fit(ABC, T_[ok_index].+Tst, lnk_[ok_index], ABC0)

# ╔═╡ d736e58c-2025-4d35-b929-8ab179c03ad9
maximum(abs.(fit_.resid))

# ╔═╡ fef97ae3-2350-4525-b8b4-148cd0613a6a
maximum(abs.(fit_.resid)) > threshold

# ╔═╡ e6526360-bb38-4d4f-b9a2-f0ec72801472
imax = findfirst(abs.(fit_.resid).==maximum(abs.(fit_.resid)))

# ╔═╡ c2c9a1d4-522a-4cff-9286-77f4595a8721
push!(flagged_index, ok_index[imax])

# ╔═╡ bf6b2c8b-c3c3-42f4-8c54-dc4279bdc0d6
ok_index_ = ok_index[findall(ok_index.!=imax)]

# ╔═╡ 75295c0c-d703-4143-b412-1e9da1bd2483
length(ok_index_) <= 3

# ╔═╡ 82daad7c-c942-4ba8-bf90-a31f4de16892
fit__ = curve_fit(ABC, T_[ok_index_].+Tst, lnk_[ok_index_], ABC0)

# ╔═╡ 967e45f9-767d-4111-87ed-d20fda066004
maximum(abs.(fit__.resid)) > threshold

# ╔═╡ 87682e01-1b20-455d-96cc-6ce4b573a685
imax_ = findfirst(abs.(fit__.resid).==maximum(abs.(fit__.resid)))

# ╔═╡ 29604d6b-5b46-41e2-aeb0-d83855cd9dd3
push!(flagged_index, ok_index_[imax_])

# ╔═╡ 3bb6eff9-7a90-4f9f-9a6b-3a53729b4d70
ok_index__ = ok_index_[findall(ok_index_.!=ok_index_[imax_])]

# ╔═╡ 4b184f2f-58f3-4437-bec9-f2f24d584be7
fit___ = curve_fit(ABC, T_[ok_index__].+Tst, lnk_[ok_index__], ABC0)

# ╔═╡ c178fa9c-1a8a-413b-9572-e4de71a523b1
maximum(abs.(fit___.resid)) > threshold

# ╔═╡ 21da172e-409c-4c50-85d1-25a9169208ff
"""
	fit_models_w(data::Array{DataFrame}, β0::Array{Float64})

Fit the ABC-model and the K-centric-model at the lnk(T) data of the `data`-array of dataframes, using LsqFit.jl 
"""
function fit_models_w(data::Array{DataFrame}, β0::Array{Float64})
	fit = Array{DataFrame}(undef, length(data))
	for i=1:length(data)
		fitABC = Array{Any}(undef, length(data[i][!,1]))
		fitKcentric = Array{Any}(undef, length(data[i][!,1]))
		T = Array{Array{Float64}}(undef, length(data[i][!,1]))
		lnk = Array{Any}(undef, length(data[i][!,1]))
		name = Array{String}(undef, length(data[i][!,1]))
		ok_i = Array{Array{Int,1}}(undef, length(data[i][!,1]))
		R²_ABC = Array{Float64}(undef, length(data[i][!,1]))
		R²_Kcentric = Array{Float64}(undef, length(data[i][!,1]))
		for j=1:length(data[i][!,1])
			T[j] = T_column_names_to_Float(data[i])
			lnk[j] = collect(Union{Missing,Float64}, data[i][j,2:end])
			ii = findall(isa.(lnk[j], Float64))
			name[j] = data[i][!, 1][j]
			ABC0 = [-100.0, 10000.0, 10.0] # for p[1] = -100.0 - log(beta0) ?
			lb_ABC = [-Inf, 0.0, 0.0]
			ub_ABC = [-log(β0[i]), Inf, Inf]
			Kcentric0 = [200.0+Tst, 30.0, 10.0]
			lb_Kcentric = [0.0, 0.0, 0.0]
			ub_Kcentric = [Inf, Inf, 500.0]
			fitABC0 = curve_fit(ABC, T[j][ii].+Tst, lnk[j][ii], ABC0, lower=lb_ABC, upper=ub_ABC)
			# estimate weights from the residuals of the ordinary least squared
			w_ABC = 1.0./fitABC0.resid.^2
			fitABC[j] = curve_fit(ABC, T[j][ii].+Tst, lnk[j][ii], w_ABC, ABC0, lower=lb_ABC, upper=ub_ABC)
			
			fitKcentric0 = curve_fit(Kcentric, T[j][ii].+Tst, lnk[j][ii], Kcentric0, lower=lb_Kcentric, upper=ub_Kcentric)
			w_Kcentric = 1.0./fitKcentric0.resid.^2
			fitKcentric[j] = curve_fit(Kcentric, T[j][ii].+Tst, lnk[j][ii], w_Kcentric, Kcentric0, lower=lb_Kcentric, upper=ub_Kcentric)
			ok_i[j] = ii
			R²_ABC[j] = coeff_of_determination(fitABC[j], lnk[j][ii])
			R²_Kcentric[j] = coeff_of_determination(fitKcentric[j], lnk[j][ii])
		end
		fit[i] = DataFrame(Name=name, T=T, lnk=lnk, fitABC=fitABC, fitKcentric=fitKcentric, ok_i=ok_i, R²_ABC=R²_ABC, R²_Kcentric=R²_Kcentric)
	end
	return fit
end

# ╔═╡ c1dc31ac-7d0b-414e-8980-a5113a68e779
"""
	fit_models_w!(meta_data::DataFrame)

Fit the ABC-model and the K-centric-model at the lnk(T) data of the `data`-array of dataframes, using LsqFit.jl and add the result in a new column (`fit`) of `meta_data 
"""
function fit_models_w!(meta_data::DataFrame)
	fit = fit_models_w(meta_data.data, meta_data.beta0)
	meta_data[!, "fit"] = fit
	return meta_data
end

# ╔═╡ fc9d9040-6033-4660-8549-7f8394bfb895
fit_w = fit_models_w(meta_data.data, meta_data.beta0)

# ╔═╡ 923c17a3-128b-4f9d-924d-066a2f6a0ee1
fit_w[2]

# ╔═╡ 27802ac2-dcca-4b0f-8027-3db6c9dc26fa
fit_w[2].fitABC[5]

# ╔═╡ 63b395c7-c3b6-466a-93ee-0493cc65f708
fit_w[2].fitKcentric[5]

# ╔═╡ 97a3bc33-954f-43ee-bdd3-82a4f7e766b0
sum(fit_w[2].fitKcentric[5].wt)

# ╔═╡ 60b5df40-2079-48de-8fa0-0fd187f59dba
fit_w[3].fitABC[5].param

# ╔═╡ 32248d5d-bc61-4b39-9f67-aa72f1505eaf
stderror(fit_w[3].fitABC[5])

# ╔═╡ b1b39ba3-9c2d-4dc2-be70-e2b099349397
"""
	extract_paramaters_from_fit(fit, β0)

Extract the parameters `A`, `B`, `C`, `Tchar`, `θchar` and `ΔCp` from the fits of lnk over T data with the `ABC`-model and the `Kcentric`-model.
"""
function extract_parameters_from_fit(fit, β0)
	Par = Array{DataFrame}(undef, length(fit))
	for i=1:length(fit)
		A = Array{Measurement{Float64}}(undef, length(fit[i].Name))
		B = Array{Measurement{Float64}}(undef, length(fit[i].Name))
		C = Array{Measurement{Float64}}(undef, length(fit[i].Name))
		Tchar = Array{Measurement{Float64}}(undef, length(fit[i].Name))
		θchar = Array{Measurement{Float64}}(undef, length(fit[i].Name))
		ΔCp = Array{Measurement{Float64}}(undef, length(fit[i].Name))
		ΔHref = Array{Measurement{Float64}}(undef, length(fit[i].Name))
		ΔSref = Array{Measurement{Float64}}(undef, length(fit[i].Name))
		for j=1:length(fit[i].Name)
			A[j] = (fit[i].fitABC[j].param[1] ± stderror(fit[i].fitABC[j])[1]) + log(β0[i])
			B[j] = fit[i].fitABC[j].param[2] ± stderror(fit[i].fitABC[j])[2]
			C[j] = fit[i].fitABC[j].param[3] ± stderror(fit[i].fitABC[j])[3]

			Tchar[j] = (fit[i].fitKcentric[j].param[1] ± stderror(fit[i].fitKcentric[j])[1]) - Tst
			θchar[j] = fit[i].fitKcentric[j].param[2] ± stderror(fit[i].fitKcentric[j])[2]
			ΔCp[j] = (fit[i].fitKcentric[j].param[3] ± stderror(fit[i].fitKcentric[j])[3]) * R
			TD = ThermodynamicData.ABC_to_TD(A[j], B[j], C[j], T0)
			ΔHref[j] = TD[1]
			ΔSref[j] = TD[2]
		end
		Par[i] = DataFrame(Name=fit[i].Name, A=A, B=B, C=C, Tchar=Tchar, thetachar=θchar, DeltaCp=ΔCp, DeltaHref=ΔHref, DeltaSref=ΔSref, R²_ABC=fit[i].R²_ABC, R²_Kcentric=fit[i].R²_Kcentric)
	end
	return Par
end			

# ╔═╡ e764b126-a0d9-42a9-a177-1bd460c05774
extract_parameters_from_fit(fit, meta_data.beta0)[1]

# ╔═╡ 882889aa-1a19-4741-8e5f-934cd60fdf49
extract_parameters_from_fit(fitt, meta_data.beta0)[1]

# ╔═╡ 1d2fd518-bdad-4ad1-8038-07229b89c447
extract_parameters_from_fit(fit_w, meta_data.beta0)[2]

# ╔═╡ 4a89410f-a205-4f8b-9af8-47fc5563e689
b[1] - Tst

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ChemicalIdentifiers = "fa4ea961-1416-484e-bda2-883ee1634ba5"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LambertW = "984bce1d-4616-540c-a9ee-88d1112d94c9"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.2"
ChemicalIdentifiers = "~0.1.5"
DataFrames = "~1.3.2"
LambertW = "~0.4.5"
LsqFit = "~0.12.1"
Measurements = "~2.7.1"
Plots = "~1.25.11"
PlutoUI = "~0.7.35"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "d49f55ff9c7ee06930b0f65b1df2bfa811418475"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.4"

[[deps.Arrow]]
deps = ["ArrowTypes", "BitIntegers", "CodecLz4", "CodecZstd", "DataAPI", "Dates", "Mmap", "PooledArrays", "SentinelArrays", "Tables", "TimeZones", "UUIDs"]
git-tree-sha1 = "85013d248b128cf13ae62c827c4bf05872e97f78"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "2.2.1"

[[deps.ArrowTypes]]
deps = ["UUIDs"]
git-tree-sha1 = "a0633b6d6efabf3f76dacd6eb1b3ec6c42ab0552"
uuid = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
version = "1.2.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitIntegers]]
deps = ["Random"]
git-tree-sha1 = "5a814467bda636f3dde5c4ef83c30dd0a19928e0"
uuid = "c3b6d118-76ef-56ca-8cc7-ebb389d030a1"
version = "0.2.6"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "9519274b50500b8029973d241d32cfbf0b127d97"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.2"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ChemicalIdentifiers]]
deps = ["Arrow", "Downloads", "Preferences", "Scratch", "UUIDs", "Unicode"]
git-tree-sha1 = "708a12392479dd484be472725a005529ca8f90a0"
uuid = "fa4ea961-1416-484e-bda2-883ee1634ba5"
version = "0.1.5"

[[deps.CodecLz4]]
deps = ["Lz4_jll", "TranscodingStreams"]
git-tree-sha1 = "59fe0cb37784288d6b9f1baebddbf75457395d40"
uuid = "5ba52731-8f18-5e0d-9241-30f10d1ec561"
version = "0.4.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.CodecZstd]]
deps = ["CEnum", "TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "849470b337d0fa8449c21061de922386f32949d9"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9d3c0c762d4666db9187f363a76b47f7346e673b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.49"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "84f04fe68a3176a583b864e492578b9466d87f1e"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "4c7d3757f3ecbcb9055870351078552b7d1dbd2d"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.0"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "ec299fdc8f49ae450807b0cb1d161c6b76fd2b60"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.10.1"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LambertW]]
git-tree-sha1 = "2d9f4009c486ef676646bca06419ac02061c088e"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.5"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a6552bfeab40de157a297d84e03ade4b8177677f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.12"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "91aa1442e63a77f101aff01dec5a821a17f43922"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.12.1"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "88cd033eb781c698e75ae0b680e5cef1553f0856"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.7.1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "7e2166042d1698b6072352c74cfd1fca2a968253"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.6"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "5c907bdee5966a9adb8a106807b7c387e51e4d6c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.11"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "85bf3e4bd279e405f91489ce518dedb1e32119cb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.35"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "65068e4b4d10f3c31aaae2e6cb92b6c6cedca610"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "2d4b6de8676b34525ac518de36006dc2e89c7e2e"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.7.2"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═ae71c5ce-9978-11ec-3f73-1197d152796a
# ╠═aaa3793e-b9c8-4742-8005-33c336787bf7
# ╠═0e80aed9-f0ac-4077-b1f4-b298349d31c8
# ╠═413fad2a-e73f-4933-9eed-b1aae8d376ce
# ╠═e2566bfc-1912-45bb-90c4-b05ee59c16af
# ╠═ffa74790-bb57-4dc9-9bfc-8309d6d58f49
# ╠═34c6c2bd-63a3-4af1-b6ac-1dccdf3ddb4c
# ╠═15047914-7463-4d22-934a-00f7c81ee44a
# ╠═baf00ccd-0a30-4cb4-9704-d5d9b23cfe69
# ╠═422b832b-0b52-4c35-995b-a904069d42c5
# ╠═29232c11-9e21-42aa-aca3-500563bbc9a9
# ╠═9ea6d953-8556-4ddc-806b-02a74a3b94d8
# ╠═fbd168b6-8d91-4e60-b2ac-342e787ccb14
# ╠═a356e7c7-c1da-4907-b8b4-e7cf7c6f5c57
# ╠═035a63b2-1d9b-417b-ae3a-f27604a50898
# ╠═b946b6a7-4453-44a6-9419-1ecd235d2956
# ╠═5655ccd8-09db-40b9-9b6d-f88b5d328d43
# ╠═77b226e1-4c36-47e7-b181-511667de9968
# ╠═e011d3c1-7070-4e4c-b804-d63424aa6f20
# ╠═9025d205-965a-4843-a0b9-eeaa8d98d997
# ╠═1c38490d-c6b6-443f-b856-2b33e8417277
# ╠═2028b40d-13a5-4c22-8335-8691bbe5a0ee
# ╠═f754b0d2-2f63-40fa-92c8-a692f52975c3
# ╠═f99640cf-58e9-4bb9-997b-e76578edebfa
# ╠═1d466a06-3865-477e-bc10-d811eb6a2732
# ╠═4a2d1ecb-02b6-44a5-953e-eedc6cef512f
# ╠═d83d45da-f852-424e-b3aa-2a5881de28fe
# ╠═049e3387-eead-4806-acdc-77ff0e292d53
# ╠═a10bec56-173e-460c-8eb5-040f3da1244c
# ╠═edda5af7-5ca1-4456-b5eb-2a2ff3ad89d4
# ╠═125495fc-1609-4b80-a1f3-ade4104f7da2
# ╠═996a3ca9-ede4-4ecf-aa24-38f3fb104cb2
# ╠═126dceaf-7e01-47a4-81ce-f83f799ac422
# ╠═13949a69-a1fc-4b81-9eb6-66269e1bf853
# ╠═4c2eed13-f7fc-4e0e-ad64-e8bf2894e98e
# ╠═d7f8407c-5ad2-4df4-a222-53f2c478fd94
# ╠═9c66eca0-9ef7-4e2e-8c37-0d2c522df5b0
# ╠═52fe47c1-d030-4936-b627-8b7886b6211b
# ╠═e7113c1b-e465-4983-a01c-52f00ed6b16f
# ╠═522d5f26-9ad5-42ce-83f2-fdabfd5d8313
# ╠═761f3c9f-7561-499b-916d-767bc3f2aafb
# ╠═45ff666f-3421-4367-bb50-4766df16efa6
# ╠═c7d96000-5a3c-4c46-85e0-ce0c101c4546
# ╠═e7dcf5dd-ab2d-4e45-ba9c-50d55896beb2
# ╠═f404511b-41d1-4cac-a260-9409184becdd
# ╠═04a0772d-d6e7-4661-88be-45f8e06d2aff
# ╠═8c2ec059-872e-4e58-851c-e3628b477b33
# ╠═3624acaf-c9f8-4548-88da-14a22e3e054f
# ╠═caddeaba-0ad6-45a0-8af0-2fa61d09663c
# ╠═69e663e4-3814-419b-9225-1574eedb9b54
# ╠═a1b16e67-4253-4e4a-98df-66671410458b
# ╠═d736e58c-2025-4d35-b929-8ab179c03ad9
# ╠═14d85174-45f2-459f-82a3-498c9aa5b1b8
# ╠═fef97ae3-2350-4525-b8b4-148cd0613a6a
# ╠═30f10519-ff73-4d46-8837-070d9c46cdbf
# ╠═e6526360-bb38-4d4f-b9a2-f0ec72801472
# ╠═c2c9a1d4-522a-4cff-9286-77f4595a8721
# ╠═bf6b2c8b-c3c3-42f4-8c54-dc4279bdc0d6
# ╠═82daad7c-c942-4ba8-bf90-a31f4de16892
# ╠═967e45f9-767d-4111-87ed-d20fda066004
# ╠═75295c0c-d703-4143-b412-1e9da1bd2483
# ╠═87682e01-1b20-455d-96cc-6ce4b573a685
# ╠═29604d6b-5b46-41e2-aeb0-d83855cd9dd3
# ╠═3bb6eff9-7a90-4f9f-9a6b-3a53729b4d70
# ╠═4b184f2f-58f3-4437-bec9-f2f24d584be7
# ╠═c178fa9c-1a8a-413b-9572-e4de71a523b1
# ╠═487afef9-9e37-48e5-bbd7-4031e6cc3a5e
# ╠═21195058-acdd-4a49-accc-fcb6d1e11c0d
# ╠═b9d90721-44f1-40cd-b663-84905787063f
# ╠═578bade5-7af7-4dae-815e-41ab1c9311d5
# ╠═0ac891c1-4ef5-4267-bfa3-13e79c316534
# ╠═be15dac4-309f-459a-b86b-c80fb12421bc
# ╠═5335f894-8b7a-41e5-a0e5-292bdde60881
# ╠═e803eef6-5b7e-4f8b-a89c-0ea99117cf35
# ╠═0a5b2bc6-9dd0-4483-af64-50e889a7a1c4
# ╠═f7e45701-3ce6-49fa-942a-90102c82c784
# ╠═9b1271c5-7686-4f23-9f31-4a7028afe1b1
# ╠═c3e1408f-9cf7-4418-8b66-fe3e502a8278
# ╠═bfbaf273-8993-4588-8822-220efc242bb9
# ╠═1d273a65-577d-42ad-9534-7acf6ca1e176
# ╠═70b59e3b-373f-455a-b58d-3083cacf3c3b
# ╠═600feb4b-bf68-48f7-a46a-fa3cc205b808
# ╠═600fd676-0677-49c8-920d-2d675ff2c64b
# ╠═ba2dfc06-37df-44df-b8f4-50a506af4151
# ╠═73511827-56ad-48a3-b514-8fa59c87366c
# ╠═17367fcc-8737-42a5-9546-4d7447420248
# ╠═cece68a2-86bf-49f1-bfb0-5dcc8a41eb5b
# ╠═bd6acff3-93b1-401f-9ec3-2d1d84485a4d
# ╠═d9494de6-248e-4731-93ba-5cce35fcf0de
# ╠═21da172e-409c-4c50-85d1-25a9169208ff
# ╠═c25f377d-cd5c-4626-a09b-f0043e478b40
# ╠═09ed6ad1-f161-4288-9152-1be342267b5e
# ╠═c1dc31ac-7d0b-414e-8980-a5113a68e779
# ╠═fc9d9040-6033-4660-8549-7f8394bfb895
# ╠═923c17a3-128b-4f9d-924d-066a2f6a0ee1
# ╠═27802ac2-dcca-4b0f-8027-3db6c9dc26fa
# ╠═94bbe840-f321-4699-9a57-21a4291cf56e
# ╠═63b395c7-c3b6-466a-93ee-0493cc65f708
# ╠═cf44ed56-ef17-4bb6-a131-82bfd405a3fe
# ╠═14bd0e24-c16d-45bc-aaa3-d0775220374b
# ╠═5ac44f86-f9c3-4a59-9764-11c3aa0a9e66
# ╠═8ff63c84-7996-47e0-87ef-d00ea1eafe93
# ╠═97a3bc33-954f-43ee-bdd3-82a4f7e766b0
# ╠═6743fce9-f000-4757-aca5-f347e4c38ef7
# ╠═252713b0-bebc-4a56-a154-44902fcee8bf
# ╠═88376f08-7ee9-43e2-a0d8-0e2dc4e63904
# ╠═60b5df40-2079-48de-8fa0-0fd187f59dba
# ╠═32248d5d-bc61-4b39-9f67-aa72f1505eaf
# ╠═5204cc96-f477-4ae4-bb27-09c94b635a07
# ╠═dd7ef4db-5843-4bbc-a2c3-bb53a23917d8
# ╠═e43a548b-a863-4ff1-ba76-7c4392175f5b
# ╠═97e9dd89-bf47-4c79-ab41-84e54bb0ec2c
# ╠═32f2b9b2-0982-41b9-b297-845f1829d301
# ╠═646b3e92-6725-4363-9c76-e94c0d84af6c
# ╠═4a1b443e-1965-432f-abb7-886edce3a7f0
# ╠═c61f2185-5c78-4cd9-81b5-2c9e472d5c37
# ╠═e270cf70-f13f-4fc5-91d4-b5181621f9da
# ╠═b1b39ba3-9c2d-4dc2-be70-e2b099349397
# ╠═9efc51b7-10ab-4b15-bdab-be7bb3f55872
# ╠═4a89410f-a205-4f8b-9af8-47fc5563e689
# ╠═ea5a5cd1-e53a-47f8-ad8d-b46aa3d61a81
# ╠═e764b126-a0d9-42a9-a177-1bd460c05774
# ╠═882889aa-1a19-4741-8e5f-934cd60fdf49
# ╠═1d2fd518-bdad-4ad1-8038-07229b89c447
# ╠═9ad6c6da-27e4-4f22-b693-1689824aad0c
# ╠═1f252318-8645-4be0-a900-a9585f1c73ff
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
