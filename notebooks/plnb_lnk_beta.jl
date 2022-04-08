### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ e3da9b92-eb94-11eb-0c86-0f1b7c98a586
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	using DrWatson
	using VGGC
	using PlutoUI
	using DataFrames, CSV, Interpolations, Plots, LsqFit, Measurements, Statistics
	plotly()
	#html"""
	#<style>
	#  main {
	#	max-width: 1000px;
	#  }
	#</style>
	#"""
end

# ╔═╡ c7679cc4-06c4-41af-8695-47436330b06f
md"""
# Investigation of the theoretical correlation of the K-centric thermodynamic parameters and the phase ratio β 
"""

# ╔═╡ 0875234e-ab67-47ad-8a86-b9121ba5e7d4
begin
	# parameters
	sp = "Rxi17SilMS"
	solute_db_path = "/Users/janleppert/Documents/GitHub/VGGC/data/exp_pro/Databases"
	solute_db = "Database_append.csv"
	gas = "He"
	T = (40.0:10.0:360.0) .+ 273.15
	β = [62.5, 125, 250, 625, 1250]
	R = 8.3145
end

# ╔═╡ b43aa25c-a502-4eda-8ef3-4f7b8b87e562
function all_solutes(sp, db)
	db_filter = filter([:Phase] => x -> x==sp, db)
	solutes = string.(db_filter.Name)
	return solutes
end

# ╔═╡ ec472f66-78ee-4cf5-a368-fe083cd67c06
begin
	db = DataFrame(CSV.File(string(solute_db_path,"/",solute_db), header=1, silencewarnings=true))
	solutes = all_solutes(sp, db)
	sub = VGGC.load_solute_database(solute_db_path, solute_db, sp, gas, solutes, zeros(length(solutes)), zeros(length(solutes)))
end

# ╔═╡ 8917fa51-228e-4b43-b3c2-f83a4aa79ce5
begin
	# calculation of the new Tchar parameter at another phase ratio
	using LambertW
	y = Array{Float64}(undef, length(sub), length(β))
	W = Array{Float64}(undef, length(sub), length(β))
	T0 = Array{Float64}(undef, length(sub), length(β))
	β₀ = 250
	for j=1:length(β)
		for i=1:length(sub)
			a = sub[i].ΔCp/R*sub[i].θchar
			expon = (log(β₀)-log(β[j]))*R/sub[i].ΔCp - sub[i].Tchar/a - 1
			y[i,j] = -(a + sub[i].Tchar)/a*exp(expon)
			W[i,j] = lambertw(y[i,j],-1)
			T0[i,j] = -sub[i].Tchar*(a+sub[i].Tchar)/(a*W[i,j])
		end
	end
	# this calculation works!
end

# ╔═╡ 8731c46e-594a-40bb-a262-ee1d0e87f1a8
begin
	# calculate the retention factors
	lnk = Array{Float64}(undef, length(sub), length(T), length(β))
	for k=1:length(β)
		for j=1:length(T)
			for i=1:length(sub)
				lnk₀ = (sub[i].ΔCp/R + sub[i].Tchar/sub[i].θchar)*(sub[i].Tchar/T[j]-1) + sub[i].ΔCp/R*log(T[j]/sub[i].Tchar)
				β₀ = 1/(4*sub[i].φ₀)
				lnk[i,j,k] = log(β₀) - log(β[k]) + lnk₀
			end
		end
	end
end

# ╔═╡ e2684d83-3251-4245-bc74-e6107888b9e3
lnk

# ╔═╡ df388e08-face-436b-8b82-0ee72e83d845
begin
	# plot ln(k) for a selected solute for all β over T
	ii = 40
	p_lnk = scatter(T, lnk[ii,:,1], xlabel="Temperature T in K", ylabel="ln(k)", label="$(sub[ii].name), β=$(β[1])")
	for j=2:length(β)
		scatter!(p_lnk, T, lnk[ii,:,j], label="$(sub[ii].name), β=$(β[j])")
	end
	p_lnk
end

# ╔═╡ 62c16f38-46bf-47e4-b602-60718d47fcaf
function Rsquare(fit, y)
		sstot = sum((y.-mean(y)).^2)
		ssres = sum(fit.resid.^2)
		R2 = 1-ssres/sstot
		return R2
	end

# ╔═╡ 397e2565-5ea4-4c97-8d43-684d0f2dee70
begin
	# fit K-centric-model to ln(k)
	model(x,p) = p[1] .*(p[2]./x .- 1) .+ p[3].*log.(x./p[2])
	p0 = [sub[ii].ΔCp/R+sub[ii].Tchar/sub[ii].θchar, sub[ii].Tchar, sub[ii].ΔCp/R]
	fit = Array{Any}(undef, length(β))
	p_lnk_fit = scatter(xlabel="Temperature T in K", ylabel="ln(k)")
	p1 = Array{Measurement{Float64}}(undef, length(β))
	p2 = Array{Measurement{Float64}}(undef, length(β))
	p3 = Array{Measurement{Float64}}(undef, length(β))
	R² = Array{Float64}(undef, length(β))
	for j=1:length(β)
		fit[j] = curve_fit(model, T, lnk[ii,:,j], p0)
		p1[j] = measurement(fit[j].param[1], standard_errors(fit[j])[1])
		p2[j] = measurement(fit[j].param[2], standard_errors(fit[j])[2])
		p3[j] = measurement(fit[j].param[3], standard_errors(fit[j])[3])
		R²[j] = Rsquare(fit[j], model(T, fit[j].param))
		scatter!(p_lnk_fit, T, lnk[ii,:,j], label="$(sub[ii].name), β=$(β[j])")
		plot!(p_lnk_fit, T, model(T, fit[j].param), label="fit β=$(β[j]), R²=$(R²[j])")	
	end
	p_lnk_fit
end

# ╔═╡ cc89b0a9-38eb-4cb9-8d66-7281258905f8
begin
	# fit K-centric-model to ln(k) for all solutes
	fit_all = Array{Any}(undef, length(sub), length(β))
	p1_all = Array{Measurement{Float64}}(undef, length(sub), length(β))
	p2_all = Array{Measurement{Float64}}(undef, length(sub), length(β))
	p3_all = Array{Measurement{Float64}}(undef, length(sub), length(β))
	R²_all = Array{Float64}(undef, length(sub), length(β))
	for j=1:length(β)
		for i=1:length(sub)
			p00 = [sub[i].ΔCp/R+sub[i].Tchar/sub[i].θchar, sub[i].Tchar, sub[i].ΔCp/R]
			fit_all[i,j] = curve_fit(model, T, lnk[i,:,j], p00)
			p1_all[i,j] = measurement(fit_all[i,j].param[1], standard_errors(fit_all[i,j])[1])
			p2_all[i,j] = measurement(fit_all[i,j].param[2], standard_errors(fit_all[i,j])[2])
			p3_all[i,j] = measurement(fit_all[i,j].param[3], standard_errors(fit_all[i,j])[3])
			R²_all[i,j] = Rsquare(fit_all[i,j], model(T, fit_all[i,j].param))
		end
	end
end

# ╔═╡ 914876b4-e3ec-4558-95fb-fabd87db0a41
begin
	d = 0.25*1e-3
	df = 1.0*1e-6
	βexakt = (d-2*df)^2/(d^2-(d-2*df)^2)
	βrund = d/(4*df)
	βexakt, βrund, log(βrund)
end

# ╔═╡ dea5b4d0-78f9-467a-905c-78f9795c4b86
# Tchar changes with ln(β) nearly linear
scatter(log.(β), p2_all[ii,:], label="$(sub[ii].name)", xlabel="ln(β)", ylabel="Tchar in K")

# ╔═╡ 04ccde67-0171-4069-9c1f-eb88bd0d2e33
scatter(log.(β), p1_all[ii,:], label="$(sub[ii].name)", xlabel="ln(β)", ylabel="p1")

# ╔═╡ 1c38f24c-657f-4860-8d6f-4e3b3e15d29e
# ΔCp/R is nearly constant
scatter(log.(β), p3_all[ii,:], label="$(sub[ii].name)", xlabel="ln(β)", ylabel="ΔCp/R")

# ╔═╡ 28e42e86-660d-4cac-8121-71adf79a7c83
R²_all

# ╔═╡ afad4f9b-c0b6-4dd7-bf66-0d23deba7e87
begin
	p_Tcharβ = scatter(xlabel="ln(β)", ylabel="Tchar in K", legend=false)
	model_lin(x,p) = p[1] .+ p[2].*x
	p_lin_0 = [1000.0, -1.01]
	fit_lin = Array{Any}(undef, length(sub))
	p1_lin = Array{Measurement{Float64}}(undef, length(sub))
	p2_lin = Array{Measurement{Float64}}(undef, length(sub))
	R²_lin = Array{Float64}(undef, length(sub))
	for i=1:length(sub)
		scatter!(p_Tcharβ, log.(β), p2_all[i,:], label="$(sub[i].name)")
		fit_lin[i] = curve_fit(model_lin, log.(β), Measurements.value.(p2_all[i,:]), p_lin_0)
		p1_lin[i] = measurement(fit_lin[i].param[1], standard_errors(fit_lin[i])[1])
		p2_lin[i] = measurement(fit_lin[i].param[2], standard_errors(fit_lin[i])[2])
		R²_lin[i] = Rsquare(fit_lin[i], model_lin(log.(β), fit_lin[i].param))
		plot!(p_Tcharβ, log.(β), model_lin(log.(β), fit_lin[i].param))	
	end
	p_Tcharβ
end

# ╔═╡ 1af6888d-200d-40da-bece-2ba53ba21aae
fitt = curve_fit(model_lin, log.(β), Measurements.value.(p2_all[ii,:]), [800.0, -0.05])

# ╔═╡ 45dfa9f0-5282-4a59-9fd0-28521a6976c2
R²_lin

# ╔═╡ 93c4de7c-b25c-4d52-9a48-33f9b6a10145
scatter(1:97, p2_lin)

# ╔═╡ cfc400d2-ee2b-4a23-b8c8-b19bd2d74b59
begin
	# fit K-centric-model to ln(k) for all solutes with ΔCp fixed
	fit_ΔCpfixed = Array{Any}(undef, length(sub), length(β))
	p1_ΔCpfixed = Array{Measurement{Float64}}(undef, length(sub), length(β))
	p2_ΔCpfixed = Array{Measurement{Float64}}(undef, length(sub), length(β))
	R²_ΔCpfixed = Array{Float64}(undef, length(sub), length(β))
	for j=1:length(β)
		for i=1:length(sub)
			model_ΔCpfixed(x,p) = (sub[i].ΔCp/R .+ p[2]./p[1]).*(p[2]./x .- 1) .+ sub[i].ΔCp/R.*log.(x./p[2])
			fit_ΔCpfixed[i,j] = curve_fit(model_ΔCpfixed, T, lnk[i,:,j], [sub[i].θchar, sub[i].Tchar])
			p1_ΔCpfixed[i,j] = measurement(fit_ΔCpfixed[i,j].param[1], standard_errors(fit_ΔCpfixed[i,j])[1])
			p2_ΔCpfixed[i,j] = measurement(fit_ΔCpfixed[i,j].param[2], standard_errors(fit_ΔCpfixed[i,j])[2])
			R²_ΔCpfixed[i,j] = Rsquare(fit_ΔCpfixed[i,j], model_ΔCpfixed(T, fit_ΔCpfixed[i,j].param))
		end
	end
end

# ╔═╡ 7301cc6c-b38a-4cee-b566-db360f132b8c
R²_ΔCpfixed

# ╔═╡ 228b97c6-4230-459d-bc19-567045201ec6
begin
	# with fixed ΔCp
	p_Tcharβ_ΔCpfixed = scatter(xlabel="ln(β)", ylabel="Tchar in K", legend=false)
	fit_Tcharβ_ΔCpfixed = Array{Any}(undef, length(sub))
	p1_Tcharβ_ΔCpfixed = Array{Measurement{Float64}}(undef, length(sub))
	p2_Tcharβ_ΔCpfixed = Array{Measurement{Float64}}(undef, length(sub))
	R²_Tcharβ_ΔCpfixed = Array{Float64}(undef, length(sub))
	for i=1:length(sub)
		scatter!(p_Tcharβ_ΔCpfixed, log.(β), p2_ΔCpfixed[i,:], label="$(sub[i].name)")
		fit_Tcharβ_ΔCpfixed[i] = curve_fit(model_lin, log.(β), Measurements.value.(p2_ΔCpfixed[i,:]), [1000.0, -30.0])
		p1_Tcharβ_ΔCpfixed[i] = measurement(fit_Tcharβ_ΔCpfixed[i].param[1], standard_errors(fit_Tcharβ_ΔCpfixed[i])[1])
		p2_Tcharβ_ΔCpfixed[i] = measurement(fit_Tcharβ_ΔCpfixed[i].param[2], standard_errors(fit_Tcharβ_ΔCpfixed[i])[2])
		R²_Tcharβ_ΔCpfixed[i] = Rsquare(fit_Tcharβ_ΔCpfixed[i], model_lin(log.(β), fit_Tcharβ_ΔCpfixed[i].param))
		plot!(p_Tcharβ_ΔCpfixed, log.(β), model_lin(log.(β), fit_Tcharβ_ΔCpfixed[i].param))	
	end
	p_Tcharβ_ΔCpfixed
	# not fully linear
end

# ╔═╡ 79a32ab3-acbd-40c5-803b-0ed47f9ef2c8
begin
	# plot of the residuen of the fit of Tchar over ln(β) for ΔCp fixed
	p_res_Tcharβ_ΔCpfixed = scatter(xlabel="ln(β)", ylabel="resid Tchar in K", legend=false)
	for i=1:length(sub)
		scatter!(p_res_Tcharβ_ΔCpfixed, log.(β),fit_Tcharβ_ΔCpfixed[i].resid, label="$(sub[i].name)")
	end
	p_res_Tcharβ_ΔCpfixed
end
# could be a quadratic part

# ╔═╡ 6af0064a-aaa8-4178-b7f5-bac6b4664d61
scatter(1:97, p2_Tcharβ_ΔCpfixed)

# ╔═╡ b15200e1-1c35-44e4-8b74-5f05646b9ad5
scatter(1:97, p1_Tcharβ_ΔCpfixed)

# ╔═╡ 2eff20d5-ced5-4143-bea6-2c980d3cb4f0
begin
	# with fixed ΔCp
	p_θcharβ_ΔCpfixed = scatter(xlabel="ln(β)", ylabel="θchar in °C", legend=false)
	fit_θcharβ_ΔCpfixed = Array{Any}(undef, length(sub))
	p1_θcharβ_ΔCpfixed = Array{Measurement{Float64}}(undef, length(sub))
	p2_θcharβ_ΔCpfixed = Array{Measurement{Float64}}(undef, length(sub))
	R²_θcharβ_ΔCpfixed = Array{Float64}(undef, length(sub))
	for i=1:length(sub)
		scatter!(p_θcharβ_ΔCpfixed, log.(β), p1_ΔCpfixed[i,:], label="$(sub[i].name)")
		fit_θcharβ_ΔCpfixed[i] = curve_fit(model_lin, log.(β), Measurements.value.(p1_ΔCpfixed[i,:]), [1000.0, -30.0])
		p1_θcharβ_ΔCpfixed[i] = measurement(fit_θcharβ_ΔCpfixed[i].param[1], standard_errors(fit_θcharβ_ΔCpfixed[i])[1])
		p2_θcharβ_ΔCpfixed[i] = measurement(fit_θcharβ_ΔCpfixed[i].param[2], standard_errors(fit_θcharβ_ΔCpfixed[i])[2])
		R²_θcharβ_ΔCpfixed[i] = Rsquare(fit_θcharβ_ΔCpfixed[i], model_lin(log.(β), fit_θcharβ_ΔCpfixed[i].param))
		plot!(p_θcharβ_ΔCpfixed, log.(β), model_lin(log.(β), fit_θcharβ_ΔCpfixed[i].param))	
	end
	p_θcharβ_ΔCpfixed
	# not fully linear
end

# ╔═╡ 9b0c1d9c-35da-4b78-8d90-33740e1cddc5
begin
	# plot of the residuen of the fit of θchar over ln(β) for ΔCp fixed
	p_res_θcharβ_ΔCpfixed = scatter(xlabel="ln(β)", ylabel="resid θchar in °C", legend=false)
	for i=1:length(sub)
		scatter!(p_res_θcharβ_ΔCpfixed, log.(β),fit_θcharβ_ΔCpfixed[i].resid, label="$(sub[i].name)")
	end
	p_res_θcharβ_ΔCpfixed
end
# could be a quadratic part

# ╔═╡ efe67771-c188-435a-8cc6-6b4b55935a0a
scatter(1:97, p2_θcharβ_ΔCpfixed)

# ╔═╡ 8df79670-e400-4ef6-92dd-0fd8e04f4467
scatter(1:97, p1_θcharβ_ΔCpfixed)

# ╔═╡ 6d78afa8-b1d4-4e06-993b-6033777a8405
log(1250)

# ╔═╡ 3d1cf5cc-d5fc-40c0-a1b4-1271709d5a07
begin
	# with fixed ΔCp, quadratic model
	model_quad(x,p) = p[1] .+ p[2].*x .+ p[3].*x.^2
	p_Tcharβ_ΔCpfixed_quad = scatter(xlabel="ln(β)", ylabel="Tchar in K", legend=false)
	fit_Tcharβ_ΔCpfixed_quad = Array{Any}(undef, length(sub))
	p1_Tcharβ_ΔCpfixed_quad = Array{Measurement{Float64}}(undef, length(sub))
	p2_Tcharβ_ΔCpfixed_quad = Array{Measurement{Float64}}(undef, length(sub))
	p3_Tcharβ_ΔCpfixed_quad = Array{Measurement{Float64}}(undef, length(sub))
	R²_Tcharβ_ΔCpfixed_quad = Array{Float64}(undef, length(sub))
	for i=1:length(sub)
		scatter!(p_Tcharβ_ΔCpfixed_quad, log.(β), p2_ΔCpfixed[i,:], label="$(sub[i].name)")
		fit_Tcharβ_ΔCpfixed_quad[i] = curve_fit(model_quad, log.(β), Measurements.value.(p2_ΔCpfixed[i,:]), [1000.0, -30.0, 1.0])
		p1_Tcharβ_ΔCpfixed_quad[i] = measurement(fit_Tcharβ_ΔCpfixed_quad[i].param[1], standard_errors(fit_Tcharβ_ΔCpfixed_quad[i])[1])
		p2_Tcharβ_ΔCpfixed_quad[i] = measurement(fit_Tcharβ_ΔCpfixed_quad[i].param[2], standard_errors(fit_Tcharβ_ΔCpfixed_quad[i])[2])
		p3_Tcharβ_ΔCpfixed_quad[i] = measurement(fit_Tcharβ_ΔCpfixed_quad[i].param[3], standard_errors(fit_Tcharβ_ΔCpfixed_quad[i])[3])
		R²_Tcharβ_ΔCpfixed_quad[i] = Rsquare(fit_Tcharβ_ΔCpfixed_quad[i], model_lin(log.(β), fit_Tcharβ_ΔCpfixed_quad[i].param))
		plot!(p_Tcharβ_ΔCpfixed_quad, log.(β), model_quad(log.(β), fit_Tcharβ_ΔCpfixed_quad[i].param))	
	end
	p_Tcharβ_ΔCpfixed_quad
end

# ╔═╡ 7104e5ee-40fe-4cca-9a05-6833a2d135a8
begin
	# plot of the residuen of the quad fit of Tchar over ln(β) for ΔCp fixed
	p_res_Tcharβ_ΔCpfixed_quad = scatter(xlabel="ln(β)", ylabel="resid Tchar in K", legend=false)
	for i=1:length(sub)
		scatter!(p_res_Tcharβ_ΔCpfixed_quad, log.(β),fit_Tcharβ_ΔCpfixed_quad[i].resid, label="$(sub[i].name)")
	end
	p_res_Tcharβ_ΔCpfixed_quad
end

# ╔═╡ 49946569-5fe4-4dc3-8e8e-dd37f5451fc6
scatter(1:97, p1_Tcharβ_ΔCpfixed_quad)

# ╔═╡ 63282744-d7fe-462a-84dc-69a1d93e3fff
scatter(1:97, p2_Tcharβ_ΔCpfixed_quad)

# ╔═╡ 1ec69ed5-bf13-4eec-b5fd-09c4a1664b8c
scatter(1:97, p3_Tcharβ_ΔCpfixed_quad)

# ╔═╡ 6ff1451b-e3af-4199-92bd-4bcec0255556
begin
	# with fixed ΔCp, quadratic model
	#model_quad(x,p) = p[1] .+ p[2].*x .+ p[3].*x.^2
	p_θcharβ_ΔCpfixed_quad = scatter(xlabel="ln(β)", ylabel="θchar in °C", legend=false)
	fit_θcharβ_ΔCpfixed_quad = Array{Any}(undef, length(sub))
	p1_θcharβ_ΔCpfixed_quad = Array{Measurement{Float64}}(undef, length(sub))
	p2_θcharβ_ΔCpfixed_quad = Array{Measurement{Float64}}(undef, length(sub))
	p3_θcharβ_ΔCpfixed_quad = Array{Measurement{Float64}}(undef, length(sub))
	R²_θcharβ_ΔCpfixed_quad = Array{Float64}(undef, length(sub))
	for i=1:length(sub)
		scatter!(p_θcharβ_ΔCpfixed_quad, log.(β), p1_ΔCpfixed[i,:], label="$(sub[i].name)")
		fit_θcharβ_ΔCpfixed_quad[i] = curve_fit(model_quad, log.(β), Measurements.value.(p1_ΔCpfixed[i,:]), [1000.0, -30.0, 1.0])
		p1_θcharβ_ΔCpfixed_quad[i] = measurement(fit_θcharβ_ΔCpfixed_quad[i].param[1], standard_errors(fit_θcharβ_ΔCpfixed_quad[i])[1])
		p2_θcharβ_ΔCpfixed_quad[i] = measurement(fit_θcharβ_ΔCpfixed_quad[i].param[2], standard_errors(fit_θcharβ_ΔCpfixed_quad[i])[2])
		p3_θcharβ_ΔCpfixed_quad[i] = measurement(fit_θcharβ_ΔCpfixed_quad[i].param[3], standard_errors(fit_θcharβ_ΔCpfixed_quad[i])[3])
		R²_θcharβ_ΔCpfixed_quad[i] = Rsquare(fit_θcharβ_ΔCpfixed_quad[i], model_lin(log.(β), fit_θcharβ_ΔCpfixed_quad[i].param))
		plot!(p_θcharβ_ΔCpfixed_quad, log.(β), model_quad(log.(β), fit_θcharβ_ΔCpfixed_quad[i].param))	
	end
	p_θcharβ_ΔCpfixed_quad
end

# ╔═╡ 23414600-154f-4394-b53c-14e0c1155225
sub[77].name

# ╔═╡ 17cdf383-a005-45b3-8388-f0cc2ba66014
begin
	# plot of the residuen of the quad fit of Tchar over ln(β) for ΔCp fixed
	p_res_θcharβ_ΔCpfixed_quad = scatter(xlabel="ln(β)", ylabel="resid Tchar in K", legend=false)
	for i=1:length(sub)
		scatter!(p_res_θcharβ_ΔCpfixed_quad, log.(β),fit_θcharβ_ΔCpfixed_quad[i].resid, label="$(sub[i].name)")
	end
	p_res_θcharβ_ΔCpfixed_quad
end

# ╔═╡ 5b17bc8a-220b-4c20-953e-4c56e6f4d67a
scatter(1:97, p1_θcharβ_ΔCpfixed_quad)

# ╔═╡ 8bf7dc39-783b-4b6a-ac64-72cb1f61e293
scatter(1:97, p2_θcharβ_ΔCpfixed_quad)

# ╔═╡ ed405c97-d0f4-4117-a22b-9c312b49ce3b
scatter(1:97, p3_θcharβ_ΔCpfixed_quad)

# ╔═╡ b09c33e0-7321-4aa6-8e3a-ec14e51db8c6
# plot parameters above over Tchar0 (Tchar at β=250, were the measurements were made

# ╔═╡ bd6b244e-3eeb-42e1-a6df-94a59626cc8c
md"""
#### Calculation of $T_{char}$

The parameters estimated with a column with a phase ratio $\beta_0$ are $T_{char,0}$, $\theta_{char,0}$ and $\Delta C_{p,0}$.

The retention factor for this phase ratio $\beta_0$ is given by:

$\ln k_0 = \left( \frac{\Delta C_{p,0}}{R} +\frac{T_{char,0}}{\theta_{char,0}} \right) \left( \frac{T_{char,0}}{T}-1 \right) + \frac{\Delta C_{p,0}}{R} \ln{\left( \frac{T}{T_{char,0}} \right)}$

The retention factor for another phase ratio $\beta$ can be calculated from $k_0$:

$k = \frac{\beta_0}{\beta}k_0$

or as 

$\ln k = \ln \beta_0 - \ln \beta + \ln k_0$

resp. as

$\ln k = \Delta(\ln\beta) + \ln k_0$ 

For a fixed temperature the $\ln k$ value shifts by the difference of the logarithm of the different phase ratios, if the phase ratio is changed.

For $\ln(k_0)$ the parameter $T_{char,0}$ is defined as the temperature, were the value of $\ln k_0=0$ (or $k_0=1$). For $\ln k$, at another phase ratio $\beta$, the parameter $T_{char,0}$ is no longer the temperature, where $\ln k=0$. The new parameter $T_{char}$ changes. It is shiftet to lower values for increasing phase ratios and increases for decreasing phase ratios.

It should be possible to calculate the new parameter $T_{char,1}$ from the old parameters and the change in phase ratio $\Delta(\ln \beta)$.

$0=\Delta (\ln \beta) +  \left( \frac{\Delta C_{p,0}}{R} +\frac{T_{char,0}}{\theta_{char,0}} \right) \left( \frac{T_{char,0}}{T}-1 \right) + \frac{\Delta C_{p,0}}{R} \ln{\left( \frac{T}{T_{char,0}} \right)}$

$0=\left( \frac{\Delta C_{p,1}}{R} +\frac{T_{char,1}}{\theta_{char,1}} \right) \left( \frac{T_{char,1}}{T}-1 \right) + \frac{\Delta C_{p,1}}{R} \ln{\left( \frac{T}{T_{char,1}} \right)}$

The parameter $T_{char,1}$ is the solution to the first equation (of the last two). WolframAlpha gives

$T(\ln k =0) = -\frac{T_{char,0}(a+T_{char,0})}{a}\frac{1}{W(y)}$

with $W$ the Lamber W function, $a=\theta_{char,0}\Delta C_{p,0}/R$ and 

$y = -\frac{a+T_{char,0}}{a}\exp{\left(\frac{\Delta(\ln \beta)R}{\Delta C_{p,0}} - \frac{T_{char,0}}{a} -1\right)}$

For the case of $\Delta (\ln \beta)=0$ $y$ becomes:

$y = -(\frac{T_{char,0}}{a}+1) \exp{-(\frac{T_{char,0}}{a}+1)}$

The Lamber W function with the argument of the form $u \exp{u}$ has $u$ as solution (this is the definition of the Lambert W function):

$W(u \exp{u}) = u$

Therfore $T(\ln k =0) = T_{char,0}$, which is the definition of $T_{char,0}$
"""

# ╔═╡ 9624d996-07e1-4f00-a3c3-21acccc84e18
function Tchar_β(β, Tchar0, θchar0, ΔCp0, β0)
	a = θchar0*ΔCp0/R
	b = Tchar0/a+1
	Δlnβ = log(β0)-log(β)
	y = -b*exp(Δlnβ*R/ΔCp0 - b)
	Tchar = -Tchar0*b/lambertw(y,-1)
	return Tchar
end

# ╔═╡ 415a4da1-b103-420f-9c42-20841cee5b4b
begin
	beta = 50:2000
	Tchar_0 = sub[ii].Tchar
	θchar_0 = sub[ii].θchar
	ΔCp_0 = sub[ii].ΔCp
	β_0 = 1/(4*sub[ii].φ₀)
	Tchar_beta = Array{Float64}(undef, length(beta))
	for i=1:length(beta)
		Tchar_beta[i] = Tchar_β(beta[i], Tchar_0, θchar_0, ΔCp_0, β_0)
	end
	plot(log.(beta), Tchar_beta, xlabel="ln(β)", ylabel="Tchar in K", label="$(sub[ii].name)")
	#plot(beta, Tchar_beta, xlabel="β", ylabel="Tchar in K", label="$(sub[ii].name)")
end

# ╔═╡ baab3335-5e2a-40fb-be02-2047db2877bf
function θchar_β(β, Tchar0, θchar0, ΔCp0, β0)
	Tchar = Tchar_β(β, Tchar0, θchar0, ΔCp0, β0)
	Δlnβ = log(β0)-log(β)
	θchar = Tchar/((Δlnβ-ΔCp0/R*log(Tchar0/Tchar))/(Tchar/Tchar0-1)-ΔCp0/R)
	return θchar
	# this calculation is wrong for β=β0!!!!
	# otherwise the values match well with the fitted values
end

# ╔═╡ 55a22483-2f58-4703-bd1d-a387a79a6eee
begin
	#beta = 50:2000
	#Tchar_0 = sub[ii].Tchar
	#θchar_0 = sub[ii].θchar
	#ΔCp_0 = sub[ii].ΔCp
	#β_0 = 1/(4*sub[ii].φ₀)
	θchar_beta = Array{Float64}(undef, length(beta))
	for i=1:length(beta)
		θchar_beta[i] = θchar_β(beta[i], Tchar_0, θchar_0, ΔCp_0, β_0)
	end
	plot(log.(beta), θchar_beta, xlabel="ln(β)", ylabel="θchar in °C", label="$(sub[ii].name)")
	#plot(beta, Tchar_beta, xlabel="β", ylabel="Tchar in K", label="$(sub[ii].name)")
end

# ╔═╡ ebe35aab-f020-4633-9ca4-6a218ad390bf
p1_ΔCpfixed[ii,:]

# ╔═╡ 44017018-a4fd-4fed-86f0-70fae565552c
β[5], θchar_β(β[5], Tchar_0, θchar_0, ΔCp_0, β_0)

# ╔═╡ a5c4e731-1cd4-45b6-8815-4239629a33e5
y

# ╔═╡ f0151f17-e908-47f1-a1a0-2fe4be40a0b9
minimum(y), maximum(y)

# ╔═╡ 7ff1054c-1a29-40d8-98d1-e065a734e059
W

# ╔═╡ 52c5ec8b-f6d5-493b-a2f6-3bceddb9709a
minimum(W), maximum(W)

# ╔═╡ 05d75b66-3686-4934-a1a6-059d096cf4da
T0[ii,:]

# ╔═╡ 73c85399-d905-416f-80fc-535efe4b22e7
Tchar_β(125, sub[ii].Tchar, sub[ii].θchar, sub[ii].ΔCp, 250.0)

# ╔═╡ cdbfc7c5-d81b-45b8-b56c-7f8f658e74c7
sub[ii].Tchar

# ╔═╡ 9fd3b35d-8e09-4271-bbb5-892f0e07b665
minimum(T0), maximum(T0)

# ╔═╡ 635c2345-0d8b-4ee7-83e1-e269e1dd7669
begin
	Tchar0 = Array{Float64}(undef, length(sub))
	for i=1:length(sub)
		Tchar0[i] = sub[i].Tchar
	end
	Tchar0
end

# ╔═╡ aa5cc982-b2f3-4749-9b8d-70be50714101
T0.-Tchar0

# ╔═╡ b9ba73ad-75f8-4037-814f-46a907fb7dfa
T0.-p2_ΔCpfixed #p2_ΔCpfixed is Tchar from the fits for different β -> nearly zero 

# ╔═╡ fb1b0a8f-c041-44d6-96c6-6f32866d387b
p1_ΔCpfixed

# ╔═╡ Cell order:
# ╠═e3da9b92-eb94-11eb-0c86-0f1b7c98a586
# ╟─c7679cc4-06c4-41af-8695-47436330b06f
# ╠═0875234e-ab67-47ad-8a86-b9121ba5e7d4
# ╠═b43aa25c-a502-4eda-8ef3-4f7b8b87e562
# ╠═ec472f66-78ee-4cf5-a368-fe083cd67c06
# ╠═8731c46e-594a-40bb-a262-ee1d0e87f1a8
# ╠═e2684d83-3251-4245-bc74-e6107888b9e3
# ╠═df388e08-face-436b-8b82-0ee72e83d845
# ╠═62c16f38-46bf-47e4-b602-60718d47fcaf
# ╠═397e2565-5ea4-4c97-8d43-684d0f2dee70
# ╠═cc89b0a9-38eb-4cb9-8d66-7281258905f8
# ╠═914876b4-e3ec-4558-95fb-fabd87db0a41
# ╠═dea5b4d0-78f9-467a-905c-78f9795c4b86
# ╠═04ccde67-0171-4069-9c1f-eb88bd0d2e33
# ╠═1c38f24c-657f-4860-8d6f-4e3b3e15d29e
# ╠═28e42e86-660d-4cac-8121-71adf79a7c83
# ╠═afad4f9b-c0b6-4dd7-bf66-0d23deba7e87
# ╠═1af6888d-200d-40da-bece-2ba53ba21aae
# ╠═45dfa9f0-5282-4a59-9fd0-28521a6976c2
# ╠═93c4de7c-b25c-4d52-9a48-33f9b6a10145
# ╠═cfc400d2-ee2b-4a23-b8c8-b19bd2d74b59
# ╠═7301cc6c-b38a-4cee-b566-db360f132b8c
# ╠═228b97c6-4230-459d-bc19-567045201ec6
# ╠═79a32ab3-acbd-40c5-803b-0ed47f9ef2c8
# ╠═6af0064a-aaa8-4178-b7f5-bac6b4664d61
# ╠═b15200e1-1c35-44e4-8b74-5f05646b9ad5
# ╠═2eff20d5-ced5-4143-bea6-2c980d3cb4f0
# ╠═9b0c1d9c-35da-4b78-8d90-33740e1cddc5
# ╠═efe67771-c188-435a-8cc6-6b4b55935a0a
# ╠═8df79670-e400-4ef6-92dd-0fd8e04f4467
# ╠═6d78afa8-b1d4-4e06-993b-6033777a8405
# ╠═3d1cf5cc-d5fc-40c0-a1b4-1271709d5a07
# ╠═7104e5ee-40fe-4cca-9a05-6833a2d135a8
# ╠═49946569-5fe4-4dc3-8e8e-dd37f5451fc6
# ╠═63282744-d7fe-462a-84dc-69a1d93e3fff
# ╠═1ec69ed5-bf13-4eec-b5fd-09c4a1664b8c
# ╠═6ff1451b-e3af-4199-92bd-4bcec0255556
# ╠═23414600-154f-4394-b53c-14e0c1155225
# ╠═17cdf383-a005-45b3-8388-f0cc2ba66014
# ╠═5b17bc8a-220b-4c20-953e-4c56e6f4d67a
# ╠═8bf7dc39-783b-4b6a-ac64-72cb1f61e293
# ╠═ed405c97-d0f4-4117-a22b-9c312b49ce3b
# ╠═b09c33e0-7321-4aa6-8e3a-ec14e51db8c6
# ╠═bd6b244e-3eeb-42e1-a6df-94a59626cc8c
# ╠═8917fa51-228e-4b43-b3c2-f83a4aa79ce5
# ╠═9624d996-07e1-4f00-a3c3-21acccc84e18
# ╠═415a4da1-b103-420f-9c42-20841cee5b4b
# ╠═baab3335-5e2a-40fb-be02-2047db2877bf
# ╠═55a22483-2f58-4703-bd1d-a387a79a6eee
# ╠═ebe35aab-f020-4633-9ca4-6a218ad390bf
# ╠═44017018-a4fd-4fed-86f0-70fae565552c
# ╠═a5c4e731-1cd4-45b6-8815-4239629a33e5
# ╠═f0151f17-e908-47f1-a1a0-2fe4be40a0b9
# ╠═7ff1054c-1a29-40d8-98d1-e065a734e059
# ╠═52c5ec8b-f6d5-493b-a2f6-3bceddb9709a
# ╠═05d75b66-3686-4934-a1a6-059d096cf4da
# ╠═73c85399-d905-416f-80fc-535efe4b22e7
# ╠═cdbfc7c5-d81b-45b8-b56c-7f8f658e74c7
# ╠═9fd3b35d-8e09-4271-bbb5-892f0e07b665
# ╠═635c2345-0d8b-4ee7-83e1-e269e1dd7669
# ╠═aa5cc982-b2f3-4749-9b8d-70be50714101
# ╠═b9ba73ad-75f8-4037-814f-46a907fb7dfa
# ╠═fb1b0a8f-c041-44d6-96c6-6f32866d387b
