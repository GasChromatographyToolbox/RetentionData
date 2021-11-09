### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 8c8f61ce-4134-11ec-166f-f122d3fd3e27
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	using DrWatson, DataFrames, CSV, Plots, LsqFit, Measurements
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 2fc3ea2c-da5b-478e-9042-c5a92932262b
md"""
# Survey of the data of thermodynamic data
"""

# ╔═╡ 95fd0d58-8a21-4117-9844-f25f9e615ae9
begin
	db = DataFrame(CSV.File(datadir("Databases", "Database_append.csv"), header=1, silencewarnings=true))
	select!(db, Not([:Cnumber, :Hnumber, :Onumber, :Nnumber, :Ringnumber]))
	sp = unique(db.Phase)
	db
end

# ╔═╡ c6ca82a6-fdda-47ad-93a3-a22926c94b10
sp

# ╔═╡ cbc4058b-8da4-4d84-9dcc-c9cacd4c4737
begin
	df_sp = Array{DataFrame}(undef, length(sp))
	for i=1:length(sp)
		df_sp[i] = filter(:Phase => Phase -> Phase==sp[i], db)
	end
	md"""
	Filter for the different stationary phases and extract Tchar, θchar and ΔCp
	
	$(embed_display(df_sp))
	"""
end

# ╔═╡ c34ffee2-b6e0-4e35-8278-cf906d3be293
md"""
## Exclusions
"""

# ╔═╡ 50cf085a-a6b4-484b-9d68-4e9763e51969
begin
	# exclude Naphthalin from ZB-PAH-CT data
	ii = findfirst(sp.=="ZB-PAH-CT")
	jj = findfirst(df_sp[8].Name.=="Naphthalin")
	delete!(df_sp[ii], jj);
	md"""
	Exclude Naphthalin from ZB-PAH-CT data.
	"""
end

# ╔═╡ f72fa197-ea0d-4bcc-99b2-6292ff4fd034
begin
	# exclude data for IL111
	iii = findfirst(sp.=="IL111")
	deleteat!(df_sp, iii)
	md"""
	Exclude the data for stationary phase IL111.
		"""
end

# ╔═╡ 9fc410e8-a9e7-4382-992c-23ed65542ab1
begin
	alkanes = String[]
	for i=1:100
		push!(alkanes, string("C",i))
	end
	PAH = ["Naphthalin", "Acenaphthylene", "Acenaphthene", "Fluorene", "Phenanthrene", "Anthracene", "Fluoranthene", "Pyrene", "Benz[a]anthracene", "Chrysene", "Benzo[b]fluoranthene", "Benzo[k]fluoranthene", "Benzo[a]pyrene", "Dibenzo[a,h]anthracene", "Indeno[1,2,3-cd]pyrene", "Benzo[ghi]perylene"]
	PAH_CAS = ["91-20-3", "208-96-8", "83-32-9", "86-73-7", "85-01-8", "120-12-7", "206-44-0", "129-00-0", "56-55-3", "218-01-9", "205-99-2", "207-08-9", "50-32-8", "53-70-3", "193-39-5", "191-24-2"]
	md"""
	## Definition of solute classes.
	
	- Alkanes: C1 to C100
	- PAH
	"""
end

# ╔═╡ 7594a514-57f8-40fe-bb92-fa8e6853f880
begin
	df_sp_alkanes = Array{DataFrame}(undef, length(df_sp))
	df_sp_PAH = Array{DataFrame}(undef, length(df_sp))
	for i=1:length(df_sp)
		df_sp_alkanes[i] = filter(:Name => Name -> in.(Name, Ref(alkanes)), df_sp[i])
		df_sp_PAH[i] = filter(:CAS => CAS -> in.(CAS, Ref(PAH_CAS)), df_sp[i])
	end
	md"""
	Filter the data for defined solute classes.
	"""
end

# ╔═╡ 1993f2bf-a421-45ea-b4a1-532cd7c44f50
md"""
## Fit function ``\theta_{char} = p_1 \left(\frac{T_{char}}{T_{st}}\right)^{p_2}``. 
"""

# ╔═╡ 5c9fe062-e4bf-48be-b9ea-4a34233b5fd1
function fit_θchar_over_Tchar(Tchar, θchar)
	# fit a model to the plot (derived from Blumberg.2010)
	model(x, p) = p[1].*(x./273.15).^p[2]
	p0 = [22, 0.7]
	fit = curve_fit(model, Tchar, θchar, p0)
	p1 = measurement(fit.param[1], standard_errors(fit)[1])
	p2 = measurement(fit.param[2], standard_errors(fit)[2])
	return p1, p2, fit
end

# ╔═╡ e5023fea-8e1a-44ec-b314-26b5bf3dc990
begin
	# thermodynamic parameters of all data set (without exclusions)
	Tchar_all = Float64[]
	θchar_all = Float64[]
	ΔCp_all = Float64[]
	for i=1:length(df_sp)
		for j=1:size(df_sp[i])[1]
			push!(Tchar_all, df_sp[i].Tchar[j])
			push!(θchar_all, df_sp[i].thetachar[j])
			push!(ΔCp_all, df_sp[i].DeltaCp[j])
		end
	end
end

# ╔═╡ b2f87308-b888-42dc-8555-d16326315704
fit_θchar_over_Tchar(Tchar_all, θchar_all)

# ╔═╡ 2b1d497d-f681-4c9d-b5a7-36f1273d31e9
begin
	p_θchar_Tchar = scatter(df_sp[1].Tchar, df_sp[1].thetachar, xlabel="Tchar in °C", ylabel="θchar in °C", label=sp[1], legend=:topleft)
	
	p_θchar_Tchar_alkanes = scatter(df_sp_alkanes[1].Tchar, df_sp_alkanes[1].thetachar, xlabel="Tchar in °C", ylabel="θchar in °C", label=sp[1], legend=:topleft, title="n-alkanes")

	p_θchar_Tchar_PAH = scatter(df_sp_PAH[1].Tchar, df_sp_PAH[1].thetachar, xlabel="Tchar in °C", ylabel="θchar in °C", label=sp[1], legend=:topleft, title="PAH")
	
	for i=2:length(df_sp)
		scatter!(p_θchar_Tchar, df_sp[i].Tchar, df_sp[i].thetachar, label=sp[i])
		
		scatter!(p_θchar_Tchar_alkanes, df_sp_alkanes[i].Tchar, df_sp_alkanes[i].thetachar, label=sp[i])

		scatter!(p_θchar_Tchar_PAH, df_sp_PAH[i].Tchar, df_sp_PAH[i].thetachar, label=sp[i])
	end

	model(x, p) = p[1].*(x./273.15).^p[2]
	p1_all, p2_all, fit_all = fit_θchar_over_Tchar(Tchar_all.+273.15, θchar_all)
	plot!(p_θchar_Tchar, sort(Tchar_all), sort(model(Tchar_all.+273.15, [Measurements.value(p1_all), Measurements.value(p2_all)])), label="fit all")
	
	md"""
	## Plot of ``\theta_{char}`` over ``T_{char}``

	$(embed_display(p_θchar_Tchar))

	$(embed_display(p_θchar_Tchar_alkanes))

	$(embed_display(p_θchar_Tchar_PAH))
	"""
end

# ╔═╡ 71de179c-efb1-47b1-b1ef-29ceb1485908
p1_all, p2_all

# ╔═╡ 445f7e73-dd48-4e19-a1ea-f3edc42ec85a
begin
	p_ΔCp_Tchar = scatter(df_sp[1].Tchar, df_sp[1].DeltaCp, xlabel="Tchar in °C", ylabel="ΔCp in J mol⁻¹ K⁻¹", label=sp[1], legend=:topleft)

	p_ΔCp_Tchar_alkanes = scatter(df_sp_alkanes[1].Tchar, df_sp_alkanes[1].DeltaCp, xlabel="Tchar in °C", ylabel="ΔCp in J mol⁻¹ K⁻¹", label=sp[1], legend=:topleft, title="n-alkanes")

	p_ΔCp_Tchar_PAH = scatter(df_sp_PAH[1].Tchar, df_sp_PAH[1].DeltaCp, xlabel="Tchar in °C", ylabel="ΔCp in J mol⁻¹ K⁻¹", label=sp[1], legend=:topleft, title="PAH")
	
	for i=2:length(df_sp)
		scatter!(p_ΔCp_Tchar, df_sp[i].Tchar, df_sp[i].DeltaCp, label=sp[i])
	
		scatter!(p_ΔCp_Tchar_alkanes, df_sp_alkanes[i].Tchar, df_sp_alkanes[i].DeltaCp, label=sp[i])
	
		scatter!(p_ΔCp_Tchar_PAH, df_sp_PAH[i].Tchar, df_sp_PAH[i].DeltaCp, label=sp[i])
	end
	
	md"""
	## Plot of ``\Delta C_p`` over ``T_{char}``

	$(embed_display(p_ΔCp_Tchar))

	$(embed_display(p_ΔCp_Tchar_alkanes))

	$(embed_display(p_ΔCp_Tchar_PAH))
	"""
end

# ╔═╡ e0f079c3-5348-4ec4-a3ce-25b8a8e2a007
fit_all

# ╔═╡ 3ef6fdfb-4e4c-4410-bcd8-f4986e55579e
model(Tchar_all, [p1_all, p2_all])

# ╔═╡ Cell order:
# ╠═8c8f61ce-4134-11ec-166f-f122d3fd3e27
# ╟─2fc3ea2c-da5b-478e-9042-c5a92932262b
# ╠═95fd0d58-8a21-4117-9844-f25f9e615ae9
# ╠═c6ca82a6-fdda-47ad-93a3-a22926c94b10
# ╠═cbc4058b-8da4-4d84-9dcc-c9cacd4c4737
# ╠═c34ffee2-b6e0-4e35-8278-cf906d3be293
# ╠═50cf085a-a6b4-484b-9d68-4e9763e51969
# ╠═f72fa197-ea0d-4bcc-99b2-6292ff4fd034
# ╠═9fc410e8-a9e7-4382-992c-23ed65542ab1
# ╠═7594a514-57f8-40fe-bb92-fa8e6853f880
# ╠═1993f2bf-a421-45ea-b4a1-532cd7c44f50
# ╠═5c9fe062-e4bf-48be-b9ea-4a34233b5fd1
# ╠═e5023fea-8e1a-44ec-b314-26b5bf3dc990
# ╠═b2f87308-b888-42dc-8555-d16326315704
# ╠═2b1d497d-f681-4c9d-b5a7-36f1273d31e9
# ╠═71de179c-efb1-47b1-b1ef-29ceb1485908
# ╠═445f7e73-dd48-4e19-a1ea-f3edc42ec85a
# ╠═e0f079c3-5348-4ec4-a3ce-25b8a8e2a007
# ╠═3ef6fdfb-4e4c-4410-bcd8-f4986e55579e
