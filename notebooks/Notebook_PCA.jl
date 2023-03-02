### A Pluto.jl notebook ###
# v0.19.14

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

# ╔═╡ 313457f3-85fd-4274-9b96-fceb277eaae4
using LaTeXStrings

# ╔═╡ 391282b0-f684-4e3d-b4ba-1eba0689fe5e
using MultivariateStats, RDatasets

# ╔═╡ c4ea2bc1-ec27-47b4-9a4b-0c03da27042b
md"""# Supplemental Materials - Investigations and PCA Analysises"""

# ╔═╡ 5861512e-a14b-11ec-3c6b-9bd2953bf909
md"""
# Load `AllParam` data
- Load all files with `AllParam`  from the folder `db_path` (and its subfolders).
- use ChemicalIdentifiers.jl
- combine all data into one DataFrame/csv-file
"""

# ╔═╡ 979d0e68-135a-4570-a835-8be8e7095326
md"""
## DataFrame with the structure of database for GasChromatographySimulator.jl v0.3 and lower
""";

# ╔═╡ 8bc510bf-b40a-4f7f-b5e9-05bdca6fe1a3
#CSV.write("../Databases/newformat_$(parset).csv", new_db)

# ╔═╡ b2b6ef1c-ca2d-4a52-a586-b3c2b13b6bac
md"""
## Plot parameters Tchar, θchar, ΔCp
"""

# ╔═╡ 93a5ed9f-64ce-4ee8-bee0-ecd2ff2dcbb8
md"""
## Plot parameters A, B, C
"""

# ╔═╡ 97d4ff70-092c-40db-90a2-6a8994ad920c
md"""
## Plot parameters ΔHref, ΔSref, ΔCp
"""

# ╔═╡ d28a4b20-af06-4613-810b-dafb9a5f0ccd
md"""
### Observations

The parameters of the Thermodynamic-model do span two planes in the parameter-space. One plane consits soly of data from the group of Harynuk (Karolat, McGinitie). The other plane consists of other sources, including ower one measurements. Interestingly, the parameters from Blumberg are in the second plane, while the data from which they are derived (Karolat) are in the first plane.

The data from the sources of Karolat2010, McGinitie2011, McGinitie2012a, McGinitie2014a and McGinitie2014b are given in the Thermodynamic-model (ΔHref, ΔSref, ΔCp).

The data from Blumberg2017 is given in the ABC and Kcentric parameter sets and with enthalpy/entropy change at Tchar as reference temperature. Therefore, the calculation of refereence entropy/entalpy change at the reference temperature of 90°C should be re-evaluated. 

#### ToDo:
- check the calculation of ΔHref and ΔSref
"""

# ╔═╡ c28cca6a-fdf3-4edf-9534-5c4f677c2889
md"""
# Flag the parameter sets
"""

# ╔═╡ 32c1889e-9da1-44e6-9af9-08af1d7f2020
md"""## K-centric Parameters"""

# ╔═╡ 17cb9134-8880-4cb2-86a1-31eb4e1ee5b5
md"""## ABC Parameters"""

# ╔═╡ 0ad0450e-199a-4be1-a6af-03eb745f08d5
md"""## Thermodynamic Parameters"""

# ╔═╡ 917e88c6-cf5d-4d8f-93e5-f50ea2bb2cdc
md"""
# Find duplicates
"""

# ╔═╡ 496b86c5-900e-4786-a254-082b8155c65b
#CSV.write("../Databases/newformat_$(parset).csv", new_db)

# ╔═╡ a319c1af-33d4-443d-9270-5a0219e25c4c
begin
	root = dirname(@__FILE__)
	project = dirname(root)
	db_path = joinpath(project, "Databases")
end

# ╔═╡ fec530e6-d675-4bb9-8b5a-6aa607574a81
begin
	using CSV, DataFrames, LambertW, Plots, LsqFit, Statistics, ChemicalIdentifiers, Measurements
	include(joinpath(project, "src", "RetentionData.jl"))
	using PlutoUI
end

# ╔═╡ ec6342ed-6f00-4c7f-93d2-21a5f39eb4db
TableOfContents()

# ╔═╡ e92ee318-9f88-40f2-81cf-e2d60bf2f45a
data = RetentionData.load_allparameter_data(db_path);

# ╔═╡ 3b9c6610-6839-4012-aa0b-219a347ca52f
data.data[end];

# ╔═╡ d26fc674-ace5-43ef-af1d-855dfc21eba5
begin
	alldata = RetentionData.dataframe_of_all(data)
	# add flags to alldata
	alldata[!, "flag"] = RetentionData.flag(alldata)
	alldata
end

# ╔═╡ 91d21ca2-5657-4f9b-ae66-629fa77236f1
fl, nfl = RetentionData.flagged_data(alldata);

# ╔═╡ 84fc7b22-2362-4b99-8fd9-3a779d706a4b
begin
	plotly()
	pKcentric_nfl_ = plot(xlabel="Tchar", ylabel="Θchar", zlabel="ΔCp")
	source = unique(nfl.Source)
	for i=1:length(source)
		nfl_filter = filter([:Source] => x -> x == source[i], nfl)
		scatter!(pKcentric_nfl_, nfl_filter.Tchar, nfl_filter.thetachar, nfl_filter.DeltaCp, label=source[i])
	end
	pKcentric_nfl_
end	

# ╔═╡ 1e912349-d26d-45ea-9167-21b0b4165060
begin 
	plotly()
	pABC_nfl_ = plot(xlabel="A", ylabel="B", zlabel="C")
	for i=1:length(source)
		nfl_filter = filter([:Source] => x -> x == source[i], nfl)
		scatter!(pABC_nfl_, nfl_filter.A, nfl_filter.B, nfl_filter.C, label=source[i])
	end
	pABC_nfl_
end

# ╔═╡ 7c95815d-11e5-4de5-8d83-a7ef8518751c
begin 
	plotly()
	pTD_nfl_ = plot(xlabel="DeltaHref", ylabel="DeltaSref", zlabel="DeltaCp")
	for i=1:length(source)
		nfl_filter = filter([:Source] => x -> x == source[i], nfl)
		scatter!(pTD_nfl_, nfl_filter.DeltaHref, nfl_filter.DeltaSref, nfl_filter.DeltaCp, label=source[i])
	end
	pTD_nfl_
end

# ╔═╡ 7eb28537-4d9e-4529-b6f7-b0407bce814c
pKcentric_nfl = scatter(nfl.Tchar, nfl.thetachar, nfl.DeltaCp, label="not flagged", xlabel="Tchar", ylabel="θchar", zlabel="ΔCp", camera=(30, 30));

# ╔═╡ d460b5e1-9060-4fe5-bcea-8c63a0ea071e
begin
	plotly()
	pKcentric = scatter(nfl.Tchar, nfl.thetachar, nfl.DeltaCp, label="not flagged", xlabel="Tchar", ylabel="thetachar", zlabel="DeltaCp")
	scatter!(pKcentric, fl.Tchar, fl.thetachar, fl.DeltaCp, label="flagged", c=:red, m=:cross)
	pKcentric
end

# ╔═╡ 17387b6a-e126-4923-ac14-229d980ed9ad
pABC_nfl = scatter(nfl.A, nfl.B, nfl.C, label="not flagged", xlabel="A", ylabel="B", zlabel="C");

# ╔═╡ 0e559fd4-5fef-4375-9e7a-1e052cf61b7a
begin
	plotly()
	pABC = scatter(nfl.A, nfl.B, nfl.C, label="not flagged", xlabel="A", ylabel="B", zlabel="C")
	scatter!(pABC, fl.A, fl.B, fl.C, label="flagged", c=:red, m=:cross)
	pABC
end

# ╔═╡ 8287590e-0e6d-4591-a616-0b64d84ca767
pTD_nfl = scatter(nfl.DeltaHref, nfl.DeltaSref, nfl.DeltaCp, label="not flagged", xlabel="DeltaHref", ylabel="DeltaSref", zlabel="DeltaCp");

# ╔═╡ ef518fe9-1920-4c95-9ef3-2b2542621fc0
begin
	plotly()
	pTD = scatter(nfl.DeltaHref, nfl.DeltaSref, nfl.DeltaCp, label="not flagged", xlabel="DeltaHref", ylabel="DeltaSref", zlabel="DeltaCp")
	scatter!(pTD, fl.DeltaHref, fl.DeltaSref, fl.DeltaCp, label="flagged", c=:red, m=:cross)
	pTD
end

# ╔═╡ fda7f746-0f7e-4dff-b8e6-7765f580e542
md"""
## Some more filters for the duplicates
"""

# ╔═╡ 9290687e-26d1-4ccf-9603-d46b8ae9dcb7
md"""
### Conclusions from duplicates
- if _Blumberg2017_ is one of the sources, keep these
- in the case of the three duplicates from _McGinitie2011_ different mobile phases where used (He, H₂ and N₂) resp. different column diameters (0.1mm, 0.2mm, 0.25mm, 0.32mm, 0.53mm; these are also duplicates with _McGinitie2014a_)-> similar results, keep the ones with He resp. 0.25mm
- duplicates from only _Marquart2020_ represent different observed isomers? related to the substance -> keep them all
- some complete duplicates (parameters have the same value) -> use only the first substance
"""

# ╔═╡ eaed9fdd-d4a2-41f7-8e05-588869e12780
#CSV.write("../Databases/oldformat_$(parset).csv", old_db)

# ╔═╡ 070c0b3b-efa4-4f4d-a567-87ba4b7c936b
unique(alldata.Phase)

# ╔═╡ ef0c1872-da54-48b5-8212-634d7d91e9ac
# filter for all alkanes
function filter_Cat(newdata, cat)
	iCat = findall(occursin.("Cat", names(newdata)))
	i_true = Int[]
	for i=1:length(newdata.CAS)
		for j=1:length(iCat)
			if ismissing(newdata[!, iCat[j]][i]) == false
				if cat == newdata[!, iCat[j]][i]
					push!(i_true, i)
				end
			end
		end
	end
	filtered_data = newdata[i_true,:]
	return filtered_data
end

# ╔═╡ 873b90e2-bc5d-4272-a90b-978bb4e1358d
md"""# Investigate Compound Influence"""

# ╔═╡ d7acb6b7-539a-4d4a-bbff-d066e6535ff6
md"""Can we found systematical influences of the compound category to the parameters?"""

# ╔═╡ d24832da-ff7d-4d4c-af5d-6662eca846cf
begin 
scatter(nfl.Tchar, nfl.thetachar, label="not flagged", xlabel="T_{char}", ylabel="θ_{char}", ylims=(0,80))
end	

# ╔═╡ a94f22a9-5332-429a-bd31-3c1799781ffe
scatter(nfl.Tchar, nfl.DeltaCp, label="not flagged", xlabel="T_{char}", ylabel="ΔC_p", camera=(30, 45))

# ╔═╡ b846be08-4850-4bbd-9366-b7e2207bb370
md"""### Filter functions"""

# ╔═╡ 7b39d781-b51c-4406-af76-ce9459b1cdb6
md"""#### Filter Data per Literature"""

# ╔═╡ c365e82d-d5e5-4be0-81c9-c5a8097d4e8c
function SourceFilter(alldata)
SourceFilter=Array{Any}(undef, size(unique(alldata.Source))[1])
	for i=1:size(unique(alldata.Source))[1]
		try
		SourceFilter[i] =filter([:Source] => x -> string(x) == unique(alldata.Source)[i], alldata)
		catch
			SourceFilter[i] =filter([:Source] => x -> string(x)== string(unique(alldata.Source)[i]), alldata)
		end	
	end	

N_SourceCompound=Array{Any}(undef, size(unique(alldata.Source))[1])
	for i=1:size(unique(alldata.Source))[1]
	N_SourceCompound[i]=size(unique(SourceFilter[i].CAS))[1]	
	end

N_Column=Array{Any}(undef, size(unique(alldata.Source))[1])
	for i=1:size(unique(alldata.Source))[1]
	N_Column[i]=size(unique(SourceFilter[i].Phase))[1]	
	end	

N_Elements=Array{Any}(undef, size(unique(alldata.Source))[1])
	for i=1:size(unique(alldata.Source))[1]
	N_Elements[i]=size(SourceFilter[i].Name)[1]	
	end
SourceFilter	

df=DataFrame(Source=unique(alldata.Source),N_Elements=N_Elements,  NumberofCompounds=N_SourceCompound, N_Column=N_Column)	
	return df, SourceFilter
end	

# ╔═╡ 512ae246-792c-45c7-ae58-406bc70d47ff
SourceFilter(alldata);

# ╔═╡ e945beeb-a9de-4a2d-bc8b-73c3f3dd37b1
SourceFilter(nfl);

# ╔═╡ 20eb935e-c66d-46fe-b309-e4f37f29e667
SourceFilter(nfl)[2][12:14];

# ╔═╡ b8c51d43-fd7d-440e-b516-a49a21cf44b8
#SV.write(string(pwd(), "\\Table1.csv"), append!(append!(SourceFilter(nfl)[2][12], SourceFilter(nfl)[2][13]), SourceFilter(nfl)[2][14]))

# ╔═╡ 015e4eb0-f60e-474f-9ab9-98f6c44b9d8f
CAS_example=
["78-70-6"
"97-53-0"
"4602-84-0"
"4602-84-0"
"106-24-1"
"54464-57-2"
"54464-57-2"
"54464-57-2"
"54464-57-2"
"76-49-3"
"97-54-1"
"102-76-1"
"60-01-5"
"621-71-6"
"621-70-5"
"538-24-9"
"538-23-8"
"37680-73-2"
"35065-28-2"
"35065-27-1"
"35065-29-3"
"7012-37-5"
"35693-99-3"
"120-12-7"
"56-55-3"
"191-24-2"
"218-01-9"
"53-70-3"
"193-39-5"];

# ╔═╡ bbaa8bd7-5b38-4fda-aab9-d97d0c14c410
#CSV.write(string(pwd(), "\\Literature.csv"), SourceFilter(alldata)[1])

# ╔═╡ 44c3be0a-3376-4b51-bd07-683e96029d69
#CSV.write(string(pwd(), "\\verfication.csv"), SourceFilter(nfl)[1])

# ╔═╡ 33ec66d1-2a25-45b6-936c-0acff05e624c
unique(nfl.CAS)

# ╔═╡ 6729a636-aa06-4695-8fd6-2078322e2ffa
md"""#### Filter for Substance category"""

# ╔═╡ 99c115a1-7bbf-4aa1-b526-6200e4d8a622
begin 
	function SubstFilter(nfl)
	SubstanceFilter=Array{Any}(undef, size(unique(nfl.Cat_1))[1])
	for i=1:size(unique(nfl.Cat_1))[1]
		try
		SubstanceFilter[i] =filter([:Cat_1] => x -> string(x) == unique(nfl.Cat_1)[i], nfl)
		catch
			SubstanceFilter[i] =filter([:Cat_1] => x -> string(x)== string(unique(nfl.Cat_1)[i]), nfl)
		end	
	end	
return SubstanceFilter	
	end	
	SubstanceFilter=SubstFilter(nfl)
end	;

# ╔═╡ 1d708017-f3d8-4126-9808-3cbc7bce0050
SubstFilter(fl);

# ╔═╡ 421e55ec-a18b-4aba-9c88-0e53cdcef75d
md""" ## ABC Model """

# ╔═╡ fb2bf9b6-e7e2-4b5f-8013-e6db5d779f85
begin
	plotly()
PlotABC=Array{Any}(undef, size(SubstanceFilter)[1])
	PlotABC=Plots.scatter(xlabel="A", ylabel="B", zlabel="C", legend=true)
	for i=1:size(SubstanceFilter)[1]
		if i==4
		Name="others"
		else
		Name=SubstanceFilter[i].Cat_1[1]	
		end
		if i < 16
		Plots.scatter!(SubstanceFilter[i].A, SubstanceFilter[i].B, SubstanceFilter[i].C, label=string(Name), legend=:outerright, m=Plots.supported_markers()[i], c=i)
		else
		Plots.scatter!(SubstanceFilter[i].A, SubstanceFilter[i].B, SubstanceFilter[i].C, label=string(Name), legend=:outerright, m=Plots.supported_markers()[i.-15], c=i)
		end
		
	end		
PlotABC
end

# ╔═╡ 21277d27-7ea7-4adc-b1b8-4429c29cdf51
md"""## K-centric Retention Parameter"""

# ╔═╡ 9b0d91eb-6672-4cf8-82e4-844becf699bb
function PaperPlot(SubstanceFilter)
plotly() 
	Plot=Array{Any}(undef, size(SubstanceFilter)[1])
	Plot=Plots.scatter(xlabel="Tchar",zlabel="thetachar", ylabel="DeltaCp")
	for i=1:size(SubstanceFilter)[1]
		Plots.scatter!(SubstanceFilter[i].Tchar, SubstanceFilter[i].DeltaCp,SubstanceFilter[i].thetachar, label=string(SubstanceFilter[i].Cat_1[1]), markers=Plots.supported_markers()[i], c=i)
	end		
return Plot
end	

# ╔═╡ ed134625-a178-48d8-856f-e8b25f57842d
begin
plotly()
PlotKcentric3D=Array{Any}(undef, size(SubstanceFilter)[1])
	PlotKcentric3D=Plots.scatter(xlabel="Tchar", ylabel=L"θchar", zlabel="ΔCp")
	for i=1:size(SubstanceFilter)[1]
		if i==4
		Name="others"
		else
		Name=SubstanceFilter[i].Cat_1[1]	
		end
		
		if i < 16
			Plots.scatter!(PlotKcentric3D, SubstanceFilter[i].Tchar, SubstanceFilter[i].thetachar, SubstanceFilter[i].DeltaCp, label=string(Name), legend=:outerright,markers=Plots.supported_markers()[i], c=[i], ylims=(0,90))
		else 
			Plots.scatter!(PlotKcentric3D, SubstanceFilter[i].Tchar, SubstanceFilter[i].thetachar, SubstanceFilter[i].DeltaCp, label=string(Name),markers=Plots.supported_markers()[i.-15], c=[i], ylims=(0,90))
		end	
		#
	end
	PlotKcentric3D
end

# ╔═╡ 4b5de04c-24ef-4f79-937a-f2816d7397d7
begin 
gr()
	PlotTcharCp=Array{Any}(undef, size(SubstanceFilter)[1])
	PlotTcharCp=Plots.scatter(xlabel=L"T_{char}", ylabel=L"ΔC_p")
	for i=1:size(SubstanceFilter)[1]
		Plots.scatter!(SubstanceFilter[i].Tchar, SubstanceFilter[i].DeltaCp, label=string(SubstanceFilter[i].Cat_1[1]), legend=:false, markers=Plots.supported_markers()[i], c=i)
	end		
	 PlotTcharCp
end

# ╔═╡ 3b0e9656-63e8-4fa3-8f06-33f8a2fcd4dc
begin
gr()
PlotTcharTheta=Array{Any}(undef, size(SubstanceFilter)[1])
	PlotTcharTheta=Plots.scatter(xlabel=L"T_{char}", ylabel=L"θ_{char}")
	for i=1:size(SubstanceFilter)[1]
		if i==4
		Name="others"
		else
		Name=SubstanceFilter[i].Cat_1[1]	
		end
		Plots.scatter!(SubstanceFilter[i].Tchar, SubstanceFilter[i].thetachar, label=string(Name), legend=:outerright,markers=Plots.supported_markers()[i], c=[i], ylims=(0,90))
		#
	end
	PlotTcharTheta
end

# ╔═╡ 2d89e6eb-4c6e-4c6f-8a82-55eaf66de976
#begin plotly()
#Plots.scatter!()
#	gui()
#end	

# ╔═╡ 09b949bf-36b5-44a6-ab13-e0074b824756
	begin 
		gr()
		Plots.plot(PlotTcharTheta, PlotTcharCp, layout=(1,2), size=(1200,600), dpi=500, minorgrid=true, silentwarning=true)
		
	end	;

# ╔═╡ 13a8c60e-979b-47e6-bcdc-7a314fd27ebb
#savefig(Plots.plot(PlotTcharTheta, PlotTcharCp, layout=(1,2), size=(1200,600), dpi=500, minorgrid=true), "Kcentric.svg" )

# ╔═╡ b5c0895f-1683-4f46-960b-9463819610a1
md""" ### Dependence between θchar and Tchar

Fit Blumberg2022/2010 Model 

$θ_{char}=(T_{char})^{0.7} \cdot (10^3 \cdot𝝋)^{0.09}$


"""

# ╔═╡ a64f174b-fac7-4566-b8b1-2c603dbd3dd1
md""" Sort data with same φ"""

# ╔═╡ 492fabbd-c84a-421e-a169-4229cf444e11
begin 
	PhiFilter=Array{Any}(undef, size(unique(nfl.beta0))[1])
	for i=1:size(unique(nfl.beta0))[1]
		try
		PhiFilter[i] =filter([:beta0] => x -> x == unique(nfl.beta0)[i], nfl)
		catch
			PhiFilter[i] =filter([:beta0] => x -> x== unique(nfl.beta0)[i], nfl)
		end	
	end	
PhiFilter	
end	

# ╔═╡ 80745a29-6a58-47c8-a176-e6f58a858d82
begin 
	Theta_calc=Array{Any}(undef, size(unique(nfl.beta0)))
	for i=1:size(unique(nfl.beta0))[1]
		
				data=Array{Any}(undef,size(PhiFilter[i].Tchar)[1])	
				for j=1:size(PhiFilter[i].Tchar)[1]		
				try
					data[j] =(PhiFilter[i].Tchar[j]).^0.7.*(1000 .*1 ./(4 .*PhiFilter[i].beta0[j])).^0.09
				catch
					data[j]=NaN
				end
				end
		Theta_calc[i]=data
	end	
	Theta_calc
end	

# ╔═╡ 333ffd2d-d4db-46fd-b3f3-5cbffeb5243d
begin
Plotti=Array{Any}(undef, size(unique(nfl.beta0)))
			Plots.scatter()
	for i=1:size(unique(nfl.beta0))[1]

		Plotti=Plots.scatter!(PhiFilter[i].Tchar, PhiFilter[i].thetachar, label=unique(nfl.beta0)[i], xlabel=L"T_{char}", ylabel=L"θ_{char}", title="Influence of Phase ratio")

	end
	Plotti
end	;

# ╔═╡ fc4a5a98-f1ab-45b5-9c2d-44693ad35313
begin
Theta_plot=Array{Any}(undef, size(unique(nfl.beta0))[1])
	for j=1:size(unique(nfl.beta0))[1]
	x=0:1:650
	y=((x.+273.15)./273.15).^0.7.*(22 .*(1000 .*1 ./(4 .* unique(nfl.beta0)[j])).^0.09)
		Theta_plot[j]=Plots.plot!(Plotti, x, y, label=unique(nfl.beta0)[j], legend=:outerright, c=j, ylims=(0,80), xlabel="Tchar", ylabel="θchar")
	end
	Theta_plot[1]
end	

# ╔═╡ bf9adfcc-5baf-486d-b606-4f1a8d935a17
unique(nfl.beta0)[5]

# ╔═╡ 156b9224-3714-4039-aa49-ecc6b04d38dc
Plots.scatter(SubstanceFilter[4].A, SubstanceFilter[4].B, SubstanceFilter[4].C, label=unique(nfl.Cat_1)[4]);

# ╔═╡ 202ab200-b358-4b3b-8bd9-da3a1c158d39
#CSV.write("$(pwd())\\Database.csv", filter([:Source]=>(x)-> occursin.(string(x), string("Marquart2020","Duong2022", "Brehmer2022")), nfl))

# ╔═╡ f687a423-bf8f-4c5f-b98e-f8afc5c4ee9e
md""" # PCA Analysis"""

# ╔═╡ ad5312d2-43bb-42f6-8327-a3de2dfc1aae
md"""## PCA ABC Model"""

# ╔═╡ 4ec1f88e-d2e4-43c0-9e43-68322ab96e1c
nfl

# ╔═╡ 534ec747-4b14-4c01-a837-7e6db2479381
ABC_Training=Matrix(nfl[1:2:end-1,5:7])'

# ╔═╡ 982446df-5021-44ca-b7a8-2f0042baf88d
 ABC_Testing=Matrix(nfl[2:2:end,5:7])'

# ╔═╡ b6fcf3bb-9b33-42ba-97fa-f0116ae0ec46
PCA_ABC=fit(PCA, ABC_Training; maxoutdim=2)

# ╔═╡ 73611573-ae96-49d5-86be-8b646cd5c8bc
ABC_Y=predict(PCA_ABC,  ABC_Testing)

# ╔═╡ 09806af8-d0d8-49a9-a2fe-4357992584b8
reconstruct(PCA_ABC, ABC_Y)

# ╔═╡ 8a4541fc-7ebb-46f0-a652-10f7ba1ea689
md"""## PCA K-centric Model"""

# ╔═╡ 52f944bf-2038-497a-ae39-a2ee95ad67d2
kcentric_Training=Matrix(nfl[1:2:end-1,8:10])'

# ╔═╡ 32471c63-a5c0-4ad9-b9e4-1701a8a9902b
kcentric_Testing=Matrix(nfl[2:2:end,8:10])'

# ╔═╡ fb7017ab-b8bb-4175-8434-5f9c1e026509
PCA_kcentric=fit(PCA, kcentric_Training; maxoutdim=2)

# ╔═╡ 0698e398-d204-448e-9460-5eb2a51a5d05
kcentric_Y= predict(PCA_kcentric, kcentric_Testing)

# ╔═╡ ff6f166d-6d18-44c1-8be3-8d6d40bcd738
reconstruct(PCA_kcentric, kcentric_Y)

# ╔═╡ 4badc718-b4a0-4642-adcb-b88a55193280
SubstanceFilter[4].Cat_1[1]

# ╔═╡ 00703bb9-21f6-4f93-95f4-45995d81ee5d
begin 
	gr()
	KcentricPCA_Plot=Array{Any}(undef, size(SubstanceFilter)[1])
		KcentricPCA_Plot=Plots.scatter(xlabel="PC1", ylabel="PC2")
	for i=1:size(SubstanceFilter)[1]
				if i==4
		Name="others"
		else
		Name=SubstanceFilter[i].Cat_1[1]	
		end	
		Predict=predict(PCA_kcentric,Matrix(SubstanceFilter[i][!,8:10])')
		Plots.scatter!(KcentricPCA_Plot, Predict[1,:],Predict[2,:], label=Name, legend=:outerright, markers=Plots.supported_markers()[i], title="Compounds",dpi=500, size=(800,400), minorgrid=true)
	end		
KcentricPCA_Plot
end

# ╔═╡ 8f2f176e-cab0-4ab8-8d9c-f635d6e95460
#savefig(KcentricPCA_Plot, "PCA.svg" )

# ╔═╡ c00ded9a-6d20-4748-9dc6-789831820427
begin 
plotly()
	KcentricPCA_Plot_beta=Array{Any}(undef, size(PhiFilter)[1])
		KcentricPCA_Plot_beta=Plots.scatter(xlabel="PC1", ylabel="PC2")
	for i=1:size(PhiFilter)[1]
		Predict=predict(PCA_kcentric,Matrix(PhiFilter[i][!,8:10])')
		Plots.scatter!(KcentricPCA_Plot_beta, Predict[1,:],Predict[2,:], label=PhiFilter[i].beta0[1], legend=:outerright, markers=Plots.supported_markers()[i], xlims=(), title="Phase Ratio")
	end		
KcentricPCA_Plot_beta
end

# ╔═╡ e2271c7c-a808-4e2e-b2e0-76393a4fcb67
predict(PCA_kcentric,Matrix(SubstanceFilter[1][!,7:9])')

# ╔═╡ b5abad47-c35c-4763-958d-d1c0b6af6fae
#File=CSV.File("E:\\Tillman\\Doktorarbeit\\Arbeit\\Labor\\Messdaten\\RXI-5SilMS05\\FAMES\\Brehmer2022_lnk_T_Rxi5SilMS_beta125.CSV", delim=',', decimal='.')

# ╔═╡ b7dff30b-541b-45b6-831e-b00ed5125a44
#CSV.write("E:\\Tillman\\Doktorarbeit\\Arbeit\\Labor\\Messdaten\\RXI-5SilMS05\\FAMES\\Brehmer2022_lnk_T_Rxi5SilMS_beta125_2.csv", File)

# ╔═╡ bc5f3e11-d5de-4343-b959-86431b8f04f1
md"""# Lambert W function"""

# ╔═╡ 23a8d4b2-dd8d-4a30-96ad-26de84f7c2a9
begin 
	gr()
Plots.plot(-1/exp(1):0.01:6, lambertw.(-1/exp(1):0.01:6, 0), label=L"W_{0}", xlabel=L"x", ylabel=L"W(x)")	
	plot!(-1 ./exp.(1):0.01:0, lambertw.(-1 ./exp.(1):0.01:0, -1), label=L"W_{-1}", dpi=500)
end

# ╔═╡ 7c23f931-9aab-4e9e-9ea8-dd6074db28ba
md"""# END"""

# ╔═╡ 76c115ae-8482-4e1d-a4ff-873beb68cb41
scatter(nfl.A, nfl.B, nfl.C, label="not flagged", xlabel="A", ylabel="B", zlabel="C", camera=(30,15));

# ╔═╡ 1c6fec9b-6719-4b11-8e5b-52d874920b8b
md"""
## ChemicalIdentifiers.jl
""";

# ╔═╡ 1e72ad83-50c4-4acc-8a35-39d6431dded8
CI = RetentionData.substance_identification(alldata);

# ╔═╡ 42aeba5d-772d-4366-8bab-51a9c0875917
rimgnumbers = RetentionData.ring_number.(CI.smiles);

# ╔═╡ a39e1661-4765-411c-a062-233a64770391
md""" 
### Extract number of rings from SMILES
""";

# ╔═╡ 4b59a4a6-4036-4563-bc06-4b6942ef4d00
search_chemical("555-44-2");

# ╔═╡ aa7bd88a-4722-493b-9c1d-2e69cf4333f3
md"""
### Extract number of elements from formula
""";

# ╔═╡ 311a09e1-7cfd-42bf-be24-0315b6727d2e
element_numbers = RetentionData.formula_to_dict.(CI.formula);

# ╔═╡ 823317eb-e0eb-465e-9bf5-a24e541445b8
# add CAS to alldata
alldata[!, "CAS"] = CI.CAS;

# ╔═╡ b54b46bb-d420-42a7-acc4-000f2177860d
RetentionData.formula_to_dict(CI.formula[end]);

# ╔═╡ bc72d1f0-e683-4fa0-aa17-701cc76e1472
alldata[!, "CAS"] = CI.CAS;

# ╔═╡ 593a3566-6df3-4b43-877c-f2524f274079
# only use substances with a CAS entry
alldata_f = filter([:CAS, :flag] => (x, y) -> ismissing(x)==false && isempty(y), alldata);

# ╔═╡ 5cf8c59b-f927-4904-bc26-21048ae3d252
filter([:CAS] => x -> ismissing(x), CI)  ;

# ╔═╡ b9e95ff4-f78a-41c9-9181-f3ca2e2c1b0d
md"""
### Filter for homologous series
""";

# ╔═╡ feb1de8a-3557-4c04-9623-490714359d05
homologous_series = DataFrame(CSV.File(joinpath(project, "data", "homologous.csv")));

# ╔═╡ cc5e8077-e4f4-4294-9ac3-3c3a1f847933
names(homologous_series);


# ╔═╡ 78c05a58-1122-45e4-9b85-45c96e2082cb
homologous_series."homologous series";

# ╔═╡ 14ca0919-6b19-41c6-a108-d2e7da3409f5
hs = unique(homologous_series[!,3]);

# ╔═╡ 5f83d618-1d7b-4c88-947a-88ffbaf74a43
function add_homologous_to_Cat!(newdata)
	hs = DataFrame(CSV.File(joinpath(project, "data", "homologous.csv")))
	CAS = Array{Array{String,1}}(undef, length(hs.CAS))
	for i=1:length(hs.CAS)
		CAS[i] = split(hs.CAS[i],',')
	end
	hs[!,"CAS"]	= CAS

	iCat = findall(occursin.("Cat", names(newdata)))
	for i=1:length(newdata.CAS)
		for j=1:length(hs.CAS)
			if newdata.CAS[i] in hs.CAS[j]
				if  ismissing(!(hs."homologous series"[j] in newdata[i, iCat])) || !(hs."homologous series"[j] in newdata[i, iCat])# group not allready in Cat
					# find the first Cat column with missing entry
					if isnothing(findfirst(ismissing.(collect(newdata[i,iCat]))))
						ii = iCat[end] + 1
					else
						ii = iCat[findfirst(ismissing.(collect(newdata[i,iCat])))]
					end
					newdata[i,ii] = hs."homologous series"[j]
				end
			end
		end
	end
	return newdata
end;

# ╔═╡ 225892b1-783f-471e-ad64-78886c26f58f
### Same categories for same substances;

# ╔═╡ 39cc1ff1-9dc3-4b24-a629-f79e11fad45c
parset = "Kcentric";

# ╔═╡ b14181e3-4edf-42a2-ac6d-c08b3ad738db
md"""
## Some Filters

- exclude flagged data
- exclude data with `CAS = missing`
""";

# ╔═╡ cc7d709d-a0ee-4917-91ea-0c85dabace86
md"""
### Filter for stationary phase and phase ratio
""";


# ╔═╡ cb032eb8-f630-46ae-9cdc-33905803b99a
filter_φ = 0.002;

# ╔═╡ 63ee1dea-332d-4922-b3bf-4aadbb922d12
filter_sp = "Rxi5ms";

# ╔═╡ 1788e068-9eaa-4909-9a2d-fdce2a83760c
# [x] filter for not identified substances -> alternative names (save in `shortnames.csv`) or they are not in the database of ChemicalIdentifiers.jl (in this case add a separate database with name, CAS, formula, MW, Smiles -> missing.csv)

# ╔═╡ 71a08e27-ba1d-4955-98ee-7a87943fe39d
# - Name, CAS, Phase, Tchar, thetachar, DeltaCp, phi0, Source, Cat_1, Cat_2, ...
# - Name, CAS, Phase, A, B, C, phi0, Source, Cat_1, Cat_2, ...

# ╔═╡ 645433a1-3e7e-4820-98b8-e398e6a35b4a
alldata_f_ = RetentionData.align_categories(alldata_f);

# ╔═╡ 50f1d5b3-2251-455a-8d7d-8bddfeed3611
dup_data, dup_entry = RetentionData.duplicated_data(alldata_f_)

# ╔═╡ 46e15b57-5192-45bf-837c-4c2b9adc5d4f
begin
	selected_dup_data = DataFrame()
	
	for i=1:length(dup_data)
		sources = unique(dup_data[i].Source)
		cols = names(dup_data[i])
		if "Blumberg2017" in sources # use Blumberg2017 data
			j = findfirst(dup_data[i].Source.=="Blumberg2017")
			push!(selected_dup_data, dup_data[i][j,cols])
		elseif ["Marquart2020"] == sources || ["Brehmer2022"] == sources # use all data, if only duplicates from Marquart (different isomers) or from Brehemer (different df)
			append!(selected_dup_data, dup_data[i][!,cols])
		elseif ["McGinitie2011"] == sources # different gases, use "He"
			j = findfirst(dup_data[i].gas.=="He")
			push!(selected_dup_data, dup_data[i][j,cols])
		elseif "McGinitie2011" in sources && "McGinitie2014a" in sources # different diameters, use d=0.25mm
			j = findfirst(dup_data[i].d.==0.25)
			push!(selected_dup_data, dup_data[i][j,cols])
		elseif length(sources) == 1 # duplicates from the same source
			# choose first
			if i == 1
				selected_dup_data =DataFrame(dup_data[i][1,cols])
			else
				push!(selected_dup_data, dup_data[i][1,cols]) 
			end
		else # not decided cases
			append!(selected_dup_data, dup_data[i][!,cols]) # choose all
			#push!(selected_dup_data, dup_data[i][1,cols]) # choose first -> error
		end
	end
	selected_dup_data

	nondup_data = alldata_f_[findall(dup_entry.==false),Not([:d, :gas, :flag])];
# 2nd filter the duplicate data dup_data according to decisions
# change the rule, if needed 

	dup_data_1, dup_entry_1 = RetentionData.duplicated_data(selected_dup_data)
# -> remaining cases of duplicates (some will stay, e.g. Marquart2020);
# 1st take only the entrys without duplicates, delete columns 'd', 'gas' and 'flag'
end;

# ╔═╡ 4a73771a-3adc-4fad-8d6a-148dcf9cc3c4
# 3rd combine both dataframes and sort for "source", "phase", "Tchar"
begin
	newdata = sort!(unique(vcat(nondup_data, selected_dup_data, cols = :union)), [:Source, :Phase, :Tchar])
	newdata
end;

# ╔═╡ 3d0b0729-c7a6-408a-8e2c-b5255347a65e
old_db = RetentionData.old_database_format(newdata);

# ╔═╡ 8397c671-c0d1-4632-ae9a-55a6dccd1002
new_db = RetentionData.new_database_format(newdata; ParSet=parset);

# ╔═╡ b28198b2-e649-4899-abb6-d4062454aa5e
new_db_filter = filter([:Phase, :phi0] => (x, y) -> x == filter_sp && y == filter_φ, new_db);

# ╔═╡ baa7d024-ec90-4744-922d-830f40683abe
dup_data_2, dup_entry_2 = RetentionData.duplicated_data(newdata)
# only duplicated data from Marquart2020 remain

# ╔═╡ 3ebdebe8-ad66-4721-b9ff-92da125bcf7c
sort(filter_Cat(newdata, "Grob"), [:Phase]);

# ╔═╡ 4cab40dc-4813-48f4-a8a5-dedb938cd25e
add_homologous_to_Cat!(newdata);
	# add CAS to alldata

# ╔═╡ 330d2dec-dd35-4469-aead-3584fa9f7a36
filter_Cat(newdata, hs[3]);

# ╔═╡ dd48df07-4bed-47ce-9799-05958e3adc7a
count(dup_entry)

# ╔═╡ 87a761c6-85ce-4342-a7d4-a13a378c6c45
md"""
### Plot of duplicates
$(@bind select_dup Slider(1:length(dup_data); show_value=true))
"""

# ╔═╡ bac7d9ef-d6bd-4955-b8e1-dfa6f6f7f038
begin
	T = 0.0:1.0:400.0
	Tst = 273.15
	R = 8.31446261815324
	lnk = Array{Float64}(undef, length(T), size(dup_data[select_dup])[1])
	plnk_dup = plot(xlabel="temperature in °C", ylabel="lnk")
	for j=1:size(dup_data[select_dup])[1]
		for i=1:length(T)
			par = [dup_data[select_dup].Tchar[j]+Tst, dup_data[select_dup].thetachar[j], dup_data[select_dup].DeltaCp[j]/R]
			lnk[i,j] = RetentionData.Kcentric(T[i]+Tst, par)
		end
		plot!(plnk_dup, T, lnk[:,j], title="duplicated data, $(dup_data[select_dup].Name[1]), $(dup_data[select_dup].Phase[1])", label=dup_data[select_dup].Source[j])
	end
end

# ╔═╡ 9a09e9cb-26cf-4576-a581-7b832fbab775
pKcentric_dup = scatter(dup_data[select_dup].Tchar, dup_data[select_dup].thetachar, dup_data[select_dup].DeltaCp, title="duplicated data, $(dup_data[select_dup].Name[1]), $(dup_data[select_dup].Phase[1])", label="", xlabel="Tchar", ylabel="θchar", zlabel="ΔCp");

# ╔═╡ fecbe30c-d7a9-41d2-8044-2ca4c4034ca3
pABC_dup = scatter(dup_data[select_dup].A, dup_data[select_dup].B, dup_data[select_dup].C, title="duplicated data, $(dup_data[select_dup].Name[1]), $(dup_data[select_dup].Phase[1])", label="", xlabel="A", ylabel="B", zlabel="C");

# ╔═╡ ba1c7f65-3e24-4b86-b2e2-521459993250
md"""
$(embed_display(plnk_dup))
$(embed_display(pABC_dup))
$(embed_display(pKcentric_dup))
"""

# ╔═╡ 610d535a-2025-4419-b35b-8d28dbaa62b8
RetentionData.add_group_to_Cat!(alldata_f_);

# ╔═╡ a60013a8-2337-4fac-a338-a0f3f78996d1
#CSV.write("../Databases/newformat_$(parset)_$(filter_sp)_beta$(1/(4*filter_φ)).csv", new_db_filter)

# ╔═╡ 845e5231-e549-42ca-b43e-8ddddc3d6446
# Name, CAS, Cnumber, Hnumber, Onumber, Nnumber, Ringnumber, Molmass, Phase, Tchar, thetachar, DeltaCp, phi0, Annotation

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ChemicalIdentifiers = "fa4ea961-1416-484e-bda2-883ee1634ba5"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LambertW = "984bce1d-4616-540c-a9ee-88d1112d94c9"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
MultivariateStats = "6f286f6a-111f-5878-ab1e-185364afe411"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RDatasets = "ce6b1742-4840-55fa-b093-852dadbb1d8b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.9"
ChemicalIdentifiers = "~0.1.7"
DataFrames = "~1.5.0"
LaTeXStrings = "~1.3.0"
LambertW = "~0.4.6"
LsqFit = "~0.13.0"
Measurements = "~2.8.0"
MultivariateStats = "~0.10.0"
Plots = "~1.38.5"
PlutoUI = "~0.7.50"
RDatasets = "~0.7.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0310e08cb19f5da31d08341c6120c047598f5b9c"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4d9946e51e24f5e509779e3e2c06281a733914c2"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.1.0"

[[deps.Arrow]]
deps = ["ArrowTypes", "BitIntegers", "CodecLz4", "CodecZstd", "DataAPI", "Dates", "LoggingExtras", "Mmap", "PooledArrays", "SentinelArrays", "Tables", "TimeZones", "UUIDs", "WorkerUtilities"]
git-tree-sha1 = "4e40f4868281b7fd702c605c764ab82a52ac3f4b"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "2.4.3"

[[deps.ArrowTypes]]
deps = ["UUIDs"]
git-tree-sha1 = "563d60f89fcb730668bd568ba3e752ee71dde023"
uuid = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
version = "2.0.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitIntegers]]
deps = ["Random"]
git-tree-sha1 = "fc54d5837033a170f3bad307f993e156eefc345f"
uuid = "c3b6d118-76ef-56ca-8cc7-ebb389d030a1"
version = "0.2.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "5084cc1a28976dd1642c9f337b28a3cb03e0f7d2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.7"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.ChemicalIdentifiers]]
deps = ["Arrow", "Downloads", "Preferences", "Scratch", "UUIDs", "Unicode"]
git-tree-sha1 = "f55715e75fb14cbe72808d43b75c4d66e5200daf"
uuid = "fa4ea961-1416-484e-bda2-883ee1634ba5"
version = "0.1.7"

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
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

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
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9258430c176319dc882efa4088e2ff882a0cb1f1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.81"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

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
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "d3ba08ab64bdfd27234d3f61956c966266757fe6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.7"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "ed1b56934a2f7a65035976985da71b6a65b4f2cf"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.18.0"

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
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

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
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "660b2ea2ec2b010bb02823c6d0ff6afd9bdc5c16"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d5e1fd17ac7f3aa4c5287a61ee28d4f8b8e98873"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

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
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "0cf92ec945125946352f3d46c96976ab972bde6f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.3.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

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
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

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
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

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
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "00f475f85c50584b12268675072663dfed5594b2"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.13.0"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "12950d646ce04fb2e89ba5bd890205882c3592d7"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.8.0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "c272302b22479a24d1cf48c114ad702933414f80"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.5"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "efe9c8ecab7a6311d4b91568bd6c88897822fabe"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

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

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

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

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "6466e524967496866901a78fca3f2e9ea445a559"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "c95373e73290cf50a8a22c3375e4625ded5c5280"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "8ac949bd0ebc46a44afb1fdca1094554a84b086e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "786efa36b7eff813723c4849c90456609cf06661"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.1"

[[deps.RData]]
deps = ["CategoricalArrays", "CodecZlib", "DataFrames", "Dates", "FileIO", "Requires", "TimeZones", "Unicode"]
git-tree-sha1 = "19e47a495dfb7240eb44dc6971d660f7e4244a72"
uuid = "df47a6cb-8c03-5eed-afd8-b6050d6c41da"
version = "0.8.3"

[[deps.RDatasets]]
deps = ["CSV", "CodecZlib", "DataFrames", "FileIO", "Printf", "RData", "Reexport"]
git-tree-sha1 = "2720e6f6afb3e562ccb70a6b62f8f308ff810333"
uuid = "ce6b1742-4840-55fa-b093-852dadbb1d8b"
version = "0.7.7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "18c35ed630d7229c5584b945641a73ca83fb5213"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.2"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "2d7d9e1ddadc8407ffd460e24218e37ef52dd9a3"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.16"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Scratch", "Unicode"]
git-tree-sha1 = "a92ec4466fc6e3dd704e2668b5e7f24add36d242"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.9.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "e4bdc63f5c6d62e80eb1c0043fcc0360d5950ff7"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.10"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

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
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6edfe154ad7b313c01aceca188c05c835c67360"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.4+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

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
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─c4ea2bc1-ec27-47b4-9a4b-0c03da27042b
# ╠═fec530e6-d675-4bb9-8b5a-6aa607574a81
# ╟─313457f3-85fd-4274-9b96-fceb277eaae4
# ╟─ec6342ed-6f00-4c7f-93d2-21a5f39eb4db
# ╟─5861512e-a14b-11ec-3c6b-9bd2953bf909
# ╟─e92ee318-9f88-40f2-81cf-e2d60bf2f45a
# ╟─3b9c6610-6839-4012-aa0b-219a347ca52f
# ╟─d26fc674-ace5-43ef-af1d-855dfc21eba5
# ╟─979d0e68-135a-4570-a835-8be8e7095326
# ╟─3d0b0729-c7a6-408a-8e2c-b5255347a65e
# ╟─8397c671-c0d1-4632-ae9a-55a6dccd1002
# ╟─8bc510bf-b40a-4f7f-b5e9-05bdca6fe1a3
# ╟─91d21ca2-5657-4f9b-ae66-629fa77236f1
# ╟─b2b6ef1c-ca2d-4a52-a586-b3c2b13b6bac
# ╟─84fc7b22-2362-4b99-8fd9-3a779d706a4b
# ╟─93a5ed9f-64ce-4ee8-bee0-ecd2ff2dcbb8
# ╟─1e912349-d26d-45ea-9167-21b0b4165060
# ╟─97d4ff70-092c-40db-90a2-6a8994ad920c
# ╟─7c95815d-11e5-4de5-8d83-a7ef8518751c
# ╟─d28a4b20-af06-4613-810b-dafb9a5f0ccd
# ╟─c28cca6a-fdf3-4edf-9534-5c4f677c2889
# ╟─32c1889e-9da1-44e6-9af9-08af1d7f2020
# ╟─7eb28537-4d9e-4529-b6f7-b0407bce814c
# ╟─d460b5e1-9060-4fe5-bcea-8c63a0ea071e
# ╟─17cb9134-8880-4cb2-86a1-31eb4e1ee5b5
# ╟─17387b6a-e126-4923-ac14-229d980ed9ad
# ╟─0e559fd4-5fef-4375-9e7a-1e052cf61b7a
# ╟─0ad0450e-199a-4be1-a6af-03eb745f08d5
# ╟─8287590e-0e6d-4591-a616-0b64d84ca767
# ╟─ef518fe9-1920-4c95-9ef3-2b2542621fc0
# ╟─46e15b57-5192-45bf-837c-4c2b9adc5d4f
# ╟─917e88c6-cf5d-4d8f-93e5-f50ea2bb2cdc
# ╟─bac7d9ef-d6bd-4955-b8e1-dfa6f6f7f038
# ╟─baa7d024-ec90-4744-922d-830f40683abe
# ╟─496b86c5-900e-4786-a254-082b8155c65b
# ╟─a319c1af-33d4-443d-9270-5a0219e25c4c
# ╟─9a09e9cb-26cf-4576-a581-7b832fbab775
# ╟─fda7f746-0f7e-4dff-b8e6-7765f580e542
# ╠═dd48df07-4bed-47ce-9799-05958e3adc7a
# ╟─50f1d5b3-2251-455a-8d7d-8bddfeed3611
# ╟─fecbe30c-d7a9-41d2-8044-2ca4c4034ca3
# ╟─87a761c6-85ce-4342-a7d4-a13a378c6c45
# ╟─ba1c7f65-3e24-4b86-b2e2-521459993250
# ╟─9290687e-26d1-4ccf-9603-d46b8ae9dcb7
# ╠═eaed9fdd-d4a2-41f7-8e05-588869e12780
# ╠═070c0b3b-efa4-4f4d-a567-87ba4b7c936b
# ╟─ef0c1872-da54-48b5-8212-634d7d91e9ac
# ╟─873b90e2-bc5d-4272-a90b-978bb4e1358d
# ╟─d7acb6b7-539a-4d4a-bbff-d066e6535ff6
# ╟─d24832da-ff7d-4d4c-af5d-6662eca846cf
# ╟─a94f22a9-5332-429a-bd31-3c1799781ffe
# ╟─b846be08-4850-4bbd-9366-b7e2207bb370
# ╟─7b39d781-b51c-4406-af76-ce9459b1cdb6
# ╟─c365e82d-d5e5-4be0-81c9-c5a8097d4e8c
# ╟─512ae246-792c-45c7-ae58-406bc70d47ff
# ╟─e945beeb-a9de-4a2d-bc8b-73c3f3dd37b1
# ╟─20eb935e-c66d-46fe-b309-e4f37f29e667
# ╟─b8c51d43-fd7d-440e-b516-a49a21cf44b8
# ╟─015e4eb0-f60e-474f-9ab9-98f6c44b9d8f
# ╟─bbaa8bd7-5b38-4fda-aab9-d97d0c14c410
# ╟─44c3be0a-3376-4b51-bd07-683e96029d69
# ╟─33ec66d1-2a25-45b6-936c-0acff05e624c
# ╟─6729a636-aa06-4695-8fd6-2078322e2ffa
# ╟─99c115a1-7bbf-4aa1-b526-6200e4d8a622
# ╟─1d708017-f3d8-4126-9808-3cbc7bce0050
# ╟─421e55ec-a18b-4aba-9c88-0e53cdcef75d
# ╟─fb2bf9b6-e7e2-4b5f-8013-e6db5d779f85
# ╟─21277d27-7ea7-4adc-b1b8-4429c29cdf51
# ╟─9b0d91eb-6672-4cf8-82e4-844becf699bb
# ╟─ed134625-a178-48d8-856f-e8b25f57842d
# ╟─4b5de04c-24ef-4f79-937a-f2816d7397d7
# ╟─3b0e9656-63e8-4fa3-8f06-33f8a2fcd4dc
# ╟─2d89e6eb-4c6e-4c6f-8a82-55eaf66de976
# ╟─09b949bf-36b5-44a6-ab13-e0074b824756
# ╟─13a8c60e-979b-47e6-bcdc-7a314fd27ebb
# ╟─b5c0895f-1683-4f46-960b-9463819610a1
# ╟─fc4a5a98-f1ab-45b5-9c2d-44693ad35313
# ╟─a64f174b-fac7-4566-b8b1-2c603dbd3dd1
# ╠═492fabbd-c84a-421e-a169-4229cf444e11
# ╟─80745a29-6a58-47c8-a176-e6f58a858d82
# ╟─333ffd2d-d4db-46fd-b3f3-5cbffeb5243d
# ╟─bf9adfcc-5baf-486d-b606-4f1a8d935a17
# ╟─156b9224-3714-4039-aa49-ecc6b04d38dc
# ╟─202ab200-b358-4b3b-8bd9-da3a1c158d39
# ╟─f687a423-bf8f-4c5f-b98e-f8afc5c4ee9e
# ╠═391282b0-f684-4e3d-b4ba-1eba0689fe5e
# ╟─ad5312d2-43bb-42f6-8327-a3de2dfc1aae
# ╠═4ec1f88e-d2e4-43c0-9e43-68322ab96e1c
# ╠═534ec747-4b14-4c01-a837-7e6db2479381
# ╠═982446df-5021-44ca-b7a8-2f0042baf88d
# ╠═b6fcf3bb-9b33-42ba-97fa-f0116ae0ec46
# ╠═73611573-ae96-49d5-86be-8b646cd5c8bc
# ╠═09806af8-d0d8-49a9-a2fe-4357992584b8
# ╟─8a4541fc-7ebb-46f0-a652-10f7ba1ea689
# ╠═52f944bf-2038-497a-ae39-a2ee95ad67d2
# ╠═32471c63-a5c0-4ad9-b9e4-1701a8a9902b
# ╠═fb7017ab-b8bb-4175-8434-5f9c1e026509
# ╠═0698e398-d204-448e-9460-5eb2a51a5d05
# ╠═ff6f166d-6d18-44c1-8be3-8d6d40bcd738
# ╟─4badc718-b4a0-4642-adcb-b88a55193280
# ╟─00703bb9-21f6-4f93-95f4-45995d81ee5d
# ╟─8f2f176e-cab0-4ab8-8d9c-f635d6e95460
# ╟─c00ded9a-6d20-4748-9dc6-789831820427
# ╟─e2271c7c-a808-4e2e-b2e0-76393a4fcb67
# ╟─b5abad47-c35c-4763-958d-d1c0b6af6fae
# ╟─b7dff30b-541b-45b6-831e-b00ed5125a44
# ╟─bc5f3e11-d5de-4343-b959-86431b8f04f1
# ╟─23a8d4b2-dd8d-4a30-96ad-26de84f7c2a9
# ╟─7c23f931-9aab-4e9e-9ea8-dd6074db28ba
# ╟─76c115ae-8482-4e1d-a4ff-873beb68cb41
# ╟─3ebdebe8-ad66-4721-b9ff-92da125bcf7c
# ╟─1c6fec9b-6719-4b11-8e5b-52d874920b8b
# ╟─1e72ad83-50c4-4acc-8a35-39d6431dded8
# ╟─610d535a-2025-4419-b35b-8d28dbaa62b8
# ╟─42aeba5d-772d-4366-8bab-51a9c0875917
# ╟─a39e1661-4765-411c-a062-233a64770391
# ╟─4b59a4a6-4036-4563-bc06-4b6942ef4d00
# ╟─aa7bd88a-4722-493b-9c1d-2e69cf4333f3
# ╟─311a09e1-7cfd-42bf-be24-0315b6727d2e
# ╟─823317eb-e0eb-465e-9bf5-a24e541445b8
# ╟─b54b46bb-d420-42a7-acc4-000f2177860d
# ╟─bc72d1f0-e683-4fa0-aa17-701cc76e1472
# ╟─593a3566-6df3-4b43-877c-f2524f274079
# ╟─4cab40dc-4813-48f4-a8a5-dedb938cd25e
# ╟─5cf8c59b-f927-4904-bc26-21048ae3d252
# ╟─b9e95ff4-f78a-41c9-9181-f3ca2e2c1b0d
# ╟─feb1de8a-3557-4c04-9623-490714359d05
# ╟─cc5e8077-e4f4-4294-9ac3-3c3a1f847933
# ╟─78c05a58-1122-45e4-9b85-45c96e2082cb
# ╟─14ca0919-6b19-41c6-a108-d2e7da3409f5
# ╟─5f83d618-1d7b-4c88-947a-88ffbaf74a43
# ╟─225892b1-783f-471e-ad64-78886c26f58f
# ╟─330d2dec-dd35-4469-aead-3584fa9f7a36
# ╟─4a73771a-3adc-4fad-8d6a-148dcf9cc3c4
# ╟─39cc1ff1-9dc3-4b24-a629-f79e11fad45c
# ╟─b14181e3-4edf-42a2-ac6d-c08b3ad738db
# ╟─cc7d709d-a0ee-4917-91ea-0c85dabace86
# ╟─cb032eb8-f630-46ae-9cdc-33905803b99a
# ╟─63ee1dea-332d-4922-b3bf-4aadbb922d12
# ╟─b28198b2-e649-4899-abb6-d4062454aa5e
# ╟─1788e068-9eaa-4909-9a2d-fdce2a83760c
# ╟─71a08e27-ba1d-4955-98ee-7a87943fe39d
# ╟─645433a1-3e7e-4820-98b8-e398e6a35b4a
# ╟─a60013a8-2337-4fac-a338-a0f3f78996d1
# ╟─845e5231-e549-42ca-b43e-8ddddc3d6446
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
