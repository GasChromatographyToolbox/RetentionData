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

# ╔═╡ 5861512e-a14b-11ec-3c6b-9bd2953bf909
md"""
# Load `AllParam` data
- Load all files with `AllParam`  from the folder `db_path` (and its subfolders).
- use ChemicalIdentifiers.jl
- combine all data into one DataFrame/csv-file
"""

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
	TableOfContents()
end

# ╔═╡ 3757151c-2244-4e45-985c-2ef869abd23d
data = RetentionData.load_allparameter_data(db_path)

# ╔═╡ 3b9c6610-6839-4012-aa0b-219a347ca52f
data.data[end]

# ╔═╡ d26fc674-ace5-43ef-af1d-855dfc21eba5
begin
	alldata = RetentionData.dataframe_of_all(data)
	# add flags to alldata
	alldata[!, "flag"] = RetentionData.flag(alldata)
	alldata
end

# ╔═╡ 070c0b3b-efa4-4f4d-a567-87ba4b7c936b
unique(alldata.Phase)

# ╔═╡ 3c860294-55cd-4b8d-8fa2-46583793fe00
md"""
## Flag the parameter sets
"""

# ╔═╡ 5d719450-b016-4cb5-b4a7-2f0e71eba5f3
fl, nfl = RetentionData.flagged_data(alldata)

# ╔═╡ 61bfdd58-cc7d-4723-9f92-7ada1d776bcc


# ╔═╡ 93a5ed9f-64ce-4ee8-bee0-ecd2ff2dcbb8
md"""
## Plot parameters A, B, C
"""

# ╔═╡ 560c76d6-9e50-461a-9815-66b40b59e580
begin
	plotly()
	pABC = scatter(nfl.A, nfl.B, nfl.C, label="not flagged", xlabel="A", ylabel="B", zlabel="C")
	scatter!(pABC, fl.A, fl.B, fl.C, label="flagged", c=:red, m=:cross)
	pABC
end

# ╔═╡ a56ed363-4cc7-4471-b0aa-34109d2dcb45
pABC_nfl = scatter(nfl.A, nfl.B, nfl.C, label="not flagged", xlabel="A", ylabel="B", zlabel="C")

# ╔═╡ 76c115ae-8482-4e1d-a4ff-873beb68cb41
scatter(nfl.A, nfl.B, nfl.C, label="not flagged", xlabel="A", ylabel="B", zlabel="C", camera=(30,15))

# ╔═╡ bfc0bd73-94c2-4742-b334-e632933d6dfe
md"""
### Observations

The parameters of the ABC-model span a plane in the ABC-space. With a estimated plane equation the divergence of single data points can be evaluated. This could also be used in the estimation of parameters.

#### ToDo:
- fit a plane equation to the points
- estimate the residua of single points
"""

# ╔═╡ 49ddd8eb-a833-43b8-9f62-18f4d5afaf10
md"""
## Plot parameters Tchar, θchar, ΔCp
"""

# ╔═╡ c99ca5d6-7f7d-4462-afab-e11154370054
begin
	plotly()
	pKcentric = scatter(nfl.Tchar, nfl.thetachar, nfl.DeltaCp, label="not flagged", xlabel="Tchar", ylabel="thetachar", zlabel="DeltaCp")
	scatter!(pKcentric, fl.Tchar, fl.thetachar, fl.DeltaCp, label="flagged", c=:red, m=:cross)
	pKcentric
end

# ╔═╡ e9f134c3-25e7-45a8-a07e-a3cfdc6c027b
pKcentric_nfl = scatter(nfl.Tchar, nfl.thetachar, nfl.DeltaCp, label="not flagged", xlabel=L"T_{char}", ylabel=L"\theta_{char}", zlabel=L"\Delta C_p", camera=(30, 30))

# ╔═╡ d24832da-ff7d-4d4c-af5d-6662eca846cf
begin gr() 
scatter(nfl.Tchar, nfl.thetachar, label="not flagged", xlabel=L"T_{char}", ylabel=L"\theta_{char}")
end	

# ╔═╡ a94f22a9-5332-429a-bd31-3c1799781ffe
scatter(nfl.Tchar, nfl.DeltaCp, label="not flagged", xlabel=L"T_{char}", ylabel=L"ΔC_p", camera=(30, 45))

# ╔═╡ 326a0c21-454d-489f-b960-be356671e1db
md"""
### Observations

The parameters of the K-centric-model do **not** span a plane in the parameter-space. But it could still be worth to fit a plane into this data cloud.

#### ToDo:
- fit a plane equation to the points
- estimate the residua of single points
"""

# ╔═╡ b4044a75-03ad-4f75-a2ac-716b0c2c628f
md"""
## Plot parameters ΔHref, ΔSref, ΔCp
"""

# ╔═╡ 80c784b2-35a6-4fc9-a1fd-7b1d6d89462f
begin
	plotly()
	pTD = scatter(nfl.DeltaHref, nfl.DeltaSref, nfl.DeltaCp, label="not flagged", xlabel="DeltaHref", ylabel="DeltaSref", zlabel="DeltaCp")
	scatter!(pTD, fl.DeltaHref, fl.DeltaSref, fl.DeltaCp, label="flagged", c=:red, m=:cross)
	pTD
end

# ╔═╡ 6260e6e0-02e4-4c95-9889-020e5c3c2d60
pTD_nfl = scatter(nfl.DeltaHref, nfl.DeltaSref, nfl.DeltaCp, label="not flagged", xlabel="DeltaHref", ylabel="DeltaSref", zlabel="DeltaCp")

# ╔═╡ 7c95815d-11e5-4de5-8d83-a7ef8518751c
begin 
	pTD_nfl_ = plot(xlabel="DeltaHref", ylabel="DeltaSref", zlabel="DeltaCp")
	source = unique(nfl.Source)
	for i=1:length(source)
		nfl_filter = filter([:Source] => x -> x == source[i], nfl)
		scatter!(pTD_nfl_, nfl_filter.DeltaHref, nfl_filter.DeltaSref, nfl_filter.DeltaCp, label=source[i])
	end
	pTD_nfl_
end

# ╔═╡ 32fd5862-e14b-4cea-b99c-0f26e0d8fbb5
md"""
### Observations

The parameters of the Thermodynamic-model do span two planes in the parameter-space. One plane consits soly of data from the group of Harynuk (Karolat, McGinitie). The other plane consists of other sources, including ower one measurements. Interestingly, the parameters from Blumberg are in the second plane, while the data from which they are derived (Karolat) are in the first plane.

The data from the sources of Karolat2010, McGinitie2011, McGinitie2012a, McGinitie2014a and McGinitie2014b are given in the Thermodynamic-model (ΔHref, ΔSref, ΔCp).

The data from Blumberg2017 is given in the ABC and Kcentric parameter sets and with enthalpy/entropy change at Tchar as reference temperature. Therefore, the calculation of refereence entropy/entalpy change at the reference temperature of 90°C should be re-evaluated. 

#### ToDo:
- check the calculation of ΔHref and ΔSref
"""

# ╔═╡ 453cbd8f-8f98-4bc7-8577-877455d8354e
md"""
## ChemicalIdentifiers.jl
"""

# ╔═╡ 5ee8d8f2-3c3c-4d4a-a14d-e9a0f23c3ef5
# [x] filter for not identified substances -> alternative names (save in `shortnames.csv`) or they are not in the database of ChemicalIdentifiers.jl (in this case add a separate database with name, CAS, formula, MW, Smiles -> missing.csv)

# ╔═╡ 4129a813-51b6-4790-9f20-a0d5e188b5c7
CI = RetentionData.substance_identification(alldata)

# ╔═╡ 5cf8c59b-f927-4904-bc26-21048ae3d252
filter([:CAS] => x -> ismissing(x), CI)  

# ╔═╡ 31a4dc3a-0b29-45b1-9876-47bd082e72bb
# add CAS to alldata
alldata[!, "CAS"] = CI.CAS

# ╔═╡ 4000c057-c75a-4bd8-93fd-dcfce907101c
search_chemical("555-44-2")

# ╔═╡ aa7bd88a-4722-493b-9c1d-2e69cf4333f3
md"""
## Extract number of elements from formula
"""

# ╔═╡ b54b46bb-d420-42a7-acc4-000f2177860d
RetentionData.formula_to_dict(CI.formula[end])

# ╔═╡ c28cca6a-fdf3-4edf-9534-5c4f677c2889
element_numbers = RetentionData.formula_to_dict.(CI.formula)

# ╔═╡ a39e1661-4765-411c-a062-233a64770391
md""" 
## Extract number of rings from SMILES
"""

# ╔═╡ e34c35f8-1ec8-47ba-af8e-77f10cf2c27b
rimgnumbers = RetentionData.ring_number.(CI.smiles)

# ╔═╡ 11de4e14-3b56-4238-bed8-c110f4de2d44
md"""
## Some Filters

- exclude flagged data
- exclude data with `CAS = missing`
"""

# ╔═╡ d619df94-5bb2-4349-a83d-ccfd13b95906
# only use substances with a CAS entry
alldata_f = filter([:CAS, :flag] => (x, y) -> ismissing(x)==false && isempty(y), alldata)

# ╔═╡ e17bf6fe-f3b6-4904-bb05-8b068fd9cf1f
md"""
## Same categories for same substances
"""

# ╔═╡ 54a2712b-d696-4097-8522-f5e1a87ecbec
alldata_f_ = RetentionData.align_categories(alldata_f)

# ╔═╡ d5381c2d-1794-4af1-9ebb-5188334fc592
md"""
## Add Category groups (e.g. BTEX, Grob, ...)
"""

# ╔═╡ 610d535a-2025-4419-b35b-8d28dbaa62b8
RetentionData.add_group_to_Cat!(alldata_f_)

# ╔═╡ 917e88c6-cf5d-4d8f-93e5-f50ea2bb2cdc
md"""
## Find duplicates
"""

# ╔═╡ e5b49869-2763-4ec9-ae7f-7b70164c0c67
dup_data, dup_entry = RetentionData.duplicated_data(alldata_f_)

# ╔═╡ dd48df07-4bed-47ce-9799-05958e3adc7a
count(dup_entry)

# ╔═╡ 47138dbf-5c60-4cdc-b484-ff168d11055f
dup_data[70]

# ╔═╡ 87a761c6-85ce-4342-a7d4-a13a378c6c45
md"""
### Plot of duplicates
$(@bind select_dup Slider(1:length(dup_data); show_value=true))
"""

# ╔═╡ 8b0b5a9f-666a-43c2-bdab-3d771334f12b
md"""
### Conclusions from duplicates
- if _Blumberg2017_ is one of the sources, keep these
- in the case of the three duplicates from _McGinitie2011_ different mobile phases where used (He, H₂ and N₂) resp. different column diameters (0.1mm, 0.2mm, 0.25mm, 0.32mm, 0.53mm; these are also duplicates with _McGinitie2014a_)-> similar results, keep the ones with He resp. 0.25mm
- duplicates from only _Marquart2020_ represent different observed isomers? related to the substance -> keep them all
- some complete duplicates (parameters have the same value) -> use only the first substance
"""

# ╔═╡ f02137a4-75a1-49ac-81a3-a1f04bab9ca2
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

# ╔═╡ c7f1acaa-d53e-427e-b646-f3920a2ce6b7
pABC_dup = scatter(dup_data[select_dup].A, dup_data[select_dup].B, dup_data[select_dup].C, title="duplicated data, $(dup_data[select_dup].Name[1]), $(dup_data[select_dup].Phase[1])", label="", xlabel="A", ylabel="B", zlabel="C");

# ╔═╡ 9a09e9cb-26cf-4576-a581-7b832fbab775
pKcentric_dup = scatter(dup_data[select_dup].Tchar, dup_data[select_dup].thetachar, dup_data[select_dup].DeltaCp, title="duplicated data, $(dup_data[select_dup].Name[1]), $(dup_data[select_dup].Phase[1])", label="", xlabel="Tchar", ylabel="θchar", zlabel="ΔCp");

# ╔═╡ ba1c7f65-3e24-4b86-b2e2-521459993250
md"""
$(embed_display(plnk_dup))
$(embed_display(pABC_dup))
$(embed_display(pKcentric_dup))
"""

# ╔═╡ fda7f746-0f7e-4dff-b8e6-7765f580e542
md"""
## Some more filters for the duplicates
"""

# ╔═╡ b7be7b57-7e69-4cc9-af54-db9465f18d05
# 1st take only the entrys without duplicates, delete columns 'd', 'gas' and 'flag'
nondup_data = alldata_f_[findall(dup_entry.==false),Not([:d, :gas, :flag])]

# ╔═╡ 81ccb3e6-dcca-4b66-9be6-a1fb63f7e056
# 2nd filter the duplicate data dup_data according to decisions
# change the rule, if needed 
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
end

# ╔═╡ a50d4f05-0f4a-4ae0-ba5d-a27bad3869e0
dup_data_1, dup_entry_1 = RetentionData.duplicated_data(selected_dup_data)
# -> remaining cases of duplicates (some will stay, e.g. Marquart2020)

# ╔═╡ 4a73771a-3adc-4fad-8d6a-148dcf9cc3c4
# 3rd combine both dataframes and sort for "source", "phase", "Tchar"
begin
	newdata = sort!(unique(vcat(nondup_data, selected_dup_data, cols = :union)), [:Source, :Phase, :Tchar])
	newdata
end

# ╔═╡ baa7d024-ec90-4744-922d-830f40683abe
dup_data_2, dup_entry_2 = RetentionData.duplicated_data(newdata)
# only duplicated data from Marquart2020 remain

# ╔═╡ 1b8f3f28-e612-40d0-9b09-136543cbb126
md"""
## DataFrame with the structure of database for GasChromatographySimulator.jl v0.3 and lower
"""

# ╔═╡ 390a66a8-e497-4a0f-b6c4-34487df787e9
# Name, CAS, Cnumber, Hnumber, Onumber, Nnumber, Ringnumber, Molmass, Phase, Tchar, thetachar, DeltaCp, phi0, Annotation

# ╔═╡ 23e0cf31-1c67-48d1-b014-26c1b44e04a8
old_db = RetentionData.old_database_format(newdata)

# ╔═╡ 3eeea8dd-de5c-4c73-8d82-4bdb4979d2b0
parset = "Kcentric"

# ╔═╡ eaed9fdd-d4a2-41f7-8e05-588869e12780
#CSV.write("../Databases/oldformat_$(parset).csv", old_db)

# ╔═╡ 9634b96c-18f1-479e-b28b-3d614893ce7c
md"""
## DataFrame with the new structure of database for GasChromatographySimulator.jl v0.4 and higher
"""

# ╔═╡ 58c8bd1d-2e69-48ac-9cfc-bb525ebe79c8
# - Name, CAS, Phase, Tchar, thetachar, DeltaCp, phi0, Source, Cat_1, Cat_2, ...
# - Name, CAS, Phase, A, B, C, phi0, Source, Cat_1, Cat_2, ...

# ╔═╡ 8397c671-c0d1-4632-ae9a-55a6dccd1002
new_db = RetentionData.new_database_format(newdata; ParSet=parset)

# ╔═╡ 496b86c5-900e-4786-a254-082b8155c65b
#CSV.write("../Databases/newformat_$(parset).csv", new_db)

# ╔═╡ fa394ad0-055e-4c8f-a3cf-931cd421cb7e
md"""
## Filter for stationary phase and phase ratio
"""

# ╔═╡ 897c8de6-c042-49aa-98ee-e5b3762756b8
filter_sp = "Rxi5ms"

# ╔═╡ abc5bab1-66c6-4ae8-95a8-3ba6d7eecbf5
filter_φ = 0.002

# ╔═╡ 1d23d9a9-996a-4c88-9bef-2fcc55fd5b44
new_db_filter = filter([:Phase, :phi0] => (x, y) -> x == filter_sp && y == filter_φ, new_db)

# ╔═╡ 4d44c3bf-910d-411e-b3b8-99c1c60d43b1
#CSV.write("../Databases/newformat_$(parset)_$(filter_sp)_beta$(1/(4*filter_φ)).csv", new_db_filter)

# ╔═╡ ccc85a17-690a-4fa4-9b14-9ca58a22e9c8
md"""
## Filter for homologous series
"""

# ╔═╡ 46e16092-d952-4c4f-a952-c5201797fcd1
homologous_series = DataFrame(CSV.File(joinpath(project, "data", "homologous.csv")))

# ╔═╡ f51eb4f5-3529-4c5a-8e2d-1902242ab7af
names(homologous_series)

# ╔═╡ 6dac8047-1170-4b65-9cee-2c8db3b62d63
homologous_series."homologous series"

# ╔═╡ e74f1990-5dd2-4062-a8d0-345a5005d0c2
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
end

# ╔═╡ ff046879-3607-47e4-afe3-b42bb1738b9f
add_homologous_to_Cat!(newdata)

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

# ╔═╡ 8d1c0954-ace2-4621-99aa-ce692936247b
hs = unique(homologous_series[!,3])

# ╔═╡ 5356fd33-c959-40b0-bdf8-d1363a0726c9
filter_Cat(newdata, hs[3])

# ╔═╡ 3ebdebe8-ad66-4721-b9ff-92da125bcf7c
sort(filter_Cat(newdata, "Grob"), [:Phase])

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
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.3"
ChemicalIdentifiers = "~0.1.5"
DataFrames = "~1.3.2"
LaTeXStrings = "~1.3.0"
LambertW = "~0.4.5"
LsqFit = "~0.12.1"
Measurements = "~2.7.1"
Plots = "~1.26.1"
PlutoUI = "~0.7.37"
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
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c46fb7dd1d8ca1d213ba25848a5ec4e47a1a1b08"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.26"

[[deps.Arrow]]
deps = ["ArrowTypes", "BitIntegers", "CodecLz4", "CodecZstd", "DataAPI", "Dates", "LoggingExtras", "Mmap", "PooledArrays", "SentinelArrays", "Tables", "TimeZones", "UUIDs", "WorkerUtilities"]
git-tree-sha1 = "3d04ab3584ece56c39397e01b55e1bd4fb8f0b30"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "2.4.1"

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
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "c5fd7cd27ac4aed0acf4b73948f0110ff2a854b2"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.7"

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
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

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
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "aaabba4ce1b7f8a9b34c015053d3b1edf60fa49c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.4.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

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
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "db2a9cb664fcea7836da4b414c3278d71dd602d2"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.6"

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
git-tree-sha1 = "a7756d098cbabec6b3ac44f369f74915e8cfd70a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.79"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

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

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "802bfc139833d2ba893dd9e62ba1767c88d708ae"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.5"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

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
git-tree-sha1 = "10fa12fe96e4d76acfa738f4df2126589a67374f"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.33"

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
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c98aea696662d09e215ef7cda5296024a9646c75"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.4"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bc9f7725571ddb4ab2c4bc74fa397c1c5ad08943"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.69.1+0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fb83fbe02fe57f2c068013aa94bcdf6760d3a7a7"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+1"

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
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

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
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

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
git-tree-sha1 = "dd8b9e6d7be9731fdaecc813acc5c3083496a251"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.7.2"

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

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

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
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

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
git-tree-sha1 = "b64719e8b4504983c7fca6cc9db3ebc8acc2a4d6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.1"

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
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "21303256d239f6b484977314674aef4bb1fe4420"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "2f041202ab4e47f4a3465e3993929538ea71bd48"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.26.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

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
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

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
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "d12e612bba40d189cead6ff857ddb67bd2e6a387"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

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
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "4e051b85454b4e4f66e6a6b7bdc452ad9da3dcf6"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.10"

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
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "13237798b407150a6d2e2bce5d793d7d9576e99e"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.13"

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
git-tree-sha1 = "d634a3641062c040fc8a7e2a3ea17661cc159688"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.9.0"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

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
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "164bf176a3bf20d037937a28f4ec9df940b26180"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.4.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

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
# ╟─5861512e-a14b-11ec-3c6b-9bd2953bf909
# ╟─a319c1af-33d4-443d-9270-5a0219e25c4c
# ╟─fec530e6-d675-4bb9-8b5a-6aa607574a81
# ╟─3757151c-2244-4e45-985c-2ef869abd23d
# ╠═3b9c6610-6839-4012-aa0b-219a347ca52f
# ╠═d26fc674-ace5-43ef-af1d-855dfc21eba5
# ╠═070c0b3b-efa4-4f4d-a567-87ba4b7c936b
# ╟─3c860294-55cd-4b8d-8fa2-46583793fe00
# ╠═5d719450-b016-4cb5-b4a7-2f0e71eba5f3
# ╠═61bfdd58-cc7d-4723-9f92-7ada1d776bcc
# ╟─93a5ed9f-64ce-4ee8-bee0-ecd2ff2dcbb8
# ╠═560c76d6-9e50-461a-9815-66b40b59e580
# ╟─a56ed363-4cc7-4471-b0aa-34109d2dcb45
# ╠═76c115ae-8482-4e1d-a4ff-873beb68cb41
# ╟─bfc0bd73-94c2-4742-b334-e632933d6dfe
# ╟─49ddd8eb-a833-43b8-9f62-18f4d5afaf10
# ╟─c99ca5d6-7f7d-4462-afab-e11154370054
# ╟─e9f134c3-25e7-45a8-a07e-a3cfdc6c027b
# ╟─313457f3-85fd-4274-9b96-fceb277eaae4
# ╟─d24832da-ff7d-4d4c-af5d-6662eca846cf
# ╟─a94f22a9-5332-429a-bd31-3c1799781ffe
# ╟─326a0c21-454d-489f-b960-be356671e1db
# ╟─b4044a75-03ad-4f75-a2ac-716b0c2c628f
# ╟─80c784b2-35a6-4fc9-a1fd-7b1d6d89462f
# ╠═6260e6e0-02e4-4c95-9889-020e5c3c2d60
# ╟─7c95815d-11e5-4de5-8d83-a7ef8518751c
# ╟─32fd5862-e14b-4cea-b99c-0f26e0d8fbb5
# ╟─453cbd8f-8f98-4bc7-8577-877455d8354e
# ╟─5ee8d8f2-3c3c-4d4a-a14d-e9a0f23c3ef5
# ╠═4129a813-51b6-4790-9f20-a0d5e188b5c7
# ╠═5cf8c59b-f927-4904-bc26-21048ae3d252
# ╠═31a4dc3a-0b29-45b1-9876-47bd082e72bb
# ╠═4000c057-c75a-4bd8-93fd-dcfce907101c
# ╠═aa7bd88a-4722-493b-9c1d-2e69cf4333f3
# ╠═b54b46bb-d420-42a7-acc4-000f2177860d
# ╠═c28cca6a-fdf3-4edf-9534-5c4f677c2889
# ╠═a39e1661-4765-411c-a062-233a64770391
# ╠═e34c35f8-1ec8-47ba-af8e-77f10cf2c27b
# ╟─11de4e14-3b56-4238-bed8-c110f4de2d44
# ╠═d619df94-5bb2-4349-a83d-ccfd13b95906
# ╠═e17bf6fe-f3b6-4904-bb05-8b068fd9cf1f
# ╠═54a2712b-d696-4097-8522-f5e1a87ecbec
# ╠═d5381c2d-1794-4af1-9ebb-5188334fc592
# ╠═610d535a-2025-4419-b35b-8d28dbaa62b8
# ╟─917e88c6-cf5d-4d8f-93e5-f50ea2bb2cdc
# ╠═e5b49869-2763-4ec9-ae7f-7b70164c0c67
# ╠═dd48df07-4bed-47ce-9799-05958e3adc7a
# ╠═47138dbf-5c60-4cdc-b484-ff168d11055f
# ╠═87a761c6-85ce-4342-a7d4-a13a378c6c45
# ╠═ba1c7f65-3e24-4b86-b2e2-521459993250
# ╠═8b0b5a9f-666a-43c2-bdab-3d771334f12b
# ╟─f02137a4-75a1-49ac-81a3-a1f04bab9ca2
# ╟─c7f1acaa-d53e-427e-b646-f3920a2ce6b7
# ╟─9a09e9cb-26cf-4576-a581-7b832fbab775
# ╠═fda7f746-0f7e-4dff-b8e6-7765f580e542
# ╠═b7be7b57-7e69-4cc9-af54-db9465f18d05
# ╠═81ccb3e6-dcca-4b66-9be6-a1fb63f7e056
# ╠═a50d4f05-0f4a-4ae0-ba5d-a27bad3869e0
# ╟─4a73771a-3adc-4fad-8d6a-148dcf9cc3c4
# ╠═baa7d024-ec90-4744-922d-830f40683abe
# ╠═1b8f3f28-e612-40d0-9b09-136543cbb126
# ╠═390a66a8-e497-4a0f-b6c4-34487df787e9
# ╠═23e0cf31-1c67-48d1-b014-26c1b44e04a8
# ╠═3eeea8dd-de5c-4c73-8d82-4bdb4979d2b0
# ╠═eaed9fdd-d4a2-41f7-8e05-588869e12780
# ╠═9634b96c-18f1-479e-b28b-3d614893ce7c
# ╠═58c8bd1d-2e69-48ac-9cfc-bb525ebe79c8
# ╠═8397c671-c0d1-4632-ae9a-55a6dccd1002
# ╠═496b86c5-900e-4786-a254-082b8155c65b
# ╠═fa394ad0-055e-4c8f-a3cf-931cd421cb7e
# ╠═897c8de6-c042-49aa-98ee-e5b3762756b8
# ╠═abc5bab1-66c6-4ae8-95a8-3ba6d7eecbf5
# ╠═1d23d9a9-996a-4c88-9bef-2fcc55fd5b44
# ╠═4d44c3bf-910d-411e-b3b8-99c1c60d43b1
# ╠═ccc85a17-690a-4fa4-9b14-9ca58a22e9c8
# ╠═46e16092-d952-4c4f-a952-c5201797fcd1
# ╠═f51eb4f5-3529-4c5a-8e2d-1902242ab7af
# ╠═6dac8047-1170-4b65-9cee-2c8db3b62d63
# ╠═e74f1990-5dd2-4062-a8d0-345a5005d0c2
# ╠═ff046879-3607-47e4-afe3-b42bb1738b9f
# ╠═ef0c1872-da54-48b5-8212-634d7d91e9ac
# ╠═8d1c0954-ace2-4621-99aa-ce692936247b
# ╠═5356fd33-c959-40b0-bdf8-d1363a0726c9
# ╠═3ebdebe8-ad66-4721-b9ff-92da125bcf7c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
