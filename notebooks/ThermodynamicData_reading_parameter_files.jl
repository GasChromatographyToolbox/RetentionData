### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 74354168-931c-11ec-26f7-0d1d90c7088b
using CSV, DataFrames

# ╔═╡ 2313df1f-6d5b-4fb0-b902-6e787b932746
using LambertW

# ╔═╡ d4ef5b75-f073-45f2-97ed-3de7c47f6732
using Plots

# ╔═╡ e395da3e-928b-48f0-9651-85cc2c5a2dd8
using ChemicalIdentifiers

# ╔═╡ 9fae926f-e423-4156-9abe-97b351fe9164
db_path = "/Users/janleppert/Documents/GitHub/ThermodynamicData/Databases"

# ╔═╡ bfd17656-0ee6-4d03-8ec0-3d6a717fc3f2
md"""
Read all files in sub-folders containing parameters and extract information from the file name:
* [x] read all objects in folder `db_path`

* [x] if object is a folder, than also read all objects in these

* [x] filter objects for filenames containing the word `Parameters` and the file extension `.csv`

* [x] read the additional information from the filename: source, stationary phase, phase ratio, optional parameters (`Tref`, `gas`, `d`)

* [ ] read the separate .csv-files (`Parameters`):

	* [x] identify the type of stored parameters and convert them to the other sets of parameters. Store these in new files in the same folder with keyword `AllParameters`.

	* [x] use ChemicalIdentifiers.jl to look up substances

	* [ ] build from all extracted informations databases (new and old format) for every .csv-file separtly AND one combined (filter out bad data, and not found substances in ChemicalIdentifiers.jl; add something for dupilcated datasets (same substance on same stationary phase))

!!! note
	Additional parameters `gas` and `d` should be attached to `source` for the exported databases

"""

# ╔═╡ fc9cbfe9-f7da-44e2-a19b-5674b6b8174a
"""
	collect_csv_paths(folder)

Collect the paths of all .csv-files in the `folder` and all its sub-folders in an array of Strings.
"""
function collect_csv_paths(folder)
	csv_obj = String[]
	
	db_obj = filter(x -> isdir.(x) || endswith.(x, ".csv"), readdir(folder, join=true))
	
	for i=1:length(db_obj)
		if endswith(db_obj[i], ".csv")
			push!(csv_obj, db_obj[i])
		elseif isdir(db_obj[i])
			subfolder_obj = collect_csv_paths(db_obj[i])
			for j=1:length(subfolder_obj)
				push!(csv_obj, subfolder_obj[j])
			end
		end
	end
	
	return csv_obj
end

# ╔═╡ bf5e6f39-70da-445e-8fa1-3685ceae6910
all_csv = collect_csv_paths(db_path)

# ╔═╡ 2fd9a066-3432-44bf-892a-34ed40ef6314
begin # filter for certain key words in the filename
	keyword = "Parameters"
	parameters_csv = filter(contains(keyword), all_csv)
end

# ╔═╡ cd66d6d6-30cc-41bc-b1d6-cf2ef9948419
"""
	extract_meta_data(csv_paths::Array{String,1})

Extract the meta data contained in the filename of a collection of paths to .csv-files. A dataframe with the columns:
- path
- filename
- source
- phase
- beta0
- Tref
- d
- gas
is returned. For the correct format of the filename see 'Note.md'.
"""
function extract_meta_data(csv_paths::Array{String,1})
	rg = r"(?<=\D)(?=\d)" # pattern: non-digits followed by digits -> pattern of beta250
	path = Array{String}(undef, length(csv_paths))
	filename = Array{String}(undef, length(csv_paths))
	source = Array{String}(undef, length(csv_paths))
	phase = Array{String}(undef, length(csv_paths))
	beta0 = Array{Float64}(undef, length(csv_paths))
	Tref = missings(Float64, length(csv_paths))
	d = missings(Float64, length(csv_paths))
	gas = missings(String, length(csv_paths))
	# split the filename from the path
	for i=1:length(csv_paths)
		path[i], filename[i] = splitdir(csv_paths[i]) # separate filename from path
		split_fn = split(split(filename[i], r".csv")[1], '_')
		source[i] = split_fn[1]
		phase[i] = split_fn[4]
		beta0[i] = parse(Float64, split(split_fn[5], rg)[2])
		for j=6:length(split_fn)
			if contains(split_fn[j], "Tref")
				Tref[i] = parse(Float64, split(split_fn[j], "Tref")[2])
			elseif contains(split_fn[j], "d")
				d[i] = parse(Float64, split(split_fn[j], "d")[2])
			elseif contains(split_fn[j], "gas")
				gas[i] = split(split_fn[j], "gas")[2]
			end
		end
	end
	meta_data = DataFrame(path=path, filename=filename, source=source, phase=phase, beta0=beta0, Tref=Tref, d=d, gas=gas)
	return meta_data
end

# ╔═╡ f3b9a58b-caaa-4413-a657-3c7a2ac5b06f
meta_data = extract_meta_data(parameters_csv)

# ╔═╡ 2175d40e-d618-4e2f-963e-446010deaba9
#="""
	load_csv_data(meta_data::DataFrame)

Load the data from the .csv-file located at `meta_data.path`, `meta_data.filename`. Data for `DeltaHref` is multiplied by 1000, if the corresponding unit is `kJ/mol`. The loaded data is returned in a dataframe.
"""
function load_csv_data(meta_data::DataFrame)
	data = Array{DataFrame}(undef, length(meta_data.filename))
	for i=1:length(meta_data.filename)
		path = joinpath(meta_data.path[i], meta_data.filename[i])
		data_ = DataFrame(CSV.File(path, header=1))
		# problem: some files have one header, some have two header lines (units in second line)
		if typeof(data_[!, 2]) == Array{Float64, 1} # no second header with units
			data[i] = data_
		else # second header with units 
			data[i] = DataFrame(CSV.File(path, header=1:2))
		end
		# if a column is DeltaHref and has a unit with kJ, than the values have to be multiplied by 1000
		col_i = findfirst(occursin.("DeltaHref", names(data[i])))
		if typeof(col_i) != Nothing
			if occursin("kJ", names(data[i])[col_i])
				data[i][!, col_i] = data[i][!, col_i].*1000.0
				rename!(data[i], names(data[i])[col_i] => "DeltaHref_J/mol")
			end
		end
	end
	return data
end=#
"""
	load_csv_data(meta_data::DataFrame)

Load the data from the .csv-file located at `meta_data.path`, `meta_data.filename`. Data for `DeltaHref` is multiplied by 1000, if the corresponding unit is `kJ/mol`. The loaded data is returned in a dataframe. If an ID-list is available in the same folder the shortname/number will be looked up in the ID-list and replaced by the full name given in the ID-List.
"""
function load_csv_data(meta_data::DataFrame)
	data = Array{DataFrame}(undef, length(meta_data.filename))
	for i=1:length(meta_data.filename)
		path = joinpath(meta_data.path[i], meta_data.filename[i])
		data_ = DataFrame(CSV.File(path, header=1, stringtype=String))
		# problem: some files have one header, some have two header lines (units in second line)
		if typeof(data_[!, 2]) == Array{Float64, 1} || typeof(data_[!, 2]) == Array{Union{Missing,Float64}, 1} # no second header with units
			data[i] = data_
		else # second header with units 
			data[i] = DataFrame(CSV.File(path, header=1:2))
		end
		# if a column is DeltaHref and has a unit with kJ, than the values have to be multiplied by 1000
		col_i = findfirst(occursin.("DeltaHref", names(data[i])))
		if typeof(col_i) != Nothing
			if occursin("kJ", names(data[i])[col_i])
				data[i][!, col_i] = data[i][!, col_i].*1000.0
				rename!(data[i], names(data[i])[col_i] => "DeltaHref_J/mol")
			end
		end
		# if the filepath has a file similar to `Source_IDList_...csv` than the names (numbers/shortnames) from the filename should be replaced with the names from the IDList (only one IDList in the same folder)
		IDfile = filter(x -> endswith.(x, ".csv") && contains.(x, "IDList"), readdir(meta_data.path[i]))
		if length(IDfile) == 1
			IDlist = DataFrame(CSV.File(joinpath(meta_data.path[i], IDfile[1]), header=1, stringtype=String))
			for j=1:length(data[i][!,1])
				j_ID = findfirst(data[i][!,1][j].==IDlist[!,1])
				if isnothing(j_ID) != true
					data[i][!,1][j] = IDlist[!,2][j_ID]
				end
			end
		end
	end
	return data
end

# ╔═╡ 5e1accb7-dbc0-43fa-96e6-657b25d3f146
data = load_csv_data(meta_data)

# ╔═╡ a5bc0bab-76f5-40a7-aa28-9643227abdb9
"""
	load_csv_data!(meta_data::DataFrame)

Load the data from the .csv-file located at `meta_data.path`, `meta_data.filename`. Data for `DeltaHref` is multiplied by 1000, if the corresponding unit is `kJ/mol`. The loaded data is returned in an additional column `data` in the `meta_data` dataframe.
"""
function load_csv_data!(meta_data::DataFrame)
	data = load_csv_data(meta_data)
	meta_data[!, "data"] = data
	return meta_data
end

# ╔═╡ ab0e8a97-bb9f-4c94-9474-cd962d43554c
load_csv_data!(meta_data)

# ╔═╡ ca41febf-48ee-4b10-8fbd-f473652eab59
"""
	identify_parameters(data::Array{DataFrame})

Extract the column indices of the parameters and identify the type of parameters sets stored in the dataframes of `data`.  

ABC-parameters:
- column name containing "A" -> A
- column name containing "B" -> B
- column name containing "C_" or column name == "C" -> C

K-centric:
- column name containing "Tchar" -> Tchar
- column name containing "thetachar" -> thetachar
- column name containing "DeltaCp" -> DeltaCp

thermodynamic parameters at reference temperature Tref
- column name containing "DeltaHref" -> DeltaHref
- column name containing "DeltaSref" -> DeltaSref
- column name containing "DeltaCp" -> DeltaCp
- additionally from meta_data a value for Tref is needed

# Output
- param_set: Array of Arrays of string, listening the types of parameter sets, possible values "ABC", "K-centric", "TD"
- index: Array of dictionaries with the colum index of the parameters (parameters name is the key, the index is the value)
"""
function indentify_parameters(data::Array{DataFrame})
	param_set = Array{Array{String,1}}(undef, length(data))
	index = Array{Dict}(undef, length(data))
	for i=1:length(data)
		# ABC model
		i_A = findfirst(occursin.("A", names(data[i])))
		i_B = findfirst(occursin.("B", names(data[i])))
		i_C = findfirst(occursin.("C_", names(data[i])))
		if typeof(i_C) == Nothing
			i_C = findfirst("C".==names(data[i]))
		end
		# K-centric model
		i_Tchar = findfirst(occursin.("Tchar", names(data[i])))
		i_thetachar = findfirst(occursin.("thetachar", names(data[i])))
		i_DeltaCp = findfirst(occursin.("DeltaCp", names(data[i])))
		# TD model
		i_DeltaHref = findfirst(occursin.("DeltaHref", names(data[i])))
		i_DeltaSref = findfirst(occursin.("DeltaSref", names(data[i])))
		# i_DeltaCp from above
		p_set = String[]
		if typeof(i_A) == Int && typeof(i_B) == Int && typeof(i_C) == Int
			push!(p_set, "ABC")
		end
		if typeof(i_Tchar) == Int && typeof(i_thetachar) == Int && typeof(i_DeltaCp) == Int
			push!(p_set, "K-centric")
		end
		if typeof(i_DeltaHref) == Int && typeof(i_DeltaSref) == Int && typeof(i_DeltaCp) == Int
			push!(p_set, "TD")
		end
		param_set[i] = p_set
		index[i] = Dict("A" => i_A, "B" => i_B, "C" => i_C,
			"Tchar" => i_Tchar, "thetachar" => i_thetachar, "DeltaCp" => i_DeltaCp,
			"DeltaHref" => i_DeltaHref, "DeltaSref" => i_DeltaSref)
	end
	return param_set, index
end

# ╔═╡ e7d29498-97b8-4458-8fa3-3690ec4ab90e
paramset, paramindex = indentify_parameters(meta_data.data)

# ╔═╡ de7c450d-f912-4715-8ae3-ec2226e024a0
paramset

# ╔═╡ ed315a81-54af-43d1-aa03-6ce6dc73f078
R = 8.31446261815324

# ╔═╡ 77303177-5328-4fc0-804e-17e91ab7ba02
Tst = 273.15

# ╔═╡ 57a9f34e-7d3b-4371-8ab4-ebc78079e597
"""
	TD_to_ABC(ΔHref, ΔSref, ΔCp, Tref)

Convert the parameters of the thermodynamic set (`ΔHref`, `ΔSref`, `ΔCp`) to the parameters of the `ABC` set (`A`, `B`, `C`). The reference temperature `Tref` is a needed additional information. The parameter should have been measured in
the following units:
- `ΔHref` in J mol⁻¹
- `ΔSref` in J mol⁻¹ K⁻¹
- `ΔCp` in J mol⁻¹ K⁻¹
- `Tref` in °C

# Output
- `A` (without unit)
- `B` in K
- `C` (without unit)
"""
function TD_to_ABC(ΔHref, ΔSref, ΔCp, Tref)
	A = (ΔSref-ΔCp*(1+log(Tref+Tst)))/R
	B = (ΔCp*(Tref+Tst) - ΔHref)/R
	C = ΔCp/R
	return A, B, C
end

# ╔═╡ 04778b8f-df60-48b4-8164-44fcf183ea5d
"""
	ABC_to_TD(A, B, C, Tref)

Convert the parameters of the `ABC` set (`A`, `B`, `C`) to the parameters of the thermodynamic set (`ΔHref`, `ΔSref`, `ΔCp`). The reference temperature `Tref` is a needed additional information. The parameter should have been measured in
the following units:
- `A` (without unit)
- `B` in K
- `C` (without unit)
- `Tref` in °C

# Output
- `ΔHref` in J mol⁻¹
- `ΔSref` in J mol⁻¹ K⁻¹
- `ΔCp` in J mol⁻¹ K⁻¹
"""
function ABC_to_TD(A, B, C, Tref)
	ΔHref = R*(C*(Tref+Tst) - B)
	ΔSref = R*(A + C + C*log(Tref+Tst))
	ΔCp = R*C
	return ΔHref, ΔSref, ΔCp
end

# ╔═╡ 1a4514a1-4c7e-4c0e-974c-81e3c19962e4
"""
	lambertw_x(A, B, C, KRef)

Calculate the argument for the Lambert W function (-1 branch) used to convert the `ABC` set to the K-centric set. The reference distribution coefficient `KRef` is a needed additional information (normally, it should be equal to the phase ratio `β` of the reference column, where the parameters where measured). The value should be -1/e ≤ x < 0.
"""
function lambertw_x(A, B, C, KRef)
	x = -B*exp(A/C)/(C*KRef^(1/C))
	return x
end

# ╔═╡ a1190a94-fc09-4727-b389-14131a6a31ab
"""
	ABC_to_Kcentric(A, B, C, KRef)

Convert the parameters of the `ABC` set (`A`, `B`, `C`) to the parameters of the K-centric set (`Tchar`, `θchar`, `ΔCp`). The reference distribution coefficient `KRef` is a needed additional information (normally, it should be equal to the phase ratio `β` of the reference column, where the parameters where measured). The parameter should have been measured in
the following units:
- `A` (without unit)
- `B` in K
- `C` (without unit)
- `KRef` (without unit)

# Output
- `Tchar` in °C
- `θchar` in °C
- `ΔCp` in J mol⁻¹ K⁻¹
"""
function ABC_to_Kcentric(A, B, C, KRef)
	x = lambertw_x(A, B, C, KRef)
	if x >= -1/exp(1) && x < 0.0
		Wx = lambertw(x,-1)
		Tchar = -B/(C*Wx) - Tst
		θchar = B/(C^2*(1 + Wx)*Wx)
	else
		Wx = lambertw(x,0)
		Tchar = -B/(C*Wx) - Tst
		θchar = B/(C^2*(1 + Wx)*Wx)
	end
	ΔCp = C*R
	return Tchar, θchar, ΔCp
end

# ╔═╡ e30f440c-d798-4092-8698-d14f512c0fed
"""
	Kcentric_to_ABC(Tchar, θchar, ΔCp, KRef)

Convert the parameters of the K-centric set (`Tchar`, `θchar`, `ΔCp`) to the parameters of the `ABC` set (`A`, `B`, `C`). The reference distribution coefficient `KRef` is a needed additional information (normally, it should be equal to the phase ratio `β` of the reference column, where the parameters where measured). The parameter should have been measured in
the following units:
- `Tchar` in °C
- `θchar` in °C
- `ΔCp` in J mol⁻¹ K⁻¹
- `KRef` (without unit)

# Output
- `A` (without unit)
- `B` in K
- `C` (without unit)
"""
function Kcentric_to_ABC(Tchar, θchar, ΔCp, KRef)
	A = log(KRef) - (Tchar+Tst)/θchar - ΔCp/R*(1 + log(Tchar+Tst))
	B = ΔCp/R*(Tchar+Tst) + (Tchar+Tst)^2/θchar
	C = ΔCp/R
	return A, B, C
end

# ╔═╡ 37f3b288-7453-44cb-803d-74ffcea945ac
meta_data

# ╔═╡ 5d57b21f-5604-416a-9eec-6df3e7014643
typeof(paramset)

# ╔═╡ 2c37624c-59f0-4ec6-84c4-c2dbcd4201b7
typeof(paramindex)

# ╔═╡ b303f4d6-eadc-4d98-921b-3b8798e2d29b
typeof(data)

# ╔═╡ 99656cd8-670b-4eb3-bf95-108abf88cfb9
ismissing(meta_data.Tref[1])

# ╔═╡ 5e7a06c3-9ce4-4080-9368-5dca65f5398c
T0 = 90.0

# ╔═╡ 1525825b-1c56-4df3-b7fb-bb80b36c8cda
data[1][!, paramindex[1]["A"]][2]

# ╔═╡ 6c679f66-0b9f-4192-8d0b-585bcaaa4b56
"""
	all_parameters(data::Array{DataFrame,1}, paramset::Array{Array{String,1},1}, paramindex::Array{Dict,1})

Calculate the other parameter sets from the given parameter sets and return them in a new array of dataframes.

# Output
A new array of dataframes. The dataframes have the following columns:
- `Name`: Name of the substance
- `A`: parameter `A` of the `ABC` set
- `B`: parameter `B` of the `ABC` set
- `C`: parameter `B` of the `ABC` set
- `Tchar`: parameter `Tchar` of the `K-centric` set
- `thetachar`: parameter `θchar` of the `K-centric` set
- `DeltaCp`: parameter `ΔCp` of the `K-centric` and `TD` set
- `DeltaHref`: parameter `ΔHref` of the `TD` set
- `DeltaSref`: parameter `ΔSref` of the `TD` set
- `Tref`: the reference temperature of the `TD` set
- `beta0`: the phase ratio β0, for which the retention was measured, `KRef=β0` of the `K-centric` set
- `lambertw_x`: the argument of the Lambert W function used in the conversion of the `ABC` set to the `K-centric` set. The value should be -1/e ≤ x < 0
"""
function all_parameters(data::Array{DataFrame,1}, paramset::Array{Array{String,1},1}, paramindex::Array{Dict,1}, β0::Array{Float64,1}, T_ref::Array{Union{Missing, Float64},1}) # calculate the not available parameter sets
	new_data = Array{DataFrame}(undef, length(data))
	for i=1:length(data)
		A = Array{Float64}(undef, length(data[i][!, 1]))
		B = Array{Float64}(undef, length(data[i][!, 1]))
		C = Array{Float64}(undef, length(data[i][!, 1]))
		Tchar = Array{Float64}(undef, length(data[i][!, 1]))
		θchar = Array{Float64}(undef, length(data[i][!, 1]))
		ΔCp = Array{Float64}(undef, length(data[i][!, 1]))
		ΔHref = Array{Float64}(undef, length(data[i][!, 1]))
		ΔSref = Array{Float64}(undef, length(data[i][!, 1]))
		#Tref = Array{Float64}(undef, length(data[i][!, 1]))
		#beta0 = Array{Float64}(undef, length(data[i][!, 1]))
		beta0 = β0[i].*ones(length(data[i][!, 1]))
		W_x = Array{Float64}(undef, length(data[i][!, 1]))
		if ismissing(T_ref[i])
			Tref = T0.*ones(length(data[i][!, 1]))
		else
			Tref = T_ref[i].*ones(length(data[i][!, 1]))
		end
		for j=1:length(data[i][!, 1])
			if "ABC" in paramset[i] && "K-centric" in paramset[i] && "TD" in paramset[i]
				A[j] = data[i][!, paramindex[i]["A"]][j]
				B[j] = data[i][!, paramindex[i]["B"]][j]
				C[j] = data[i][!, paramindex[i]["C"]][j]
				Tchar[j] = data[i][!, paramindex[i]["Tchar"]][j]
				θchar[j] = data[i][!, paramindex[i]["thetachar"]][j]
				ΔCp[j] = data[i][!, paramindex[i]["DeltaCp"]][j]
				ΔHref[j] = data[i][!, paramindex[i]["DeltaHref"]][j]
				ΔSref[j] = data[i][!, paramindex[i]["DeltaSref"]][j]
			elseif "ABC" in paramset[i] && "K-centric" in paramset[i] &&!("TD" in paramset[i])
				A[j] = data[i][!, paramindex[i]["A"]][j]
				B[j] = data[i][!, paramindex[i]["B"]][j]
				C[j] = data[i][!, paramindex[i]["C"]][j]
				Tchar[j] = data[i][!, paramindex[i]["Tchar"]][j]
				θchar[j] = data[i][!, paramindex[i]["thetachar"]][j]
				ΔCp[j] = data[i][!, paramindex[i]["DeltaCp"]][j]
				TD = ABC_to_TD(A[j], B[j], C[j], Tref[j])
				ΔHref[j] = TD[1]
				ΔSref[j] = TD[2]
			elseif "ABC" in paramset[i] && !("K-centric" in paramset[i]) && "TD" in paramset[i]
				A[j] = data[i][!, paramindex[i]["A"]][j]
				B[j] = data[i][!, paramindex[i]["B"]][j]
				C[j] = data[i][!, paramindex[i]["C"]][j]
				ΔCp[j] = data[i][!, paramindex[i]["DeltaCp"]][j]
				ΔHref[j] = data[i][!, paramindex[i]["DeltaHref"]][j]
				ΔSref[j] = data[i][!, paramindex[i]["DeltaSref"]][j]
				Kcentric = ABC_to_Kcentric(A[j], B[j], C[j], beta0[j])
				Tchar[j] = Kcentric[1]
				θchar[j] = Kcentric[2]
			elseif "ABC" in paramset[i] && !("K-centric" in paramset[i]) && !("TD" in paramset[i])
				A[j] = data[i][!, paramindex[i]["A"]][j]
				B[j] = data[i][!, paramindex[i]["B"]][j]
				C[j] = data[i][!, paramindex[i]["C"]][j]
				Kcentric = ABC_to_Kcentric(A[j], B[j], C[j], beta0[j])
				Tchar[j] = Kcentric[1]
				θchar[j] = Kcentric[2]
				ΔCp[j] = Kcentric[3]
				TD = ABC_to_TD(A[j], B[j], C[j], Tref[j])
				ΔHref[j] = TD[1]
				ΔSref[j] = TD[2]
			elseif !("ABC" in paramset[i]) && "K-centric" in paramset[i] && "TD" in paramset[i]
				Tchar[j] = data[i][!, paramindex[i]["Tchar"]][j]
				θchar[j] = data[i][!, paramindex[i]["thetachar"]][j]
				ΔCp[j] = data[i][!, paramindex[i]["DeltaCp"]][j]
				ΔHref[j] = data[i][!, paramindex[i]["DeltaHref"]][j]
				ΔSref[j] = data[i][!, paramindex[i]["DeltaSref"]][j]
				ABC = Kcentric_to_ABC(Tchar[j], θchar[j], ΔCp[j], beta0[j])
				A[j] = ABC[1]
				B[j] = ABC[2]
				C[j] = ABC[3]
			elseif !("ABC" in paramset[i]) && "K-centric" in paramset[i] && !("TD" in paramset[i])
				Tchar[j] = data[i][!, paramindex[i]["Tchar"]][j]
				θchar[j] = data[i][!, paramindex[i]["thetachar"]][j]
				ΔCp[j] = data[i][!, paramindex[i]["DeltaCp"]][j]
				ABC = Kcentric_to_ABC(Tchar[j], θchar[j], ΔCp[j], beta0[j])
				A[j] = ABC[1]
				B[j] = ABC[2]
				C[j] = ABC[3]
				TD = ABC_to_TD(A[j], B[j], C[j], Tref[j])
				ΔHref[j] = TD[1]
				ΔSref[j] = TD[2]
			elseif !("ABC" in paramset[i]) && !("K-centric" in paramset[i]) && "TD" in paramset[i]
				ΔCp[j] = data[i][!, paramindex[i]["DeltaCp"]][j]
				ΔHref[j] = data[i][!, paramindex[i]["DeltaHref"]][j]
				ΔSref[j] = data[i][!, paramindex[i]["DeltaSref"]][j]
				ABC = TD_to_ABC(ΔHref[j], ΔSref[j], ΔCp[j], Tref[j])
				A[j] = ABC[1]
				B[j] = ABC[2]
				C[j] = ABC[3]
				Kcentric = ABC_to_Kcentric(A[j], B[j], C[j], beta0[j])
				Tchar[j] = Kcentric[1]
				θchar[j] = Kcentric[2]
			end
			W_x[j] = lambertw_x(A[j], B[j], C[j], beta0[j])
		end
		new_data[i] = DataFrame(Name=data[i][!, 1], A=A, B=B, C=C, Tchar=Tchar, thetachar=θchar, DeltaCp=ΔCp, DeltaHref=ΔHref, DeltaSref=ΔSref, Tref=Tref, beta0=beta0, lambertw_x=W_x)
	end
	return new_data
end

# ╔═╡ 8c1fcf33-091a-43a2-ac6f-23991bf733e2
"""
	all_parameters!(meta_data)

Calculate the other parameter sets from the given parameter sets and replace the `data` in meta_data with the following:

A new array of dataframes. The dataframes have the following columns:
- `Name`: Name of the substance
- `A`: parameter `A` of the `ABC` set
- `B`: parameter `B` of the `ABC` set
- `C`: parameter `B` of the `ABC` set
- `Tchar`: parameter `Tchar` of the `K-centric` set
- `thetachar`: parameter `θchar` of the `K-centric` set
- `DeltaCp`: parameter `ΔCp` of the `K-centric` and `TD` set
- `DeltaHref`: parameter `ΔHref` of the `TD` set
- `DeltaSref`: parameter `ΔSref` of the `TD` set
- `Tref`: the reference temperature of the `TD` set
- `beta0`: the phase ratio β0, for which the retention was measured, `KRef=β0` of the `K-centric` set
- `lambertw_x`: the argument of the Lambert W function used in the conversion of the `ABC` set to the `K-centric` set. The value should be -1/e ≤ x < 0
"""
function all_parameters!(meta_data)
	paramset, paramindex = indentify_parameters(meta_data.data)
	new_data = all_parameters(meta_data.data, paramset, paramindex, meta_data.beta0, meta_data.Tref)
	meta_data[!, "data"] = new_data
	return meta_data
end

# ╔═╡ 4d9d44c7-c7d8-4375-b8d8-429930cea772
all_parameters!(meta_data)

# ╔═╡ 8fda6700-4765-4118-a9f5-cb390aaf7736
meta_data.filename[30:33] # problematic data, also some data in 34 and 35

# ╔═╡ 6ce59875-9209-4bd7-a042-7dbc04a5d8f2
# combine all data into one big dataframe

# ╔═╡ bfb32370-3e19-4534-b809-3bb8279a196f
begin
	plotly()
	pABC = scatter(meta_data.data[1].A, meta_data.data[1].B, meta_data.data[1].C, label=1, xlabel="A", ylabel="B", zlabel="C")
	for i=2:length(meta_data.data)
		scatter!(pABC, meta_data.data[i].A, meta_data.data[i].B, meta_data.data[i].C, label=i)
	end
	pABC
end

# ╔═╡ fee366a1-8374-474d-86db-e987399b97ef
begin
	plotly()
	pKcentric = scatter(meta_data.data[1].Tchar, meta_data.data[1].thetachar, meta_data.data[1].DeltaCp, label=1, xlabel="Tchar", ylabel="thetachar", zlabel="DeltaCp")
	for i=2:length(meta_data.data)
		scatter!(pKcentric, meta_data.data[i].Tchar, meta_data.data[i].thetachar, meta_data.data[i].DeltaCp, label=i)
	end
	pKcentric
end

# ╔═╡ 196fcd36-432f-4e8d-b00a-33ad788d328e
begin
	plotly()
	pTD = scatter(meta_data.data[1].DeltaHref, meta_data.data[1].DeltaSref, meta_data.data[1].DeltaCp, label=1, xlabel="DeltaHref", ylabel="DeltaSref", zlabel="DeltaCp")
	for i=2:length(meta_data.data)
		scatter!(pTD, meta_data.data[i].DeltaHref, meta_data.data[i].DeltaSref, meta_data.data[i].DeltaCp, label=i)
	end
	pTD
end

# ╔═╡ 779f25cd-d107-45b5-9cc3-4177d0f61339
names(meta_data)

# ╔═╡ 0133ad3a-7d7b-47bd-9299-5cd5a39e8893
"""
	load_parameter_data(db_path)

Load the data files (.csv format) with the `Parameters` data from the folder `db_path` including all subfolders. Based on the loaded parameters, the parameters of the not included parameter sets are calculated. Additional information from the filenames are also saved.

# Output
A dataframes with the following columns:
- `path`: Path of the folder from where the data was loaded.
- `filename`: Name of the file from where the data was loaded.
- `source`: Short name for the source from where the data original is taken.
- `phase`: Name of the stationary phase corresponding to the data.
- `beta0`: The phase ratio corresponding to the data.
- `Tref`: The reference temperature used for the thermodynamic parameter set. Optional parameter, if not available it has the value `missing`.
- `d`: The column diameter. Optional parameter, if not available it has the value `missing`.
- `gas`: The used gas for the mobile phase. Optional parameter, if not available it has the value `missing`.
- `data`: Dataframes with the parameters of the three different parameter sets. See function all_parameters().

# Note
For the naming convention of the filenames see Note.md.
"""
function load_parameter_data(db_path)
	all_csv = collect_csv_paths(db_path)
	keyword = "Parameters"
	parameters_csv = filter(contains(keyword), all_csv)
	meta_data = extract_meta_data(parameters_csv)
	load_csv_data!(meta_data)
	all_parameters!(meta_data)
	return meta_data
end

# ╔═╡ 549c8d8a-f6f3-4771-bc68-f0a8d1d3cca8
md = load_parameter_data(db_path)

# ╔═╡ 1db01b46-98a3-40be-9ad8-96eb878cc1c4
md.data[7]

# ╔═╡ 43a09a5d-111a-4378-8791-163012ec83ce
# export the separate all_parameter data to the original folders, use keyword "AllParam", add stationary phase in the filename

# ╔═╡ 7142ad38-5fc7-435f-b2d9-fece141235b1
# Filename: Source_AllParam_statPhase.csv

# ╔═╡ 4335b7cb-65dc-467c-b957-f6579da3bf6d
"""
	save_all_parameter_data(meta_data::DataFrame)

Save the `data` in the `meta_data` dataframe in new .csv-files in the same folder as the original data using the new filename `Source_AllParam_Tablename_statPhase(_d)(_gas).csv`. `Source` is the name of the source of the original data, `statPhase` is the stationary phase corresponding to the data and `Tablename` the name of the table of the original data. The optional entrys `d` and `gas` stand for the column diameter and the gas of the mobile phase. 
"""
function save_all_parameter_data(meta_data::DataFrame)
	for i=1:length(meta_data.data)
		tablename = split(meta_data.filename[i], "_")[3]
		if ismissing(meta_data.d[i]) && ismissing(meta_data.gas[i])
			new_filename = string(meta_data.source[i], "_AllParam_", tablename, "_", meta_data.phase[i], ".csv")
		elseif ismissing(meta_data.d[i])
			new_filename = string(meta_data.source[i], "_AllParam_", tablename, "_", meta_data.phase[i], "_gas", meta_data.gas[i], ".csv")
		elseif ismissing(meta_data.gas[i])
			new_filename = string(meta_data.source[i], "_AllParam_", tablename, "_", meta_data.phase[i], "_d", meta_data.d[i], ".csv")
		else
			new_filename = string(meta_data.source[i], "_AllParam_", tablename, "_", meta_data.phase[i], "_d", meta_data.d[i], "_gas", meta_data.gas[i], ".csv")
		end
		# add a round for significant digits
		CSV.write(joinpath(meta_data.path[i], new_filename), meta_data.data[i])
	end
	
end

# ╔═╡ 3b43115f-1836-4403-b4cb-ed74f6c8ba11
names(meta_data)

# ╔═╡ f2368c38-aa15-4c71-8634-509e92544550
names(meta_data.data[1])

# ╔═╡ 4ce0159b-56d6-40d9-9c25-71a0121cb4d1
length(meta_data.data)

# ╔═╡ 82e1a3f8-71b9-4d57-a6e2-919375de0a67
length(meta_data.data[1].Name)

# ╔═╡ 4d9cf1bc-5903-42e8-82cf-5bee2dbd75e5
"""
	datafram_of_all(meta_data::DataFrame)

Combine the separate dataframes with the parameter sets of the different entrys of the meta_data dataframe into one big dataframe.
"""
function dataframe_of_all(meta_data)
	# number of all data entrys
	Nall = 0
	for i=1:length(meta_data.data)
		Nall = Nall + length(meta_data.data[i].Name)
	end

	Name = String[]
	Phase = String[]
	Source = String[]
	A = Float64[]
	B = Float64[]
	C = Float64[]
	Tchar = Float64[]
	thetachar = Float64[]
	DeltaCp = Float64[]
	DeltaHref = Float64[]
	DeltaSref = Float64[]
	Tref = Float64[]
	beta0 = Float64[]
	lambertw_x = Float64[]
	d = Any[]
	gas = Any[]

	for i=1:length(meta_data.data)
		for j=1:length(meta_data.data[i].Name)
			push!(Name, meta_data.data[i].Name[j])
			push!(Phase, meta_data.phase[i])
			push!(Source, meta_data.source[i])
			push!(A, meta_data.data[i].A[j])
			push!(B, meta_data.data[i].B[j])
			push!(C, meta_data.data[i].C[j])
			push!(Tchar, meta_data.data[i].Tchar[j])
			push!(thetachar, meta_data.data[i].thetachar[j])
			push!(DeltaCp, meta_data.data[i].DeltaCp[j])
			push!(DeltaHref, meta_data.data[i].DeltaHref[j])
			push!(DeltaSref, meta_data.data[i].DeltaSref[j])
			push!(Tref, meta_data.data[i].Tref[j])
			push!(beta0, meta_data.data[i].beta0[j])
			push!(lambertw_x, meta_data.data[i].lambertw_x[j])
			push!(d, meta_data.d[i])
			push!(gas, meta_data.gas[i])
		end
	end
	dfall = DataFrame(Name=Name, Phase=Phase, Source=Source,
						A=A, B=B, C=C,
						Tchar=Tchar, thetachar=thetachar, DeltaCp=DeltaCp,
						DeltaHref=DeltaHref, DeltaSref=DeltaSref, Tref=Tref,
						beta0=beta0, lambertw_x=lambertw_x, d=d, gas=gas)
	return dfall
end

# ╔═╡ d39b95ae-cd18-4eea-9cce-c27ecf85f6ed
dataframe_of_all(meta_data)

# ╔═╡ fe625929-7f02-4e1f-b91c-f518bdc7dc7f
-1/exp(1)

# ╔═╡ fac1ea92-ebad-43c3-a945-e68d0e727a7c
"""
	flagged_data(alldata::DataFrame)

Flag the substance data for which a certain criteria of the parameters is not fullfilled, by filter `alldata` for these substances. The criterias are:
- value of `lambertw_x` is < -1/e or > 0
- value of `DeltaCp` < 0
"""
function flagged_data(alldata::DataFrame)
	# 1st flag reason: not -1/e ≤ lambertw_x < 0
	df = filter([:lambertw_x] => x -> x < -1/exp(1) || x > 0.0, alldata)
	# additional flag reasons
	df1 = filter([:DeltaCp] => x -> x < 0.0, alldata)
	# combine the filtered dataframes and delete duplicates
	df_final = unique(vcat(df, df1))
	return df_final
end

# ╔═╡ 579cb2f2-0970-48d5-a846-5115cf1c307c
flagged_data(dataframe_of_all(meta_data))

# ╔═╡ 059c0d49-748e-49eb-a334-acf93eefc3a2
ccii = search_chemical("n-tricosane")

# ╔═╡ 9d563f19-2262-4c4d-b094-7583bf50d3ca
typeof(ccii)

# ╔═╡ 039ce767-0be1-4326-bd7e-842bc2e75de3
ismissing(ccii)

# ╔═╡ 914913e3-b174-4963-8e74-ab700b6f606b
alldata = dataframe_of_all(meta_data);

# ╔═╡ 9658b138-5c09-440e-861c-42e9f016d032
"""
	substance_identification(data::DataFrame)

Look up the substance name from the `data` dataframe with ChemicalIdentifiers.jl to find the `CAS`-number, the `formula`, the molecular weight `MW` and the `smiles`-identifier. If the name is not found in the database of ChemicalIdentifiers.jl a list with alternative names (`shortnames.csv`) is used. If there are still no matches, `missing` is used.
"""
function substance_identification(data::DataFrame)
	shortnames = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/ThermodynamicData/data/shortnames.csv"))
	
	CAS = Array{Union{Missing,String}}(missing, length(data.Name))
	formula = Array{Union{Missing,String}}(missing, length(data.Name))
	MW = Array{Union{Missing,Float64}}(missing, length(data.Name))
	smiles = Array{Union{Missing,String}}(missing, length(data.Name))
	for i=1:length(data.Name)
		if data.Name[i] in shortnames.shortname
			j = findfirst(data.Name[i].==shortnames.shortname)
			ci = search_chemical(string(shortnames.name[j]))
		else
			ci = search_chemical(data.Name[i])
		end
		if ismissing(ci)
			CAS[i] = missing
			formula[i] = missing
			MW[i] = missing
			smiles[i] = missing
		else
			CAS[i] = string(ci.CAS[1], "-", ci.CAS[2], "-", ci.CAS[3])
			formula[i] = ci.formula
			MW[i] = ci.MW
			smiles[i] = ci.smiles
		end
	end
	id = DataFrame(Name=data.Name, CAS=CAS, formula=formula, MW=MW, smiles=smiles)
	#no_id = filter([:CAS] => x -> ismissing(x), id)
	return id
end

# ╔═╡ a835efcd-7ed2-4e25-a050-36591afe2fd4
id = substance_identification(alldata)

# ╔═╡ 0ef7ac98-b7ad-47e6-af7c-3a6174d62eaa
#not found substances:
filter([:CAS] => x -> ismissing(x), id)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ChemicalIdentifiers = "fa4ea961-1416-484e-bda2-883ee1634ba5"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LambertW = "984bce1d-4616-540c-a9ee-88d1112d94c9"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
CSV = "~0.10.2"
ChemicalIdentifiers = "~0.1.5"
DataFrames = "~1.3.2"
LambertW = "~0.4.5"
Plots = "~1.25.11"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Arrow]]
deps = ["ArrowTypes", "BitIntegers", "CodecLz4", "CodecZstd", "DataAPI", "Dates", "Mmap", "PooledArrays", "SentinelArrays", "Tables", "TimeZones", "UUIDs"]
git-tree-sha1 = "d4a35c773dd7b07ddeeba36f3520aefe517a70f2"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "2.2.0"

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

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "7dd38532a1115a215de51775f9891f0f3e1bac6a"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.1"

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

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

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

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

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

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

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

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

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
git-tree-sha1 = "37c1631cb3cc36a535105e6d5557864c82cd8c2b"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.0"

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

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "6354dfaf95d398a1a70e0b28238321d5d17b2530"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.0"

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

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

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
git-tree-sha1 = "0f1017f68dc25f1a0cb99f4988f78fe4f2e7955f"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.7.1"

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
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

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
# ╠═74354168-931c-11ec-26f7-0d1d90c7088b
# ╠═9fae926f-e423-4156-9abe-97b351fe9164
# ╠═bfd17656-0ee6-4d03-8ec0-3d6a717fc3f2
# ╠═fc9cbfe9-f7da-44e2-a19b-5674b6b8174a
# ╠═bf5e6f39-70da-445e-8fa1-3685ceae6910
# ╠═2fd9a066-3432-44bf-892a-34ed40ef6314
# ╠═cd66d6d6-30cc-41bc-b1d6-cf2ef9948419
# ╠═f3b9a58b-caaa-4413-a657-3c7a2ac5b06f
# ╠═2175d40e-d618-4e2f-963e-446010deaba9
# ╠═5e1accb7-dbc0-43fa-96e6-657b25d3f146
# ╠═a5bc0bab-76f5-40a7-aa28-9643227abdb9
# ╠═ab0e8a97-bb9f-4c94-9474-cd962d43554c
# ╠═ca41febf-48ee-4b10-8fbd-f473652eab59
# ╠═e7d29498-97b8-4458-8fa3-3690ec4ab90e
# ╠═de7c450d-f912-4715-8ae3-ec2226e024a0
# ╠═57a9f34e-7d3b-4371-8ab4-ebc78079e597
# ╠═ed315a81-54af-43d1-aa03-6ce6dc73f078
# ╠═77303177-5328-4fc0-804e-17e91ab7ba02
# ╠═04778b8f-df60-48b4-8164-44fcf183ea5d
# ╠═2313df1f-6d5b-4fb0-b902-6e787b932746
# ╠═1a4514a1-4c7e-4c0e-974c-81e3c19962e4
# ╠═a1190a94-fc09-4727-b389-14131a6a31ab
# ╠═e30f440c-d798-4092-8698-d14f512c0fed
# ╠═37f3b288-7453-44cb-803d-74ffcea945ac
# ╠═5d57b21f-5604-416a-9eec-6df3e7014643
# ╠═2c37624c-59f0-4ec6-84c4-c2dbcd4201b7
# ╠═b303f4d6-eadc-4d98-921b-3b8798e2d29b
# ╠═99656cd8-670b-4eb3-bf95-108abf88cfb9
# ╠═5e7a06c3-9ce4-4080-9368-5dca65f5398c
# ╠═1525825b-1c56-4df3-b7fb-bb80b36c8cda
# ╠═6c679f66-0b9f-4192-8d0b-585bcaaa4b56
# ╠═8c1fcf33-091a-43a2-ac6f-23991bf733e2
# ╠═4d9d44c7-c7d8-4375-b8d8-429930cea772
# ╠═8fda6700-4765-4118-a9f5-cb390aaf7736
# ╠═6ce59875-9209-4bd7-a042-7dbc04a5d8f2
# ╠═d4ef5b75-f073-45f2-97ed-3de7c47f6732
# ╠═bfb32370-3e19-4534-b809-3bb8279a196f
# ╠═fee366a1-8374-474d-86db-e987399b97ef
# ╠═196fcd36-432f-4e8d-b00a-33ad788d328e
# ╠═779f25cd-d107-45b5-9cc3-4177d0f61339
# ╠═0133ad3a-7d7b-47bd-9299-5cd5a39e8893
# ╠═549c8d8a-f6f3-4771-bc68-f0a8d1d3cca8
# ╠═1db01b46-98a3-40be-9ad8-96eb878cc1c4
# ╠═43a09a5d-111a-4378-8791-163012ec83ce
# ╠═7142ad38-5fc7-435f-b2d9-fece141235b1
# ╠═4335b7cb-65dc-467c-b957-f6579da3bf6d
# ╠═3b43115f-1836-4403-b4cb-ed74f6c8ba11
# ╠═f2368c38-aa15-4c71-8634-509e92544550
# ╠═4ce0159b-56d6-40d9-9c25-71a0121cb4d1
# ╠═82e1a3f8-71b9-4d57-a6e2-919375de0a67
# ╠═4d9cf1bc-5903-42e8-82cf-5bee2dbd75e5
# ╠═d39b95ae-cd18-4eea-9cce-c27ecf85f6ed
# ╠═fe625929-7f02-4e1f-b91c-f518bdc7dc7f
# ╟─fac1ea92-ebad-43c3-a945-e68d0e727a7c
# ╠═579cb2f2-0970-48d5-a846-5115cf1c307c
# ╠═e395da3e-928b-48f0-9651-85cc2c5a2dd8
# ╠═059c0d49-748e-49eb-a334-acf93eefc3a2
# ╠═9d563f19-2262-4c4d-b094-7583bf50d3ca
# ╠═039ce767-0be1-4326-bd7e-842bc2e75de3
# ╠═914913e3-b174-4963-8e74-ab700b6f606b
# ╠═9658b138-5c09-440e-861c-42e9f016d032
# ╠═a835efcd-7ed2-4e25-a050-36591afe2fd4
# ╠═0ef7ac98-b7ad-47e6-af7c-3a6174d62eaa
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
