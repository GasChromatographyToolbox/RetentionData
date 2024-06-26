# script for loading all lnk_T data and calculating the parameters of the thermodynamic retention models
root = dirname(@__FILE__)
project = dirname(root)
db_path = joinpath(project, "Databases")
	
using CSV, DataFrames, LambertW, Plots, LsqFit, Statistics, ChemicalIdentifiers, Measurements, RAFF
include(joinpath(project, "src", "RetentionData.jl"))

# steps from Fit_lnk-T_Outliers.jl
	# load the lnk_T data
	data = RetentionData.load_lnkT_data(db_path)

	# load the saved settings for excluding found outliers in the lnk_T data (`true`) or ignore outlier test (`false`)
	# if no file with settings exist exclude found outliers
	CheckBase = Array{Array{Bool,1}}(undef, size(data.filename)[1])
	for i=1:size(data.filename)[1]
		if isfile(joinpath(project, "OutlierCheck", string("Check_",data.filename[i])))
			CheckBase[i] = DataFrame(CSV.File(joinpath(project, "OutlierCheck", string("Check_",data.filename[i])); types=[String, Bool])).OutlierTest
		else 
			CheckBase[i] = trues(size(data.data[i])[1])	
		end
	end	

	# make the fitting (with outlier check using RAFF.jl)
	fitting = RetentionData.fit_models(data.data, CheckBase)
	data[!, "fitting"] = fitting
	# extract the fitting results in a dataframe format for saving
	RetentionData.extract_parameters_from_fit!(data)

	# save the estimated parameters:
	RetentionData.save_all_parameter_data(data; rounding=true, sigdigits=5, errors=true)

# steps from Convert_Parameters.jl
	# load the parameter data, automatic convertion of parameter sets
	data_p = RetentionData.load_parameter_data(db_path)

	# save the converted parameters:
	RetentionData.save_all_parameter_data(data_p)

# database_all
db_all = unique(RetentionData.database(db_path; filter_CAS=false, filter_flag=false, db_format="all")[1])
CSV.write(joinpath(db_path, "database_all.csv"), db_all) # without duplicated rows

# database_nonflag
db_nonflag = unique(RetentionData.database(db_path; filter_CAS=false, filter_flag=true, db_format="all")[1])
CSV.write(joinpath(db_path, "database_nonflag.csv"), db_nonflag) # without duplicated rows

# GCSim_database_nonflag
db_GCSim_nonflag = unique(RetentionData.database(db_path; filter_CAS=false, filter_flag=true, db_format="GasChromatographySimulator")[1])
CSV.write(joinpath(db_path, "GCSim_database_nonflag.csv"), db_GCSim_nonflag) # without duplicated rows

# GCSim_database_all
db_GCSim_all = unique(RetentionData.database(db_path; filter_CAS=false, filter_flag=false, db_format="GasChromatographySimulator")[1])
select!(db_GCSim_all, Not(:flag))
CSV.write(joinpath(db_path, "GCSim_database_all.csv"), db_GCSim_all)

# unique solutes in GCSim_database_nonflag
unique(db_GCSim_nonflag.CAS)
unique(db_GCSim_nonflag.Name)
# unique stationary phases
unique(db_GCSim_nonflag.Phase)
# substances without CAS
noCAS = unique(db_GCSim_nonflag.Name[findall(ismissing.(db_GCSim_nonflag.CAS).==1)])
noCAS_all = unique(db_GCSim_all.Name[findall(ismissing.(db_GCSim_all.CAS).==1)])
unique(db_GCSim_all.CAS)
unique(db_GCSim_all.Phase)