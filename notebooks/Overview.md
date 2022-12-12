# Overview of Notebooks of the Project _ThermodynamicData

Several notebooks, written in the Julia programming language (https://julialang.org/) and using the Pluto.jl package (https://github.com/fonsp/Pluto.jl), are used to
calculate the retention parameters from retention data (`Fit_lnk-T_data.jl`) or convert available sets of retention parameters into the other sets (`Convert_Parameters.jl`). These notebooks can be found in the folder `notebooks`. Functions, which are used in these notebooks can be found in the file `ThermodynamicData.jl` in the `src` folder. 

- `Convert_Parameters.jl`:
  - Load all files with the keyword `Parameters` from a folder
  - convert the parameter set from the file into the other parameter sets
  - ToDo: if no column for `CAS` number exists in the data, add one with the `CAS`-numbers from ChemicalIdentifiers.jl (also by using the file `shortname.csv`) resp. from the file `missing.csv` (chemical substances, which are not found in the database of ChemicalIdentifiers.jl)
  - optional: save the data of all three parameter sets into a file with the keyword `AllParam`

- `Fit_lnk-T_data.jl`
  - Load all files with the keyword `lnk-T` resp. `log10k-T` from a folder
  - fit of the ABC-model and the Kcentric model to the data and convert them to the third parameter set
  - export the parameter sets together with some statistics of the fits
  - ToDo: if no column for `CAS` number exists in the data, add one with the `CAS`-numbers from ChemicalIdentifiers.jl (also by using the file `shortname.csv`) resp. from the file `missing.csv` (chemical substances, which are not found in the database of ChemicalIdentifiers.jl)
  - optional: save the data of all three parameter sets into a file with the keyword `AllParam`

- `Database.jl`:
  - Load all files with keyword `AllParam`
  - combine the data of all these files in one dataframe
  - using the options:
    - filter for substances with a CAS number
    - filter out flagged substances (according to different criteria the retention parameters are problematic)
    - format (all parameter sets or only `ABC`, `Kcentric` or `TD` (thermodynamic) parameter sets)
  - additional filter can be applied for
    - stationary phases
    - dimensionless film thickness
    - source
    - categories
  - save/download option for the filtered database
- `AllParam.jl`: 
  - Load all files with keyword `AllParam` from a folder
  - combine the data of all these files in one dataframe
  - flag substances with parameters which do not fulfill certain conditions
  - plots of the parameter sets
  - use of ChemicalIdentifiers.jl to identify CAS, formula, MW and smiles
  - optional: save the dataset of the Kcentric model of all not-flagged substances into one database file, which could be used for the simulation of GC with the Julia package `GasChromatographySimulator.jl`.

- `AllParam_FitPlane.jl`: 
  - addition to `AllParam.jl`
  - fit a 2-D plane to the 3D data of (A,B,C) resp. (Tchar,thetachar,DeltaCp) with different methods