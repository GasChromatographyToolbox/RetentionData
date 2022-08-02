# Overview of Notebooks of the Project _ThermodynamicData

- `AllParam.jl`: 
  - Load all files with keyword `AllParam` from a folder
  - combine the data of all these files in one dataframe
  - flag substances with parameters which do not fulfill certain conditions
  - plots of the parameter sets
  - use of ChemicalIdentifiers.jl to identify CAS, formula, MW and smiles

- `AllParam_FitPlane.jl`: 
  - addition to `AllParam.jl`
  - fit a 2-D plane to the 3D data of (A,B,C) resp. (Tchar,thetachar,DeltaCp) with different methods

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

- `Eval_Telu_Tchar.jl`
  - Load simulation results from the folder `/data/sims/Telu_Tchar/` (results from the script `/scripts/Telu_Tchar.jl`)
  - plots of Telu over Tchar and Tinit for different heating rates and inlet pressures (heating rates amd inlet pressure can be combined in dimensionless heating rate), identification of outliers
  - definition of four models to fit the data of Telu(Tchar)
  - fits for different heating rates, inlet pressures resp. different dimensionless heating rates
    - general trends for parameters:
      - the constant temperature T1 is a linear function of the start temperature Tinit, offset is depending on rT
      - the slope is near 1 (but seems to also depend on theta_char, especially for low heating rates)
      - curvature b could be constant around 0.05?
      - intercept of increasing line depends on rT, could be independent of Tinit