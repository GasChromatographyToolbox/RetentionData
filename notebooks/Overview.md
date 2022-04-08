# Overview of Notebooks of the Project _ThermodynamicData

- nb_data_survey.jl:

- nb_heating_rate.jl:

- plnb_lnk_beta.jl:

- plnb_Telu_Tchar_sim_phase_ratio_corr.jl: (**TODO: rename the notebook**)
    - Simulation of all solutes available for different stationary phases
    - variation of heating rate (temperature program constist only of on heating
      ramp, no holding times, no thermal gradients)
    - variation of inlet pressure (constant pressure)
    - the difference ``\Delta T_{elu}`` between ``T_{char}``(corrected for the phase ratio) and
      elution temperature ``T_{elu}`` is determined (**TODO: change the**
      **algorithmen for choosing the allowed solutes for this estimation, the**
      **usage of a constant mobility and only allowing a certain deviation has**
      **some problems and leads to erronous estimations of ``\Delta T_{elu}``**
      **and to jumps in the data**)
    - these differences ``\Delta T_{elu}`` are plotted against the dimensionless
      heating rate (hold-up time at 150°C as reference time) and three different
      models are fitted against the data (**TODO: derive the heating rate for**
      **which ``\Delta T_{char}=0``°C in relation to the phase ratio ``\beta``**) 


- `AllParam.jl`: 
  - Load all files with keyword `AllParam` from a folder
  - combine the data of all these files in one dataframe
  - flag substances with parameters which do not fulfill certain conditions
  - plots of the parameter sets
  - use of ChemicalIdentifiers.jl to identify CAS, formula, MW and smiles

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