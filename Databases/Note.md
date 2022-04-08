# Note

## Structure of filenames:

The general structure of filenames in this Database is the following:

Source\_Datatype\_(Substances)or(Table)\_StationaryPhase\_PhaseRatio\_(AddParameters).csv

with:
- **Source**: Name of the source, identical with the key of the Bibtex file `references.bib` (name of first author and year). For not yet publicized datasets the name consists of the name of the primary person which measured the data and the year.
- **Datatype**: Which data is in the file
  - `Parameters`: One of the following sets of thermodynamic parameters (or combinations):
    - `A, B, C`
    - `DeltaHref`, `DeltaSref`, `DeltaCp`
    - `Tchar`, `thetachar`, `DeltaCp`
  - `lnk-T`: The measured `lnk`-values (natural logarithm of the retention factor) of the substances for isothermal GC measurements at defined temperatures `T`. The thermodynamic parameters can be estimated from these measurements by fitting of the model.
  - `log10k-T`: The measured `log10(k)`-values (decadic logarithm of the retention factor) of the substances for isothermal GC measurements at defined temperatures `T`. The thermodynamic parameters can be estimated from these measurements by fitting of the model.
  - `tR-T`: The measured `tR`-values (retention time) of the substances for isothermal GC measurements at defined temperatures `T`. Also, the hold-up time `tM` must be included (measured or calculated), to calculate `lnk` from this data and estimate the thermodynamic parameters by model fitting. 
  - `AllParam`: thermodynamic parameters of all three sets are in the data, also a reference temperature `Tref`, the phase ratio $\beta$ (see below), and the argument for the Lambert W function are included 
  - `IDList`: A list of at least two columns. First column has the shortname or number (ID), the second column has a full name or alternative name. Some data (`Parameters`, `lnk-T` or `tR-T`) use the shortname/ID-number instead of a full name. This file can be used to convert the ID to full name.
- **(Substances)**: Optional. A label for the substances in the file, e.g. `Mix1` or `Alkanes`. 
- **(Table)**: Table + Number of the table in the source.
- **StationaryPhase**: The name of the stationary phase. Spaces or hyphens are removed from the name. The name is follows in general the name given by the manufacturer.
- **PhaseRatio**: The nominal phase ratio of the used column according the the data of the manufacturer. Calculated by the approximation $\beta \approx 1/4 d/d_f$
- **(AddParameters)**: Additional parameters (the value follows the parameters without spaces):
  - `Tref`: Reference temperature for the thermodynamic parameters `DeltaHref` and `DeltaSref`
  - `gas`: The type of gas used as stationary phase for the measurement of isothermal chromatograms or thermodynamic parameters, e.g. `He` or `H2`.
  - `d`: The diameter of the column used for the measurements of isothermal chromatograms or thermodynamic parameters, e.g. `0.10` or `0.32`.

## Structure of the files:

- comma separated
- first column with the names of the substance (header can be different)
- following columns depend on the **Datatype**
  - `Parameters`: three columns of the different thermodynamic parameter sets (or multiple for multiple datasets)
  - `lnk_T`: several columns with the `lnk`-values at the several temperatures `T`. Every columns stands for one temperature.
  - `tR_T`: similar to `lnk_T` plus an additional column with the hold-up time `tM`
- optional columns containing information about affiliations of the substances to different classes or categories, name with `Cat` or `Category` or `Cat_1` (important is the string "Cat") 
- first row contains the name of the columns, resp. in case of `lnk_T`or `tR_T` the name of the column is the value of the temperature
- (Optional) second row contains the units corresponding to the values in the columns
- if no row with units is present, it is assumed, that the values have their corresponding SI-units
- temperatures are assumed to be measured in °C