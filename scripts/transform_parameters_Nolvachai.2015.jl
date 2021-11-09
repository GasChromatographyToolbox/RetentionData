# transform from the parameters (A, B, C) to (Tchar, θchar, ΔCp) under knowledge of β 
using DrWatson
@quickactivate "VGGC"
#DrWatson.greet()
using CSV, DataFrames
using LambertW
using Plots
using LsqFit

# first database

#db_file = "Nolvachai.2015-SuppMat-SLB-IL111-ID010mm.csv"
db_file = "Nolvachai.2015-SuppMat-SLB-IL111-ID025mm.csv"
db_path = string(datadir("exp_pro", "Databases"),"/",db_file)
db = DataFrame(CSV.File(db_path, header=1, silencewarnings=true))
#β = 250 # to do: exact calc from d and df
d = 0.25e-3
df = 0.20e-6
β = (d-2*df)^2/(d^2-(d-2*df)^2)
φ₀ = 1/(4*β)# 0.0008
R = 8.3145

# according to Table S11
# beta(T) = -9.377e-5*T^2+8.559e-2*T-1.858e+1 # 0.1mm ID
beta(T) = -1.252e-4*T^2+1.147e-1*T-2.519e+1 # 0.25mm ID
# !!! several orders of magnitude lower than β (0.5 -> 311)
# according to Nolvachai.2015-SuppMat p. S14, eq. (S11)
lnk(T,p,β) = p[1] - p[2]/T + p[3]*log(T) -log(β)

# plot lnk with the data
T = (100:10:230).+273.15
lnk_array = Array{Float64}(undef, length(T), size(db)[1])
lnk_array_βconst = Array{Float64}(undef, length(T), size(db)[1])
for j=1:size(db)[1]
    for i=1:length(T)
        lnk_array[i,j] = lnk(T[i], [db.A[j], db.B[j], db.C[j]], beta(T[i]))
        lnk_array_βconst[i,j] = lnk(T[i], [db.A[j], db.B[j], db.C[j]], β)
    end
end
plot(T, lnk_array)
plot(T, lnk_array_βconst) # these values would not make much sense

#=
x = Array{Float64}(undef, size(db)[1])
W = Array{Float64}(undef, size(db)[1])
Tchar = Array{Float64}(undef, size(db)[1])
θchar = Array{Float64}(undef, size(db)[1])
ΔCp = Array{Float64}(undef, size(db)[1])
for i=1:size(db)[1]
    x[i] = -db.B[i]/(db.C[i]*β^(1/db.C[i]))*exp(db.A[i]/db.C[i])
    if x[i] > -1/exp(1) && x[i] < 0
        W[i] = lambertw(x[i],-1)
    else
        W[i] = NaN
    end
    Tchar[i] = -db.B[i]/(db.C[i]*W[i])
    θchar[i] = db.B[i]/(db.C[i]^2*(1+W[i])*W[i])
    ΔCp[i] = db.C[i]*R
    println("i=$(i), x=$(x[i])")
    # Error for i=24, C12 FAME, x=-1.409 -> DomainError for lambertw
    # values <-1/e (-0.367879) are not possible
end

scatter(Tchar, θchar)
scatter(Tchar, ΔCp)
scatter(db.B, db.A)
scatter(db.B, db.C)

T=300:10:600
lnk_ABC = Array{Float64}(undef, length(T), size(db)[1])
lnk_Kcent = Array{Float64}(undef, length(T), size(db)[1])
for j=1:size(db)[1]
    for i=1:length(T)
        lnk_ABC[i,j] = db.A[j] + db.B[j]/T[i] + db.C[j]*log(T[i]) - log(β)
        lnk_Kcent[i,j] = (ΔCp[j]/R + Tchar[j]/θchar[j])*(Tchar[j]/T[i]-1) + ΔCp[j]/R*log(T[i]/Tchar[j])
    end
end
ii = 24
scatter(T, lnk_ABC[:,ii])
plot!(T, lnk_Kcent[:,ii])
=#
# test the calculation of the K-centric parameters
@. model(x,p) = (p[3]/R + p[1]/p[2])*(p[1]/x-1) + p[3]/R*log(x/p[1])
Tchar_fit = Array{Float64}(undef, size(db)[1])
θchar_fit = Array{Float64}(undef, size(db)[1])
ΔCp_fit = Array{Float64}(undef, size(db)[1])
for i=1:size(db)[1]
    p0 = [160+273.15, 30.0, 100.0]
    fit = curve_fit(model, T, lnk_array[:,i], p0)
    Tchar_fit[i] = fit.param[1]
    θchar_fit[i] = fit.param[2]
    ΔCp_fit[i] = fit.param[3]
end

lnk_Kcent_fit = Array{Float64}(undef, length(T), size(db)[1])
res = Array{Float64}(undef, length(T), size(db)[1])
for j=1:size(db)[1]
    for i=1:length(T)
        lnk_Kcent_fit[i,j] = (ΔCp_fit[j]/R + Tchar_fit[j]/θchar_fit[j])*(Tchar_fit[j]/T[i]-1) + ΔCp_fit[j]/R*log(T[i]/Tchar_fit[j])
        res[i,j] = lnk_array[i,j] - lnk_Kcent_fit[i,j]
    end
end
scatter!(T, lnk_Kcent_fit, shape=:x)
# plot of the residuen between the models
plot(T, res)


db[!,:Tchar] = Tchar_fit.-273.15
db[!,:thetachar] = θchar_fit
db[!,:DeltaCp] = ΔCp_fit
db[!,:phi0] = φ₀.*ones(length(Tchar_fit))
#rename!(db, :Column5 => :Tchar, :Column6 => :thetachar, :Column7 => :DeltaCp)

# export converted parameters
export_db = false#true
if export_db==true
    CSV.write(datadir("exp_pro", "Databases", string("Convert_Nolvachai.2015_025.csv")), db)
end
