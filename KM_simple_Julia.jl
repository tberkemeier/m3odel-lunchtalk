## KM-MEMO Generated Model Code ##
# Developed by: Berkemeier Group
# Multiphase Chemistry Department
# Max Planck Institute for Chemistry, Mainz

using DifferentialEquations

# define a custom struct for ODE params
struct ODEParams

    O3_a0::Float64
    O3_Hcp::Float64
    OL_Db::Float64
    O3_Td::Float64
    O3_sigma::Float64
    O3_KHcc::Float64
    O3_MTV::Float64
    A::Vector{Float64}
    V::Vector{Float64}
    kg1::Vector{Float64}
    kd1::Vector{Float64}
    kbs1::Vector{Float64}
    kbs2::Vector{Float64}
    ksb2::Vector{Float64}
    kb::Matrix{Float64}
    kbr::Float64
    c1::Float64
    c2::Float64

end

function OLEIC_F1(dydt, y, ode_params, t)

    ### Unpack ODEparams ###

    O3_a0 = ode_params.O3_a0
    O3_Hcp = ode_params.O3_Hcp
    OL_Db = ode_params.OL_Db
    OL_Db = ode_params.OL_Db
    O3_Hcp = ode_params.O3_Hcp
    O3_Td = ode_params.O3_Td
    O3_sigma = ode_params.O3_sigma
    O3_KHcc = ode_params.O3_KHcc
    O3_MTV = ode_params.O3_MTV
    A = ode_params.A
    V = ode_params.V
    kg1 = ode_params.kg1
    kd1 = ode_params.kd1
    kbs1 = ode_params.kbs1
    kbs2 = ode_params.kbs2
    ksb2 = ode_params.ksb2
    kb = ode_params.kb
    kbr = ode_params.kbr
    c1 = ode_params.c1
    c2 = ode_params.c2

    ### ODEs ###

    dydt[1] = 0
    dydt[2] = kg1[1] * y[1] * A[1] / V[2] - kg1[1] * y[2] * A[1] / V[2] - O3_a0 * O3_MTV / 4 * (1 - (O3_sigma * y[3])) * y[2] * A[2] / V[2] + kd1[1] * y[3] * A[2] / V[2]
    dydt[3] = -c1 * 1e-12 * y[3] * y[4] + O3_a0 * O3_MTV / 4 * (1 - (O3_sigma * y[3])) * y[2] - kd1[1] * y[3] - O3_KHcc * kbs1[1] / O3_Td / O3_MTV / (O3_a0 * (1 - (+O3_sigma * y[3]))) * 4 * y[3] + kbs1[1] * y[5]
    dydt[4] = -c1 * 1e-12 * y[3] * y[4] - ksb2[2] * y[4] + kbs2[2] * y[6]
    dydt[5] = -c2 * 1e-15 * kbr * y[5] * y[6] + O3_KHcc * kbs1[1] / O3_Td / O3_MTV / (O3_a0 * (1 - (+O3_sigma * y[3]))) * 4 * y[3] * (A[3] + A[4]) / 2 / V[5] - kbs1[1] * y[5] * (A[3] + A[4]) / 2 / V[5] - kb[1, 1] * y[5] * A[5] / V[5] + kb[1, 1] * y[7] * A[5] / V[5]
    dydt[6] = -c2 * 1e-15 * kbr * y[5] * y[6] + ksb2[2] * y[4] * A[4] / V[5] - kbs2[2] * y[6] * A[4] / V[5] - kb[2, 1] * y[6] * A[5] / V[5] + kb[2, 1] * y[8] * A[5] / V[5]
    dydt[7] = -c2 * 1e-15 * kbr * y[7] * y[8] + kb[1, 1] * y[5] * A[5] / V[6] - kb[1, 1] * y[7] * A[5] / V[6]
    dydt[8] = -c2 * 1e-15 * kbr * y[7] * y[8] + kb[2, 1] * y[6] * A[5] / V[6] - kb[2, 1] * y[8] * A[5] / V[6]

end

function simple(inpt::Vector{Float64})

    #inpt::Vector{Float64} = [10, 1, 1E-05, 1.9E-07, 4.8E-04, 8.09E-08, 4.63e-08, 0.1, 1.52E-15, 0.00225, 1, 1.5e3, 1]

    subset::Int = 1

    ### OPT ###
    bulkL::Int = 2
    opt_T::Float64 = 296.0
    opt_start::Float64 = -4.0
    opt_stop::Float64 = 4.8603
    abstolvar::Float64 = 0.1
    reltolvar::Float64 = 0.001

    ### DATABASE ###

    ## O3 ##
    O3_M::Float64 = 48.0
    O3_Dg::Float64 = 0.14
    O3_delta::Float64 = 3.89e-08

    ## OL ##
    OL_delta::Float64 = 8.09e-08

    ### INPUT ###

    O3_kbsxfac::Float64 = 10
    OL_kssbyfac::Float64 = 1
    O3_Db::Float64 = 1e-5
    OL_Db::Float64 = inpt[1]
    O3_Hcp::Float64 = inpt[2]
    OL_Ksurface_statbulk::Float64 = 8.09E-08
    O3_Td::Float64 = inpt[3]
    O3_a0::Float64 = 0.1
    O3_sigma::Float64 = 1.52E-15
    rout::Float64 = inpt[5]
    kbr::Float64 = 1
    c1::Float64 = inpt[4]
    c2::Float64 = 1

    ### INITIAL ###

    y0_raw::Matrix{Float64} = zeros(8, 1)
    y0_raw[1, 1] = (inpt[6])
    y0_raw[2, 1] = (inpt[6])
    y0_raw[4, 1] = (1 / (1 / (OL_Ksurface_statbulk * 1.89E+21) + (((OL_delta^3)^(1 / 3))^2)))
    y0_raw[6, 1] = (1.89E+21)
    y0_raw[8, 1] = (1.89E+21)
    y0::Vector{Float64} = y0_raw[:, subset]

    ### PRE PROCESSING ###

    R::Float64 = 8.314e7 #gas constant [g cm2 s-2 K-1 mol-1]
    R2::Float64 = 82.0578
    O3_KHcc::Float64 = O3_Hcp * R2 * opt_T #Henry"s law coefficient of X in Y []
    O3_MTV::Float64 = sqrt(8 * R * opt_T / pi / O3_M)#thermal velocity [cm s-1]
    mfp::Float64 = 1.7 * O3_Dg / O3_MTV     #mean free path in air [cm]

    ##################GEOMETRY#######################

    #LAYER SPACINGS
    Dcoarse::Float64 = rout - OL_delta
    Lsize::Float64 = Dcoarse / bulkL

    #Ltot is the total number of layers in the model, including gas, surface and bulk
    Ltot::Int = 4 + bulkL

    #SURFACE AREA AND VOLUME
    A::Vector{Float64} = zeros(Ltot - 1)
    V::Vector{Float64} = zeros(Ltot)
    D::Vector{Float64} = zeros(Ltot)
    r::Vector{Float64} = zeros(bulkL)

    #r is the position of the "top" of all bulk layers in the model, ranging from 1 to opt["fine_bulkL"]+opt["bulkL"]
    for i = 1:bulkL
        r[i] = rout - OL_delta - (i - 1) * Lsize
    end

    for i = 2:Ltot-4
        #surface area of bulk layers 1:n-1 [cm2] 
        A[i+2] = 4 * pi * r[i]^2
    end
    for i = 1:Ltot-5
        #volume of bulk layers 1:n-1 [cm3]
        V[i+4] = 4 / 3 * pi * (r[i]^3 - r[1+i]^3)
    end

    #volume of bulk layer n [cm3]
    V[Ltot] = 4 / 3 * pi * ((r[Ltot-4])^3)

    #surface area of bulk layer n [cm2]
    A[Ltot-1] = 4 * pi * (r[Ltot-4])^2

    #f.s. gas diffusion shell 
    V[1] = 4 / 3 * pi * ((2 * mfp + rout + O3_delta)^3 - (mfp + rout + O3_delta)^3)

    #n.s. gas diffusion shell
    A[1] = 4 * pi * (rout + mfp + O3_delta)^2
    V[2] = 4 / 3 * pi * ((mfp + rout + O3_delta)^3 - (rout + O3_delta)^3)

    #sorption layer
    A[2] = 4 * pi * (rout + O3_delta)^2
    V[3] = 4 / 3 * pi * ((rout + O3_delta)^3 - rout^3)

    #static surface layer
    A[3] = 4 * pi * rout^2
    V[4] = 4 / 3 * pi * ((rout)^3 - (rout - OL_delta)^3)

    #diameter of each layer
    D[1] = mfp                                #thickness gas layer 1 [cm]
    D[2] = mfp                                #thickness gas layer 1 [cm]
    D[3] = O3_delta                           #thickness sorption layer [cm]
    D[4] = OL_delta                           #thickness static surface layer [cm] 
    D[5:bulkL+4] = ones(1, bulkL) * Lsize     #thickness bulk layer [cm]

    ### TRANSPORT ###

    #initialize
    kg1 = zeros(2)
    kd1 = zeros(2)
    kbs1 = zeros(2)
    kbs2 = zeros(2)
    ksb2 = zeros(2)
    kb = zeros(2, 1)

    #set values
    kg1[1] = O3_Dg / min(D[1], D[2])
    kd1[1] = 1 / O3_Td
    kbs1[1] = 2 * O3_Db / (D[3] + 2 * D[4] + D[5]) * O3_kbsxfac
    kbs2[2] = 2 * OL_Db / (D[4] + D[5]) * OL_kssbyfac
    ksb2[2] = 2 * OL_Db / OL_Ksurface_statbulk / (D[4] + D[5]) * OL_kssbyfac
    kb[1, 1] = 2 * O3_Db / (D[5] + D[6])
    kb[2, 1] = 2 * OL_Db / (D[5] + D[6])

    ### SOLVER OPTIONS ###

    #Output time step
    tspan::Tuple{Float64,Float64} = (10.0^(opt_start), 10.0^(opt_stop))

    #assign ODEparams
    ode_params = ODEParams(
        O3_a0,
        O3_Hcp,
        OL_Db,
        O3_Td,
        O3_sigma,
        O3_KHcc,
        O3_MTV,
        A,
        V,
        kg1,
        kd1,
        kbs1,
        kbs2,
        ksb2,
        kb,
        kbr,
        c1,
        c2
    )

    ### SOLVER CALL ###
    sol = solve(ODEProblem(OLEIC_F1, y0, tspan, ode_params), KenCarp3(), abstol=abstolvar, reltol=reltolvar)

    ### POST PROCESSING ###

    #get Z matrix
    Z = reduce(hcat, sol.u)
    Z = Z'

    #OL over time (normalized)
    Y::Vector{Float64} = A[4] * @views Z[:, 4]
    for j = 1:bulkL
        Y += V[4+j] * @views Z[:, 4+2*j]
    end
    Y ./= Y[1]
    
    return sol.t, Y

end