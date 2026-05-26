## setup
using DifferentialEquations
#using GLMakie
using GLMakie



## define model
#= lotka-volterra 2 species competition =#
function LV!(dn, n, p, t)
    #= dn[1] = pop1*rate1(carry_cap1 - comp_coeff_1invades2*pop1 - comp_coeff_2inv1*pop2)
    dn[2] = pop2*rate2(carry_cap2 - comp_coeff_2invades1*pop2 - comp_coeff_1inv2*pop1) =#
    (; r1, r2, α11, α12, α21, α22) = p # the ; makes it a destructured tuple. for cleaner ordering
    n1, n2 = n
    # dn[1] = n1 * r1 * (K1 - α11*n1 - α12*n2) / K1
    # dn[2] = n2 * r2 * (K2 - α22*n2 - α21*n1) / K2
    dn[1] = n1 * (r1 - α11 * n1 - α12 * n2)
    dn[2] = n2 * (r2 - α21 * n1 - α22 * n2)

    nothing
end
## run
n0 = [0.1, 0.12] # initial conditions
tspan = (0.0, 40.0)
# p  = (r1=1.2, r2=1.0, K1=1.0, K2=1.0, α11=1.0, α12=1.0, α21=0.5, α22=1.0)
p = (
    r1=1.2, 
    r2=1.0, 
    α11=1.0, 
    α12=1.0, 
    α21=0.5, 
    α22=1.0)

prob = ODEProblem{true}(LV!, n0, tspan, p)
sol = solve(prob, Tsit5(); reltol = 1e-6, abstol = 1e-6)
##plot
fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Time",
    ylabel = "Population",
    xautolimitmargin = (0.0, 0.0)  # remove padding
)

lines!(ax, sol.t, sol[1, :], label="Species 1", color=:blue)
lines!(ax, sol.t, sol[2, :], label="Species 2", color=:red)
axislegend(ax)

param_str = join(["$(k) = $(v)" for (k, v) in pairs(p)], "    ")
Label(fig[2,1], param_str, tellwidth = false)
# this figure shows the dynamics of a two species lotka-volterra competition model from t = [ 0.0, 40.0 ]
fig

# save("lotka_volterra.png", fig) 
