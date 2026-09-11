using ParameterizedFunctions, BenchmarkTools

const SUITE = BenchmarkGroup()

# @ode_def builds an ODEFunction with analytic Jacobian/etc from DSL
lv = @ode_def LotkaVolterra begin
    dx = a * x - b * x * y
    dy = -c * y + d * x * y
end a b c d

f = lv  # ODEFunction-like callable
u0 = [1.0, 1.0]
p = [1.5, 1.0, 3.0, 1.0]
du = similar(u0)
J = zeros(2, 2)
t = 0.0

# =============================================================================
# Function evaluation
# =============================================================================

SUITE["call"] = BenchmarkGroup()

SUITE["call"]["f"] = @benchmarkable $f($du, $u0, $p, $t)
SUITE["call"]["f_oop"] = @benchmarkable $f($u0, $p, $t)
SUITE["call"]["jac"] = @benchmarkable $(f.jac)($J, $u0, $p, $t)

# =============================================================================
# Problem construction through the generated function
# =============================================================================

SUITE["problem"] = BenchmarkGroup()

SUITE["problem"]["odeproblem"] = @benchmarkable ODEProblem($f, $u0, (0.0, 10.0), $p)

# Larger generated system (codegen cost at macro expansion is load-time;
# benchmark steady-state evaluation of a 20-component system)
big = @ode_def BigSys begin
    dx1 = -x1 + a
    dx2 = x1 - x2
    dx3 = x2 - x3
    dx4 = x3 - x4
    dx5 = x4 - x5
end a
u0b = ones(5)
pb = [1.0]
SUITE["call"]["big_sys"] = @benchmarkable $big($(zeros(5)), $u0b, $pb, $t)
