
using QuadGK
using Polynomials
using Plots
using Roots

# funkcja obliczająca odcięte punktów Gaussa oraz ich współczynniki 
(xp,a)=gauss(Float64,7)

f(x)=x^2
sum(a .* f.(xp)) 

function legendre(n)
    p0 = Poly([1.])
    p1 = Poly([0., 1.])  
    if n == 0
        return p0
    end
    for k = 1:(n - 1)
        p = Poly([0., (2 * k + 1.)]) * p1
        p2 = (p - (k * p0)) / (k + 1.)
        p0 = p1
        p1 = p2
    end
    return p1
end

x = -1.:0.05:1.
for i = 2:5
    y = legendre(i)(x)
    p = plot!(x, y, label = "P_$i")
end
plot!(size = (950, 500), legend = :bottomright)

for i = 2:4
    println("degree: ", i)
    r = find_zeros(legendre(i), -1., 1.)
    println("\troots: ", r)
    (xp, a) = gauss(Float64, i)
    println("\tgauss: ", xp)
end

function integrate(f, k)
    (xp, w) = gauss(Float64, k)
    sum(w .* f.(xp))  
end

eps = 1.0e-9

for qdeg = 1:5     
    exact = true
    pdeg = 0
    p = Poly([rand()])
    
    while exact   
        p *= Poly([rand(), rand()])
        pdeg += 1
        
        res = polyint(p, -1, 1)
        qres = integrate(p, qdeg)
        diff = abs(qres - res)
        
        print("qdeg = ", qdeg, ", pdeg = ", pdeg)
        println(" | res = ", res, ", qres = ", qres, ", diff = ", diff)
        
        if diff > eps
            exact = false
        end
    end
    
    println("qdeg = ", qdeg, " ==> exact up to ", pdeg - 1, "-degree poly\n")
end

function integrate(f, k, a, b)
    fn = x -> f((b + a)/2 + (b - a)/2 * x)
    return (b - a)/2 * integrate(fn, k)
end

p1 = Poly([1., 2., 3.]) # 3x^2 + 2x + 1 => x^3 + x^2 + x
println("polyint:   ", polyint(p1, -0.4, 3))
println("integrate: ", integrate(p1, 4, -0.4, 3))

println("polyint:   ", exp(3.3) - exp(-1.1))
println("integrate: ", integrate(exp, 5, -1.1, 3.3))

p1 = Poly([1., 2., 3.]) # 3x^2 + 2x + 1 => x^3 + x^2 + x
println("polyint: ", polyint(p1, -5, 7))
println("quadgk: ", quadgk(p1, -5, 7))

gaussian = x -> (1 / (sqrt(2*pi))) * exp(-(x*x)/2)
println("quadqk for gaussian: ", quadgk(gaussian, -Inf, Inf))
plot(gaussian, -4, 4, label = "gaussian distribution", size = (900, 350))

function trapezoid(f, a, b, n)
    sum = 0
    dx = (b - a) / n
    xs = range(a, stop = b, length = n + 1)
    for i = 1:n
        sum = sum + (f(xs[i]) + f(xs[i+1])) * dx / 2
    end
    sum
end

p1 = (1/20) * Poly([64., 10., -27., -11., 3., 1.])
println(p1)
a = -4
b = 2
plot(x -> p1(x), a-0.1, b+0.1, label = "p1 - sample poly")

for n = 1:2:7
    xs = range(a, stop = b, length = n + 1)
    p = plot!(xs, p1(xs), label = "trapezoid n = $n")
end
plot!(size = (900, 350), legend = :bottomleft)

exact = polyint(p1, a, b)
ns = 1:2:200
err = [abs(exact - trapezoid(p1, a, b, n)) for n in ns]
scatter(ns, err, xlabel = "subintervals", ylabel = "absolute error",
    label = "abs_error(subintervals)", yscale = :log10, size = (900, 400))
