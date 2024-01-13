function entropy_decomposition(H_theta, H_theta0, xHist, maxInner, maxOuter, maxFunEvals; decompose=true)

    # Assuming H_theta, H_theta0, and xH are defined and have the same shape
    phi_omega = clamp.(cumsum(H_theta / sum(H_theta)), 1e-99, 1-1e-99)
    inv_phi_omega = clamp.(1 .- phi_omega, 1e-99, 1-1e-99)
    psi = -sum(inv_phi_omega .* log.(inv_phi_omega))

    if !decompose
        return Dict("psi" => psi, "csi" => 0, "zeta" => 0, "converged" => 0, "penalty" => 0)
    else
        # Computing proxy for CSI, ZETA1, ZETA2
        pi_tea = clamp.(H_theta0 / sum(H_theta0),1e-99, 1-1e-99)
        U = H_theta .- H_theta0
        med = findfirst(xHist .== median(xHist))
        
        U1 = U[1:(med-1)]
        U2 = U[(med+1):end]
        lambda1 = clamp.(U1 / (sum(U1) + 1e-99), 1e-99, 1-1e-99)
        lambda2 = clamp.(U2 / (sum(U1) + 1e-99), 1e-99, 1-1e-99)

        # Set parameters for the minimum KL-procedure
        K = [length(pi_tea), length(lambda1), length(lambda2)]
        q = [pi_tea; lambda1; lambda2]
        p0 = [ones(K[1]) .* 1 / K[1]; ones(K[2]) .* 1 / K[2]; ones(K[3]) .* 1 / K[2]]
        lb = ones(sum(K)) .* 1e-99
        ub = ones(sum(K))
        
        optimizer_opts = Optim.Options(f_calls_limit = maxFunEvals, iterations = maxInner, outer_iterations = maxOuter, show_trace = false)  # Adjust as needed

        results = optimize(p -> objective(p, q, psi, K, pi, lambda1, lambda2), lb, ub, p0, Fminbox(LBFGS()), optimizer_opts)
        converged = Optim.converged(results)
        
        optimized_p = Optim.minimizer(results)
        
        penalty, csi, zeta1, zeta2, tau, tau_star, u1, u1_star, u2, u2_star = H_cons(optimized_p, psi, K, pi_tea, lambda1, lambda2)

        return Dict("psi" => psi, "tau" => tau, "tau_star" => tau_star, 
                    "u1" => u1, "u1_star" => u1_star, "u2" => u2, "u2_star" => u2_star, 
                    "csi" => csi, "zeta1" => zeta1, "zeta2" => zeta2, "zeta" => zeta1 + zeta2,
                    "H_theta" => H_theta, "H_theta0" => H_theta0, "xHist" => xHist, 
                    "U" => U, "U1" => U1, "U2" => U2, "pi_tea" => pi_tea, "lambda1" => lambda1, 
                    "lambda2" => lambda2, "converged" => converged, "penalty" => penalty)
    end
end