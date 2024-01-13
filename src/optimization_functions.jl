# KL-divergence function
function H_Function(p, q) 
    return sum(p .* NaNMath.log.(p ./ clamp.(q, 1e-99, Inf)))
end

# Constraint function
function H_cons(p, psi, K, pi_tea, lambda1, lambda2)
    ϵ = 0
    
    pi_tea = sum(pi_tea)
    lambda1 = sum(lambda1)
    lambda2 = sum(lambda2)
        
    
    tau = p[1:K[1]]; tau_star = clamp.(cumsum(tau), 1e-99, 1 - 1e-99)
    u1 = p[K[1]+1:K[1]+K[2]]; u1_star = clamp.(cumsum(u1), 1e-99, 1 - 1e-99)
    u2 = p[K[1]+K[2]+1:end]; u2_star = clamp.(cumsum(u2), 1e-99, 1 - 1e-99)
    
    

    # Check if any values are outside the [0, 1) range
    if any(t -> t < 0 || t > 1, tau_star)
        println("tau out of bounds")
        # Handle the situation, maybe by setting a flag or adjusting values
    end
    
    if any(t -> t < 0 || t > 1, u1_star)
        println("u1 out of bounds")
        # Handle the situation, maybe by setting a flag or adjusting values
    end
    
    if any(t -> t < 0 || t > 1, u2_star)
        println("u2 out of bounds")
        # Handle the situation, maybe by setting a flag or adjusting values
    end
    
    inv_tau_star = clamp.(1 .- tau_star, 1e-99, 1 - 1e-99)
    inv_u1_star = clamp.(1 .- u1_star, 1e-99, 1 - 1e-99)
    inv_u2_star = clamp.(1 .- u2_star, 1e-99, 1 - 1e-99)

    csi = -sum(inv_tau_star .* log.(inv_tau_star))
    zeta1 = -sum(inv_u1_star .* log.(inv_u1_star))
    zeta2 = -sum(inv_u2_star .* log.(inv_u2_star))

    ceq = [real(pi_tea - sum(tau)), real(lambda1 - sum(u1)), real(lambda2 - sum(u2)), real(psi - csi - zeta1 - zeta2)]

    penalty = 0
    for val in ceq
        if abs(val) > 1e-6  # Adjust tolerance as needed
            penalty += 1e6 * abs(val)  # Large penalty for constraint violation
        end
    end

    return penalty, csi, zeta1, zeta2, tau, tau_star, u1, u1_star, u2, u2_star
end

# Objective function
function objective(p, q, psi, K, pi_tea, lambda1, lambda2)
    penalty, _, _, _, _, _, _, _, _, _,  = H_cons(p, psi, K, pi_tea, lambda1, lambda2)
    return H_Function(p, q) + penalty
end