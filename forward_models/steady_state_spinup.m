function Css = steady_state_spinup(N0, eps_erosion, kappa, lambda, z, targetC, Pnow, tol, maxit)
    % Choose a pseudo-timestep that respects near-surface CFL
    dz = gradient(z); dz(1)=z(2)-z(1); dz(end)=z(end)-z(end-1);
    if abs(eps_erosion) > 0
        dt_ps = 0.5 * min(dz) / max(abs(eps_erosion), eps);  % keeps local C<=0.5 per substep
    else
        dt_ps = 1.0;  % any positive value; diffusion+decay is implicit anyway
    end

    Cold = N0;
    for it = 1:maxit
        Cnew = step_forward_with_prod(Cold, eps_erosion, kappa, lambda, z, dt_ps, targetC, Pnow);
        % relative change
        denom = max(1e-16, norm(Cnew,2));
        rel = norm(Cnew-Cold,2) / denom;
        Cold = Cnew;
        if rel < tol
            break
        end
    end
    Css = Cold;
end
