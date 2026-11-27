"""
1. Inspect the UNIVERSAL continuum function:
     Psi_{D, tau}(W) = -log(W) + int_{-inf}^{inf} log(1 - W * exp(-2*pi*i*D*x - tau*x*(1+x))) dx
   at W = -Lambda (Lambda > 0) for D = 0, 1, 2, 3 and tau in [0.2, 15.0]!
   Compute:
     - Saddle radius Lambda_*(D, tau) solving W * Psi'(W) = 0 on W = -Lambda (or global max of Re Psi on |W|=Lambda)!
     - Curvature H_cont(D, tau) = -d^2/dtheta^2 Re Psi(-Lambda_* e^{i theta}) |_{theta=0}
     - Does H_cont(D, tau) stay STRICTLY POSITIVE for ALL tau > 0 and ALL D >= 1?!
     - Does Lambda_*(D, tau) stay strictly separated from exp(-tau/4) for ALL tau > 0 when D is odd?!
2. Test our exact 2-step Abel summation bound on |B_r - A_r| / r at lambda = 1!
"""
import numpy as np

def eval_psi_cont_derivatives(D, tau, lam, M_x=4000):
    # Integrate x in [-L, L] centered at x = -0.5: x = -0.5 + u
    L = np.sqrt(45.0 / max(tau, 1e-4)) + 1.0
    us = np.linspace(-L, L, M_x)
    du = us[1] - us[0]
    # W = -lam => -W * exp(-2pi i D x - tau x(1+x)) = +lam * (-1)^D * exp(-2pi i D u - tau*(u^2 - 0.25))
    re_z = np.log(lam) + tau * (0.25 - us * us)
    im_z = -2.0 * np.pi * D * us + (np.pi if D % 2 == 1 else 0.0)
    # Wait: in F(W) = W^{-c} prod (1 - a_m W), at W = -lam:
    # 1 - a_m W = 1 + lam * (-1)^D * exp(-2pi i D u - tau*(u^2-0.25))
    # Let z = -a_m W = +lam * (-1)^D * exp(...) so factor is (1 + z)!
    # Wait: check sign! Earlier a_m = exp(-2pi i D m/c - alpha m(c+m)).
    # At x = m/c = -0.5 + u:
    #   -2pi i D (-0.5 + u) = +i pi D - 2pi i D u => exp = (-1)^D exp(-2pi i D u).
    #   So -a_m W = +lam * (-1)^D * exp(-2pi i D u - tau(u^2 - 0.25)).
    # Factor is (1 + z)!
    # d/dtheta log(1 + z e^{i theta}) = i * z / (1 + z)
    # d^2/dtheta^2 log(1 + z e^{i theta}) = -z / (1 + z)^2
    # Wait: if |z| > 1: z/(1+z) = 1 - 1/(1+z) = 1 - (1/z)/(1 + 1/z)
    #       z/(1+z)^2 = (1/z)/(1 + 1/z)^2
    big = re_z > 0.0
    g1 = -1.0 + 0.0j
    H = 0.0j
    if np.any(big):
        inv_z = np.exp(-re_z[big] - 1j * im_z[big])
        g1 += du * np.sum(1.0 - inv_z / (1.0 + inv_z))
        H += du * np.sum(inv_z / ((1.0 + inv_z) ** 2))
    sm = ~big
    if np.any(sm):
        z = np.exp(re_z[sm] + 1j * im_z[sm])
        g1 += du * np.sum(z / (1.0 + z))
        H += du * np.sum(z / ((1.0 + z) ** 2))
    return g1.real, H.real

print("=== UNIVERSAL CONTINUUM SADDLE CURVE Lambda_*(D, tau) & CURVATURE H_cont(D, tau) ===")
for D in [0, 1, 2, 3]:
    print(f"\n--- D = {D} ---")
    print(f"{'tau':>6} {'e^{-tau/4}':>10} {'g1@1':>9} | {'lam_*':>9} {'lam*/e^-t/4':>11} {'H_cont':>9}")
    for tau in [0.4, 0.8, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0]:
        em = np.exp(-tau / 4.0)
        g1_at_1, _ = eval_psi_cont_derivatives(D, tau, 1.0)
        # Solve g1(lam) = 0 on (em*1.005 if D%2==1 else 0.01, 5.0)
        lo = em * 1.002 if (D % 2 == 1) else 0.01
        hi = 5.0
        g_lo, _ = eval_psi_cont_derivatives(D, tau, lo)
        g_hi, _ = eval_psi_cont_derivatives(D, tau, hi)
        if g_lo * g_hi <= 0:
            for _ in range(35):
                mid = 0.5 * (lo + hi)
                gm, _ = eval_psi_cont_derivatives(D, tau, mid)
                if gm * g_lo <= 0:
                    hi = mid
                else:
                    lo = mid; g_lo = gm
            lam_s = 0.5 * (lo + hi)
            _, Hs = eval_psi_cont_derivatives(D, tau, lam_s)
            print(f"{tau:6.2f} {em:10.5f} {g1_at_1:9.4f} | {lam_s:9.5f} {lam_s/em:11.5f} {Hs:9.4f}")
        else:
            # Scan minimum of |g1| on [lo, hi]
            ls = np.linspace(lo, 2.0, 100)
            gs = [eval_psi_cont_derivatives(D, tau, x)[0] for x in ls]
            idx = int(np.argmin(np.abs(gs)))
            _, Hs = eval_psi_cont_derivatives(D, tau, ls[idx])
            print(f"{tau:6.2f} {em:10.5f} {g1_at_1:9.4f} | NO_ROOT (min g1={gs[idx]:.4f} at lam={ls[idx]:.4f}, H={Hs:.4f})")
