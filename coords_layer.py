import torch
import numpy as np

# Constants for BDS (CGCS2000)
MU = 3.986004418e14  # Earth's gravitational constant
OMEGA_E = 7.2921150e-5  # Earth rotation rate

def bds_ephemeris_to_ecef(params, t_obs):
    """
    Differentiable conversion from BDS Ephemeris parameters to ECEF Coordinates (X, Y, Z).
    
    Args:
        params: Tensor (Batch, 22/23). Our model output.
                Expects specific order of encoded features.
                We need to decode Sin/Cos back to angles first if encoded.
        t_obs:  Time of Observation (seconds within week).
        
    Returns:
        pos: Tensor (Batch, 3) -> [X, Y, Z] in meters
        clk: Tensor (Batch, 1) -> Clock Bias in seconds
    """
    
    # 1. Extract Parameters (Assuming un-normalized, decoded format for calculation)
    # We assume the input 'params' here are already De-Normalized and Angles Decoded
    # Order based on our parser:
    # 0: ClockBias, 1: Drift, 2: DriftRate
    # 3: Crs, 4: DeltaN, 5: M0
    # 6: Cuc, 7: e, 8: Cus, 9: SqrtA
    # 10: Toe, 11: Cic, 12: Omega0, 13: Cis, 14: i0, 15: Crc, 16: omega, 17: OmegaDot, 18: IDOT
    
    # Clock Correction
    a0 = params[:, 0]
    a1 = params[:, 1]
    a2 = params[:, 2]
    toc = params[:, 10] # Using Toe as Toc for simplicity in broadcast
    
    dt = t_obs - toc
    # Check for week crossover (simplification)
    dt = dt - torch.round(dt / 604800.0) * 604800.0
    
    sat_clk_err = a0 + a1 * dt + a2 * (dt ** 2)
    
    # Orbital Calculation
    sqrt_A = params[:, 9]
    A = sqrt_A ** 2
    n0 = torch.sqrt(MU / (A ** 3))
    delta_n = params[:, 4]
    n = n0 + delta_n
    
    M0 = params[:, 5]
    tk = dt # Time since ephemeris reference epoch
    M_k = M0 + n * tk
    
    # Solve Kepler's Eq for E (Eccentric Anomaly)
    # Iterative solution is not differentiable easily in loop, 
    # but for small e, M approx E. Or we use the Model's Predicted E if available!
    # Let's use the Model's latent E output if valid, otherwise approx.
    # For robust loss, we usually perform 2-3 fixed iterations (unrolled loop).
    
    e = params[:, 7]
    E_k = M_k
    # 3 Iterations of Newton-Raphson (sufficient for GPS/BDS)
    for _ in range(3):
        E_k = M_k + e * torch.sin(E_k)
        
    # True Anomaly v_k
    sin_v = (torch.sqrt(1 - e**2) * torch.sin(E_k)) / (1 - e * torch.cos(E_k))
    cos_v = (torch.cos(E_k) - e) / (1 - e * torch.cos(E_k))
    v_k = torch.atan2(sin_v, cos_v)
    
    # Argument of Latitude u_k
    omega = params[:, 16]
    phi_k = v_k + omega
    
    # Second Harmonic Perturbations
    cus = params[:, 8]; cuc = params[:, 6]
    crs = params[:, 3]; crc = params[:, 15]
    cis = params[:, 13]; cic = params[:, 11]
    
    sin_2phi = torch.sin(2 * phi_k)
    cos_2phi = torch.cos(2 * phi_k)
    
    du_k = cus * sin_2phi + cuc * cos_2phi
    dr_k = crs * sin_2phi + crc * cos_2phi
    di_k = cis * sin_2phi + cic * cos_2phi
    
    u_k = phi_k + du_k
    r_k = A * (1 - e * torch.cos(E_k)) + dr_k
    i_k = params[:, 14] + params[:, 18] * tk + di_k # i0 + IDOT*tk
    
    # Position in Orbital Plane
    x_k_prime = r_k * torch.cos(u_k)
    y_k_prime = r_k * torch.sin(u_k)
    
    # Corrected Longitude of Ascending Node
    Omega0 = params[:, 12]
    OmegaDot = params[:, 17]
    Omega_k = Omega0 + (OmegaDot - OMEGA_E) * tk - OMEGA_E * toc
    
    # Earth Fixed Coordinates (ECEF)
    x = x_k_prime * torch.cos(Omega_k) - y_k_prime * torch.cos(i_k) * torch.sin(Omega_k)
    y = x_k_prime * torch.sin(Omega_k) + y_k_prime * torch.cos(i_k) * torch.cos(Omega_k)
    z = y_k_prime * torch.sin(i_k)
    
    # Stack position
    pos = torch.stack([x, y, z], dim=1)
    
    return pos, sat_clk_err
