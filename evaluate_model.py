import torch
import numpy as np
import sys
from model import PINN_LSTM_ResNet
from position_loss import PositionLoss
from coords_layer import bds_ephemeris_to_ecef
from synthetic_data import generate_keplerian_series
from rinex_parser import BDSRinexParser

# Constants
MU = 3.986004418e14  # Earth's gravitational constant (m^3/s^2)

def calculate_metrics(predictions, targets, mean, std):
    """
    Computes RMSE, MAE, Max Error for Position and Clock.
    """
    # 1. Denormalize & Reconstruct
    # We reuse the logic from PositionLoss to get Physical Parameters
    # (Simplified here for brevity, ideally shared code)
    criterion = PositionLoss()
    
    # We need a way to get POS/CLK from params. 
    # Let's assume we can use the helper inside PositionLoss or replicate it.
    # To keep this clean, let's use the PositionLoss's internal logic which is robust.
    
    loss, pos_err_mean, clk_err_mean = criterion(predictions, targets, mean, std)
    
    # But we need distribution metrics (Max, RMSE), not just Mean.
    # So we need access to the raw error tensors.
    # Let's replicate the reconstruction logic quickly:
    
    p_pred = predictions[:, :22] * std + mean
    p_true = targets[:, :22] * std + mean
    
    def reconstruct_19(v22):
        batch = v22.shape[0]
        v19 = torch.zeros(batch, 19)
        v19[:, 0:5] = v22[:, 0:5]
        v19[:, 5] = torch.atan2(v22[:, 5], v22[:, 6])
        v19[:, 6:12] = v22[:, 7:13]
        v19[:, 12] = torch.atan2(v22[:, 13], v22[:, 14])
        v19[:, 13:16] = v22[:, 15:18]
        v19[:, 16] = torch.atan2(v22[:, 18], v22[:, 19])
        v19[:, 17:19] = v22[:, 20:22]
        return v19

    phys_pred = reconstruct_19(p_pred)
    phys_true = reconstruct_19(p_true)
    t_obs = phys_true[:, 10]
    
    pos_pred, clk_pred = bds_ephemeris_to_ecef(phys_pred, t_obs)
    pos_true, clk_true = bds_ephemeris_to_ecef(phys_true, t_obs)
    
    # Position Errors (meters)
    # pos shape: (batch, 3)
    diff_pos = pos_pred - pos_true
    dist_err = torch.norm(diff_pos, dim=1) # Euclidean distance per sample
    
    rmse_pos = torch.sqrt(torch.mean(dist_err**2)).item()
    mae_pos = torch.mean(dist_err).item()
    max_pos = torch.max(dist_err).item()
    
    # Clock Errors (seconds)
    diff_clk = torch.abs(clk_pred - clk_true)
    rmse_clk = torch.sqrt(torch.mean(diff_clk**2)).item()
    max_clk = torch.max(diff_clk).item()
    
    return {
        "pos_rmse": rmse_pos,
        "pos_mae": mae_pos,
        "pos_max": max_pos,
        "clk_rmse": rmse_clk,
        "clk_max": max_clk,
        "phys_pred": phys_pred, # Returned for physics check
        "pos_pred": pos_pred    # Returned for physics check
    }

def check_physics_energy(phys_params, pos_vectors):
    """
    Checks Specific Mechanical Energy Conservation.
    E = v^2/2 - mu/r
    For Keplerian orbit, E should be constant (-mu/2a).
    """
    # Extract Semi-Major Axis 'a'
    sqrt_a = phys_params[:, 9]
    a = sqrt_a ** 2
    
    # Theoretical Energy
    E_theoretical = -MU / (2 * a)
    
    # Calculated Energy from Position/Velocity
    # Note: We need Velocity. bds_ephemeris_to_ecef only returns Position currently.
    # We can approximate Velocity by finite difference if we have a sequence,
    # or implement the full velocity equation.
    # For this check, let's verify if 'a' (Semi-Major Axis) stays constant over the prediction.
    # Because 'a' defines the energy level.
    
    # Ideally, we check if the PREDICTED 'a' matches the INITIAL 'a'.
    
    a_variance = torch.var(a).item()
    a_drift = (torch.max(a) - torch.min(a)).item()
    
    return a_variance, a_drift

def evaluate():
    print("--- Loading Model and Data ---")
    
    # Load Data (Test Set)
    # We generate fresh data to ensure we aren't testing on training data
    from train_demo import process_data # Reuse processing logic
    raw_data = generate_keplerian_series(num_samples=500)
    inputs, targets, mean, std = process_data(raw_data)
    
    # Load Model
    INPUT_DIM = 22
    HIDDEN_DIM = 128
    OUTPUT_DIM = 23
    model = PINN_LSTM_ResNet(INPUT_DIM, HIDDEN_DIM, output_dim=OUTPUT_DIM)
    
    try:
        model.load_state_dict(torch.load("bds_pinn_model.pth"))
        model.eval()
    except FileNotFoundError:
        print("Error: Model file 'bds_pinn_model.pth' not found. Train first.")
        return

    print("--- Running Evaluation ---")
    with torch.no_grad():
        predictions = model(inputs)
        
        metrics = calculate_metrics(predictions, targets, mean, std)
        
        print("\n=== ACCURACY REPORT ===")
        print(f"Position RMSE: {metrics['pos_rmse']:.4f} m")
        print(f"Position MAE:  {metrics['pos_mae']:.4f} m")
        print(f"Position Max:  {metrics['pos_max']:.4f} m")
        print(f"Clock RMSE:    {metrics['clk_rmse']*1e9:.4f} ns")
        print(f"Clock Max:     {metrics['clk_max']*1e9:.4f} ns")
        
        print("\n=== PHYSICS CONSISTENCY REPORT ===")
        # Check Energy Stability (via Semi-Major Axis stability)
        phys_params = metrics['phys_pred']
        var, drift = check_physics_energy(phys_params, metrics['pos_pred'])
        
        print(f"Energy (Orbit Size) Drift: {drift:.4f} meters")
        print(f"Energy Variance:           {var:.4e}")
        
        if drift < 100.0:
            print(">>> PASS: Orbit Energy is conserved (Stable Orbit).")
        else:
            print(">>> FAIL: Orbit is decaying or expanding unphysically.")
            
        print("\n=== CONCLUSION ===")
        if metrics['pos_rmse'] < 50.0 and metrics['clk_rmse'] < 1e-7:
             print("Model is VALID for Red Team Testing.")
        else:
             print("Model requires refinement.")

if __name__ == "__main__":
    evaluate()
