import torch
import torch.nn as nn
from coords_layer import bds_ephemeris_to_ecef

class PositionLoss(nn.Module):
    def __init__(self, weight_pos=1.0, weight_time=1e9):
        super(PositionLoss, self).__init__()
        self.mse = nn.MSELoss()
        self.w_pos = weight_pos
        self.w_time = weight_time

    def forward(self, pred_params, true_params, mean, std):
        # 1. Denormalize
        p_pred = pred_params[:, :22] * std + mean
        p_true = true_params[:, :22] * std + mean
        
        # 2. Reconstruct 19-dim physical params from 22-dim encoded params
        def reconstruct_19(v22):
            batch = v22.shape[0]
            v19 = torch.zeros(batch, 19, device=v22.device)
            
            # Linear 0-4 (ClockBias, Drift, DriftRate, Crs, DeltaN)
            v19[:, 0:5] = v22[:, 0:5]
            # M0 (5,6)
            v19[:, 5] = torch.atan2(v22[:, 5], v22[:, 6])
            # Linear 6-11 in 19-dim (Cuc, e, Cus, SqrtA, Toe, Cic) 
            # These are indices 7,8,9,10,11,12 in 22-dim
            v19[:, 6:12] = v22[:, 7:13]
            # Omega0 (13,14)
            v19[:, 12] = torch.atan2(v22[:, 13], v22[:, 14])
            # Linear 13-15 (Cis, i0, Crc) -> indices 15,16,17 in 22-dim
            v19[:, 13:16] = v22[:, 15:18]
            # omega (18,19)
            v19[:, 16] = torch.atan2(v22[:, 18], v22[:, 19])
            # Linear 17-18 (OmegaDot, IDOT) -> indices 20,21 in 22-dim
            v19[:, 17:19] = v22[:, 20:22]
            
            return v19

        phys_pred = reconstruct_19(p_pred)
        phys_true = reconstruct_19(p_true)
        
        # 3. Position and Time calculation
        # Observation time t_obs = Toe (index 10)
        t_obs = phys_true[:, 10]
        
        pos_pred, clk_pred = bds_ephemeris_to_ecef(phys_pred, t_obs)
        pos_true, clk_true = bds_ephemeris_to_ecef(phys_true, t_obs)
        
        # 4. Losses
        pos_err = torch.norm(pos_pred - pos_true, dim=1)
        loss_pos = torch.mean(pos_err ** 2)
        
        # Clock error is in seconds
        clk_err = torch.abs(clk_pred - clk_true)
        loss_clk = torch.mean(clk_err ** 2)
        
        total_loss = (self.w_pos * loss_pos) + (self.w_time * loss_clk)
        
        return total_loss, torch.mean(pos_err), torch.mean(clk_err)