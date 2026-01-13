import torch
import torch.nn as nn

# --- 1. ResNet Block with BatchNorm (Improved Stability) ---
class ResNetBlock(nn.Module):
    def __init__(self, hidden_dim):
        super(ResNetBlock, self).__init__()
        self.fc1 = nn.Linear(hidden_dim, hidden_dim)
        self.bn1 = nn.BatchNorm1d(hidden_dim) # Key addition for deep networks
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_dim, hidden_dim)
        self.bn2 = nn.BatchNorm1d(hidden_dim)

    def forward(self, x):
        residual = x
        out = self.fc1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.fc2(out)
        out = self.bn2(out)
        out += residual  # The "Skip Connection"
        return self.relu(out)

# --- 2. Main Hybrid Model ---
class PINN_LSTM_ResNet(nn.Module):
    def __init__(self, input_dim=22, hidden_dim=128, num_res_blocks=2, output_dim=23):
        """
        Args:
            input_dim: 22 (Encoded Ephemeris features)
            output_dim: 23 (22 Params + 1 Latent 'E')
        """
        super(PINN_LSTM_ResNet, self).__init__()
        
        # LSTM Encoder: Captures Time-Series Trends
        self.lstm = nn.LSTM(input_dim, hidden_dim, batch_first=True)
        
        # ResNet Backend: Refines the features
        self.res_blocks = nn.ModuleList([
            ResNetBlock(hidden_dim) for _ in range(num_res_blocks)
        ])
        
        # --- specialized Heads ---
        # Clock Dynamics are simple (bias, drift), Orbit is complex.
        # We split the output to let heads specialize.
        
        # Head 1: Clock Parameters (Indices 0, 1, 2 in decoded 19-dim / 0, 1, 2 in encoded 22-dim)
        # Output: 3 values [ClockBias, Drift, DriftRate]
        self.clk_head = nn.Linear(hidden_dim, 3) 
        
        # Head 2: Orbital Parameters (The rest)
        # Output: 19 values (Original 19 params - 3 clock + 1 latent E + extra sin/cos dims)
        # We need to sum up to 'output_dim'.
        # Total Output = 23.
        # Clock = 3.
        # Orbit = 20.
        self.orbit_head = nn.Linear(hidden_dim, 20)

    def forward(self, x):
        # x shape: [batch, seq_len, input_dim]
        
        # 1. LSTM Pass
        lstm_out, _ = self.lstm(x) 
        # Take the last time step's feature
        feature = lstm_out[:, -1, :] 
        
        # 2. ResNet Pass
        for block in self.res_blocks:
            feature = block(feature)
            
        # 3. Prediction Heads
        clk_params = self.clk_head(feature)
        orbit_params = self.orbit_head(feature)
        
        # Recombine for loss calculation
        # [Clock(3) | Orbit(20)] = 23 output values
        return torch.cat([clk_params, orbit_params], dim=1)