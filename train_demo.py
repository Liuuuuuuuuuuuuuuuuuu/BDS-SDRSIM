import torch
import torch.optim as optim
import numpy as np
import sys
from torch.utils.data import TensorDataset, DataLoader
from model import PINN_LSTM_ResNet # Updated Class Name
from position_loss import PositionLoss
from rinex_parser import BDSRinexParser
from synthetic_data import generate_keplerian_series

# ... (rest of imports and helper functions encode_data/load_data/process_data/prepare_sequences remain same)

def train(data_path=None):
    INPUT_DIM = 22
    HIDDEN_DIM = 128
    LAYERS = 2
    OUTPUT_DIM = 23
    LR = 0.002
    EPOCHS = 20

    inputs, targets, mean, std = load_data(data_path)
    if inputs is None: 
        print("Error: No data to train on.")
        return
        
    split_idx = int(len(inputs) * 0.8)
    train_in, test_in = inputs[:split_idx], inputs[split_idx:]
    train_tgt, test_tgt = targets[:split_idx], targets[split_idx:]
    
    train_loader = DataLoader(TensorDataset(train_in, train_tgt), batch_size=32, shuffle=True)
    
    # Use the new Architecture
    model = PINN_LSTM_ResNet(input_dim=INPUT_DIM, hidden_dim=HIDDEN_DIM, output_dim=OUTPUT_DIM)
    
    criterion = PositionLoss(weight_pos=1.0, weight_time=1e9) 
    optimizer = optim.Adam(model.parameters(), lr=LR)


    print(f"--- Starting Mini-Batch Training with Coordinate Loss ---")
    
    for epoch in range(EPOCHS):
        model.train()
        total_pos = 0; total_clk = 0
        
        for batch_in, batch_tgt in train_loader:
            optimizer.zero_grad()
            predictions = model(batch_in)
            loss, pos_err, clk_err = criterion(predictions, batch_tgt, mean, std)
            loss.backward()
            optimizer.step()
            total_pos += pos_err.item(); total_clk += clk_err.item()
        
        if epoch % 5 == 0:
            print(f"Epoch {epoch} | Pos: {total_pos/len(train_loader):.1f}m | Clk: {(total_clk/len(train_loader))*1e9:.1f}ns")

    print("\n--- Final Validation ---")
    model.eval()
    with torch.no_grad():
        test_preds = model(test_in)
        _, final_pos, final_clk = criterion(test_preds, test_tgt, mean, std)
        print(f"Test Set Error -> Pos: {final_pos.item():.2f}m | Clk: {final_clk.item()*1e9:.2f}ns")

    torch.save(model.state_dict(), "bds_pinn_model.pth")
    print("SUCCESS: Model saved.")

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else None
    train(path)