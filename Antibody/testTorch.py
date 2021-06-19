import torch

if __name__ == "__main__":
    print("CUDA available: " + torch.cuda.is_available().__str__())
    print("GPU count: " + torch.cuda.device_count().__str__())