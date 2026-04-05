import numpy as np
from typing import Optional

def random_skew_symmetric_F2(m):
    upper = np.triu(np.random.randint(0, 2, size=(m, m), dtype=np.uint8), k=1)
    
    A = upper + upper.T
    
    return A % 2  

def f2_rank(arr: np.ndarray) -> int:
    """Compute rank of a matrix over F2 using Gaussian elimination."""
    M = arr.copy().astype(int) % 2
    rows, cols = M.shape
    rank = 0
    
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if M[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        M[[rank, pivot]] = M[[pivot, rank]]
        for row in range(rows):
            if row != rank and M[row, col] == 1:
                M[row] = (M[row] + M[rank]) % 2
        rank += 1
    return rank


class LinearSpaceF2:

    def __init__(self, generators: list[np.ndarray]):
        self.generators = [g for g in generators if np.any(g != 0)]
        self.m = generators[0].size if generators else 0
    
    def is_element(self, v: np.ndarray):
        if v.size != self.m:
            return False
        if not self.generators:
            return np.all(v == 0)
        arr_with = np.array(self.generators + [v], dtype=np.uint8)
        rank_with_v = f2_rank(arr_with)
    
        arr_gen = np.array(self.generators, dtype=np.uint8)
        rank_without_v = f2_rank(arr_gen)
        
        return rank_with_v == rank_without_v        
        
    def add_generator(self, v: np.ndarray):
        if np.all(v == 0):
            return
        if not self.is_element(v):
            self.generators.append(v)

    def random_element(self) -> np.ndarray:
        if not self.generators:
            return np.zeros(self.m, dtype=np.uint8)
        coeffs = np.random.randint(0, 2, len(self.generators))
        return sum(c * g for c, g in zip(coeffs, self.generators)) % 2

    def dim(self) -> int:
        return len(self.generators)
        
    def show_generators(self, labels: bool = True):
        for i, v in enumerate(self.generators):
            if labels:
                print(f"v{i+1}: {v}")
            else:
                print(v)


def f2_nullspace(arr: np.ndarray) -> list[np.ndarray]:
    M = arr.copy().astype(int) % 2
    rows, cols = M.shape
    
    # Row reduce to RREF
    pivot_cols = []
    row = 0
    
    for col in range(cols):
        if row >= rows:
            break
        # Find pivot
        pivot_row = None
        for r in range(row, rows):
            if M[r, col] == 1:
                pivot_row = r
                break
        if pivot_row is None:
            continue
        # Swap
        M[[row, pivot_row]] = M[[pivot_row, row]]
        # Eliminate
        for r in range(rows):
            if r != row and M[r, col] == 1:
                M[r] = (M[r] + M[row]) % 2
        pivot_cols.append((row, col))
        row += 1
    
    pivot_set = set(c for _, c in pivot_cols)
    free_vars = [c for c in range(cols) if c not in pivot_set]

    basis = []
    for free in free_vars:
        v = np.zeros(cols, dtype=np.uint8)
        v[free] = 1
        for r, pc in pivot_cols:
            v[pc] = M[r, free]
        basis.append(v)
    
    return basis


class LinearMapF2:
    def __init__(self, matrix: np.ndarray):
        self.matrix = matrix.astype(np.uint8) % 2
        self.m = matrix.shape[0]
        self.n = matrix.shape[1]
    
    def apply(self, v: np.ndarray) -> np.ndarray:
        return (self.matrix @ v.astype(np.uint8)) % 2
    
    def kernel(self) -> "LinearSpaceF2":
        basis = f2_nullspace(self.matrix)
        return LinearSpaceF2(basis)


class skewSymetricFormF2:

    def __init__(self, B: list[np.ndarray]):
        self.B = B
        self.m = B[0].shape[0]
        self.n = len(B)

    def show(self):
        for i in range(self.n):
            print("B" + str(i))
            print(self.B[i])
        
    def prod(self,x: np.ndarray, y: np.ndarray):
        result = np.zeros(self.n, dtype = np.uint8)
        for i, B in enumerate(self.B):
            result[i] = (x @ B @ y.T).item() % 2
        return result

    def myrank(self):
        num_pairs = self.m * (self.m - 1) // 2
        beta_tilde = np.zeros((num_pairs, self.n), dtype=np.int8)
        
        row = 0
        for i in range(self.m):
            for j in range(i + 1, self.m):
                for α, B in enumerate(self.B):
                    beta_tilde[row, α] = B[i, j]
                row += 1
        
        rank = f2_rank(beta_tilde)
        return rank

    def findK_RSBSM(self, depth: int = 10):
        """ 
        RSBSM - Random Step By Step Method:
            Take random elemnt from W and construct step by step bigger K.
            depth: controls how many random vectors to try before giving up
        """
        m = self.m
        n = self.n
        W_gen = []
        for i in range(m):
            e = np.array([0 for i in range(m)], dtype = np.uint8)
            e[i] = 1
            W_gen.append(e)

        W = LinearSpaceF2(W_gen)
        x0 = W.random_element()
        while np.all(x0 == 0):
            x0 = W.random_element()
        K = LinearSpaceF2([x0])

        iterator = 0

        while True:
            
            if iterator > 20:
                break

            r = K.dim()
            Psi_matrix = np.zeros((r*n, m), dtype = np.uint8)
            for i, xi in enumerate(K.generators):
                for α, B in enumerate(self.B):
                    row_idx = i * self.n + α
                    Psi_matrix[row_idx] = (B @ xi) % 2
            
            Psi = LinearMapF2(Psi_matrix)
            Psi_kernel = Psi.kernel()
            
            random_control = 0
            x = W.random_element()
            while not Psi_kernel.is_element(x) or np.all(x == 0):
                if random_control > depth:
                    break
                x = W.random_element()
                random_control += 1
            
            if random_control > depth:
                break

            K.add_generator(x)
            iterator += 1
    
        return K
    
    def verifyK(self, K: LinearSpaceF2):
        K_gen = K.generators
        for i, xi in enumerate(K_gen):
            for j, xj in enumerate(K_gen):
                print(self.prod(xi,xj))

import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Find maximal isotropic subspace K for skew-symmetric form over F2")
    parser.add_argument("-m", type=int, default=10, help="Dimension of W (default: 10)")
    parser.add_argument("-n", type=int, default=2, help="Dimension of V (default: 2)")
    parser.add_argument("-t", "--tries", type=int, default=100, help="Number of iterations (default: 100)")
    parser.add_argument("-d", "--depth", type=int, default=10, help="Search depth (random control limit, default: 10)")
    parser.add_argument("-i", "--input", type=str, default=None, help="Input file with beta definitions")
    return parser.parse_args()

def load_beta_from_file(filepath: str) -> list[np.ndarray]:
    betas = []
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    idx = 0
    while idx < len(lines):
        m, n = map(int, lines[idx].split())
        idx += 1
        beta_matrices = []
        for _ in range(n):
            matrix = []
            for _ in range(m):
                row = list(map(int, lines[idx].split()))
                matrix.append(row)
                idx += 1
            beta_matrices.append(np.array(matrix, dtype=np.uint8))
        betas.append(beta_matrices)
    
    return betas

def main():
    args = parse_args()
    m = args.m
    n = args.n
    iterations = args.tries
    depth = args.depth
    input_file = args.input

    import os
    os.makedirs("logs", exist_ok=True)

    if input_file:
        beta_list = load_beta_from_file(input_file)
        for i, B in enumerate(beta_list):
            beta = skewSymetricFormF2(B)
            log_file = f"logs/beta_{i+1}.txt"
            
            print(f"Processing beta #{i+1}/{len(beta_list)}...")
            print("========== BETA ==========")
            for j, mat in enumerate(beta.B):
                print(f"B{j}\n{mat}")
            print("==========================\n")
            
            max_dim = 0
            best_K = None
            
            for j in range(iterations):
                K = beta.findK_RSBSM(depth)
                curr_dim = K.dim()
                if curr_dim > max_dim:
                    max_dim = curr_dim
                    best_K = K
                print(f"\r--- #{j+1} --- dim: {K.dim()}", end="", flush=True)
            
            with open(log_file, mode='w') as f:
                f.write("========== BETA ==========\n")
                for j, mat in enumerate(beta.B):
                    f.write(f"B{j}\n{mat}\n")
                f.write("==========================\n\n")
                f.write(f"Max dim: {max_dim}\n")
                if best_K:
                    for k, v in enumerate(best_K.generators):
                        f.write(f"v{k+1}: {v}\n")
                f.write("=" * 40 + "\n\n")
                
                # Now run and log each attempt
                for j in range(iterations):
                    K = beta.findK_RSBSM(depth)
                    curr_dim = K.dim()
                    f.write(f"--- #{j+1} ---\n")
                    f.write(f"dim: {curr_dim}\n")
                    for k, v in enumerate(K.generators):
                        f.write(f"v{k+1}: {v}\n")
                    f.write("---\n\n")
            
            print(f"\n\nMax dim found: {max_dim}")
            if best_K:
                print(f"Best K generators:")
                for k, v in enumerate(best_K.generators):
                    print(f"v{k+1}: {v}")
            print(f"Log saved to {log_file}")
        return

    # Default: generate random beta
    B = []
    for i in range(n):
        B.append(random_skew_symmetric_F2(m))

    beta = skewSymetricFormF2(B)

    print("========== BETA ==========")
    print(f"B0\n{beta.B[0]}")
    print(f"B1\n{beta.B[1]}")
    print("==========================\n")

    max_dim = 0
    best_K = None
    
    for i in range(iterations):
        K = beta.findK_RSBSM(depth)
        curr_dim = K.dim()
        if curr_dim > max_dim:
            max_dim = curr_dim
            best_K = K
        print(f"\r--- #{i+1} --- dim: {K.dim()}", end="", flush=True)
    
    with open("log.txt", mode = 'w') as f:
        f.write("========== BETA ==========\n")
        f.write(f"B0\n{beta.B[0]}\n")
        f.write(f"B1\n{beta.B[1]}\n")
        f.write("==========================\n\n")
        f.write(f"Max dim: {max_dim}\n")
        if best_K:
            for j, v in enumerate(best_K.generators):
                f.write(f"v{j+1}: {v}\n")
        f.write("=" * 40 + "\n\n")
        
        # Now run and log each attempt
        for i in range(iterations):
            K = beta.findK_RSBSM(depth)
            curr_dim = K.dim()
            f.write(f"--- #{i+1} ---\n")
            f.write(f"dim: {curr_dim}\n")
            for j, v in enumerate(K.generators):
                f.write(f"v{j+1}: {v}\n")
            f.write("---\n\n")

    print(f"\n\nMax dim found: {max_dim}")
    if best_K:
        print(f"Best K generators:")
        for j, v in enumerate(best_K.generators):
            print(f"v{j+1}: {v}")

if __name__ == "__main__":
    main()

