import matplotlib.pyplot as plt

def plot_eigenstates(x, energies, wavefuncs, n_states=5):
    """
    Plot first n_states eigenfunctions with corresponding energy levels.
    """
    plt.figure(figsize=(8,6))
    for i in range(n_states):
        plt.plot(x, wavefuncs[:,i] + energies[i], label=f"n={i}")
    plt.xlabel("x")
    plt.ylabel("Energy + ψ(x)")
    plt.title("Eigenstates of the Schrödinger Equation")
    plt.legend()
    plt.show()
