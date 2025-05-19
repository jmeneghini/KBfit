import numpy as np
import os
from scipy.linalg import sqrtm



def getCB(B):
    """
    Compute the CB matrix from the B matrix.
    (1 + i B)(1 - iB)^{-1}
    """
    I = np.eye(B.shape[0], dtype=complex)
    return (I + 1j * B) @ np.linalg.inv(I - 1j * B)

def getStilde(Ktilde):
    """
    Compute the Stilde matrix from the Ktilde matrix.
    (1 + i Ktilde)(1 - i Ktilde)^{-1}
    """
    I = np.eye(Ktilde.shape[0], dtype=complex)
    return  (I + 1j * Ktilde) @ np.linalg.inv(I - 1j * Ktilde)

def getStildeFromInv(Ktilde_inv):
    """
    Compute the Stilde_inv matrix from the Ktilde_inv matrix.
    -(1 - i Ktilde_inv)(1 + i Ktilde_inv)^{-1}
    """
    I = np.eye(Ktilde_inv.shape[0], dtype=complex)
    return  (-I + 1j * Ktilde_inv) @ np.linalg.inv(I + 1j * Ktilde_inv)

def getOmega(mu, A):
    """
    Compute Omega(mu, A) = det(A)/(det(mu^2 + A A^dagger)^{1/2})
    """
    A_dagger = np.conjugate(A.T)
    det_A = np.linalg.det(A)
    det_term = np.linalg.det(sqrtm(mu**2 * np.eye(A.shape[0], dtype=complex) + A @ A_dagger))
    return det_A / det_term

def getQuantConds(Ktilde, B, Ktilde_inv, B_from_inv):
    """
    Compute the quantization conditions for the given Ktilde, B, and Ktilde_inv matrices.
    """
    Stilde = getStilde(Ktilde)
    StildeFromInv = getStildeFromInv(Ktilde_inv)
    CB = getCB(B)
    CBfromInv = getCB(B_from_inv)
    mu = 5.0

    qc_mats = [np.eye(Stilde.shape[0], dtype=complex) + Stilde@CB,
               np.eye(Stilde.shape[0], dtype=complex) - Ktilde@B]
    qc_mats_inv = [np.conj(StildeFromInv).T + CBfromInv,
                   Ktilde_inv - B_from_inv]

    qc_evs = [sorted(list(np.linalg.eigvals(mat)), key=lambda x: np.real(x)) for mat in qc_mats]
    qc_evs_inv = [sorted(list(np.linalg.eigvals(mat)), key=lambda x: np.real(x)) for mat in qc_mats_inv]
    qc_omegas = [getOmega(mu, mat) for mat in [qc_mats[0], qc_mats_inv[0], qc_mats[1], qc_mats_inv[1]]]
    return qc_evs, qc_evs_inv, qc_omegas


def main():
    # set pwd to the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    # load in the arrays
    B = np.load('B.npy')
    B_from_inv = np.load('B_from_inv.npy')
    Ktilde = np.load('Ktilde.npy')
    Ktilde_inv = np.load('Ktilde_inv.npy')
    np.set_printoptions(precision=3)
    
    # get their Cayley transforms
    CB = getCB(B)
    CBfromInv = getCB(B_from_inv)
    Stilde = getStilde(Ktilde)
    Stilde_inv = getStildeFromInv(Ktilde_inv)
    # write them to files
    np.save('CB.npy', CB)
    np.save('Stilde.npy', Stilde)
    np.save('Stilde_inv.npy', Stilde_inv)
    
    # get the quantization conditions
    qc_evs, qc_evs_inv, qc_omegas = getQuantConds(Ktilde, B, Ktilde_inv, B_from_inv)

    # write them to files
    np.save('qc_evs.npy', qc_evs)
    np.save('qc_evs_inv.npy', qc_evs_inv)
    np.save('qc_omegas.npy', qc_omegas)
    
    return 0
    
if __name__ == "__main__":
    main()