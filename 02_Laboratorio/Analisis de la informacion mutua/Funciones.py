def graficar_probabilidades(probs, 
                            titulo="Amplitudes de probabilidad", 
                            tipo_etiqueta="bin", 
                            etiquetas_personalizadas=None,
                            rango=None,
                            ys_exitosos=None):

    probs = np.array(probs, dtype=float)
    if not np.isclose(np.sum(probs), 1):
        probs = probs / np.sum(probs)

    n_total = len(probs)
    x_total = np.arange(n_total)

    # --- Generación de etiquetas 
    if tipo_etiqueta == "num":
        etiquetas_total = [str(i) for i in range(n_total)]

    elif tipo_etiqueta == "bin":
        k = max(1, math.ceil(math.log2(n_total)))
        etiquetas_total = [format(i, f"0{k}b") for i in range(n_total)]

    else:
        raise ValueError("tipo_etiqueta debe ser 'num' o 'bin'")

    # --- Aplicar rango ---
    if rango is not None:
        inicio, fin = rango
        if not (0 <= inicio < fin <= n_total):
            raise ValueError(f"Rango inválido: debe estar dentro de [0, {n_total})")
        probs = probs[inicio:fin]
        x = np.arange(inicio, fin)
        etiquetas = etiquetas_total[inicio:fin]
    else:
        x = x_total
        etiquetas = etiquetas_total

    n = len(probs)

    # --- Colores personalizados si ys_exitosos está definido ---
    if ys_exitosos is not None:
        colors = ["red" if i in ys_exitosos else "steelblue" for i in x]
    else:
        colors = "steelblue"

    # --- Ajustes automáticos ---
    ancho_figura = max(8, n * 0.25)
    rotacion_etiquetas = 90 if n > 16 else 45 if n > 8 else 0
    if n > 512:
        fontsize_etiquetas = 8
    elif n > 128:
        fontsize_etiquetas = 10
    elif n > 32:
        fontsize_etiquetas = 12
    else:
        fontsize_etiquetas = 14

    # --- Graficar ---
    plt.figure(figsize=(ancho_figura, 5))
    plt.bar(x, probs, width=0.8, color=colors)
    plt.xticks(x, etiquetas, rotation=rotacion_etiquetas, fontsize=fontsize_etiquetas)
    plt.ylabel("Probabilidad", fontsize=14)
    plt.title(titulo, fontsize=16)
    plt.ylim(0, max(probs) * 1.1)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

def probs_teoricas_QPE_l(r, l):
    Q = 2**l  # tamaño del registro ancilla
    q = Q // r
    c = Q % r
    Q0 = r * q
    
    probs = np.zeros(Q, dtype=float)
    
    for y in range(Q):
        if (r * y) % Q != 0:  # caso general
            num = (
                c * np.sin(np.pi * r * y / Q * (Q0 / r + 1))**2
                + (r - c) * np.sin(np.pi * r * y / Q * (Q0 / r))**2
            )
            den = (Q**2) * (np.sin(np.pi * r * y / Q)**2)
            probs[y] = num / den if den != 0 else 0.0
        else:  # caso múltiplo exacto
            probs[y] = (
                c * (Q0 + r)**2 + (r - c) * (Q0**2)
            ) / (Q**2 * r**2)
    
    return probs

def U_swap(n, a, b):

    if not (0 <= a < n) or not (0 <= b < n):
        raise ValueError("Los índices de qubits deben cumplir 0 <= a,b < n.")
    if a == b:
        raise ValueError("a y b deben ser qubits distintos.")

    dim = 1 << n  # 2**n
    U = np.zeros((dim, dim), dtype=complex)

    # Convertimos de índice físico a posición en el entero (LSB indexing)
    j1 = n - 1 - a
    j2 = n - 1 - b

    for i in range(dim):
        # extraer bit j1 y j2
        b1 = (i >> j1) & 1
        b2 = (i >> j2) & 1

        if b1 == b2:
            j = i  # no cambia nada si los bits son iguales
        else:
            # flip de ambos bits
            j = i ^ (1 << j1) ^ (1 << j2)

        U[j, i] = 1.0

    return U



def U_swap_circular(l):
    U = np.eye(2**l, dtype=complex)
    for k in range(l-1):
        U = U_swap(l, k, k+1) @ U
    return U


def U_swap_circular_k(l, k):
    k = k % l  # rotaciones equivalentes mod l
    U = np.eye(2**l, dtype=complex)
    if k == 0:
        return U
    U1 = U_swap_circular(l)
    for _ in range(k):
        U = U1 @ U
    return U


def U_swap_all_pairs(n):
    """
    Aplica SWAP en todas las parejas de qubits adyacentes:
    (0<->1), (2<->3), (4<->5), ...
    """
    if n % 2 != 0:
        raise ValueError("Para intercambiar todas las parejas, n debe ser par.")
    
    U = np.eye(2**n, dtype=complex)
    for k in range(0, n, 2):
        U = U_swap_adjacent(n, k) @ U
    return U

def U_swap_adjacent(n, k):
    if not (0 <= k < n-1):
        raise ValueError("k debe cumplir 0 <= k < n-1.")

    dim = 1 << n  # 2**n
    U = np.zeros((dim, dim), dtype=complex)

    # En esta convención, el bit físico k (MSB indexing) está en posición LSB j1 = n-1-k.
    j1 = n - 1 - k
    j2 = n - 1 - (k + 1)

    for i in range(dim):
        # extrae bits
        b1 = (i >> j1) & 1
        b2 = (i >> j2) & 1
        if b1 == b2:
            j = i
        else:
            # flip de ambos bits
            j = i ^ (1 << j1) ^ (1 << j2)
        U[j, i] = 1.0
    # return U

def y_que_recuperan_r(r_real, l):
    Q = 2**l
    ys_exitosos = []

    for y in range(Q):
        x = y / Q
        convergentes = fracciones_continuas(x, M=Q)

        for (p, q) in convergentes:
            if q != 0:
                approx = Fraction(p, q).limit_denominator()
                if approx.denominator == r_real:
                    ys_exitosos.append(y)
                    break  

    return ys_exitosos

def fracciones_continuas(x, M):
    convergentes = []
    a0 = int(np.floor(x))
    p_menos2, q_menos2 = 0, 1
    p_menos1, q_menos1 = 1, 0
    
    xi = x
    k = 0
    while True:
        ak = int(np.floor(xi))
        pk = ak * p_menos1 + p_menos2
        qk = ak * q_menos1 + q_menos2
        
        convergentes.append((pk, qk))
        
        if qk > M or abs(xi - ak) < 1e-12:
            break
        
        p_menos2, q_menos2 = p_menos1, q_menos1
        p_menos1, q_menos1 = pk, qk
        xi = 1/(xi - ak)
        k += 1
    
    return convergentes

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Statevector, partial_trace
import numpy as np
from math import pi, gcd

def U_mod(a, N, n):
    """Construye la matriz unitaria modular U|y> = |a*y mod N>."""
    dim = 2**n
    if dim < N:
        raise ValueError("2**n debe ser >= N")
    if gcd(a, N) != 1:
        raise ValueError("a y N deben ser coprimos")

    U = np.zeros((dim, dim), dtype=complex)
    for y in range(N):
        U[(a * y) % N, y] = 1.0
    for y in range(N, dim):
        U[y, y] = 1.0
    return U


def QPE_qiskit(l, n, a, N):
    """
    Simula el algoritmo de Quantum Phase Estimation (QPE)
    con QFT inversa manual, usando una unitaria modular U_mod(a, N, n).

    Parámetros:
        l (int): número de qubits en el registro ancilla
        n (int): número de qubits en el registro principal
        a (int): base modular (coprima con N)
        N (int): módulo usado en U_mod

    Devuelve:
        rho_ancilla (np.ndarray): matriz densidad del registro ancilla
    """

    ancilla = QuantumRegister(l, "ancilla")
    principal = QuantumRegister(n, "principal")
    c_ancilla = ClassicalRegister(l, "c_ancilla")
    qc = QuantumCircuit(ancilla, principal, c_ancilla)

    qc.h(ancilla)             # superposición uniforme en ancilla
    qc.x(principal[0])        # inicializa |1> en el registro principal
    qc.barrier()

    U = U_mod(a, N, n)
    assert np.allclose(U.conj().T @ U, np.eye(U.shape[0])), "La matriz U no es unitaria"

    # Cadena de controlled-U^(2^j)
    for j in range(l):
        U_power = np.linalg.matrix_power(U, 2**j)
        CU_gate = UnitaryGate(U_power, label=f"U^{2**j}").control(1)
        qc.append(CU_gate, [ancilla[j]] + list(principal))

    qc.barrier()

    
    # Paso 1: swap del orden de qubits
    for i in range(l // 2):
        qc.swap(ancilla[i], ancilla[l - i - 1])

    # Paso 2: rotaciones controladas + Hadamards
    for j in range(l):
        for k in range(j):
            qc.cp(-pi / (2 ** (j - k)), ancilla[k], ancilla[j])
        qc.h(ancilla[j])

    qc.barrier()


    qc.save_density_matrix(qubits=ancilla)

    # (Opcional) medir el registro ancilla
    qc.measure(ancilla, c_ancilla)

    # ------------------------------------------------------------
    # 6. Simulación
    # ------------------------------------------------------------
    simulator = AerSimulator()
    qc_t = transpile(qc, simulator)
    job = simulator.run(qc_t)
    result = job.result()

    rho_ancilla = result.data()['density_matrix']
    return rho_ancilla


def probabilidad_exito_QPE(r, l):
    """
    Devuelve la probabilidad total de éxito del algoritmo QPE
    para periodo r y registro de l qubits.
    """

    probs = probs_teoricas_QPE_l(r, l)
    ys_exitosos = y_que_recuperan_r(r, l)

    prob_exito = sum(probs[y] for y in ys_exitosos)

    return prob_exito

def graficar_pendientes_swaps_circulares(r, l, n_puntos=30):

    ps = np.linspace(0, 1, n_puntos)

    probs = probs_teoricas_QPE_l(r, l)
    ys_exitosos = y_que_recuperan_r(r, l)

    pendientes = []

    for k in range(1, l):

        prob_exito_vs_p = []
        U = U_swap_circular_k(l, k)

        for p in ps:
            errprobs = (1 - p) * probs + p * (U @ probs)
            prob_exito = sum(errprobs[y] for y in ys_exitosos)
            prob_exito_vs_p.append(prob_exito)

        pendiente, _ = np.polyfit(ps, prob_exito_vs_p, 1)
        pendientes.append(pendiente)

    ks = np.arange(1, l)

    plt.figure(figsize=(8,5))
    plt.bar(ks, pendientes)
    plt.xlabel("Rotación circular k")
    plt.ylabel("Pendiente")
    plt.title(f"Pendientes SWAP circular (r={r}, l={l})")
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()


def graficar_pendientes_swap_fijo(r, l, qubit_fijo, n_puntos=30):

    ps = np.linspace(0, 1, n_puntos)

    probs = probs_teoricas_QPE_l(r, l)
    ys_exitosos = y_que_recuperan_r(r, l)

    pendientes = []
    qubits = []

    for j in range(l):

        if j == qubit_fijo:
            continue

        prob_exito_vs_p = []
        U = U_swap(l, qubit_fijo, j)

        for p in ps:
            errprobs = (1 - p) * probs + p * (U @ probs)
            prob_exito = sum(errprobs[y] for y in ys_exitosos)
            prob_exito_vs_p.append(prob_exito)

        pendiente, _ = np.polyfit(ps, prob_exito_vs_p, 1)

        pendientes.append(pendiente)
        qubits.append(j)

    plt.figure(figsize=(8,5))
    plt.bar(qubits, pendientes)
    plt.xlabel("Qubit intercambiado")
    plt.ylabel("Pendiente")
    plt.title(f"Pendientes SWAP({qubit_fijo}, j) (r={r}, l={l})")
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

def estimar_r(a, N, l):
    r_real = orden_modulo(a, N)
    probs = probs_teoricas_QPE_l(r_real, l)
    y, Q = medir_y_from_probs(probs, n_samples=1)
    x = y / Q

    # convergentes
    convergentes = fracciones_continuas(x, M=Q)

    for (p, q) in convergentes:
        if q != 0:
            approx = Fraction(p, q).limit_denominator()
            if approx.denominator == r_real:
                return {
                    "a": a, "N": N,
                    "orden_real": r_real,
                    "y": y, "x": x,
                    "resultado": approx,
                    "exito": True
                }

    return {
        "a": a, "N": N,
        "orden_real": r_real,
        "y": y, "x": x,
        "resultado": convergentes[-1],
        "exito": False
    }

def y_que_recuperan_r(r_real, l):
    Q = 2**l
    ys_exitosos = []

    for y in range(Q):
        x = y / Q
        convergentes = fracciones_continuas(x, M=Q)

        for (p, q) in convergentes:
            if q != 0:
                approx = Fraction(p, q).limit_denominator()
                if approx.denominator == r_real:
                    ys_exitosos.append(y)
                    break  

    return ys_exitosos



def probabilidad_exito_teorica(probs, ys_exitosos):
    prob_total = sum(probs[y] for y in ys_exitosos)
    return prob_total

def graficar_probabilidades_comparadas(probs1, probs2,
                                       titulo="Comparación de distribuciones de probabilidad",
                                       tipo_etiqueta="bin",
                                       etiquetas_personalizadas=None,
                                       rango=None,
                                       ys_exitosos=None,
                                       label1="Distribución 1",
                                       label2="Distribución 2"):
    # --- Normalizar probabilidades ---
    probs1 = np.array(probs1, dtype=float)
    probs2 = np.array(probs2, dtype=float)

    if not np.isclose(np.sum(probs1), 1):
        probs1 = probs1 / np.sum(probs1)
    if not np.isclose(np.sum(probs2), 1):
        probs2 = probs2 / np.sum(probs2)

    if len(probs1) != len(probs2):
        raise ValueError("Ambos vectores deben tener el mismo tamaño.")

    n_total = len(probs1)
    x_total = np.arange(n_total)

    # --- Generación de etiquetas ---
    if tipo_etiqueta == "num":
        etiquetas_total = [str(i) for i in range(n_total)]
    elif tipo_etiqueta == "bin":
        k = max(1, math.ceil(math.log2(n_total)))
        etiquetas_total = [format(i, f"0{k}b") for i in range(n_total)]
    else:
        raise ValueError("tipo_etiqueta debe ser 'num' o 'bin'")

    # --- Aplicar rango ---
    if rango is not None:
        inicio, fin = rango
        if not (0 <= inicio < fin <= n_total):
            raise ValueError(f"Rango inválido: debe estar dentro de [0, {n_total})")
        probs1 = probs1[inicio:fin]
        probs2 = probs2[inicio:fin]
        x = np.arange(inicio, fin)
        etiquetas = etiquetas_total[inicio:fin]
    else:
        x = x_total
        etiquetas = etiquetas_total

    n = len(probs1)

    # --- Colores si hay ys_exitosos ---
    if ys_exitosos is not None:
        colors1 = ["red" if i in ys_exitosos else "steelblue" for i in x]
        colors2 = ["orange" if i in ys_exitosos else "green" for i in x]
    else:
        colors1 = "steelblue"
        colors2 = "green"

    # --- Ajustes automáticos ---
    ancho_figura = max(8, n * 0.25)
    rotacion_etiquetas = 90 if n > 16 else 45 if n > 8 else 0
    fontsize_etiquetas = 8 if n > 512 else 10 if n > 128 else 12 if n > 32 else 14

    # --- Gráfica ---
    plt.figure(figsize=(ancho_figura, 6))
    ancho_barra = 0.4
    plt.bar(x - ancho_barra/2, probs1, width=ancho_barra, color=colors1, label=label1)
    plt.bar(x + ancho_barra/2, probs2, width=ancho_barra, color=colors2, label=label2)

    plt.xticks(x, etiquetas, rotation=rotacion_etiquetas, fontsize=fontsize_etiquetas)
    plt.ylabel("Probabilidad", fontsize=14)
    plt.title(titulo, fontsize=16)
    plt.ylim(0, max(max(probs1), max(probs2)) * 1.1)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()

def qft_matrix(k):  # QFT de k qubits
    N = 2**k
    omega = np.exp(2j * np.pi / N)
    QFT = np.array([[omega**(j*k2) for k2 in range(N)] for j in range(N)], dtype=complex)
    return QFT / np.sqrt(N)

def qft_operator(l, n):  # QFT en ancilla
    QFT = qft_matrix(l)
    I = np.eye(2**n)
    return np.kron(QFT, I)

def inverse_qft_operator(l, n):  # QFT inversa en ancilla
    QFT = qft_operator(l, n)
    return QFT.conj().T

def controlled_unitary(U, l, n):  # Control ancilla -> U^(2^(l-k)) en main
    C = np.eye(2**(l + n), dtype=complex)
    for k in range(l, 0, -1):
        Id1 = np.eye(2**(k-1)) if k > 1 else 1
        Id2 = np.eye(2**(l - k)) if k < l else 1

        P0 = np.kron(np.kron(Id1, np.array([[1, 0], [0, 0]])), Id2)
        P1 = np.kron(np.kron(Id1, np.array([[0, 0], [0, 1]])), Id2)

        C_U = np.kron(P0, np.eye(2**n)) + np.kron(P1, np.linalg.matrix_power(U, 2**(l - k)))
        C = C_U @ C
    return C

def initialize_state(l, n, vector_index=1):  # |0..0> ancilla ⊗ |vector_index> main
    ancilla = np.zeros(2**l, dtype=complex); ancilla[0] = 1.0
    main = np.zeros(2**n, dtype=complex); main[vector_index] = 1.0
    return np.kron(ancilla, main)

def qpe(U, l, n, vector_index=1):  # Quantum Phase Estimation
    u = initialize_state(l, n, vector_index)
    QFT = qft_operator(l, n)
    QFT_dag = inverse_qft_operator(l, n)
    C = controlled_unitary(U, l, n)
    return QFT_dag @ C @ QFT @ u


def orden_modulo(a, N):
    # Reducimos a módulo N
    a = a % N
    
    # Comprobamos que sean coprimos
    if np.gcd(a, N) != 1:
        raise ValueError("a y N no son coprimos, el orden no está definido.")
    
    # Buscamos el menor r > 0 tal que a^r ≡ 1 (mod N)
    r = 1
    valor = a % N
    while valor != 1:
        valor = (valor * a) % N
        r += 1
        # Por seguridad, límite en caso de ciclo inesperado
        if r > N:
            raise RuntimeError("No se encontró orden en rango esperado.")
    
    return r

import numpy as np
from math import gcd

def U_mod(a, N, n):
    dim = 2**n
    if dim < N:
        raise ValueError("2**n debe ser >= N")
    if gcd(a, N) != 1:
        raise ValueError("a y N deben ser coprimos")

    U = np.zeros((dim, dim), dtype=complex)
    for y in range(N):
        U[(a * y) % N, y] = 1.0
    for y in range(N, dim):
        U[y, y] = 1.0
    return U

def amplitudes_sistema(estado, n1, k, n2):
    total_qubits = n1 + k + n2
    if len(estado) != 2**total_qubits:
        raise ValueError("Dimensiones no coinciden con el tamaño del vector de estado")

    probs_central = np.zeros(2**k, dtype=float)
    for idx, amp in enumerate(estado):
        bits = format(idx, f"0{total_qubits}b") if total_qubits > 0 else ""
        central_bits = bits[n1 : n1 + k] if k > 0 else ""
        central_index = int(central_bits, 2) if k > 0 else 0
        probs_central[central_index] += abs(amp)**2
    return probs_central

import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.linalg import eigvals

from scipy.linalg import sqrtm

def fidelity(rho, sigma):
    sqrt_rho = sqrtm(rho)
    inner = sqrt_rho @ sigma @ sqrt_rho
    F = np.trace(sqrtm(inner))
    return np.real(F)

def estado_puro(psi):
    psi = np.asarray(psi, dtype=complex).reshape(-1, 1)  # fuerza vector columna
    return psi @ psi.conj().T


def operacion_unitaria(U,rho): #Aplica una operación unitaria a una matriz de densidad
    return U@rho@U.conj().T

def partial_trace_2systems(rho, n, m, keep="A"):
    """
    Supone orden de base |i>_A ⊗ |j>_B con indexación combinada i*m + j.
    """
    rho = np.asarray(rho, dtype=complex)
    assert rho.shape == (n*m, n*m), "ρ debe ser de tamaño (n*m)×(n*m)"

    rho_t = rho.reshape(n, m, n, m)  # ρ[i, j, k, l]

    if keep in ("A", 0):
        # ρ_A[i,k] = sum_j ρ[i, j, k, j]
        return np.trace(rho_t, axis1=1, axis2=3)
    elif keep in ("B", 1):
        # ρ_B[j,l] = sum_i ρ[i, j, i, l]
        return np.trace(rho_t, axis1=0, axis2=2)
    else:
        raise ValueError("keep debe ser 'A','B', 0 o 1")



def U_swap_adjacent(n, k):
    if not (0 <= k < n-1):
        raise ValueError("k debe cumplir 0 <= k < n-1.")

    dim = 1 << n  # 2**n
    U = np.zeros((dim, dim), dtype=complex)

    # En esta convención, el bit físico k (MSB indexing) está en posición LSB j1 = n-1-k.
    j1 = n - 1 - k
    j2 = n - 1 - (k + 1)

    for i in range(dim):
        # extrae bits
        b1 = (i >> j1) & 1
        b2 = (i >> j2) & 1
        if b1 == b2:
            j = i
        else:
            # flip de ambos bits
            j = i ^ (1 << j1) ^ (1 << j2)
        U[j, i] = 1.0
    return U

import numpy as np

def U_swap(n, a, b):

    if not (0 <= a < n) or not (0 <= b < n):
        raise ValueError("Los índices de qubits deben cumplir 0 <= a,b < n.")
    if a == b:
        raise ValueError("a y b deben ser qubits distintos.")

    dim = 1 << n  # 2**n
    U = np.zeros((dim, dim), dtype=complex)

    # Convertimos de índice físico a posición en el entero (LSB indexing)
    j1 = n - 1 - a
    j2 = n - 1 - b

    for i in range(dim):
        # extraer bit j1 y j2
        b1 = (i >> j1) & 1
        b2 = (i >> j2) & 1

        if b1 == b2:
            j = i  # no cambia nada si los bits son iguales
        else:
            # flip de ambos bits
            j = i ^ (1 << j1) ^ (1 << j2)

        U[j, i] = 1.0

    return U



def U_swap_circular(l):
    U = np.eye(2**l, dtype=complex)
    for k in range(l-1):
        U = U_swap(l, k, k+1) @ U
    return U


def U_swap_circular_k(l, k):
    k = k % l  # rotaciones equivalentes mod l
    U = np.eye(2**l, dtype=complex)
    if k == 0:
        return U
    U1 = U_swap_circular(l)
    for _ in range(k):
        U = U1 @ U
    return U


def U_swap_all_pairs(n):
    """
    Aplica SWAP en todas las parejas de qubits adyacentes:
    (0<->1), (2<->3), (4<->5), ...
    """
    if n % 2 != 0:
        raise ValueError("Para intercambiar todas las parejas, n debe ser par.")
    
    U = np.eye(2**n, dtype=complex)
    for k in range(0, n, 2):
        U = U_swap_adjacent(n, k) @ U
    return U

def swap_channel_adjacent(rho, n, k, p):
    if not (0 <= p <= 1):
        raise ValueError("p debe estar en [0,1].")

    dim = 1 << n
    rho = np.asarray(rho, dtype=complex)
    if rho.shape != (dim, dim):
        raise ValueError(f"rho debe ser de tamaño {(dim, dim)} para n={n}.")

    U = U_swap_adjacent(n, k)
    rho_sw = operacion_unitaria(U, rho)
    return (1 - p) * rho + p * rho_sw


def plot_probabilities(rho):
    rho = np.asarray(rho, dtype=complex)
    dim = rho.shape[0]

    # verificar cuadratura
    if rho.shape[0] != rho.shape[1]:
        raise ValueError("rho debe ser cuadrada.")

    # verificar que dim sea potencia de 2
    n = int(np.log2(dim))
    if 2**n != dim:
        raise ValueError("La dimensión de rho no es potencia de 2.")

    # probabilidades = diagonal
    probs = np.real(np.diag(rho))
    probs = np.maximum(probs, 0)   # corrige posibles errores numéricos
    probs = probs / probs.sum()    # normaliza

    # etiquetas binarias
    labels = [format(i, f'0{n}b') for i in range(dim)]

    # gráfico
    plt.figure(figsize=(1.5*n, 4))
    plt.bar(labels, probs)
    plt.xlabel("Estados (base computacional)")
    plt.ylabel("Probabilidad")
    plt.title(f"Distribución de medición ({n} qubits)")
    plt.xticks(rotation=90)  # <-- etiquetas en vertical
    plt.show()

    return probs, labels


def graficar_probabilidades_comparadas_rho(rho1, rho2,
                                           titulo="Comparación de distribuciones (matrices de densidad)",
                                           etiquetas=None,
                                           label1="ρ1", label2="ρ2"):

    # Asegurar arrays complejos
    rho1 = np.asarray(rho1, dtype=complex)
    rho2 = np.asarray(rho2, dtype=complex)

    # Dimensiones
    if rho1.shape[0] != rho1.shape[1] or rho2.shape[0] != rho2.shape[1]:
        raise ValueError("Las matrices deben ser cuadradas")
    if rho1.shape != rho2.shape:
        raise ValueError("Ambas matrices deben tener la misma dimensión")

    dim = rho1.shape[0]
    n = int(np.log2(dim))
    if 2**n != dim:
        raise ValueError("La dimensión no corresponde a un número entero de qubits")

    # Probabilidades = diagonal normalizada
    probs1 = np.real(np.diag(rho1))
    probs1 = np.maximum(probs1, 0)
    probs1 = probs1 / probs1.sum()

    probs2 = np.real(np.diag(rho2))
    probs2 = np.maximum(probs2, 0)
    probs2 = probs2 / probs2.sum()

    # Etiquetas binarias
    if etiquetas is None:
        etiquetas = [format(i, f"0{n}b") for i in range(dim)]

    # Gráfico comparado
    x = np.arange(dim)
    ancho = 0.4
    plt.figure(figsize=(max(15, 6), 4))
    plt.bar(x - ancho/2, probs1, width=ancho, label=label1)
    plt.bar(x + ancho/2, probs2, width=ancho, label=label2)

    plt.xticks(x, etiquetas, rotation=90)
    plt.ylabel("Probabilidad")
    plt.title(titulo)
    plt.ylim(0, max(max(probs1), max(probs2)) * 1.1)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return probs1, probs2, etiquetas

import numpy as np

def diagonal_probs(rho, normalize=True):
    rho = np.asarray(rho, dtype=complex)
    if rho.shape[0] != rho.shape[1]:
        raise ValueError("rho debe ser una matriz cuadrada")

    diag = np.real(np.diag(rho))
    diag = np.maximum(diag, 0)  # elimina pequeños negativos numéricos

    if normalize:
        diag = diag / diag.sum()

    return diag

ket0 = np.array([[1],
                 [0]], dtype=complex)

ket1 = np.array([[0],
                 [1]], dtype=complex)

X = np.array([[0, 1],
              [1, 0]], dtype=complex)

Y = np.array([[0, -1j],
              [1j, 0]], dtype=complex)

Z = np.array([[1, 0],
              [0, -1]], dtype=complex)

I = np.array([[1, 0],
              [0, 1]], dtype=complex)

H = (1/np.sqrt(2)) * np.array([[1, 1],
                                [1, -1]], dtype=complex)



import numpy as np
import math
import matplotlib.pyplot as plt

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Statevector, partial_trace
import numpy as np
from math import pi, gcd

def U_mod(a, N, n):
    """Construye la matriz unitaria modular U|y> = |a*y mod N>."""
    dim = 2**n
    if dim < N:
        raise ValueError("2**n debe ser >= N")
    if gcd(a, N) != 1:
        raise ValueError("a y N deben ser coprimos")

    U = np.zeros((dim, dim), dtype=complex)
    for y in range(N):
        U[(a * y) % N, y] = 1.0
    for y in range(N, dim):
        U[y, y] = 1.0
    return U


def QPE_qiskit(l, n, a, N):
    """
    Simula el algoritmo de Quantum Phase Estimation (QPE)
    con QFT inversa manual, usando una unitaria modular U_mod(a, N, n).

    Parámetros:
        l (int): número de qubits en el registro ancilla
        n (int): número de qubits en el registro principal
        a (int): base modular (coprima con N)
        N (int): módulo usado en U_mod

    Devuelve:
        rho_ancilla (np.ndarray): matriz densidad del registro ancilla
    """

    # ------------------------------------------------------------
    # 1. Registros cuánticos
    # ------------------------------------------------------------
    ancilla = QuantumRegister(l, "ancilla")
    principal = QuantumRegister(n, "principal")
    c_ancilla = ClassicalRegister(l, "c_ancilla")
    qc = QuantumCircuit(ancilla, principal, c_ancilla)

    # ------------------------------------------------------------
    # 2. Inicialización
    # ------------------------------------------------------------
    qc.h(ancilla)             # superposición uniforme en ancilla
    qc.x(principal[0])        # inicializa |1> en el registro principal
    qc.barrier()

    # ------------------------------------------------------------
    # 3. Construcción de la unitaria U(a, N, n)
    # ------------------------------------------------------------
    U = U_mod(a, N, n)
    assert np.allclose(U.conj().T @ U, np.eye(U.shape[0])), "La matriz U no es unitaria"

    # Cadena de controlled-U^(2^j)
    for j in range(l):
        U_power = np.linalg.matrix_power(U, 2**j)
        CU_gate = UnitaryGate(U_power, label=f"U^{2**j}").control(1)
        qc.append(CU_gate, [ancilla[j]] + list(principal))

    qc.barrier()

    # ------------------------------------------------------------
    # 4. QFT inversa manual
    # ------------------------------------------------------------
    # Paso 1: swap del orden de qubits
    for i in range(l // 2):
        qc.swap(ancilla[i], ancilla[l - i - 1])

    # Paso 2: rotaciones controladas + Hadamards
    for j in range(l):
        for k in range(j):
            qc.cp(-pi / (2 ** (j - k)), ancilla[k], ancilla[j])
        qc.h(ancilla[j])

    qc.barrier()

    # ------------------------------------------------------------
    # 5. Guardar matriz densidad del registro ancilla
    # ------------------------------------------------------------
    qc.save_density_matrix(qubits=ancilla)

    # (Opcional) medir el registro ancilla
    qc.measure(ancilla, c_ancilla)

    # ------------------------------------------------------------
    # 6. Simulación
    # ------------------------------------------------------------
    simulator = AerSimulator()
    qc_t = transpile(qc, simulator)
    job = simulator.run(qc_t)
    result = job.result()

    rho_ancilla = result.data()['density_matrix']
    return rho_ancilla

def QPE_qiskit_vec(l, n, a, N):
    """
    Simula el algoritmo de Quantum Phase Estimation (QPE)
    con QFT inversa manual, usando una unitaria modular U_mod(a, N, n).
    
    Ahora devuelve el STATEVECTOR COMPLETO.
    """

    # ------------------------------------------------------------
    # 1. Registros cuánticos
    # ------------------------------------------------------------
    ancilla = QuantumRegister(l, "ancilla")
    principal = QuantumRegister(n, "principal")
    qc = QuantumCircuit(ancilla, principal)

    # ------------------------------------------------------------
    # 2. Inicialización
    # ------------------------------------------------------------
    qc.h(ancilla)             # superposición uniforme en ancilla
    qc.x(principal[0])        # inicializa |1> en el registro principal
    qc.barrier()

    # ------------------------------------------------------------
    # 3. Construcción de la unitaria U(a, N, n)
    # ------------------------------------------------------------
    U = U_mod(a, N, n)
    assert np.allclose(U.conj().T @ U, np.eye(U.shape[0])), "La matriz U no es unitaria"

    # Cadena de controlled-U^(2^j)
    for j in range(l):
        U_power = np.linalg.matrix_power(U, 2**j)
        CU_gate = UnitaryGate(U_power, label=f"U^{2**j}").control(1)
        qc.append(CU_gate, [ancilla[j]] + list(principal))

    qc.barrier()

    # ------------------------------------------------------------
    # 4. QFT inversa manual en el registro ancilla
    # ------------------------------------------------------------
    # Paso 1: swap del orden de qubits
    for i in range(l // 2):
        qc.swap(ancilla[i], ancilla[l - i - 1])

    # Paso 2: rotaciones controladas + Hadamards
    for j in range(l):
        for k in range(j):
            qc.cp(-pi / (2 ** (j - k)), ancilla[k], ancilla[j])
        qc.h(ancilla[j])

    qc.barrier()

    # ------------------------------------------------------------
    # 5. Guardar el STATEVECTOR completo
    # ------------------------------------------------------------
    qc.save_statevector()

    # ------------------------------------------------------------
    # 6. Simulación
    # ------------------------------------------------------------
    simulator = AerSimulator(method="statevector")
    qc_t = transpile(qc, simulator)
    job = simulator.run(qc_t)
    result = job.result()

    statevector = result.get_statevector()

    return statevector

#Definiremos algunas funciones extra, donde el nucleo sera la función probs_teoricas_QPE_l

import numpy as np
from fractions import Fraction


def medir_y_from_probs(probs, n_samples=1):  #Mide alguno de los valores usando una distribución de probabilidad pre dada 
    Q = len(probs)  # dimensión del sistema
    probs = np.array(probs, dtype=float)
    probs /= np.sum(probs)  # normalizamos, pues la suma de las probabilidades deben de ser 1
    
    results = np.random.choice(Q, size=n_samples, p=probs) #medimos
    
    if n_samples == 1:
        return results[0], Q
    else:
        return results, Q


#Algoritmo de fracciones continuas

def fracciones_continuas(x, M):
    convergentes = []
    a0 = int(np.floor(x))
    p_menos2, q_menos2 = 0, 1
    p_menos1, q_menos1 = 1, 0
    
    xi = x
    k = 0
    while True:
        ak = int(np.floor(xi))
        pk = ak * p_menos1 + p_menos2
        qk = ak * q_menos1 + q_menos2
        
        convergentes.append((pk, qk))
        
        if qk > M or abs(xi - ak) < 1e-12:
            break
        
        p_menos2, q_menos2 = p_menos1, q_menos1
        p_menos1, q_menos1 = pk, qk
        xi = 1/(xi - ak)
        k += 1
    
    return convergentes




def joint_distribution(P, U, p):
 
    
    dim = len(P)
    joint = np.zeros((dim, dim))


    perm = np.argmax(U, axis=0)

    for x in range(dim):
        y_swap = perm[x]

        joint[x, x] += (1 - p) * P[x]
        joint[x, y_swap] += p * P[x]

    return joint

def mutual_information(Pxy):

    Px = np.sum(Pxy, axis=1)
    Py = np.sum(Pxy, axis=0)

    I = 0.0

    for x in range(len(Px)):
        for y in range(len(Py)):

            pxy = Pxy[x, y]

            if pxy > 0:
                I += pxy * np.log2(pxy / (Px[x] * Py[y]))

    return I

def plot_swap_analysis(r, l, a):
    
    P = probs_teoricas_QPE_l(r, l)
    ys_exitosos = y_que_recuperan_r(r, l)

    ps_I = np.linspace(0, 1, 100)
    ps_P = np.linspace(0, 1, 30)

    fig, axs = plt.subplots(1, 2, figsize=(14,5))

    for b in range(l):
        if b == a:
            continue

        U = U_swap(l, a, b)

        # ---- Información mutua ----
        axs[0].plot(
            ps_I,
            [mutual_information(joint_distribution(P, U, p)) for p in ps_I],
            label=f"swap({a},{b})"
        )

        # ---- Probabilidad de éxito ----
        axs[1].plot(
            ps_P,
            [
                sum(((1-p)*P + p*(U @ P))[y] for y in ys_exitosos)
                for p in ps_P
            ],
            marker='o',
            label=f"swap({a},{b})"
        )

    # Configuración gráfica
    axs[0].set(
        xlabel="Probabilidad de error p",
        ylabel="Información mutua I(X;Y)",
        title=f"Información mutua (qubit fijo {a})"
    )
    axs[0].legend()
    axs[0].grid(True)

    axs[1].set(
        xlabel="Probabilidad de error p",
        ylabel="Probabilidad de éxito",
        title=f"Probabilidad de éxito (qubit fijo {a})"
    )
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()

def plot_circular_swap_analysis(r, l):

    P = probs_teoricas_QPE_l(r, l)
    ys_exitosos = y_que_recuperan_r(r, l)

    ps = np.linspace(0, 1, 100)

    fig, axs = plt.subplots(1, 2, figsize=(14,5))

    for k in range(1, l):

        U = U_swap_circular_k(l, k)

        # ---- Información mutua ----
        axs[0].plot(
            ps,
            [mutual_information(joint_distribution(P, U, p)) for p in ps],
            label=f"k={k}"
        )

        # ---- Probabilidad de éxito ----
        axs[1].plot(
            ps,
            [
                sum(((1-p)*P + p*(U @ P))[y] for y in ys_exitosos)
                for p in ps
            ],
            label=f"k={k}"
        )

    axs[0].set(
        xlabel="Probabilidad de error p",
        ylabel="Información mutua I(X;Y)",
        title="Información mutua (rotación circular)"
    )
    axs[0].legend()
    axs[0].grid(True)

    axs[1].set(
        xlabel="Probabilidad de error p",
        ylabel="Probabilidad de éxito",
        title="Probabilidad de éxito (rotación circular)"
    )
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()

def gaussiana(n,mu,sigma):
  x = np.arange(n)
  g = np.exp(-(x-mu)**2/(2*sigma**2))
  g = g/ np.linalg.norm(g)
  return g

def graficas_prob_exito_circular(r, l):

    P = probs_teoricas_QPE_l(r, l)
    ys_exitosos = y_que_recuperan_r(r, l)

    ps = np.linspace(0, 1, 100)

    plt.figure(figsize=(7,5))

    for k in range(1, l):

        U = U_swap_circular_k(l, k)

        probs_exito = [
            sum(((1-p)*P + p*(U @ P))[y] for y in ys_exitosos)
            for p in ps
        ]

        plt.plot(ps, probs_exito, label=f"k={k}")

    plt.xlabel("Probabilidad de error p")
    plt.ylabel("Probabilidad de éxito")
    plt.title("Probabilidad de éxito (rotación circular)")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

def graficar_informacion_mutua_circular(r, l):

    P = probs_teoricas_QPE_l(r, l)
    ps = np.linspace(0, 1, 100)

    plt.figure(figsize=(7,5))

    for k in range(1, l):

        U = U_swap_circular_k(l, k)

        infos = [
            mutual_information(joint_distribution(P, U, p))
            for p in ps
        ]

        plt.plot(ps, infos, label=f"k={k}")

    plt.xlabel("Probabilidad de error p")
    plt.ylabel("Información mutua I(X;Y)")
    plt.title("Información mutua (rotación circular)")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

def plot_pairwise_swap_analysis(r, l, fixed_qubit):

    P = probs_teoricas_QPE_l(r, l)
    ys_exitosos = y_que_recuperan_r(r, l)

    ps = np.linspace(0, 1, 100)

    fig, axs = plt.subplots(1, 2, figsize=(14,5))

    for q in range(l):

        if q == fixed_qubit:
            continue

        U = U_swap(l, fixed_qubit, q)

        # ---- Información mutua ----
        axs[0].plot(
            ps,
            [mutual_information(joint_distribution(P, U, p)) for p in ps],
            label=f"{fixed_qubit}↔{q}"
        )

        # ---- Probabilidad de éxito ----
        axs[1].plot(
            ps,
            [
                sum(((1-p)*P + p*(U @ P))[y] for y in ys_exitosos)
                for p in ps
            ],
            label=f"{fixed_qubit}↔{q}"
        )

    axs[0].set(
        xlabel="Probabilidad de error p",
        ylabel="Información mutua I(X;Y)",
        title=f"Información mutua (SWAP qubit {fixed_qubit})"
    )
    axs[0].legend()
    axs[0].grid(True)

    axs[1].set(
        xlabel="Probabilidad de error p",
        ylabel="Probabilidad de éxito",
        title=f"Probabilidad de éxito (SWAP qubit {fixed_qubit})"
    )
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()