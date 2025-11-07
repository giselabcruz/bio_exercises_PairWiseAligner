import argparse
import random
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner

BASES_ADN = "mouse_actb.fasta"
random.seed(42)

def generar_adn_aleatorio(longitud):
    return "".join(random.choice(BASES_ADN) for _ in range(longitud))

def calcular_identidad_y_huecos(secuencia1_alineada, secuencia2_alineada):
    coincidencias = sum(1 for a, b in zip(secuencia1_alineada, secuencia2_alineada) if a == b and a != "-" and b != "-")
    posiciones_alineadas = sum(1 for a, b in zip(secuencia1_alineada, secuencia2_alineada) if a != "-" and b != "-")
    identidad = (coincidencias / posiciones_alineadas) * 100 if posiciones_alineadas else 0.0
    huecos = secuencia1_alineada.count("-") + secuencia2_alineada.count("-")
    return round(identidad, 2), coincidencias, huecos

def reconstruir_con_huecos(alignment, seq_a, seq_b):
    A, B = [], []
    coords = alignment.coordinates
    for j in range(coords.shape[1] - 1):
        a0, a1 = coords[0, j], coords[0, j+1]
        b0, b1 = coords[1, j], coords[1, j+1]
        while a0 < a1 and b0 < b1:
            A.append(seq_a[a0]); B.append(seq_b[b0])
            a0 += 1; b0 += 1
        while a0 < a1:
            A.append(seq_a[a0]); B.append('-')
            a0 += 1
        while b0 < b1:
            A.append('-'); B.append(seq_b[b0])
            b0 += 1
    return "".join(A), "".join(B)

def formatear_bloques(alA, alB, etiqueta_arriba="Secuencia A", etiqueta_abajo="Secuencia B", ancho=60, off_a=0, off_b=0):
    i = 0
    out = []
    while i < len(alA):
        sA = alA[i:i+ancho]
        sB = alB[i:i+ancho]
        mid = "".join('|' if (a == b and a != '-' and b != '-') else ('.' if a != '-' and b != '-' else ' ')
                      for a,b in zip(sA,sB))
        startA = off_a + sum(1 for c in alA[:i] if c != '-') + 1
        endA   = startA + sum(1 for c in sA      if c != '-') - 1
        startB = off_b + sum(1 for c in alB[:i] if c != '-') + 1
        endB   = startB + sum(1 for c in sB      if c != '-') - 1

        out.append(f"{etiqueta_arriba:<12} {startA:>6} {sA} {endA}")
        out.append(f"{'':<12} {'':>6} {mid}")
        out.append(f"{etiqueta_abajo:<12} {startB:>6} {sB} {endB}\n")
        i += ancho
    return "\n".join(out).rstrip()

def ejecutar_alineamiento(secuencia_a, secuencia_b, descripcion, modo,
                          coincidencia=2, diferencia=-1, hueco_lineal=None,
                          apertura_hueco=-10, extension_hueco=-0.5,
                          etiqueta_arriba="Secuencia A", etiqueta_abajo="Secuencia B"):
    alineador = PairwiseAligner()
    alineador.mode = modo
    alineador.match_score = coincidencia
    alineador.mismatch_score = diferencia
    if hueco_lineal is not None:
        alineador.open_gap_score = hueco_lineal
        alineador.extend_gap_score = hueco_lineal
    else:
        alineador.open_gap_score = apertura_hueco
        alineador.extend_gap_score = extension_hueco

    mejor_alineamiento = alineador.align(secuencia_a, secuencia_b)[0]

    alA, alB = reconstruir_con_huecos(mejor_alineamiento, secuencia_a, secuencia_b)
    identidad, num_coincidencias, num_huecos = calcular_identidad_y_huecos(alA, alB)

    a_ini = mejor_alineamiento.aligned[0][0][0] if mejor_alineamiento.aligned[0].size else 0
    b_ini = mejor_alineamiento.aligned[1][0][0] if mejor_alineamiento.aligned[1].size else 0

    print("\n" + "═" * 70)
    print(f"Configuración: {descripcion}")
    print(f"Modo de alineamiento: {modo.upper()}")
    print("─" * 70)
    print(formatear_bloques(alA, alB, etiqueta_arriba=etiqueta_arriba, etiqueta_abajo=etiqueta_abajo,
                            ancho=60, off_a=a_ini, off_b=b_ini))
    print("─" * 70)
    print(f"Puntuación total         : {mejor_alineamiento.score:.3f}")
    print(f"Porcentaje de identidad  : {identidad}%")
    print(f"Número de coincidencias  : {num_coincidencias}")
    print(f"Número total de huecos   : {num_huecos}")
    print("═" * 70)

def comparar_alineamientos(secuencia_a, secuencia_b, titulo):
    print("\n" + "=" * 70)
    print(f"{titulo.center(70)}")
    print("=" * 70)
    configuraciones = [
        ("Global (afín)   m=2  mm=-1  open=-10  ext=-0.5", "global", None, -10, -0.5),
        ("Global (lineal) m=2  mm=-1  gap=-2",             "global", -2, None, None),
        ("Local (afín)    m=2  mm=-1  open=-10  ext=-0.5", "local",  None, -10, -0.5),
        ("Local (lineal)  m=2  mm=-1  gap=-2",             "local",  -2, None, None),
    ]
    for descripcion, modo, hueco_lineal, apertura_hueco, extension_hueco in configuraciones:
        ejecutar_alineamiento(
            secuencia_a, secuencia_b, descripcion, modo,
            coincidencia=2, diferencia=-1,
            hueco_lineal=hueco_lineal,
            apertura_hueco=apertura_hueco if apertura_hueco is not None else -10,
            extension_hueco=extension_hueco if extension_hueco is not None else -0.5,
            etiqueta_arriba="Secuencia A",   # ← aquí personalizas “target”
            etiqueta_abajo="Secuencia B"     # ← aquí personalizas “query”
        )

def leer_secuencia_fasta(ruta_fichero):
    registro = next(SeqIO.parse(str(ruta_fichero), "fasta"))
    return str(registro.seq).upper().replace("U", "T")

def main():
    parser = argparse.ArgumentParser(description="Ejercicio 2: Alineamiento de secuencias con PairwiseAligner (ADN)")
    parser.add_argument("--fasta", nargs=2, metavar=("SEQ1.fasta", "SEQ2.fasta"),
                        help="Alinear dos secuencias obtenidas desde archivos FASTA")
    parser.add_argument("--lenA", type=int, default=40, help="Longitud de la primera secuencia aleatoria")
    parser.add_argument("--lenB", type=int, default=32, help="Longitud de la segunda secuencia aleatoria")
    args = parser.parse_args()

    if args.fasta:
        fasta1, fasta2 = map(Path, args.fasta)
    else:
        fasta1 = Path("./fasta/mouse_actb.fasta")
        fasta2 = Path("./fasta/lizard_actb.fasta")

    secuencia_a = leer_secuencia_fasta(fasta1)
    secuencia_b = leer_secuencia_fasta(fasta2)

    print("\nArchivos FASTA cargados correctamente:")
    print(f" - {fasta1.name} (longitud: {len(secuencia_a)})")
    print(f" - {fasta2.name} (longitud: {len(secuencia_b)})")

    comparar_alineamientos(secuencia_a, secuencia_b, titulo="(b) Secuencias obtenidas de bases de datos")

if __name__ == "__main__":
    main()