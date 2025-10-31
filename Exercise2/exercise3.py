import argparse
import random
from pathlib import Path
from typing import Optional, Tuple, List

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices as sm

AA20 = "ACDEFGHIKLMNPQRSTVWY"
random.seed(42)


def generar_proteina_aleatoria(longitud: int) -> str:
    return "".join(random.choice(AA20) for _ in range(longitud))


def limpiar_proteina(seq: str, reemplazo: str = "X") -> str:
    seq = seq.upper()
    return "".join(c if c in AA20 else reemplazo for c in seq)


def calcular_identidad_y_huecos(secuencia1_alineada: str, secuencia2_alineada: str) -> Tuple[float, int, int]:
    coincidencias = sum(
        1 for a, b in zip(secuencia1_alineada, secuencia2_alineada)
        if a == b and a != "-" and b != "-"
    )
    posiciones_alineadas = sum(
        1 for a, b in zip(secuencia1_alineada, secuencia2_alineada)
        if a != "-" and b != "-"
    )
    identidad = (coincidencias / posiciones_alineadas) * 100 if posiciones_alineadas else 0.0
    huecos = secuencia1_alineada.count("-") + secuencia2_alineada.count("-")
    return round(identidad, 2), coincidencias, huecos


def reconstruir_con_huecos(alignment, seq_a: str, seq_b: str) -> Tuple[str, str]:
    A, B = [], []
    coords = alignment.coordinates
    for j in range(coords.shape[1] - 1):
        a0, a1 = coords[0, j], coords[0, j + 1]
        b0, b1 = coords[1, j], coords[1, j + 1]
        while a0 < a1 and b0 < b1:
            A.append(seq_a[a0]); B.append(seq_b[b0])
            a0 += 1; b0 += 1
        while a0 < a1:
            A.append(seq_a[a0]); B.append('-'); a0 += 1
        while b0 < b1:
            A.append('-'); B.append(seq_b[b0]); b0 += 1
    return "".join(A), "".join(B)


def formatear_bloques(alA: str, alB: str, etiqueta_arriba="Proteína A", etiqueta_abajo="Proteína B",
                      ancho=60, off_a=0, off_b=0) -> str:
    i = 0
    out: List[str] = []
    while i < len(alA):
        sA = alA[i:i + ancho]
        sB = alB[i:i + ancho]
        mid = "".join(
            '|' if (a == b and a != '-' and b != '-') else ('.' if a != '-' and b != '-' else ' ')
            for a, b in zip(sA, sB)
        )
        startA = off_a + sum(1 for c in alA[:i] if c != '-') + 1
        endA   = startA + sum(1 for c in sA      if c != '-') - 1
        startB = off_b + sum(1 for c in alB[:i] if c != '-') + 1
        endB   = startB + sum(1 for c in sB      if c != '-') - 1

        out.append(f"{etiqueta_arriba:<12} {startA:>6} {sA} {endA}")
        out.append(f"{'':<12} {'':>6} {mid}")
        out.append(f"{etiqueta_abajo:<12} {startB:>6} {sB} {endB}\n")
        i += ancho
    return "\n".join(out).rstrip()


def cargar_matriz(nombre_o_path: Optional[str]):
    if nombre_o_path is None:
        return sm.load("BLOSUM62")

    try:
        return sm.load(nombre_o_path)
    except Exception:
        pass

    path = Path(nombre_o_path)
    if not path.exists():
        raise FileNotFoundError(f"No se encuentra la matriz: {nombre_o_path}")

    with path.open() as fh:
        rows = [line.strip() for line in fh if line.strip() and not line.strip().startswith("#")]
    header = rows[0].split()
    alphabet = "".join(header)
    n = len(header)
    arr = [[0] * n for _ in range(n)]
    for i, line in enumerate(rows[1:1 + n]):
        parts = line.split()
        if parts[0] != header[i]:
            raise ValueError("Fila/columna no alineadas en la matriz personalizada.")
        scores = list(map(float, parts[1:1 + n]))
        for j, val in enumerate(scores):
            arr[i][j] = val

    mat = sm.Array(alphabet=alphabet, dims=2)
    for i, a in enumerate(alphabet):
        for j, b in enumerate(alphabet):
            mat[(a, b)] = arr[i][j]
    return mat

def ejecutar_alineamiento(
    secuencia_a: str,
    secuencia_b: str,
    descripcion: str,
    modo: str,
    matriz: str,
    gap_lineal: Optional[float],
    open_gap: float,
    ext_gap: float,
    etiqueta_arriba="Proteína A",
    etiqueta_abajo="Proteína B",
):
    alineador = PairwiseAligner()
    alineador.mode = modo
    subm = cargar_matriz(matriz)
    alineador.substitution_matrix = subm
    if gap_lineal is not None:
        alineador.open_gap_score = gap_lineal
        alineador.extend_gap_score = gap_lineal
    else:
        alineador.open_gap_score = open_gap
        alineador.extend_gap_score = ext_gap

    mejor = alineador.align(secuencia_a, secuencia_b)[0]
    alA, alB = reconstruir_con_huecos(mejor, secuencia_a, secuencia_b)
    identidad, nmatch, ngaps = calcular_identidad_y_huecos(alA, alB)

    a_ini = mejor.aligned[0][0][0] if mejor.aligned[0].size else 0
    b_ini = mejor.aligned[1][0][0] if mejor.aligned[1].size else 0

    print("\n" + "═" * 78)
    print(f"Configuración: {descripcion}")
    print(f"Modo: {modo.upper()} | Matriz: {matriz}")
    print("─" * 78)
    print(
        formatear_bloques(
            alA, alB,
            etiqueta_arriba=etiqueta_arriba,
            etiqueta_abajo=etiqueta_abajo,
            ancho=60, off_a=a_ini, off_b=b_ini
        )
    )
    print("─" * 78)
    print(f"Puntuación total         : {mejor.score:.3f}")
    print(f"Porcentaje de identidad  : {identidad}%")
    print(f"Número de coincidencias  : {nmatch}")
    print(f"Número total de huecos   : {ngaps}")
    print("═" * 78)


def comparar_matrices(secuencia_a: str, secuencia_b: str, titulo: str,
                      matrices: List[str], modo: str):
    print("\n" + "=" * 78)
    print(f"{titulo.center(78)}")
    print("=" * 78)
    configuraciones = []
    for M in matrices:
        configuraciones.append((f"{modo.capitalize()} (afín)   {M}  open=-10  ext=-0.5", modo, M, None, -10, -0.5))
        configuraciones.append((f"{modo.capitalize()} (lineal) {M}  gap=-2",             modo, M, -2,  None, None))

    for desc, md, M, gap_lineal, og, eg in configuraciones:
        ejecutar_alineamiento(
            secuencia_a, secuencia_b, desc, md,
            matriz=M,
            gap_lineal=gap_lineal,
            open_gap=og if og is not None else -10,
            ext_gap=eg if eg is not None else -0.5,
            etiqueta_arriba="Proteína A",
            etiqueta_abajo="Proteína B",
        )

def leer_secuencia_fasta(ruta_fichero: Path) -> str:
    registro = next(SeqIO.parse(str(ruta_fichero), "fasta"))
    seq = str(registro.seq)
    return limpiar_proteina(seq)


def main():
    parser = argparse.ArgumentParser(
        description="Alineamiento de proteínas con PairwiseAligner y matrices de sustitución"
    )
    parser.add_argument("--fasta", nargs=2, metavar=("SEQ1.fasta", "SEQ2.fasta"),
                        help="Alinear dos secuencias desde archivos FASTA (p. ej., descargados de UniProt)")
    parser.add_argument("--lenA", type=int, default=60, help="Longitud de la proteína aleatoria A")
    parser.add_argument("--lenB", type=int, default=48, help="Longitud de la proteína aleatoria B")
    parser.add_argument("--matrices", nargs="+", default=["BLOSUM62", "PAM250"],
                        help='Lista de matrices (incluye nombres de Biopython como "BLOSUM62", "PAM250" o rutas a fichero personalizado)')
    parser.add_argument("--modo", choices=["global", "local"], default="global",
                        help="Modo de alineamiento: global o local (por defecto: global)")
    args = parser.parse_args()

    if args.fasta:
        f1, f2 = map(Path, args.fasta)
        prot_a, prot_b = leer_secuencia_fasta(f1), leer_secuencia_fasta(f2)
        print("\nArchivos FASTA cargados correctamente:")
        print(f" - {f1.name} (longitud: {len(prot_a)})")
        print(f" - {f2.name} (longitud: {len(prot_b)})")
        comparar_matrices(prot_a, prot_b, titulo="(b) Secuencias de bases de datos (p. ej., UniProt)",
                          matrices=args.matrices, modo=args.modo)
        print("\nCómo se obtuvieron: Descarga los FASTA desde UniProt (https://www.uniprot.org/) buscando las proteínas de interés, abre la entrada y usa el botón “Download” → FASTA.")
    else:
        prot_a = generar_proteina_aleatoria(args.lenA)
        prot_b = generar_proteina_aleatoria(args.lenB)
        print("\nProteínas aleatorias generadas:")
        print(f" - Proteína A (longitud {len(prot_a)}): {prot_a}")
        print(f" - Proteína B (longitud {len(prot_b)}): {prot_b}")
        comparar_matrices(prot_a, prot_b, titulo="(a) Proteínas aleatorias",
                          matrices=args.matrices, modo=args.modo)


if __name__ == "__main__":
    main()