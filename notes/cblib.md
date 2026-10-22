# CBLIB como fuente de instancias: notas para el informe

Notas técnicas de cómo se incorporan las instancias de CBLIB al conjunto de
instancias SDP de la app. Sirven de base para la sección de fuentes de
instancias de la tesis.

## 1. Qué es CBLIB

- *Conic Benchmark Library*, mantenida en ZIB: <https://cblib.zib.de>.
- Referencia: H. A. Friberg, "CBLIB 2014: a benchmark library for conic
  mixed-integer and continuous optimization", *Mathematical Programming
  Computation* 8(3), 2016. doi:10.1007/s12532-015-0092-4.
- Formato de archivo: *Conic Benchmark Format* (CBF), versión 4. El manual
  técnico es <https://cblib.zib.de/doc/format4.pdf>.
- Contenido: problemas cónicos mixto-enteros y continuos con conos lineales,
  de segundo orden (Q, QR), exponenciales, de potencia y semidefinidos. La
  mayor parte **no** es SDP.
- Las instancias con cono PSD son, en su mayoría, **SDP mixto-enteros
  (MISDP)** de la colección de Gally, Pfetsch y Ulbrich (SCIP-SDP):
  k-particionamiento mínimo (`kpart_*`), clique (`clique_*`), mínimos
  cuadrados con restricción de cardinalidad (`clsq_*`, `coloncancer_*`),
  diseño de reticulados (*truss topology*) robusto o con barras enteras
  (`rt_*`, `bridge_*`) y diseño de experimentos (`expdesign_*`).
  *Referencia sugerida: T. Gally, M. E. Pfetsch, S. Ulbrich, "A framework for
  solving mixed-integer semidefinite programs", Optimization Methods and
  Software 33(3), 2018. Verificar los datos bibliográficos antes de citar.*
- Volumen del índice `download/all/` (octubre de 2026): 3346 archivos
  `.cbf.gz`, unos 7,2 GB comprimidos.

## 2. Pipeline en la app

Es el mismo pipeline que usan las demás fuentes (`download → extract → sort →
prepare`). CBLIB está registrada en `app/config/instances/sources.json` como
fuente de tipo `index`.

| Etapa | Qué pasa con CBLIB | Archivo |
|---|---|---|
| Descarga | Se lista el índice y se bajan todos los `.cbf.gz` | `src/renderer/renderer.ts` (`INDEX_FILE_SUFFIXES`), `scripts/instances/list_files.py`, `download.py` |
| Extracción | Los `.cbf.gz` se descomprimen a `.cbf` | `scripts/instances/extract.py` |
| Clasificación | `.cbf` se reconoce como instancia | `scripts/instances/sort.py` |
| Preparación | Filtro, relajación y conversión a SDPA (`.dat-s`) | `scripts/instances/prepare.py`, `instance_grouping.py` |

El filtrado ocurre en la etapa de preparación, no en la descarga. Por eso se
descarga la librería completa, aunque solo una fracción termina en
`data/instances/prepared/CBLIB/`.

El conversor está en `app/scripts/instances/format_converters/`:

- `cbf_parser.py`: lee el formato CBF.
- `cbf_to_sdpa.py`: hace el mapeo matemático a SDPA (sección 4).
- `cbf_to_dats.py`: aplica filtro y relajación, escribe el `.dat-s` y expone
  la interfaz de línea de comandos.

## 3. Filtro y relajación

**Filtro.** Una instancia se conserva solo si cumple las dos condiciones
siguientes:

1. Todos sus conos escalares, en `VAR` y en `CON`, son lineales: `F` (libre),
   `L+`, `L-` o `L=`.
2. Tiene al menos una variable PSD (`PSDVAR`) o una restricción LMI
   (`PSDCON`).

Se descartan las instancias con conos de segundo orden, exponenciales, de
potencia, paramétricos o `SVECPSD`. **No se reformulan.** Un cono de segundo
orden podría escribirse como bloque PSD de tipo "flecha", pero eso cambiaría
la estructura de la instancia y, con ella, las *features* que se extraen. El
criterio es que la estructura SDP que ve el análisis sea la que escribieron
los autores de la instancia.

También se descarta una instancia cuando, tras la conversión, queda una
restricción vacía con lado derecho distinto de cero, porque sería
trivialmente infactible o no acotada. Las restricciones vacías con lado
derecho cero se eliminan sin más; CSDP no acepta matrices de restricción
vacías.

Al preparar, la app informa cuántas instancias se descartaron por cada
motivo; por ejemplo: `Filtered out N CBF instance(s): non-linear cones QR`.

**Relajación.** Se ignora la sección `INT`, es decir, se elimina la
integralidad. Una instancia MISDP se transforma así en su **relajación SDP
continua**. Las cotas de las variables, por ejemplo `0 ≤ x ≤ 1` en variables
binarias, se mantienen porque vienen como restricciones lineales explícitas.

**Advertencia para el informe:** estas instancias son relajaciones de
problemas MISDP, no SDP nativos. Hay que decirlo explícitamente.

## 4. Mapeo matemático CBF → SDPA

SDPA trabaja con el par primal-dual

$$
\text{(P)}\quad \min_{x\in\mathbb{R}^m} c^\top x \;\;\text{s.a.}\;\; \sum_{k=1}^m F_k x_k - F_0 \succeq 0,
\qquad
\text{(D)}\quad \max_{Y\succeq 0} \langle F_0, Y\rangle \;\;\text{s.a.}\;\; \langle F_k, Y\rangle = c_k,
$$

donde las matrices son diagonales por bloques: hay bloques SDP y un bloque
diagonal "LP" para las restricciones lineales.

En CBF el problema es $\min$ (o $\max$) de $g_{obj}$ sujeto a:

- restricciones escalares $g_i = \sum_j a_{ij}x_j + \sum_j\langle F_{ij},X_j\rangle + b_i$, cada una en un cono lineal;
- LMI $G_i = \sum_j x_j H_{ij} + D_i \succeq 0$;
- variables escalares $x_j$ en conos lineales, y variables matriciales $X_j \succeq 0$.

Hay dos rutas de conversión. Se elige automáticamente según la forma de la
instancia.

### 4.1 Ruta LMI → primal de SDPA (`lmi_to_sdpa`)

Se usa cuando la instancia no tiene `PSDVAR` ni restricciones escalares de
igualdad (`L=`). Es el caso de casi todas las MISDP de CBLIB (kpart, clique,
clsq, rt_bridge, coloncancer).

- Las variables escalares de CBF pasan a ser las $x$ libres de SDPA.
- Las variables en `L=` son cero; se eliminan.
- Cada LMI $G_i$ es un bloque SDP, con $F_k = H_{ik}$ y $F_0 = -D_i$.
- Cada $x_j \ge 0$ (o $\le 0$) y cada restricción escalar `L+` o `L-` es una
  entrada del bloque LP. Para $g_i \ge 0$ se usa $F_k = a_{ik}$ y
  $F_0 = -b_i$; para `L-` se cambia el signo.
- Restricciones `F` (libres): se ignoran.
- Objetivo: $c = a_{obj}$ si el problema es MIN, y $c = -a_{obj}$ si es MAX.

Esta ruta es exacta y no agrega variables ni restricciones: el tamaño de la
instancia no cambia.

### 4.2 Ruta estándar → dual de SDPA (`standard_to_sdpa`)

Se usa en el resto de los casos: hay `PSDVAR`, hay igualdades escalares, o
hay mezcla de formas. Un ejemplo es `expdesign_*`. Todo se lleva a partes no
negativas de $Y$:

- Cada variable PSD $X_j$ es un bloque SDP de $Y$.
- Las variables escalares se escriben con columnas del bloque LP: si
  $x \in$ `L+`, entonces $x = p$; si $x \in$ `L-`, entonces $x = -p$; si $x$
  es libre, entonces $x = p - q$, con $p, q \ge 0$.
- Restricciones escalares, con lado derecho $c_k = -b_i$:
  - `L=` queda como igualdad;
  - `L+` se escribe $g - s = 0$ con $s\ge0$;
  - `L-` se escribe $g + s = 0$ con $s \ge 0$;
  - `F` se ignora.
- Cada LMI $G_i(x)\succeq 0$, de tamaño $n$, recibe un bloque de holgura
  $S_i \succeq 0$ con $G_i(x) - S_i = 0$. Esto agrega una igualdad por cada
  entrada del triángulo inferior, es decir, $n(n+1)/2$ restricciones.
  - Detalle: en SDPA una entrada fuera de la diagonal cuenta dos veces en
    $\langle F_k, Y\rangle$. Por eso el coeficiente de $S_{rc}$ es $-1/2$ si
    $r\ne c$ y $-1$ en la diagonal.
- Objetivo: $F_0 = -g_{obj}$ si el problema es MIN, porque SDPA maximiza en
  (D), y $F_0 = g_{obj}$ si es MAX.

### 4.3 Consideraciones

- **Constante del objetivo:** el término `OBJBCOORD` se descarta porque SDPA
  no admite constantes. El óptimo de la instancia convertida difiere del de
  CBF en esa constante.
- **Signo del óptimo:** según la ruta y el sentido del problema, el valor
  óptimo de SDPA es $+$ o $-$ el valor de CBF:
  - ruta LMI, problema MIN: mismo signo;
  - ruta LMI, problema MAX: signo opuesto;
  - ruta estándar, problema MIN: signo opuesto;
  - ruta estándar, problema MAX: mismo signo.
- **Relajaciones no acotadas:** si una variable entera no tiene cota
  explícita, la relajación podría quedar no acotada. En las muestras
  revisadas no ocurrió.
- **Hotstart:** la palabra clave `CHANGE`, que define secuencias de
  instancias, se trata como fin de archivo. Esto lo permite la
  especificación.
- **Trazabilidad:** el catálogo (`config/instances/catalog.json`) no marca
  qué instancias son relajadas; esa información queda en el log de
  preparación (`instances`). Si el informe necesita distinguirlas, por
  ahora todas las de `source: "CBLIB"` con integrales son relajaciones.

## 5. Validación realizada (2026-10-22)

**Problemas construidos a mano con óptimo conocido**, resueltos con CSDP 6.2:

| Prueba | Ruta | Óptimo esperado | Resultado SDPA |
|---|---|---|---|
| $\min x$ s.a. $\begin{psmallmatrix}x&1\\1&1\end{psmallmatrix}\succeq0$ | LMI | 1 | 1.0000 |
| `PSDVAR` + igualdades + libre + `L+` (MIN) | estándar | 6 | −6 (signo opuesto, correcto) |
| MAX con `L+`, `L-`, LMI y una variable entera | LMI | 3 | −3 (signo opuesto, correcto) |

**Consistencia entre rutas sobre instancias reales:** cada instancia se
convirtió por las dos rutas y se resolvió con CSDP. Los valores coinciden en
magnitud con signo opuesto, como se espera:

| Instancia | Ruta LMI | Ruta estándar |
|---|---|---|
| `kpart_2g_4_164_k3_5_6` | −7.0714784e5 | 7.0714784e5 |
| `clique_20_k3_6_7` | 143.93877 | −143.93877 |
| `rt_bridge_2x10_2bars_2scen` | 6.2394798 | −6.2394799 |

**Muestra de 9 instancias:**

- 6 convertidas: `clique`, `clsq_random`, `coloncancer`, `expdesign_A`,
  `kpart` y `rt_bridge`. Todas, salvo `expdesign_A`, fueron por la ruta LMI.
- 3 filtradas por conos QR: `achtziger_stolpe`, `ck_n` y `wiener`.
- `expdesign_A_12_6` resolvió con CSDP (−7.1765628) después de eliminar sus
  restricciones vacías.

## 6. Pendientes y posibles mejoras

- Filtrar ya en la descarga, para no bajar 7,2 GB: la estructura de conos
  (`VAR`, `CON`, `PSDVAR`, `PSDCON`) viene antes que los datos, así que
  bastaría leer en streaming el comienzo de cada `.cbf.gz`.
- Agregar al catálogo un campo `relaxed_from_misdp`.
- Revisar duplicados con Mittelmann y SDPLIB.
- Después de agregar las instancias, verificar en la proyección 2D del
  espacio de instancias que CBLIB efectivamente aporte cobertura nueva.
