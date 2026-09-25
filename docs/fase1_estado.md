# Fase 1 — Estado: **parcialmente bloqueada**

Fecha: 2026-09-25

## Resumen en una línea

No he podido extraer los metadatos de NCBI: **el proxy de red de este entorno
deniega `ncbi.nlm.nih.gov` y `ebi.ac.uk`**. Sí he podido verificar cosas
sustanciales sobre PRJNA1248021 a través del repositorio de araclab, que sí es
accesible. No hay ni una cifra inventada en este documento.

---

## 1. Qué está bloqueado, exactamente

Probé las dos rutas de salida del contenedor (curl a través del proxy, y la
herramienta de fetch). Ambas fallan en CONNECT con 403 para todo host de NCBI y
EBI:

| Host | Resultado |
|---|---|
| `eutils.ncbi.nlm.nih.gov` | bloqueado (403 en CONNECT) |
| `www.ncbi.nlm.nih.gov` | bloqueado |
| `trace.ncbi.nlm.nih.gov` | bloqueado |
| `ftp.ncbi.nlm.nih.gov` | bloqueado |
| `sra-download.ncbi.nlm.nih.gov` | bloqueado |
| `ftp.sra.ebi.ac.uk` | bloqueado |
| `www.ebi.ac.uk` | bloqueado |
| `github.com` / `raw.githubusercontent.com` | accesible |

**Cómo desbloquearlo:** en los ajustes del entorno (menú del entorno en la barra
de título de la sesión → Edit → Network access), subir el nivel de acceso o
añadir a los dominios permitidos:

- `ncbi.nlm.nih.gov` — imprescindible para Fase 1 (E-utilities, BioSample) y
  para Fase 2 (descarga de lecturas).
- `ebi.ac.uk` — opcional pero recomendable: ENA sirve los mismos FASTQ por
  HTTPS directo, bastante más rápido que `sra-tools` para la Fase 2.

Niveles de acceso descritos en
<https://code.claude.com/docs/en/claude-code-on-the-web>.

En cuanto esté permitido, la Fase 1 completa sale con dos comandos (ver §5).

---

## 2. Qué SÍ he podido verificar (PRJNA1248021)

Fuente: clon de `github.com/araclab/general` en HEAD `95610e8`
(`microbiome/mb_analysis/`). Esto es código y ficheros reales del laboratorio,
no recuerdo mío.

### 2.1 El pipeline

DADA2 v1.16 + QIIME 1.9.1 + RDP Classifier 2.12, orquestado en SLURM por
`run_codes/dada2_pipeline_vP_v11_cleaned.sh`. Cadena: Trimmomatic → Cutadapt
(retirada de cebadores) → `dada2_filter.R` (maxEE 2,2) → `dada2_inference.R`
(minOverlap 12) → quimeras + taxonomía (minBoot 80) → BIOM → RDP a género y a
especie con un clasificador propio ("V6-species").

Relevante para la Fase 2: la taxonomía a especie no es genérica, es un
clasificador entrenado a medida sobre una lista cerrada de 707 entradas
(`database_files/16S_Sequencing_QC_16SAcceptedSpecies.txt`). **Comprobado que
esa lista cubre los taxones que te interesan**: `Anaerococcus_vaginalis`,
`A. tetradius`, `A. lactolyticus`, `A. prevotii`, `A. octavius` y ~15
*Anaerococcus* más, además de `Atopobium_vaginae`. Los ASV fuera de ese
ámbito se caracterizan por BLAST, es decir, manualmente.

### 2.2 Región del 16S — y por qué importa más de lo que parece

El repositorio asigna cebadores por tipo de muestra, y lo dice explícitamente
en su README:

- **Muestras genitales** → `amplicon_primers_bactquant.fasta` (cita PMID 22510143)
- Muestras nasales → `amplicon_primers_fadrosh.fasta` (cita PMID 24558975)

Secuencias reales de los dos ficheros:

| Set | Forward | Reverse |
|---|---|---|
| BactQuant (genital) | `CCTACGGGDGGCWGCA` (16 nt) | `GGACTACHVGGGTMTCTAATC` (21 nt) |
| Fadrosh (nasal) | `ACTCCTACGGGAGGCAGCAG` (20 nt) | `GGACTACHVGGGTWTCTAAT` (20 nt) |

El propio pipeline declara la región en el parámetro de Trimmomatic:
`--minlen 225  # Drop the read if it is below 225 bp (for V3-V4 regions)`.

Identificación mía, corroborando: el forward cae en la posición 338/341 y el
reverse es 806R, de modo que el amplicón abarca **V3–V4**. Consistente con lo
que declara el repo.

**El punto no obvio:** los dos sets son V3–V4 y aun así **no son
intercambiables**. Difieren en longitud y en bases degeneradas (`DGGCWGCA` vs
`AGGCAGCAG`; `TMTCTAATC` vs `TWTCTAAT`), luego producen amplicones con
extremos distintos. Un ASV es una secuencia exacta: dos amplicones que
empiezan o acaban en bases diferentes generan ASV **distintos para el mismo
organismo**. Es decir, "misma región" no basta para fusionar tablas de ASV —
ni siquiera dentro de un mismo laboratorio. Esto endurece tu criterio inicial:
región distinta descarta la combinación, pero región igual **no** la autoriza.

### 2.3 Variables del estudio SECS (Park)

El repo incluye el script de análisis del propio paper
(`paper_scripts/Post-Coital_Dynamics_.../Toronto_SECS_Git.sas`, 4.010 líneas,
autor Dan Park). Las variables que maneja, leídas del código:

| Eje | Variables |
|---|---|
| Identificadores | `ID`, `PID` (identificador de pareja), `Gender` |
| Diseño | `site` (valores `Cervical`, `CorSul` = surco coronal), `time` (1, 2, 3, 4), `DayPostBln` |
| **Carga bacteriana** | `std_qty`, `log_std_qty` |
| **Estado de VB** | `BVStatusA`, `BVStatusE`, `NugentSA`, `NugentSE`, y `NugentBV1…4` con sus `NugentBV1_Cat…4_Cat` |
| **Ciclo menstrual** | `Cyclephase` |
| **Circuncisión** | `circumcision` |
| Conducta sexual | `condom`, `sex_type`, `SexParN` |
| Otras | `Age`, `FemaleEthnicity`, `MaleEthnicity`, `pH_high`, datos inmunes (`immune_raw`) |

**Esto contradice tres de tus cuatro suposiciones de partida**, al menos para
este estudio:

- Carga bacteriana: **existe**. `std_qty` es cuantificación bacteriana total —
  coherente con que el set de cebadores sea el de BactQuant, que es
  precisamente un ensayo de qPCR de carga bacteriana total (PMID 22510143).
  Los cebadores del 16S y los de la qPCR son los mismos.
- Nugent por punto temporal: **existe**, y en cuatro puntos (`NugentBV1` a
  `NugentBV4`), cada uno con versión continua y categorizada.
- Fase del ciclo: **existe** (`Cyclephase`).
- Circuncisión: **existe** (`circumcision`), y además se usa como variable de
  agrupación en los gráficos de carga bacteriana.

### 2.4 Advertencia importante sobre §2.3

Esas variables están en el **fichero analítico interno de los autores**, que el
script importa con `PROC IMPORT DATAFILE="path/file"` — una ruta local que no
está en el repositorio. **No he verificado que ninguna de ellas esté depositada
en el BioSample público.** Lo habitual en SRA es que el BioSample lleve cuatro
atributos mínimos y que el resto viva en las tablas suplementarias del paper o
sólo en el ordenador del grupo.

Por eso el script de Fase 2 (`02_parse_metadata.py`) hace esa distinción
explícita en su tabla de auditoría: separa "no está en el registro público" de
"no existe". Son cosas distintas y la diferencia decide si el reanálisis es
viable o si hay que escribir a los autores.

---

## 3. Qué sigue siendo desconocido

Todo lo que pediste para **PRJNA398590 y PRJNA735440**, sin excepción: número
de muestras, de sujetos, de parejas emparejadas, puntos temporales, sitios,
región del 16S, plataforma y variables clínicas. No tengo forma de verificarlo
sin NCBI y no voy a rellenarlo de memoria.

De PRJNA1248021 tampoco sé aún: número de runs, de muestras, de sujetos, ni qué
subconjunto de §2.3 llegó al depósito público.

**En consecuencia, no puedo darte todavía el veredicto de combinabilidad.**
Lo único que puedo adelantar es el criterio, ya afinado en §2.2: aunque los
tres coincidieran en V3–V4, la fusión a nivel de ASV seguiría sin estar
justificada salvo que compartan cebadores exactos; lo defendible sería un
meta-análisis a nivel de género con el proyecto como efecto aleatorio.

---

## 4. Sobre la pregunta (a), la principal

Un apunte de viabilidad que sí puedo hacer ahora: tu pregunta sobre
*Corynebacterium*, *Staphylococcus*, *Anaerococcus* y *Finegoldia* tras
clindamicina es, en parte, una pregunta sobre **abundancia absoluta, no
relativa**. Si los anaerobios asociados a VB se desploman con el tratamiento y
los comensales cutáneos no, la proporción de estos últimos sube aunque su
carga real no se mueva. Con datos composicionales de 16S a secas, "reaparición"
y "no haber bajado nunca" son indistinguibles.

Eso hace que `std_qty` de PRJNA1248021 sea más valioso de lo que parecía, y
convierte en pregunta crítica de Fase 1 si los Plummer depositaron alguna
medida de carga. Si no la hay, la pregunta (a) es respondible sólo de forma
cualitativa, y conviene saberlo antes de descargar nada.

---

## 5. Qué queda listo para ejecutar

```bash
bash scripts/01_fetch_sra_metadata.sh data/raw_metadata
python3 scripts/02_parse_metadata.py data/raw_metadata docs/fase1_tablas.md
```

- `01_fetch_sra_metadata.sh` — runinfo, BioSample XML, experiment XML y
  BioProject XML de los tres proyectos. Reintentos con backoff, respeta el
  límite de peticiones de NCBI, admite `NCBI_API_KEY`. **Hace preflight**: si
  la red sigue bloqueada, aborta con el mensaje de §1 en vez de dejar ficheros
  vacíos que luego parezcan datos.
- `02_parse_metadata.py` — sólo stdlib (este contenedor no tiene pandas ni R).
  Genera la tabla comparativa, el inventario completo de atributos con tasa de
  relleno, la estructura del diseño (sujeto / punto temporal / sitio / pareja),
  la evidencia textual de la región 16S y la auditoría de ausencias. Cada
  detección va con el nombre del atributo y sus valores, para que sea
  verificable y no un veredicto a ciegas.

Probado contra un fixture sintético: detecta correctamente sujeto, punto
temporal, sitio, emparejamiento de pareja y región, y marca como AUSENTE lo que
no está. El fixture es de prueba y no está en `data/`.

---

## 6. Decisión que te toca a ti

1. Desbloquear `ncbi.nlm.nih.gov` en los ajustes de red del entorno y me dices,
   y termino la Fase 1 entera en una sola pasada; o
2. Si prefieres no tocar la red, ejecutas tú los dos comandos de §5 donde
   tengas salida a internet y me pasas `docs/fase1_tablas.md`, y sigo desde
   ahí.

No he pasado a la Fase 2, que es lo correcto: descargar y procesar antes de
saber si los metadatos sostienen las preguntas sería gastar cómputo a ciegas.
