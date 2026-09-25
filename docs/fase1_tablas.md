# Fase 1 — Inventario de metadatos (generado automáticamente)

## 1. Tabla comparativa

|  | PRJNA398590 | PRJNA735440 | PRJNA1248021 |
|---|---|---|---|
| Estudio | Plummer 2018 — parejas, tratamiento antibiótico, 4 semanas | Plummer 2021 — 34 parejas, 12 semanas, 3 sitios | Park 2026 — dinámica post-coital (SECS, Toronto) |
| Fuente de los datos | `runinfo.csv` | `runinfo.csv` | `runinfo.csv` |
| Runs (SRR) | 140 | 418 | 284 |
| BioSamples | — | — | — |
| Sujetos (inferido) | no inferible | no inferible | no inferible |
| Atributo sujeto | — | — | — |
| Región 16S | **341F, V3-V4** | **341F, V3-V4** | **V3-V4** |
| Plataforma | ILLUMINA (140) | ILLUMINA (418) | ILLUMINA (284) |
| Modelo | Illumina MiSeq (140) | Illumina MiSeq (418) | Illumina MiSeq (284) |
| Layout | PAIRED (140) | PAIRED (418) | PAIRED (284) |
| Estrategia | AMPLICON (140) | AMPLICON (418) | AMPLICON (284) |
| Profundidad mediana/run | 35,373 spots | 35,712 spots | 23,950 spots |

## 2. Auditoría de ausencias

Marcado **AUSENTE** = ningún atributo de BioSample con ese nombre o valor. Un atributo puede existir en el dataset analítico de los autores y aun así no estar depositado en BioSample: esta tabla sólo habla del registro público.

| Variable | PRJNA398590 | PRJNA735440 | PRJNA1248021 |
|---|---|---|---|
| Carga bacteriana (qPCR / BactQuant) | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Estado de VB (Nugent / Amsel) | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Fase del ciclo menstrual | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Estado de circuncisión | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Tratamiento / brazo del ensayo | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Uso de preservativo / actividad sexual | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| pH vaginal | **AUSENTE** | **AUSENTE** | **AUSENTE** |

## 3. Estructura del diseño detectada

### PRJNA398590 — Plummer 2018 — parejas, tratamiento antibiótico, 4 semanas

- **sujeto**: no detectado
- **punto_temporal**: no detectado
- **sitio_anatomico**: no detectado
- **pareja**: no detectado
- **sexo**: no detectado

Evidencia de región 16S:

  - DESIGN_DESCRIPTION: Dual indexed universal primers 341F (CCTACGGGNGGCWGCAG) and 805R (GACTACHVGGGTATCTAATCC) were used for PCR amplification of the V3-V4 hypervariable regions of the 16S rRNA gene, as previously described in Shipitsyna et al 2013;8(4) PloseOne and Fadrosh et al 2014;2(1):6 Microbiome.

### PRJNA735440 — Plummer 2021 — 34 parejas, 12 semanas, 3 sitios

- **sujeto**: no detectado
- **punto_temporal**: no detectado
- **sitio_anatomico**: no detectado
- **pareja**: no detectado
- **sexo**: no detectado

Evidencia de región 16S:

  - DESIGN_DESCRIPTION: Universal primers 341F (CCTACGGGNGGCWGCAG) and 805R (GACTACHVGGGTATCTAATCC) were used for PCR amplification of the V3-V4 hypervariable regions of the 16S rRNA gene

### PRJNA1248021 — Park 2026 — dinámica post-coital (SECS, Toronto)

- **sujeto**: no detectado
- **punto_temporal**: no detectado
- **sitio_anatomico**: no detectado
- **pareja**: no detectado
- **sexo**: no detectado

Evidencia de región 16S:

  - DESIGN_DESCRIPTION: 16SrRNA amplicon sequencing of the V3-V4 region using bactquant primers with heterogeneity spacers

## 4. Inventario completo de atributos BioSample

### PRJNA398590

_Sin atributos._

### PRJNA735440

_Sin atributos._

### PRJNA1248021

_Sin atributos._

## 5. Combinabilidad

- PRJNA398590: 341F, V3-V4
- PRJNA735440: 341F, V3-V4
- PRJNA1248021: V3-V4

**NO combinables**: las regiones del 16S difieren. Los ASV no son comparables entre proyectos y los sesgos de amplificación por taxón tampoco. Analizar por separado.
