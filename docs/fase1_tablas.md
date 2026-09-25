# Fase 1 — Inventario de metadatos (generado automáticamente)

## 1. Tabla comparativa

|  | PRJNA398590 | PRJNA735440 | PRJNA1248021 |
|---|---|---|---|
| Estudio | Plummer 2018 — parejas, tratamiento antibiótico, 4 semanas | Plummer 2021 — 34 parejas, 12 semanas, 3 sitios | Park 2026 — dinámica post-coital (SECS, Toronto) |
| Fuente de los datos | `runinfo.csv`, `biosample.xml`, `experiment.xml` | `runinfo.csv`, `biosample.xml`, `experiment.xml` | `runinfo.csv`, `biosample.xml`, `experiment.xml` |
| Runs (SRR) | 140 | 418 | 284 |
| BioSamples | 140 | 418 | 284 |
| Sujetos (inferido) | no inferible | 54 | 284 |
| Atributo sujeto | — | `host_subject_id` | `host_subject_id` |
| Región 16S | **V3-V4** | **V3-V4** | **V3-V4** |
| Plataforma | ILLUMINA (140) | ILLUMINA (418) | ILLUMINA (284) |
| Modelo | Illumina MiSeq (140) | Illumina MiSeq (418) | Illumina MiSeq (284) |
| Layout | PAIRED (140) | PAIRED (418) | PAIRED (284) |
| Estrategia | AMPLICON (140) | AMPLICON (418) | AMPLICON (284) |
| Profundidad mediana/run | 35,373 spots | 35,712 spots | 23,950 spots |

## 2. Auditoría de ausencias

**AUSENTE** = hay registros BioSample y ninguno declara ese atributo. **SIN DATOS** = no se han podido leer registros BioSample, así que no hay nada que auditar; no es un hallazgo, es una laguna de la descarga.

Un atributo puede existir en el dataset analítico de los autores y aun así no estar depositado en BioSample: esta tabla sólo habla del registro público.

| Variable | PRJNA398590 | PRJNA735440 | PRJNA1248021 |
|---|---|---|---|
| Carga bacteriana (qPCR / BactQuant) | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Estado de VB (Nugent / Amsel) | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Fase del ciclo menstrual | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Estado de circuncisión | **AUSENTE** | **AUSENTE** | `host_phenotype` (100%) |
| Tratamiento / brazo del ensayo | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| Uso de preservativo / actividad sexual | **AUSENTE** | **AUSENTE** | **AUSENTE** |
| pH vaginal | **AUSENTE** | **AUSENTE** | **AUSENTE** |

## 3. Estructura del diseño detectada

### PRJNA398590 — Plummer 2018 — parejas, tratamiento antibiótico, 4 semanas

- **sujeto**: no detectado
- **punto_temporal**: `collection_date` (140 muestras, 100%) → 2015-08-18, 2015-08-19, 2015-08-21, 2015-08-24, 2015-08-25 …; `day` (140 muestras, 100%) → D0, D14, D21, D28, D8
- **sitio_anatomico**: no detectado
- **pareja**: no detectado
- **sexo**: no detectado

Nombres de muestra (140 runs, separador `.`, reparto de campos {2: 62, 3: 78}):

  - ejemplos: `3002.D8`, `3003.D0`, `3003.D28`, `3003.D8`, `3009.D0`, `3009.D14`, `3009.D28`, `3009.D8`
  - campo 1: 21 valores → 4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008, 4009, 4010, 4012, 4013 …
  - campo 2: 6 valores → 28, D0, D14, D21, D28, D8
  - campo 3: 2 valores → S, U

Evidencia de región 16S:

  - DESIGN_DESCRIPTION: Dual indexed universal primers 341F (CCTACGGGNGGCWGCAG) and 805R (GACTACHVGGGTATCTAATCC) were used for PCR amplification of the V3-V4 hypervariable regions of the 16S rRNA gene, as previously described in Shipitsyna et al 2013;8(4) PloseOne and Fadrosh et al 2014;2(1):6 Microbiome.

### PRJNA735440 — Plummer 2021 — 34 parejas, 12 semanas, 3 sitios

- **sujeto**: `host_subject_id` (375 muestras, 90%) → 12F, 12M, 15F, 15M, 18F …
- **punto_temporal**: `collection_date` (404 muestras, 97%) → 2018-03-28, 2018-03-29, 2018-04-06, 2018-04-10, 2018-04-14 …
- **sitio_anatomico**: `isolation_source` (418 muestras, 100%) → Community standard - ZymoBIOMICS Microbial Community Standard Catalog No. D6300, Community standard -BEI resources HM-276D Genomic DNA from Microbial Mock Community B, First pass urine, Genelock AssayAssure, PBS …
- **pareja**: no detectado
- **sexo**: `host_sex` (375 muestras, 90%) → female, male
- **muestras por sujeto**: 4 muestras×5 sujetos, 5 muestras×21 sujetos, 6 muestras×5 sujetos, 8 muestras×5 sujetos, 10 muestras×18 sujetos

Nombres de muestra (418 runs, separador `_`, reparto de campos {2: 375, 1: 43}):

  - ejemplos: `29V_End`, `29V_W4`, `29V_W8`, `30V_D0`, `30V_D8`, `30V_End`, `30V_W4`, `30V_W8`
  - campo 1: 81 valores → 12P, 12U, 12V, 15P, 15U, 15V, 18P, 18U, 18V, 20P, 20U, 20V …
  - campo 2: 6 valores → D0, D8, End, R, W4, W8

Evidencia de región 16S:

  - DESIGN_DESCRIPTION: Universal primers 341F (CCTACGGGNGGCWGCAG) and 805R (GACTACHVGGGTATCTAATCC) were used for PCR amplification of the V3-V4 hypervariable regions of the 16S rRNA gene

### PRJNA1248021 — Park 2026 — dinámica post-coital (SECS, Toronto)

- **sujeto**: `host_subject_id` (284 muestras, 100%) → #142_A, #142_E, #142_H, #142_P, #154_A …
- **punto_temporal**: no detectado
- **sitio_anatomico**: no detectado
- **pareja**: no detectado
- **sexo**: `host_sex` (284 muestras, 100%) → female, male
- **muestras por sujeto**: 1 muestras×284 sujetos

Nombres de muestra (284 runs, separador `-`, reparto de campos {3: 160, 2: 110, 5: 13, 4: 1}):

  - ejemplos: `SECS-179-A`, `SECS-179-E`, `SECS-183-H`, `SECS-415-H`, `SECS-415-P`, `SECS-436-A`, `SECS-436-E`, `SECS-436-H`
  - campo 1: 1 valores → SECS
  - campo 2: 49 valores → 142, 154, 179, 183, 190, 213, 214, 224, 415, 416, 418, 419 …
  - campo 3: 8 valores → A, A2, E, E2, H, H2, P, P2

Evidencia de región 16S:

  - DESIGN_DESCRIPTION: 16SrRNA amplicon sequencing of the V3-V4 region using bactquant primers with heterogeneity spacers

## 4. Inventario completo de atributos BioSample

### PRJNA398590

| Atributo | n | % relleno | n distintos | ejemplos |
|---|---|---|---|---|
| collection_date | 140 | 100% | 74 | 2015-08-18, 2015-08-19, 2015-08-21, 2015-08-24 … |
| day | 140 | 100% | 5 | D0, D14, D21, D28 … |
| env_broad_scale | 140 | 100% | 3 | human penile skin, human urine, human vagina |
| env_local_scale | 140 | 100% | 3 | penile skin, urine, vagina |
| env_medium | 140 | 100% | 3 | penile skin swab, urine swab, vaginal swab |
| geo_loc_name | 140 | 100% | 1 | Australia: Melbourne |
| host | 140 | 100% | 1 | Homo sapiens |
| id | 140 | 100% | 42 | 3001, 3002, 3003, 3004 … |
| lat_lon | 140 | 100% | 1 | 37.81 S 144.96 E |

### PRJNA735440

| Atributo | n | % relleno | n distintos | ejemplos |
|---|---|---|---|---|
| isolation_source | 418 | 100% | 11 | Community standard - ZymoBIOMICS Microbial Community Standard Catalog No. D6300, Community standard -BEI resources HM-276D Genomic DNA from Microbial Mock Community B, First pass urine, Genelock AssayAssure … |
| collection_date | 404 | 97% | 159 | 2018-03-28, 2018-03-29, 2018-04-06, 2018-04-10 … |
| geo_loc_name | 375 | 90% | 1 | Australia: Melbourne |
| host | 375 | 90% | 1 | Homo sapiens |
| host_sex | 375 | 90% | 2 | female, male |
| host_subject_id | 375 | 90% | 54 | 12F, 12M, 15F, 15M … |
| samp_collect_device | 375 | 90% | 2 | Copan flocked swab, Urine collected into 15ml tube with 830µl of AssayAssure® Genelock |
| Aliquot | 43 | 10% | 14 | 1, 10, 11, 12 … |
| env_broad_scale | 0 | 0% | 0 |  |
| env_local_scale | 0 | 0% | 0 |  |
| env_medium | 0 | 0% | 0 |  |
| lat_lon | 0 | 0% | 0 |  |

### PRJNA1248021

| Atributo | n | % relleno | n distintos | ejemplos |
|---|---|---|---|---|
| env_broad_scale | 284 | 100% | 1 | Genital |
| env_medium | 284 | 100% | 1 | Swab |
| ethnicity | 284 | 100% | 4 | ACB, Asian, Mixed, White |
| host | 284 | 100% | 1 | Homo sapiens |
| host_age | 284 | 100% | 18 | 18, 19, 20, 21 … |
| host_phenotype | 284 | 100% | 2 | circumcised, uncircumcised |
| host_sex | 284 | 100% | 2 | female, male |
| host_subject_id | 284 | 100% | 284 | #142_A, #142_E, #142_H, #142_P … |
| env_local_scale | 283 | 100% | 2 | Cervical, CorSul |
| collection_date | 0 | 0% | 0 |  |
| geo_loc_name | 0 | 0% | 0 |  |
| isol_growth_condt | 0 | 0% | 0 |  |
| lat_lon | 0 | 0% | 0 |  |
| strain | 0 | 0% | 0 |  |

## 5. Combinabilidad

| Proyecto | Región | Cebadores (nombres) | Cebadores (secuencia) |
|---|---|---|---|
| PRJNA398590 | V3-V4 | 341F, 805R | CCTACGGGNGGCWGCAG<br>GACTACHVGGGTATCTAATCC |
| PRJNA735440 | V3-V4 | 341F, 805R | CCTACGGGNGGCWGCAG<br>GACTACHVGGGTATCTAATCC |
| PRJNA1248021 | V3-V4 | — | **no declarada** |

Veredicto por pares (la combinabilidad no es global: dos proyectos pueden ser fusionables entre sí y ninguno de los dos con el tercero).

| Par | Veredicto |
|---|---|
| PRJNA398590 ↔ PRJNA735440 | **SÍ a nivel de ASV** — misma región y cebadores idénticos. Sigue haciendo falta procesar junto (un solo `learnErrors` por run de secuenciación) y modelar el estudio como efecto de lote |
| PRJNA398590 ↔ PRJNA1248021 | Misma región, cebadores no declarados en al menos uno. Sólo a nivel de género, y con cautela |
| PRJNA735440 ↔ PRJNA1248021 | Misma región, cebadores no declarados en al menos uno. Sólo a nivel de género, y con cautela |

Nota sobre el criterio: la región variable por sí sola no decide. Lo que fija si dos proyectos producen el mismo ASV para el mismo organismo son las secuencias exactas de los cebadores, porque son las que determinan dónde empieza y acaba el amplicón una vez recortados. Misma región con cebadores distintos obliga a agregar a género antes de comparar.

(Los tres proyectos declaran V3-V4.)
