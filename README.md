# Reanálisis 16S — microbioma genital en parejas

Evaluación de viabilidad y generación de hipótesis sobre tres BioProjects del
NCBI SRA:

| BioProject | Estudio |
|---|---|
| PRJNA398590 | Plummer 2018 — parejas, tratamiento antibiótico, 4 semanas |
| PRJNA735440 | Plummer 2021 — 34 parejas, 12 semanas, 3 sitios anatómicos |
| PRJNA1248021 | Park 2026 — dinámica post-coital ([pipeline araclab](https://github.com/araclab/general/tree/main/microbiome/mb_analysis)) |

## Estado

**Fase 1: parcialmente bloqueada.** El proxy de red del entorno deniega
`ncbi.nlm.nih.gov` y `ebi.ac.uk`, así que no se han podido extraer los
metadatos. Ver [`docs/fase1_estado.md`](docs/fase1_estado.md) para el detalle,
lo que sí se ha verificado por otras vías y cómo desbloquearlo.

Fases 2 y 3: no iniciadas, condicionadas al resultado de la Fase 1.

## Uso

```bash
bash scripts/01_fetch_sra_metadata.sh data/raw_metadata
python3 scripts/02_parse_metadata.py data/raw_metadata docs/fase1_tablas.md
```

Opcional: `export NCBI_API_KEY=...` para subir el límite de NCBI de 3 a 10
peticiones/s.

El primer script hace preflight de red y aborta con instrucciones si NCBI
sigue bloqueado, en vez de dejar ficheros vacíos. El segundo sólo usa la
stdlib de Python 3.

## Principio

Ningún número de este repositorio procede de recuerdo o inferencia: o sale de
un fichero descargado, o se marca explícitamente como desconocido.
