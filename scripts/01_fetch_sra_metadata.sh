#!/usr/bin/env bash
# Fase 1 — descarga de metadatos (NO descarga secuencias).
#
# Para cada BioProject obtiene:
#   - runinfo.csv   : una fila por Run (SRR), vía esearch+efetch sobre db=sra
#   - biosample.xml : registro BioSample completo, con todos los <Attribute>
#   - bioproject.xml: descripción del proyecto
#
# Requiere salida HTTPS hacia eutils.ncbi.nlm.nih.gov.
# Si el proxy de red del entorno bloquea ese host, el preflight aborta con
# un mensaje explícito en vez de generar ficheros vacíos.
#
# Uso:  bash scripts/01_fetch_sra_metadata.sh [directorio_salida]
# Opcional: export NCBI_API_KEY=... para subir el límite a 10 req/s.

set -euo pipefail

OUTDIR="${1:-data/raw_metadata}"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PROJECTS=(PRJNA398590 PRJNA735440 PRJNA1248021)

# Sin API key NCBI limita a 3 peticiones/s; con key, a 10.
if [[ -n "${NCBI_API_KEY:-}" ]]; then
  KEYARG="&api_key=${NCBI_API_KEY}"
  THROTTLE=0.11
else
  KEYARG=""
  THROTTLE=0.35
fi

mkdir -p "$OUTDIR"

# curl con reintentos y backoff; falla si el cuerpo viene vacío.
fetch() {
  local url="$1" dest="$2" label="$3"
  local attempt
  for attempt in 1 2 3 4; do
    if curl -sS --fail --max-time 120 "$url" -o "$dest" 2>"$dest.err"; then
      if [[ -s "$dest" ]]; then
        rm -f "$dest.err"
        sleep "$THROTTLE"
        return 0
      fi
      echo "  ! $label: respuesta vacía (intento $attempt)" >&2
    else
      echo "  ! $label: curl falló (intento $attempt): $(tr -d '\n' < "$dest.err")" >&2
    fi
    sleep $((2 ** attempt))
  done
  echo "  X $label: agotados los reintentos" >&2
  return 1
}

preflight() {
  echo "== Preflight: comprobando acceso a NCBI E-utilities =="
  local code
  # curl escribe "000" y además sale con error cuando el CONNECT es rechazado,
  # así que no se puede encadenar un `|| echo 000`: duplicaría el código.
  code=$(curl -sS -o /dev/null -w '%{http_code}' --max-time 30 \
         "${EUTILS}/einfo.fcgi?db=sra&retmode=json${KEYARG}" 2>/dev/null) || true
  code="${code:-000}"
  if [[ "$code" != "200" ]]; then
    cat >&2 <<EOF

ABORTADO: no hay salida de red hacia eutils.ncbi.nlm.nih.gov (HTTP ${code}).

El proxy de egreso de este entorno deniega el host. Para desbloquearlo hay que
permitir, en la configuración de red del entorno (Network access → allowed
domains):

    ncbi.nlm.nih.gov        (E-utilities, BioSample, BioProject, descarga SRA)
    ebi.ac.uk               (opcional: ENA, espejo de los FASTQ, suele ser más
                             rápido que sra-tools para la Fase 2)

Documentación: https://code.claude.com/docs/en/claude-code-on-the-web

No se ha escrito ningún fichero. Este script no inventa metadatos.
EOF
    exit 1
  fi
  echo "   OK (HTTP 200)"
  echo
}

# Devuelve "WebEnv QueryKey Count" para un término en una db dada.
history_handle() {
  local db="$1" term="$2" tmp
  tmp=$(mktemp)
  fetch "${EUTILS}/esearch.fcgi?db=${db}&term=${term}&usehistory=y&retmode=json${KEYARG}" \
        "$tmp" "esearch ${db} ${term}"
  python3 - "$tmp" <<'PY'
import json, sys
with open(sys.argv[1]) as fh:
    r = json.load(fh)["esearchresult"]
print(r.get("webenv", ""), r.get("querykey", ""), r.get("count", "0"))
PY
  rm -f "$tmp"
}

preflight

for prj in "${PROJECTS[@]}"; do
  echo "== ${prj} =="
  pdir="${OUTDIR}/${prj}"
  mkdir -p "$pdir"

  # --- SRA runinfo ---
  read -r webenv qkey count <<<"$(history_handle sra "${prj}%5BBioProject%5D")"
  echo "   SRA: ${count} runs"
  if [[ "$count" == "0" ]]; then
    echo "   ! Sin resultados en db=sra para ${prj}. ¿Proyecto embargado o con otro accession?" >&2
  else
    fetch "${EUTILS}/efetch.fcgi?db=sra&WebEnv=${webenv}&query_key=${qkey}&rettype=runinfo&retmode=text${KEYARG}" \
          "${pdir}/runinfo.csv" "runinfo ${prj}"
    # El XML completo del experimento lleva los campos de librería/primers que
    # runinfo recorta (LIBRARY_CONSTRUCTION_PROTOCOL, target_subfragment...).
    fetch "${EUTILS}/efetch.fcgi?db=sra&WebEnv=${webenv}&query_key=${qkey}&retmode=xml${KEYARG}" \
          "${pdir}/experiment.xml" "experiment xml ${prj}" || true
  fi

  # --- BioSample ---
  read -r bwebenv bqkey bcount <<<"$(history_handle biosample "${prj}%5BBioProject%5D")"
  echo "   BioSample: ${bcount} muestras"
  if [[ "$bcount" != "0" ]]; then
    fetch "${EUTILS}/efetch.fcgi?db=biosample&WebEnv=${bwebenv}&query_key=${bqkey}&rettype=full&retmode=xml${KEYARG}" \
          "${pdir}/biosample.xml" "biosample ${prj}"
  fi

  # --- BioProject ---
  read -r pwebenv pqkey pcount <<<"$(history_handle bioproject "${prj}%5BProject%20Accession%5D")"
  if [[ "$pcount" != "0" ]]; then
    fetch "${EUTILS}/efetch.fcgi?db=bioproject&WebEnv=${pwebenv}&query_key=${pqkey}&retmode=xml${KEYARG}" \
          "${pdir}/bioproject.xml" "bioproject ${prj}" || true
  fi

  echo "   -> ${pdir}"
  echo
done

echo "Descarga de metadatos completa. Siguiente paso:"
echo "  python3 scripts/02_parse_metadata.py ${OUTDIR} docs/fase1_tablas.md"
