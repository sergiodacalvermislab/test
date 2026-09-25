# Comandos para ejecutar tú (Fase 1)

Este contenedor no tiene salida a `ncbi.nlm.nih.gov`. Estos comandos hacen la
extracción desde una máquina que sí la tenga. Elige **una** de las tres rutas.

Requisitos de A y B: `bash`, `curl` y `python3` (sólo stdlib — no hace falta
pandas, ni R, ni sra-tools, ni conda). Cualquier Mac o Linux los tiene.

---

## Ruta A — clonar y ejecutar (recomendada)

Es la que produce exactamente el formato que espero, y me devuelve los datos
por git sin que tengas que adjuntar nada.

```bash
git clone -b claude/microbioma-16s-parejas-mo2mxg \
  https://github.com/sergiodacalvermislab/test.git microbioma
cd microbioma

# 1) Descarga de metadatos (NO descarga secuencias: son unos pocos MB)
bash scripts/01_fetch_sra_metadata.sh data/raw_metadata

# 2) Inventario
python3 scripts/02_parse_metadata.py data/raw_metadata docs/fase1_tablas.md

# 3) Me lo devuelves por la rama
git add -f data/raw_metadata docs/fase1_tablas.md
git commit -m "Fase 1: metadatos SRA descargados"
git push
```

El `-f` del paso 3 es necesario: `data/raw_metadata/*/` está en `.gitignore`
para no versionar descargas por accidente, y aquí sí las queremos.

Cuando termine el paso 1 verás, por proyecto, cuántos runs y cuántas muestras
ha encontrado. Si algún proyecto sale con 0, dímelo: significa que está
embargado o que el accession no es el que pensamos, y es información útil de
por sí.

### Si tienes clave de NCBI (opcional)

Sube el límite de 3 a 10 peticiones/s. Se saca gratis en la cuenta de NCBI:

```bash
export NCBI_API_KEY=tu_clave
```

---

## Ruta B — sin clonar el repo

Un solo bloque, pégalo tal cual. Deja los ficheros en `./mb_meta/`.

```bash
mkdir -p mb_meta && cd mb_meta
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils

for P in PRJNA398590 PRJNA735440 PRJNA1248021; do
  echo "== $P =="
  for DB in sra biosample; do
    curl -s "$E/esearch.fcgi?db=$DB&term=${P}%5BBioProject%5D&usehistory=y&retmode=json" \
      -o "$P.$DB.search.json"
    read -r WE QK N <<<"$(python3 -c "
import json,sys
r=json.load(open('$P.$DB.search.json'))['esearchresult']
print(r.get('webenv',''), r.get('querykey',''), r.get('count','0'))")"
    echo "   $DB: $N registros"
    [ "$N" = "0" ] && continue
    if [ "$DB" = "sra" ]; then
      curl -s "$E/efetch.fcgi?db=sra&WebEnv=$WE&query_key=$QK&rettype=runinfo&retmode=text" \
        -o "$P.runinfo.csv"
      curl -s "$E/efetch.fcgi?db=sra&WebEnv=$WE&query_key=$QK&retmode=xml" \
        -o "$P.experiment.xml"
    else
      curl -s "$E/efetch.fcgi?db=biosample&WebEnv=$WE&query_key=$QK&rettype=full&retmode=xml" \
        -o "$P.biosample.xml"
    fi
    sleep 1
  done
done

ls -la && echo "---" && wc -l *.csv
```

Después me pasas el contenido de `mb_meta/` (o lo comprimes:
`tar czf mb_meta.tgz mb_meta/`).

---

## Ruta C — sin terminal, sólo navegador

La más rápida si no quieres tocar la línea de comandos. El `SraRunTable.csv`
del Run Selector trae los atributos de BioSample ya fusionados como columnas,
que es justo lo que necesito. El parser ya lo admite.

Para cada uno de los tres proyectos:

1. <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA398590>
2. <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA735440>
3. <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1248021>

En cada página: en el recuadro **Select**, fila *Total*, columna **Metadata**
→ botón **Download**. Baja un `SraRunTable.csv`.

Renómbralos y colócalos así (el nombre del directorio importa, el del fichero
también):

```
data/raw_metadata/PRJNA398590/SraRunTable.csv
data/raw_metadata/PRJNA735440/SraRunTable.csv
data/raw_metadata/PRJNA1248021/SraRunTable.csv
```

Y luego, o me los pasas tal cual, o si has clonado el repo:

```bash
python3 scripts/02_parse_metadata.py data/raw_metadata docs/fase1_tablas.md
```

Limitación de esta ruta: no trae `experiment.xml`, donde vive el
`LIBRARY_CONSTRUCTION_PROTOCOL` — el campo de texto libre donde los autores
suelen declarar los cebadores exactos. Si la región del 16S sale
`INDETERMINADA` en alguno de los tres, tendremos que ir a buscarla ahí o en el
paper, y es precisamente el dato que decide la combinabilidad. Por eso A o B
son preferibles si puedes.

---

## Qué necesito de vuelta, en orden de preferencia

1. Push a la rama con `data/raw_metadata/` y `docs/fase1_tablas.md` (ruta A).
2. El `.tgz` de `mb_meta/` (ruta B).
3. Los tres `SraRunTable.csv` (ruta C).

Con cualquiera de los tres cierro la Fase 1: tabla comparativa, auditoría de
ausencias y veredicto de combinabilidad.
