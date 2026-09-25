#!/usr/bin/env python3
"""Fase 1 — inventario de metadatos de los tres BioProjects.

Lee lo que descargó 01_fetch_sra_metadata.sh y produce:
  1. Resumen por proyecto (runs, muestras, plataforma, librería).
  2. Inventario completo de atributos BioSample con tasa de relleno.
  3. Inferencia de estructura del diseño: sujeto, punto temporal, sitio
     anatómico, emparejamiento de parejas.
  4. Región del 16S, a partir de la evidencia textual cruda del registro.
  5. Auditoría explícita de AUSENCIAS (carga qPCR, Nugent/Amsel por punto
     temporal, fase del ciclo, circuncisión).
  6. Veredicto de combinabilidad.

Sólo stdlib. No infiere nada que no esté en los ficheros: cada afirmación
va acompañada del nombre de atributo y los valores que la sostienen.

Uso: python3 scripts/02_parse_metadata.py [dir_metadatos] [salida.md]
"""
from __future__ import annotations

import csv
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

PROJECTS = ["PRJNA398590", "PRJNA735440", "PRJNA1248021"]

ETIQUETAS = {
    "PRJNA398590": "Plummer 2018 — parejas, tratamiento antibiótico, 4 semanas",
    "PRJNA735440": "Plummer 2021 — 34 parejas, 12 semanas, 3 sitios",
    "PRJNA1248021": "Park 2026 — dinámica post-coital (SECS, Toronto)",
}

# --- Heurísticas de detección -------------------------------------------------
# (clave legible, patrones sobre el NOMBRE del atributo, patrones sobre el VALOR)
# Se informa siempre del atributo concreto que disparó la detección, para que
# el resultado sea verificable y no una afirmación a ciegas.

ESTRUCTURA = {
    "sujeto": ([r"subject", r"host_subject", r"patient", r"participant", r"\bpid\b",
                r"individual", r"donor"], []),
    "punto_temporal": ([r"time", r"visit", r"week", r"day", r"timepoint", r"collection_date",
                        r"month", r"follow", r"baseline"], []),
    "sitio_anatomico": ([r"body_site", r"isolation_source", r"tissue", r"source_ma",
                         r"body_?habitat", r"\bsite\b", r"organism_part"], []),
    "pareja": ([r"partner", r"couple", r"dyad", r"\bpair\b", r"\bpid\b", r"spouse"], []),
    "sexo": ([r"\bsex\b", r"gender", r"host_sex"], []),
}

AUDITORIA = {
    "Carga bacteriana (qPCR / BactQuant)": (
        [r"qpcr", r"bacterial_?load", r"total_?bacteri", r"copy_?number", r"copies",
         r"16s_?quant", r"bactquant", r"std_?qty", r"abundance_?absolute", r"quantif"],
        [r"copies/?m[lL]", r"gene copies"],
    ),
    "Estado de VB (Nugent / Amsel)": (
        [r"nugent", r"amsel", r"\bbv\b", r"bv_?status", r"vaginos", r"bvstatus"],
        [r"\bnugent\b", r"bacterial vaginosis"],
    ),
    "Fase del ciclo menstrual": (
        [r"menstrua", r"cycle", r"cyclephase", r"\blmp\b", r"last_?menstrual", r"menses"],
        [r"follicular", r"luteal", r"menses"],
    ),
    "Estado de circuncisión": (
        [r"circumcis", r"\bcirc\b"],
        [r"circumcised", r"uncircumcised"],
    ),
    "Tratamiento / brazo del ensayo": (
        [r"treatment", r"antibiotic", r"clindamycin", r"metronidazol", r"\barm\b",
         r"therapy", r"drug", r"placebo", r"intervention"],
        [r"clindamycin", r"metronidazole", r"placebo"],
    ),
    "Uso de preservativo / actividad sexual": (
        [r"condom", r"sexual", r"coitus", r"intercourse", r"sex_?type", r"abstin"],
        [],
    ),
    "pH vaginal": ([r"\bph\b", r"vaginal_?ph"], []),
}

REGION_ATTRS = [r"target_?gene", r"target_?subfragment", r"pcr_?primers", r"primer",
                r"amplicon", r"\bregion\b", r"variable_?region"]

# Tres cosas distintas, que antes se mezclaban en una sola cadena y producían
# falsos "las regiones difieren" cuando lo único que cambiaba era si el
# registro citaba o no el nombre del cebador:
#   - la región variable      -> "V3-V4"
#   - el nombre del cebador   -> "341F", "805R"
#   - la secuencia del cebador-> "CCTACGGGNGGCWGCAG"
# La secuencia es la que manda: es lo que fija los extremos del amplicón y,
# por tanto, si dos proyectos producen o no los mismos ASV.
REGION_V = re.compile(r"\bV[1-9](?:\s*[-–]\s*V?[1-9])?\b", re.IGNORECASE)
PRIMER_NOMBRE = re.compile(r"\b\d{2,4}[FR]\b")
PRIMER_SEQ = re.compile(r"\b[ACGTRYSWKMBDHVN]{15,30}\b")


def normaliza_region(txt: str) -> str:
    """'V3 - V4', 'v3-4', 'V3-V4' -> 'V3-V4'."""
    t = txt.upper().replace(" ", "").replace("–", "-")
    m = re.fullmatch(r"V([1-9])-V?([1-9])", t)
    return f"V{m.group(1)}-V{m.group(2)}" if m else t


def busca(patrones: list[str], texto: str) -> bool:
    return any(re.search(p, texto, re.IGNORECASE) for p in patrones)


# --- Lectura ------------------------------------------------------------------

def lee_runinfo(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        return [fila for fila in csv.DictReader(fh) if fila.get("Run")]


def lee_biosample(path: Path) -> list[dict]:
    """Devuelve una lista de dicts: accession, title, organism + atributos."""
    if not path.exists():
        return []
    try:
        raiz = ET.parse(path).getroot()
    except ET.ParseError as exc:
        print(f"  ! {path}: XML ilegible ({exc})", file=sys.stderr)
        return []
    muestras = []
    for bs in raiz.iter("BioSample"):
        reg = {"_accession": bs.get("accession", ""), "_id": bs.get("id", "")}
        titulo = bs.find("./Description/Title")
        reg["_title"] = (titulo.text or "").strip() if titulo is not None else ""
        org = bs.find("./Description/Organism/OrganismName")
        reg["_organism"] = (org.text or "").strip() if org is not None else ""
        for at in bs.iter("Attribute"):
            nombre = (at.get("harmonized_name") or at.get("attribute_name") or "").strip()
            if nombre:
                reg[nombre] = (at.text or "").strip()
        muestras.append(reg)
    return muestras


def lee_texto(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""


# Columnas de SraRunTable.csv que describen el run/la secuenciación, no la
# biología de la muestra. Se excluyen del inventario de atributos para que la
# tabla no se llene de ruido técnico, pero se siguen usando para plataforma,
# modelo y layout.
COLS_RUN = {
    "run", "releasedate", "loaddate", "createdate", "spots", "bases", "spots_with_mates",
    "avglength", "size_mb", "assemblyname", "download_path", "experiment", "libraryname",
    "librarystrategy", "libraryselection", "librarysource", "librarylayout", "insertsize",
    "insertdev", "platform", "model", "srastudy", "bioproject", "study_pubmed_id",
    "projectid", "sample", "biosample", "samplename", "submission", "consent",
    "runhash", "readhash", "center_name", "instrument", "assay_type", "datastore_filetype",
    "datastore_provider", "datastore_region", "bytes", "avgspotlen", "sra_study",
    "sample_name", "experiment_title", "library_name", "sra_accession", "version",
}


def lee_srarun_table(path: Path) -> tuple[list[dict], list[dict]]:
    """Lee un SraRunTable.csv del SRA Run Selector.

    Ese CSV trae una fila por run con los atributos de BioSample ya fusionados
    como columnas, que es justo lo que necesita la Fase 1. Devuelve
    (runs, muestras) deduplicando muestras por BioSample.
    """
    if not path.exists():
        return [], []
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        filas = [f for f in csv.DictReader(fh) if any((v or "").strip() for v in f.values())]
    if not filas:
        return [], []

    # Normaliza los nombres de columna que el Run Selector escribe distinto
    # que runinfo, para que el resto del script no tenga que saberlo.
    alias = {"Assay Type": "LibraryStrategy", "Instrument": "Model",
             "Platform": "Platform", "LibraryLayout": "LibraryLayout",
             "AvgSpotLen": "avgspotlen", "Bases": "bases", "Bytes": "bytes"}
    runs = []
    for f in filas:
        r = dict(f)
        for origen, destino in alias.items():
            if origen in f and destino not in r:
                r[destino] = f[origen]
        r.setdefault("Run", f.get("Run") or f.get("run") or "")
        runs.append(r)

    vistos, muestras = set(), []
    for f in filas:
        acc = (f.get("BioSample") or f.get("biosample") or f.get("Sample Name")
               or f.get("Run") or "").strip()
        if not acc or acc in vistos:
            continue
        vistos.add(acc)
        reg = {"_accession": acc, "_id": "", "_title": (f.get("Sample Name") or "").strip(),
               "_organism": (f.get("Organism") or f.get("organism") or "").strip()}
        for k, v in f.items():
            if k and k.strip().lower() not in COLS_RUN and (v or "").strip():
                reg[k.strip()] = v.strip()
        muestras.append(reg)
    return runs, muestras


# --- Análisis -----------------------------------------------------------------

def inventario_atributos(muestras: list[dict]) -> list[tuple]:
    """(atributo, n_relleno, pct, n_distintos, ejemplos) ordenado por relleno."""
    total = len(muestras) or 1
    claves = sorted({k for m in muestras for k in m if not k.startswith("_")})
    filas = []
    for k in claves:
        vals = [m[k] for m in muestras if m.get(k) not in (None, "", "not applicable",
                                                           "missing", "not collected", "NA")]
        distintos = sorted(set(vals))
        ejemplos = ", ".join(distintos[:4]) + (" …" if len(distintos) > 4 else "")
        filas.append((k, len(vals), 100.0 * len(vals) / total, len(distintos), ejemplos))
    filas.sort(key=lambda f: (-f[1], f[0]))
    return filas


def detecta(muestras: list[dict], patrones_n: list[str], patrones_v: list[str]) -> list[tuple]:
    """Atributos que casan por nombre o por valor. -> (nombre, n_relleno, ejemplos)"""
    total = len(muestras) or 1
    hits = []
    claves = sorted({k for m in muestras for k in m if not k.startswith("_")})
    for k in claves:
        vals = [m[k] for m in muestras if m.get(k) not in (None, "", "not applicable",
                                                           "missing", "not collected", "NA")]
        if not vals:
            continue
        por_nombre = busca(patrones_n, k)
        por_valor = bool(patrones_v) and any(busca(patrones_v, v) for v in vals[:200])
        if por_nombre or por_valor:
            distintos = sorted(set(vals))
            hits.append((k, len(vals), 100.0 * len(vals) / total,
                         ", ".join(distintos[:5]) + (" …" if len(distintos) > 5 else "")))
    return hits


def evidencia_region(muestras: list[dict], exp_xml: str) -> list[str]:
    """Cadenas crudas que documentan la región 16S. Sin adivinar."""
    ev = []
    for m in muestras[:200]:
        for k, v in m.items():
            if not k.startswith("_") and busca(REGION_ATTRS, k) and v:
                ev.append(f"BioSample `{k}` = {v}")
    for etiqueta in ("DESIGN_DESCRIPTION", "LIBRARY_CONSTRUCTION_PROTOCOL"):
        for trozo in re.findall(rf"<{etiqueta}>(.*?)</{etiqueta}>", exp_xml, re.S):
            texto = re.sub(r"\s+", " ", trozo).strip()
            if texto and (REGION_V.search(texto) or PRIMER_SEQ.search(texto)):
                ev.append(f"{etiqueta}: {texto[:400]}")
    # deduplicar conservando orden
    vistos, out = set(), []
    for e in ev:
        if e not in vistos:
            vistos.add(e)
            out.append(e)
    return out[:8]


def region_resumida(ev: list[str]) -> str:
    regiones = {normaliza_region(m.group(0)) for e in ev for m in REGION_V.finditer(e)}
    return ", ".join(sorted(regiones)) if regiones else "INDETERMINADA"


def primers_resumidos(ev: list[str]) -> tuple[set[str], set[str]]:
    """(secuencias, nombres) de cebador halladas en la evidencia."""
    seqs, nombres = set(), set()
    for e in ev:
        # Recorta la etiqueta del campo para no capturar texto del prefijo.
        cuerpo = e.split(":", 1)[1] if ":" in e else e
        for m in PRIMER_SEQ.finditer(cuerpo):
            s = m.group(0).upper()
            # Un tramo de sólo ACGT largo suele ser secuencia real; exigir al
            # menos una base degenerada o una longitud típica de cebador evita
            # capturar palabras en mayúsculas o identificadores.
            if 15 <= len(s) <= 30:
                seqs.add(s)
        nombres.update(PRIMER_NOMBRE.findall(cuerpo))
    return seqs, nombres


def nombres_de_muestra(runs: list[dict]) -> dict:
    """Radiografía de SampleName / LibraryName.

    Cuando el BioSample no trae atributos estructurados, el nombre de muestra
    es lo único que queda para reconstruir el diseño. No se adivina qué
    significa cada campo: se reporta el separador dominante, cuántos campos
    hay, y el repertorio de valores de cada campo, que es lo que permite
    reconocer a ojo cuál es el sujeto, cuál la visita y cuál el sitio.
    """
    nombres = [(r.get("SampleName") or r.get("Sample Name") or "").strip() for r in runs]
    nombres = [n for n in nombres if n]
    if not nombres:
        nombres = [(r.get("LibraryName") or "").strip() for r in runs]
        nombres = [n for n in nombres if n]
    if not nombres:
        return {}

    sep = max("-_.", key=lambda s: sum(n.count(s) for n in nombres))
    partes = [n.split(sep) for n in nombres]
    n_campos = Counter(len(p) for p in partes)
    modal = n_campos.most_common(1)[0][0]
    campos = []
    for i in range(modal):
        vals = [p[i] for p in partes if len(p) == modal]
        distintos = sorted(set(vals))
        campos.append({
            "pos": i + 1,
            "n_distintos": len(distintos),
            "ejemplos": ", ".join(distintos[:12]) + (" …" if len(distintos) > 12 else ""),
        })
    return {
        "n": len(nombres),
        "ejemplos": nombres[:8],
        "separador": sep,
        "reparto_campos": dict(n_campos),
        "campos": campos,
    }


def analiza(prj: str, base: Path) -> dict:
    d = base / prj
    runs = lee_runinfo(d / "runinfo.csv")
    muestras = lee_biosample(d / "biosample.xml")
    exp_xml = lee_texto(d / "experiment.xml")

    # Alternativa sin E-utilities: SraRunTable.csv descargado del Run Selector,
    # que ya trae los atributos de BioSample fusionados como columnas. Rellena
    # sólo lo que falte, para poder mezclar ambas fuentes.
    fuentes = ["runinfo.csv" if runs else None, "biosample.xml" if muestras else None,
               "experiment.xml" if exp_xml else None]
    if not runs or not muestras:
        srt_runs, srt_muestras = lee_srarun_table(d / "SraRunTable.csv")
        if srt_runs or srt_muestras:
            fuentes.append("SraRunTable.csv")
        runs = runs or srt_runs
        muestras = muestras or srt_muestras

    res = {
        "prj": prj,
        "etiqueta": ETIQUETAS[prj],
        "fuentes": [f for f in fuentes if f],
        "n_runs": len(runs),
        "n_biosamples": len(muestras),
        "descargado": bool(runs or muestras),
        "plataforma": Counter(r.get("Platform", "?") for r in runs),
        "modelo": Counter(r.get("Model", "?") for r in runs),
        "layout": Counter(r.get("LibraryLayout", "?") for r in runs),
        "strategy": Counter(r.get("LibraryStrategy", "?") for r in runs),
        # runinfo trae "spots"; SraRunTable no, así que se cae a "bases".
        "spots": [int(r["spots"]) for r in runs if (r.get("spots") or "").isdigit()],
        "bases": [int(r["bases"]) for r in runs if (r.get("bases") or "").isdigit()],
        "inventario": inventario_atributos(muestras),
        "estructura": {k: detecta(muestras, pn, pv) for k, (pn, pv) in ESTRUCTURA.items()},
        "auditoria": {k: detecta(muestras, pn, pv) for k, (pn, pv) in AUDITORIA.items()},
        "ev_region": evidencia_region(muestras, exp_xml),
    }
    res["region"] = region_resumida(res["ev_region"])
    res["primer_seqs"], res["primer_nombres"] = primers_resumidos(res["ev_region"])

    # Los nombres de muestra suelen codificar el diseño (sujeto, visita, sitio)
    # aunque el BioSample no traiga atributos estructurados. No se interpreta el
    # código: se muestran ejemplos y el reparto de campos para que lo descifre
    # una persona.
    res["nombres_muestra"] = nombres_de_muestra(runs)

    # Nº de sujetos y puntos temporales, si hay un atributo de sujeto usable.
    hits_suj = res["estructura"]["sujeto"]
    if hits_suj:
        clave = max(hits_suj, key=lambda h: h[1])[0]
        por_sujeto = defaultdict(int)
        for m in muestras:
            if m.get(clave):
                por_sujeto[m[clave]] += 1
        res["clave_sujeto"] = clave
        res["n_sujetos"] = len(por_sujeto)
        res["muestras_por_sujeto"] = Counter(por_sujeto.values())
    else:
        res["clave_sujeto"] = None
        res["n_sujetos"] = None
        res["muestras_por_sujeto"] = Counter()
    return res


# --- Informe ------------------------------------------------------------------

def tabla(cabeceras: list[str], filas: list[list[str]]) -> str:
    out = ["| " + " | ".join(cabeceras) + " |",
           "|" + "|".join("---" for _ in cabeceras) + "|"]
    for f in filas:
        out.append("| " + " | ".join(str(c).replace("|", "\\|") for c in f) + " |")
    return "\n".join(out)


def top(c: Counter, n: int = 3) -> str:
    return ", ".join(f"{k} ({v})" for k, v in c.most_common(n)) if c else "—"


def informe(resultados: list[dict]) -> str:
    L = ["# Fase 1 — Inventario de metadatos (generado automáticamente)", ""]
    faltan = [r["prj"] for r in resultados if not r["descargado"]]
    if faltan:
        L += [f"> **Sin datos descargados para: {', '.join(faltan)}.** "
              "Ejecuta `scripts/01_fetch_sra_metadata.sh` con salida de red hacia NCBI.", ""]

    L += ["## 1. Tabla comparativa", ""]
    L.append(tabla(
        ["", "PRJNA398590", "PRJNA735440", "PRJNA1248021"],
        [[campo] + [fn(r) for r in resultados] for campo, fn in [
            ("Estudio", lambda r: r["etiqueta"]),
            ("Fuente de los datos", lambda r: ", ".join(f"`{f}`" for f in r["fuentes"]) or "—"),
            ("Runs (SRR)", lambda r: r["n_runs"] or "—"),
            ("BioSamples", lambda r: r["n_biosamples"] or "—"),
            ("Sujetos (inferido)", lambda r: r["n_sujetos"] if r["n_sujetos"] else "no inferible"),
            ("Atributo sujeto", lambda r: f"`{r['clave_sujeto']}`" if r["clave_sujeto"] else "—"),
            ("Región 16S", lambda r: f"**{r['region']}**"),
            ("Plataforma", lambda r: top(r["plataforma"])),
            ("Modelo", lambda r: top(r["modelo"])),
            ("Layout", lambda r: top(r["layout"])),
            ("Estrategia", lambda r: top(r["strategy"])),
            ("Profundidad mediana/run", lambda r: (
                f"{sorted(r['spots'])[len(r['spots']) // 2]:,} spots" if r["spots"]
                else f"{sorted(r['bases'])[len(r['bases']) // 2] / 1e6:,.0f} Mbases"
                if r["bases"] else "—")),
        ]]))
    L.append("")

    L += ["## 2. Auditoría de ausencias", "",
          "**AUSENTE** = hay registros BioSample y ninguno declara ese atributo. "
          "**SIN DATOS** = no se han podido leer registros BioSample, así que no hay "
          "nada que auditar; no es un hallazgo, es una laguna de la descarga.", "",
          "Un atributo puede existir en el dataset analítico de los autores y aun así "
          "no estar depositado en BioSample: esta tabla sólo habla del registro público.", ""]
    # Sin BioSamples no se puede afirmar ausencia: sería confundir "no lo he
    # mirado" con "no está". Esa confusión es exactamente lo que invalidaría
    # todo el juicio de viabilidad.
    sin_bs = [r["prj"] for r in resultados if r["descargado"] and not r["n_biosamples"]]
    if sin_bs:
        L += [f"> Sin registros BioSample legibles para **{', '.join(sin_bs)}**. "
              "Su columna queda en SIN DATOS: no se puede concluir ausencia de nada.", ""]
    claves_aud = list(AUDITORIA)
    L.append(tabla(
        ["Variable"] + [r["prj"] for r in resultados],
        [[k] + [("SIN DATOS" if not r["n_biosamples"]
                 else "**AUSENTE**" if not r["auditoria"][k]
                 else ", ".join(f"`{h[0]}` ({h[2]:.0f}%)" for h in r["auditoria"][k]))
                for r in resultados] for k in claves_aud]))
    L.append("")

    L += ["## 3. Estructura del diseño detectada", ""]
    for r in resultados:
        L += [f"### {r['prj']} — {r['etiqueta']}", ""]
        if not r["descargado"]:
            L += ["_Sin metadatos descargados._", ""]
            continue
        for eje, hits in r["estructura"].items():
            if hits:
                L.append(f"- **{eje}**: " + "; ".join(
                    f"`{h[0]}` ({h[1]} muestras, {h[2]:.0f}%) → {h[3]}" for h in hits))
            else:
                L.append(f"- **{eje}**: no detectado")
        if r["muestras_por_sujeto"]:
            reparto = ", ".join(f"{n} muestras×{c} sujetos"
                                for n, c in sorted(r["muestras_por_sujeto"].items()))
            L.append(f"- **muestras por sujeto**: {reparto}")
        nm = r.get("nombres_muestra") or {}
        if nm:
            L += ["", f"Nombres de muestra ({nm['n']} runs, separador `{nm['separador']}`, "
                      f"reparto de campos {nm['reparto_campos']}):", "",
                  "  - ejemplos: " + ", ".join(f"`{e}`" for e in nm["ejemplos"])]
            for c in nm["campos"]:
                L.append(f"  - campo {c['pos']}: {c['n_distintos']} valores → {c['ejemplos']}")
        L += ["", "Evidencia de región 16S:", ""]
        L += [f"  - {e}" for e in r["ev_region"]] or ["  - _ninguna evidencia textual_"]
        L.append("")

    L += ["## 4. Inventario completo de atributos BioSample", ""]
    for r in resultados:
        L += [f"### {r['prj']}", ""]
        if not r["inventario"]:
            L += ["_Sin atributos._", ""]
            continue
        L.append(tabla(["Atributo", "n", "% relleno", "n distintos", "ejemplos"],
                       [[k, n, f"{p:.0f}%", d, ej] for k, n, p, d, ej in r["inventario"]]))
        L.append("")

    L += ["## 5. Combinabilidad", ""]
    vivos = [r for r in resultados if r["descargado"]]
    if not vivos:
        L += ["_No evaluable: faltan metadatos._", ""]
        return "\n".join(L)

    L.append(tabla(["Proyecto", "Región", "Cebadores (nombres)", "Cebadores (secuencia)"],
                   [[r["prj"], r["region"],
                     ", ".join(sorted(r["primer_nombres"])) or "—",
                     "<br>".join(sorted(r["primer_seqs"])) or "**no declarada**"]
                    for r in vivos]))
    L.append("")

    regiones = {r["region"] for r in vivos}
    # Se comparan por pares: la combinabilidad no es una propiedad del conjunto.
    L += ["Veredicto por pares (la combinabilidad no es global: dos proyectos "
          "pueden ser fusionables entre sí y ninguno de los dos con el tercero).", ""]
    filas = []
    for i, a in enumerate(vivos):
        for b in vivos[i + 1:]:
            if "INDETERMINADA" in (a["region"], b["region"]):
                v = "No concluyente — región sin declarar"
            elif a["region"] != b["region"]:
                v = f"**NO** — regiones distintas ({a['region']} vs {b['region']})"
            elif not a["primer_seqs"] or not b["primer_seqs"]:
                v = ("Misma región, cebadores no declarados en al menos uno. "
                     "Sólo a nivel de género, y con cautela")
            elif a["primer_seqs"] == b["primer_seqs"]:
                v = ("**SÍ a nivel de ASV** — misma región y cebadores idénticos. "
                     "Sigue haciendo falta procesar junto (un solo `learnErrors` "
                     "por run de secuenciación) y modelar el estudio como efecto "
                     "de lote")
            else:
                v = ("**NO a nivel de ASV**, sí a nivel de género — misma región "
                     "pero cebadores distintos: los amplicones tienen extremos "
                     "distintos, luego el mismo organismo da ASV distintos")
            filas.append([f"{a['prj']} ↔ {b['prj']}", v])
    L.append(tabla(["Par", "Veredicto"], filas))
    L += ["", "Nota sobre el criterio: la región variable por sí sola no decide. "
          "Lo que fija si dos proyectos producen el mismo ASV para el mismo "
          "organismo son las secuencias exactas de los cebadores, porque son "
          "las que determinan dónde empieza y acaba el amplicón una vez "
          "recortados. Misma región con cebadores distintos obliga a agregar a "
          "género antes de comparar.", ""]
    if len(regiones) == 1 and "INDETERMINADA" not in regiones:
        L += [f"(Los tres proyectos declaran {regiones.pop()}.)", ""]
    return "\n".join(L)


def main() -> int:
    base = Path(sys.argv[1] if len(sys.argv) > 1 else "data/raw_metadata")
    destino = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if not base.exists():
        print(f"No existe {base}. Ejecuta antes scripts/01_fetch_sra_metadata.sh",
              file=sys.stderr)
        return 1
    resultados = [analiza(p, base) for p in PROJECTS]
    texto = informe(resultados)
    if destino:
        destino.parent.mkdir(parents=True, exist_ok=True)
        destino.write_text(texto, encoding="utf-8")
        print(f"Escrito: {destino}")
    else:
        print(texto)
    return 0


if __name__ == "__main__":
    sys.exit(main())
