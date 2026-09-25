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
REGION_EN_TEXTO = re.compile(
    r"\bV[1-9](?:\s*[-–to]+\s*V?[1-9])?\b|\b(?:27F|338F|341F|515F|534R|785R|806R|926R|1492R)\b",
    re.IGNORECASE)


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
            if texto and REGION_EN_TEXTO.search(texto):
                ev.append(f"{etiqueta}: {texto[:300]}")
    # deduplicar conservando orden
    vistos, out = set(), []
    for e in ev:
        if e not in vistos:
            vistos.add(e)
            out.append(e)
    return out[:8]


def region_resumida(ev: list[str]) -> str:
    regiones = {m.group(0).upper().replace(" ", "")
                for e in ev for m in REGION_EN_TEXTO.finditer(e)}
    return ", ".join(sorted(regiones)) if regiones else "INDETERMINADA"


def analiza(prj: str, base: Path) -> dict:
    d = base / prj
    runs = lee_runinfo(d / "runinfo.csv")
    muestras = lee_biosample(d / "biosample.xml")
    exp_xml = lee_texto(d / "experiment.xml")

    res = {
        "prj": prj,
        "etiqueta": ETIQUETAS[prj],
        "n_runs": len(runs),
        "n_biosamples": len(muestras),
        "descargado": bool(runs or muestras),
        "plataforma": Counter(r.get("Platform", "?") for r in runs),
        "modelo": Counter(r.get("Model", "?") for r in runs),
        "layout": Counter(r.get("LibraryLayout", "?") for r in runs),
        "strategy": Counter(r.get("LibraryStrategy", "?") for r in runs),
        "spots": [int(r["spots"]) for r in runs if r.get("spots", "").isdigit()],
        "inventario": inventario_atributos(muestras),
        "estructura": {k: detecta(muestras, pn, pv) for k, (pn, pv) in ESTRUCTURA.items()},
        "auditoria": {k: detecta(muestras, pn, pv) for k, (pn, pv) in AUDITORIA.items()},
        "ev_region": evidencia_region(muestras, exp_xml),
    }
    res["region"] = region_resumida(res["ev_region"])

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
            ("Runs (SRR)", lambda r: r["n_runs"] or "—"),
            ("BioSamples", lambda r: r["n_biosamples"] or "—"),
            ("Sujetos (inferido)", lambda r: r["n_sujetos"] if r["n_sujetos"] else "no inferible"),
            ("Atributo sujeto", lambda r: f"`{r['clave_sujeto']}`" if r["clave_sujeto"] else "—"),
            ("Región 16S", lambda r: f"**{r['region']}**"),
            ("Plataforma", lambda r: top(r["plataforma"])),
            ("Modelo", lambda r: top(r["modelo"])),
            ("Layout", lambda r: top(r["layout"])),
            ("Estrategia", lambda r: top(r["strategy"])),
            ("Mediana spots/run", lambda r: (sorted(r["spots"])[len(r["spots"]) // 2]
                                             if r["spots"] else "—")),
        ]]))
    L.append("")

    L += ["## 2. Auditoría de ausencias", "",
          "Marcado **AUSENTE** = ningún atributo de BioSample con ese nombre o valor. "
          "Un atributo puede existir en el dataset analítico de los autores y aun así "
          "no estar depositado en BioSample: esta tabla sólo habla del registro público.", ""]
    claves_aud = list(AUDITORIA)
    L.append(tabla(
        ["Variable"] + [r["prj"] for r in resultados],
        [[k] + [("**AUSENTE**" if not r["auditoria"][k]
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
    regiones = {r["prj"]: r["region"] for r in resultados if r["descargado"]}
    if not regiones:
        L.append("_No evaluable: faltan metadatos._")
    else:
        distintas = set(regiones.values())
        for prj, reg in regiones.items():
            L.append(f"- {prj}: {reg}")
        L.append("")
        if "INDETERMINADA" in distintas:
            L.append("**No concluyente**: al menos un proyecto no declara la región. "
                     "Hay que confirmarla con la publicación o los primers antes de decidir.")
        elif len(distintas) > 1:
            L.append("**NO combinables**: las regiones del 16S difieren. "
                     "Los ASV no son comparables entre proyectos y los sesgos de "
                     "amplificación por taxón tampoco. Analizar por separado.")
        else:
            L.append("**Misma región declarada.** Condición necesaria pero no suficiente: "
                     "cebadores concretos, longitud de lectura y química pueden seguir "
                     "difiriendo, y eso rompe la equivalencia a nivel de ASV. "
                     "Un meta-análisis a nivel de género, con proyecto como efecto "
                     "aleatorio, es lo defendible; fusionar tablas de ASV no lo es.")
    L.append("")
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
