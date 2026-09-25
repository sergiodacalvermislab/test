#!/usr/bin/env python3
"""Fase 1 — reconstrucción del diseño experimental.

El inventario (02) dice qué atributos hay. Este script responde a lo que
decide la viabilidad: quién es quién, quién es pareja de quién, y qué
puntos temporales tiene cada cual.

Ninguno de los tres proyectos deposita un identificador de pareja. En dos de
ellos el emparejamiento se puede reconstruir del identificador de sujeto; en
el tercero no, y eso se reporta como tal en vez de forzarlo.

Salida: un TSV tidy por proyecto (base para la Fase 2) y un resumen por
consola.

Uso: python3 scripts/03_diseno.py [dir_metadatos] [dir_salida]
"""
from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

SITIO_735 = {"V": "vaginal", "P": "piel_peneana", "U": "orina"}
ORDEN_735 = ["D0", "D8", "W4", "W8", "End", "R"]
ORDEN_398 = ["D0", "D8", "D14", "D21", "D28"]


def carga(path: Path) -> list[dict]:
    if not path.exists():
        return []
    raiz = ET.parse(path).getroot()
    out = []
    for bs in raiz.iter("BioSample"):
        d = {"biosample": bs.get("accession", "")}
        for i in bs.iter("Id"):
            if i.get("db_label") == "Sample name":
                d["sample_name"] = (i.text or "").strip()
            elif i.get("db") == "SRA":
                d["sra_sample"] = (i.text or "").strip()
        for a in bs.iter("Attribute"):
            k = a.get("harmonized_name") or a.get("attribute_name")
            if k:
                d[k] = (a.text or "").strip()
        out.append(d)
    return out


def diseno_398590(ms: list[dict]) -> tuple[list[dict], list[str]]:
    """Serie 3xxx = mujeres (vagina); 4xxx = hombres (piel peneana, orina).

    La pareja es el sufijo común de las dos series. Se comprueba que las dos
    series tengan exactamente el mismo juego de sufijos antes de afirmarlo.
    """
    filas, notas = [], []
    series = defaultdict(set)
    for m in ms:
        if m.get("id"):
            series[m["id"][0]].add(m["id"][1:])
    if len(series) == 2:
        a, b = (series[k] for k in sorted(series))
        notas.append(f"series {sorted(series)}: {len(a)} y {len(b)} sujetos, "
                     f"{len(a & b)} sufijos comunes"
                     + ("" if a == b else f" — DESCUADRE: {sorted(a ^ b)}"))
    for m in ms:
        sid = m.get("id", "")
        filas.append({
            "biosample": m["biosample"], "sample_name": m.get("sample_name", ""),
            "sujeto": sid, "pareja": sid[1:] if sid else "",
            "sexo": {"3": "female", "4": "male"}.get(sid[:1], ""),
            "sitio": m.get("env_local_scale", ""), "punto": m.get("day", ""),
            "fecha": m.get("collection_date", ""), "control": "",
        })
    return filas, notas


def diseno_735440(ms: list[dict]) -> tuple[list[dict], list[str]]:
    """host_subject_id = <nº pareja><F|M>; sample_name = <nº><V|P|U>_<punto>."""
    filas, notas = [], []
    for m in ms:
        sid = m.get("host_subject_id", "")
        g = re.fullmatch(r"(\d+)([FM])", sid)
        n = re.fullmatch(r"(\d+)([VPU])_(\w+)", m.get("sample_name", "") or "")
        es_control = not sid
        filas.append({
            "biosample": m["biosample"], "sample_name": m.get("sample_name", ""),
            "sujeto": sid, "pareja": g.group(1) if g else "",
            "sexo": m.get("host_sex", ""),
            "sitio": SITIO_735.get(n.group(2), "") if n else "",
            "punto": n.group(3) if n else "",
            "fecha": m.get("collection_date", ""),
            "control": m.get("isolation_source", "") if es_control else "",
        })
    ctrl = [f for f in filas if f["control"]]
    notas.append(f"{len(ctrl)} muestras de control/sin sujeto: "
                 + ", ".join(f"{k}×{v}" for k, v in
                             Counter(f["control"][:45] for f in ctrl).most_common(6)))
    comp = defaultdict(set)
    for f in filas:
        if f["pareja"]:
            comp[f["pareja"]].add(f["sexo"][:1].upper())
    completas = sum(1 for v in comp.values() if v == {"F", "M"})
    notas.append(f"{len(comp)} números de pareja; {completas} con ambos miembros")
    return filas, notas


def diseno_1248021(ms: list[dict]) -> tuple[list[dict], list[str]]:
    """host_subject_id = #<nº persona>_<código de punto>.

    El número identifica a una PERSONA, no a una pareja: cada número aparece
    con un solo sexo. El identificador de pareja (PID en el script de análisis
    de los autores) no está depositado, así que no se puede emparejar.
    """
    filas, notas = [], []
    por_num = defaultdict(set)
    for m in ms:
        g = re.match(r"#?(\d+)_([A-Z]\d?)", m.get("host_subject_id", "") or "")
        if g:
            por_num[g.group(1)].add(m.get("host_sex", ""))
        filas.append({
            "biosample": m["biosample"], "sample_name": m.get("sample_name", ""),
            "sujeto": g.group(1) if g else "", "pareja": "",
            "sexo": m.get("host_sex", ""),
            "sitio": m.get("env_local_scale", ""),
            "punto": g.group(2) if g else "",
            "fecha": m.get("collection_date", ""), "control": "",
            "circuncision": m.get("host_phenotype", ""),
            "edad": m.get("host_age", ""), "etnia": m.get("ethnicity", ""),
        })
    mixtos = sum(1 for v in por_num.values() if len(v) > 1)
    notas.append(f"{len(por_num)} números de sujeto; {mixtos} con más de un sexo")
    notas.append("Emparejamiento de parejas NO recuperable: cada número es una "
                 "persona y el identificador de pareja no está depositado.")
    return filas, notas


HANDLERS = {"PRJNA398590": (diseno_398590, ORDEN_398),
            "PRJNA735440": (diseno_735440, ORDEN_735),
            "PRJNA1248021": (diseno_1248021, None)}


def matriz(filas: list[dict], orden: list[str] | None) -> str:
    reales = [f for f in filas if not f["control"]]
    sitios = sorted({f["sitio"] for f in reales if f["sitio"]})
    puntos = orden or sorted({f["punto"] for f in reales if f["punto"]})
    c = Counter((f["sitio"], f["punto"]) for f in reales)
    out = [f"{'sitio':<15}" + "".join(f"{p:>7}" for p in puntos) + f"{'total':>8}"]
    for s in sitios:
        tot = sum(v for (ss, _), v in c.items() if ss == s)
        out.append(f"{s:<15}" + "".join(f"{c[(s, p)]:>7}" for p in puntos) + f"{tot:>8}")
    return "\n".join(out)


def cobertura(filas: list[dict], primero: str, ultimo: str) -> list[str]:
    """Cuántos sujetos y parejas cubren el intervalo completo."""
    out = []
    por_sitio = defaultdict(lambda: defaultdict(set))
    for f in filas:
        if f["sujeto"] and f["sitio"]:
            por_sitio[f["sitio"]][f["sujeto"]].add(f["punto"])
    for s, d in sorted(por_sitio.items()):
        completos = sum(1 for v in d.values() if {primero, ultimo} <= v)
        out.append(f"  {s:<15} {len(d):>3} sujetos, {completos:>3} con {primero} y {ultimo}"
                   f"  (nº de puntos: {dict(sorted(Counter(len(v) for v in d.values()).items()))})")
    parejas = defaultdict(lambda: defaultdict(set))
    for f in filas:
        if f["pareja"]:
            parejas[f["pareja"]][f["sexo"][:1].upper()].add(f["punto"])
    ok = sum(1 for v in parejas.values()
             if len(v) == 2 and all({primero, ultimo} <= p for p in v.values()))
    if parejas:
        out.append(f"  parejas con {primero}+{ultimo} en AMBOS miembros: {ok}/{len(parejas)}")
    return out


def main() -> int:
    base = Path(sys.argv[1] if len(sys.argv) > 1 else "data/raw_metadata")
    out = Path(sys.argv[2] if len(sys.argv) > 2 else "data/diseno")
    out.mkdir(parents=True, exist_ok=True)
    for prj, (fn, orden) in HANDLERS.items():
        ms = carga(base / prj / "biosample.xml")
        print("=" * 72); print(prj); print("=" * 72)
        if not ms:
            print("  sin biosample.xml\n"); continue
        filas, notas = fn(ms)
        cols = sorted({k for f in filas for k in f})
        dest = out / f"{prj}_diseno.tsv"
        with dest.open("w", encoding="utf-8") as fh:
            fh.write("\t".join(cols) + "\n")
            for f in filas:
                fh.write("\t".join(str(f.get(c, "")) for c in cols) + "\n")
        print(f"  {len(filas)} muestras -> {dest}")
        for n in notas:
            print(f"  · {n}")
        print(); print(matriz(filas, orden))
        if orden:
            print("\nCobertura longitudinal:")
            print("\n".join(cobertura(filas, orden[0], orden[4] if len(orden) > 4 else orden[-1])))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
