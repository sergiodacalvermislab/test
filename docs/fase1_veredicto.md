# Fase 1 — Veredicto

Basado en los metadatos reales de los tres BioProjects (140 + 418 + 284
registros BioSample, descargados y en `data/raw_metadata/`). Todo lo que sigue
sale de esos ficheros; lo que no se puede sostener con ellos está marcado como
desconocido.

## Veredicto en tres líneas

- **PRJNA735440 es un dataset excelente para tu pregunta principal.** Diseño
  completo y equilibrado: 27 parejas, las 27 con los dos miembros y con los
  tres sitios, serie temporal completa D0→End.
- **No hay una sola variable clínica en ninguno de los tres**, salvo
  circuncisión en el de Park. Ni Nugent, ni carga bacteriana, ni fase del
  ciclo. Eso no impide la Fase 2, pero cambia qué preguntas son respondibles.
- **Los dos Plummer son combinables entre sí a nivel de ASV**; ninguno lo es
  con el de Park.

---

## 1. Tabla comparativa

| | PRJNA398590 | PRJNA735440 | PRJNA1248021 |
|---|---|---|---|
| Estudio | Tratamiento oral+tópico a la pareja masculina (aceptabilidad/tolerabilidad) | Piloto de tratamiento concurrente a la pareja, 12 semanas | Dinámica post-coital (SECS) |
| Centro | Murdoch Childrens RI, Melbourne | Murdoch Childrens RI, Melbourne | George Washington Univ. |
| Runs / BioSamples | 140 / 140 | 418 / 418 | 284 / 284 |
| Muestras reales | 140 | **375** (+43 controles) | 284 |
| Sujetos | 42 | 54 | 73 (37 H, 36 M) |
| **Parejas** | **21, reconstruibles** | **27, reconstruibles** | **no reconstruibles** |
| Puntos temporales | D0, D8, D14, D21, D28 | D0, D8, W4, W8, End (+R) | A, E, H, P (+repeticiones) |
| Sitios | vagina, piel peneana, orina | vaginal, piel peneana, orina | Cervical, CorSul |
| Región 16S | V3–V4 | V3–V4 | V3–V4 |
| Cebadores | 341F/805R | 341F/805R | BactQuant + espaciadores |
| Plataforma | MiSeq, PAIRED | MiSeq, PAIRED | MiSeq, PAIRED |
| Profundidad mediana | 35.373 spots | 35.712 | 23.950 |
| Controles negativos | ninguno depositado | **43** | ninguno depositado |

### Cómo se reconstruyen las parejas

Ninguno deposita un identificador de pareja. En dos se deduce:

- **PRJNA398590**: el atributo `id` tiene dos series, `3xxx` (21 sujetos, todos
  vaginales) y `4xxx` (21 sujetos, piel peneana y orina). Las dos series
  comparten **exactamente** los 21 sufijos, sin descuadre. La pareja es el
  sufijo; la serie es el sexo.
- **PRJNA735440**: `host_subject_id` es `<nº><F|M>` (`12F`, `12M`). 27 números,
  los 27 con ambos miembros.
- **PRJNA1248021**: `host_subject_id` es `#<nº>_<código>`, y cada número
  aparece con **un solo sexo** (37 hombres, 36 mujeres, 0 mixtos). El número es
  una persona, no una pareja. El `PID` que usa el script SAS de los autores no
  está depositado. **El emparejamiento no es recuperable de los datos
  públicos.**

---

## 2. Matrices de diseño

### PRJNA735440 — el bueno

```
sitio               D0     D8     W4     W8    End      R   total
orina               27     27     23     18     27      0     122
piel_peneana        27     27     23     18     27      0     122
vaginal             27     27     26     21     27      3     131
```

Cobertura longitudinal: **27/27 hombres con D0 y End en piel peneana y en
orina. 27/27 parejas con D0+End en ambos miembros.** 18 de 27 tienen los cinco
puntos. La pérdida se concentra en W4 y W8, no en los extremos — es decir, el
intervalo completo está intacto y lo que falta es resolución intermedia.

### PRJNA398590 — bastante más irregular

```
sitio               D0     D8    D14    D21    D28   total
penile skin         21     16      4      2     15      58
urine                8      7      0      0      5      20
vagina              20     18      4      3     17      62
```

D14 y D21 están prácticamente vacíos (4 y 2-3 muestras). Sólo **14 de 21
parejas** tienen D0+D28 en ambos miembros, y la orina es casi inutilizable (20
muestras en total, 2 sujetos con inicio y fin). Sirve como réplica del de 2021,
no como dataset principal.

### PRJNA1248021

```
sitio                A     A2      E     E2      H     H2      P     P2   total
Cervical            36      1     35      1     34      1     31      1     140
CorSul              36      1     35      1     36      1     32      1     143
```

Muy equilibrado por persona (4 puntos × 2 sitios), pero **qué intervalo
representa cada código A/E/H/P no está documentado en los metadatos**. El
script de los autores usa `time` 1–4 y `DayPostBln` (días desde basal), así que
el mapeo existe pero hay que sacarlo del paper o pedírselo a ellos.

---

## 3. Lo que NO está

Ahora sí es una auditoría válida: hay 140/418/284 registros y se han leído
todos. "Ausente" significa ausente.

| Variable | 398590 | 735440 | 1248021 |
|---|---|---|---|
| Carga bacteriana (qPCR) | ausente | ausente | ausente |
| **Nugent / Amsel / estado de VB** | **ausente** | **ausente** | **ausente** |
| Fase del ciclo menstrual | ausente | ausente | ausente |
| Circuncisión | ausente | ausente | **presente** (`host_phenotype`, 100%) |
| Brazo de tratamiento | ausente | ausente | n/a |
| Preservativo / actividad sexual | ausente | ausente | ausente |
| pH vaginal | ausente | ausente | ausente |
| Edad | ausente | ausente | presente (`host_age`) |
| Etnia | ausente | ausente | presente (`ethnicity`) |

El inventario completo por proyecto está en `docs/fase1_tablas.md` §4. Los
BioSample de los Plummer llevan poco más que el mínimo MIMARKS: fecha, sitio,
sexo, sujeto, dispositivo de recogida y coordenadas de Melbourne.

Dos matices importantes:

- **El brazo de tratamiento probablemente no falta, es que no hay.** Los dos
  son estudios de un solo brazo: "tratamiento concurrente a la pareja" y
  "aceptabilidad y tolerabilidad", no ensayos aleatorizados con control. Si
  todos los participantes recibieron tratamiento, no hay variable de brazo que
  depositar. Conviene confirmarlo en los papers antes de darlo por bueno.
- **Lo de Park es una pérdida real, no una ausencia estructural.** Su script de
  análisis maneja `std_qty` (carga), `NugentBV1`–`NugentBV4`, `Cyclephase`,
  `condom`, `pH_high` y `PID`. Todo eso existe en el fichero interno de los
  autores y nada de ello llegó al depósito, salvo circuncisión, edad y etnia.
  Es material que se consigue escribiéndoles, no regenerándolo.

---

## 4. Combinabilidad

| Par | Veredicto |
|---|---|
| 398590 ↔ 735440 | **Sí, a nivel de ASV** |
| 398590 ↔ 1248021 | No a nivel de ASV; sí a género |
| 735440 ↔ 1248021 | No a nivel de ASV; sí a género |

Los dos Plummer declaran los mismos cebadores byte a byte —
`CCTACGGGNGGCWGCAG` (341F) y `GACTACHVGGGTATCTAATCC` (805R)— y salen del mismo
grupo, la misma ciudad y el mismo tipo de diseño. Es el mejor escenario
posible para fusionar.

El de Park usa BactQuant (`CCTACGGGDGGCWGCA` / `GGACTACHVGGGTMTCTAATC`, del
repositorio de araclab) con espaciadores de heterogeneidad. Misma región
V3–V4, pero forward de 16 nt frente a 17 y reverse que empieza `GGACTAC`
frente a `GACTAC`: el amplicón tiene otros extremos, así que el mismo organismo
produce ASV distintos. Comparable sólo tras agregar a género.

**Tu criterio inicial era correcto pero incompleto:** región distinta descarta;
región igual no autoriza. Lo que decide son las secuencias exactas de los
cebadores.

Aun siendo ASV-combinables, los dos Plummer deben procesarse aprendiendo el
modelo de error por corrida de secuenciación (`learnErrors` separado) y
modelando el estudio como efecto de lote. Son de 2015 y 2018, con química
distinta muy probablemente.

---

## 5. Viabilidad de tus tres preguntas

### (a) Recolonización masculina en PRJNA735440 — **VIABLE**, es la pregunta bien planteada

El diseño da exactamente lo que hace falta: 27 hombres con serie completa en
piel peneana **y** en orina, de D0 a End (12 semanas), con D8, W4 y W8 en
medio. Los 43 controles (PBS, agua ultrapura, negativos de PCR, y dos
comunidades mock, Zymo y BEI) permiten descontaminar en serio, que a baja
biomasa como la piel peneana no es un lujo.

**La limitación, y es seria:** sin carga bacteriana, el 16S sólo da
proporciones. Si los anaerobios asociados a VB se desploman con el
tratamiento, la fracción de *Corynebacterium*, *Staphylococcus*,
*Anaerococcus* y *Finegoldia* sube aunque su número absoluto no se mueva.
"Reaparecen" y "nunca se fueron, sólo dejaron de estar diluidos" son
indistinguibles con estos datos. Se puede acotar razonando sobre ratios entre
taxones y sobre el orden temporal de reaparición —que es justo lo que
preguntas— pero no resolver. Hay que decirlo en los resultados, no esconderlo.

Un apunte a favor: la pregunta por el **orden** de reaparición es más robusta a
la composicionalidad que la pregunta por la magnitud. Si los anaerobios vuelven
en una secuencia reproducible entre hombres, eso es señal aunque las
proporciones estén distorsionadas.

### (b) Reaparición en él ↔ recurrencia en ella — **PARCIALMENTE VIABLE**

La mitad buena: el emparejamiento está completo, 27/27 parejas con ambos
miembros y con inicio y fin. La comparación entre miembros es posible.

La mitad mala: **no hay Nugent, ni Amsel, ni ningún estado de VB por punto
temporal.** "Recurrencia" es una definición clínica y no está en los datos. Lo
único que se puede hacer es un sustituto microbiológico —pérdida de dominancia
de *Lactobacillus*, cambio de tipo de comunidad, aparición de *Gardnerella*
y otros anaerobios— y **eso no es recurrencia de VB, es un proxy**. Hay que
llamarlo por su nombre en cualquier cosa que se escriba.

Tu sub-pregunta ("¿hay hombres que se repueblan y cuya pareja no recae?") sigue
siendo abordable, pero se convierte en "¿hay hombres que se repueblan y cuya
pareja no pierde la dominancia de *Lactobacillus*?". Es una pregunta distinta y
más débil. El dato clínico hay que pedírselo a los autores.

### (c) Efecto del coito en PRJNA1248021 — **VIABLE CON RESERVAS**

A favor: 4 puntos por persona, dos sitios, 73 personas, circuncisión al 100%,
edad y etnia. Suficiente para describir qué cambia en la piel peneana y qué
vuelve primero, que es lo que preguntas.

En contra, dos cosas:

1. **Los códigos A/E/H/P no están documentados.** Sin saber a qué intervalo
   corresponde cada uno, "cuánto dura" no tiene respuesta. Es lo primero que
   hay que resolver y no se resuelve con los datos: o el paper o los autores.
2. **Sin emparejamiento**, no se puede estudiar transferencia entre miembros,
   sólo la dinámica dentro de cada persona. Para tu pregunta (c) tal como está
   formulada basta; para conectarla con (a) y (b), no.

---

## 6. Recomendación para la Fase 2

Priorizar así:

1. **PRJNA735440 solo.** Es donde está tu pregunta principal y es el dataset
   mejor construido de los tres. Descargar, procesar con DADA2 y responder (a),
   con el proxy microbiológico para (b) claramente etiquetado como tal.
2. **Añadir PRJNA398590 como réplica**, no como ampliación de n. Mismos
   cebadores, así que se puede procesar junto, pero su cobertura es irregular
   (14/21 parejas con inicio y fin, orina casi vacía). Vale para comprobar si
   el patrón de (a) se repite en una cohorte independiente del mismo grupo.
3. **PRJNA1248021 aparte**, a nivel de género, y sólo después de aclarar el
   mapeo temporal.

Volumen a descargar en el paso 1: 375 muestras × ~36k pares de lecturas.
Del orden de 10–15 GB de FASTQ comprimido. Conviene traerlos de ENA por HTTPS
en vez de `sra-tools`.

**Antes de empezar la Fase 2 hay dos decisiones tuyas** (ver el mensaje que
acompaña a este documento).
