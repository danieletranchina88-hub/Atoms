# Atlante Orbitale

Un laboratorio di chimica quantistica che gira nel browser: atomi, orbitali, molecole, legami e reazioni calcolati **dai principi primi**, cioè risolvendo l'equazione di Schrödinger. Nessuna forma è disegnata a mano e nessun valore è inventato: ogni grandezza viene da un calcolo numerico, verificato contro programmi e dati di riferimento.

## Le cinque modalità

| Modalità | Cosa fa |
|---|---|
| **Atomo** | I 118 elementi con un calcolo autoconsistente DFT-LDA. Nuvola elettronica, gusci, configurazione (con le eccezioni alla regola di Madelung), caselle ↑↓ secondo Hund, energie degli orbitali, Z<sub>eff</sub> di Slater, energia di ionizzazione calcolata contro quella misurata. |
| **Orbitale** | Qualunque orbitale (n, l, m) fino a n = 7: p<sub>x</sub>, d<sub>z²</sub>, i sette f… Nuvola Monte Carlo, isosuperficie colorata per fase, nodi radiali e angolari, sezione piana di ψ. |
| **Legami σ π** | Orbitali molecolari LCAO σ, σ\*, π, π\* di molecole biatomiche e orbitali ibridi sp…sp³d², costruiti dagli orbitali atomici calcolati. |
| **Molecole** | Sandbox molecolare con Hartree–Fock *ab initio*. Si parte da una libreria di 67 molecole, da una stringa SMILES o modificando la molecola con un clic sugli atomi. Calcola: legami singoli, doppi e tripli (ordine di Lewis e ordine di Mayer quantistico), elettronegatività (Pauling, Allen), cariche formali, numeri di ossidazione, cariche parziali (Mulliken, Löwdin), VSEPR e ibridazione, risonanza, momento di dipolo, orbitali molecolari in 3D, densità elettronica, potenziale elettrostatico, geometria ottimizzata, vibrazioni con spettro IR, termodinamica statistica, curva di dissociazione. |
| **Reazioni** | Termochimica e cinetica: ΔE, ΔH, ΔS, ΔG e K<sub>p</sub> in funzione della temperatura, confronto con i ΔH° sperimentali. Stati di transizione, energie di attivazione e costanti di velocità (Eyring con correzione tunnel), grafico di Arrhenius, animazione della coordinata di reazione. |

## La fisica

### Atomi: campo centrale autoconsistente (DFT-LDA)
- Equazione radiale −½u″ + [l(l+1)/2r² + V(r)]u = εu, risolta con il metodo di Numerov su griglia logaritmica. Autovalori per shooting con il numero di nodi n − l − 1.
- V(r) = −Z/r + V<sub>H</sub>[ρ] + v<sub>xc</sub>[ρ]: scambio di Dirac e correlazione di Vosko–Wilk–Nusair. Si parte da Thomas–Fermi e si mescola con il metodo di Anderson (Herman–Skillman / Kohn–Sham).
- Parte angolare: armoniche sferiche reali. Nuvole campionate esattamente da |ψ|², isosuperfici con marching cubes.

### Molecole: Hartree–Fock *ab initio*
- Basi gaussiane contratte STO-3G, 3-21G, 6-31G, 6-31G\*, 6-31G\*\* dal [Basis Set Exchange](https://www.basissetexchange.org), elementi H–Kr.
- Integrali di sovrapposizione, cinetici, di attrazione nucleare, di dipolo e bielettronici con lo schema di McMurchie–Davidson e la funzione di Boys. Screening di Schwarz e simmetria a 8 permutazioni.
- Equazioni di Roothaan–Hall (RHF) e Pople–Nesbet (UHF) con ortogonalizzazione di Löwdin, stima iniziale di Hückel generalizzata e convergenza DIIS.
- Correlazione elettronica MP2 e UMP2 (Møller–Plesset al secondo ordine, nucleo congelato).
- Gradiente analitico dell'energia (derivate degli integrali), ottimizzazione quasi-Newton con Hessiana modello di Lindh, aggiornamento BFGS e passi RFO.
- Frequenze armoniche dall'Hessiana (differenze finite dei gradienti analitici), proiezione di traslazioni e rotazioni, intensità IR da ∂μ/∂Q.
- Stati di transizione con il partitioned-RFO di Baker e aggiornamento di Bofill.
- Termodinamica statistica (gas ideale, rotore rigido, oscillatore armonico): funzioni di partizione traslazionale, rotazionale e vibrazionale. Numero di simmetria dal gruppo puntuale, riconosciuto automaticamente.
- Analisi della funzione d'onda: popolazioni di Mulliken e Löwdin, ordini di legame di Mayer, ⟨S²⟩, potenziale elettrostatico esatto sulla superficie di densità 0,002 e/bohr³.

### Chimica "classica" dalla struttura
- Lettore e scrittore SMILES, con aromaticità e strutture di Kekulé.
- Ricerca delle strutture di Lewis (ottetti, cariche formali, ottetti espansi) e delle strutture di risonanza equivalenti.
- Numeri di ossidazione: gli elettroni di ogni legame vanno all'atomo più elettronegativo.
- VSEPR (AX<sub>n</sub>E<sub>m</sub>), ibridazione, carattere ionico di Pauling 1 − e<sup>−(Δχ)²/4</sup>.
- Geometria 3D di partenza: geometria delle distanze più repulsione dei domini elettronici (modello di Gillespie).

## Verifiche

`npm test` esegue tutte le verifiche numeriche:

| Verifica | Riferimento | Accordo |
|---|---|---|
| Energie atomiche LDA (H, He, Ne, Ar) | NIST | 10⁻⁶ Ha |
| Orbitali dell'idrogeno (energie, ⟨r⟩) | soluzione esatta | 10⁻⁷ |
| Energie HF (H₂O, NH₃, O₂ tripletto; STO-3G…6-31G\*\*) | PySCF 2.14 | 10⁻¹² Ha |
| Gradienti analitici HF e UHF | PySCF | 10⁻⁸ Ha/bohr |
| Cariche di Mulliken, dipolo, MP2, UMP2 | PySCF | 10⁻⁷ |
| Frequenze HF/6-31G\* dell'acqua (1827, 4071, 4189 cm⁻¹) | letteratura | 1 cm⁻¹ |
| Strutture di Lewis, risonanza, numeri di ossidazione, VSEPR, gruppi puntuali | casi da manuale | esatto |

`npm run test:all` calcola tutti i 118 atomi.

## Limiti

- Gli atomi sono trattati senza relatività e senza polarizzazione di spin: per l'oro l'energia di ionizzazione calcolata è 7,6 eV contro 9,2 eV misurati.
- Hartree–Fock non contiene la correlazione elettronica. MP2 ne recupera gran parte, ma con basi piccole restano errori di decine di kJ/mol sulle energie di reazione.
- Le molecole sono isolate, in fase gassosa (niente solvente). Le basi disponibili arrivano fino al kripton.
- Le frequenze sono armoniche e scalate con i fattori empirici di Scott e Radom.

## Avvio

Serve un server statico (i moduli ES e i Web Worker non funzionano aprendo il file direttamente):

```bash
npm start                                  # poi http://localhost:8080
npm test                                   # verifiche numeriche
node tools/precompute.mjs --parallel 4     # ricalcola la libreria di molecole (ore)
```

Funziona anche su GitHub Pages. Collegamenti diretti: `index.html#Fe`, `index.html#molecole`, `index.html#reazioni`.

## Struttura

```
index.html, css/style.css
js/physics/   atomi: elements, configuration (Madelung, Slater, Hund), scf (Numerov + LDA),
              harmonics, wavefunction, bonds (LCAO, ibridi), atomStore, worker
js/chem/      molecole: basisData, boys, integrals, hf, gradient, optimize, vibrations, symmetry,
              properties (Mulliken, Mayer, MP2), smiles, structure (Lewis, VSEPR, ossidazione),
              embed, elementData, library, libraryData (precalcolata), moleculeWorker, moleculeClient
js/render/    viewer (three.js), marching (marching cubes), moleculeView
js/ui/        periodicTable, charts, chemCharts, moleculeMode, reactionMode
tools/        precompute.mjs (libreria di molecole e stati di transizione)
tests/        verifiche numeriche
```
