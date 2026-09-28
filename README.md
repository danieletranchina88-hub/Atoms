# Atlante Orbitale

Un laboratorio di chimica quantistica che gira nel browser: atomi, orbitali, molecole, legami e reazioni calcolati **dai principi primi**, cioè risolvendo l'equazione di Schrödinger. Nessuna forma è disegnata a mano e nessun valore è inventato: ogni grandezza viene da un calcolo numerico, verificato contro programmi e dati di riferimento.

## Le modalità

| Modalità | Cosa fa |
|---|---|
| **Atomo** | I 118 elementi con un calcolo autoconsistente DFT-LDA. Tavola periodica colorabile per elettronegatività, ionizzazione, affinità elettronica e raggio calcolato. Nuvola elettronica, gusci, configurazione (con le eccezioni alla regola di Madelung), caselle ↑↓ secondo Hund, energie degli orbitali, Z<sub>eff</sub> di Slater, energia di ionizzazione calcolata contro quella misurata. |
| **Orbitale** | Qualunque orbitale (n, l, m) fino a n = 7: p<sub>x</sub>, d<sub>z²</sub>, i sette f… Nuvola Monte Carlo, isosuperficie colorata per fase, nodi radiali e angolari, sezione piana di ψ. |
| **Legami σ π** | Orbitali molecolari LCAO σ, σ\*, π, π\* di molecole biatomiche e orbitali ibridi sp…sp³d², costruiti dagli orbitali atomici calcolati. |
| **Molecole** | Sandbox molecolare con Hartree–Fock *ab initio*. Si parte da una libreria di 67 molecole, da una stringa SMILES o modificando la molecola con un clic sugli atomi. Calcola: legami singoli, doppi e tripli (ordine di Lewis e ordine di Mayer quantistico), elettronegatività (Pauling, Allen), cariche formali, numeri di ossidazione, cariche parziali (Mulliken, Löwdin), VSEPR e ibridazione, risonanza, momento di dipolo, orbitali molecolari in 3D (con carattere σ/π e legante/antilegante), densità elettronica, potenziale elettrostatico, geometria ottimizzata, vibrazioni con spettro IR, spettro UV-visibile (CIS), spettro fotoelettronico (Koopmans), termodinamica statistica, curva di dissociazione RHF/UHF. |
| **Reazioni** | Termochimica e cinetica: ΔE, ΔH, ΔS, ΔG e K<sub>p</sub> in funzione della temperatura, confronto con i ΔH° sperimentali. Stati di transizione, energie di attivazione e costanti di velocità (Eyring con correzione tunnel), grafico di Arrhenius, animazione della coordinata di reazione. |
| **Sandbox** | Dinamica molecolare in tempo reale: si mettono atomi e molecole in una scatola (24 elementi, 20 molecole pronte, qualunque SMILES) e si osservano legami che si formano e si rompono. Si controllano temperatura (termostato o sistema isolato), volume, luce di lunghezza d'onda scelta (fotolisi), scintille; si possono afferrare gli atomi con una "pinzetta". Riconosce da sola specie e reazioni (registro delle reazioni, curve di concentrazione), misura la pressione sulle pareti e la distribuzione di Maxwell–Boltzmann, verifica il primo principio passo per passo. Per sistemi fino a 8 atomi le forze possono venire direttamente da Hartree–Fock (dinamica *ab initio*). |
| **Becher** | Chimica in soluzione con reagenti liberi: si compone qualunque sale, acido o base da 21 cationi × 28 anioni (in soluzione a qualunque concentrazione o come solido), si aggiungono molecole (NH₃, etilendiammina, EDTA, KHP, indicatori), metalli e minerali, si scalda con una piastra termostatata. Per ogni aggiunta calcola l'equilibrio completo: pH, precipitati (K<sub>sp</sub>), complessi, idrossidi anfoteri, reazioni redox con i metalli (Nernst), gas (H₂, CO₂), calore e temperatura, colore della soluzione dallo spettro, equazione ionica netta, curva di titolazione goccia a goccia. |
| **Cinetica** | Meccanismi di reazione scritti liberamente (costanti fisse, Arrhenius, forma k₀(T/300)ⁿ, tempi di dimezzamento, reazioni reversibili, ordini non stechiometrici, specie costanti), integrati con un metodo implicito per sistemi rigidi. Curve c(t), ordine apparente dalle leggi integrate, tempo di dimezzamento, grafico di Arrhenius, parametri di Eyring, velocità dei singoli passi, equilibri. Esperienze con costanti di letteratura: N₂O₅, ciclopropano, SN1 e SN2 (Atkins), serie radioattiva ²¹⁰Bi → ²¹⁰Po (NNDC), ciclo di Chapman dell'ozono (IUPAC), reazione oscillante di Belousov–Zhabotinsky (Oregonator). |
| **Laboratorio** | Le leggi della chimica fisica calcolate punto per punto: gas reali (van der Waals, costruzione di Maxwell, punto critico, fattore Z, Clausius–Clapeyron), acidi e basi (pH esatto dal bilancio di carica, titolazioni anche di acidi poliprotici, diagrammi di distribuzione, capacità tampone, indicatori), elettrochimica (pile, Nernst, ΔG = −nFE, K), cinetica (leggi integrate, reazioni consecutive, Arrhenius). |

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
- Stati di transizione con il partitioned-RFO di Baker e aggiornamento di Bofill. Ogni minimo della libreria è controllato: se compare una frequenza immaginaria, la struttura viene spostata lungo quel modo e riottimizzata.
- Stati eccitati di singoletto con CIS (Tamm–Dancoff), forze dell'oscillatore e spettro UV-visibile.
- Termodinamica statistica (gas ideale, rotore rigido, oscillatore armonico): funzioni di partizione traslazionale, rotazionale e vibrazionale. Numero di simmetria dal gruppo puntuale, riconosciuto automaticamente.
- Analisi della funzione d'onda: popolazioni di Mulliken e Löwdin, ordini di legame di Mayer, ⟨S²⟩, potenziale elettrostatico esatto sulla superficie di densità 0,002 e/bohr³.

### Chimica "classica" dalla struttura
- Lettore e scrittore SMILES, con aromaticità e strutture di Kekulé.
- Ricerca delle strutture di Lewis (ottetti, cariche formali, ottetti espansi) e delle strutture di risonanza equivalenti.
- Numeri di ossidazione: gli elettroni di ogni legame vanno all'atomo più elettronegativo.
- VSEPR (AX<sub>n</sub>E<sub>m</sub>), ibridazione, carattere ionico di Pauling 1 − e<sup>−(Δχ)²/4</sup>.
- Geometria 3D di partenza: geometria delle distanze più repulsione dei domini elettronici (modello di Gillespie).

### Sandbox: dinamica molecolare reattiva
- Equazioni di Newton integrate con velocity Verlet; termostato stocastico di Bussi–Donadio–Parrinello (insieme canonico) oppure sistema isolato (NVE); pressione come forza media sulle pareti.
- Campo di forze reattivo a ordine di legame (forma di Abell–Tersoff–Brenner). L'ordine di legame nasce dalla valenza libera degli atomi, con la conservazione della valenza σ e la riduzione dei π coniugati (benzene: 1,5). Lunghezza ed energia seguono la relazione di Pauling r(n) = r₁ − c ln n, D(n) = D₁nᵖ, con parametri ricavati dalle energie medie di legame tabulate.
- Cariche parziali per equalizzazione dell'elettronegatività (split-charge equilibration, parametri QEq/UFF); angoli dal teorema di Coulson (cos θ = −1/λ) con la correzione VSEPR per le coppie solitarie; van der Waals UFF; legame a idrogeno con il termine di DREIDING.
- Fotoni di energia E = hc/λ assorbiti dai legami; riconoscimento automatico di specie e reazioni con isteresi sull'ordine di legame.
- Dinamica *ab initio*: UHF a ogni passo con gradiente analitico, stato di spin più stabile (O₂ tripletto), soluzione a spin rotto per la rottura dei legami.
- Chimica fisica sulla traiettoria: funzione di distribuzione radiale g(r) normalizzata sulla distribuzione esatta delle distanze nel cubo (senza condizioni periodiche), spostamento quadratico medio con esponente α e coefficiente di diffusione dalla relazione di Einstein (solo nel regime diffusivo), capacità termica C<sub>V</sub> dalle fluttuazioni dell'energia nell'insieme canonico (con errore a blocchi e controllo di stazionarietà), barostato di Berendsen per lavorare a pressione costante.
- Esperimenti personali salvati nel browser (posizioni e velocità) ed esportazione delle coordinate in formato XYZ.

### Becher: equilibri in soluzione acquosa
- Database termodinamico **MINTEQ v4** (U.S. EPA, distribuito con [PHREEQC dell'USGS](https://www.usgs.gov/software/phreeqc-version-3)): 668 specie in soluzione e 74 solidi, ognuno con la sua fonte (NIST Critical Stability Constants, Bard 1985…). Il file è convertito da `tools/buildAqueousDB.mjs`. Aggiunte con fonte: ammino-complessi di Cu, Zn, Ni (database LLNL), FeSCN²⁺ (Inorg. Chim. Acta 2018), indicatori (Harris, *Quantitative Chemical Analysis*), entalpie di soluzione (CRC).
- Bilanci di massa per ogni componente, compresi H⁺ ed elettroni (pe, Eh): legge di azione di massa risolta con Newton–Raphson in ln c. Coefficienti di attività di Debye–Hückel esteso (parametri WATEQ) o di Davies, con A e B calcolati dalla costante dielettrica (Malmberg–Maryott) e dalla densità dell'acqua alla temperatura del becher.
- Dipendenza dalla temperatura con van 't Hoff dalle ΔrH del database: il pH neutro a 60 °C è 6,51.
- Precipitazione e dissoluzione con gli indici di saturazione (regola delle fasi); redox con i metalli fino ad annullare la f.e.m. di Nernst, con la sovratensione dell'idrogeno; CO₂ oltre la solubilità di Henry. L'Eh è mostrato solo quando c'è una coppia redox con entrambe le forme presenti (altrimenti un elettrodo non misurerebbe un potenziale stabile).
- Calore con la legge di Hess dalle ΔrH delle reazioni e ΔT = q/(m c<sub>p</sub>); colore dalla legge di Lambert–Beer e dalle funzioni colorimetriche CIE 1931.

### Cinetica: meccanismi e sistemi rigidi
- Legge di azione di massa per ogni passo, r = k Π cᵢ^oᵢ, e k(T) = A (T/300)ⁿ exp(−E<sub>a</sub>/RT).
- Integratore di Rosenbrock ROS2 (Verwer, Spee, Blom e Hundsdorfer, SIAM J. Sci. Comput. 20, 1456, 1999): L-stabile, del secondo ordine, jacobiana analitica, passo adattivo.
- Analisi: leggi integrate di ordine 0, 1, 2 con R², tempo di dimezzamento, Arrhenius, Eyring (ΔH‡ = E<sub>a</sub> − RT), bilancio di elementi e carica di ogni passo.
- Fonti: Atkins, *Physical Chemistry*, tabelle 22.1 e 22.4; IUPAC (Atkinson et al., Atmos. Chem. Phys. 4, 1461, 2004); NNDC/ENSDF; Field e Noyes (1974) con le costanti di Field e Försterling (1986); atmosfera standard USA 1976.

## Verifiche

`npm test` esegue tutte le verifiche numeriche:

| Verifica | Riferimento | Accordo |
|---|---|---|
| Energie atomiche LDA (H, He, Ne, Ar) | NIST | 10⁻⁶ Ha |
| Orbitali dell'idrogeno (energie, ⟨r⟩) | soluzione esatta | 10⁻⁷ |
| Energie HF (H₂O, NH₃, O₂ tripletto; STO-3G…6-31G\*\*) | PySCF 2.14 | 10⁻¹² Ha |
| Gradienti analitici HF e UHF | PySCF | 10⁻⁸ Ha/bohr |
| Cariche di Mulliken, dipolo, MP2, UMP2, potenziale elettrostatico | PySCF | 10⁻⁷ |
| Energie di eccitazione CIS e forze dell'oscillatore | PySCF (TDA) | 10⁻³ eV |
| Frequenze HF/6-31G\* dell'acqua (1827, 4071, 4189 cm⁻¹) | letteratura | 1 cm⁻¹ |
| Strutture di Lewis, risonanza, numeri di ossidazione, VSEPR, gruppi puntuali | casi da manuale | esatto |
| Campo reattivo: forze = −∇E | differenze finite | 10⁻⁸ |
| Campo reattivo: geometrie (H₂O 105°, CH₄ 109,47°, C–C/C=C/C≡C) ed energie di reazione (2 H₂ + O₂: −434 contro −484 kJ/mol) | dati sperimentali | pochi % |
| Campo reattivo: entalpie di atomizzazione di 23 molecole contro la legge di Hess sulle ΔfH° (CODATA; NBS in Atkins, tab. 2.5) | sperimentale | scarto medio 2,7 % su 20 molecole con legami ordinari (tutte entro il 7 %); N₂ −11 %, CO₂ −18 %, CO −37 % |
| Dinamica molecolare: ⟨r⟩ nel cubo = 0,6617 L (costante di Robbins), g(r) = 1 per il gas ideale, C<sub>V</sub> dell'argon gassoso = 3/2 Nk<sub>B</sub>, volume a pressione costante uguale a quello del gas ideale | analitico | entro l'errore statistico; volume entro il 3 % |
| Cinetica: soluzioni esatte (primo e secondo ordine, Bateman), problema rigido di Robertson (riferimento di Hairer e Wanner), equilibrio Q = k/k<sub>r</sub>, stato stazionario di Chapman, k del ciclopropano a 500 °C | analitico / Atkins | 10⁻⁵–10⁻⁷ relativo; k entro il 3 % |
| Becher: pH di acidi e basi deboli, tamponi, K<sub>sp</sub>, calore di neutralizzazione (ΔT = 6,67 K), pK<sub>w</sub>(T) (Bandura–Lvov 2006), Nernst per Fe³⁺/Fe²⁺, redox con i metalli; convergenza di tutte le 588 combinazioni catione × anione | valori da manuale | 0,01–0,05 unità di pH |

`npm run test:all` calcola tutti i 118 atomi.

## Limiti

- Gli atomi sono trattati senza relatività e senza polarizzazione di spin: per l'oro l'energia di ionizzazione calcolata è 7,6 eV contro 9,2 eV misurati.
- Hartree–Fock non contiene la correlazione elettronica. MP2 ne recupera gran parte, ma con basi piccole restano errori di decine di kJ/mol sulle energie di reazione.
- Le molecole sono isolate, in fase gassosa (niente solvente). Le basi disponibili arrivano fino al kripton.
- Le frequenze sono armoniche e scalate con i fattori empirici di Scott e Radom.
- Sandbox: il campo reattivo è classico (niente stati di spin, niente ipervalenza, barriere qualitative); la dinamica ab initio è limitata a pochi atomi e Hartree–Fock con base minima sovrastima le barriere (H + H₂: 1,0 eV contro 0,42 eV).
- Sandbox: il campo reattivo sottostima le energie di atomizzazione di qualche percento (la regolarizzazione degli ordini di legame dà n ≈ 0,97 per un legame singolo) e molto di più per i legami multipli corti o con cariche formali (N₂, CO₂, CO).
- Cinetica: simulazioni isoterme a volume costante; le costanti di fotolisi del ciclo di Chapman sono ordini di grandezza; il quinto passo dell'Oregonator è fenomenologico.
- Becher: il modello di attività vale fino a I ≈ 0,5–1 M; le ΔrH mancanti nel database vengono trascurate (segnalato con *); nessuna cinetica: ogni aggiunta arriva subito all'equilibrio (i solidi che si formano solo in tempi geologici sono esclusi, come l'ossidazione e la riduzione dell'acqua).

## Avvio

Serve un server statico (i moduli ES e i Web Worker non funzionano aprendo il file direttamente):

```bash
npm start                                  # poi http://localhost:8080
npm test                                   # verifiche numeriche
node tools/precompute.mjs --parallel 4     # ricalcola la libreria di molecole (ore)
```

Funziona anche su GitHub Pages. Collegamenti diretti: `index.html#Fe`, `index.html#molecole`, `index.html#reazioni`, `index.html#laboratorio`, `index.html#sandbox`, `index.html#becher`, `index.html#cinetica`.

## Struttura

```
index.html, css/style.css
js/physics/   atomi: elements, configuration (Madelung, Slater, Hund), scf (Numerov + LDA),
              harmonics, wavefunction, bonds (LCAO, ibridi), atomStore, worker
js/chem/      molecole: basisData, boys, integrals, hf, gradient, optimize, vibrations, symmetry,
              properties (Mulliken, Mayer, MP2), smiles, structure (Lewis, VSEPR, ossidazione),
              embed, elementData, library, libraryData (precalcolata), moleculeWorker, moleculeClient
js/render/    viewer (three.js), marching (marching cubes), moleculeView
js/ui/        periodicTable, charts, chemCharts, moleculeMode, reactionMode, labMode
js/chem/      (laboratorio) labData, labPhysics
js/chem/      (sandbox) reactive (campo di forze reattivo), reactiveData, md (dinamica), mdAnalysis (g(r), MSD, C_V), aimd (ab initio),
              sandboxData, sandboxWorker
js/chem/      (becher) aqueous (equilibri, colore), aqueousDB (generato da MINTEQ v4), aqueousExtra, beakerReagents
js/chem/      (cinetica) kinetics (meccanismi, integratore rigido, esperienze)
js/ui/        sandboxMode, beakerMode, kineticsMode
tools/        precompute.mjs (libreria di molecole e stati di transizione), atomSummary.mjs, buildAqueousDB.mjs,
              build-artifact.py (versione a pagina singola)
tests/        verifiche numeriche
```
