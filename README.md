# Atoms

Simulatore didattico nel browser focalizzato su due strumenti:

1. **Tavola periodica + visualizzazione di atomi e orbitali**  
   Tutti i 118 elementi con calcolo autoconsistente DFT-LDA relativistico scalare (Koelling–Harmon).  
   Densità |ψ|², orbitali individuali, sezioni, livelli energetici, configurazione elettronica.

2. **Sandbox di dinamica molecolare reattiva**  
   Atomi e molecole in una scatola. Motori: Lennard–Jones, GFN2-xTB (fino a 120 atomi, anche con solvente implicito ALPB), MINDO/3, Hartree–Fock.  
   Viste elettroniche dalla funzione d’onda, reazioni emergenti, microscopio di reazione, grafici di energia/temperatura/g(r).

Le altre modalità (molecole complete, reazioni, becher, cinetica, laboratorio, fasi, formazione del legame) sono state rimosse per concentrare l’interfaccia su questi due nuclei.


Un simulatore didattico di chimica nel browser. Calcola orbitali, molecole e grandezze termodinamiche con **modelli espliciti e approssimati**. L’accordo con un altro programma verifica l’implementazione dello stesso modello, non rende esatti i risultati fisici.

La vista dà priorità a esperimenti, misure e grafici. Le spiegazioni estese sono richiudibili; ogni modalità mostra il modello impiegato. Vedi [SCIENTIFIC_MODEL.md](SCIENTIFIC_MODEL.md) per dominio di validità, fonti e verifiche.

Novità: **viste degli elettroni** nella sandbox, tutte calcolate dalla funzione d'onda: *elettroni in movimento* (decine di migliaia di posizioni estratte con probabilità |ψ|² che seguono i nuclei), *coppie e ibridi* (orbitali localizzati di Pipek–Mezey: legami σ, π e doppietti solitari, ciascuno col suo colore, con l'ibridazione spⁿ di ogni atomo calcolata dalla composizione s/p/d: CH₄ sp³, C₂H₄ sp², C₂H₂ sp, benzene sp²), *densità di legame* (ρ meno gli atomi sferici) e *flusso degli elettroni* fra due istanti, più il diagramma cliccabile dei livelli orbitali. GFN2-xTB è da 2 a 3 volte più veloce (ortogonalizzazione di Cholesky, diagonalizzazione contigua in memoria, estrapolazione delle cariche) e converge anche nei sistemi quasi metallici. Prima: **GFN2-xTB** in JavaScript (Bannwarth, Ehlert, Grimme 2019), il tight binding quantistico di Grimme per tutti gli elementi da H a Rn, con **solvente implicito ALPB** (acqua): energie e forze uguali a quelle di tblite entro 10⁻¹⁰ hartree e 10⁻⁹ hartree/bohr. Nella sandbox è il motore predefinito fino a 120 atomi: con l'acqua implicita HCl si dissocia da solo in H₃O⁺ + Cl⁻ in un decimo di picosecondo e Na⁺ e Cl⁻ restano separati, mentre nel vuoto HCl resta molecolare e gli ioni si riattraggono. Prima ancora: atomi **relativistici** (equazione di Koelling–Harmon, verificata su tutti i 92 atomi delle tabelle NIST ScRLDA: per l'oro l'ionizzazione passa da 7,6 a 9,7 eV e il 6s si contrae del 14 %). Sandbox: **ioni** e specie cariche, **spin** fissato o libero, **campo elettrico** uniforme nella funzione d'onda MINDO/3, nessun trasferimento di elettroni fra molecole lontane; esperimenti di neutralizzazione H₃O⁺ + OH⁻ e di nitrazione del benzene (intermedio di Wheland → nitrobenzene). La vista predefinita è la densità di probabilità, non il modello a palline. L’atomo è |ψ|² degli orbitali calcolati (regola di Born); l’orbitale mostra il segno di ψ, che è una fase e non un colore osservato; la molecola, dopo Hartree–Fock, è un’isosuperficie di ρ e il legame è l’addensamento fra i nuclei. Il nucleo è fuori scala. Lewis, CPK e bastoncini restano una convenzione. Fedeltà automatica in sandbox: Lennard–Jones pubblicato per i nobili, GFN2-xTB fino a 120 atomi (MINDO/3 e Hartree–Fock a scelta), oltre il limite la modalità automatica rifiuta l’aggiunta; il campo classico qualitativo richiede una scelta esplicita.

## Le modalità

| Modalità | Cosa fa |
|---|---|
| **Atomo** | I 118 elementi con un calcolo autoconsistente DFT-LDA relativistico scalare (o non relativistico, a scelta), con la scheda dell'effetto relativistico. Vista predefinita: nuvola |ψ|² della densità totale. Si può scomporre negli orbitali occupati. Configurazione (con le eccezioni alla regola di Madelung), caselle ↑↓ secondo Hund, energie degli orbitali, Z<sub>eff</sub> di Slater, energia di ionizzazione calcolata contro quella misurata. |
| **Orbitale** | Qualunque orbitale (n, l, m) fino a n = 7: p<sub>x</sub>, d<sub>z²</sub>, i sette f… Nuvola Monte Carlo, isosuperficie colorata per fase, nodi radiali e angolari, sezione piana di ψ. |
| **Legami σ π** | Orbitali molecolari LCAO σ, σ\*, π, π\* di molecole biatomiche e orbitali ibridi sp…sp³d², costruiti dagli orbitali atomici calcolati. |
| **Come nasce un legame** | Due atomi si avvicinano e a ogni distanza l'equazione di Schrödinger viene risolta davvero: Hartree–Fock ristretto e non ristretto, interazione di configurazioni completa (soluzione esatta nella base) per H₂, GFN2-xTB. Per 12 molecole biatomiche (H₂⁺, H₂, He₂, LiH, Li₂, N₂, O₂, F₂, HF, CO, HCl, Cl₂) si vedono: l'orbitale molecolare che nasce dagli orbitali atomici, la densità di legame Δρ (dove gli elettroni si accumulano), la densità ρ, gli spin che si accoppiano (UHF); la curva di energia confrontata con le costanti spettroscopiche sperimentali (r<sub>e</sub>, D<sub>e</sub>, ω<sub>e</sub> di Huber e Herzberg, NIST); il diagramma degli orbitali con la composizione calcolata in 1s, 2s, 2p…; l'energia cinetica e potenziale (Ruedenberg) dal teorema del viriale; il punto di Coulson–Fischer; animazione "avvicina gli atomi". |
| **Molecole** | Sandbox molecolare con Hartree–Fock *ab initio*. Si parte da una libreria di 67 molecole, da una stringa SMILES o modificando la molecola con un clic sugli atomi. Calcola: legami singoli, doppi e tripli (ordine di Lewis e ordine di Mayer quantistico), elettronegatività (Pauling, Allen), cariche formali, numeri di ossidazione, cariche parziali (Mulliken, Löwdin), VSEPR e ibridazione, risonanza, momento di dipolo, orbitali molecolari in 3D (con carattere σ/π e legante/antilegante), densità elettronica, potenziale elettrostatico, geometria ottimizzata, vibrazioni con spettro IR, spettro UV-visibile (CIS), spettro fotoelettronico (Koopmans), termodinamica statistica, curva di dissociazione RHF/UHF. |
| **Reazioni** | Termochimica e cinetica: ΔE, ΔH, ΔS, ΔG e K<sub>p</sub> in funzione della temperatura, confronto con i ΔH° sperimentali. Stati di transizione, energie di attivazione e costanti di velocità (Eyring con correzione tunnel), grafico di Arrhenius, animazione della coordinata di reazione. |
| **Sandbox** | Dinamica molecolare in tempo reale, con reazioni e fasi visibili sulla scena: si mettono atomi e molecole in una scatola (tutti gli elementi da H a Rn con GFN2-xTB, 24 anche con il campo classico, 20 molecole pronte, qualunque SMILES) e si osservano cambi di connettività previsti dal modello. Si controllano temperatura (termostato o sistema isolato), volume, impulsi meccanici di energia hc/λ (non fotolisi), scintille; si possono afferrare gli atomi con una "pinzetta". Riconosce da sola specie e reazioni (registro delle reazioni, curve di concentrazione), misura la pressione sulle pareti e la distribuzione di Maxwell–Boltzmann, misura la deriva del bilancio energetico. Le forze possono venire da un calcolo quantistico a ogni passo: GFN2-xTB fino a 120 atomi, nel vuoto o in acqua implicita (ALPB), MINDO/3 semiempirico fino a 90 atomi (reazioni e spin dalla funzione d'onda) o Hartree–Fock *ab initio* fino a 8 atomi. Microscopio con replay della traiettoria, rilevamento dei singoli legami a ogni passo, ordini di legame, popolazioni elettroniche, geometria ed energia nel tempo. Viste ρ(r), Δρ, densità di spin e orbitali di frontiera selezionabili; urti controllati fra due reagenti da SMILES; g(r), diffusione, C<sub>V</sub>, barostato. |
| **Fasi · MD** | 256 atomi Lennard–Jones, frontiere periodiche, NVE/NVT, riscaldamento e raffreddamento. Energia, calore, lavoro, pressione, g(r), ordine e MSD. Modello di gas nobile; nessuna estensione arbitraria a sali o metalli. |
| **Becher** | Chimica in soluzione con reagenti liberi: si compone qualunque sale, acido o base da 21 cationi × 28 anioni (in soluzione a qualunque concentrazione o come solido), si aggiungono molecole (NH₃, etilendiammina, EDTA, KHP, indicatori), metalli e minerali, si scalda con una piastra termostatata. Per ogni aggiunta calcola l'equilibrio nel database disponibile: pH, precipitati (K<sub>sp</sub>), complessi, idrossidi anfoteri, reazioni redox con i metalli (Nernst), gas (H₂, CO₂), calore e temperatura, colore della soluzione dallo spettro, variazioni in mmol, curva di titolazione goccia a goccia. |
| **Cinetica** | Meccanismi di reazione scritti liberamente (costanti fisse, Arrhenius, forma k₀(T/300)ⁿ, tempi di dimezzamento, reazioni reversibili, ordini non stechiometrici, specie costanti), integrati con un metodo implicito per sistemi rigidi. Curve c(t), ordine apparente dalle leggi integrate, tempo di dimezzamento, grafico di Arrhenius, parametri di Eyring, velocità dei singoli passi, equilibri. Esperienze con costanti di letteratura: N₂O₅, ciclopropano, SN1 e SN2 (Atkins), serie radioattiva ²¹⁰Bi → ²¹⁰Po (NNDC), ciclo di Chapman dell'ozono (IUPAC), reazione oscillante di Belousov–Zhabotinsky (Oregonator). |
| **Laboratorio** | Le leggi della chimica fisica calcolate punto per punto: gas reali (van der Waals, costruzione di Maxwell, punto critico, fattore Z, Clausius–Clapeyron), acidi e basi (pH esatto dal bilancio di carica, titolazioni anche di acidi poliprotici, diagrammi di distribuzione, capacità tampone, indicatori), elettrochimica (pile, Nernst, ΔG = −nFE, K), cinetica (leggi integrate, reazioni consecutive, Arrhenius). |

## Seguire una reazione nel sandbox

Apri **Sandbox → Ab initio: H + H₂ oltre la barriera**, poi **Avvia**. Nel **Microscopio di reazione** vedi H1–H2 formarsi e H2–H3 rompersi, con ordini di legame, variazioni delle popolazioni elettroniche e bilancio energetico. Puoi fermarti automaticamente a un cambio di legame, spostare il cursore del replay o premere **Rivedi** accanto a un evento. Il replay non modifica la simulazione; le modifiche ai reagenti o alle condizioni tornano allo stato attuale.

Seleziona due, tre o quattro nuclei per misurare distanza, angolo e torsione. Attiva etichette, tracce dei nuclei o vettori di forza/velocità. **Fissa riferimento Δρ** confronta due densità elettroniche salvate; **densità di spin** mostra ρα−ρβ. Una sezione 2D della stessa griglia permette di guardare all’interno delle isosuperfici. L’elenco degli orbitali disponibili riporta energia, spin e occupazione.

**Urto controllato fra due reagenti** accetta molecole e ioni da SMILES: geometrie iniziali rilassate con GFN2-xTB, energia relativa, separazione dei centri di massa, parametro d’urto e rotazione di una molecola. Le velocità iniziali conservano la quantità di moto totale e riproducono l’energia richiesta; le forze determinano il risultato. Sono disponibili coppie iniziali per scambio H/H₂, incontro di radicali, associazione ionica, attacco nucleofilo, acido-base e coordinazione. Non sono sequenze di prodotti preprogrammate.

La registrazione salva fino a **600 campioni / 32 MiB**, rimuovendo i più vecchi. Automatico: un campione ogni passo fino a 8 atomi, ogni 5 fino a 40, ogni 20 sopra 40; ogni cambio di legame conserva comunque un campione se la registrazione è attiva. Le funzioni d’onda si archiviano fino a 40 atomi; sopra questo numero le viste elettroniche restano disponibili nello stato attuale. **Traiettoria JSON** esporta geometrie, velocità, forze, cariche, spin, legami, energie, condizioni e registro degli eventi, con unità dichiarate.

**Limite fisico:** queste sono dinamiche di Born–Oppenheimer e modelli approssimati, non una simulazione di ogni processo chimico. Δρ comprende anche il moto dei nuclei; Δn è una variazione di popolazione atomica. Nessuna delle due è una corrente elettronica. Non vengono simulate traiettorie di singoli elettroni, salti non adiabatici, stati eccitati o tunneling nucleare. Un massimo di energia nel tempo non è uno stato di transizione verificato.

## Rame, bicarbonato e cristalli

Apri **Becher → CuCl₂ + bicarbonato**, aggiungi NaHCO₃ e seleziona **Reticolo atomico · malachite**. La vista usa coordinate neutroniche sperimentali e permette rotazione, sezione e ispezione degli atomi. Acido, concentrazione e diluizione cambiano l’equilibrio. Puoi annullare le aggiunte ed esportare lo stato.

**Cu + HCl senza ossidante** è un controllo separato: non si impone la formazione di CuCl₂. Il carbonato basico è Cu₂CO₃(OH)₂, non un bicarbonato basico. La crescita del cristallo da soluzione non è simulata. La dinamica di fusione e raffreddamento in **Fasi · MD** opera entro il modello LJ.

Fonti, unità, test e limiti: [SCIENTIFIC_MODEL.md](SCIENTIFIC_MODEL.md).

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
- Potenziale classico didattico a ordine di legame, ispirato ad Abell–Tersoff–Brenner; non è una parametrizzazione REBO/ReaxFF validata. L'ordine di legame nasce dalla valenza libera degli atomi, con la conservazione della valenza σ e la riduzione dei π coniugati (benzene: 1,5). Lunghezza ed energia seguono la relazione di Pauling r(n) = r₁ − c ln n, D(n) = D₁nᵖ, con parametri ricavati dalle energie medie di legame tabulate.
- Cariche parziali per equalizzazione dell'elettronegatività (split-charge equilibration, parametri QEq/UFF); angoli dal teorema di Coulson (cos θ = −1/λ) con la correzione VSEPR per le coppie solitarie; van der Waals UFF; legame a idrogeno con il termine di DREIDING.
- Impulsi meccanici di energia E = hc/λ, conservando la quantità di moto; nessun modello di assorbimento; riconoscimento automatico di specie e reazioni con isteresi sull'ordine di legame.
- Dinamica *ab initio*: UHF a ogni passo con gradiente analitico, minimo energetico tra le due molteplicità più basse provate (non ricerca completa degli stati di spin), soluzione a spin rotto per la rottura dei legami.
- **Viste degli elettroni** (sandbox, *Vista microscopica*): gli **elettroni in movimento** sono punti estratti con probabilità ρ(r)·dV dalla densità di valenza calcolata (regola di Born) e agganciati al nucleo più vicino, così la nuvola segue la dinamica fra un campionamento e l'altro; le **coppie e gli ibridi** sono gli orbitali localizzati di Pipek–Mezey (rotazione degli orbitali occupati che massimizza le popolazioni di Mulliken atomiche, senza cambiare densità né energia), classificati in legami σ, π e doppietti solitari; per ogni centro l'ibrido spⁿ ha n = p/s della parte di orbitale su quell'atomo; la **densità di legame** è ρ − ρ⁰ con ρ⁰ gli atomi neutri sferici nella stessa base GFN2; il **flusso** è la variazione della densità di legame fra due istanti successivi (toglie gli atomi che si spostano con i nuclei e lascia la ridistribuzione fra legami e doppietti). Il pannello *Elettroni* elenca le coppie con gli ibridi e l'ibridazione di ogni atomo, oppure il diagramma dei livelli orbitali con le occupazioni; un clic isola un orbitale.
- **GFN2-xTB** (Bannwarth, Ehlert, Grimme, J. Chem. Theory Comput. 15, 1652, 2019), scritto da zero in JavaScript seguendo tblite: base STO-nG di valenza con armoniche sferiche, hamiltoniana di core dipendente dal numero di coordinazione, elettrostatica isotropa di secondo e terzo ordine per shell, elettrostatica anisotropa con dipoli e quadrupoli atomici (AES2), dispersione D4 autoconsistente (pesi dipendenti dalle cariche) più il termine a tre corpi di Axilrod–Teller–Muto, repulsione efficace, smearing di Fermi a 300 K con spin fissato, miscelamento di Broyden, campo elettrico uniforme e gradiente analitico completo (Hellmann–Feynman più Pulay). Parametri esportati da dxtb/tad-dftd4 (`tools/xtb/dump_gfn2.py`). Su 12 sistemi (neutri, ioni, O₂ tripletto, Zn, I, Fe quintetto, Pb, un cluster HCl·8H₂O) energie entro 10⁻¹⁰ hartree, cariche entro 10⁻⁸ e gradienti entro 10⁻⁹ hartree/bohr da tblite 0.7.0.
- **Solvente implicito ALPB** (Ehlert, Stahn, Spicher, Grimme, JCTC 17, 4250, 2021) per l'acqua: Born generalizzato con kernel P16, raggi di Born GBOBC e correzione ALPB di forma, area accessibile al solvente su griglia di Lebedev a 230 punti per le tensioni superficiali e i legami a idrogeno con il solvente, spostamento costante. Parametri presi dai sorgenti di tblite (`tools/xtb/dump_alpb.py`); energie e gradienti uguali a tblite entro 10⁻¹⁰ hartree e 10⁻⁹ hartree/bohr.
- **Dinamica quantistica MINDO/3** (Bingham, Dewar, Lo, JACS 97, 1285, 1975): a ogni passo una SCF UHF per tutti gli elettroni di valenza (orbitali di Slater STO-6G), forze dal gradiente analitico, fino a 90 atomi di H, B, C, N, O, F, P, S, Cl. Parametri e formule verificati contro l'implementazione di riferimento di PySCF (calori di formazione ed energie totali entro 10⁻⁶). Lo spin è libero (livello di Fermi comune a α e β, DIIS, smearing di Fermi a 300 K): O₂ esce tripletto, CH₃ doppietto, e due radicali con spin paralleli non si legano. Le reazioni non sono programmate: nella combustione dell'idrogeno compaiono da sole H₂ + O → OH + H, H + O₂ → HO₂, OH + H₂ → H₂O + H.
- **Come appare davvero la materia**: stile "nuvola elettronica" (densità ρ(r) dalla funzione d'onda: valenza MINDO/3 con la base ortogonalizzata alla Löwdin più il core degli atomi calcolato in DFT-LDA, oppure Hartree–Fock completo; con il campo classico la densità promolecolare degli atomi isolati) e stile "orbitali" (HOMO e LUMO calcolati, con il segno della funzione d'onda). Le griglie si calcolano in un Web Worker separato e si aggiornano durante la dinamica.
- Chimica fisica sulla traiettoria: funzione di distribuzione radiale g(r) normalizzata sulla distribuzione esatta delle distanze nel cubo (senza condizioni periodiche), spostamento quadratico medio con esponente α e coefficiente di diffusione dalla relazione di Einstein (solo nel regime diffusivo), capacità termica C<sub>V</sub> dalle fluttuazioni dell'energia nell'insieme canonico (con errore a blocchi e controllo di stazionarietà), barostato di Berendsen per lavorare a pressione costante.
- Esperimenti personali salvati nel browser (posizioni e velocità) ed esportazione delle coordinate in formato XYZ.

### Becher: equilibri in soluzione acquosa
- Database termodinamico **MINTEQ v4** (U.S. EPA, distribuito con [PHREEQC dell'USGS](https://www.usgs.gov/software/phreeqc-version-3)): 668 specie in soluzione e 74 solidi, ognuno con la sua fonte (NIST Critical Stability Constants, Bard 1985…). Il file è convertito da `tools/buildAqueousDB.mjs`. Aggiunte con fonte: ammino-complessi di Cu, Zn, Ni (database LLNL), FeSCN²⁺ (Inorg. Chim. Acta 2018), indicatori (Harris, *Quantitative Chemical Analysis*), entalpie di soluzione (CRC).
- Bilanci di massa per ogni componente, compresi H⁺ ed elettroni (pe, Eh): legge di azione di massa risolta con Newton–Raphson in ln c. Coefficienti di attività di Debye–Hückel esteso (parametri WATEQ) o di Davies, con A e B calcolati dalla costante dielettrica (Malmberg–Maryott) e dalla densità dell'acqua alla temperatura del becher.
- Dipendenza dalla temperatura con van 't Hoff dalle ΔrH del database: il pH neutro a 60 °C è 6,51.
- Precipitazione e dissoluzione con gli indici di saturazione (regola delle fasi); redox con i metalli fino ad annullare la f.e.m. di Nernst, senza sovratensioni empiriche (fugacità H₂ unitaria); CO₂ oltre la solubilità di Henry. L'Eh è mostrato solo quando c'è una coppia redox con entrambe le forme presenti (altrimenti un elettrodo non misurerebbe un potenziale stabile).
- Calore con la legge di Hess dalle ΔrH delle reazioni e ΔT = q/(m c<sub>p</sub>); colore dalla legge di Lambert–Beer e dalle funzioni colorimetriche CIE 1931.

### Come nasce un legame
- A ogni distanza R (37 punti, fitti vicino all'equilibrio) si risolvono RHF (Roothaan 1951) e UHF (Pople e Nesbet 1954) nella base scelta; ogni distanza parte dalla densità della distanza vicina, così si segue con continuità lo stesso stato. UHF si cerca in due passate (dagli atomi separati con spin opposti, e dall'equilibrio rompendo la simmetria HOMO/LUMO) e si tiene la soluzione più bassa.
- H₂: interazione di configurazioni completa nei due elettroni (Szabo e Ostlund, cap. 4), con occupazioni naturali e pesi delle configurazioni: a grande distanza è (σ<sub>g</sub>)² − (σ<sub>u</sub>)², la funzione di Heitler e London.
- Orbitali molecolari riscritti esattamente nella base degli orbitali degli atomi liberi (calcolati nella stessa base, medie sferiche delle densità di spin) con pesi di Mulliken; simmetria σ/π e g/u dalla funzione d'onda, carattere legante/antilegante dalla popolazione di sovrapposizione, ordine di legame di Mayer, integrali di sovrapposizione fra orbitali atomici.
- Energia cinetica e potenziale rispetto agli atomi dal teorema del viriale molecolare (Slater 1933), ΔT = −ΔE − R dE/dR e ΔV = 2ΔE + R dE/dR, sia sulla curva calcolata sia su quella sperimentale (Ruedenberg 1962; Bacskay e Nordholm 2013).
- Esperimento: r<sub>e</sub>, ω<sub>e</sub>, ω<sub>e</sub>x<sub>e</sub> di Huber e Herzberg (1979) dal NIST WebBook; D<sub>0</sub> da Huber e Herzberg, Liu et al. 2009 (H₂), Le Roy et al. 2009 (Li₂), di Lonardo e Douglas 1973 (HF), Le Roy 1973 (Cl₂), Aziz et al. 1995 (He₂); D<sub>e</sub> = D<sub>0</sub> + G(0). La curva sperimentale è il potenziale di Morse con queste costanti.

### Cinetica: meccanismi e sistemi rigidi
- Legge di azione di massa per ogni passo, r = k Π cᵢ^oᵢ, e k(T) = A (T/300)ⁿ exp(−E<sub>a</sub>/RT).
- Integratore di Rosenbrock ROS2 (Verwer, Spee, Blom e Hundsdorfer, SIAM J. Sci. Comput. 20, 1456, 1999): L-stabile, del secondo ordine, jacobiana analitica, passo adattivo.
- Analisi: leggi integrate di ordine 0, 1, 2 con R², tempo di dimezzamento, Arrhenius, Eyring (ΔH‡ = E<sub>a</sub> − RT), bilancio di elementi e carica di ogni passo.
- Fonti: Atkins, *Physical Chemistry*, tabelle 22.1 e 22.4; IUPAC (Atkinson et al., Atmos. Chem. Phys. 4, 1461, 2004); NNDC/ENSDF; Field e Noyes (1974) con le costanti di Field e Försterling (1986); atmosfera standard USA 1976.

## Verifiche

`npm test` esegue tutte le verifiche numeriche:

| Verifica | Riferimento | Accordo |
|---|---|---|
| Energie atomiche LDA (H, He, Ne, Ar) | NIST | tolleranza del test: 2·10⁻⁵ Ha |
| Orbitali dell'idrogeno (energie, ⟨r⟩) | soluzione esatta | 10⁻⁷ |
| Energie HF (H₂O, NH₃, O₂ tripletto; STO-3G…6-31G\*\*) | PySCF 2.14 | 10⁻¹² Ha |
| Gradienti analitici HF e UHF | PySCF | 10⁻⁸ Ha/bohr |
| Cariche di Mulliken, dipolo, MP2, UMP2, potenziale elettrostatico | PySCF | 10⁻⁷ |
| Energie di eccitazione CIS e forze dell'oscillatore | PySCF (TDA) | 10⁻³ eV |
| Frequenze HF/6-31G\* dell'acqua (1827, 4071, 4189 cm⁻¹) | letteratura | 1 cm⁻¹ |
| Strutture di Lewis, risonanza, numeri di ossidazione, VSEPR, gruppi puntuali | casi da manuale | esatto |
| Campo reattivo: forze = −∇E | differenze finite | 10⁻⁸ |
| Campo reattivo: geometrie (H₂O 105°, CH₄ 109,47°, C–C/C=C/C≡C) ed energie di reazione (2 H₂ + O₂: −434 contro −484 kJ/mol) | dati sperimentali | verifica qualitativa, scarto circa 10% per H₂/O₂ |
| Orbitali localizzati: numero di legami σ, π e doppietti e ibridazione di CH₄, NH₃, H₂O (sp³), C₂H₄, C₆H₆ (sp²), C₂H₂, CO₂ (sp); densità invariata e ortonormalità dopo la localizzazione; 8 elettroni di valenza campionati per H₂O; densità di legame a integrale nullo | chimica di riferimento; proprietà esatte | ibridi entro gli intervalli attesi; densità 10⁻¹⁵ |
| GFN2-xTB: energia, cariche, dipolo e gradiente di 12 sistemi (ioni, tripletto, Zn, I, Fe, Pb, HCl·8H₂O), nel vuoto, in campo elettrico e in acqua ALPB | tblite 0.7.0 | 10⁻⁹ hartree; gradienti 10⁻⁸ hartree/bohr (misurato ≤ 1,3·10⁻⁹) |
| MINDO/3: ΔfH di H₂O, CH₄, OH ed energie di H₂O singoletto e tripletto contro PySCF; forze contro differenze finite; 2 CH₃ → C₂H₆ (C–C 1,53 Å) | riferimento PySCF / sperimentale | 10⁻⁶ kcal/mol; forze 10⁻⁶ eV/Å |
| Campo reattivo: entalpie di atomizzazione di 23 molecole contro la legge di Hess sulle ΔfH° (CODATA; NBS in Atkins, tab. 2.5) | sperimentale | scarto medio 2,7 % su 20 molecole con legami ordinari (tutte entro il 7 %); N₂ −11 %, CO₂ −18 %, CO −37 % |
| Dinamica molecolare: ⟨r⟩ nel cubo = 0,6617 L (costante di Robbins), g(r) = 1 per il gas ideale, C<sub>V</sub> dell'argon gassoso = 3/2 Nk<sub>B</sub>, volume a pressione costante uguale a quello del gas ideale | analitico | entro l'errore statistico; volume entro il 3 % |
| Come nasce un legame: H₂ STO-3G a 1,4 bohr (RHF −1,1167, FCI −1,1373 Eh); D<sub>e</sub> sperimentali di H₂ e H₂⁺ contro le soluzioni esatte (Kołos–Wolniewicz, Bates et al.); FCI e UHF dissociano in due atomi H, RHF 6,3 eV sopra; H₂ FCI/6-31G\*\* r<sub>e</sub> 0,739 Å, D<sub>e</sub> 4,59 eV (esp. 0,741 Å, 4,75 eV); N₂ RHF più corto, più debole e più rigido; F₂ non legato in RHF; configurazione di N₂, polarità di LiH, dissociazione UHF di HCl; viriale | Szabo e Ostlund; teoria esatta; Huber e Herzberg | 2·10⁻⁴ Eh; 2 meV |
| Cinetica: soluzioni esatte (primo e secondo ordine, Bateman), problema rigido di Robertson (riferimento di Hairer e Wanner), equilibrio Q = k/k<sub>r</sub>, stato stazionario di Chapman, k del ciclopropano a 500 °C | analitico / Atkins | 10⁻⁵–10⁻⁷ relativo; k entro il 3 % |
| Becher: pH di acidi e basi deboli, tamponi, K<sub>sp</sub>, calore di neutralizzazione (ΔT = 6,67 K), pK<sub>w</sub>(T) (Bandura–Lvov 2006), Nernst per Fe³⁺/Fe²⁺, redox con i metalli; 587 ricette convergenti su 588, una redox fuori database rifiutata | valori da manuale | 0,01–0,05 unità di pH |

`npm run test:all` calcola tutti i 118 atomi.

## Limiti

- Gli atomi sono trattati senza relatività e senza polarizzazione di spin: per l'oro l'energia di ionizzazione calcolata è 7,6 eV contro 9,2 eV misurati.
- Hartree–Fock non contiene la correlazione elettronica. MP2 ne recupera gran parte, ma con basi piccole restano errori di decine di kJ/mol sulle energie di reazione.
- Le molecole sono isolate, in fase gassosa (niente solvente). Le basi disponibili arrivano fino al kripton.
- Le frequenze sono armoniche e scalate con i fattori empirici di Scott e Radom.
- Sandbox: il campo reattivo è classico (niente stati di spin, niente ipervalenza, barriere qualitative); la dinamica ab initio è limitata a pochi atomi e Hartree–Fock con base minima sovrastima le barriere (H + H₂: 1,0 eV contro 0,42 eV).
- GFN2-xTB è semiempirico: errori tipici di alcune kcal/mol sulle energie di reazione, barriere spesso sottostimate, base minima. ALPB rappresenta il solvente come continuo: niente struttura della seconda sfera di solvatazione né trasporto del protone lontano dallo ione. La dissoluzione di un cristallo (nanosecondi) è fuori portata di una dinamica quantistica nel browser.
- Sandbox quantistica: MINDO/3 è un metodo semiempirico (errori tipici sui calori di formazione di circa 11 kcal/mol, legami a idrogeno sottostimati); il cambio di spin fra superfici non è modellato (lo spin segue la funzione d'onda). Oltre 60–90 atomi il calcolo diventa lento.
- Sandbox: il campo reattivo sottostima le energie di atomizzazione di qualche percento (la regolarizzazione degli ordini di legame dà n ≈ 0,97 per un legame singolo) e molto di più per i legami multipli corti o con cariche formali (N₂, CO₂, CO).
- Come nasce un legame: nuclei fermi (si confronta D<sub>e</sub>, non D<sub>0</sub>), nessun effetto relativistico, basi finite (in He₂ un minimo spurio di pochi μeV per l'errore di sovrapposizione). Hartree–Fock sottostima D<sub>e</sub> di 1–5 eV e non lega F₂; UHF dissocia bene ma non è autostato dello spin; per O₂ un solo determinante con M<sub>S</sub> = 1 non descrive due atomi O(³P) e la curva UHF resta circa 1 eV sopra gli atomi. GFN2-xTB non è parametrizzato per la dissociazione (D<sub>e</sub> non affidabile). Il potenziale di Morse è esatto solo vicino al minimo e al limite di dissociazione.
- Cinetica: simulazioni isoterme a volume costante; le costanti di fotolisi del ciclo di Chapman sono ordini di grandezza; il quinto passo dell'Oregonator è fenomenologico.
- Becher: il modello di attività vale fino a I ≈ 0,5–1 M; le ΔrH mancanti nel database vengono trascurate (segnalato con n.d. nel registro); nessuna cinetica: ogni aggiunta arriva subito all'equilibrio (i solidi che si formano solo in tempi geologici sono esclusi, come l'ossidazione e la riduzione dell'acqua).

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
              mindo3 (dinamica quantistica semiempirica), densityWorker (nuvola elettronica e orbitali),
              sandboxData, sandboxWorker
js/chem/      (becher) aqueous (equilibri, colore), aqueousDB (generato da MINTEQ v4), aqueousExtra, beakerReagents
js/chem/      (cinetica) kinetics (meccanismi, integratore rigido, esperienze)
js/chem/      (legame) diatomicData (costanti sperimentali e fonti), bondFormation (scansioni, FCI, analisi), bondFormationWorker
js/ui/        sandboxMode, beakerMode, kineticsMode, bondFormationMode
tools/        precompute.mjs (libreria di molecole e stati di transizione), atomSummary.mjs, buildAqueousDB.mjs,
              build-artifact.py (versione a pagina singola)
tests/        verifiche numeriche
```