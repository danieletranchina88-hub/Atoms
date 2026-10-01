# Modelli e verifiche scientifiche

Questo progetto è un simulatore didattico, non un predittore universale della chimica. Un algoritmo numericamente corretto può avere un modello fisico inadeguato. L'interfaccia distingue dati misurati, calcoli e rappresentazioni illustrative.

## Atomi e orbitali

DFT-LDA in campo centrale sferico, senza relatività né polarizzazione di spin. Le configurazioni sono assegnate dal dataset. Le energie degli orbitali Kohn–Sham non sono tutte energie di ionizzazione sperimentali. Gli orbitali virtuali utilizzano una correzione di coda di Latter. Le combinazioni reali degli armonici sferici non sono in generale autostati di Lz.

La sonda radiale integra la densità normalizzata sulla griglia: P(r′ ≤ r) = ∫₀ʳ |u(r′)|² dr′. Interpolazione e campionamento hanno errore numerico. Le nuvole sono probabilità, non traiettorie elettroniche. Gli ibridi sp³d/sp³d² sono costruzioni illustrative e non dimostrano il coinvolgimento degli orbitali d nelle molecole ipervalenti.

Riferimento: [NIST SRD 141](https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations), energie **teoriche** LDA. La suite verifica H, He, Ne, Ar con tolleranza 2·10⁻⁵ Ha, non tutti i 118 elementi contro misure sperimentali.

## Molecole e reattore ideale

HF/UHF con basi gaussiane finite; MP2 opzionale. Nuclei classici, Born–Oppenheimer, molecole isolate. Correlazione, stati multiriferimento, spin e basi incomplete limitano i risultati. [NIST CCCBDB](https://cccbdb.nist.gov/) distingue confronti tra calcoli e confronti con esperimenti.

Termodinamica: gas ideale, rotore rigido, oscillatore armonico, stato standard p° = **101325 Pa (1 atm)**, coerente con `thermochemistry`. Le frequenze scalate sono una correzione empirica opzionale.

Il reattore è un bilancio per una sola reazione bilanciata a T,V fissi:

- nᵢ(ξ) = nᵢ₀ + νᵢ ξ; quantità iniziali dei reagenti modificabili e prodotti inizialmente assenti.
- pᵢ = nᵢRT/V; Q = ∏(pᵢ/p°)^νᵢ, quindi Q e K sono adimensionali.
- ΔᵣG = RT(ln Q − ln K); l'equilibrio è trovato per bisezione sul dominio nᵢ ≥ 0.
- Non è una simulazione del tempo: nessuna velocità è dedotta da ΔG. Reazioni concorrenti, condensazione e non idealità escluse.
- Per K estremi l'equilibrio può essere indistinguibile numericamente dal limite stechiometrico. Il programma dichiara questo limite.

## Elettrochimica

[IUPAC: equazione di Nernst](https://goldbook.iupac.org/terms/view/09068): E = E° − RT/(nF) Σνᵢ ln aᵢ. Ogni semireazione possiede stechiometria esplicita: include H⁺ per H₂, O₂ e MnO₄⁻, le specie ridotte solubili e le pressioni dei gas.

Attività adimensionali. Gas ideali: a ≈ p/(1 bar); solidi e liquidi puri: a = 1. L'utente può approssimare a ≈ c/(1 mol/L) solo nel limite ideale diluito. L'interfaccia mantiene 298,15 K perché non possiede E°(T). I potenziali del dataset sono arrotondati e non hanno precisione metrologica.

La reazione scritta ossida la semicella sinistra e riduce la destra. Il potenziale ha segno: se E < 0 il verso spontaneo è inverso. Lo schema inverte anodo, catodo e frecce; a E ≈ 0 non assegna corrente. Le semireazioni vengono moltiplicate per il minimo comune multiplo degli elettroni; **i potenziali non vengono moltiplicati**. Il diagramma è termodinamico: non garantisce realizzabilità pratica (metalli alcalini in acqua, passivazione, cinetica, giunzioni e reazioni concorrenti escluse). La curva E(pH) mantiene fisse le specie, anche dove la loro stabilità chimica reale cambierebbe.

## Sandbox

### UHF Born–Oppenheimer

Fino a 8 atomi e 40 funzioni di base, UHF/STO-3G, forze da gradiente analitico. All'inizio vengono confrontate le due molteplicità più basse consentite, poi si mantiene quella selezionata: non è una ricerca completa dello stato fondamentale. SCF non convergente o forze non finite arrestano il moto. Cambi di soluzione UHF possono introdurre discontinuità: la deriva energetica va controllata.

### Potenziale classico qualitativo

Il vecchio potenziale è mantenuto per esplorazione, con etichetta esplicita. È un assemblaggio specifico del progetto ispirato a forme Abell–Tersoff–Brenner, UFF, QEq e DREIDING, **non** un'implementazione parametrizzata e validata di REBO/ReaxFF. Include costanti empiriche, medie di legame e stime di Pauling. Non predice quantitativamente meccanismi, barriere, fasi o costanti cinetiche. Le cariche parziali e i frammenti neutri non simulano la redox in soluzione. La formula di un frammento non distingue isomeri; gli ordini frazionari non assegnano numeri di ossidazione o ibridazione.

### Integrazione e bilanci

Unità interne: Å, fs, u, eV, K. Velocity Verlet a passo fisso selezionabile. Un passo che prevede uno spostamento oltre 0,25 Å viene rifiutato prima di modificare la traiettoria: **non si tagliano velocità o energie per nascondere un'instabilità**. Il limite è un controllo numerico, non una legge fisica.

NVT: termostato CSVR di [Bussi, Donadio e Parrinello (2007)](https://doi.org/10.1063/1.2408420), inclusi segno del fattore e contabilizzazione del riscaldamento da energia cinetica nulla. NVE: termostato disabilitato. Con pareti esterne e nessuna rimozione del moto del centro di massa si usano 3N gradi di libertà per la temperatura cinetica. Per pochi atomi o urti diretti, T non implica equilibrio termico.

Energia totale: potenziale interatomico + cinetica + pareti + pinzetta. Spostare pareti o ancoraggio della pinzetta contabilizza il cambiamento di potenziale come lavoro esterno. Aggiunta/rimozione di atomi è uno scambio di energia con la materia. La diagnostica mostra la variazione di U − Q − W − E_materia dalla misura iniziale. Azzerare le misure cambia il riferimento; non corregge la traiettoria.

Il comando hc/λ è un **deposito meccanico**: risolve il cambiamento di velocità relativa includendo quella iniziale, aggiungendo esattamente l'energia richiesta e conservando la quantità di moto della coppia. Non calcola sezioni d'urto, regole di selezione, stati eccitati o rese fotochimiche. Il riscaldamento locale non è un modello di plasma o scintilla elettrica.

Pressione: media della forza sulle pareti morbide / 6L². Il volume geometrico differisce dal volume accessibile effettivo; fluttuazioni e transitori possono essere grandi. Le curve Maxwell–Boltzmann sono confronti di equilibrio, non un vincolo imposto ai dati.

## Verifica automatica

`npm test` comprende test originali (SCF, orbitali, HF, gradienti, proprietà, Lewis, potenziale, soluzione) e nuovi controlli:

- energia hc/λ e quantità di moto, per moto iniziale concorde e opposto;
- convergenza NVE al dimezzamento del passo, rifiuto senza mutazione di un passo instabile;
- lavoro di pareti e pinzetta, scambi di materia, termalizzazione e azzeramento;
- riproducibilità del seed; probabilità cumulativa 1s contro la soluzione coulombiana;
- Daniell, pila di concentrazione, inversione, dipendenza dal pH, Fe³⁺/Fe²⁺, alogeni, MnO₄⁻ e rifiuto di attività non positive;
- equilibrio analitico per H₂ + I₂ ⇌ 2HI con K fornito, bilancio degli elementi e risposta alla compressione per una dissociazione;
- protocollo Worker reale: avvio UHF in pausa, passi, svuotamento, aggiunta molecolare e censimento immediato.

Questi controlli verificano casi specifici e invarianti. Non certificano tutte le reazioni o tutti i parametri esistenti nel repository.

## Soluzioni, precipitati e reticoli atomici

Il becher calcola equilibri nel database MINTEQ v4, non traiettorie di reazione. Newton nei logaritmi delle concentrazioni, attività Debye–Hückel estese / Davies e insieme attivo delle fasi solide. Le costanti cambiano con T secondo van ’t Hoff dove ΔrH è disponibile. Il modello usa concentrazioni molari come approssimazione delle molalità nel limite diluito. Prima di accettare un'aggiunta si verificano bilanci delle componenti, consistenza della forza ionica e complementarità: SI ≤ 0 per le fasi assenti, SI = 0 per quelle presenti. Se un insieme attivo degenera si provano altri insiemi iniziali; nessun risultato non convergente viene accettato. Il precedente stato e il registro vengono ripristinati in caso di errore. Per I > 0,5 M compare un avviso: convergenza numerica non implica accuratezza fisica.

Le costanti provengono dal database **MINTEQ v4**, distribuito con PHREEQC, consultabile nel [mirror di PhreePlot](https://www.phreeplot.org/ppihtml/minteq.v4.dat.html). Le specie rame–cloruro e rame–carbonato riportano le costanti cumulative, non quelle successive. Le costanti dei solidi seguono le reazioni di dissoluzione sotto indicate; non vanno interpretate tutte come Ksp nel solo anione:

| Fase | Reazione che definisce K (acqua pura: attività 1) | log₁₀K a 25 °C |
|---|---|---:|
| Malachite | Cu₂CO₃(OH)₂ + 2 H⁺ ⇌ 2 Cu²⁺ + CO₃²⁻ + 2 H₂O | −5,306 |
| Azzurrite | Cu₃(CO₃)₂(OH)₂ + 2 H⁺ ⇌ 3 Cu²⁺ + 2 CO₃²⁻ + 2 H₂O | −16,906 |
| Atacamite | Cu₂(OH)₃Cl + 3 H⁺ ⇌ 2 Cu²⁺ + Cl⁻ + 3 H₂O | 7,391 |
| Smithsonite | ZnCO₃ ⇌ Zn²⁺ + CO₃²⁻ | −10 |
| Siderite | FeCO₃ ⇌ Fe²⁺ + CO₃²⁻ | −10,24 |

Il rame metallico non viene trasformato automaticamente in CuCl₂ da HCl diluito privo di ossidanti. La corrosione con O₂ atmosferico, la passivazione e i potenziali misti non sono implementati. Cu(I) e i suoi complessi sono inclusi: l’equilibrio formale può dare piccole quantità di Cu(I) stabilizzato dal cloruro, senza prevedere una velocità o formazione macroscopica di Cu(II). Il becher rifiuta esplicitamente Cu(II)/ioduro e metalli in nitrato acido, che richiedono specie redox mancanti; questa lista non certifica tutte le altre combinazioni. Il dataset è finito: mancano altri equilibri e fasi, inclusi precursori amorfi. Il risultato è l'equilibrio **entro le specie ammesse**, non una previsione universale dei prodotti reali.

La soglia empirica di 0,40 V per H₂ è stata rimossa: il calcolo redox usa Nernst con fugacità H₂ unitaria e nessuna previsione della velocità. La fuga di CO₂ conserva il criterio precedente di Henry con pCO₂ = 1 atm: rappresenta l'uscita di bolle, **non** un equilibrio con l'atmosfera terrestre né un recipiente sigillato con headspace finito. La stima calorimetrica adiabatica viene accoppiata all’equilibrio variando T, mentre la piastra impone T. Le ΔrH mancanti sono trascurate nel calcolo: il registro mostra n.d. quando mancano dati per le specie principali. Temperatura adiabatica e calore totale restano approssimazioni; usare la piastra per confronti a T nota. Le bande ottiche disponibili non coprono tutte le specie: i colori di alcune miscele sono incompleti.

Il registro mostra variazioni in mmol, senza costruire equazioni stechiometriche arrotondate e potenzialmente sbilanciate. Il JSON contiene composizione, totali, diagnostica e variazioni disponibili. Le bolle e i grani del becher sono illustrazioni. Il microscopio campiona le **specie** con il metodo dei maggiori resti (120 simboli, acqua esclusa); non pretende che tali simboli siano atomi individuali o che le posizioni siano traiettorie.

La malachite ha invece una vista atomica da coordinate sperimentali: [Zigan, Joswig, Schuster e Mason (1977), Z. Kristallogr. 145, 412–426](https://rruff.info/uploads/ZK145_421.pdf), **AMCSD 0010795**, coordinate depositate anche nel [catalogo Larixite](https://larixite.seescience.org/cifs/10795). Cella monoclinica a = 9,502 Å, b = 11,974 Å, c = 3,240 Å, β = 98,75°, P2₁/a, Z = 4. Le quattro operazioni di simmetria producono 40 atomi per cella: Cu₈C₄O₂₀H₈. Supercelle e sezione sono disponibili. I segmenti sono contatti geometrici (C–O < 1,55 Å, O–H < 1,15 Å, Cu–O < 2,55 Å), non ordini di legame calcolati. Il reticolo è statico e non simula nucleazione, crescita o forma del precipitato. Quando la fase non è presente nel becher, viene indicata come struttura di riferimento.

## Fasi · dinamica molecolare

Nuovo motore indipendente dal vecchio potenziale reattivo: Lennard–Jones monocomponente 12–6, con **forza traslata linearmente a rc = 2,5σ**, senza correzioni di coda. U_FS(r) = U_LJ(r) − U_LJ(rc) + (r − rc)F_LJ(rc), F_FS(r) = F_LJ(r) − F_LJ(rc); entrambi nulli oltre rc. Questa scelta è una specifica variante del potenziale, documentata nei [riferimenti di simulazione NIST](https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-properties); non si trasferiscono automaticamente i punti di transizione del LJ non troncato.

256 atomi, condizioni periodiche cubiche, convenzione della minima immagine con L > 2rc. Coordinate non avvolte per MSD; avvolte solo per la vista. Velocity Verlet NVE; splitting Langevin **BAOAB** NVT, propagazione esatta del passo Ornstein–Uhlenbeck e proiezione del moto del centro di massa, 3N−3 gradi di libertà. Riferimento per lo splitting: [Leimkuhler e Matthews, 2013](https://doi.org/10.1093/amrx/abs010). γ* = 1, Δt* = 0,002. Un passo instabile viene annullato, mai corretto tagliando le forze.

Il termostato registra ΔK come calore Q; il ridimensionamento affine della scatola registra ΔU come lavoro W. Diagnostica E−E₀−Q−W, pressione viriale P = (2K + Σrᵢⱼ·Fᵢⱼ)/(3V), g(r) con normalizzazione N(N−1), MSD e intensità di Bragg (200) rispetto all'FCC iniziale. Quest'ultima **non riconosce tutti i cristalli**, né distingue da sola vetro, liquido e gas. Nessuna etichetta di fase viene imposta da una soglia di temperatura. Il raffreddamento può produrre un liquido sovraraffreddato o un vetro nel tempo accessibile; non si forza la cristallizzazione.

Unità ridotte; mappatura **approssimata** all'argon con σ = 3,405 Å, ε/kB = 119,8 K, m = 39,948 u, τ ≈ 2,156 ps. Questi parametri sono comuni ma non un potenziale universale o di precisione: [Méndez-Bermúdez et al., Argon force field revisited](https://arxiv.org/abs/2201.08155). Non applicare questo modello a sali, acqua, metalli o reazioni chimiche. Il volume fissato può portare a pressione negativa in un cristallo in tensione; non è un errore di segno.

`tests/materials.mjs` verifica il gradiente del potenziale, continuità al cutoff, conservazione NVE, bilancio NVT e lavoro del volume, fusione emergente, simmetria/stechiometria/volume della cella sperimentale, precipitazione e dissoluzione della malachite, controllo Cu/HCl, rollback e 99 miscele diluite tra reagenti. La suite controlla implementazione e invarianti; non costituisce validazione sperimentale di tutti i prodotti e di tutte le condizioni.

La suite MINTEQ esamina inoltre le 588 ricette catione × anione: 587 convergono; rame/ioduro viene rifiutato esplicitamente perché manca la redox di I₂. Questi numeri descrivono copertura numerica, non certificazione sperimentale di tutte le ricette. Sono preservate le modalità MINDO/3, cinetica e analisi MD introdotte separatamente su main.
