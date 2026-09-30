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
