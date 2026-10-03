# Modelli e verifiche scientifiche

Questo progetto è un simulatore didattico, non un predittore universale della chimica. Un algoritmo numericamente corretto può avere un modello fisico inadeguato. L'interfaccia distingue dati misurati, calcoli e rappresentazioni illustrative.

## Atomi e orbitali

DFT-LDA in campo centrale sferico, senza polarizzazione di spin. Predefinito: **relativistico scalare**, lo schema ScRLDA del NIST: equazione radiale di Koelling–Harmon
u″ = [l(l+1)/r² + 2M(V−ε)]u + (M′/M)(u′ − u/r), M = 1 + (ε − V)/2c², che contiene la massa relativistica e il termine di Darwin ma non lo spin–orbita (p₁/₂ e p₃/₂ mediati); correzione relativistica dello scambio di MacDonald–Vosko (J. Phys. C 12, 2977, 1979); densità dalla sola grande componente, normalizzata. Numericamente: con u = √M·w il termine in u′ sparisce e si usa lo stesso Numerov; il termine −Z/r di V′ e V″ è derivato analiticamente; la griglia parte dentro la regione r < Z/2c², dove u ~ r^γ con γ = √(l(l+1) + 1 − Z²/c²). Confronto con le tabelle NIST ScRLDA di tutti i 92 elementi: |ΔE_tot| ≤ 5·10⁻⁴ Ha (≤ 4·10⁻⁷ relativo), autovalori entro 10⁻³ Ha (Sm escluso: la tabella NIST ha le colonne sfalsate). Effetto sull'oro: ionizzazione ΔSCF da 7,61 eV (non relativistica) a 9,71 eV (sperimentale 9,23 eV), 6s contratto del 14 %; sull'intera tavola l'errore medio della ionizzazione resta 0,5–0,8 eV per l'approssimazione LDA, non per la relatività. La scelta "non relativistica" riproduce il NIST LDA. Le configurazioni sono assegnate dal dataset. Le energie degli orbitali Kohn–Sham non sono tutte energie di ionizzazione sperimentali. Gli orbitali virtuali utilizzano una correzione di coda di Latter. Le combinazioni reali degli armonici sferici non sono in generale autostati di Lz.

La sonda radiale integra la densità normalizzata sulla griglia: P(r′ ≤ r) = ∫₀ʳ |u(r′)|² dr′. Interpolazione e campionamento hanno errore numerico. Le nuvole sono probabilità, non traiettorie elettroniche. Gli ibridi sp³d/sp³d² sono costruzioni illustrative e non dimostrano il coinvolgimento degli orbitali d nelle molecole ipervalenti.

Riferimento: [NIST SRD 141](https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations), energie **teoriche** LDA e ScRLDA. `tests/scf-reference.mjs` verifica H, He, Ne, Ar non relativistici con tolleranza 2·10⁻⁵ Ha; `tests/scrlda-reference.mjs` dodici atomi da H a U relativistici (tolleranza 10⁻³ Ha). Non sono confronti con misure sperimentali.

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

### Scelta automatica del modello

In modalità fedeltà automatica la scatola non usa un solo potenziale. He, Ne, Ar, Kr e Xe usano Lennard–Jones con σ e ε pubblicati, mixing di Lorentz–Berthelot e taglio a 2,5σ con correzione lineare della forza. La fase è letta dalla struttura locale (parametri q₆ di Steinhardt entro 1,5σ, criterio dei legami solidi di ten Wolde–Frenkel; vedi Fasi · dinamica molecolare), non da soglie su T* e ρ*: con le pareti e la coesistenza liquido–vapore la densità media della scatola non è quella delle fasi, e i punti triplo e critico del LJ completo (0,69 e 1,31) non valgono per il potenziale troncato. Fino a 36 atomi coperti da MINDO/3 le forze vengono dalla SCF semiempirica. Fino a 6 atomi trattabili, Hartree–Fock/STO-3G. Oltre, resta il campo classico qualitativo. L’utente può bloccare un modello.

Il passo Δt è quello scelto (o del preset) e non viene mai aumentato automaticamente. Una versione precedente lo adattava alla forza massima: con atomi di idrogeno arrivava a 0,5–1 fs e la deriva dell'energia saliva a 10–40 eV in pochi ps (fino a 10⁵ eV nel preset Na + Cl₂), contro 0,003–0,8 eV con il passo fisso. Se un passo è instabile (spostamento > 0,25 Å) viene rifiutato prima di modificare la traiettoria, Δt si dimezza e il passo si ripete; la sandbox lo segnala. `tests/sandbox-worker.mjs` controlla che Δt non cresca e che l'acqua conservi l'energia.

### Ioni, carica totale, spin e campo elettrico

Specie cariche (dalla lista degli ioni o da SMILES come [NH4+], C[O-], O=[N+]=O) richiedono un motore quantistico: MINDO/3 (fino a 90 atomi) o Hartree–Fock, che ricevono la carica totale della scatola (somma delle cariche formali delle specie inserite). Il campo classico ha frammenti neutri e rifiuta una carica netta. La carica mostrata per ogni molecola è la somma delle cariche atomiche calcolate, non un'etichetta: se un protone passa da H₃O⁺ a OH⁻, passa anche la carica. Calori di formazione di NH₄⁺, OH⁻, H₃O⁺, NO₂⁺, CH₃⁺ e H₂O⁺ uguali a quelli di PySCF (pyscf-semiempirical 0.1.1) entro 4·10⁻⁵ kcal/mol; le nostre forze coincidono con le differenze finite (5·10⁻⁷ eV/Å), quelle di PySCF per H₃O⁺ e NO₂⁺ no (errore 10⁻³–10⁻² eV/Å).

Molteplicità: libera (livello di Fermi comune ai due spin, lo spin esce dal calcolo) o fissata 2S+1; una molteplicità incompatibile con il numero di elettroni viene rifiutata.

Campo elettrico uniforme E (V/Å), solo MINDO/3: l'hamiltoniana di core riceve +E·R_A sugli orbitali dell'atomo A e +E·D_A fra s e p dello stesso atomo, con D_A = ⟨ns|x|np_x⟩ = (2n+1)/√3·(4ζ_sζ_p)^(n+½)/(ζ_s+ζ_p)^(2n+2) per orbitali di Slater; i core contribuiscono −Z_A E·R_A. La forza su ogni atomo riceve q_A·E. Verificato: momento di dipolo uguale a −∂E/∂E entro 10⁻⁸ e·Å (acqua: 2,10 D; sperimentale 1,85 D), forze entro 2·10⁻⁷ eV/Å dalle differenze finite, forza netta su OH⁻ esattamente −E. Accendere o cambiare il campo, o fissare lo spin, è lavoro esterno contato nel bilancio dell'energia. Con base minima di valenza non c'è ionizzazione per effetto tunnel.

Conservazione degli elettroni per gruppi: la Fock fra atomi senza sovrapposizione (oltre 7 Å) è nulla, quindi un elettrone non può passare fra molecole lontane. Con un livello di Fermi unico, invece, l'aufbau lo spostava lo stesso: due acque a 15 Å in 2 V/Å diventavano H₂O⁻ e H₂O⁺, e CH₃⁻ con NO₂⁺ a 12 Å diventavano due radicali. Ora ogni gruppo di atomi accoppiati conserva il suo numero di elettroni (dalle popolazioni del passo precedente o, all'inizio, dalle cariche formali) con un proprio livello di Fermi; quando i gruppi si avvicinano e si accoppiano, il trasferimento torna possibile. Le forze restano il gradiente esatto (10⁻⁶ eV/Å). Una SCF non convergente viene ripresa dall'ultima densità convergente con smorzamento forte.

Limite: la dinamica segue un solo stato elettronico (Born–Oppenheimer). Quando due stati quasi degeneri si scambiano (per esempio H₃O⁺ e OH⁻ a 3,6 Å, gap HOMO–LUMO da 2,5 a 4,6 eV) il bilancio mostra un salto (0,34 eV in quel caso). Esperimenti: neutralizzazione H₃O⁺ + OH⁻ → 2 H₂O; deprotonazione dell'intermedio di Wheland da parte di NH₃ (nitrobenzene + NH₄⁺ in 0,1 ps), con l'avvertenza che MINDO/3 lega il complesso σ benzene–NO₂⁺ di soli 3 kcal/mol; allineamento di sei molecole d'acqua in 0,5 V/Å (dipolo totale ≈ 12 D, quasi 6 × 2,1 D).

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

Il registro mostra un'equazione netta solo quando coefficienti interi piccoli (≤ 12, scarto ≤ 4 % dalle variazioni calcolate) la bilanciano **esattamente** in tutte le componenti del database, nelle cariche e nell'acqua; altrimenti mostra le variazioni in mmol. H₂CO₃* del database è scritto come CO₂(aq), con l'acqua corrispondente. Esempio: 2 Cu²⁺ + 4 HCO₃⁻ → Cu₂(OH)₂CO₃(s) + 3 CO₂(aq) + H₂O. Il JSON contiene composizione, totali, diagnostica e variazioni disponibili. Le bolle e i grani del becher sono illustrazioni. Il microscopio campiona le **specie** con il metodo dei maggiori resti (120 simboli, acqua esclusa); non pretende che tali simboli siano atomi individuali o che le posizioni siano traiettorie.

La malachite ha invece una vista atomica da coordinate sperimentali: [Zigan, Joswig, Schuster e Mason (1977), Z. Kristallogr. 145, 412–426](https://rruff.info/uploads/ZK145_421.pdf), **AMCSD 0010795**, coordinate depositate anche nel [catalogo Larixite](https://larixite.seescience.org/cifs/10795). Cella monoclinica a = 9,502 Å, b = 11,974 Å, c = 3,240 Å, β = 98,75°, P2₁/a, Z = 4. Le quattro operazioni di simmetria producono 40 atomi per cella: Cu₈C₄O₂₀H₈. Supercelle e sezione sono disponibili. I segmenti sono contatti geometrici (C–O < 1,55 Å, O–H < 1,15 Å, Cu–O < 2,55 Å), non ordini di legame calcolati. Il reticolo è statico e non simula nucleazione, crescita o forma del precipitato. Quando la fase non è presente nel becher, viene indicata come struttura di riferimento.

## Fasi · dinamica molecolare

Nuovo motore indipendente dal vecchio potenziale reattivo: Lennard–Jones monocomponente 12–6, con **forza traslata linearmente a rc = 2,5σ**, senza correzioni di coda. U_FS(r) = U_LJ(r) − U_LJ(rc) + (r − rc)F_LJ(rc), F_FS(r) = F_LJ(r) − F_LJ(rc); entrambi nulli oltre rc. Questa scelta è una specifica variante del potenziale, documentata nei [riferimenti di simulazione NIST](https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-properties); non si trasferiscono automaticamente i punti di transizione del LJ non troncato.

256 atomi, condizioni periodiche cubiche, convenzione della minima immagine con L > 2rc. Coordinate non avvolte per MSD; avvolte solo per la vista. Velocity Verlet NVE; splitting Langevin **BAOAB** NVT, propagazione esatta del passo Ornstein–Uhlenbeck e proiezione del moto del centro di massa, 3N−3 gradi di libertà. Riferimento per lo splitting: [Leimkuhler e Matthews, 2013](https://doi.org/10.1093/amrx/abs010). γ* = 1, Δt* = 0,002. Un passo instabile viene annullato, mai corretto tagliando le forze.

Il termostato registra ΔK come calore Q; il ridimensionamento affine della scatola registra ΔU come lavoro W. Diagnostica E−E₀−Q−W, pressione viriale P = (2K + Σrᵢⱼ·Fᵢⱼ)/(3V), g(r) con normalizzazione N(N−1), MSD e intensità di Bragg (200) rispetto all'FCC iniziale. Quest'ultima **non riconosce un cristallo con orientazione diversa**, per esempio uno nato dal liquido. Per questo la fase è letta dalla struttura locale: per ogni atomo i parametri di Steinhardt q₆ₘ sui vicini entro 1,5σ (primo minimo di g(r)); due vicini con prodotto scalare normalizzato dei q₆ > 0,7 formano un legame solido e un atomo con almeno 7 legami solidi è cristallino (ten Wolde, Ruiz-Montero, Frenkel, J. Chem. Phys. 104, 9932, 1996); con meno di 3 vicini è vapore. Per separare un liquido da un solido disordinato (vetro, aggregato non cristallino) si misura la mobilità senza risentire di traslazioni e rotazioni d'insieme: la frazione di coppie di primi vicini (r < 1,35σ) che in 3τ si separano oltre 1,65σ. Misurata: 0,5–0,6 nell'argon liquido a 90 K, 0,29 nel liquido di bulk a T* = 0,75, 0,05–0,12 nell'aggregato di 160 atomi a 40 K (si muove solo la superficie), 0 nel cristallo; soglia 0,2. Il criterio q₆ è invariante per rotazione; Q₆ globale vale 0,5745 per l'FCC perfetto (verificato) e ≈ 0,03 nel liquido. Non distingue FCC, HCP e BCC fra loro e, sopra il punto critico, "liquido + vapore" indica solo un fluido non uniforme. Nessuna etichetta di fase viene imposta da una soglia di temperatura. Il raffreddamento può produrre un liquido sovraraffreddato o un vetro nel tempo accessibile; non si forza la cristallizzazione.

Unità ridotte; mappatura **approssimata** all'argon con σ = 3,405 Å, ε/kB = 119,8 K, m = 39,948 u, τ ≈ 2,156 ps. Questi parametri sono comuni ma non un potenziale universale o di precisione: [Méndez-Bermúdez et al., Argon force field revisited](https://arxiv.org/abs/2201.08155). Non applicare questo modello a sali, acqua, metalli o reazioni chimiche. Il volume fissato può portare a pressione negativa in un cristallo in tensione; non è un errore di segno.

`tests/structure.mjs` verifica Q₆ dell'FCC, invarianza per rotazione, cristallo, liquido e gas, e le equazioni nette bilanciate del becher. `tests/materials.mjs` verifica il gradiente del potenziale, continuità al cutoff, conservazione NVE, bilancio NVT e lavoro del volume, fusione emergente, simmetria/stechiometria/volume della cella sperimentale, precipitazione e dissoluzione della malachite, controllo Cu/HCl, rollback e 99 miscele diluite tra reagenti. La suite controlla implementazione e invarianti; non costituisce validazione sperimentale di tutti i prodotti e di tutte le condizioni.

La suite MINTEQ esamina inoltre le 588 ricette catione × anione: 587 convergono; rame/ioduro viene rifiutato esplicitamente perché manca la redox di I₂. Questi numeri descrivono copertura numerica, non certificazione sperimentale di tutte le ricette. Sono preservate le modalità MINDO/3, cinetica e analisi MD introdotte separatamente su main.
