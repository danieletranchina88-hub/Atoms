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

## Come nasce un legame

Per ogni molecola biatomica la distanza R fra i nuclei viene variata da circa 2,4 r<sub>e</sub> (o r<sub>e</sub> + 3,4 Å) fino a 0,6 r<sub>e</sub>; a ogni R l'hamiltoniano elettronico non relativistico con nuclei fissi viene risolto nella base gaussiana scelta (STO-3G, 6-31G, 6-31G\*\*):

- **RHF** e **UHF** (Hartree–Fock ristretto e non ristretto, DIIS). Ogni distanza parte dalla densità della distanza vicina già risolta. Per UHF di singoletto si fanno due passate e si tiene a ogni R l'energia più bassa: una dagli atomi separati verso l'interno, con le densità UHF degli atomi liberi e lo spin dell'atomo B invertito (più una breve SCF con un termine che spinge α verso A e β verso B, che non entra nell'energia); l'altra dall'equilibrio verso l'esterno, ruotando l'HOMO α verso il LUMO e l'HOMO β in verso opposto. La distanza più grande a cui UHF coincide con RHF (⟨S²⟩ < 0,01) è il punto di Coulson–Fischer.
- **FCI per due elettroni** (H₂): Ψ = Σ c<sub>pq</sub> φ<sub>p</sub>(1) φ<sub>q</sub>(2) con c simmetrica negli orbitali RHF; (Hc)<sub>pq</sub> = Σ h<sub>pr</sub> c<sub>rq</sub> + Σ c<sub>ps</sub> h<sub>sq</sub> + Σ (pr|qs) c<sub>rs</sub>, diagonalizzata nello spazio delle coppie p ≥ q. Matrice densità γ = 2 c c, occupazioni naturali dai suoi autovalori.
- **GFN2-xTB** con atomi liberi calcolati con lo stesso metodo (anche H⁺).

Analisi della funzione d'onda:

- *Orbitali atomici di riferimento*: per ogni atomo libero (UHF ad alto spin nella stessa base) si costruisce la Fock degli elettroni α con le densità di spin mediate sfericamente, F = H + J[P̄<sub>α</sub> + P̄<sub>β</sub>] − K[P̄<sub>α</sub>], e la si diagonalizza separatamente nelle classi di parità (x, y, z) delle funzioni cartesiane: si ottengono 1s, 2s, 2p<sub>x</sub>, 2p<sub>y</sub>, 2p<sub>z</sub>… con ε ≈ −energia di ionizzazione (Koopmans). Insieme coprono esattamente la base dell'atomo.
- *Composizione*: un orbitale molecolare c viene riscritto in questa base, d<sub>a</sub> = u<sub>a</sub>ᵀ S<sub>AA</sub> c<sub>A</sub>, e il peso dell'orbitale atomico a è la popolazione di Mulliken d<sub>a</sub>(Uᵀ S c)<sub>a</sub>.
- *Simmetria*: σ, π, δ dal peso delle componenti cartesiane dispari in x e y; g/u dal segno di ⟨φ|î φ⟩ (solo se |⟨φ|î φ⟩| > 0,8: nelle soluzioni UHF a simmetria rotta gli orbitali sono localizzati).
- *Carattere*: popolazione di sovrapposizione 2 Σ<sub>μ∈A,ν∈B</sub> c<sub>μ</sub> c<sub>ν</sub> S<sub>μν</sub> (> 0,03 legante, < −0,03 antilegante). Ordine di legame di Mayer, cariche e spin di Mulliken.
- *Densità di legame*: Δρ = ρ − ρ<sub>pro</sub>, con ρ<sub>pro</sub> la somma delle densità sferiche degli atomi liberi alla stessa distanza.
- *Energia cinetica e potenziale*: teorema del viriale molecolare 2T + V = −R dE/dR (Slater 1933), che con E = T + V dà ΔT = −ΔE − R dE/dR e ΔV = 2ΔE + R dE/dR rispetto agli atomi (per i quali T = −E). dE/dR con differenze finite di Lagrange a tre punti sulla griglia non uniforme. Si usa la curva e non ⟨T⟩ della funzione d'onda perché, con gaussiane non riscalate, −⟨V⟩/⟨T⟩ differisce da 2 di qualche millesimo: un errore di alcuni eV sull'energia cinetica totale di N₂, più grande delle variazioni dovute al legame. Il rapporto della funzione d'onda è comunque mostrato come controllo.
- *Costanti dalla curva*: cubica per i quattro punti attorno al minimo → r<sub>e</sub>, curvatura k, ω<sub>e</sub> = √(k/μ)/(2πc) con le masse isotopiche; D<sub>e</sub> rispetto agli atomi liberi. Un minimo meno profondo di un decimo di D<sub>e</sub> sperimentale è indicato come "non legato".

Esperimento (js/chem/diatomicData.js): ω<sub>e</sub>, ω<sub>e</sub>x<sub>e</sub>, B<sub>e</sub>, r<sub>e</sub> dello stato X da Huber e Herzberg (1979) tramite il NIST Chemistry WebBook; D<sub>0</sub> dalle fonti citate per ogni molecola, D<sub>e</sub> = D<sub>0</sub> + ω<sub>e</sub>/2 − ω<sub>e</sub>x<sub>e</sub>/4. Verifica: per H₂ e H₂⁺ questo D<sub>e</sub> coincide entro 1 meV con i valori teorici esatti. La curva sperimentale è il potenziale di Morse V = D<sub>e</sub>(1 − e<sup>−β(R−r<sub>e</sub>)</sup>)² − D<sub>e</sub>, β = ω<sub>e</sub>√(μ/2D<sub>e</sub>).

## Elettrochimica

[IUPAC: equazione di Nernst](https://goldbook.iupac.org/terms/view/09068): E = E° − RT/(nF) Σνᵢ ln aᵢ. Ogni semireazione possiede stechiometria esplicita: include H⁺ per H₂, O₂ e MnO₄⁻, le specie ridotte solubili e le pressioni dei gas.

Attività adimensionali. Gas ideali: a ≈ p/(1 bar); solidi e liquidi puri: a = 1. L'utente può approssimare a ≈ c/(1 mol/L) solo nel limite ideale diluito. L'interfaccia mantiene 298,15 K perché non possiede E°(T). I potenziali del dataset sono arrotondati e non hanno precisione metrologica.

La reazione scritta ossida la semicella sinistra e riduce la destra. Il potenziale ha segno: se E < 0 il verso spontaneo è inverso. Lo schema inverte anodo, catodo e frecce; a E ≈ 0 non assegna corrente. Le semireazioni vengono moltiplicate per il minimo comune multiplo degli elettroni; **i potenziali non vengono moltiplicati**. Il diagramma è termodinamico: non garantisce realizzabilità pratica (metalli alcalini in acqua, passivazione, cinetica, giunzioni e reazioni concorrenti escluse). La curva E(pH) mantiene fisse le specie, anche dove la loro stabilità chimica reale cambierebbe.

## Sandbox

### UHF Born–Oppenheimer

Fino a 8 atomi e 40 funzioni di base, UHF/STO-3G, forze da gradiente analitico. All'inizio vengono confrontate le due molteplicità più basse consentite, poi si mantiene quella selezionata: non è una ricerca completa dello stato fondamentale. SCF non convergente o forze non finite arrestano il moto. Cambi di soluzione UHF possono introdurre discontinuità: la deriva energetica va controllata.

### Scelta automatica del modello

In modalità fedeltà automatica la scatola non usa un solo potenziale. He, Ne, Ar, Kr e Xe usano Lennard–Jones con σ e ε pubblicati, mixing di Lorentz–Berthelot e taglio a 2,5σ con correzione lineare della forza. La fase è letta dalla struttura locale (parametri q₆ di Steinhardt entro 1,5σ, criterio dei legami solidi di ten Wolde–Frenkel; vedi Fasi · dinamica molecolare), non da soglie su T* e ρ*: con le pareti e la coesistenza liquido–vapore la densità media della scatola non è quella delle fasi, e i punti triplo e critico del LJ completo (0,69 e 1,31) non valgono per il potenziale troncato. Fino a 40 atomi le forze vengono da GFN2-xTB, che copre tutti gli elementi da H a Rn; MINDO/3 (fino a 90 atomi di H, B, C, N, O, F, P, S, Cl) e Hartree–Fock/STO-3G (fino a 8 atomi) si scelgono a mano. Oltre, resta il campo classico qualitativo. Elementi senza parametri classici (oltre H–Ca, Br, Kr, I, Xe) richiedono GFN2-xTB. L’utente può bloccare un modello.

Il passo Δt è quello scelto (o del preset) e non viene mai aumentato automaticamente. Una versione precedente lo adattava alla forza massima: con atomi di idrogeno arrivava a 0,5–1 fs e la deriva dell'energia saliva a 10–40 eV in pochi ps (fino a 10⁵ eV nel preset Na + Cl₂), contro 0,003–0,8 eV con il passo fisso. Se un passo è instabile (spostamento > 0,25 Å) viene rifiutato prima di modificare la traiettoria, Δt si dimezza e il passo si ripete; la sandbox lo segnala. `tests/sandbox-worker.mjs` controlla che Δt non cresca e che l'acqua conservi l'energia.

### Ioni, carica totale, spin e campo elettrico

Specie cariche (dalla lista degli ioni o da SMILES come [NH4+], C[O-], O=[N+]=O) richiedono un motore quantistico: GFN2-xTB (predefinito), MINDO/3 o Hartree–Fock, che ricevono la carica totale della scatola (somma delle cariche formali delle specie inserite). Il campo classico ha frammenti neutri e rifiuta una carica netta. La carica mostrata per ogni molecola è la somma delle cariche atomiche calcolate, non un'etichetta: se un protone passa da H₃O⁺ a OH⁻, passa anche la carica. Calori di formazione di NH₄⁺, OH⁻, H₃O⁺, NO₂⁺, CH₃⁺ e H₂O⁺ uguali a quelli di PySCF (pyscf-semiempirical 0.1.1) entro 4·10⁻⁵ kcal/mol; le nostre forze coincidono con le differenze finite (5·10⁻⁷ eV/Å), quelle di PySCF per H₃O⁺ e NO₂⁺ no (errore 10⁻³–10⁻² eV/Å).

Molteplicità: libera (livello di Fermi comune ai due spin, lo spin esce dal calcolo) o fissata 2S+1; una molteplicità incompatibile con il numero di elettroni viene rifiutata.

Campo elettrico uniforme E (V/Å), con GFN2-xTB o MINDO/3. In GFN2-xTB il potenziale −E·r agisce sulle cariche e sui dipoli atomici (energia −E·μ, μ = Σ q_A R_A + d_A): energia e dipolo uguali a tblite entro 10⁻¹¹ hartree; il gradiente coincide con le differenze finite (10⁻¹⁰ hartree/bohr), mentre quello di tblite 0.7.0 con il campo non è la derivata della sua energia. In MINDO/3 l'hamiltoniana di core riceve +E·R_A sugli orbitali dell'atomo A e +E·D_A fra s e p dello stesso atomo, con D_A = ⟨ns|x|np_x⟩ = (2n+1)/√3·(4ζ_sζ_p)^(n+½)/(ζ_s+ζ_p)^(2n+2) per orbitali di Slater; i core contribuiscono −Z_A E·R_A. La forza su ogni atomo riceve q_A·E. Verificato: momento di dipolo uguale a −∂E/∂E entro 10⁻⁸ e·Å (acqua: 2,10 D; sperimentale 1,85 D), forze entro 2·10⁻⁷ eV/Å dalle differenze finite, forza netta su OH⁻ esattamente −E. Accendere o cambiare il campo, o fissare lo spin, è lavoro esterno contato nel bilancio dell'energia. Con base minima di valenza non c'è ionizzazione per effetto tunnel.

Conservazione degli elettroni per gruppi: la Fock fra atomi senza sovrapposizione (oltre 7 Å) è nulla, quindi un elettrone non può passare fra molecole lontane. Con un livello di Fermi unico, invece, l'aufbau lo spostava lo stesso: due acque a 15 Å in 2 V/Å diventavano H₂O⁻ e H₂O⁺, e CH₃⁻ con NO₂⁺ a 12 Å diventavano due radicali. Ora ogni gruppo di atomi accoppiati conserva il suo numero di elettroni (dalle popolazioni del passo precedente o, all'inizio, dalle cariche formali) con un proprio livello di Fermi; quando i gruppi si avvicinano e si accoppiano, il trasferimento torna possibile. Le forze restano il gradiente esatto (10⁻⁶ eV/Å). Una SCF non convergente viene ripresa dall'ultima densità convergente con smorzamento forte.


### GFN2-xTB e solvente implicito

Implementazione JavaScript di GFN2-xTB (Bannwarth, Ehlert, Grimme, JCTC 15, 1652, 2019) in `js/chem/xtb/`, validata contro tblite 0.7.0 (`tests/gfn2.mjs`, riferimenti generati da `tools/xtb/reference_gfn2.py`). Energia: E = Σ P·H⁰ + E_ES2 + E_ES3 + E_AES2 + E_D4(q) − TS + E_rep + E_ATM.
- Base STO-nG di valenza (armoniche sferiche reali, normalizzazione sull'autosovrapposizione); H⁰ fuori sito = ½(h_A + h_B)·S·π(R)·k_ll'·(2√(ζζ')/(ζ+ζ'))^½·(1 + 0,02 ΔEN²), con energie di sito h = h⁰ − k_CN·CN (CN a doppia esponenziale, raggi D3).
- Elettrostatica isotropa per shell: γ = (R² + η̄⁻²)^(−½) con η̄ media aritmetica; terzo ordine Σ Γ_l q³/3. Anisotropa (AES2): dipoli e quadrupoli atomici a traccia nulla dalla matrice densità, interazioni carica–dipolo, dipolo–dipolo e carica–quadrupolo smorzate con raggi dipendenti dal CN, più termini di nucleo.
- Dispersione D4 autoconsistente: C₆ dai riferimenti di Casimir–Polder pesati con gaussiane nel CN (erf, pesato con ΔEN) e con la funzione di carica ζ(q); il potenziale ∂E/∂q entra nella Fock. ATM con pesi a carica nulla.
- Occupazioni di Fermi a 300 K con livelli separati per α e β (spin fissato, orbitali comuni); l'energia è l'energia libera di Mermin, il gradiente usa W = Σ f ε c cᵀ.
- Costante eV→hartree 1/27,21138505 (quella di tblite per i livelli GFN2): con il valore CODATA 2018 l'energia si sposta di circa 6·10⁻⁸ hartree per atomo.

Solvente implicito ALPB per l'acqua (Ehlert et al., JCTC 17, 4250, 2021; ε = 80,2): E = ½ qᵀJq + Σ σ_A S_A + Σ h_A S_A q_A² + ΔG_shift, con J dal Born generalizzato P16 (raggi di Born GBOBC con parametri di descreening per elemento) più la correzione ALPB α/(ε·A_det) che dipende dalla forma della molecola, S_A l'area accessibile al solvente (sonda 1,13 Å, griglia di Lebedev a 230 punti, bordo smussato di 0,3 Å), σ le tensioni superficiali e h i legami a idrogeno con il solvente. Gradiente completo, comprese le derivate dei raggi di Born e delle aree. ΔG di solvatazione dell'acqua: −12,1 kcal/mol (sperimentale −6,3: è un modello parametrizzato).

Che cosa mostra, verificato in dinamica (300 K, termostato): HCl con quattro acque nel vuoto resta molecolare per 1 ps (Cl–H 1,3 Å, q(Cl) ≈ −0,4); con ALPB il protone passa all'acqua entro 0,15 ps e resta su H₃O⁺ (q(Cl) ≈ −0,9, Cl···H 1,8–2,0 Å). Na⁺ e Cl⁻ partiti a 6,4 Å con dieci acque: nel vuoto tornano a coppia ionica (2,3 Å) in 1,3 ps; con ALPB si allontanano (7,7 Å a 1,5 ps), con 4–6 ossigeni attorno al sodio. La coppia NaCl a contatto, anche in ALPB, resta unita per 3 ps: nell'acqua reale è un minimo metastabile che vive decine di picosecondi, e un cristallo si scioglie in tempi molto più lunghi.

Limite: la dinamica segue un solo stato elettronico (Born–Oppenheimer). Quando due stati quasi degeneri si scambiano (per esempio H₃O⁺ e OH⁻ a 3,6 Å, gap HOMO–LUMO da 2,5 a 4,6 eV) il bilancio mostra un salto (0,34 eV in quel caso). Esperimenti: neutralizzazione H₃O⁺ + OH⁻ → 2 H₂O; deprotonazione dell'intermedio di Wheland da parte di NH₃ (nitrobenzene + NH₄⁺ in 0,1 ps), con l'avvertenza che MINDO/3 lega il complesso σ benzene–NO₂⁺ di soli 3 kcal/mol; allineamento di sei molecole d'acqua in 0,5 V/Å (dipolo totale ≈ 12 D, quasi 6 × 2,1 D).

### Viste degli elettroni

Tutte calcolate dalla funzione d'onda corrente (GFN2-xTB salvo dove indicato) in un worker separato, su una griglia limitata alla zona delle molecole (margine 3,2 Å, passo 0,24 Å) e ricalcolate con un intervallo che si adatta al tempo di calcolo (mai più di una volta ogni 1,3 volte quel tempo).
- **Elettroni in movimento.** K punti (900 per elettrone di valenza, fra 9 000 e 45 000) estratti con probabilità ρ(r)·dV dalla densità di valenza (regola di Born); ogni punto è legato al nucleo più vicino e lo segue fra due campionamenti. Un punto non è un elettrone: dove i punti sono fitti è più probabile trovarlo. Verifica: la densità campionata integra a 8,000 elettroni di valenza per H₂O.
- **Coppie e ibridi.** Orbitali localizzati di Pipek–Mezey (J. Chem. Phys. 90, 4916, 1989) sugli orbitali doppiamente occupati: rotazioni di Jacobi a coppie che massimizzano Σ_i Σ_A (Q_A^ii)² con Q le popolazioni di Mulliken; a differenza di Foster–Boys non mescolano σ e π. Un orbitale è un doppietto se ≥ 85 % sta su un atomo, un legame a due centri se ≥ 80 % sta su due; σ o π secondo il carattere p lungo l'asse del legame o perpendicolare; per ogni centro l'ibrido è spⁿ con n = p/s delle popolazioni per momento angolare. L'ibridazione dell'atomo è la media pesata sui suoi legami σ e doppietti. Verifiche: CH₄, NH₃, H₂O sp³ (n = 3,0–3,2), C₂H₄ e benzene sp² (n = 1,8), C₂H₂ e CO₂ sp (n = 0,8–0,9); la densità 2Σcc^T non cambia (10⁻¹⁵) e gli orbitali restano ortonormali. Limiti: per l'acqua Pipek–Mezey dà un doppietto ricco di s e uno p puro (non due "orecchie da coniglio" equivalenti); nel benzene i tre π sono parzialmente delocalizzati e, con l'anello che vibra, possono passare da una struttura di Kekulé all'altra.
- **Densità di legame.** Δρ = ρ − ρ⁰ con ρ⁰ gli atomi neutri sferici nella stessa base (occupazioni di riferimento della shell divise fra le sue 2l+1 funzioni): integrale nullo, positiva dove gli elettroni si accumulano rispetto agli atomi isolati. Nei legami fra atomi molto elettronegativi (O–H, N≡N) il valore a metà legame rispetto ad atomi sferici è piccolo o negativo: è il comportamento noto delle mappe di deformazione (gli elettroni vanno nei doppietti e sull'atomo più elettronegativo); in H₂ e C–H è positivo.
- **Flusso degli elettroni.** Differenza fra le densità di legame in due istanti successivi della dinamica (con gli altri modelli: densità di valenza); soglia proporzionale all'intervallo, 7·10⁻⁴ e/bohr³ per fs, scelta perché su HCl in acqua a 300 K racchiuda circa l'1 % della griglia. Non è la corrente quantistica j(r) né la traiettoria di singoli elettroni.
- **Diagramma dei livelli.** Energie orbitali e occupazioni dell'ultimo calcolo (HOMO−8 … LUMO+5); un clic disegna l'orbitale canonico corrispondente.

### Prestazioni del calcolo GFN2-xTB

Ortogonalizzazione con la fattorizzazione di Cholesky S = LLᵀ (F' = L⁻¹FL⁻ᵀ, metà dei prodotti di matrice rispetto a Löwdin; Löwdin canonico se S è quasi singolare); diagonalizzazione di Householder + QL sulla matrice trasposta, con i cicli interni contigui in memoria e gli autovettori già per righe; punto di partenza dell'SCF estrapolato linearmente dagli ultimi due passi; tolleranza SCF di 10⁻⁵ in dinamica (la deriva NVE resta quella dell'integratore: 2,0 meV in 240 fs per otto acque, come con 10⁻⁶); ordini di legame di Mayer solo per le coppie entro 4,5 Å; termine ATM senza allocazioni. Risultato: un passo con 40 molecole d'acqua (120 atomi, 240 funzioni di base) passa da 1,6 a 0,6 s, 10 acque da 26 a 37 fotogrammi al secondo. Nei sistemi quasi metallici (Cu₂ + 4 Cl₂, gap 0,04 eV) gli elettroni oscillano fra gli atomi e il miscelamento di Broyden con smorzamento 0,4 non converge (nemmeno tblite in 250 cicli): se dopo 80 iterazioni non converge, l'SCF riparte dalla soluzione migliore con smorzamento 0,2 e poi 0,1, sempre a 300 K.

### Microscopio di reazione e traiettorie

`reactionTrace.js` osserva ogni **passo MD accettato** senza imporre reazioni. Il rilevatore di legami mantiene l’identità dei nuclei: passa a connesso per B_AB > 0,55 e torna disconnesso per B_AB ≤ 0,35, inclusi cambi interni a uno stesso frammento. Sono soglie operative di visualizzazione, non criteri fisici universali di legame. Si registrano tempo, coppia di nuclei, ordine prima/dopo, distanza, posizione e variazioni di popolazione. Confrontare Δt e Δt/2 resta necessario per eventi molto rapidi. Il comando manuale di avanzamento usa lo stesso controllo e rigetto del passo della dinamica continua.

La modalità automatica usa GFN2-xTB fino al limite di 120 atomi per la chimica e LJ pubblicato per soli gas nobili. Superato il limite, rifiuta l’aggiunta senza degradare automaticamente il modello. I vecchi esperimenti classici sono esplicitamente classici. Mescolare reagenti a un gas nobile seleziona il motore quantistico **prima** del calcolo delle forze sulla nuova composizione.

Il registratore conserva campioni reali, senza interpolazioni: coordinate, velocità, forze totali (interatomiche, pareti e pinzetta), cariche, popolazioni di spin, ordini di legame ed energia con il bilancio di calore/lavoro. I nuclei hanno indici stabili dentro un segmento; inventario, modello, carica totale, campo, solvente o spin cambiati iniziano un segmento nuovo. Buffer circolare: ≤600 campioni o 32 MiB; onde archiviate fino a 40 atomi. Il replay legge i campioni e le matrici SCF originali, senza reintegrare, ricalcolare gli elettroni o modificare le coordinate della dinamica. Le matrici vengono copiate all’archiviazione, evitando la sovrascrittura dal passo successivo. Il cursore non ripristina la simulazione: tornare alla vista attuale mostra il sistema rimasto in pausa. Il JSON esporta tutte le osservabili e le unità, senza matrici d’onda.

- **Popolazione atomica:** Δn_A = −[q_A(t)−q_A(t_ref)]. La somma è zero entro l’errore numerico a inventario/carica fissi. HF e GFN2 usano Mulliken; MINDO/3 popolazioni nella base ortogonalizzata del modello. Queste quantità dipendono dalla partizione e non assegnano un trasferimento da una coppia specifica di atomi o numeri di ossidazione. Nei grafici il riferimento è il primo campione ancora nel buffer.
- **Spin:** s_A = Σ_{μ∈A}(P^αS−P^βS)_μμ per HF/GFN2; MINDO/3 usa la base ortogonale. GFN2 condivide gli orbitali spaziali tra i due canali di Fermi: non è UHF, non introduce rottura della simmetria di spin o spin-flip dinamici. La modalità automatica GFN2 sceglie il numero minimo compatibile di elettroni spaiati; il tripletto O₂ richiede selezione esplicita. La griglia mostra ρ_α−ρ_β = Σ(P^α−P^β)_μν φ_μ φ_ν. Il core atomico appaiato non viene aggiunto allo spin. MINDO/3 trasforma anche P^α−P^β con S^−1/2, come la densità totale.
- **Differenza di densità:** Δρ(r) = ρ(r;t)−ρ(r;t_ref), con stesso inventario e metodo, su una griglia cartesiana comune e stessa risoluzione. Confronta **configurazioni diverse**: contiene sia ridistribuzione sia traslazione/rotazione dei nuclei e del core atomico. È volutamente diversa da una densità di deformazione a nuclei fissi. Non si deduce una corrente dalla sola variazione di popolazioni o densità. Le isosuperfici usano ±0,005 e/bohr³; la sezione mostra i valori della stessa griglia con scala cromatica logaritmica simmetrica. I tempi delle griglie sono riportati: il calcolo visivo asincrono può essere in ritardo rispetto ai nuclei nello stato attuale.
- **Orbitali:** lista finita di orbitali disponibili vicino ai livelli di frontiera, con energia, canale di spin e occupazione. Gli indici si riferiscono alla SCF corrente: non tracciano automaticamente un orbitale lungo crossing o degenerazioni. Il segno di ψ è una fase; non distingue carica positiva/negativa. Distanze in Å, angoli e torsioni in gradi; geometrie degeneri non ricevono angoli inventati. Tracce: posizioni realmente registrate, segmenti rettilinei di disegno tra campioni. Vettori: direzione di forze/velocità calcolate, riscalati per leggibilità; i moduli e le unità sono nelle etichette.

Perché questo non è un flusso elettronico: la dinamica segue una superficie elettronica adiabatica. Una corrente elettronica richiede ulteriori informazioni sulla funzione d’onda e la sua dipendenza temporale; non è ricostruita collegando arbitrariamente gli atomi che perdono e guadagnano popolazione. Riferimenti: [Bannwarth et al., GFN2-xTB (2019)](https://doi.org/10.1021/acs.jctc.8b01176); [Patchkovskii, Electronic currents and Born–Oppenheimer molecular dynamics (2012)](https://doi.org/10.1063/1.4747540).

### Urti controllati

`collision.js` prepara due reagenti con masse M_A e M_B, massa ridotta μ = M_AM_B/(M_A+M_B) e velocità relativa v = √[2E_rel/(μ·MV2)]. Nel riferimento del centro di massa: v_A = M_B v/(M_A+M_B), v_B = −M_A v/(M_A+M_B), quindi P_tot = 0 ed E_cin,trasl = E_rel. Separazione iniziale e parametro d’urto sono fra **centri di massa**, non distanze minime fra nuclei. Una rotazione rigida di B attorno a z cambia l’orientazione senza alterare la geometria interna. Sovrapposizioni inter-reagenti sotto 1,2 Å vengono rifiutate come controllo di preparazione, non come legge della chimica.

Le geometrie dei reagenti vengono rilassate con GFN2-xTB nel vuoto; la procedura è una minimizzazione locale con tolleranza di forza, non una ricerca globale dei conformeri. Il modello dinamico selezionabile è GFN2-xTB oppure UHF/STO-3G (≤8 atomi); per HF l’iniziale rilassata GFN2 non è necessariamente un minimo HF. Le velocità interne iniziali sono nulle; il sistema è NVE e il termostato parte spento. Non è un campionamento termico rovibrazionale, una sezione d’urto, una costante cinetica o un meccanismo garantito. Il solvente e lo spin possono poi essere cambiati usando i controlli, contabilizzando il lavoro.

Verifiche in `tests/reaction-trace.mjs`: energia traslazionale e quantità di moto per masse diverse, più energie e parametri d’urto; geometria rigida; normalizzazione della densità di spin per un H doppietto in tre metodi e spin nullo per H₂ singoletto; Δρ identica zero; scambio reale H1 + H2–H3 → H1–H2 + H3 con UHF, eventi ai singoli legami, conservazione dei campioni/SCF durante replay, buffer limitato, esportazione, passaggio LJ→GFN2 quando si aggiunge H₂ a 40 Ar, rifiuto di aggiunta oltre 120 atomi senza modificare l’inventario.

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


## Metadinamica Well-Tempered

Estensione della dinamica molecolare per l'esplorazione accelerata di reazioni chimiche con barriere energetiche elevate. Implementazione dell'algoritmo di Barducci, Bussi e Parrinello (Physical Review Letters 100, 020603, 2008).

**Variabili Collettive (CV)**: grandezze macroscopiche che descrivono il progresso della reazione. Il modulo `js/chem/metadynamics.js` implementa tre tipi di CV con gradienti analitici:

- **Distanza**: s = |R_i − R_j| per monitorare rottura/formazione di legami chimici
- **Angolo**: s = ∠(R_i, R_j, R_k) per cambiamenti conformazionali e riarrangiamenti
- **Numero di coordinazione**: s = Σ_j (1 − (r_ij/r₀)^n) / (1 − (r_ij/r₀)^m) con n = 6, m = 12 per reazioni complesse con cambiamento di coordinazione

**Potenziale di bias**: somma di gaussiane depositate nello spazio delle CV durante la simulazione:

V(s,t) = Σ_{t'<t} W₀ exp(−V(s(t'),t')/(kB ΔT)) exp(−Σ_k (s_k − s_k(t'))²/(2σ_k²))

dove W₀ è l'altezza iniziale della gaussiana, σ_k la larghezza lungo la k-esima CV, e ΔT = (γ − 1)T controlla la convergenza asintotica con γ = T_bias/T (bias factor).

**Forze di bias**: calcolate analiticamente come gradiente del potenziale di bias rispetto alle coordinate atomiche:

F^bias_i = −∇_{R_i} V(s,t) = −Σ_k (∂V/∂s_k)(∂s_k/∂R_i)

Le forze di bias vengono sommate vettorialmente alle forze del campo reattivo (`ReactiveFF`) a ogni passo di integrazione Velocity Verlet, accelerando termodinamicamente l'esplorazione dello spazio delle fasi.

**Superficie di Energia Libera (FES)**: a convergenza, il bias accumulato fornisce direttamente l'energia libera di Gibbs della reazione:

F(s) ≈ −(1 + T/ΔT) V(s) + costante = −(γ/(γ − 1)) V(s) + costante

Il metodo `reconstructFES(gridSize)` ricostruisce la FES su una griglia regolare nello spazio delle CV, permettendo di identificare stati di transizione (massimi locali) e intermedi di reazione (minimi locali).

**Parametri tipici**: W₀ = 0.001–0.01 Hartree, γ = 10–20, σ = 0.1–0.2 Å per distanze, 0.1–0.3 rad per angoli. Il deposition stride (tipicamente 50–100 step MD) controlla la frequenza di deposizione delle gaussiane.

**Test di verifica** (`tests/metadynamics.mjs`): gradienti analitici verificati contro differenze finite (tolleranza 1e-7), esplorazione di doppio pozzo armonico, ricostruzione FES con corretta topologia (minimi e barriere).

**Esempio applicativo** (`examples/metadynamics-water-dissociation.js`): dissociazione H₂O → OH + H con CV = distanza O-H, dimostrazione del superamento della barriera di attivazione e ricostruzione della FES.

**Limitazioni**: la convergenza della FES richiede simulazioni sufficientemente lunghe (tipicamente > 10000 step MD). Per sistemi con più di 2-3 CV, la ricostruzione della FES diventa computazionalmente proibitiva (maledizione della dimensionalità). Il metodo non fornisce direttamente le velocità di reazione, ma solo le energie libere.
