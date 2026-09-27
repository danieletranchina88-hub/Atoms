# Atlante Orbitale

Simulazione quantistica degli atomi: le forme reali degli orbitali **s, p, d, f** di tutti i **118 elementi** della tavola periodica, calcolate risolvendo l'equazione di Schrödinger, più gli orbitali di legame **σ e π** e gli orbitali **ibridi**.

Non ci sono immagini disegnate a mano: ogni forma viene da un calcolo numerico fatto nel browser.

## Cosa si può fare

| Modalità | Contenuto |
|---|---|
| **Atomo** | Nuvola elettronica dell'atomo intero: ogni elettrone nel suo orbitale reale, riempito con la regola di Hund, oppure la densità totale con i gusci K, L, M… Tabella degli orbitali con caselle ↑↓, energie calcolate e Z<sub>eff</sub> di Slater. |
| **Orbitale** | Qualunque orbitale (n, l, m) fino a n = 7, occupato o eccitato: p<sub>x</sub>, d<sub>xy</sub>, d<sub>z²</sub>, i sette f… Nuvola Monte Carlo, isosuperficie colorata per fase, funzione radiale con i nodi, sezione piana di ψ. |
| **Legami σ π** | Orbitali molecolari LCAO (σ, σ*, π, π*) delle molecole biatomiche omonucleari alla distanza sperimentale, con diagramma degli OM, ordine di legame e magnetismo. Ibridi sp, sp², sp³, sp³d, sp³d² e la loro energia. |

## La fisica

1. **Campo centrale.** Ogni elettrone si muove nel potenziale del nucleo più il campo medio degli altri elettroni. La funzione d'onda si separa in ψ = R<sub>nl</sub>(r)·Y<sub>lm</sub>(θ,φ).
2. **Equazione radiale.** −½u″ + [l(l+1)/2r² + V(r)]u = εu, risolta con il metodo di Numerov su griglia logaritmica. L'autovalore si trova per shooting, imponendo n − l − 1 nodi radiali.
3. **Autoconsistenza (DFT-LDA).** V(r) = −Z/r + V<sub>Hartree</sub>[ρ] + v<sub>xc</sub>[ρ], con scambio di Dirac e correlazione di Vosko–Wilk–Nusair. Si itera partendo da Thomas–Fermi, con mescolamento di Anderson, fino alla convergenza (metodo di Herman–Skillman / Kohn–Sham).
4. **Parte angolare.** Armoniche sferiche reali, calcolate con i polinomi associati di Legendre.
5. **Visualizzazione.** Punti estratti esattamente da |ψ|² (inversione della distribuzione radiale e metodo del rigetto per quella angolare). Isosuperfici che racchiudono una frazione scelta della probabilità, estratte con marching cubes.
6. **Legami.** Combinazioni lineari degli orbitali atomici calcolati, con integrale di sovrapposizione S calcolato numericamente.

### Verifica

`npm test` confronta i risultati con i dati di riferimento del NIST (LDA non relativistica) e con le soluzioni esatte dell'idrogeno:

| Atomo | Energia totale calcolata (Ha) | NIST (Ha) |
|---|---|---|
| H | −0,445671 | −0,445671 |
| He | −2,834836 | −2,834836 |
| Ne | −128,233482 | −128,233481 |
| Ar | −525,946199 | −525,946195 |

`npm run test:all` esegue il calcolo per tutti i 118 elementi (tutti convergono, ciascuno in meno di 0,6 s) e confronta l'energia di ionizzazione ΔSCF con quella misurata.

### Limiti

Il calcolo è non relativistico, senza polarizzazione di spin e con occupazioni mediate sfericamente. Per gli elementi pesanti gli effetti relativistici diventano importanti (Au: 7,6 eV calcolati contro 9,2 eV misurati). Per Z > 103 le configurazioni sono previsioni teoriche.

## Avvio

Serve un semplice server statico (i moduli ES non funzionano aprendo il file direttamente):

```bash
npm start          # poi apri http://localhost:8080
```

Funziona anche con GitHub Pages. three.js viene caricato da CDN. Si può aprire direttamente un elemento con l'ancora, es. `index.html#Fe`.

## Struttura

```
index.html            pagina
css/style.css         stile (tema scuro e chiaro)
js/physics/           elements.js, configuration.js (Madelung, Slater, Hund),
                      scf.js (Numerov + DFT-LDA), harmonics.js, wavefunction.js,
                      bonds.js (LCAO, ibridi), atom.js, atomStore.js, worker.js
js/render/            viewer.js (three.js), marching.js (marching cubes)
js/ui/                periodicTable.js, charts.js
js/main.js            interfaccia
tests/                verifiche numeriche
```
