// Argomento 1 — Materia: proprietà e misura.
// Riferimento: R. H. Petrucci, F. G. Herring, J. D. Madura, C. Bissonnette, Chimica generale. Principi ed applicazioni
// moderne, Piccin, cap. 1. Definizioni delle unità: BIPM, Il Sistema Internazionale di Unità (SI), 9ª ed. (2019).
// Densità: CRC Handbook of Chemistry and Physics.

import { parseNumber, formatSig, formatPlain } from '../numeri.js';
import { sorter, quiz, practice, steppers, checklist, rand } from '../ui.js';
import { UNITA, NOTE_UNITA, convert, convertTemperature, SCALE_T, GENERATORI, digitsExplanation } from './materia-calcoli.js';

export const SEZIONI = [
  ['obiettivi', 'Obiettivi'],
  ['metodo', 'Metodo scientifico'],
  ['proprieta', 'Proprietà'],
  ['stati', 'Stati della materia'],
  ['classificazione', 'Classificazione'],
  ['si', 'Unità SI'],
  ['temperatura', 'Temperatura'],
  ['densita', 'Densità e %'],
  ['incertezza', 'Incertezza'],
  ['cifre', 'Cifre significative'],
  ['conversioni', 'Analisi dimensionale'],
  ['svolti', 'Esercizi svolti'],
  ['palestra', 'Palestra'],
  ['trappole', 'Errori tipici'],
  ['orale', 'Domande d\'orale'],
  ['quiz', 'Quiz finale'],
];

const TRASFORMAZIONI = [
  { name: 'Fusione del ghiaccio', answer: 'f', why: 'Cambia lo stato di aggregazione, non la sostanza: prima e dopo è H₂O.' },
  { name: 'Arrugginimento del ferro', answer: 'c', why: 'Il ferro reagisce con ossigeno e acqua e forma ossidi idrati di ferro(III): una sostanza nuova, con proprietà diverse.' },
  { name: 'Combustione del metano', answer: 'c', why: 'CH₄ + 2 O₂ → CO₂ + 2 H₂O: i reagenti spariscono e compaiono sostanze nuove.' },
  { name: 'Sublimazione dello iodio', answer: 'f', why: 'I₂ solido diventa I₂ gassoso: le molecole restano le stesse, cambia solo il modo in cui sono disposte.' },
  { name: 'Zucchero che si scioglie nel tè', answer: 'f', why: 'Le molecole di saccarosio si disperdono fra quelle dell\'acqua senza cambiare: facendo evaporare l\'acqua si recupera lo zucchero.' },
  { name: 'Cottura di un uovo', answer: 'c', why: 'Le proteine si denaturano e formano nuovi legami fra loro: raffreddando l\'uovo non torna crudo.' },
  { name: 'Elettrolisi dell\'acqua', answer: 'c', why: '2 H₂O → 2 H₂ + O₂: il composto viene scomposto nei suoi elementi.' },
  { name: 'Ebollizione dell\'acqua', answer: 'f', why: 'Il vapore è ancora H₂O; condensandolo si riottiene acqua liquida.' },
  { name: 'Annerimento dell\'argenteria', answer: 'c', why: 'L\'argento reagisce con i composti dello zolfo presenti nell\'aria e forma solfuro d\'argento, Ag₂S, nero.' },
  { name: 'Macinare il sale grosso', answer: 'f', why: 'Cambiano forma e dimensione dei granuli, non la composizione: resta NaCl.' },
  { name: 'Fermentazione del mosto', answer: 'c', why: 'I lieviti trasformano gli zuccheri in etanolo e diossido di carbonio.' },
  { name: 'Magnetizzazione di un chiodo', answer: 'f', why: 'Si orientano i domini magnetici del ferro, ma il materiale resta ferro.' },
];

const INTENSIVE = [
  { name: 'Massa', answer: 'e', why: 'Raddoppia se raddoppi il campione.' },
  { name: 'Volume', answer: 'e', why: 'Dipende dalla quantità di materia.' },
  { name: 'Densità', answer: 'i', why: 'È il rapporto fra due grandezze estensive (m/V): raddoppiando il campione raddoppiano entrambe e il rapporto non cambia.' },
  { name: 'Temperatura', answer: 'i', why: 'Due bicchieri d\'acqua a 20 °C versati insieme danno acqua a 20 °C, non a 40 °C.' },
  { name: 'Temperatura di ebollizione', answer: 'i', why: 'È una caratteristica della sostanza (a pressione fissata), utile per riconoscerla.' },
  { name: 'Colore', answer: 'i', why: 'Non dipende da quanto campione hai.' },
  { name: 'Capacità termica (calore per scaldare il campione di 1 °C)', answer: 'e', why: 'Serve più calore per scaldare più materia. Il calore specifico, riferito a 1 g, è invece intensivo.' },
  { name: 'Calore specifico', answer: 'i', why: 'È riferito all\'unità di massa (J g⁻¹ °C⁻¹): non dipende dalla quantità.' },
  { name: 'Pressione di un gas', answer: 'i', why: 'Forza per unità di superficie: non si somma unendo due porzioni dello stesso gas alle stesse condizioni.' },
  { name: 'Lunghezza di una barra', answer: 'e', why: 'Due barre messe in fila sono lunghe il doppio.' },
  { name: 'Concentrazione di una soluzione', answer: 'i', why: 'Un cucchiaio e un litro della stessa soluzione hanno la stessa concentrazione.' },
];

const CLASSI = [
  { name: 'Diamante', answer: 'el', why: 'È formato solo da atomi di carbonio, C.' },
  { name: 'Ozono, O₃', answer: 'el', why: 'Molecole fatte di un solo elemento: è una forma (allotropo) dell\'ossigeno, non un composto.' },
  { name: 'Acqua distillata', answer: 'co', why: 'H₂O: idrogeno e ossigeno combinati in rapporto fisso (11,19% di H in massa), separabili solo con una reazione.' },
  { name: 'Cloruro di sodio puro', answer: 'co', why: 'NaCl: due elementi in rapporto fisso.' },
  { name: 'Aria secca filtrata', answer: 'om', why: 'N₂, O₂, Ar, CO₂… in un\'unica fase gassosa; la composizione può variare (per esempio l\'umidità).' },
  { name: 'Ottone', answer: 'om', why: 'Lega di rame e zinco: una soluzione solida, uniforme, a composizione variabile.' },
  { name: 'Acqua di mare filtrata', answer: 'om', why: 'Sali disciolti in acqua: una soluzione, una sola fase.' },
  { name: 'Granito', answer: 'et', why: 'Si distinguono a occhio cristalli diversi (quarzo, feldspati, miche).' },
  { name: 'Sabbia in acqua', answer: 'et', why: 'Due fasi visibili, separabili per filtrazione o decantazione.' },
  { name: 'Latte', answer: 'et', why: 'Sembra uniforme, ma al microscopio si vedono goccioline di grasso disperse: è un colloide (emulsione), classificato fra le miscele eterogenee.' },
  { name: 'Ghiaccio che galleggia in acqua', answer: 'co', why: 'Due fasi, ma una sola sostanza, H₂O: è una sostanza pura. Più fasi non bastano per avere una miscela.' },
  { name: 'Mercurio', answer: 'el', why: 'Elemento metallico, liquido a temperatura ambiente.' },
  { name: 'Saccarosio', answer: 'co', why: 'C₁₂H₂₂O₁₁: composto di carbonio, idrogeno e ossigeno.' },
  { name: 'Aceto', answer: 'om', why: 'Soluzione di acido acetico in acqua (circa 4–6%) con altre sostanze disciolte.' },
  { name: 'Sangue', answer: 'et', why: 'Cellule in sospensione nel plasma: si separano per centrifugazione.' },
];

export const DOMANDE = [
  { q: 'Quante cifre significative ha 0,004050?', options: ['3', '4', '6', '7'], correct: 1, why: 'Gli zeri iniziali (0,00) non contano; 4, 0 intermedio, 5 e lo zero finale dopo la virgola sì: 4 c.s.' },
  { q: 'Il risultato di 4,51 × 3,6666 scritto correttamente è:', options: ['16,536', '16,54', '16,5', '17'], correct: 2, why: 'Nel prodotto vale il dato con meno c.s.: 4,51 ne ha 3, quindi 16,536… → 16,5.' },
  { q: '12,11 + 18,0 + 1,013 = ?', options: ['31,123', '31,12', '31,1', '31'], correct: 2, why: 'Nella somma conta la posizione: 18,0 si ferma ai decimi, quindi anche il risultato: 31,1.' },
  { q: 'Quale di queste è una proprietà intensiva?', options: ['Massa', 'Volume', 'Densità', 'Capacità termica'], correct: 2, why: 'La densità non dipende dalla quantità di campione; massa, volume e capacità termica sì.' },
  { q: 'Quale di questi processi è una trasformazione chimica?', options: ['Sublimazione dello iodio', 'Arrugginimento del ferro', 'Fusione della cera', 'Dissoluzione del sale in acqua'], correct: 1, why: 'Solo l\'arrugginimento produce una sostanza nuova (ossidi idrati di ferro).' },
  { q: 'L\'ottone (lega rame-zinco) è:', options: ['un elemento', 'un composto', 'una miscela omogenea', 'una miscela eterogenea'], correct: 2, why: 'È una soluzione solida: una sola fase, composizione variabile.' },
  { q: 'Un ghiacciolo d\'acqua pura che galleggia in acqua pura è:', options: ['una miscela eterogenea, perché ci sono due fasi', 'una sostanza pura', 'una miscela omogenea', 'un elemento'], correct: 1, why: 'C\'è una sola sostanza, H₂O, in due stati di aggregazione. Le fasi non decidono se è una miscela: decide il numero di sostanze.' },
  { q: '1 m³ corrisponde a:', options: ['10 L', '100 L', '1000 L', '10⁶ L'], correct: 2, why: '1 m³ = (10 dm)³ = 10³ dm³ e 1 dm³ = 1 L.' },
  { q: '1 cm³ in m³ è:', options: ['10⁻² m³', '10⁻³ m³', '10⁻⁶ m³', '10⁻⁹ m³'], correct: 2, why: '1 cm = 10⁻² m, quindi 1 cm³ = (10⁻² m)³ = 10⁻⁶ m³. L\'esponente si moltiplica per 3.' },
  { q: '37,0 °C corrispondono a:', options: ['310 K', '310,2 K', '310,15 K', '236,2 K'], correct: 1, why: 'T = 37,0 + 273,15 = 310,15; il dato si ferma ai decimi, quindi 310,2 K (273,15 è esatto).' },
  { q: 'Una bilancia fornisce 10,52 g, 10,51 g e 10,53 g per una massa campione certificata di 10,00 g. Le misure sono:', options: ['precise e accurate', 'precise ma non accurate', 'accurate ma non precise', 'né precise né accurate'], correct: 1, why: 'Sono molto vicine fra loro (precise) ma tutte lontane dal valore vero: c\'è un errore sistematico, per esempio la bilancia non tarata.' },
  { q: 'Dal 20 maggio 2019 il chilogrammo è definito fissando il valore numerico:', options: ['della massa del prototipo internazionale di Sèvres', 'della costante di Planck h', 'del numero di Avogadro', 'della massa di 1 L d\'acqua a 4 °C'], correct: 1, why: 'h = 6,626 070 15 × 10⁻³⁴ J s esatti. Il prototipo di platino-iridio non è più la definizione.' },
  { q: 'Quante cifre significative ha "1500 m" scritto così?', options: ['2', '4', 'È ambiguo: 2, 3 o 4', 'Infinite: è un numero intero'], correct: 2, why: 'Gli zeri finali di un intero senza virgola sono ambigui. Si toglie l\'ambiguità con la notazione scientifica: 1,50 × 10³ m ha 3 c.s.' },
  { q: 'Un cubetto di 54,0 g occupa 20,0 cm³. Probabilmente è di:', options: ['magnesio (1,74 g/cm³)', 'alluminio (2,70 g/cm³)', 'ferro (7,87 g/cm³)', 'piombo (11,3 g/cm³)'], correct: 1, why: 'd = 54,0 / 20,0 = 2,70 g/cm³.' },
  { q: 'Quale di questi numeri è esatto, cioè senza incertezza?', options: ['La massa di una moneta, 2,30 g', 'Il numero di studenti in aula, 32', 'La temperatura dell\'aula, 21 °C', 'Il volume di una goccia, 0,05 mL'], correct: 1, why: 'I numeri ottenuti contando oggetti (e quelli fissati per definizione, come 1 in = 2,54 cm) sono esatti: hanno infinite cifre significative.' },
  { q: 'Per separare il sale dall\'acqua di mare si usa:', options: ['la filtrazione', 'la decantazione', 'l\'evaporazione o la distillazione', 'la separazione magnetica'], correct: 2, why: 'Il sale è disciolto (miscela omogenea): passa attraverso il filtro. Si sfrutta la diversa volatilità: l\'acqua evapora, il sale resta.' },
  { q: 'Durante l\'ebollizione a pressione costante, una sostanza pura:', options: ['si scalda in modo continuo', 'resta a temperatura costante finché tutto il liquido è evaporato', 'si raffredda', 'si decompone'], correct: 1, why: 'Il calore fornito serve a vaporizzare il liquido. Una soluzione invece bolle in un intervallo di temperatura, perché la sua composizione cambia mentre il solvente evapora.' },
  { q: 'Su un\'unità di misura, che cosa cambia fra la Terra e la Luna per uno stesso oggetto?', options: ['La massa', 'Il peso', 'Sia massa sia peso', 'Né massa né peso'], correct: 1, why: 'La massa misura la quantità di materia e non cambia; il peso è la forza m·g e g sulla Luna è circa 1,62 m/s² invece di 9,81 m/s².' },
];

const html = String.raw;

const PAGE = html`
<section id="obiettivi" class="sec">
  <h2>Cosa devi saper fare</h2>
  <p class="lead">Questo è il capitolo degli strumenti: unità, misure e numeri. Tutti gli esercizi del corso (stechiometria, gas, equilibri, pH) richiedono conversioni e cifre significative corrette. Spunta un obiettivo quando sai farlo senza guardare gli appunti.</p>
  <ul class="checklist" id="chk">
    <li><label><input type="checkbox" id="o1"> Distinguere legge, ipotesi e teoria</label></li>
    <li><label><input type="checkbox" id="o2"> Distinguere proprietà fisiche e chimiche, intensive ed estensive</label></li>
    <li><label><input type="checkbox" id="o3"> Classificare la materia: elementi, composti, miscele omogenee ed eterogenee</label></li>
    <li><label><input type="checkbox" id="o4"> Conoscere le 7 unità SI di base e i prefissi da tera a femto</label></li>
    <li><label><input type="checkbox" id="o5"> Convertire fra °C, K e °F</label></li>
    <li><label><input type="checkbox" id="o6"> Usare densità e percentuali come fattori di conversione</label></li>
    <li><label><input type="checkbox" id="o7"> Spiegare precisione e accuratezza, errori casuali e sistematici</label></li>
    <li><label><input type="checkbox" id="o8"> Contare le cifre significative e applicarle nei calcoli</label></li>
    <li><label><input type="checkbox" id="o9"> Risolvere problemi con l'analisi dimensionale, anche con unità al quadrato e al cubo</label></li>
  </ul>
</section>

<section id="metodo" class="sec">
  <h2>Il metodo scientifico</h2>
  <p>La chimica studia la composizione, la struttura e le proprietà della materia e le sue trasformazioni. Come le altre scienze procede per un ciclo di osservazioni ed esperimenti.</p>
  <dl class="defs">
    <dt>Osservazione e dato</dt><dd>Ciò che si misura o si osserva. Può essere <b>qualitativo</b> (il liquido è blu) o <b>quantitativo</b> (la massa è 2,35 g).</dd>
    <dt>Ipotesi</dt><dd>Una spiegazione provvisoria delle osservazioni, che si deve poter mettere alla prova con esperimenti.</dd>
    <dt>Legge scientifica</dt><dd>Un enunciato sintetico, spesso matematico, che <b>descrive</b> un comportamento osservato sempre allo stesso modo. Esempio: nelle reazioni chimiche la massa totale si conserva (Lavoisier).</dd>
    <dt>Teoria</dt><dd>Un modello che <b>spiega</b> le leggi e permette nuove previsioni. Esempio: la teoria atomica spiega perché la massa si conserva (gli atomi si riorganizzano ma non si creano né si distruggono).</dd>
  </dl>
  <p>Una teoria resta valida finché le sue previsioni sono confermate. Basta un esperimento riproducibile in contraddizione per doverla modificare o abbandonare. Un esperimento è utile se è <b>controllato</b>: si cambia una variabile alla volta.</p>
  <aside class="trap"><b>All'orale.</b> "Legge" e "teoria" non indicano un grado di certezza diverso. La legge dice <i>che cosa</i> succede, la teoria <i>perché</i>.</aside>
</section>

<section id="proprieta" class="sec">
  <h2>Proprietà e trasformazioni della materia</h2>
  <p>La <b>materia</b> è tutto ciò che ha massa e occupa spazio. Si descrive con la sua <b>composizione</b> (quali parti e in che proporzioni: l'acqua è 11,19% di idrogeno e 88,81% di ossigeno in massa) e con la sua <b>struttura</b> (come sono disposti gli atomi).</p>
  <div class="two">
    <div>
      <h3>Proprietà fisiche</h3>
      <p>Si osservano senza cambiare la composizione: colore, stato fisico, densità, temperatura di fusione e di ebollizione, durezza, conducibilità elettrica.</p>
      <p>Una <b>trasformazione fisica</b> cambia l'aspetto o lo stato, non la sostanza: il ghiaccio fonde ma resta H₂O.</p>
    </div>
    <div>
      <h3>Proprietà chimiche</h3>
      <p>Descrivono la capacità di una sostanza di trasformarsi in altre: il ferro arrugginisce, il metano brucia, il sodio reagisce con l'acqua.</p>
      <p>Una <b>trasformazione chimica</b> (reazione) produce sostanze nuove con proprietà diverse.</p>
    </div>
  </div>
  <h3>Fisica o chimica?</h3>
  <div class="widget" id="w-trasf"></div>

  <h3>Proprietà intensive ed estensive</h3>
  <p>Una proprietà <b>estensiva</b> dipende dalla quantità di materia (massa, volume, capacità termica): si somma unendo due campioni. Una proprietà <b>intensiva</b> non dipende dalla quantità (densità, temperatura, colore, temperatura di ebollizione): serve a riconoscere una sostanza.</p>
  <p>Il rapporto di due grandezze estensive è intensivo: densità = massa/volume, calore specifico = capacità termica/massa.</p>
  <div class="widget" id="w-int"></div>
</section>

<section id="stati" class="sec">
  <h2>Stati di aggregazione</h2>
  <div class="table-wrap"><table>
    <thead><tr><th>Stato</th><th>Forma</th><th>Volume</th><th>Particelle</th></tr></thead>
    <tbody>
      <tr><td>Solido</td><td>propria</td><td>proprio</td><td>a contatto, in posizioni fisse; vibrano attorno a esse</td></tr>
      <tr><td>Liquido</td><td>del recipiente</td><td>proprio</td><td>a contatto ma libere di scorrere le une sulle altre</td></tr>
      <tr><td>Gas</td><td>del recipiente</td><td>del recipiente</td><td>lontane fra loro rispetto alle loro dimensioni, in moto rapido e disordinato; il gas è comprimibile</td></tr>
    </tbody>
  </table></div>
  <figure class="fig">
    <svg viewBox="0 0 600 340" role="img" aria-labelledby="stati-title">
      <title id="stati-title">I passaggi di stato fra solido, liquido e gas</title>
      <defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="svg-ink"/></marker></defs>
      <g class="svg-node"><rect x="10" y="250" width="150" height="56" rx="6"/><text x="85" y="285">Solido</text></g>
      <g class="svg-node"><rect x="440" y="250" width="150" height="56" rx="6"/><text x="515" y="285">Liquido</text></g>
      <g class="svg-node"><rect x="225" y="16" width="150" height="56" rx="6"/><text x="300" y="51">Gas</text></g>
      <g class="svg-arrows">
        <path d="M165 265 L435 265" marker-end="url(#arr)"/><text x="300" y="254" class="svg-lbl hot">fusione</text>
        <path d="M435 292 L165 292" marker-end="url(#arr)"/><text x="300" y="322" class="svg-lbl cold">solidificazione</text>
        <path d="M500 245 L360 77" marker-end="url(#arr)"/><text x="452" y="140" class="svg-lbl hot" text-anchor="start">vaporizzazione</text>
        <path d="M335 80 L470 245" marker-end="url(#arr)"/><text x="392" y="215" class="svg-lbl cold" text-anchor="end">condensazione</text>
        <path d="M100 245 L240 77" marker-end="url(#arr)"/><text x="148" y="140" class="svg-lbl hot" text-anchor="end">sublimazione</text>
        <path d="M265 80 L130 245" marker-end="url(#arr)"/><text x="208" y="185" class="svg-lbl cold" text-anchor="start">brinamento</text>
      </g>
    </svg>
    <figcaption>In arancio i passaggi che assorbono calore, in verde quelli che lo cedono. La vaporizzazione può avvenire in superficie a qualunque temperatura (evaporazione) o in tutto il liquido (ebollizione). Il passaggio gas → solido si chiama anche deposizione.</figcaption>
  </figure>
  <p>Durante un passaggio di stato di una sostanza pura, a pressione costante, la temperatura resta costante finché il passaggio non è completo: il calore fornito serve a separare le particelle, non ad aumentarne l'agitazione.</p>
</section>

<section id="classificazione" class="sec">
  <h2>Classificazione della materia</h2>
  <div class="tree" role="img" aria-label="Schema: la materia si divide in sostanze pure (elementi e composti) e miscele (omogenee ed eterogenee)">
    <div class="tree-root">Materia</div>
    <div class="tree-branches">
      <div class="tree-col">
        <div class="tree-node strong">Sostanze pure<small>composizione fissa, proprietà costanti</small></div>
        <div class="tree-pair">
          <div class="tree-node">Elementi<small>un solo tipo di atomo: Fe, O₂, C (diamante)</small></div>
          <div class="tree-node">Composti<small>più elementi in rapporto fisso: H₂O, NaCl</small></div>
        </div>
        <p class="tree-note">composto ⇄ elementi: solo con trasformazioni chimiche</p>
      </div>
      <div class="tree-col">
        <div class="tree-node strong">Miscele<small>composizione variabile</small></div>
        <div class="tree-pair">
          <div class="tree-node">Omogenee (soluzioni)<small>una sola fase: aria, acqua salata, ottone</small></div>
          <div class="tree-node">Eterogenee<small>più fasi distinguibili: granito, sabbia in acqua</small></div>
        </div>
        <p class="tree-note">miscela ⇄ sostanze pure: con metodi fisici</p>
      </div>
    </div>
  </div>
  <p>Un <b>elemento</b> non si può scomporre in sostanze più semplici con mezzi chimici. Un <b>composto</b> contiene due o più elementi combinati in proporzioni fisse (legge delle proporzioni definite, Proust), e ha proprietà diverse da quelle degli elementi che lo formano: il sodio è un metallo che reagisce con l'acqua, il cloro un gas tossico, il cloruro di sodio è il sale da cucina.</p>
  <p>Una <b>miscela</b> ha composizione variabile e conserva in parte le proprietà dei componenti. È <b>omogenea</b> se ha la stessa composizione in ogni punto (una sola fase), <b>eterogenea</b> se si distinguono zone diverse. Una <b>fase</b> è una porzione di materia uniforme, separata dalle altre da una superficie netta.</p>
  <aside class="trap"><b>Trappola.</b> Il numero di fasi non decide se c'è una sostanza pura: ghiaccio in acqua sono due fasi ma una sola sostanza. Una sostanza pura fonde e bolle a temperatura costante; una soluzione bolle in un intervallo di temperature.</aside>
  <h3>Metodi di separazione delle miscele</h3>
  <div class="table-wrap"><table>
    <thead><tr><th>Metodo</th><th>Proprietà sfruttata</th><th>Esempio</th></tr></thead>
    <tbody>
      <tr><td>Filtrazione</td><td>dimensione delle particelle (solido non disciolto)</td><td>sabbia da acqua</td></tr>
      <tr><td>Decantazione, centrifugazione</td><td>densità</td><td>globuli rossi dal plasma</td></tr>
      <tr><td>Distillazione</td><td>volatilità (temperatura di ebollizione)</td><td>acqua dolce da acqua di mare; etanolo dal vino</td></tr>
      <tr><td>Cristallizzazione</td><td>solubilità che cambia con la temperatura</td><td>purificazione di un sale</td></tr>
      <tr><td>Cromatografia</td><td>affinità diversa per una fase fissa e una mobile</td><td>pigmenti di un inchiostro</td></tr>
      <tr><td>Estrazione con solvente</td><td>solubilità in solventi non miscibili</td><td>caffeina dal caffè</td></tr>
      <tr><td>Separazione magnetica</td><td>proprietà magnetiche</td><td>limatura di ferro da zolfo</td></tr>
    </tbody>
  </table></div>
  <h3>Classifica tu</h3>
  <div class="widget" id="w-classi"></div>
</section>

<section id="si" class="sec">
  <h2>Il Sistema Internazionale (SI)</h2>
  <p>Una misura è un numero <b>e</b> un'unità: "2,5" non significa nulla, "2,5 m" sì. Dal 20 maggio 2019 le sette unità di base sono definite fissando il valore esatto di sette costanti della natura.</p>
  <div class="table-wrap"><table class="si">
    <thead><tr><th>Grandezza</th><th>Unità</th><th>Simbolo</th><th>Costante fissata</th></tr></thead>
    <tbody>
      <tr><td>Lunghezza</td><td>metro</td><td class="sym">m</td><td>velocità della luce, c = 299 792 458 m s⁻¹</td></tr>
      <tr><td>Massa</td><td>chilogrammo</td><td class="sym">kg</td><td>costante di Planck, h = 6,626 070 15 × 10⁻³⁴ J s</td></tr>
      <tr><td>Tempo</td><td>secondo</td><td class="sym">s</td><td>frequenza iperfine del cesio-133, Δν<sub>Cs</sub> = 9 192 631 770 Hz</td></tr>
      <tr><td>Temperatura</td><td>kelvin</td><td class="sym">K</td><td>costante di Boltzmann, k = 1,380 649 × 10⁻²³ J K⁻¹</td></tr>
      <tr><td>Quantità di sostanza</td><td>mole</td><td class="sym">mol</td><td>costante di Avogadro, N<sub>A</sub> = 6,022 140 76 × 10²³ mol⁻¹</td></tr>
      <tr><td>Corrente elettrica</td><td>ampere</td><td class="sym">A</td><td>carica elementare, e = 1,602 176 634 × 10⁻¹⁹ C</td></tr>
      <tr><td>Intensità luminosa</td><td>candela</td><td class="sym">cd</td><td>efficacia luminosa K<sub>cd</sub> = 683 lm W⁻¹ a 540 × 10¹² Hz</td></tr>
    </tbody>
  </table></div>
  <p class="note">Il chilogrammo è l'unica unità di base con un prefisso nel nome. Prima del 2019 era la massa di un cilindro di platino-iridio conservato a Sèvres; la mole erano gli atomi contenuti in 12 g di carbonio-12. Alcune edizioni dei libri riportano ancora le definizioni precedenti.</p>
  <h3>Prefissi</h3>
  <div class="prefixes">
    ${[['tera', 'T', 12], ['giga', 'G', 9], ['mega', 'M', 6], ['kilo', 'k', 3], ['etto', 'h', 2], ['deca', 'da', 1], ['deci', 'd', -1], ['centi', 'c', -2], ['milli', 'm', -3], ['micro', 'µ', -6], ['nano', 'n', -9], ['pico', 'p', -12], ['femto', 'f', -15]]
      .map(([n, s, e]) => `<div class="pfx"><span class="pfx-s">${s}</span><span class="pfx-n">${n}</span><span class="pfx-e">10${String(e).replace('-', '⁻').replace(/\d/g, d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d])}</span></div>`).join('')}
  </div>
  <h3>Unità derivate e unità d'uso comune</h3>
  <ul>
    <li><b>Volume</b>: m³. In chimica si usano il <b>litro</b>, 1 L = 1 dm³ (esatto), e il millilitro, 1 mL = 1 cm³.</li>
    <li><b>Densità</b>: kg/m³; in laboratorio g/cm³ = g/mL per solidi e liquidi, g/L per i gas.</li>
    <li><b>Lunghezze atomiche</b>: pm, nm o ångström (1 Å = 10⁻¹⁰ m = 100 pm, non SI).</li>
    <li><b>Pressione</b> e <b>energia</b> si vedranno con i gas e la termodinamica: Pa = N/m², J = N·m = kg m² s⁻².</li>
  </ul>
  <aside class="trap"><b>Massa e peso.</b> La massa misura la quantità di materia e non dipende dal luogo. Il peso è la forza con cui la gravità attrae un corpo, P = m·g: sulla Terra g ≈ 9,81 m/s², sulla Luna circa 1,62 m/s². In laboratorio la bilancia confronta masse, per questo si parla di "massa" anche se si dice "pesare".</aside>
</section>

<section id="temperatura" class="sec">
  <h2>Scale di temperatura</h2>
  <div class="formula-row">
    <div class="formula">T(K) = t(°C) + 273,15</div>
    <div class="formula">t(°F) = <span class="frac"><span>9</span><span>5</span></span> t(°C) + 32</div>
  </div>
  <p>Lo zero della scala Kelvin è lo <b>zero assoluto</b> (0 K = −273,15 °C = −459,67 °F): una temperatura più bassa non esiste, per questo non ci sono valori negativi in kelvin. Il simbolo è K, senza il segno di grado. Una <b>differenza</b> di 1 K è uguale a una differenza di 1 °C, ma vale 1,8 °F.</p>
  <div class="table-wrap"><table class="num-table">
    <thead><tr><th>Punto di riferimento</th><th>K</th><th>°C</th><th>°F</th></tr></thead>
    <tbody>
      <tr><td>Zero assoluto</td><td>0</td><td>−273,15</td><td>−459,67</td></tr>
      <tr><td>Fusione del ghiaccio (1 atm)</td><td>273,15</td><td>0</td><td>32</td></tr>
      <tr><td>Temperatura corporea</td><td>310,2</td><td>37,0</td><td>98,6</td></tr>
      <tr><td>Ebollizione dell'acqua (1 atm)</td><td>373,15</td><td>100</td><td>212</td></tr>
    </tbody>
  </table></div>
  <p class="note">273,15, 9/5 e 32 sono esatti per definizione: nelle conversioni non limitano le cifre significative. Aggiungere 273,15 è una somma, quindi si conserva la posizione dell'ultima cifra del dato (25,0 °C → 298,2 K).</p>
  <h3>Convertitore con i passaggi</h3>
  <div class="widget conv" id="w-temp">
    <div class="conv-row">
      <label>Valore<input id="t-v" value="-195,79" inputmode="decimal" autocomplete="off"></label>
      <label>Da<select id="t-from">${SCALE_T.map(s => `<option${s === '°C' ? ' selected' : ''}>${s}</option>`).join('')}</select></label>
      <label>A<select id="t-to">${SCALE_T.map(s => `<option${s === 'K' ? ' selected' : ''}>${s}</option>`).join('')}</select></label>
    </div>
    <div class="conv-out" id="t-out" aria-live="polite"></div>
  </div>
</section>

<section id="densita" class="sec">
  <h2>Densità e composizione percentuale</h2>
  <div class="formula-row"><div class="formula">d = <span class="frac"><span>m</span><span>V</span></span></div></div>
  <p>La densità è la massa per unità di volume. È una proprietà intensiva e caratteristica di una sostanza, ma <b>dipende dalla temperatura</b> (quasi tutte le sostanze si dilatano scaldandole): l'acqua ha la densità massima, 0,99997 g/cm³, a circa 4 °C, e 0,99705 g/cm³ a 25 °C. Il ghiaccio (0,917 g/cm³) è meno denso dell'acqua liquida e galleggia.</p>
  <p>Nei problemi la densità è un <b>fattore di conversione</b> fra massa e volume: si scrive come <span class="frac"><span>0,789 g</span><span>1 mL</span></span> per passare da volume a massa, o capovolta per il passaggio inverso.</p>
  <div class="table-wrap"><table class="num-table dens">
    <caption>Densità a 20–25 °C (CRC Handbook)</caption>
    <thead><tr><th>Sostanza</th><th>g/cm³</th></tr></thead>
    <tbody>
      <tr><td>Aria (20 °C, 1 atm)</td><td>0,00120</td></tr><tr><td>Etanolo</td><td>0,789</td></tr><tr><td>Acqua (20 °C)</td><td>0,998</td></tr>
      <tr><td>Ghiaccio (0 °C)</td><td>0,917</td></tr><tr><td>Magnesio</td><td>1,74</td></tr><tr><td>Alluminio</td><td>2,70</td></tr><tr><td>Titanio</td><td>4,51</td></tr>
      <tr><td>Ferro</td><td>7,87</td></tr><tr><td>Rame</td><td>8,96</td></tr><tr><td>Argento</td><td>10,5</td></tr><tr><td>Piombo</td><td>11,3</td></tr>
      <tr><td>Mercurio</td><td>13,53</td></tr><tr><td>Oro</td><td>19,3</td></tr><tr><td>Osmio (il più denso)</td><td>22,6</td></tr>
    </tbody>
  </table></div>
  <h3>Composizione percentuale in massa</h3>
  <div class="formula-row"><div class="formula">% in massa di A = <span class="frac"><span>massa di A</span><span>massa totale</span></span> × 100</div></div>
  <p>Anche la percentuale è un fattore di conversione: "ottone al 67,0% di rame" significa <span class="frac"><span>67,0 g rame</span><span>100 g ottone</span></span>. Il 100 è esatto.</p>
  <p>Per un oggetto irregolare il volume si misura per <b>spostamento d'acqua</b> in un cilindro graduato: V = V<sub>finale</sub> − V<sub>iniziale</sub>. Attenzione: questa sottrazione riduce le cifre significative (vedi gli esercizi svolti).</p>
</section>

<section id="incertezza" class="sec">
  <h2>Incertezza: precisione e accuratezza</h2>
  <p>Ogni misura ha un'incertezza. Si scrivono tutte le cifre certe più <b>una cifra stimata</b>: con un righello graduato in millimetri si legge 12,35 cm (il 5 è stimato fra due tacche); con una buretta graduata ogni 0,1 mL si legge 23,46 mL.</p>
  <div class="two">
    <div>
      <h3>Precisione</h3>
      <p>Quanto sono vicine fra loro misure ripetute. Dipende dagli <b>errori casuali</b>, che cambiano segno da una misura all'altra e si riducono facendo la media di molte misure. Si quantifica con la deviazione standard s.</p>
    </div>
    <div>
      <h3>Accuratezza</h3>
      <p>Quanto il risultato è vicino al valore vero. È peggiorata dagli <b>errori sistematici</b>, che spostano tutte le misure nello stesso verso (strumento non tarato, metodo sbagliato) e <b>non</b> si riducono ripetendo la misura.</p>
    </div>
  </div>
  <div class="formula-row">
    <div class="formula">x̄ = <span class="frac"><span>Σ x<sub>i</sub></span><span>n</span></span></div>
    <div class="formula">s = √<span class="over"><span class="frac"><span>Σ (x<sub>i</sub> − x̄)²</span><span>n − 1</span></span></span></div>
  </div>
  <h3>Prova: una bilancia e una massa certificata di 10,000 g</h3>
  <div class="widget" id="w-prec">
    <div class="prec-presets" role="group" aria-label="Situazioni tipiche">
      <button type="button" class="chip" data-p="0,0.012">Precisa e accurata</button>
      <button type="button" class="chip" data-p="0.22,0.012">Precisa, non accurata</button>
      <button type="button" class="chip" data-p="0,0.13">Accurata, non precisa</button>
      <button type="button" class="chip" data-p="0.25,0.14">Né l'una né l'altra</button>
    </div>
    <div class="prec-controls">
      <label for="p-sys">Errore sistematico <output id="p-sys-o"></output></label>
      <input type="range" id="p-sys" min="-0.3" max="0.3" step="0.01" value="0.22">
      <label for="p-rnd">Errore casuale (dispersione) <output id="p-rnd-o"></output></label>
      <input type="range" id="p-rnd" min="0.005" max="0.2" step="0.005" value="0.015">
      <button type="button" class="btn-quiet" id="p-new">Ripeti le 20 misure</button>
    </div>
    <svg id="p-svg" viewBox="0 0 640 190" role="img" aria-label="Venti misure della massa sulla scala da 9,4 a 10,6 g"></svg>
    <div class="prec-stats" id="p-stats" aria-live="polite"></div>
  </div>
</section>

<section id="cifre" class="sec">
  <h2>Cifre significative</h2>
  <p>Le cifre significative (c.s.) sono le cifre certe di una misura più la prima incerta. Scrivere "2,0 g" o "2,00 g" non è la stessa cosa: il secondo dato è dieci volte più preciso.</p>
  <ol class="rules">
    <li>Le cifre diverse da zero sono sempre significative: 4,56 → 3 c.s.</li>
    <li>Gli zeri fra cifre non nulle sono significativi: 4006 → 4 c.s.; 2,05 → 3 c.s.</li>
    <li>Gli zeri iniziali non sono significativi: indicano solo la posizione della virgola. 0,00453 → 3 c.s. (= 4,53 × 10⁻³)</li>
    <li>Gli zeri finali dopo la virgola sono significativi: 0,0450 → 3 c.s.; 2,00 → 3 c.s.</li>
    <li>Gli zeri finali di un intero senza virgola sono ambigui: 1500 può avere 2, 3 o 4 c.s. Si scrive 1,5 × 10³, 1,50 × 10³ o 1,500 × 10³.</li>
  </ol>
  <p><b>Numeri esatti</b>: i conteggi (12 uova, 3 atomi in H₂O) e le definizioni (1 m = 100 cm, 1 in = 2,54 cm, 1 L = 1000 mL) hanno infinite cifre significative e non limitano mai il risultato.</p>
  <h3>Contatore</h3>
  <div class="widget" id="w-sig">
    <label for="sig-in" class="lbl">Scrivi un numero (virgola o punto; notazione scientifica come 1,50e3 o 1,50 × 10^3)</label>
    <input id="sig-in" value="0,004050" inputmode="decimal" autocomplete="off" spellcheck="false">
    <div class="chips" role="group" aria-label="Esempi">
      ${['0,004050', '4006', '1500', '1500,', '1,50 × 10^3', '100,0', '0,1020', '6,022 × 10^23'].map(x => `<button type="button" class="chip" data-x="${x}">${x}</button>`).join('')}
    </div>
    <div id="sig-out" class="sig-out" aria-live="polite"></div>
  </div>
  <h3>Cifre significative nei calcoli</h3>
  <div class="two">
    <div class="rule-card">
      <h4>Moltiplicazione e divisione</h4>
      <p>Il risultato ha tante c.s. quante il dato che ne ha <b>meno</b>.</p>
      <p class="num">14,79 × 12,11 × 5,05 = 904,48… → <b>904</b> <span class="note">(5,05 ha 3 c.s.)</span></p>
    </div>
    <div class="rule-card">
      <h4>Addizione e sottrazione</h4>
      <p>Il risultato si ferma alla <b>stessa posizione decimale</b> del dato meno preciso, indipendentemente dal numero di c.s.</p>
      <p class="num">15,02 + 9986,0 + 3,518 = 10 004,538 → <b>10 004,5</b> <span class="note">(9986,0 si ferma ai decimi)</span></p>
    </div>
  </div>
  <p><b>Arrotondamento.</b> Se la prima cifra da eliminare è 5 o più si aumenta di 1 l'ultima cifra tenuta, altrimenti si lascia (2,346 → 2,35; 2,344 → 2,34). Alcuni testi, per un 5 esatto seguito solo da zeri, arrotondano alla cifra pari: segui la convenzione del tuo docente.</p>
  <aside class="trap"><b>Si arrotonda solo alla fine.</b> Nei calcoli a più passaggi tieni almeno una cifra in più (cifra di guardia) e arrotonda il risultato finale. Arrotondare a ogni passaggio accumula errori: è uno degli sbagli più frequenti negli scritti.</aside>
</section>

<section id="conversioni" class="sec">
  <h2>Analisi dimensionale (fattori di conversione)</h2>
  <p>Un fattore di conversione è una frazione che vale 1, perché numeratore e denominatore sono la stessa quantità in unità diverse: <span class="frac"><span>100 cm</span><span>1 m</span></span> = 1. Moltiplicando per un fattore il valore non cambia, cambiano solo le unità. Si sceglie il verso del fattore in modo che l'unità da eliminare stia al denominatore e si semplifichi.</p>
  <div class="formula-row"><div class="formula">grandezza cercata = grandezza data × fattore di conversione</div></div>
  <p><b>Unità al quadrato o al cubo</b>: il fattore va elevato alla stessa potenza. 1 m = 10² cm, quindi 1 m² = 10⁴ cm² e 1 m³ = 10⁶ cm³. Dimenticare l'esponente è l'errore più comune.</p>
  <p><b>Metodo</b>: 1) scrivi il dato con la sua unità; 2) scrivi l'unità del risultato; 3) metti in fila i fattori finché le unità si semplificano; 4) controlla che l'ordine di grandezza sia ragionevole; 5) arrotonda alle c.s. dei dati misurati.</p>
  <h3>Convertitore con i fattori</h3>
  <div class="widget conv" id="w-conv">
    <div class="conv-row">
      <label>Grandezza<select id="c-q">${Object.entries(UNITA).map(([k, g]) => `<option value="${k}">${g.nome}</option>`).join('')}</select></label>
      <label>Valore<input id="c-v" value="2,50" inputmode="decimal" autocomplete="off"></label>
      <label>Da<select id="c-from"></select></label>
      <label>A<select id="c-to"></select></label>
    </div>
    <div class="conv-out" id="c-out" aria-live="polite"></div>
  </div>
</section>

<section id="svolti" class="sec">
  <h2>Esercizi svolti</h2>
  <p class="lead">Prova a risolverli da solo, poi scopri un passaggio alla volta.</p>

  <article class="svolto">
    <h3>1 · Densità per spostamento d'acqua</h3>
    <p class="testo">Un oggetto metallico di 25,46 g viene immerso in un cilindro graduato: l'acqua sale da 20,0 mL a 23,2 mL. Calcola la densità del metallo. Che metallo può essere?</p>
    <ol class="passi">
      <li>Volume dell'oggetto: V = 23,2 mL − 20,0 mL = 3,2 mL = 3,2 cm³. È una sottrazione: il risultato si ferma ai decimi e ha solo <b>2 c.s.</b></li>
      <li>d = m / V = 25,46 g / 3,2 cm³ = 7,956… g/cm³</li>
      <li>Nella divisione vale il dato con meno c.s. (3,2 → 2 c.s.): <b>d = 8,0 g/cm³</b>.</li>
      <li>Il valore è compatibile con il ferro (7,87 g/cm³), ma con 2 c.s. non si può escludere l'ottone (8,4–8,7 g/cm³). La misura del volume limita la precisione: servirebbe un cilindro più fine o un oggetto più grande.</li>
    </ol>
  </article>

  <article class="svolto">
    <h3>2 · Dalla massa al volume</h3>
    <p class="testo">Un serbatoio contiene 1,20 × 10³ kg di etanolo (d = 0,789 g/cm³). Qual è il volume in litri?</p>
    <ol class="passi">
      <li>Catena: kg → g → cm³ → L. La densità va capovolta perché si parte dalla massa.</li>
      <li class="chain">1,20 × 10³ <span class="cancel">kg</span> × <span class="frac"><span>10³ <span class="cancel">g</span></span><span>1 <span class="cancel">kg</span></span></span> × <span class="frac"><span>1 <span class="cancel">cm³</span></span><span>0,789 <span class="cancel">g</span></span></span> × <span class="frac"><span>1 L</span><span>10³ <span class="cancel">cm³</span></span></span> = 1,5209… × 10³ L</li>
      <li>Dati misurati: 1,20 × 10³ (3 c.s.) e 0,789 (3 c.s.); i fattori 10³ sono esatti. <b>V = 1,52 × 10³ L</b>.</li>
      <li>Controllo: l'etanolo è meno denso dell'acqua, quindi 1200 kg devono occupare più di 1200 L. ✓</li>
    </ol>
  </article>

  <article class="svolto">
    <h3>3 · Temperatura dell'azoto liquido</h3>
    <p class="testo">L'azoto bolle a −195,79 °C. Esprimi la temperatura in kelvin e in gradi Fahrenheit.</p>
    <ol class="passi">
      <li>T = −195,79 + 273,15 = <b>77,36 K</b> (somma: si conservano i centesimi del dato; 273,15 è esatto).</li>
      <li>9/5 × (−195,79) = −352,422. Il prodotto per un numero esatto conserva le 5 c.s. del dato: −352,42.</li>
      <li>−352,42 + 32 = <b>−320,42 °F</b> (32 è esatto).</li>
    </ol>
  </article>

  <article class="svolto">
    <h3>4 · Operazioni miste e cifra di guardia</h3>
    <p class="testo">Calcola (12,5 + 3,456) × 2,1 con il numero corretto di cifre significative.</p>
    <ol class="passi">
      <li>Somma: 12,5 + 3,456 = 15,956. Il risultato si fermerebbe ai decimi (16,0, cioè 3 c.s.), ma <b>non arrotondi ancora</b>: tieni 15,956 e ricorda che vale 3 c.s.</li>
      <li>Prodotto: 15,956 × 2,1 = 33,5076.</li>
      <li>Fattori: 3 c.s. (la somma) e 2 c.s. (2,1) → 2 c.s.: <b>34</b>.</li>
    </ol>
  </article>

  <article class="svolto">
    <h3>5 · Volume di un cubo e unità al cubo</h3>
    <p class="testo">Un cubo di alluminio ha il lato di 2,50 cm (d = 2,70 g/cm³). Calcola la massa, poi esprimi il volume in m³.</p>
    <ol class="passi">
      <li>V = (2,50 cm)³ = 15,625 cm³ (3 c.s.: 15,6 cm³, ma si tiene la cifra di guardia).</li>
      <li>m = d × V = 2,70 g/cm³ × 15,625 cm³ = 42,1875 g → <b>42,2 g</b>.</li>
      <li>1 cm = 10⁻² m, quindi 1 cm³ = (10⁻² m)³ = 10⁻⁶ m³: V = 15,6 × 10⁻⁶ m³ = <b>1,56 × 10⁻⁵ m³</b>.</li>
    </ol>
  </article>

  <article class="svolto">
    <h3>6 · Percentuale come fattore</h3>
    <p class="testo">L'acqua di mare contiene in media il 3,5% in massa di sali disciolti. Quanti kg di acqua di mare bisogna evaporare per ottenere 1,00 kg di sali?</p>
    <ol class="passi">
      <li>Fattore: <span class="frac"><span>100 kg acqua di mare</span><span>3,5 kg sali</span></span> (capovolto, perché si parte dai sali).</li>
      <li>1,00 kg × 100 / 3,5 = 28,57… kg</li>
      <li>3,5 ha 2 c.s.: <b>29 kg</b>.</li>
    </ol>
  </article>

  <article class="svolto">
    <h3>7 · Unità composte</h3>
    <p class="testo">La velocità del suono nell'aria a 20 °C è 343 m/s. Quanto vale in km/h?</p>
    <ol class="passi">
      <li class="chain">343 <span class="frac"><span><span class="cancel">m</span></span><span><span class="cancel">s</span></span></span> × <span class="frac"><span>1 km</span><span>10³ <span class="cancel">m</span></span></span> × <span class="frac"><span>3600 <span class="cancel">s</span></span><span>1 h</span></span> = 1234,8 km/h</li>
      <li>3600 s/h è esatto: restano le 3 c.s. di 343 → <b>1,23 × 10³ km/h</b>. Scrivere "1230 km/h" sarebbe ambiguo.</li>
    </ol>
  </article>
</section>

<section id="palestra" class="sec">
  <h2>Palestra</h2>
  <p class="lead">Esercizi sempre nuovi, generati a caso. La correzione controlla sia il valore sia le cifre significative: se il numero è giusto ma arrotondato male te lo segnala. Puoi scrivere 1,23e3, 1,23 × 10^3 o 1230.</p>
  <div class="widget" id="w-palestra"></div>
</section>

<section id="trappole" class="sec">
  <h2>Errori tipici all'esame</h2>
  <ul class="traps">
    <li><b>Dimenticare l'esponente nelle unità al cubo.</b> 1 m³ non è 100 cm³ ma 10⁶ cm³.</li>
    <li><b>Arrotondare a metà calcolo.</b> Tieni una cifra di guardia e arrotonda solo il risultato finale.</li>
    <li><b>Applicare la regola sbagliata.</b> Somme e sottrazioni: posizione decimale. Prodotti e divisioni: numero di c.s.</li>
    <li><b>Contare le c.s. dei numeri esatti.</b> 273,15, 100 nelle percentuali, i fattori fra prefissi e i conteggi non limitano il risultato.</li>
    <li><b>Scrivere zeri finali ambigui.</b> "1200 g" non dice quante c.s. hai: scrivi 1,20 × 10³ g.</li>
    <li><b>Capovolgere il fattore nel verso sbagliato.</b> Controlla che l'unità da eliminare stia al denominatore.</li>
    <li><b>Confondere omogeneo con puro.</b> L'acqua salata è omogenea ma è una miscela; ghiaccio in acqua ha due fasi ma è una sostanza pura.</li>
    <li><b>Scrivere "°K".</b> Il kelvin non ha il simbolo di grado.</li>
    <li><b>Dimenticare l'unità.</b> Un risultato senza unità di misura è considerato sbagliato.</li>
  </ul>
</section>

<section id="orale" class="sec">
  <h2>Domande d'orale</h2>
  <p class="lead">Rispondi a voce, poi apri la risposta modello e confronta.</p>
  <div class="oral">
    <details><summary>Che differenza c'è fra una legge e una teoria scientifica?</summary>
      <p>La legge descrive in modo sintetico un comportamento che si osserva sempre (la massa si conserva nelle reazioni). La teoria è un modello che spiega le leggi e fa nuove previsioni (gli atomi si riorganizzano senza crearsi né distruggersi). Una teoria resta valida finché nessun esperimento riproducibile la contraddice.</p></details>
    <details><summary>Perché la densità è una proprietà intensiva se massa e volume sono estensive?</summary>
      <p>Perché è il loro rapporto: se raddoppio il campione raddoppiano sia la massa sia il volume e il rapporto resta uguale. Per questo la densità caratterizza una sostanza (a temperatura e pressione date) e non la quantità che ne ho.</p></details>
    <details><summary>Come distingui sperimentalmente una sostanza pura da una soluzione?</summary>
      <p>Una sostanza pura ha composizione fissa e proprietà costanti: fonde e bolle a una temperatura costante (a pressione fissata). Una soluzione bolle in un intervallo di temperature, perché mentre il solvente evapora la composizione del liquido cambia. Inoltre i componenti di una miscela si separano con metodi fisici (distillazione, cristallizzazione), mentre un composto si scompone solo con una reazione.</p></details>
    <details><summary>Qual è la differenza fra elemento e composto?</summary>
      <p>Un elemento è formato da un solo tipo di atomo e non si scompone con mezzi chimici. Un composto è formato da due o più elementi combinati in proporzioni fisse (legge di Proust), e ha proprietà diverse da quelle degli elementi che lo formano. Esempio: Na (metallo reattivo) e Cl₂ (gas tossico) formano NaCl, il sale da cucina.</p></details>
    <details><summary>Spiega precisione e accuratezza con un esempio.</summary>
      <p>La precisione è l'accordo fra misure ripetute; l'accuratezza è la vicinanza al valore vero. Se una bilancia non tarata dà 10,52, 10,51 e 10,53 g per una massa vera di 10,00 g, le misure sono precise ma non accurate: c'è un errore sistematico, che non si elimina ripetendo la misura ma solo tarando lo strumento. Gli errori casuali invece peggiorano la precisione e si riducono mediando molte misure.</p></details>
    <details><summary>Che cosa sono le cifre significative e perché sono importanti?</summary>
      <p>Sono le cifre certe di una misura più la prima incerta: dicono con che precisione è stata fatta la misura. Un risultato calcolato non può essere più preciso dei dati: nei prodotti e nei quozienti si tengono le c.s. del dato che ne ha meno, nelle somme e nelle differenze la posizione decimale del dato meno preciso. I numeri esatti (conteggi, definizioni) non limitano il risultato.</p></details>
    <details><summary>Perché non esistono temperature negative nella scala Kelvin?</summary>
      <p>Lo zero della scala Kelvin è lo zero assoluto, −273,15 °C, la temperatura più bassa possibile. La temperatura è legata all'energia di agitazione delle particelle, che non può scendere sotto il minimo. Una variazione di 1 K è uguale a una variazione di 1 °C; T(K) = t(°C) + 273,15.</p></details>
    <details><summary>Quali sono le unità di base del SI e come sono definite oggi?</summary>
      <p>Metro, chilogrammo, secondo, kelvin, mole, ampere e candela. Dal 2019 ciascuna è definita fissando il valore esatto di una costante: c per il metro, h per il chilogrammo, la frequenza del cesio per il secondo, k di Boltzmann per il kelvin, N<sub>A</sub> per la mole, la carica elementare e per l'ampere, K<sub>cd</sub> per la candela. Le altre unità (newton, joule, pascal, litro) derivano da queste.</p></details>
  </div>
</section>

<section id="quiz" class="sec">
  <h2>Quiz finale</h2>
  <p class="lead">${DOMANDE.length} domande in ordine casuale. Il miglior risultato resta salvato in questo browser e compare nel programma.</p>
  <div id="w-quiz"></div>
</section>
`;

export function render(root, { onScore } = {}) {
  root.innerHTML = PAGE;
  checklist(root.querySelector('#chk'), 'materia');
  sorter(root.querySelector('#w-trasf'), { items: TRASFORMAZIONI, options: [{ id: 'f', label: 'Fisica' }, { id: 'c', label: 'Chimica' }], label: 'Trasformazione' });
  sorter(root.querySelector('#w-int'), { items: INTENSIVE, options: [{ id: 'i', label: 'Intensiva' }, { id: 'e', label: 'Estensiva' }], label: 'Proprietà' });
  sorter(root.querySelector('#w-classi'), { items: CLASSI, options: [{ id: 'el', label: 'Elemento' }, { id: 'co', label: 'Composto' }, { id: 'om', label: 'Miscela omogenea' }, { id: 'et', label: 'Miscela eterogenea' }], label: 'Classe' });
  temperatureWidget(root);
  converterWidget(root);
  sigWidget(root);
  precisionWidget(root);
  steppers(root);
  practice(root.querySelector('#w-palestra'), { kinds: GENERATORI });
  quiz(root.querySelector('#w-quiz'), DOMANDE, { topic: 'materia', onScore });
}

function temperatureWidget(root) {
  const $ = (s) => root.querySelector(s);
  const draw = () => {
    const r = convertTemperature($('#t-v').value, $('#t-from').value, $('#t-to').value);
    const out = $('#t-out');
    if (!r) { out.innerHTML = '<p class="warn">Scrivi un numero, per esempio −195,79.</p>'; return; }
    const below = ($('#t-from').value === 'K' && parseNumber($('#t-v').value).value < 0);
    out.innerHTML = `${below ? '<p class="warn">Una temperatura in kelvin non può essere negativa.</p>' : ''}<ol class="tight">${r.steps.map(s => `<li>${s}</li>`).join('')}</ol>`;
  };
  ['#t-v', '#t-from', '#t-to'].forEach(s => $(s).addEventListener('input', draw));
  draw();
}

function converterWidget(root) {
  const $ = (s) => root.querySelector(s);
  const fill = () => {
    const units = Object.keys(UNITA[$('#c-q').value].u);
    const defaults = { lunghezza: ['km', 'cm'], massa: ['mg', 'kg'], volume: ['cm³', 'm³'], densita: ['g/cm³', 'kg/m³'] }[$('#c-q').value];
    $('#c-from').innerHTML = units.map(u => `<option${u === defaults[0] ? ' selected' : ''}>${u}</option>`).join('');
    $('#c-to').innerHTML = units.map(u => `<option${u === defaults[1] ? ' selected' : ''}>${u}</option>`).join('');
  };
  const draw = () => {
    const group = $('#c-q').value, from = $('#c-from').value, to = $('#c-to').value;
    const c = convert($('#c-v').value, from, to, group);
    const out = $('#c-out');
    if (!c) { out.innerHTML = '<p class="warn">Scrivi un numero, per esempio 2,50 oppure 3,1e-4.</p>'; return; }
    const notes = [...new Set([from, to])].map(u => NOTE_UNITA[u]).filter(Boolean);
    out.innerHTML = `<p class="chain">${c.chain} = <b>${c.text} ${to}</b></p>
      <p class="note">${c.ambiguous ? `Il dato ha zeri finali ambigui (da ${c.sigMin} a ${c.sig} c.s.): qui li considero misurati. Meglio scriverlo in notazione scientifica. ` : `Il risultato conserva le ${c.sig} c.s. del dato: i fattori sono esatti. `}${notes.join(' ')}</p>`;
  };
  $('#c-q').addEventListener('change', () => { fill(); draw(); });
  ['#c-v', '#c-from', '#c-to'].forEach(s => $(s).addEventListener('input', draw));
  fill(); draw();
}

function sigWidget(root) {
  const input = root.querySelector('#sig-in'), out = root.querySelector('#sig-out');
  const draw = () => {
    const raw = input.value;
    const p = parseNumber(raw);
    if (!p) { out.innerHTML = '<p class="warn">Non riesco a leggerlo come numero. Esempi validi: 0,0450 · 1500 · 1,50e3 · 6,022 × 10^23.</p>'; return; }
    if (p.sig === null) { out.innerHTML = '<p class="warn">Lo zero da solo non ha un numero definito di cifre significative.</p>'; return; }
    let k = 0;
    const digit = (c) => `<span class="dg mk-${p.marks[k++]}">${c}</span>`;
    const mant = `${p.sign === '-' ? '−' : ''}${[...p.intPart].map(digit).join('') || '<span class="dg">0</span>'}${p.hasPoint ? ',' : ''}${[...p.fracPart].map(digit).join('')}`;
    const expo = p.sci ? ` × 10<sup>${String(p.exp).replace('-', '−')}</sup>` : '';
    const count = p.ambiguous ? `ambiguo: da <b>${p.sig}</b> a <b>${p.sigMax}</b> c.s.` : `<b>${p.sig}</b> ${p.sig === 1 ? 'cifra significativa' : 'cifre significative'}`;
    const sci = p.ambiguous
      ? [...new Set([p.sig, p.sigMax])].map(n => formatSig(p.value, n, { forceSci: true })).join(' oppure ')
      : formatSig(p.value, p.sig, { forceSci: true });
    out.innerHTML = `<p class="sig-num">${mant}${expo}</p><p class="sig-count">${count}</p>${digitsExplanation(p)}<p class="note">In notazione scientifica: ${sci}${p.ambiguous ? ', a seconda di quante cifre sono state misurate' : ''}.</p>`;
  };
  input.addEventListener('input', draw);
  root.querySelector('#w-sig .chips').addEventListener('click', (ev) => {
    const b = ev.target.closest('[data-x]');
    if (b) { input.value = b.dataset.x; draw(); }
  });
  draw();
}

function precisionWidget(root) {
  const $ = (s) => root.querySelector(s);
  const TRUE = 10, N = 20, LO = 9.4, HI = 10.6;
  let shots = [];
  const gauss = () => { let u = 0, v = 0; while (!u) u = Math.random(); while (!v) v = Math.random(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
  const sample = () => {
    const sys = +$('#p-sys').value, rnd = +$('#p-rnd').value;
    shots = Array.from({ length: N }, () => ({ x: TRUE + sys + rnd * gauss(), y: rand(0.15, 0.85) }));
    draw();
  };
  const X = (v) => 30 + (Math.min(HI, Math.max(LO, v)) - LO) / (HI - LO) * 580;
  const draw = () => {
    $('#p-sys-o').textContent = `${(+$('#p-sys').value >= 0 ? '+' : '−')}${formatPlain(Math.abs(+$('#p-sys').value))} g`.replace('.', ',');
    $('#p-rnd-o').textContent = `± ${String(+$('#p-rnd').value).replace('.', ',')} g`;
    const xs = shots.map(s => s.x);
    const mean = xs.reduce((a, b) => a + b, 0) / N;
    const s = Math.sqrt(xs.reduce((a, b) => a + (b - mean) ** 2, 0) / (N - 1));
    const ticks = [];
    for (let v = 9.4; v <= 10.61; v += 0.4) ticks.push(v);
    $('#p-svg').innerHTML = `
      <rect x="30" y="20" width="580" height="100" class="svg-band"/>
      ${ticks.map(v => `<line x1="${X(v)}" x2="${X(v)}" y1="120" y2="126" class="svg-axis"/><text x="${X(v)}" y="148" class="svg-tick">${v.toFixed(1).replace('.', ',')}</text>`).join('')}
      <line x1="30" x2="610" y1="120" y2="120" class="svg-axis"/>
      <line x1="${X(TRUE)}" x2="${X(TRUE)}" y1="20" y2="120" class="svg-true"/>
      <text x="${X(TRUE)}" y="14" class="svg-tick">valore vero</text>
      ${shots.map(p => `<circle cx="${X(p.x)}" cy="${20 + p.y * 100}" r="7" class="svg-shot"/>`).join('')}
      <line x1="${X(mean)}" x2="${X(mean)}" y1="20" y2="120" class="svg-mean"/>
      <text x="${X(mean)}" y="176" class="svg-tick svg-mean-t">media</text>`;
    const err = mean - TRUE;
    const precise = s <= 0.05, accurate = Math.abs(err) <= 0.05;
    const verdict = precise && accurate ? 'precise e accurate' : precise ? 'precise ma non accurate' : accurate ? 'accurate ma poco precise' : 'né precise né accurate';
    $('#p-stats').innerHTML = `
      <dl class="stats">
        <div><dt>Media x̄</dt><dd>${formatPlain(mean, 5)} g</dd></div>
        <div><dt>Deviazione standard s</dt><dd>${formatPlain(s, 2)} g</dd></div>
        <div><dt>x̄ − valore vero</dt><dd>${err >= 0 ? '+' : ''}${formatPlain(err, 2)} g</dd></div>
      </dl>
      <p>Le misure sono <b>${verdict}</b>. ${!accurate && precise ? 'Ripetere la misura non serve: tutte sono spostate dallo stesso errore sistematico. Bisogna tarare la bilancia.' : !precise && accurate ? 'La media è vicina al valore vero, ma le singole misure sono disperse: con più misure la media diventa più affidabile.' : !precise && !accurate ? 'Ci sono sia un errore sistematico sia una grande dispersione.' : 'Media vicina al valore vero e misure raccolte.'} <span class="note">(Qui "precise" vuol dire s ≤ 0,05 g, "accurate" |x̄ − 10,000| ≤ 0,05 g.)</span></p>`;
  };
  $('#p-sys').addEventListener('input', sample);
  $('#p-rnd').addEventListener('input', sample);
  $('#p-new').addEventListener('click', sample);
  root.querySelector('.prec-presets').addEventListener('click', (ev) => {
    const b = ev.target.closest('[data-p]');
    if (!b) return;
    const [a, c] = b.dataset.p.split(',');
    $('#p-sys').value = a; $('#p-rnd').value = c;
    sample();
  });
  sample();
}
