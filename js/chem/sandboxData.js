// Molecole, atomi ed esperimenti pronti per la sandbox di dinamica molecolare reattiva.

export const SANDBOX_MOLECULES = [
  { id: 'H2', smiles: '[H][H]', name: 'idrogeno' },
  { id: 'O2', smiles: 'O=O', name: 'ossigeno' },
  { id: 'N2', smiles: 'N#N', name: 'azoto' },
  { id: 'H2O', smiles: 'O', name: 'acqua' },
  { id: 'CH4', smiles: 'C', name: 'metano' },
  { id: 'CO2', smiles: 'O=C=O', name: 'anidride carbonica' },
  { id: 'NH3', smiles: 'N', name: 'ammoniaca' },
  { id: 'Cl2', smiles: 'ClCl', name: 'cloro' },
  { id: 'HCl', smiles: 'Cl', name: 'acido cloridrico' },
  { id: 'F2', smiles: 'FF', name: 'fluoro' },
  { id: 'HF', smiles: 'F', name: 'acido fluoridrico' },
  { id: 'H2O2', smiles: 'OO', name: 'perossido di idrogeno' },
  { id: 'C2H6', smiles: 'CC', name: 'etano' },
  { id: 'C2H4', smiles: 'C=C', name: 'etilene' },
  { id: 'C2H2', smiles: 'C#C', name: 'acetilene' },
  { id: 'C6H6', smiles: 'c1ccccc1', name: 'benzene' },
  { id: 'CH3OH', smiles: 'CO', name: 'metanolo' },
  { id: 'HCN', smiles: 'C#N', name: 'acido cianidrico' },
  { id: 'NaCl', smiles: '[Na]Cl', name: 'cloruro di sodio (gas)' },
  { id: 'CH3', smiles: '[CH3]', name: 'radicale metile' },
  { id: 'OH', smiles: '[OH]', name: 'radicale ossidrile' },
];

// Ioni e specie cariche: richiedono un motore quantistico (MINDO/3 o Hartree–Fock), che riceve la carica totale.
// MINDO/3 non ha i parametri della coppia Cl–O: Cl⁻ non si può mettere con acqua o ossigeno.
export const SANDBOX_IONS = [
  { id: 'H3O+', smiles: '[OH3+]', name: 'ione ossonio (idronio)' },
  { id: 'OH-', smiles: '[OH-]', name: 'ione idrossido' },
  { id: 'NH4+', smiles: '[NH4+]', name: 'ione ammonio' },
  { id: 'NH2-', smiles: '[NH2-]', name: 'ione ammiduro' },
  { id: 'H+', smiles: '[H+]', name: 'protone nudo' },
  { id: 'H-', smiles: '[H-]', name: 'ione idruro' },
  { id: 'CH3+', smiles: '[CH3+]', name: 'catione metile (carbocatione)' },
  { id: 'CH3-', smiles: '[CH3-]', name: 'anione metile (carbanione)' },
  { id: 'NO2+', smiles: 'O=[N+]=O', name: 'ione nitronio, l\'elettrofilo della nitrazione' },
  { id: 'NO3-', smiles: '[O-][N+](=O)[O-]', name: 'ione nitrato' },
  { id: 'CH3O-', smiles: 'C[O-]', name: 'ione metossido' },
  { id: 'F-', smiles: '[F-]', name: 'ione fluoruro' },
  { id: 'Cl-', smiles: '[Cl-]', name: 'ione cloruro (non con O: manca la coppia Cl–O)' },
];

export const SANDBOX_ELEMENTS = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 35, 36, 53, 54];

// Nomi delle specie riconosciute durante la simulazione (formula come la scrive hillFormula).
export const SPECIES_NAMES = {
  // ioni (la carica è letta dalla funzione d'onda dei motori quantistici)
  'H₃O⁺': 'ione ossonio', 'HO⁻': 'ione idrossido', 'NH₄⁺': 'ione ammonio', 'H₂N⁻': 'ione ammiduro', 'H⁺': 'protone',
  'H⁻': 'ione idruro', 'CH₃⁺': 'catione metile', 'CH₃⁻': 'anione metile', 'NO₂⁺': 'ione nitronio', 'NO₃⁻': 'ione nitrato',
  'CH₃O⁻': 'ione metossido', 'F⁻': 'ione fluoruro', 'Cl⁻': 'ione cloruro', 'H₂O⁺': 'catione radicale dell\'acqua',
  'C₆H₆NO₂⁺': 'complesso benzene–nitronio (σ o π)', 'C₆H₅NO₂': 'nitrobenzene', 'C₆H₆': 'benzene', 'NO₂': 'diossido di azoto',
  'H': 'idrogeno atomico', 'H₂': 'idrogeno', 'O': 'ossigeno atomico', 'O₂': 'ossigeno', 'O₃': 'ozono',
  'HO': 'radicale ossidrile ·OH', 'HO₂': 'radicale idroperossile', 'H₂O': 'acqua', 'H₂O₂': 'perossido di idrogeno',
  'N': 'azoto atomico', 'N₂': 'azoto', 'NH₃': 'ammoniaca', 'H₂N': 'radicale amminico', 'HN': 'imidogeno', 'N₂H₄': 'idrazina',
  'Cl': 'cloro atomico', 'Cl₂': 'cloro', 'HCl': 'acido cloridrico', 'F': 'fluoro atomico', 'F₂': 'fluoro', 'HF': 'acido fluoridrico',
  'CH₄': 'metano', 'CH₃': 'radicale metile', 'CH₂': 'metilene', 'C₂H₆': 'etano', 'C₂H₄': 'etilene', 'C₂H₂': 'acetilene',
  'C₂H₅': 'radicale etile', 'CO₂': 'anidride carbonica', 'CO': 'monossido di carbonio', 'CH₂O': 'formaldeide',
  'CH₄O': 'metanolo', 'CHO': 'radicale formile', 'CH₂O₂': 'acido formico', 'C₆H₆': 'benzene', 'CHN': 'acido cianidrico',
  'ClNa': 'cloruro di sodio', 'NaCl': 'cloruro di sodio', 'Na': 'sodio', 'Na₂': 'sodio biatomico',
  'He': 'elio', 'Ne': 'neon', 'Ar': 'argon', 'Kr': 'kripton', 'Xe': 'xeno', 'C': 'carbonio atomico',
};

/** 8 molecole d'acqua su un reticolo 2×2×2 (O–O = 3,1 Å, densità 1,0 g/cm³ in una scatola di 6,2 Å), orientazioni diverse. */
function waterLattice() {
  const out = [];
  const r = 0.957, th = 104.5 * Math.PI / 180;
  let k = 0;
  for (const x of [-1.55, 1.55]) for (const y of [-1.55, 1.55]) for (const z of [-1.55, 1.55]) {
    const a = 0.9 * k, b = 1.7 * k + 0.4;
    // due vettori O–H nel piano ruotato di (a, b)
    const u = [Math.cos(a), Math.sin(a) * Math.cos(b), Math.sin(a) * Math.sin(b)];
    const w0 = [-Math.sin(a), Math.cos(a) * Math.cos(b), Math.cos(a) * Math.sin(b)];
    const h1 = [0, 1, 2].map(c => r * (Math.cos(th / 2) * u[c] + Math.sin(th / 2) * w0[c]));
    const h2 = [0, 1, 2].map(c => r * (Math.cos(th / 2) * u[c] - Math.sin(th / 2) * w0[c]));
    out.push({ Z: 8, pos: [x, y, z] }, { Z: 1, pos: [x + h1[0], y + h1[1], z + h1[2]] }, { Z: 1, pos: [x + h2[0], y + h2[1], z + h2[2]] });
    k++;
  }
  return out;
}

// Esperimenti: ogni riga di "add" è [SMILES o simbolo, numero di copie].
export const PRESETS = [
  // --- dinamica quantistica: a ogni passo una SCF MINDO/3 (elettroni trattati con la meccanica quantistica) ---
  {
    id: 'q-h2-o2', name: 'Combustione dell\'idrogeno', forceField: 'mindo3', quantum: true,
    box: 10, T: 3000, add: [['[H][H]', 8], ['O=O', 4], ['[O]', 2]], style: 'cloud',
    text: 'Dinamica quantistica: a ogni passo si risolvono le equazioni di Schrödinger (SCF MINDO/3) per tutti gli elettroni di valenza e le forze sui nuclei vengono dalla funzione d\'onda. Due atomi di ossigeno innescano la catena radicalica: H₂ + O → OH + H, H + O₂ → HO₂, OH + H₂ → H₂O + H. Nessuna di queste reazioni è programmata: emergono dagli elettroni.',
    tips: ['Stile "nuvola elettronica": vedi la densità degli elettroni spostarsi quando un legame si rompe o si forma.', 'Stile "orbitali": l\'HOMO dei radicali (OH, H, O) è l\'orbitale spaiato che attacca le altre molecole.', 'Lo spin totale (radicali, O₂ tripletto) è calcolato, non imposto.'],
  },
  {
    id: 'q-h2-cl2', name: 'Idrogeno e cloro alla luce', forceField: 'mindo3', quantum: true,
    box: 10, T: 900, add: [['[H][H]', 6], ['ClCl', 6]], light: { lambda: 300, on: true, rate: 3 }, style: 'cloud',
    text: 'Impulsi meccanici di energia hc/λ (4,1 eV a 300 nm) depositati nei legami rompono alcune molecole di Cl₂: è un deposito di energia, non un modello di assorbimento della luce. Gli atomi di cloro avviano la catena Cl + H₂ → HCl + H, H + Cl₂ → HCl + Cl, con le forze calcolate dalla meccanica quantistica a ogni passo.',
    tips: ['Colora per carica: nel prodotto HCl il cloro prende carica negativa.'],
  },
  {
    id: 'q-methyl', name: 'Radicali metile che si uniscono', forceField: 'mindo3', quantum: true,
    box: 11, T: 600, add: [['[CH3]', 8]], style: 'orbital',
    text: 'Ogni radicale CH₃ ha un elettrone spaiato in un orbitale p del carbonio (è il suo HOMO). Quando due radicali con spin opposti si incontrano, i due orbitali si sovrappongono e formano il legame σ C–C dell\'etano: 2 CH₃ → C₂H₆. Se gli spin sono paralleli (stato di tripletto) il principio di Pauli impedisce il legame e i radicali si respingono: un effetto puramente quantistico, per cui in media solo un urto su quattro può formare il legame.',
    tips: ['Stile "orbitali": guarda l\'HOMO passare da un orbitale p isolato a un legame fra due carboni.'],
  },
  {
    id: 'q-water', name: 'Acqua liquida quantistica', forceField: 'mindo3', quantum: true,
    box: 6.2, T: 300, atoms: waterLattice(), thermalize: 300, thermostat: true, speed: 20, style: 'cloud',
    text: 'Otto molecole d\'acqua alla densità del liquido (1,0 g/cm³), con la struttura elettronica ricalcolata a ogni passo: le cariche, la polarizzazione delle molecole e i legami a idrogeno cambiano con la geometria. Limite noto: MINDO/3 sottostima i legami a idrogeno.',
    tips: ['Colora per carica: O negativo, H positivo, calcolati dalla funzione d\'onda.'],
  },
  {
    id: 'q-neutral', name: 'Neutralizzazione: H₃O⁺ + OH⁻', forceField: 'mindo3', quantum: true,
    box: 8, T: 300, dt: 0.1, add: [['[OH3+]', 1], ['[OH-]', 1]], thermostat: false, speed: 30,
    text: 'Uno ione ossonio e uno ione idrossido in fase gassosa, con la struttura elettronica calcolata a ogni passo e la carica di ogni specie letta dalla funzione d\'onda. L\'attrazione elettrostatica li avvicina, il protone passa e si formano due molecole d\'acqua: H₃O⁺ + OH⁻ → 2 H₂O. Senza solvente l\'energia liberata (circa 10 eV, l\'affinità protonica di OH⁻ meno quella di H₂O) resta nelle due molecole, che si scaldano a migliaia di kelvin: in acqua il solvente la disperderebbe come calore di neutralizzazione (57 kJ/mol).',
    tips: ['Sistema isolato: guarda salire la temperatura quando il protone passa.', 'Colora per carica: la carica positiva si sposta con il protone.', 'Il bilancio dell\'energia può mostrare un salto di qualche decimo di eV quando, avvicinandosi, gli ioni passano a un altro stato elettronico (il gap HOMO–LUMO cambia di colpo): la dinamica segue un solo stato, quello più basso, e non descrive transizioni fra stati.'],
  },
  {
    id: 'q-wheland', name: 'Nitrazione del benzene: l\'intermedio di Wheland perde H⁺', forceField: 'mindo3', quantum: true,
    box: 14, T: 300, thermalize: 300, thermostat: true, speed: 20, style: 'ball',
    atoms: [[8, -3.134, -0.1045, -0.7237], [7, -2.243, -0.117, 0.0904], [8, -2.1615, -0.2046, 1.2976], [6, -0.8523, 0.0053, -0.5204], [6, -0.0185, -1.2344, -0.3517], [6, 1.2613, -1.2034, 0.1618], [6, 1.8356, 0.0405, 0.5469], [6, 1.1429, 1.2723, 0.3787], [6, -0.1385, 1.2703, -0.1319], [1, -0.9033, 0.1026, -1.6497], [1, -0.458, -2.1935, -0.6856], [1, 1.8456, -2.1326, 0.2655], [1, 2.8564, 0.052, 0.9729], [1, 1.6374, 2.2205, 0.6466], [1, -0.6701, 2.2265, -0.2975], [7, -0.9888, 0.2656, -3.5408], [1, -1.0058, -0.6383, -3.9997], [1, -1.8191, 0.7634, -3.8422], [1, -0.1926, 0.7697, -3.9151]].map(([Z, x, y, z], i) => ({ Z, pos: [x, y, z], formal: i === 3 ? 1 : 0 })),
    text: 'Sostituzione elettrofila aromatica, secondo passo. Si parte dal complesso σ (intermedio di Wheland) C₆H₆NO₂⁺: l\'elettrofilo NO₂⁺ è legato a un carbonio diventato sp³, l\'anello ha perso l\'aromaticità e porta la carica positiva. Un\'ammoniaca vicina strappa il protone di quel carbonio: si ottengono nitrobenzene e NH₄⁺ e l\'anello torna aromatico. Geometria iniziale ottimizzata con MINDO/3, dinamica quantistica a ogni passo. Nel laboratorio la base è HSO₄⁻ o l\'acqua.',
    tips: ['Colora per carica: la carica positiva passa dall\'anello allo ione ammonio.', 'Stile "orbitali": l\'HOMO del nitrobenzene è di nuovo un orbitale π delocalizzato sull\'anello.', 'Primo passo: in "Scatola vuota" aggiungi benzene e NO₂⁺ (lista degli ioni). MINDO/3 lega il complesso σ solo di circa 3 kcal/mol (la chimica reale molto di più): a 300 K resta spesso un complesso di incontro.'],
  },
  {
    id: 'q-field', name: 'Molecole d\'acqua in un campo elettrico', forceField: 'mindo3', quantum: true,
    box: 10, T: 150, add: [['O', 6]], field: [0, 0, 0.5], thermostat: true, speed: 20,
    text: 'Un campo elettrico uniforme di 0,5 V/Å lungo z entra nella funzione d\'onda di MINDO/3: polarizza le molecole e ne orienta i dipoli. L\'energia −μ·E è minima con il dipolo parallelo al campo, cioè con gli idrogeni (estremo positivo) verso +z e l\'ossigeno verso −z. L\'agitazione termica si oppone all\'allineamento. Cambia intensità e direzione del campo nei controlli: il lavoro fatto dal campo è contato nel bilancio dell\'energia.',
    tips: ['Confronta il dipolo totale (pannello) a campo acceso e spento.', 'Abbassa la temperatura: l\'allineamento aumenta (legge di Langevin–Debye).'],
  },
  {
    id: 'argon', name: 'Gas di argon: Maxwell–Boltzmann',
    box: 42, T: 300, dt: 2, add: [['[Ar]', 120]],
    text: "Confronta la distribuzione delle velocità e la pressione alle pareti con Maxwell–Boltzmann e nkT/V. Con il modello automatico l'argon usa Lennard–Jones (σ = 3,405 Å, ε/k = 119,8 K): attrazioni e pareti morbide producono piccoli scarti dal gas ideale.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
  },
  {
    id: 'argon-liquid', name: 'Argon: condensazione e cristallo',
    box: 22, T: 60, dt: 2, add: [['[Ar]', 160]],
    text: "Esplora aggregazione e ordinamento con Lennard–Jones. La fase è letta dalla struttura locale (q₆ di Steinhardt): goccia, coesistenza con il vapore o cristallo. Il condensarsi libera calore: con il termostato spento la temperatura sale. Taglio, pareti e numero finito di atomi spostano le temperature di transizione rispetto all'argon reale.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
  },
  {
    id: 'water', name: 'Acqua: legami a idrogeno',
    box: 15.5, T: 300, add: [['O', 40]],
    text: "Osserva geometria e interazioni tra molecole nel potenziale classico. Non è un modello validato delle proprietà dell’acqua liquida.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
    color: 'charge',
  },
  {
    id: 'h2-o2', name: 'Combustione dell\'idrogeno',
    box: 16, T: 2800, add: [['[H][H]', 24], ['O=O', 12]], photons: { lambda: 200, count: 4, pair: [8, 8] },
    text: "Miscela H₂/O₂ con impulsi meccanici iniziali. Registra le variazioni di connettività; prodotti e velocità dipendono dal potenziale qualitativo, non da una cinetica di combustione validata.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
  },
  {
    id: 'methane', name: 'Combustione del metano',
    box: 16, T: 3200, add: [['C', 8], ['O=O', 16]], spark: { radius: 6, T: 9000 },
    text: "Miscela CH₄/O₂ riscaldata localmente. Confronta energia e frammenti: il modello non predice quantitativamente il meccanismo di combustione.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
  },
  {
    id: 'h2-cl2', name: 'Idrogeno e cloro alla luce',
    box: 18, T: 1100, add: [['[H][H]', 20], ['ClCl', 20]], light: { lambda: 400, on: true },
    text: "Miscela H₂/Cl₂ con deposito meccanico di energia hc/λ. Non vengono calcolati assorbimento ottico, selettività o rese fotochimiche.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
  },
  {
    id: 'na-cl2', name: 'Sodio e cloro: si forma il sale',
    box: 16, T: 900, add: [['[Na]', 16], ['ClCl', 8]],
    text: "Interazioni Na/Cl nel potenziale classico. Le cariche sono parziali e i frammenti restano neutri; non è un modello del reticolo ionico del sale.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
    color: 'charge',
  },
  {
    id: 'carbon', name: 'Vapore di carbonio: catene e anelli',
    box: 14, T: 3500, add: [['[C]', 44]],
    text: "Aggregazione del carbonio nel potenziale del progetto. La comparsa di anelli non identifica fullereni o nanotubi stabili reali.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
  },
  {
    id: 'methyl', name: 'Radicali metile → etano',
    box: 18, T: 500, add: [['[CH3]', 24]],
    text: "Esplora gli urti di frammenti CH₃ nel modello classico. Stabilizzazione, barriere e velocità di ricombinazione non sono validate.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
  },
  {
    id: 'haber', name: 'Azoto e idrogeno: perché serve un catalizzatore',
    box: 16, T: 1500, add: [['N#N', 10], ['[H][H]', 30]],
    text: "Confronta la resistenza dei legami nel modello. Nessun catalizzatore o meccanismo Haber–Bosch è implementato; assenza di reazioni in pochi ps non prova inerzia chimica.",
    tips: ['Confronta NVE e NVT e dimezza Δt per verificare la stabilità numerica.'],
  },
  {
    id: 'aimd-h3', name: 'Ab initio: H + H₂ oltre la barriera', forceField: 'hf', box: 10, T: 300, dt: 0.2, speed: 2,
    atoms: [
      { Z: 1, pos: [-3.2, 0, 0], vel: [0.1433, 0, 0] },
      { Z: 1, pos: [-0.355, 0, 0], vel: [-0.0717, 0, 0] },
      { Z: 1, pos: [0.355, 0, 0], vel: [-0.0717, 0, 0] },
    ],
    text: 'Dinamica ab initio: a ogni passo si risolvono le equazioni di Hartree–Fock e le forze sono il gradiente dell\'energia. L\'atomo di idrogeno arriva con 1,6 eV di energia relativa, più della barriera (1,0 eV in HF/STO-3G; 0,42 eV sperimentale): si forma una nuova molecola H₂ e l\'altro atomo se ne va.',
    tips: ['Guarda gli ordini di legame di Mayer: al passaggio per lo stato di transizione valgono circa 0,5 e 0,5.', 'Prova il preset "sotto la barriera": l\'atomo rimbalza.'],
  },
  {
    id: 'aimd-h3-slow', name: 'Ab initio: H + H₂ sotto la barriera', forceField: 'hf', box: 10, T: 300, dt: 0.2, speed: 2,
    atoms: [
      { Z: 1, pos: [-3.2, 0, 0], vel: [0.0878, 0, 0] },
      { Z: 1, pos: [-0.355, 0, 0], vel: [-0.0439, 0, 0] },
      { Z: 1, pos: [0.355, 0, 0], vel: [-0.0439, 0, 0] },
    ],
    text: 'Stessa collisione con 0,6 eV: l\'energia non basta a superare la barriera di attivazione, l\'atomo rimbalza e la molecola resta vibrando. È l\'origine microscopica dell\'energia di attivazione dell\'equazione di Arrhenius.',
    tips: [],
  },
  {
    id: 'aimd-h2', name: 'Ab initio: rottura del legame di H₂', forceField: 'hf', box: 10, T: 300, dt: 0.2, speed: 2,
    atoms: [
      { Z: 1, pos: [-0.37, 0, 0], vel: [-0.245, 0, 0] },
      { Z: 1, pos: [0.37, 0, 0], vel: [0.245, 0, 0] },
    ],
    text: 'H₂ riceve 6,2 eV di energia vibrazionale, più dell\'energia di dissociazione: il legame si allunga e si rompe. Per descrivere due atomi separati serve la funzione d\'onda a spin non ristretto (UHF): gli elettroni α e β si localizzano su atomi diversi.',
    tips: ['Osserva ⟨S²⟩ nel pannello: da 0 (singoletto puro) sale verso 1 quando il legame si rompe.'],
  },
  {
    id: 'aimd-water', name: 'Ab initio: vibrazioni dell\'acqua', forceField: 'hf', box: 10, T: 1500, dt: 0.25, thermalize: 1500, speed: 4,
    atoms: [
      { Z: 8, pos: [0, 0, 0.07] },
      { Z: 1, pos: [0.76, 0, -0.52] },
      { Z: 1, pos: [-0.76, 0, -0.52] },
    ],
    text: 'Una molecola d\'acqua isolata con energia termica di 1500 K. Stiramenti e piegamento si mescolano; le cariche di Mulliken e gli ordini di legame di Mayer oscillano con la geometria.',
    tips: ['Colora per carica parziale.'],
    color: 'charge',
  },
  {
    id: 'aimd-o2', name: 'Ab initio: O₂ è un tripletto',
    forceField: 'hf', box: 10, T: 300, dt: 0.25, thermalize: 600, speed: 4,
    atoms: [{ Z: 8, pos: [-0.6, 0, 0] }, { Z: 8, pos: [0.6, 0, 0] }],
    text: 'Il calcolo confronta singoletto e tripletto e sceglie il più stabile: per O₂ è il tripletto, con due elettroni spaiati negli orbitali π* (regola di Hund). È il motivo per cui l\'ossigeno liquido è attratto da un magnete.',
    tips: [],
  },
  { id: 'empty', name: 'Scatola vuota', box: 20, T: 300, add: [], text: 'Aggiungi atomi e molecole dal pannello a sinistra, oppure con lo strumento "aggiungi" cliccando nella scatola.', tips: [] },
];
