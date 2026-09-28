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

export const SANDBOX_ELEMENTS = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 35, 36, 53, 54];

// Nomi delle specie riconosciute durante la simulazione (formula come la scrive hillFormula).
export const SPECIES_NAMES = {
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

// Esperimenti: ogni riga di "add" è [SMILES o simbolo, numero di copie].
export const PRESETS = [
  {
    id: 'argon', name: 'Gas di argon: Maxwell–Boltzmann',
    box: 42, T: 300, dt: 2, add: [['[Ar]', 120]],
    text: 'Un gas quasi ideale. Le velocità degli atomi si distribuiscono secondo Maxwell–Boltzmann (grafico centrale). La pressione, misurata come forza media sulle pareti, rispetta PV = nRT: il fattore Z resta vicino a 1.',
    tips: ['Cambia la temperatura: la curva si allarga e la pressione cresce in proporzione.', 'Riduci il volume: la pressione sale (legge di Boyle).'],
  },
  {
    id: 'argon-liquid', name: 'Argon: condensazione e cristallo',
    box: 22, T: 60, dt: 2, add: [['[Ar]', 160]],
    text: 'Sotto la temperatura di ebollizione (87 K) le attrazioni di van der Waals vincono l\'agitazione termica: gli atomi si raccolgono in una goccia di liquido. Z scende molto sotto 1.',
    tips: ['Porta la temperatura a 30 K: la goccia solidifica in un cristallo.', 'Riscalda a 150 K: il liquido evapora.'],
  },
  {
    id: 'water', name: 'Acqua: legami a idrogeno',
    box: 15.5, T: 300, add: [['O', 40]],
    text: 'Le cariche parziali (O negativo, H positivo, dall\'elettronegatività) fanno attrarre le molecole d\'acqua tra loro. I legami a idrogeno O–H···O (linee tratteggiate) si formano e si rompono in continuazione, ogni pochi picosecondi.',
    tips: ['Colora per carica per vedere δ− e δ+.', 'Scalda a 3000 K: l\'acqua si dissocia in H e OH.'],
    color: 'charge',
  },
  {
    id: 'h2-o2', name: 'Combustione dell\'idrogeno',
    box: 16, T: 2800, add: [['[H][H]', 24], ['O=O', 12]], photons: { lambda: 200, count: 4, pair: [8, 8] },
    text: 'Un lampo ultravioletto spezza alcune molecole di O₂. Gli atomi di ossigeno avviano la reazione a catena radicalica: O + H₂ → OH + H, H + O₂ → HO₂, OH + H₂ → H₂O + H. L\'acqua si accumula.',
    tips: ['Spegni il termostato (sistema isolato): il calore di reazione fa salire la temperatura.', 'Segui il registro delle reazioni e le curve delle specie.'],
  },
  {
    id: 'methane', name: 'Combustione del metano',
    box: 16, T: 3200, add: [['C', 8], ['O=O', 16]], spark: { radius: 6, T: 9000 },
    text: 'CH₄ + 2 O₂ → CO₂ + 2 H₂O passa per decine di passaggi radicalici: CH₃, OH, HO₂, CH₂O, CO… Una scintilla innesca la catena.',
    tips: ['Usa lo strumento "scintilla" per riaccendere la reazione.'],
  },
  {
    id: 'h2-cl2', name: 'Idrogeno e cloro alla luce',
    box: 18, T: 1100, add: [['[H][H]', 20], ['ClCl', 20]], light: { lambda: 400, on: true },
    text: 'La luce viola (400 nm, 3,1 eV = 299 kJ/mol) basta per rompere Cl–Cl (242 kJ/mol) ma non H–H (436 kJ/mol). Gli atomi di cloro avviano la catena Cl + H₂ → HCl + H, H + Cl₂ → HCl + Cl.',
    tips: ['Porta la lunghezza d\'onda a 600 nm: i fotoni non hanno più l\'energia per rompere Cl₂.'],
  },
  {
    id: 'na-cl2', name: 'Sodio e cloro: si forma il sale',
    box: 16, T: 900, add: [['[Na]', 16], ['ClCl', 8]],
    text: '2 Na + Cl₂ → 2 NaCl. La grande differenza di elettronegatività (Δχ = 2,2) rende il legame fortemente ionico: le cariche calcolate si avvicinano a ±1 e le coppie ioniche si aggregano.',
    tips: ['Colora per carica.'],
    color: 'charge',
  },
  {
    id: 'carbon', name: 'Vapore di carbonio: catene e anelli',
    box: 14, T: 3500, add: [['[C]', 44]],
    text: 'Atomi di carbonio isolati a 3500 K formano spontaneamente catene, anelli e reticoli con legami doppi e tripli: è così che nascono fullereni e nanotubi nelle scariche ad arco.',
    tips: ['Raffredda lentamente: le strutture si riorganizzano.'],
  },
  {
    id: 'methyl', name: 'Radicali metile → etano',
    box: 18, T: 500, add: [['[CH3]', 24]],
    text: 'Due radicali si uniscono senza barriera: 2 CH₃· → C₂H₆. Ogni atomo di carbonio completa i suoi quattro legami (ottetto).',
    tips: [],
  },
  {
    id: 'haber', name: 'Azoto e idrogeno: perché serve un catalizzatore',
    box: 16, T: 1500, add: [['N#N', 10], ['[H][H]', 30]],
    text: 'N₂ + 3 H₂ ⇌ 2 NH₃ è esotermica, ma il triplo legame N≡N (945 kJ/mol) non si rompe negli urti: senza catalizzatore non succede nulla, anche a 1500 K. Nel processo Haber–Bosch il ferro abbassa la barriera.',
    tips: ['Prova a rompere N₂ con un fotone UV da 120 nm (10 eV).'],
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
